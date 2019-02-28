'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


import numpy as np
import scipy.interpolate
import utm
import pyproj


class BigLakeModel(object):
    
    def __init__(self, height_model):
        self.height_model = height_model
        
    def fit_biglake_model(self, pixc, pixc_mask, grid_spacing, grid_resolution, plot=False):
        """
        Compute new point heights using a 'big lake' model.

        First, the model creates a grid of spacing 'grid_spacing', centered
        at the lake centroid.
        Then, for each point of this grid we consider pixels whose noisy
        ground positions are inside a square box of side length
        grid_resolution around the current grid cell center.
        Then, we compute a mean height, range index and azimuth index for each grid point
        using those pixels.
        Finally, new heights are interpolated for each pixel cloud coordinate (range, azimuth).

        Args:
            pixc: pixel cloud
            pixc_mask: indices of pixels of current water body
            grid_spacing spacing between averaging points in biglake model
            grid_resolution: side length of square averaging box in biglake model
            plot: debugging plots

        returns:
            1D array of new smooth heights for each pixel
        """

        if grid_resolution < grid_spacing:
            raise ValueError("grid resolution is less than grid spacing")

        # Center of lake in UTM
        x_c, y_c, zone_number, zone_letter = utm.from_latlon(
                np.mean(pixc.latitude[pixc_mask]),
                np.mean(pixc.longitude[pixc_mask]))

        # Convert pixel cloud to UTM (zone of the centroid)
        latlon = pyproj.Proj(init="epsg:4326")
        utm_proj = pyproj.Proj("+proj=utm +zone={}{} +ellps=WGS84 +datum=WGS84 +units=m +no_defs".format(zone_number, zone_letter))
        X, Y = pyproj.transform(latlon, utm_proj, pixc.longitude[pixc_mask], pixc.latitude[pixc_mask])
        XY = np.vstack((X, Y)).T # column-variables format

        # Lake bounding box to limit grid extent
        min_x, max_x = np.min(X), np.max(X)
        min_y, max_y = np.min(Y), np.max(Y)

        # The lake center and grid_spacing define our local grid
        def local_to_utm(xi, yi):
            return (x_c + xi*grid_spacing), (y_c + yi*grid_spacing)

        def utm_to_local(x, y):
            return (x - x_c) / grid_spacing, (y - y_c) / grid_spacing

        # Setup grid bounds in grid coordinates
        index_min_x, index_max_x = int((min_x - x_c)/grid_spacing - 0.5), int((max_x - x_c)/grid_spacing + 0.5)
        index_min_y, index_max_y = int((min_y - y_c)/grid_spacing - 0.5), int((max_y - y_c)/grid_spacing + 0.5)

        if plot:
            import matplotlib.pyplot as plt
            f, ax = plt.subplots()

            ax.plot(X[::100], Y[::100], color="k", linestyle="", marker=".")
            ax.plot(x_c, y_c, markersize=5, marker="o", color="r", linestyle="")
            ax.plot((min_x, min_x, max_x, max_x),
                    (min_y, max_y, min_y, max_y),
                    color="r", linestyle="", marker="o")

            ax.plot(*local_to_utm(index_min_x, index_min_y),
                    color="b", linestyle="", marker="o")

            ax.plot(*local_to_utm(index_min_x, index_max_y),
                    color="b", linestyle="", marker="o")

            ax.plot(*local_to_utm(index_max_x, index_min_y),
                    color="b", linestyle="", marker="o")

            ax.plot(*local_to_utm(index_max_x, index_max_y),
                    color="b", linestyle="", marker="o")

        # Points and values to be interpolated for the height model
        points = [] # mean range_index and azimuth index for each grid point
        values = [] # mean height for each grid point

        # For each grid point
        for index_x in np.arange(index_min_x, index_max_x + 1):
            for index_y in np.arange(index_min_y, index_max_y + 1):
                utm_index = np.array(local_to_utm(index_x, index_y))

                # Infinity norm is a square box around utm_index (center of grid point)
                norm = np.linalg.norm(XY - utm_index, ord=np.Inf, axis=1)

                # indices of points within the grid cell (for both grid_spacing and grid_resolution)
                points_in_square_resolution = np.where(norm < grid_resolution / 2.0)[0]
                points_in_square_spacing = np.where(norm < grid_spacing / 2.0)[0]

                # If there are points in the grid square, average their range, azimuth and height
                # TODO improve this criterion?
                if len(np.isfinite(points_in_square_spacing)) > 0:
                    h_mean = np.nanmean(pixc.height[pixc_mask][points_in_square_resolution])
                    range_index_mean = np.nanmean(pixc.range_idx[pixc_mask][points_in_square_spacing])
                    azimuth_index_mean = np.nanmean(pixc.azimuth_idx[pixc_mask][points_in_square_spacing])

                    points.append([range_index_mean, azimuth_index_mean])
                    values.append(h_mean)

                if plot:
                    r, g, b = np.random.rand(3)
                    ax.plot(X[points_in_square_resolution], Y[points_in_square_resolution], marker=".", linestyle="", color=(r, g, b, 0.2))

                    ax.plot(utm_index[0], utm_index[1], marker="x", color="k", linestyle="")

        if plot:
            plt.show(block=False)

        # interpolate new heights for each pixel in slant plane
        xi = np.vstack((pixc.range_idx[pixc_mask], pixc.azimuth_idx[pixc_mask])).T

        # Cubic interpolation and extrapolate using nearest interpolation on the border
        res_linear = scipy.interpolate.griddata(np.array(points), np.array(values), xi, method="cubic", fill_value='nan')

        outside = np.isnan(res_linear)

        # TODO could speed this up by only interpolating where necessary (where xi[outside])
        res_nearest = scipy.interpolate.griddata(np.array(points), np.array(values), xi, method="nearest")

        # replace nans from the cubic interpolation with the result of nearest interpolation
        res = np.where(outside, res_nearest, res_linear)

        if plot:
            image_noisy = np.full((pixc.nb_pix_azimuth, pixc.nb_pix_range), 0.0)
            image_clean = np.full((pixc.nb_pix_azimuth, pixc.nb_pix_range), 0.0)

            for i in range(len(res)):
                if i % 200 == 0:
                    image_noisy[pixc.azimuth_idx[pixc_mask][i], pixc.range_idx[pixc_mask][i]] = pixc.height[pixc_mask][i]
                    image_clean[pixc.azimuth_idx[pixc_mask][i], pixc.range_idx[pixc_mask][i]] = res[i]

            plt.figure()
            plt.imshow(image_noisy, vmin=-3, vmax=-1, origin="lower")

            plt.title("noisy")

            plt.figure()
            plt.imshow(image_clean, vmin=-3, vmax=-1, origin="lower")
            plt.colorbar()

            plt.plot(np.array(points)[:, 0], np.array(points)[:, 1],
                    color="r", linestyle="", marker="+")

            plt.title("clean")

            plt.show(block=True)

        return res

    def fit_biglake_model_polyfit(self, pixc, pixc_mask):


        # Center of lake in UTM and find the zone_number, zone_letter
        
        latitude = pixc.latitude[pixc_mask]
        longitude = pixc.longitude[pixc_mask]
        Z = pixc.height[pixc_mask]

        good_ind = np.where(np.isnan(latitude) == False)

        ind = np.where(longitude > 180.)
        if ind is not None:
            longitude[ind]+= -360.
        
        x_c, y_c, zone_number, zone_letter = utm.from_latlon(
                np.nanmean(latitude),
                np.nanmean(longitude))
                
        # Convert pixel cloud to UTM (zone of the centroid)
        latlon = pyproj.Proj(init="epsg:4326")
        utm_proj = pyproj.Proj("+proj=utm +zone={}{} +ellps=WGS84 +datum=WGS84 +units=m +no_defs".format(zone_number, zone_letter))
        X, Y = pyproj.transform(latlon, utm_proj, longitude, latitude)

        # 2D polynomial fitting of heights on the ground grid 
        # For now, degre 2, may be change
        A = np.array([X[good_ind]*0+1, X[good_ind], Y[good_ind], X[good_ind]**2, X[good_ind]**2*Y[good_ind], X[good_ind]**2*Y[good_ind]**2, Y[good_ind]**2, X[good_ind]*Y[good_ind]**2, X[good_ind]*Y[good_ind]]).T
        B = Z[good_ind].flatten()

        coeff, r, rank, s = np.linalg.lstsq(A, B)
        
        res = coeff[0]*(X*0+1) + coeff[1]*X + coeff[2]*Y + coeff[3]*(X**2) + coeff[4]*(X**2*Y) + coeff[5]*(X**2*Y**2) + coeff[6]*(Y**2) + coeff[7]*(X*Y**2) + coeff[8]*(X*Y)
        
        return res
    
    
    
    
    
