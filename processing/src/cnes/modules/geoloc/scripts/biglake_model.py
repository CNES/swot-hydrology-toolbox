# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2019/05/17:version initiale.
# FIN-HISTORIQUE
# ======================================================
'''
.. module:: biglake_model.py

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


import numpy as np
import scipy.interpolate
import utm
import pyproj
import os


class BigLakeModel(object):
    """
        class BigLakeModel
    """
    def __init__(self, height_model):
        """
        Constructor of BigLakeModel
         :param height_model: height model
         :type height_model: string
        Variables of the object:
        - height_model / string: biglake height model
        """
        self.height_model = height_model
        
    def fit_biglake_model(self, pixc, pixc_mask, grid_spacing, grid_resolution):
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


        :param pixc: pixel cloud from which to compute lake products
        :type in_obj_pixc: proc_pixc_sp.PixelCloudSP object
        :param pixc_mask: selected indices
        :type pixc_mask: numpy 1D array of integer
        :param grid_spacing spacing between averaging points in biglake model
        :type grid_spacing: float
        :param grid_resolution: side length of square averaging box in biglake model
        :type grid_resolution: float

        return: new smooth heights for each
        :rtype: numpy 1D array of float
        """

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

        # Setup grid bounds in grid coordinates
        index_min_x, index_max_x = int((min_x - x_c)/grid_spacing - 0.5), int((max_x - x_c)/grid_spacing + 0.5)
        index_min_y, index_max_y = int((min_y - y_c)/grid_spacing - 0.5), int((max_y - y_c)/grid_spacing + 0.5)


        # Points and values to be interpolated for the height model
        points = [] # mean range_index and azimuth index for each grid point
        values = [] # mean height for each grid point

        # For each grid point
        for index_x in np.arange(index_min_x, index_max_x + 1):
            for index_y in np.arange(index_min_y, index_max_y + 1):
                # The lake center and grid_spacing define our local grid
                utm_index = np.array([(x_c + index_x*grid_spacing), (y_c + index_y*grid_spacing)])

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

        # interpolate new heights for each pixel in slant plane
        xi = np.vstack((pixc.range_idx[pixc_mask], pixc.azimuth_idx[pixc_mask])).T

        # Cubic interpolation and extrapolate using nearest interpolation on the border
        res_linear = scipy.interpolate.griddata(np.array(points), np.array(values), xi, method="cubic", fill_value='nan')

        outside = np.isnan(res_linear)

        # TODO could speed this up by only interpolating where necessary (where xi[outside])
        res_nearest = scipy.interpolate.griddata(np.array(points), np.array(values), xi, method="nearest")

        # replace nans from the cubic interpolation with the result of nearest interpolation
        res = np.where(outside, res_nearest, res_linear)


        return res

    def fit_biglake_model_polyfit(self, pixc, pixc_mask, classif):
        """
        Compute new point heights using a 'big lake' model.

        First, the model creates a grid of spacing 'grid_spacing', centered
        at the lake centroid.
        Then we perform a 2D polynomial fit of the heights on the grid on the ground
        Finally, new heights are interpolated for each pixel cloud coordinate (range, azimuth).

        :param pixc: pixel cloud from which to compute lake products
        :type pixc: proc_pixc.PixelCloud
        :param pixc_mask: selected indices
        :type pixc_mask: numpy 1D array of integer
        :param classif: classification flags
        :type classif: dictionnary

        return: new smooth heights for each pixel
        :rtype : 1D array
        """

        # Center of lake in UTM and find the zone_number, zone_letter
        latitude = pixc.latitude[pixc_mask]
        longitude = pixc.longitude[pixc_mask]
        Z = pixc.height[pixc_mask]

        # TBD : replace by a threshold od dark water %
        if classif["water"] is not None:
            good_ind = np.where(np.isnan(latitude[classif["water"]]) == False)

        else:
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
        A = np.array([X[good_ind]*0+1, X[good_ind], Y[good_ind], X[good_ind]**2, X[good_ind]**2*Y[good_ind], X[good_ind]**2*Y[good_ind]**2,
                      Y[good_ind]**2, X[good_ind]*Y[good_ind]**2, X[good_ind]*Y[good_ind]]).T
        B = Z[good_ind].flatten()
        coeff, r, rank, s = np.linalg.lstsq(A, B)
        
        res = coeff[0]*(X*0+1) + coeff[1]*X + coeff[2]*Y + coeff[3]*(X**2) + coeff[4]*(X**2*Y) + coeff[5]*(X**2*Y**2) + coeff[6]*(Y**2) + \
              coeff[7]*(X*Y**2) + coeff[8]*(X*Y)
        
        return res
    
    
    
    
    
