#!/usr/bin/env python
'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''



from cnes.modules.geoloc.lib.interface import Sensor, RiverTile, PixcvecRiver, PixelCloud
import cnes.modules.geoloc.lib.geoloc as geoloc
from cnes.modules.geoloc.lib.pixc_to_shp import pixc_to_shp
from cnes.modules.geoloc.lib.geoloc import GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE

import numpy as np
import argparse
import os
import os.path
from collections import defaultdict
from numpy import polyfit
import shutil
from netCDF4 import Dataset
from scipy import interpolate


# used to filter out fill value and clearly wrong heights, dead sea is about -410 (for now)
MIN_HEIGHT_VALUE = -800

class GeolocRiver(object):
    def __init__(self, pixc, pixcvec, sensor, rivertile):
        
        self.pixc = pixc       
        self.pixcvec = pixcvec
        self.sensor = sensor
        self.rivertile = rivertile

        # Build a dict of (reach id, node id) to node index in the node list
        # Used to find node index (in the node file) corresponding to the (reach index, nodex index) tuple from the pixel cloud
        self.node_reach_dict = {}
        for rivertile_index in range(self.rivertile.node_indx.size):
            rivertile_reach_indx = self.rivertile.reach_indx[rivertile_index]
            rivertile_node_indx = self.rivertile.node_indx[rivertile_index]
            self.node_reach_dict[(rivertile_reach_indx, rivertile_node_indx)] = rivertile_index

        # Build a dict of reach_id -> list of nodes ids
        # Used to interpolate node heights per reach
        self.reach_dict = defaultdict(list)
        for rivertile_index in range(self.rivertile.node_indx.size):
            rivertile_reach_indx = self.rivertile.reach_indx[rivertile_index]
            rivertile_node_indx = self.rivertile.node_indx[rivertile_index]
            self.reach_dict[rivertile_reach_indx].append(rivertile_node_indx)

    def fit_node_heights(self, plot=False):
        """Linear fit on node heights per reach"""

        # dict of reach id -> linear fit coefficients
        reach_fit = {}

        # For each reach, compute the fit
        for reach_id in self.reach_dict.keys():

            # x variable are the node ids. This assumes that nodes are equidistant
            # to be improved using the along reach coordinate of the node instead
            x = self.reach_dict[reach_id]

            # y variable is height
            # get the node heights for the current reach from the reach dict
            y = np.array([self.rivertile.h_n_ave[self.node_reach_dict[(reach_id, node_indx)]] for node_indx in x])

            # filter bad data in y (nan, -inf, fill values)
            good = (y >= MIN_HEIGHT_VALUE)
            x = np.array(x)[good]
            y = np.array(y)[good]

            # Linear fitting (if enough points)
            reach_str = ""
            if np.sum(good) >= 2:
                a, b = polyfit(x, y, 1)
                reach_fit[reach_id] = (a, b)
                reach_str = "reach %d"%(reach_id)
            else:
                reach_fit[reach_id] = None

            if plot: # for debug/visu
                import matplotlib.pyplot as plt
                f, ax = plt.subplots()
                ax.plot(x, y)

                X = np.linspace(np.min(x), np.max(x))
                Y = a*X+b

                ax.plot(X, Y)
                plt.title(reach_str)
                plt.show(block=True)

        # Reconstruct node heights using the computed coefficients
        for rivertile_index in range(self.rivertile.node_indx.size):
            rivertile_reach_indx = self.rivertile.reach_indx[rivertile_index]
            rivertile_node_indx = self.rivertile.node_indx[rivertile_index]
            if reach_fit[rivertile_reach_indx] is not None:
                a, b = reach_fit[rivertile_reach_indx]
                self.rivertile.h_n_ave_fit[rivertile_index] = a*rivertile_node_indx + b

    def estimate_pixvec_height_for_geoloc(self, interpolate_pixc_between_nodes):
        """
        Compute new heights for each pixcvec point,
        either the associated node height (interpolate_pixc_between_nodes == False),
        or interpolating between the associated node and the next/previous closest (interpolate_pixc_between_nodes == True)
        """

        self.new_height = np.copy(self.pixcvec.height)

        unfound_keys = 0

        # For each pixcvec point
        for point_index in range(self.pixcvec.height.size):
            reach_index = self.pixcvec.reach_index[point_index]
            node_index = self.pixcvec.node_index[point_index]

            # Check if the reach and node of the point exist in the node file
            key = (reach_index, node_index)
            if key not in self.node_reach_dict:
                unfound_keys += 1
            else:
                rivertile_index = self.node_reach_dict[key]

                # If possible, interpolate new point height from the associated node and the next/previous closest
                if interpolate_pixc_between_nodes and (reach_index, node_index-1) in self.node_reach_dict and (reach_index,node_index+1) in self.node_reach_dict:
                    rivertile_index_less_one = self.node_reach_dict[(reach_index,node_index-1)]
                    rivertile_index_plus_one = self.node_reach_dict[(reach_index,node_index+1)]

                    # 's' means curvilinear coordinate (distance along reach)
                    s_pixc = self.pixcvec.along_reach[point_index]
                    s_node = self.rivertile.s[rivertile_index]
                    s_node_previous = self.rivertile.s[rivertile_index_less_one]
                    s_node_after = self.rivertile.s[rivertile_index_plus_one]

                    # logic to handle if point is before/after the associated node
                    if s_node < s_pixc:
                        d1 = s_pixc - s_node
                        d2 = s_node_after- s_pixc
                        h1 = self.rivertile.h_n_ave_fit[rivertile_index]
                        h2 = self.rivertile.h_n_ave_fit[rivertile_index_plus_one]
                    else:
                        d1 = s_pixc - s_node_previous
                        d2 = s_node - s_pixc
                        h1 = self.rivertile.h_n_ave_fit[rivertile_index_less_one]
                        h2 = self.rivertile.h_n_ave_fit[rivertile_index]

                    # Linear interpolation of height
                    if h1 >= MIN_HEIGHT_VALUE and h2 >= MIN_HEIGHT_VALUE:
                        self.new_height[point_index] = h1 + (h2-h1)/(d1+d2)*d1
                    if reach_index == 5382:
                        print ("h1:",h1,' h2:',h2,' d1:',d1,' d2:',d2, 
                               's_prev',s_node_previous,' s_node',s_node,'s_aft',s_node_after,
                               ' new_height:',self.new_height[point_index])
                        
                else:
                    # Simple case: just use the node height to estimate improved_height
                    new_height = self.rivertile.h_n_ave[rivertile_index]
                    # avoid bad data
                    if new_height >= MIN_HEIGHT_VALUE:
                        self.new_height[point_index] = new_height
                if reach_index == 5382:
                    print ("new_height:",self.new_height[point_index],'node_height')

        if unfound_keys > 0:
            print("Warning: {} points' (reach, node) were not found in the node file".format(unfound_keys))

    def apply_improved_geoloc(self,method='taylor'):
        """
        Compute the new lat, lon, height using the new heights
        """
        if method == 'taylor':
            self.taylor_improved_geoloc()
        else:
            # swotlib
            self.swotlib_improved_geoloc()

    def swotlib_improved_geoloc(self):
        """
        wrapper for functions that use the actual geolocation in swotlib
        """
        # import python/swotlib-only stuff
        import swot.proc.geolocate
        import swot.proc.base_classes
        import swot.multilook.multilook as ml
        import swot.refloc.refloc
        #print ("got here", self.sensor)
        # make swot.proc sensor file from damiens version
        sensor = swot.proc.base_classes.TVP()#Sensor()
        sensor.time = np.array(self.sensor.time, dtype=np.double)
        sensor.x = np.array(self.sensor.nadir_x, dtype=np.double)
        sensor.y = np.array(self.sensor.nadir_y, dtype=np.double)
        sensor.z = np.array(self.sensor.nadir_z, dtype=np.double)
        sensor.vx = np.array(self.sensor.nadir_vx, dtype=np.double)
        sensor.vy = np.array(self.sensor.nadir_vy, dtype=np.double)
        sensor.vz = np.array(self.sensor.nadir_vz, dtype=np.double)
        sensor.ref_leverarm_x = np.array(self.sensor.ref_leverarm_x, dtype=np.double)
        sensor.ref_leverarm_y = np.array(self.sensor.ref_leverarm_y, dtype=np.double)
        sensor.ref_leverarm_z = np.array(self.sensor.ref_leverarm_z, dtype=np.double)
        sensor.sec_leverarm_x = np.array(self.sensor.sec_leverarm_x, dtype=np.double)
        sensor.sec_leverarm_y = np.array(self.sensor.sec_leverarm_y, dtype=np.double)
        sensor.sec_leverarm_z = np.array(self.sensor.sec_leverarm_z, dtype=np.double)
        #
        #sensor.baseline_left_x = np.array(self.sensor.baseline_left_x,dtype=np.double)
        #sensor.baseline_left_y = np.array(self.sensor.baseline_left_y,dtype=np.double)
        #sensor.baseline_left_z = np.array(self.sensor.baseline_left_z, dtype=np.double)
        #sensor.baseline_right_x = np.array(self.sensor.baseline_right_x, dtype=np.double)
        #sensor.baseline_right_y = np.array(self.sensor.baseline_right_y, dtype=np.double)
        #sensor.baseline_right_z = np.array(self.sensor.baseline_right_z ,dtype=np.double)
        #print ("#########sensor:",sensor.x)
        # get the azimuth looks
        azimuth_looks = int(np.nanmedian(self.pixc.num_rare_looks))
        #print ('#########az_looks:',azimuth_looks)
        line_to_sensor = np.arange(len(sensor.x))
        # fake these so we dont have to read them in and mess with the interface
        #sensor.time = np.array(line_to_sensor,dtype=np.double)
        #lat, lon, height = geoloc.convert_ecef2llh(
        #    sensor.x, sensor.y, sensor.z, GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE)
        #sensor.latitude = np.zeros(np.shape(sensor.x))
        #sensor.longitude = np.zeros(np.shape(sensor.x))
        #sensor.altitude =  np.zeros(np.shape(sensor.x))
        #sensor.heading = np.zeros(np.shape(sensor.x))
        
        rare_line_to_sensor = ml.boxcar_downsample(
            line_to_sensor[np.newaxis].T, azimuth_looks, 1,
            agg_type='sample')
        
        # get the swath side
        swath_side = 'Left'
        swath_side_int = 1
        if self.pixc.tile_ref.endswith('R'):
            swath_side = 'Right'
            swath_side_int = -1
        # make a refloc object
        attributes = {}#pixc.attributes
        attributes['nr_pixels'] = self.pixc.nr_pixels
        attributes['nr_lines'] = self.pixc.nr_lines
        attributes['wavelength'] = self.pixc.wavelength
        attributes['range_spacing'] = self.pixc.range_spacing
        attributes['near_range'] = self.pixc.near_range
        attributes['swath_side_flag'] = swath_side_int
        refloc = swot.refloc.refloc.refloc_from_sensor(
            sensor, attributes, rare_line_to_sensor)
        
        # set xyz image
        # Convert geodetic coordinates (lat, lon, height) to cartesian coordinates (x, y, z)
        x, y, z = geoloc.convert_llh2ecef(
            self.pixcvec.latitude, 
            self.pixcvec.longitude, 
            self.pixcvec.height, 
            GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE)
        #
        X = np.zeros((attributes['nr_lines'],attributes['nr_pixels']))
        Y = np.zeros((attributes['nr_lines'],attributes['nr_pixels']))
        Z = np.zeros((attributes['nr_lines'],attributes['nr_pixels']))
        X[self.pixcvec.azimuth_idx,self.pixcvec.range_idx] = x
        Y[self.pixcvec.azimuth_idx,self.pixcvec.range_idx] = y
        Z[self.pixcvec.azimuth_idx,self.pixcvec.range_idx] = z
        refloc.set_xyz(X,Y,Z)
        
        # set s_image
        use_illumination_time = True
        if use_illumination_time:
            # use the s_image in the pixc file
            s_image = np.zeros((attributes['nr_lines'],attributes['nr_pixels']))
                                    
            #### OLD WAY
            #~ sensor_s = [np.argmin(abs(self.sensor.time 
                                      #~ - self.pixc.illumination_time[k])) 
                        #~ for k in range(len(self.pixc.illumination_time))]
                        
            f = interpolate.interp1d(self.sensor.time,range(len(self.sensor.time)))
            sensor_s = (np.rint(f(self.pixc.illumination_time))).astype(np.double)
                        
            s_image[self.pixc.azimuth_index,self.pixc.range_index] = sensor_s[:]
            refloc.set_s_prof(s_image)
        else:
            # recompute s-image using pixc geoloc as refloc
            refloc.compute_s_prof(0, attributes['nr_lines'])
        #refloc.compute_s_prof(0, attributes['nr_lines'])
        # make H_new as a 2d image
        H_new = np.zeros((attributes['nr_lines'],attributes['nr_pixels']))
        H_new[self.pixcvec.azimuth_idx,self.pixcvec.range_idx] = self.new_height[:]
        
        # finally call the function that does the work 
        #llh = swot.proc.geolocate.height_constrained_geoloc_grad_search (
        #    refloc,
        #    sensor,
        #    H_new,
        #    self.pixc.wavelength,
        #    swath_side_int,
        #    maxiter=10)
        xyz = swot.proc.geolocate.fixed_height_geolocation (
            refloc,
            sensor,
            H_new,
            self.pixc.wavelength)
        llh = swot.proc.geolocate.xyz2llh(xyz)
        # map to outputs
        lat = np.squeeze(llh[0,:,:])
        lon = np.squeeze(llh[1,:,:])
        height = np.squeeze(llh[2,:,:])
        
        msk = np.zeros((attributes['nr_lines'],attributes['nr_pixels']))
        msk[self.pixcvec.azimuth_idx,self.pixcvec.range_idx] = 1
        self.OUT_lat_corr = lat[msk==1]
        self.OUT_lon_corr = lon[msk==1]
        self.OUT_height_corr = height[msk==1]

    def taylor_improved_geoloc(self):
        nb_pix = self.pixcvec.height.size
        #print ("got here", self.sensor)
        # Convert geodetic coordinates (lat, lon, height) to cartesian coordinates (x, y, z)
        x, y, z = geoloc.convert_llh2ecef(self.pixcvec.latitude, self.pixcvec.longitude, self.pixcvec.height, GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE)

        # Get position of associated along-track pixels (in cartesian coordinates)
        nadir_x = self.sensor.nadir_x
        nadir_y = self.sensor.nadir_y
        nadir_z = self.sensor.nadir_z

        # Get velocity of associated along-track pixels (in cartesian coordinates)
        nadir_vx = self.sensor.nadir_vx
        nadir_vy = self.sensor.nadir_vy
        nadir_vz = self.sensor.nadir_vz

        # Get distance from satellite to target point
        
        ri = self.pixcvec.near_range + self.pixcvec.range_idx * self.pixcvec.range_spacing

        # Init output vectors
        self.OUT_lat_corr = np.zeros((nb_pix)) # Improved latitudes
        self.OUT_lon_corr = np.zeros((nb_pix)) # Improved longitudes
        self.OUT_height_corr = np.zeros((nb_pix)) # Improved heights

        # need to remap illumnation time to nearest sensor index
        # TODO replace this by a call to a get_sensor_index or equivalent function
        # that either interpolates the sensor or does something more efficient
        
        #### OLD WAY
        #~ sensor_s = [np.argmin(abs(self.sensor.time 
                                  #~ - self.pixc.illumination_time[k])) 
                    #~ for k in range(len(self.pixc.illumination_time))]
                    
                    
        f = interpolate.interp1d(self.sensor.time,range(len(self.sensor.time)))
        sensor_s = (np.rint(f(self.pixc.illumination_time))).astype(int)
        
        # could call swot.proc.geolocate.height_constrained_geoloc directly here
        # Loop over each pixel (could be vectorized)
        
        # test with vectorisation ##
        h_noisy = self.pixcvec.height 
        nadir_x_vect = np.zeros(nb_pix)
        nadir_y_vect = np.zeros(nb_pix)
        nadir_z_vect = np.zeros(nb_pix)
        nadir_vx_vect = np.zeros(nb_pix)
        nadir_vy_vect = np.zeros(nb_pix)
        nadir_vz_vect = np.zeros(nb_pix)
        
        for i in np.arange(nb_pix):
            ind_sensor = sensor_s[i]
            nadir_x_vect[i] = nadir_x[ind_sensor]
            nadir_y_vect[i] = nadir_y[ind_sensor]
            nadir_z_vect[i] = nadir_z[ind_sensor]
            nadir_vx_vect[i] = nadir_vx[ind_sensor]
            nadir_vy_vect[i] = nadir_vy[ind_sensor]
            nadir_vz_vect[i] = nadir_vz[ind_sensor]
         
        p_final, p_final_llh, h_mu, (iter_grad,nfev_minimize_scalar) = geoloc.pointcloud_height_geoloc_vect(np.transpose(np.array([x, y, z])), h_noisy,
                                                                   np.transpose(np.array([nadir_x_vect, nadir_y_vect, nadir_z_vect])),
                                                                   np.transpose(np.array([nadir_vx_vect, nadir_vy_vect, nadir_vz_vect])),
                                                                   ri, self.new_height, 
                                                                   recompute_Doppler=True, recompute_R=True, verbose = False, 
                                                                   max_iter_grad=1, height_goal = 1.e-3, safe_flag=True)
                                                                   
        self.OUT_lat_corr, self.OUT_lon_corr, self.OUT_height_corr = np.transpose(p_final_llh)

        
def geoloc_river(pixc, pixcvec, sensor, rivertile, fit_heights_per_reach, interpolate_pixc_between_nodes,method='taylor'):
        geolocriver = GeolocRiver(pixc, pixcvec, sensor, rivertile)

        # Do the improved river geolocation
        print("Improved geolocation")
        if fit_heights_per_reach:
            geolocriver.fit_node_heights()
        geolocriver.estimate_pixvec_height_for_geoloc(interpolate_pixc_between_nodes)
        geolocriver.apply_improved_geoloc(method=method)
        
        return (geolocriver.OUT_lat_corr, geolocriver.OUT_lon_corr, geolocriver.OUT_height_corr)

def get_pixc_index(pixc_azimuth_idx,pixc_range_idx,pixcvec, pixc):
    # make a 2d slant image array
    arr2d = np.zeros((pixc.nr_lines,pixc.nr_pixels))
    # assign different value to the pixels that are in the pixcvec vs only pixc
    arr2d[pixc_azimuth_idx,pixc_range_idx] = 1
    arr2d[pixcvec.azimuth_idx,pixcvec.range_idx] = 2
    # keep only the good pixels in pixc 
    arr1d = arr2d[arr2d>0]
    N = len(arr1d)
    indx = np.arange(N,dtype='int32')
    pixc_index = indx[arr1d==2]
    return pixc_index


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pixc', help='full pixel cloud file', type=str)
    parser.add_argument('pixcvec', help='pixcvec file', type=str)
    parser.add_argument('sensor', help='Sensor file', type=str)
    parser.add_argument('rivertile', help='river tile file', type=str)
    parser.add_argument('output_dir', help='output directory', type=str)
    parser.add_argument("--write_shapefile", help="True to write a shapefile version of pixc input and output", nargs='?', type=bool, default=False, const=True)
    parser.add_argument("--fit_heights_per_reach", help="fit and interpolate node heights in each reach", nargs='?', type=bool, default=False, const=True)
    parser.add_argument("--interpolate_pixc_between_nodes", help="interpolate pixc heights between the two closest nodes", nargs='?', type=bool, default=False, const=True)
    args = parser.parse_args()

    # Read inputs
    objPixc = PixcvecRiver(args.pixc)
    objPixcvec = PixcvecRiver(args.pixcvec)
    objSens = Sensor().from_file(args.sensor)
    objRiver = RiverTile().from_file(args.rivertile)

    pixc_dataset = Dataset(args.pixc)

    # Do the improved river geolocation
    (lat_corr, lon_corr, height_corr) = geoloc_river(
            objPixcvec, objSens, objRiver, 
            args.fit_heights_per_reach, 
            args.interpolate_pixc_between_nodes)

    # get the pixc_index
    pixc_index = get_pixc_index(
        objPixc.azimuth_idx,objPixc.range_idx,objPixcvec, objPixc)
    
    # Write outputs

    # format cycle_number and pass_number with leading zeros
    cycle_number = "{:03d}".format(pixc_dataset.getncattr("cycle_number"))
    pass_number = "{:03d}".format(pixc_dataset.getncattr("pass_number"))
    tile_ref = pixc_dataset.getncattr("tile_ref")

    # Symlinks to pixc and sensor with the right name for locnes input
    def swot_symlink(target, name):
        outname = name.format(cycle_number, pass_number, tile_ref)
        outpath = os.path.join(os.path.abspath(args.output_dir), outname)

        # Overwrite the link is existing
        if os.path.islink(outpath):
            print("Overwritting existing {}".format(outpath))
            os.remove(outpath)

        os.symlink(target, outpath)

    swot_symlink(os.path.abspath(args.pixc), "SWOT_L2_HR_PIXC_{}_{}_{}_main.nc")
    swot_symlink(os.path.abspath(args.sensor), "SWOT_L2_HR_PIXC_{}_{}_{}_sensor.nc")

    # write output pixcvec
    pixcvec_outname = "SWOT_L2_HR_PIXC_VEC_RIVER_{}_{}_{}.nc".format(cycle_number, pass_number, tile_ref)
    pixcvec_outpath = os.path.join(os.path.abspath(args.output_dir), pixcvec_outname)
    print("Writing {}".format(pixcvec_outpath))
    shutil.copyfile(args.pixcvec, pixcvec_outpath)

    pixcvec_improved = PixcvecRiver(pixcvec_outpath)
    pixcvec_improved.set_index_file(lat_corr, lon_corr, height_corr, pixc_index)

    if args.write_shapefile:
        print("Writing output shapefile")
        root, ext = os.path.splitext(pixcvec_outpath)
        pixc_to_shp(pixcvec_outpath, root + ".shp", "latitude_vectorproc", "longitude_vectorproc", ["height_vectorproc"])


if __name__ == "__main__":
    main()
