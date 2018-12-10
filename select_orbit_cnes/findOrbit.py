#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

import os
from netCDF4 import Dataset
import numpy as np
from shapely.geometry import box, Polygon

from ressources.utils.inversion_algo import inversionCore
import ressources.utils.vincenty_direct_formula as vincenty

GEN_RAD_EARTH = 6378137.0
GEN_RAD_EARTH_POLE = 6356752.31425

RECORD_MARGIN = 3 # Define a margin between 2 orbits points
SWATH_MARGIN = 0 # Define (in meters) a margin to add at the end of the swath


class findOrbit(object):

    def __init__(self, north, south, east, west, near_range, swath_length):
        print("[findOrbit] == INIT ==")
        
        self.north_lat = north
        self.south_lat = south
        self.east_lon = east
        self.west_lon = west
        self.near_range = near_range
        self.swath= swath_length

    def orbit_over_dem(self, orbit_directory, file_prefix, start_mission_time = "0000_00_00"):
        print("[findOrbit] == orbit_over_dem ==")
        
        orbit_over_dem = []
        
        # DEM reference polygon
        polygon_ref = box(self.south_lat, self.west_lon, self.north_lat, self.east_lon)
        
        # Find all files in this directory
        for orbit_file in os.listdir(orbit_directory):
            
            if ~os.path.isdir(orbit_file):
                index_over_dem = []

                # Open orbit file
                data_orbit = Dataset(os.path.join(orbit_directory, orbit_file))
                lat= data_orbit.variables['latitude'][:]
                lon= data_orbit.variables['longitude'][:]

                for i in range(lat[:].size - RECORD_MARGIN):
                    
                    # Calculate angle between range and latitude axe - invert phi_left with 2016 orbites
                    if (lat[i+RECORD_MARGIN] > lat[i] and lon[i+RECORD_MARGIN] > lon[i]) or (lat[i+RECORD_MARGIN] < lat[i] and lon[i+RECORD_MARGIN] < lon[i]):
                        phi_left = np.rad2deg(np.arccos(np.abs(lat[i] - lat[i+RECORD_MARGIN]) / np.sqrt(pow(lat[i] - lat[i+RECORD_MARGIN], 2) + pow(lon[i] - lon[i+RECORD_MARGIN], 2)))) - 90
                    else:
                        phi_left = 90 - np.rad2deg(np.arccos(np.abs(lat[i] - lat[i+RECORD_MARGIN]) / np.sqrt(pow(lat[i] - lat[i+RECORD_MARGIN], 2) + pow(lon[i] - lon[i+RECORD_MARGIN], 2))))

                    if phi_left > 0:
                        phi_right = phi_left - 180
                    else:
                        phi_right = phi_left + 180

                    # Swath calculation
                    lat_left_nr_first, lon_left_nr_first, left_deg_nr_first = vincenty.dest_vincenty(lat[i], lon[i], phi_left, self.near_range)
                    lat_left_fr_first, lon_left_fr_first, left_deg_fr_first = vincenty.dest_vincenty(lat[i], lon[i], phi_left, self.swath/2 + SWATH_MARGIN)
                    lat_left_nr_second, lon_left_nr_second, left_deg_nr_second = vincenty.dest_vincenty(lat[i+RECORD_MARGIN], lon[i+RECORD_MARGIN], phi_left, self.near_range)
                    lat_left_fr_second, lon_left_fr_second, left_deg_fr_second = vincenty.dest_vincenty(lat[i+RECORD_MARGIN], lon[i+RECORD_MARGIN], phi_left, self.swath/2 + SWATH_MARGIN)

                    lat_right_nr_first, lon_right_nr_first, right_deg_nr_first = vincenty.dest_vincenty(lat[i], lon[i], phi_right, self.near_range)
                    lat_right_fr_first, lon_right_fr_first, right_deg_fr_first = vincenty.dest_vincenty(lat[i], lon[i], phi_right, self.swath/2 + SWATH_MARGIN)
                    lat_right_nr_second, lon_right_nr_second, right_deg_nr_second = vincenty.dest_vincenty(lat[i+RECORD_MARGIN], lon[i+RECORD_MARGIN], phi_right, self.near_range)
                    lat_right_fr_second, lon_right_fr_second, right_deg_fr_second = vincenty.dest_vincenty(lat[i+RECORD_MARGIN], lon[i+RECORD_MARGIN], phi_right, self.swath/2 + SWATH_MARGIN)

                    polygon_data_left = Polygon([[lat_left_nr_first, lon_left_nr_first], [lat_left_fr_first, lon_left_fr_first], [lat_left_fr_second, lon_left_fr_second], [lat_left_nr_second, lon_left_nr_second]])
                    polygon_data_right= Polygon([[lat_right_nr_first, lon_right_nr_first], [lat_right_fr_first, lon_right_fr_first], [lat_right_fr_second, lon_right_fr_second], [lat_right_nr_second, lon_right_nr_second]])

                    # Save file if intersection with DEM > 0
                    if ((polygon_data_left.intersection(polygon_ref).area > 0 or polygon_data_right.intersection(polygon_ref).area > 0) and (-10 < (lat[i] - self.south_lat) < 10 and -10 < (lon[i] - self.east_lon) < 10)):
                        if i not in index_over_dem:
                            index_over_dem.append(i)
                        if i+RECORD_MARGIN < lat[:].size:
                            index_over_dem.append(i+RECORD_MARGIN)

                if len(index_over_dem) > 1:
                    print("> Orbit file = %s" % orbit_file)
                    
                    # Data sampling
                    nb_sampling_points = int(vincenty.dist_vincenty(lon[index_over_dem[0]], lat[index_over_dem[0]], lon[index_over_dem[-1]], lat[index_over_dem[-1]])/21.875)
                    print("  Number of sampling point = %d" % nb_sampling_points)

                    # Cut valid files and save in new files
                    pass_num = int(orbit_file.split('.')[0].split("_")[-1]) + 1  # Compute pass number wrt SWOT KMLs available on AVISO+
                    if pass_num > 584:
                        pass_num = pass_num % 585 + 1
                    out_filename = file_prefix + "_cycle_0000_pass_%04d.nc" % pass_num
                    print("  Save as %s" % out_filename)
                    output_orbit_file = Dataset(out_filename, "w", format="NETCDF4")
                
                    # Dimensions
                    output_orbit_file.createDimension('record', nb_sampling_points)

                    # Variables
                    for v_name, varin in iter(data_orbit.variables.items()):
                        outVar = output_orbit_file.createVariable(v_name, varin.datatype, 'record')
                        outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
                        # Linear regression of variable
                        lin_reg = np.polyfit(index_over_dem[:], varin[index_over_dem], 1)
                        give_output = np.poly1d(lin_reg)
                        output_scale = np.arange(index_over_dem[0], index_over_dem[-1], float(index_over_dem[-1] - index_over_dem[0]) / float(nb_sampling_points))
                        #print(give_output(output_scale))
                        outVar[:] = give_output(output_scale)
                    
                    # Creating x, y and z variables
                    x, y, z = inversionCore.convert_llh2ecef(output_orbit_file.variables['latitude'][:], output_orbit_file.variables['longitude'][:], output_orbit_file.variables['altitude'][:], GEN_RAD_EARTH, GEN_RAD_EARTH_POLE)
                    outVar = output_orbit_file.createVariable('x', np.float64, 'record')
                    outVar[:] = x[:]
                    
                    outVar = output_orbit_file.createVariable('y', np.float64, 'record')
                    outVar[:] = y[:]

                    outVar = output_orbit_file.createVariable('z', np.float64, 'record')
                    outVar[:] = z[:]
                    # Global attributes
                    output_orbit_file.setncattr('repeat_cycle_period', 1802697.12)
                    output_orbit_file.setncattr('pass_number', orbit_file[-7:]) 
                    output_orbit_file.setncattr('cycle_number', 1) 
                    output_orbit_file.setncattr('beginning_of_mission_time', 0.)
                    output_orbit_file.setncattr('azimuth_spacing', 21.875)
                    output_orbit_file.setncattr('swath_width', 120000.)
                    output_orbit_file.setncattr('release', "Build 1051 made by user bawillia on 10/20/2016 12:32:57")
                    output_orbit_file.setncattr('mission start time', "2014-01-01")
                    output_orbit_file.setncattr('cycle_duration', 1802645.8059698)
                    output_orbit_file.setncattr('dem south latitude', self.south_lat)
                    output_orbit_file.setncattr('dem north katitude', self.north_lat)
                    output_orbit_file.setncattr('dem west longitude', self.west_lon)
                    output_orbit_file.setncattr('dem east longitude', self.east_lon)

                    data_orbit.close()
                    output_orbit_file.close()
                    
                else:
                    print("> NOT KEPT: orbit file = %s" % orbit_file)

        return orbit_over_dem
    