#!/usr/bin/python2.7
#-*- coding: utf-8 -*-
"""
.. module find_orbit.py
    :synopsis: handle orbit files
    Created on 21 sept. 2015

.. moduleauthor: Capgemini

 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
 
"""

from netCDF4 import Dataset
import numpy as np
import os
from shapely.geometry import box, Polygon

from ressources.utils.inversion_algo import inversionCore
import ressources.utils.vincenty_direct_formula as vincenty


GEN_RAD_EARTH = 6378137.0
GEN_RAD_EARTH_POLE = 6356752.31425

RECORD_MARGIN = 3  # Margin between 2 orbit points
SWATH_MARGIN = 0  # Margin (in meters) to add at the end of the swath


class findOrbit(object):

    def __init__(self, in_north, in_south, in_east, in_west, in_near_range, in_swath_width):
        """
        Init orbit caracteristics
        
        :param in_south: Southern latitude of the studied area
        :type in_south: float
        :param in_north: Northern latitude of the studied area
        :type in_north: float
        :param in_west: Western latitude of the studied area
        :type in_west: float
        :param in_east: Eastern longitude of the studied area
        :type in_east: float
        :param in_swath_width: Swath width
        :type in_swath_width: float
        :param in_near_range: NR cross track
        :type in_near_range: float
        """
        print("[findOrbit] == INIT ==")
        
        # Studied area
        self.north_lat = in_north
        self.south_lat = in_south
        self.east_lon = in_east
        self.west_lon = in_west
        
        # Simulation caracteristics
        self.near_range = in_near_range
        self.swath_width = in_swath_width
    
    #----------------------------------

    def orbit_over_dem(self, in_orbit_directory, in_file_prefix, in_azimuth_spacing, in_swath_width, in_mission_name="SWOT", in_mission_start_time="0000_00_00"):
        """
        Extract parts of orbit from input files, that cover the studied area.
        Output these parts with the sampling specified in configuration file.
        
        :param in_orbit_directory: directory of input orbit files
        :type in_orbit_directory: str
        :param in_file_prefix: prefix for output files (include full path)
        :type in_file_prefix: str
        :param in_azimuth_spacing: azimuth spacing for output file, used to interpolate input orbit files
        :type in_azimuth_spacing: float
        :param in_swath_width: swath width
        :type in_swath_width: float
        :param in_mission_name: mission name (default=SWOT), used in case of specific processing
        :type in_mission_name: str
        :param in_mission_start_time: mission start time
        :type in_mission_start_time: str
        
        :return: out_cycle_duration = cycle duration, read from input orbit files
        :rtype: float
        """
        print("[findOrbit] == orbit_over_dem ==")
        
        # DEM reference polygon as a shapely.geometry.box
        polygon_ref = box(self.south_lat, self.west_lon, self.north_lat, self.east_lon)
        
        # Find all orbit files in the input directory
        for orbit_file in os.listdir(os.path.expandvars(in_orbit_directory)):
            
            if ~os.path.isdir(orbit_file):  # Don't go down the file tree
                
                index_over_dem = []  # Init list of indices of nadir points corresponding to part of orbit overfliying the studied area

                # Open orbit file and get some variables
                data_orbit = Dataset(os.path.join(os.path.expandvars(in_orbit_directory), orbit_file))
                lat = data_orbit.variables['latitude'][:]
                lon = data_orbit.variables['longitude'][:]
                out_cycle_duration = data_orbit.getncattr('repeat_cycle_period')

                for ind_pt in range(lat[:].size - RECORD_MARGIN):
                    
                    # Calculate angle between range and latitude axe - invert phi_left with 2016 orbites
                    if (lat[ind_pt+RECORD_MARGIN] > lat[ind_pt] and lon[ind_pt+RECORD_MARGIN] > lon[ind_pt]) or (lat[ind_pt+RECORD_MARGIN] < lat[ind_pt] and lon[ind_pt+RECORD_MARGIN] < lon[ind_pt]):
                        phi_left = np.rad2deg(np.arccos(np.abs(lat[ind_pt] - lat[ind_pt+RECORD_MARGIN]) / np.sqrt(pow(lat[ind_pt] - lat[ind_pt+RECORD_MARGIN], 2) + pow(lon[ind_pt] - lon[ind_pt+RECORD_MARGIN], 2)))) - 90
                    else:
                        phi_left = 90 - np.rad2deg(np.arccos(np.abs(lat[ind_pt] - lat[ind_pt+RECORD_MARGIN]) / np.sqrt(pow(lat[ind_pt] - lat[ind_pt+RECORD_MARGIN], 2) + pow(lon[ind_pt] - lon[ind_pt+RECORD_MARGIN], 2))))

                    if phi_left > 0:
                        phi_right = phi_left - 180
                    else:
                        phi_right = phi_left + 180

                    # Swath calculation
                    lat_left_nr_first, lon_left_nr_first, left_deg_nr_first = vincenty.dest_vincenty(lat[ind_pt], lon[ind_pt], phi_left, self.near_range)
                    lat_left_fr_first, lon_left_fr_first, left_deg_fr_first = vincenty.dest_vincenty(lat[ind_pt], lon[ind_pt], phi_left, self.swath_width/2 + SWATH_MARGIN)
                    lat_left_nr_second, lon_left_nr_second, left_deg_nr_second = vincenty.dest_vincenty(lat[ind_pt+RECORD_MARGIN], lon[ind_pt+RECORD_MARGIN], phi_left, self.near_range)
                    lat_left_fr_second, lon_left_fr_second, left_deg_fr_second = vincenty.dest_vincenty(lat[ind_pt+RECORD_MARGIN], lon[ind_pt+RECORD_MARGIN], phi_left, self.swath_width/2 + SWATH_MARGIN)

                    lat_right_nr_first, lon_right_nr_first, right_deg_nr_first = vincenty.dest_vincenty(lat[ind_pt], lon[ind_pt], phi_right, self.near_range)
                    lat_right_fr_first, lon_right_fr_first, right_deg_fr_first = vincenty.dest_vincenty(lat[ind_pt], lon[ind_pt], phi_right, self.swath_width/2 + SWATH_MARGIN)
                    lat_right_nr_second, lon_right_nr_second, right_deg_nr_second = vincenty.dest_vincenty(lat[ind_pt+RECORD_MARGIN], lon[ind_pt+RECORD_MARGIN], phi_right, self.near_range)
                    lat_right_fr_second, lon_right_fr_second, right_deg_fr_second = vincenty.dest_vincenty(lat[ind_pt+RECORD_MARGIN], lon[ind_pt+RECORD_MARGIN], phi_right, self.swath_width/2 + SWATH_MARGIN)

                    polygon_data_left = Polygon([[lat_left_nr_first, lon_left_nr_first], [lat_left_fr_first, lon_left_fr_first], [lat_left_fr_second, lon_left_fr_second], [lat_left_nr_second, lon_left_nr_second]])
                    polygon_data_right= Polygon([[lat_right_nr_first, lon_right_nr_first], [lat_right_fr_first, lon_right_fr_first], [lat_right_fr_second, lon_right_fr_second], [lat_right_nr_second, lon_right_nr_second]])

                    # Save file if intersection with DEM > 0
                    if ((polygon_data_left.intersection(polygon_ref).area > 0 or polygon_data_right.intersection(polygon_ref).area > 0) and (-10 < (lat[ind_pt] - self.south_lat) < 10 and -10 < (lon[ind_pt] - self.east_lon) < 10)):
                        if ind_pt not in index_over_dem:
                            index_over_dem.append(ind_pt)
                        if ind_pt+RECORD_MARGIN < lat[:].size:
                            index_over_dem.append(ind_pt+RECORD_MARGIN)

                if len(index_over_dem) > 1:
                    print("> Orbit file = %s" % orbit_file)
                    
                    # Data sampling
                    nb_sampling_points = int(vincenty.dist_vincenty(lat[index_over_dem[0]], lon[index_over_dem[0]], lat[index_over_dem[-1]], lon[index_over_dem[-1]])/in_azimuth_spacing)
                    print("  Number of sampling points = %d" % nb_sampling_points)

                    # Cut valid files and save in new files
                    if in_mission_name == "SWOT":
                        pass_num = int(orbit_file.split('.')[0].split("_")[-1]) + 332  # Compute pass number wrt SWOT KMLs available on AVISO+ (sept2015-v2)
                        if pass_num > 584:
                            pass_num -= 584
                    else:
                        pass_num = int(orbit_file.split('.')[0].split("_")[-1])
                    out_filename = in_file_prefix + "_cycle_0001_pass_%04d.nc" % pass_num
                    print("  Save as %s" % out_filename)
                    output_orbit_file = Dataset(out_filename, "w", format="NETCDF4")
                    
                    # SWOT only: update time vector to be coherent with new pass number
                    tmp_time = data_orbit.variables['time'][:]
                    if in_mission_name == "SWOT":
                        tmp_time += 1024820.9861689  # = 332/2 (orbit number) * 6173.62​0398608 (nodal period)
                        tmp_ind = np.where(tmp_time > out_cycle_duration)[0]
                        if len(tmp_ind) > 0:
                            tmp_time[tmp_ind] -= out_cycle_duration
                
                    # Dimensions
                    output_orbit_file.createDimension('record', nb_sampling_points)

                    # Variables
                    for v_name, varin in iter(data_orbit.variables.items()):
                        outVar = output_orbit_file.createVariable(v_name, varin.datatype, 'record')
                        outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
                        # Linear regression of variable
                        if v_name == "time":  # Specific consideration of time variable
                            lin_reg = np.polyfit(index_over_dem[:], tmp_time[index_over_dem], 1)
                        else:
                            lin_reg = np.polyfit(index_over_dem[:], varin[index_over_dem], 1)
                        give_output = np.poly1d(lin_reg)
                        output_scale = np.linspace(index_over_dem[0], index_over_dem[-1], nb_sampling_points)
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
                    output_orbit_file.setncattr('repeat_cycle_period', out_cycle_duration)
                    output_orbit_file.setncattr('pass_number', pass_num) 
                    output_orbit_file.setncattr('cycle_number', 1) 
                    output_orbit_file.setncattr('beginning_of_mission_time', 0.)
                    output_orbit_file.setncattr('azimuth_spacing', in_azimuth_spacing)
                    output_orbit_file.setncattr('swath_width', in_swath_width)
                    output_orbit_file.setncattr('release', "select_orbit_cnes")
                    output_orbit_file.setncattr('mission start time', in_mission_start_time)
                    output_orbit_file.setncattr('cycle_duration', out_cycle_duration)
                    output_orbit_file.setncattr('dem south latitude', self.south_lat)
                    output_orbit_file.setncattr('dem north latitude', self.north_lat)
                    output_orbit_file.setncattr('dem west longitude', self.west_lon)
                    output_orbit_file.setncattr('dem east longitude', self.east_lon)

                    # Close output orbit file
                    output_orbit_file.close()
                    
                else:
                    print("> NOT KEPT: orbit file = %s" % orbit_file)
                
                # Close input orbit file
                data_orbit.close()
    
        # Return cycle duration
        return out_cycle_duration
