# -*- coding: utf8 -*-
"""
.. module proc_real_pixc.py
    :synopsis: Deal with official pixel cloud (L2_HR_PIXC) files
    Created on 08/24/2017
    2018/11/30 (D. Desroches, V. Poughon - CNES): change variables names wrt to new PixC naming convention

.. module author: Claire POTTIER - CNES DSO/SI/TR

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from __future__ import absolute_import, division, print_function, unicode_literals

from datetime import datetime, timedelta
import numpy as np
import os
from osgeo import ogr, osr

import lib.my_api as my_api
import lib.my_netcdf_file as my_nc
import lib.my_variables as my_var
import mathematical_function as math_fct
from lib.my_variables import GEN_APPROX_RAD_EARTH


def fill_vector_param(variable, variable_name, ref_size, data_param, group=None):
    """
    Fill variable data field
    
    :param variable: input data
    :type variable: array
    :param variable_name: name for the variable in the NetCDF file
    :type variable_name: string
    :param ref_size: number of elements associated to the variable in the NetCDF file
    :type ref_size: int
    :param data_param: pointer to the NetCDF writer (=NcWrite(OUT_file))
    :type data_param: -
    :param group:
    :type group:
    """
    if variable is not None:
        xsize = len(variable)
        if xsize != ref_size:
            exc = '[proc_realPixC/fill_vector_param] ERROR = There is a problem with the size of ' + variable_name
            exit(exc)
        else:
            data_param.fill_variable(variable_name, variable, group=group)
            
            
#################################################


class l2_hr_pixc(object):

    def __init__(self, IN_azimuth_index, IN_range_index, IN_classification, IN_pixel_area, IN_latitude, IN_longitude, IN_height, IN_crosstrack,
                 IN_nadir_time, IN_nadir_latitude, IN_nadir_longitude, IN_nadir_altitude, IN_nadir_heading, IN_nadir_x, IN_nadir_y, IN_nadir_z, IN_nadir_vx, IN_nadir_vy, IN_nadir_vz, IN_nadir_near_range,
                 IN_mission_start_time, IN_cycle_duration, IN_cycle_num, IN_pass_num, IN_tile_ref, IN_nb_pix_range, IN_nb_pix_azimuth, IN_azimuth_spacing, IN_range_spacing, IN_near_range):
        """
        Constructor of the pixel cloud product

        :param IN_azimuth_index: azimuth indices
        :type IN_azimuth_index: 1D-array of int
        :param IN_range_index: range indices
        :type IN_range_index: 1D-array of int
        :param IN_classification: classification values
        :type IN_classification: 1D-array of int
        :param IN_pixel_area: surface area
        :type IN_pixel_area: 1D-array of float
        :param IN_latitude: latitude values
        :type IN_latitude: 1D-array of float
        :param IN_longitude: longitude values
        :type IN_longitude: 1D-array of float
        :param IN_height: height
        :type IN_height: 1D-array of float
        :param IN_crosstrack: crosstrack distance from nadir track
        :type IN_crosstrack: 1D-array of float
            
        :param IN_nadir_time: time tags for each nadir points ; provided as UTC seconds since begin of current cycle
        :type IN_nadir_time: 1D-array of float
        :param IN_nadir_latitude: latitude values
        :type IN_nadir_latitude: 1D-array of float
        :param IN_nadir_longitude: longitude values
        :type IN_nadir_longitude: 1D-array of float
        :param IN_nadir_altitude: altitude values
        :type IN_nadir_altitude: 1D-array of float
        :param IN_nadir_heading: heading values
        :type IN_nadir_heading: 1D-array of float
        :param IN_nadir_x|y|z: x|y|z cartesian coordinate values
        :type IN_nadir_x|y|z: 1D-array of float
        :param IN_nadir_vx|vy|vz: vx|vy|vz cartesian velocity values
        :type IN_nadir_vx|vy|vz: 1D-array of float
        :param IN_nadir_near_range: near range distance for each time tag
        :type IN_nadir_near_range: 1D-array of float
        
        :param IN_mission_start_time: mission start time
        :type IN_mission_start_time: string (yyyy-mm-dd)
        :param IN_cycle_duration: number of seconds in a cycle
        :type IN_cycle_duration: int
        :param IN_cycle_num: cycle number
        :type IN_cycle_num: int
        :param IN_pass_num: pass number
        :type IN_pass_num: int
        :param IN_tile_ref: tile reference
        :type IN_tile_ref: string
        :param IN_nb_pix_range: number of pixels in range of the interferogram
        :type IN_nb_pix_range: int
        :param IN_nb_pix_azimuth: number of pixels in azimuth of the interferogram
        :type IN_nb_pix_azimuth: int
        :param IN_azimuth_spacing: azimuth spacing
        :type IN_azimuth_spacing: float
        :param IN_range_spacing: range spacing
        :type IN_range_spacing: float
        :param IN_near_range: range distance at the near range
        :type IN_near_range: float
            
        + nb_water_pix(int) : number of water pixels, i.e. pixels in azimuth_index, ..., crosstrack vectors
        + nb_nadir_pix(int) : number of pixels on the nadir track, i.e. pixels in time, ..., near_range vectors
        + pattern(str): filename pattern
        + tile_boundaries (list of tuple of float) : tile boundaries [inner_first, outer_first, inner_last, outer_last]
            inner_first : lonlat(rg_min, az_min)
            outer_first : lonlat(rg_max, az_min)
            inner_last : lonlat(rg_min, az_max)
            outer_last : lonlat(rg_max, az_max)
        """
        my_api.printInfo("[proc_real_pixc] == INIT ==") 
        
        self.azimuth_index = IN_azimuth_index
        self.range_index = IN_range_index
        self.classification = IN_classification
        self.pixel_area = IN_pixel_area
        self.latitude = IN_latitude
        self.longitude = IN_longitude
        self.height = IN_height
        self.crosstrack = IN_crosstrack
        self.nb_water_pix = IN_azimuth_index.size
        
        # Modification to have sensor_s (sensor azimuth position for each pixel) to be compatible with HR simulator. It is a duplication of azimuth_index in the large scale simulator
        self.sensor_s = IN_azimuth_index
        self.nadir_time = IN_nadir_time
        
        self.illumination_time = np.zeros(len(IN_azimuth_index))
        
        for i in range(self.illumination_time.size):
            self.illumination_time[i] = self.nadir_time[self.sensor_s[i]]
            
        self.nadir_latitude = IN_nadir_latitude
        self.nadir_longitude = IN_nadir_longitude
        self.nadir_altitude = IN_nadir_altitude
        self.nadir_heading = IN_nadir_heading
        self.nadir_x = IN_nadir_x
        self.nadir_y = IN_nadir_y
        self.nadir_z = IN_nadir_z
        self.nadir_vx = IN_nadir_vx
        self.nadir_vy = IN_nadir_vy
        self.nadir_vz = IN_nadir_vz
        self.nadir_near_range = IN_nadir_near_range
        self.nb_nadir_pix = IN_nadir_time.size
        
        self.mission_start_time = IN_mission_start_time
        self.cycle_duration = IN_cycle_duration
        self.cycle_num = IN_cycle_num
        self.pass_num = IN_pass_num
        self.tile_ref = IN_tile_ref
        self.nb_pix_range = IN_nb_pix_range
        self.nb_pix_azimuth = IN_nb_pix_azimuth
        self.azimuth_spacing = IN_azimuth_spacing
        self.range_spacing = IN_range_spacing
        self.near_range = IN_near_range 

        self.tile_boundaries = []
    #----------------------------------

    def set_tile_boundaries(self, IN_attributes, min_az):
        """
        this function set tile boundaries in order to fill inner/outer first/last longitude/latitude of tile.

        :param IN_attributes:
        :type IN_attributes:
        :param min_az: min azimuth at the begining of the tile
        :type min_az: float
        """
        az = np.array([self.nb_pix_azimuth, self.nb_pix_azimuth, 0, 0]) + min_az
        rg = np.array([0, self.nb_pix_range, 0, self.nb_pix_range])

        if self.tile_ref[-1] == "L" :
            swath = "Left"
        else :
            swath = "Right"

        r0 = np.sqrt((IN_attributes.alt + (IN_attributes.nr_cross_track ** 2) / ( 2 * GEN_APPROX_RAD_EARTH)) ** 2 + IN_attributes.nr_cross_track ** 2)  # Radar to ground distance for near range pixels
        ri = r0[az] + rg * IN_attributes.range_sampling  # Radar-to-ground distance

        lon_list, lat_list = math_fct.lonlat_from_azy(az, ri, IN_attributes, swath, 0, "deg")

        for i, lon in enumerate(lon_list):
            lat = lat_list[i]
            self.tile_boundaries.append((lon, lat))


    def write_pixc_file(self, IN_output_file, compress=False):
        """
        Write the main file of real pixel cloud product (L2_HR_PIXC product, main file)

        :param IN_output_file: output full path
        :type IN_output_file: string
        :param compress: parameter the define to compress or not the file
        :type compress: boolean
        """
        my_api.printInfo("[proc_real_pixc] == write_pixc_file : %s ==" % IN_output_file)
    
        # 1 - Open NetCDF file in writing mode
        data = my_nc.myNcWriter(IN_output_file)
        
        # Global attributes
        data.add_global_attribute('Conventions', 'CF-1.7')
        data.add_global_attribute('title', 'Level 2 KaRIn High Rate Water Mask Pixel Clould Data Product')
        data.add_global_attribute('institution', 'CNES - Large scale simulator')
        data.add_global_attribute('source', 'Ka-band radar interferometer')
        data.add_global_attribute('history', "%sZ: Creation" % datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
        data.add_global_attribute('mission_name', "SWOT")
        data.add_global_attribute('references', 'Large scale simulator')
        data.add_global_attribute('reference_document', 'JPL D-56411 - Initial release - February 11, 2019')
        data.add_global_attribute('contact', 'None')
        data.add_global_attribute('cycle_number', self.cycle_num)
        data.add_global_attribute('pass_number', np.int(self.pass_num))
        data.add_global_attribute('tile_number', int(self.tile_ref[0:-1]))
        data.add_global_attribute('swath_side', self.tile_ref[-1])
        data.add_global_attribute('tile_name', "%03d_%03d%s" % (np.int(self.pass_num), int(self.tile_ref[0:-1]), self.tile_ref[-1]))
        data.add_global_attribute("wavelength", 0.008385803020979)
        data.add_global_attribute('near_range', np.min(self.near_range))  # TODO: improve
        data.add_global_attribute('nominal_slant_range_spacing', self.range_spacing)
        data.add_global_attribute('start_time', self.computeDate(self.nadir_time[0]))    
        data.add_global_attribute('stop_time', self.computeDate(self.nadir_time[-1]))  
        data.add_global_attribute('polarization', 'None')         
        data.add_global_attribute('transmit_antenna', 'None')
        data.add_global_attribute('processing_beamwidth', 'None')

        data.add_global_attribute("inner_first_longitude", self.tile_boundaries[0][0])
        data.add_global_attribute("inner_first_latitude", self.tile_boundaries[0][1])
        data.add_global_attribute("inner_last_longitude", self.tile_boundaries[1][0])
        data.add_global_attribute("inner_last_latitude", self.tile_boundaries[1][1])
        data.add_global_attribute("outer_first_longitude", self.tile_boundaries[2][0])
        data.add_global_attribute("outer_first_latitude", self.tile_boundaries[2][1])
        data.add_global_attribute("outer_last_longitude", self.tile_boundaries[3][0])
        data.add_global_attribute("outer_last_latitude", self.tile_boundaries[3][1])

        data.add_global_attribute("slc_first_line_index_in_tvp", 'None')
        data.add_global_attribute("slc_last_line_index_in_tvp", 'None')
        data.add_global_attribute("xref_input_l1b_hr_slc_file", 'None')
        data.add_global_attribute("xref_input_static_karin_cal_file", 'None')
        data.add_global_attribute("xref_input_ref_dem_file", 'None')
        data.add_global_attribute("xref_input_water_mask_file", 'None')
        data.add_global_attribute("xref_input_static_geophys_file", 'None')
        data.add_global_attribute("xref_input_dynamic_geophys_file", 'None')
        data.add_global_attribute("xref_input_int_lr_xover_cal_file", 'None')
        data.add_global_attribute("xref_l2_hr_pixc_config_parameters_file", 'None')
        data.add_global_attribute("ellipsoid_semi_major_axis", 'None')
        data.add_global_attribute("ellipsoid_flattening", 'None')

        # =======================
        # == Group pixel_cloud ==
        # =======================
        pixc = data.add_group("pixel_cloud")
        
        # Group attributes
        data.add_global_attribute('description', 'cloud of geolocated interferogram pixels', group=pixc)     
        data.add_global_attribute('interferogram_size_azimuth', self.nb_pix_azimuth, group=pixc) 
        data.add_global_attribute('interferogram_size_range', self.nb_pix_range, group=pixc)      
        data.add_global_attribute('looks_to_efflooks', 1.75, group=pixc)   
  
        # Group dimensions
        data.add_dimension('points', self.nb_water_pix, group=pixc)
        data.add_dimension('depth', 2, group=pixc)
        
        # Group variables
        data.add_variable('azimuth_index', np.int32, 'points', my_var.FV_NETCDF["int32"], compress, group=pixc)
        fill_vector_param(self.azimuth_index, 'azimuth_index', self.nb_water_pix, data, group=pixc)
        data.add_variable('range_index', np.int32, 'points', my_var.FV_NETCDF["int32"], compress, group=pixc)
        fill_vector_param(self.range_index, 'range_index', self.nb_water_pix, data, group=pixc)
        
        data.add_variable('interferogram', np.float32, ('points', 'depth'), my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros([self.nb_water_pix, 2]), 'interferogram', self.nb_water_pix, data, group=pixc)
        data.add_variable('power_plus_y', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'power_plus_y', self.nb_water_pix, data, group=pixc)
        data.add_variable('power_minus_y', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'power_minus_y', self.nb_water_pix, data, group=pixc)
        data.add_variable('coherent_power', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'coherent_power', self.nb_water_pix, data, group=pixc)
        data.add_variable('x_factor_plus_y', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'x_factor_plus_y', self.nb_water_pix, data, group=pixc)
        data.add_variable('x_factor_minus_y', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'x_factor_minus_y', self.nb_water_pix, data, group=pixc)  
        
        data.add_variable('water_frac', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.ones(self.nb_water_pix), 'water_frac', self.nb_water_pix, data, group=pixc)       
        data.add_variable('water_frac_uncert', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'water_frac_uncert', self.nb_water_pix, data, group=pixc)              
        data.add_variable('classification', np.int8, 'points', my_var.FV_NETCDF["int8"], compress, group=pixc)
        fill_vector_param(self.classification, 'classification', self.nb_water_pix, data, group=pixc) 
        data.add_variable('false_detection_rate', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'false_detection_rate', self.nb_water_pix, data, group=pixc)
        data.add_variable('missed_detection_rate', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'missed_detection_rate', self.nb_water_pix, data, group=pixc)
        data.add_variable('prior_water_prob', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'prior_water_prob', self.nb_water_pix, data, group=pixc)
        data.add_variable('bright_land_flag', np.int8, 'points', my_var.FV_NETCDF["int8"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'bright_land_flag', self.nb_water_pix, data, group=pixc)          
        data.add_variable('layover_impact', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'layover_impact', self.nb_water_pix, data, group=pixc)
        
        data.add_variable('num_rare_looks', np.int8, 'points', my_var.FV_NETCDF["int8"], compress, group=pixc)
        fill_vector_param(np.full(self.nb_water_pix, 7.), 'num_rare_looks', self.nb_water_pix, data, group=pixc) 
        
        data.add_variable('latitude', np.float64, 'points', my_var.FV_NETCDF["float64"], compress, group=pixc)
        data.add_variable_attribute('latitude', 'units', 'degrees_north', group=pixc)
        fill_vector_param(self.latitude, 'latitude', self.nb_water_pix, data, group=pixc)
        data.add_variable('longitude', np.float64, 'points', my_var.FV_NETCDF["float64"], compress, group=pixc)
        data.add_variable_attribute('longitude', 'units', 'degrees_east', group=pixc)
        fill_vector_param(self.longitude, 'longitude', self.nb_water_pix, data, group=pixc)
        data.add_variable('height', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        data.add_variable_attribute('height', 'units', 'm', group=pixc)
        fill_vector_param(self.height, 'height', self.nb_water_pix, data, group=pixc)
        
        data.add_variable('cross_track', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(self.crosstrack, 'cross_track', self.nb_water_pix, data, group=pixc)
        data.add_variable('pixel_area', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(self.pixel_area, 'pixel_area', self.nb_water_pix, data, group=pixc)        
        data.add_variable('inc', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'inc', self.nb_water_pix, data, group=pixc)
        
        data.add_variable('phase_noise_std', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'phase_noise_std', self.nb_water_pix, data, group=pixc)        
        data.add_variable('dlatitude_dphase', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dlatitude_dphase', self.nb_water_pix, data, group=pixc)
        data.add_variable('dlongitude_dphase', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dlongitude_dphase', self.nb_water_pix, data, group=pixc)                  
        data.add_variable('dheight_dphase', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dheight_dphase', self.nb_water_pix, data, group=pixc)  
        data.add_variable('dheight_droll', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dheight_droll', self.nb_water_pix, data, group=pixc)
        data.add_variable('dheight_dbaseline', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dheight_dbaseline', self.nb_water_pix, data, group=pixc)                  
        data.add_variable('dheight_drange', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dheight_drange', self.nb_water_pix, data, group=pixc)       
        data.add_variable('darea_dheight', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'darea_dheight', self.nb_water_pix, data, group=pixc)
        
        data.add_variable('illumination_time', np.float64, 'points', my_var.FV_NETCDF["float64"], compress, group=pixc)
        fill_vector_param(self.computeTime_UTC(self.illumination_time), 'illumination_time', self.nb_water_pix, data, group=pixc)
        data.add_variable('illumination_time_tai', np.float64, 'points', my_var.FV_NETCDF["float64"], compress, group=pixc)
        fill_vector_param(self.computeTime_TAI(self.illumination_time), 'illumination_time_tai', self.nb_water_pix, data, group=pixc)  # TODO: to improve
        
        data.add_variable('num_med_looks', np.int32, 'points', my_var.FV_NETCDF["int32"], compress, group=pixc)
        fill_vector_param(np.full(self.nb_water_pix, 63.), 'num_med_looks', self.nb_water_pix, data, group=pixc)
        data.add_variable('sig0', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'sig0', self.nb_water_pix, data, group=pixc)
        data.add_variable('phase_unwrapping_region', np.int32, 'points', my_var.FV_NETCDF["int32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'phase_unwrapping_region', self.nb_water_pix, data, group=pixc)
        
        data.add_variable('instrument_range_cor', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'instrument_range_cor', self.nb_water_pix, data, group=pixc)
        data.add_variable('instrument_phase_cor', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'instrument_phase_cor', self.nb_water_pix, data, group=pixc)
        data.add_variable('instrument_baseline_cor', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'instrument_baseline_cor', self.nb_water_pix, data, group=pixc)
        data.add_variable('instrument_attitude_cor', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'instrument_attitude_cor', self.nb_water_pix, data, group=pixc)

        data.add_variable('model_dry_tropo_cor', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'model_dry_tropo_cor', self.nb_water_pix, data, group=pixc)
        data.add_variable('model_wet_tropo_cor', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'model_wet_tropo_cor', self.nb_water_pix, data, group=pixc)
        data.add_variable('iono_cor_gim_ka', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'iono_cor_gim_ka', self.nb_water_pix, data, group=pixc)     
        data.add_variable('xover_height_cor', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'xover_height_cor', self.nb_water_pix, data, group=pixc)
        data.add_variable('load_tide_sol1', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'load_tide_sol1', self.nb_water_pix, data, group=pixc)
        data.add_variable('load_tide_sol2', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'load_tide_sol2', self.nb_water_pix, data, group=pixc)
        data.add_variable('pole_tide', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'pole_tide', self.nb_water_pix, data, group=pixc)
        data.add_variable('solid_earth_tide', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'solid_earth_tide', self.nb_water_pix, data, group=pixc)
        data.add_variable('geoid', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'geoid', self.nb_water_pix, data, group=pixc)
        data.add_variable('surface_type_flag', np.int8, 'points', my_var.FV_NETCDF["int8"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'surface_type_flag', self.nb_water_pix, data, group=pixc)
        data.add_variable('pixc_qual', np.int8, 'points', my_var.FV_NETCDF["int8"], compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'pixc_qual', self.nb_water_pix, data, group=pixc) 
        
        # ===============
        # == Group TVP ==
        # ===============
        sensor = data.add_group("tvp")
        
        # Group attributes
        data.add_global_attribute('description', 'Time varying parameters group including spacecraft attitude, position, velocity, and antenna position information', group=sensor)

        # Group dimension
        data.add_dimension('num_tvps', self.nb_nadir_pix, group=sensor)

        # Group variables
        data.add_variable('time', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(self.computeTime_UTC(self.nadir_time), 'time', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('time_tai', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(self.computeTime_TAI(self.nadir_time), 'time_tai', self.nb_nadir_pix, data, group=sensor)
        
        data.add_variable('latitude', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        data.add_variable_attribute('latitude', 'units', 'degrees_north', group=sensor)
        fill_vector_param(self.nadir_latitude, 'latitude', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('longitude', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        data.add_variable_attribute('longitude', 'units', 'degrees_east', group=sensor)
        fill_vector_param(self.nadir_longitude, 'longitude', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('altitude', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(self.nadir_altitude, 'altitude', self.nb_nadir_pix, data, group=sensor)
        
        data.add_variable('roll', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        data.add_variable_attribute('roll', 'units', 'degrees', group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'roll', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('pitch', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        data.add_variable_attribute('pitch', 'units', 'degrees', group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'pitch', self.nb_nadir_pix, data, group=sensor)  
        data.add_variable('yaw', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        data.add_variable_attribute('yaw', 'units', 'degrees', group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'yaw', self.nb_nadir_pix, data, group=sensor)   
        data.add_variable('velocity_heading', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(self.nadir_heading, 'velocity_heading', self.nb_nadir_pix, data, group=sensor)
        data.add_variable_attribute('velocity_heading', 'units', 'degrees', group=sensor)
        
        data.add_variable('x', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(self.nadir_x, 'x', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('y', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(self.nadir_y, 'y', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('z', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(self.nadir_z, 'z', self.nb_nadir_pix, data, group=sensor)
    
        data.add_variable('vx', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(self.nadir_vx, 'vx', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('vy', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(self.nadir_vy, 'vy', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('vz', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(self.nadir_vz, 'vz', self.nb_nadir_pix, data, group=sensor)
        
        data.add_variable('plus_y_antenna_x', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'plus_y_antenna_x', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('plus_y_antenna_y', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'plus_y_antenna_y', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('plus_y_antenna_z', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'plus_y_antenna_z', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('minus_y_antenna_x', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'minus_y_antenna_x', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('minus_y_antenna_y', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'minus_y_antenna_y', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('minus_y_antenna_z', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'minus_y_antenna_z', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('record_counter', np.int32, 'num_tvps', my_var.FV_NETCDF["int32"], compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'record_counter', self.nb_nadir_pix, data, group=sensor) 
        data.add_variable('sc_event_flag', np.int8, 'num_tvps', my_var.FV_NETCDF["int8"], compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'sc_event_flag', self.nb_nadir_pix, data, group=sensor) 
        data.add_variable('tvp_qual', np.int8, 'num_tvps', my_var.FV_NETCDF["int8"], compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'tvp_qual', self.nb_nadir_pix, data, group=sensor) 
                
        # =================
        # == Group Noise ==
        # =================
        noise = data.add_group("noise")
        
        # Group attributes
        data.add_global_attribute('description', 'Measured noise power for each recieve echo of the plus_y and minus_y SLC channels', group=noise)
 
        # Group dimension        
        data.add_dimension('num_lines', self.nb_nadir_pix, group=noise)

        # Group variables
        data.add_variable('noise_plus_y', np.float32, 'num_lines', my_var.FV_NETCDF["float32"], compress, group=noise)
        fill_vector_param(np.full(self.nb_nadir_pix, -116.845780895788), 'noise_plus_y', self.nb_nadir_pix, data, group=noise)
        data.add_variable('noise_minus_y', np.float32, 'num_lines', my_var.FV_NETCDF["float32"], compress, group=noise)
        fill_vector_param(np.full(self.nb_nadir_pix, -116.845780895788), 'noise_minus_y', self.nb_nadir_pix, data, group=noise)

        # Close NetCDF file
        data.close()
    
    #----------------------------------
 
    def write_annotation_file(self, IN_output_file, IN_pixc_file):
        """
        write the river-annotation.rdf file so that lake processor can run
        
        :param IN_output_file: output full path
        :type IN_output_file: string
        :param IN_pixc_file: PIXC full path
        :type IN_pixc_file: string
        """
        my_api.printInfo("[proc_real_pixc] == write_annotation_file : %s ==" % IN_output_file)
        
        f = open(IN_output_file, 'w')
        f.write("l2pixc file = %s\n" % IN_pixc_file)
        
        f.close()
    
    #----------------------------------

    def write_pixc_asShp(self, IN_output_file):
        """
        Write some of the pixel cloud attributes in a shapefile

        :param IN_output_file: output full path
        :type IN_output_file: string
        """
        my_api.printInfo("[proc_real_pixc] == write_pixc_asShp : %s ==" % IN_output_file) 
        
        # 1 - Initialisation du fichier de sortie
        # 1.1 - Driver
        shpDriver = ogr.GetDriverByName(str("ESRI Shapefile"))
        # 1.2 - Creation du fichier
        if os.path.exists(IN_output_file):
            shpDriver.DeleteDataSource(IN_output_file)
        outDataSource = shpDriver.CreateDataSource(IN_output_file)
        # 1.3 - Creation de la couche
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)  # WGS84
        outLayer = outDataSource.CreateLayer(str(os.path.basename(IN_output_file).split('.')[0]+"_pixc"), srs, geom_type=ogr.wkbPoint)
        # 1.4 - Creation des attributs
        outLayer.CreateField(ogr.FieldDefn(str('az_index'), ogr.OFTInteger))  # Azimuth index
        outLayer.CreateField(ogr.FieldDefn(str('r_index'), ogr.OFTInteger))  # Range index
        outLayer.CreateField(ogr.FieldDefn(str('classif'), ogr.OFTInteger))  # Classification
        tmpField = ogr.FieldDefn(str('pix_area'), ogr.OFTReal)  # Pixel area
        tmpField.SetWidth(15)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('lat'), ogr.OFTReal)  # Latitude
        tmpField.SetWidth(15)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('long'), ogr.OFTReal)  # Longitude
        tmpField.SetWidth(15)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('wse'), ogr.OFTReal)  # Hauteur
        tmpField.SetWidth(15)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('cr_track'), ogr.OFTReal)  # Distance dans la fauchee
        tmpField.SetWidth(15)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        # 1.5 - On recupere la definition de la couche
        outLayerDefn = outLayer.GetLayerDefn()
        
        # 2 - On traite point par point
        for az_ind, range_index, classif, pixel_area, lat, lng, height, crosstrack in zip(self.azimuth_index, self.range_index, self.classification, self.pixel_area, self.latitude, self.longitude, self.height, self.crosstrack):
            # 2.1 - On cree l'objet dans le format de la couche de sortie
            outFeature = ogr.Feature(outLayerDefn)
            # 2.2 - On lui assigne le point
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(lng, lat)
            outFeature.SetGeometry(point)
            # 2.3 - On lui assigne les attributs
            outFeature.SetField(str('az_index'), float(az_ind))
            outFeature.SetField(str('r_index'), float(range_index))
            outFeature.SetField(str('classif'), float(classif))
            outFeature.SetField(str('pix_area'), float(pixel_area))
            outFeature.SetField(str('lat'), float(lat))
            outFeature.SetField(str('long'), float(lng))
            outFeature.SetField(str('wse'), float(height))
            outFeature.SetField(str('cr_track'), float(crosstrack))
            # 2.4 - On ajoute l'objet dans la couche de sortie
            outLayer.CreateFeature(outFeature)
            
        # 3 - Destroy the data sources to free resources
        outDataSource.Destroy()
        
    def write_tvp_asShp(self, IN_output_file):
        """
        Write some of the TVP attributes in a shapefile

        :param IN_output_file: output full path
        :type IN_output_file: string
        """
        my_api.printInfo("[proc_real_pixc] == write_tvp_asShp : %s ==" % IN_output_file) 
    
        # 1 - Initialisation du fichier de sortie
        # 1.1 - Driver
        shpDriver = ogr.GetDriverByName(str("ESRI Shapefile"))
        # 1.2 - Creation du fichier
        if os.path.exists(IN_output_file):
            shpDriver.DeleteDataSource(IN_output_file)
        outDataSource = shpDriver.CreateDataSource(IN_output_file)
        # 1.3 - Creation de la couche
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)  # WGS84
        outLayer = outDataSource.CreateLayer(str(os.path.basename(IN_output_file).split('.')[0]+"_tvp"), srs, geom_type=ogr.wkbPoint)
        # 1.4 - Creation des attributs
        tmpField = ogr.FieldDefn(str('time'), ogr.OFTReal)  # Time
        tmpField.SetWidth(20)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('lat'), ogr.OFTReal)  # Latitude
        tmpField.SetWidth(15)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('long'), ogr.OFTReal)   # Longitude
        tmpField.SetWidth(15)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('altitude'), ogr.OFTReal)  # Altitude
        tmpField.SetWidth(15)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('heading'), ogr.OFTReal)  # Heading
        tmpField.SetWidth(15)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        # 1.5 - On recupere la definition de la couche
        outLayerDefn = outLayer.GetLayerDefn()
        
        # 2 - On traite point par point
        for lng, lat, t, heading, alt in zip(self.nadir_longitude, self.nadir_latitude, self.nadir_time, self.nadir_heading, self.nadir_altitude):
            # 2.1 - On cree l'objet dans le format de la couche de sortie

            outFeature = ogr.Feature(outLayerDefn)
            # 2.2 - On lui assigne le point
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(lng, lat)
            outFeature.SetGeometry(point)
            # 2.3 - On lui assigne les attributs
            outFeature.SetField(str('time'), float(t)) 
            outFeature.SetField(str('lat'), float(lat)) 
            outFeature.SetField(str('long'), float(lng)) 
            outFeature.SetField(str('altitude'), float(alt)) 
            outFeature.SetField(str('heading'), float(heading))
            # 2.4 - On ajoute l'objet dans la couche de sortie
            outLayer.CreateFeature(outFeature)
            
        # 3 - Destroy the data sources to free resources
        outDataSource.Destroy()
    
    #----------------------------------
        
    def computeDate(self, IN_sec_from_start):
        """
        Compute date
        
        :param IN_sec_from_start: number of seconds from mission start time
        :type IN_sec_from_start: int
        
        :return: date in UTC
        :rtype: string YYYYMMDDThhmmss
        """
        
        # Computation
        tmp_time_split = self.mission_start_time.split("-")
        date_in_sec = datetime(int(tmp_time_split[0]), int(tmp_time_split[1]), int(tmp_time_split[2])) + timedelta(seconds=IN_sec_from_start)
        
        # Format
        return datetime.strftime(date_in_sec, '%Y%m%dT%H%M%S')
        
    def computeTime_UTC(self, IN_sec_from_start):
        """
        Compute time in seconds from 01/01/2000 00:00:00
        
        :param IN_sec_from_start: number of seconds from mission start time
        :type IN_sec_from_start: int
        
        :return: time in seconds in UTC time scale
        :rtype: float
        """
        
        # Convert mission start time to datetime
        tmp_time_split = self.mission_start_time.split("-")
        mission_start_time = datetime(int(tmp_time_split[0]), int(tmp_time_split[1]), int(tmp_time_split[2]))
        
        # Convert reference to datetime
        ref_time = datetime(2000,1,1)
        
        # Compute difference
        diff = mission_start_time - ref_time
        
        # Return number of seconds of difference
        return IN_sec_from_start + diff.total_seconds()
        
    def computeTime_TAI(self, IN_sec_from_start):
        """
        Compute time in seconds from 01/01/2000 00:00:32
        
        :param IN_sec_from_start: number of seconds from mission start time
        :type IN_sec_from_start: int
        
        :return: time in seconds in TAI time scale
        :rtype: float
        """
        
        # Convert mission start time to datetime
        tmp_time_split = self.mission_start_time.split("-")
        mission_start_time = datetime(int(tmp_time_split[0]), int(tmp_time_split[1]), int(tmp_time_split[2]))
        
        # Convert reference to datetime
        ref_time = datetime(2000,1,1,0,0,32)
        
        # Compute difference
        diff = mission_start_time - ref_time
        
        # Return number of seconds of difference
        return IN_sec_from_start + diff.total_seconds()

