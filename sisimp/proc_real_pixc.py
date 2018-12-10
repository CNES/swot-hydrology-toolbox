# -*- coding: utf8 -*-
"""
.. module proc_real_pixc.py
    :synopsis: Deal with official pixel cloud (L2_HR_PIXC) files
    Created on 08/24/2017
    2018/11/30 (D. Desroches, V. Poughon - CNES): change variables names wrt to new PixC naming convention

.. module author: Claire POTTIER - CNES DSO/SI/TR

Copyright (c) 2017 CNES. All rights reserved.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import os
from osgeo import ogr, osr
from datetime import datetime, timedelta

import lib.my_api as my_api
import lib.my_netcdf_file as my_nc


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

    def __init__(self, IN_azimuth_index, IN_range_index, IN_classification, IN_pixel_area, IN_latitude, IN_longitude, IN_height, IN_crosstrack, \
                 IN_nadir_time, IN_nadir_latitude, IN_nadir_longitude, IN_nadir_altitude, IN_nadir_heading, IN_nadir_x, IN_nadir_y, IN_nadir_z, IN_nadir_vx, IN_nadir_vy, IN_nadir_vz, IN_nadir_near_range, \
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
            
        :param IN_nadir_time: time tags for each nadir record ; provided as UTC seconds since begin of current cycle
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
        """
        my_api.printInfo("[l2_hr_pixc] == INIT ==") 
        
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

        #self.illumination_time = [self.nadir_time[self.sensor_s[i]] for i in range((self.illumination_time).size)]
            
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

    def write_pixc_file(self, IN_output_file, noval, compress=False):
        """
        Write the main file of real pixel cloud product (L2_HR_PIXC product, main file)

        :param IN_output_file: output full path
        :type IN_output_file: string
        :param noval: No data value
        :type noval: float
        :param compress: parameter the define to compress or not the file
        :type compress: boolean
        """
        my_api.printInfo("[l2_hr_pixc] == write_pixc_file : %s ==" % IN_output_file) 
    
        data = my_nc.myNcWriter(IN_output_file)
        
        if noval is None:
            noval = -999000000.

        # Group pixel_cloud

        pixc = data.add_group("pixel_cloud")
  
        data.add_dimension('record', self.nb_water_pix, group=pixc)
        data.add_variable('azimuth_index', np.int64, 'record', np.int(noval), compress, group=pixc)
        fill_vector_param(self.azimuth_index, 'azimuth_index', self.nb_water_pix, data, group=pixc)
        data.add_variable('range_index', np.int64, 'record', np.int(noval), compress, group=pixc)
        fill_vector_param(self.range_index, 'range_index', self.nb_water_pix, data, group=pixc)
    
        data.add_variable('classification', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(self.classification, 'classification', self.nb_water_pix, data, group=pixc)
        data.add_variable('continuous_classification', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(self.classification, 'continuous_classification', self.nb_water_pix, data, group=pixc)
        data.add_variable('reference_layover_flag', np.byte, 'record', None, compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'reference_layover_flag', self.nb_water_pix, data, group=pixc)
        data.add_variable('reference_layover_height_error', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'reference_layover_height_error', self.nb_water_pix, data, group=pixc)
    
        data.add_variable('ice_flag', np.byte, 'record', None, compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'ice_flag', self.nb_water_pix, data, group=pixc)
        data.add_variable('dark_water_flag', np.byte, 'record', None, compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dark_water_flag', self.nb_water_pix, data, group=pixc)
    
        data.add_variable('pixel_area', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(self.pixel_area, 'pixel_area', self.nb_water_pix, data, group=pixc)
    
        data.add_variable('latitude', np.float64, 'record', np.float(noval), compress, group=pixc)
        data.add_variable_attribute('latitude', 'units', 'degrees_north', group=pixc)
        fill_vector_param(self.latitude, 'latitude', self.nb_water_pix, data, group=pixc)
        data.add_variable('longitude', np.float64, 'record', np.float(noval), compress, group=pixc)
        data.add_variable_attribute('longitude', 'units', 'degrees_east', group=pixc)
        fill_vector_param(self.longitude, 'longitude', self.nb_water_pix, data, group=pixc)
        data.add_variable('height', np.float64, 'record', np.float(noval), compress, group=pixc)
        data.add_variable_attribute('height', 'units', 'm', group=pixc)
        fill_vector_param(self.height, 'height', self.nb_water_pix, data, group=pixc)
        data.add_variable('cross_track', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(self.crosstrack, 'cross_track', self.nb_water_pix, data, group=pixc)
        data.add_variable('look_unit_x', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'look_unit_x', self.nb_water_pix, data, group=pixc)
        data.add_variable('look_unit_y', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'look_unit_y', self.nb_water_pix, data, group=pixc)
        data.add_variable('look_unit_z', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'look_unit_z', self.nb_water_pix, data, group=pixc)
        data.add_variable('dlook_dphase_x', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dlook_dphase_x', self.nb_water_pix, data, group=pixc)
        data.add_variable('dlook_dphase_y', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dlook_dphase_y', self.nb_water_pix, data, group=pixc)
        data.add_variable('dlook_dphase_z', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dlook_dphase_z', self.nb_water_pix, data, group=pixc)
        data.add_variable('dheight_dphase', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dheight_dphase', self.nb_water_pix, data, group=pixc)
        data.add_variable('dheight_droll', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dheight_droll', self.nb_water_pix, data, group=pixc)
        data.add_variable('dheight_dbaseline', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dheight_dbaseline', self.nb_water_pix, data, group=pixc)
        data.add_variable('dheight_drange', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dheight_drange', self.nb_water_pix, data, group=pixc)
        data.add_variable('dphase', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dphase', self.nb_water_pix, data, group=pixc)
        data.add_variable('delta_x', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'delta_x', self.nb_water_pix, data, group=pixc)
        data.add_variable('delta_y', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'delta_y', self.nb_water_pix, data, group=pixc)
        data.add_variable('delta_z', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'delta_z', self.nb_water_pix, data, group=pixc)
        
        data.add_variable('ifgram_real', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'ifgram_real', self.nb_water_pix, data, group=pixc)
        data.add_variable('ifgram_imag', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'ifgram_imag', self.nb_water_pix, data, group=pixc)
        data.add_variable('power_left', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'power_left', self.nb_water_pix, data, group=pixc)
        data.add_variable('power_right', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'power_right', self.nb_water_pix, data, group=pixc)
        data.add_variable('coherent_power', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'coherent_power', self.nb_water_pix, data, group=pixc)
        data.add_variable('num_rare_looks', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.full(self.nb_water_pix, 7.), 'num_rare_looks', self.nb_water_pix, data, group=pixc)
    
        data.add_variable('instrument_range_correction', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'instrument_range_correction', self.nb_water_pix, data, group=pixc)
        data.add_variable('instrument_phase_correction', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'instrument_phase_correction', self.nb_water_pix, data, group=pixc)
        data.add_variable('instrument_baseline_correction', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'instrument_baseline_correction', self.nb_water_pix, data, group=pixc)
        data.add_variable('instrument_attitude_correction', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'instrument_attitude_correction', self.nb_water_pix, data, group=pixc)
        data.add_variable('xover_roll_correction', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'xover_roll_correction', self.nb_water_pix, data, group=pixc)
        data.add_variable('ionosphere_range_correction', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'ionosphere_range_correction', self.nb_water_pix, data, group=pixc)
        data.add_variable('wet_tropo_range_correction', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'wet_tropo_range_correction', self.nb_water_pix, data, group=pixc)
        data.add_variable('dry_tropo_range_correction', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dry_tropo_range_correction', self.nb_water_pix, data, group=pixc)
        data.add_variable('solid_earth_tide_height_correction', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'solid_earth_tide_height_correction', self.nb_water_pix, data, group=pixc)
        data.add_variable('phase_screen', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'phase_screen', self.nb_water_pix, data, group=pixc)
 
        data.add_variable('sensor_s', np.int64, 'record', np.int(noval), compress, group=pixc)
        fill_vector_param(self.sensor_s, 'sensor_s', self.nb_water_pix, data, group=pixc)
        data.add_variable('illumination_time', np.float64, 'record', np.float(noval), compress, group=pixc)
        fill_vector_param(self.illumination_time, 'illumination_time', self.nb_water_pix, data, group=pixc)

        # Group TVP
        
        sensor = data.add_group("tvp")
        data.add_dimension('record', self.nb_nadir_pix, group=sensor)

        data.add_variable('time', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(self.nadir_time, 'time', self.nb_nadir_pix, data, group=sensor)
    
        data.add_variable('latitude', np.float64, 'record', np.float(noval), compress, group=sensor)
        data.add_variable_attribute('latitude', 'units', 'degrees_north', group=sensor)
        fill_vector_param(self.nadir_latitude, 'latitude', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('longitude', np.float64, 'record', np.float(noval), compress, group=sensor)
        data.add_variable_attribute('longitude', 'units', 'degrees_east', group=sensor)
        fill_vector_param(self.nadir_longitude, 'longitude', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('height', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(self.nadir_altitude, 'height', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('heading', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(self.nadir_heading, 'heading', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('x', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(self.nadir_x, 'x', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('y', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(self.nadir_y, 'y', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('z', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(self.nadir_z, 'z', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('near_range', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(self.nadir_near_range, 'near_range', self.nb_nadir_pix, data, group=sensor)
    
        data.add_variable('baseline_left_x', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'baseline_left_x', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('baseline_left_y', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'baseline_left_y', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('baseline_left_z', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'baseline_left_z', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('baseline_right_x', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'baseline_right_x', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('baseline_right_y', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'baseline_right_y', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('baseline_right_z', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'baseline_right_z', self.nb_nadir_pix, data, group=sensor)
    
        data.add_variable('vx', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(self.nadir_vx, 'vx', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('vy', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(self.nadir_vy, 'vy', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('vz', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(self.nadir_vz, 'vz', self.nb_nadir_pix, data, group=sensor)
        
        data.add_variable('ref_leverarm_x', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'ref_leverarm_x', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('ref_leverarm_y', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'ref_leverarm_y', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('ref_leverarm_z', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'ref_leverarm_z', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('sec_leverarm_x', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'sec_leverarm_x', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('sec_leverarm_y', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'sec_leverarm_y', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('sec_leverarm_z', np.float64, 'record', np.float(noval), compress, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'sec_leverarm_z', self.nb_nadir_pix, data, group=sensor)
        
        # Global attributes
                   
        data.add_global_attribute('description', 'L2_HR_PIXC product obtained by CNES/LEGOS Large Scale Simulator', group=pixc)
        data.add_global_attribute('mission_start_time', self.mission_start_time, group=pixc)
        data.add_global_attribute('repeat_cycle_period', self.cycle_duration, group=pixc)
        data.add_global_attribute('start_time', np.min(self.nadir_time), group=pixc)
        data.add_global_attribute('stop_time', np.max(self.nadir_time), group=pixc)
        data.add_global_attribute('cycle_number', self.cycle_num, group=pixc)
        data.add_global_attribute('pass_number', np.int(self.pass_num), group=pixc)
        data.add_global_attribute('tile_ref', self.tile_ref, group=pixc)
        data.add_global_attribute('nr_pixels', self.nb_pix_range, group=pixc)
        data.add_global_attribute('nr_lines', self.nb_pix_azimuth, group=pixc)

        data.add_global_attribute('range_spacing', self.range_spacing, group=pixc)
        data.add_global_attribute('azimuth_spacing', self.azimuth_spacing, group=pixc)
        data.add_global_attribute('near_range', self.near_range, group=pixc)
        
        data.add_global_attribute("inner_first_lat", self.latitude[np.argmin(self.latitude)], group=pixc)
        data.add_global_attribute("outer_first_lat", self.latitude[np.argmax(self.latitude)], group=pixc)
        data.add_global_attribute("inner_last_lat", self.latitude[np.argmin(self.longitude)], group=pixc)
        data.add_global_attribute("outer_last_lat", self.latitude[np.argmax(self.longitude)], group=pixc)
    
        data.add_global_attribute("inner_first_lon", self.longitude[np.argmin(self.latitude)], group=pixc)
        data.add_global_attribute("outer_first_lon", self.longitude[np.argmax(self.latitude)], group=pixc)
        data.add_global_attribute("inner_last_lon", self.longitude[np.argmin(self.longitude)], group=pixc)
        data.add_global_attribute("outer_last_lon", self.longitude[np.argmax(self.longitude)], group=pixc)

        # Needed by some processing code, default values here
        data.add_global_attribute('noise_power', 0.0, group=pixc)
        data.add_global_attribute('nb_az_only_looks', 4, group=pixc)
        data.add_global_attribute("noise_power_left", -116.845780895788, group=pixc)
        data.add_global_attribute("noise_power_right", -116.845780895788, group=pixc)
        data.add_global_attribute("wavelength", 0.008385803020979, group=pixc)
        
        data.close()
 
    def write_annotation_file(self, IN_output_file, IN_pixc_file):
        """
        write the river-annotation.rdf file so that lake processor can run
        
        :param IN_output_file: output full path
        :type IN_output_file: string
        :param IN_pixc_file: PIXC full path
        :type IN_pixc_file: string
        """
        my_api.printInfo("[l2_hr_pixc] == write_annotation_file : %s ==" % IN_output_file)
        
        f = open(IN_output_file, 'w')
        f.write("l2pixc file = %s\n" % IN_pixc_file)
        f.write("true GDEM file = /work/ALT/swot/swotdev/desrochesd/swot-sds/run/camargue/output_simu/gdem-dem-camargue/cycle_0001_pass_0004_LeftSwath_nlcd50_darkWater_25m_CBE/gdem_truth.LeftSwath.nc")
        
        f.close()

    def write_pixc_asShp(self, IN_output_file):
        """
        Write some of the pixel cloud attributes in a shapefile

        :param IN_output_file: output full path
        :type IN_output_file: string
        """
        my_api.printInfo("[l2_hr_pixc] == write_pixc_asShp : %s ==" % IN_output_file) 
        
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
        outLayer.CreateField(ogr.FieldDefn(str('AZ_INDEX'), ogr.OFTInteger))  # Azimuth index
        outLayer.CreateField(ogr.FieldDefn(str('R_INDEX'), ogr.OFTInteger))  # Range index
        outLayer.CreateField(ogr.FieldDefn(str('CLASSIF'), ogr.OFTInteger))  # Classification
        tmpField = ogr.FieldDefn(str('PIX_AREA'), ogr.OFTReal)  # Pixel area
        tmpField.SetWidth(10)
        tmpField.SetPrecision(5)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('LAT'), ogr.OFTReal)  # Latitude
        tmpField.SetWidth(10)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('LONG'), ogr.OFTReal)  # Longitude
        tmpField.SetWidth(10)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('HEIGHT'), ogr.OFTReal)  # Hauteur
        tmpField.SetWidth(10)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('CR_TRACK'), ogr.OFTReal) # Distance dans la fauchee
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
            outFeature.SetField(str('AZ_INDEX'), float(az_ind))
            outFeature.SetField(str('R_INDEX'), float(range_index))
            outFeature.SetField(str('CLASSIF'), float(classif))
            outFeature.SetField(str('PIX_AREA'), float(pixel_area))
            outFeature.SetField(str('LAT'), float(lat))
            outFeature.SetField(str('LONG'), float(lng))
            outFeature.SetField(str('HEIGHT'), float(height))
            outFeature.SetField(str('CR_TRACK'), float(crosstrack))
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
        my_api.printInfo("[l2_hr_pixc] == write_tvp_asShp : %s ==" % IN_output_file) 
    
        # 1 - Initialisation du fichier de sortie
        # 1.1 - Driver
        shpDriver = ogr.GetDriverByName(str("ESRI Shapefile"))
        # 1.2 - Creation du fichier
        if os.path.exists(IN_output_file):
            shpDriver.DeleteDataSource(IN_output_file)
        outDataSource = shpDriver.CreateDataSource(IN_output_file)
        # 1.3 - Creation de la couche
        #print '> Creating layer pixel_cloud'
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)  # WGS84
        outLayer = outDataSource.CreateLayer(str(os.path.basename(IN_output_file).split('.')[0]+"_tvp"), srs, geom_type=ogr.wkbPoint)
        # 1.4 - Creation des attributs
        tmpField = ogr.FieldDefn(str('TIME'), ogr.OFTReal)  # Time
        tmpField.SetWidth(10)
        tmpField.SetPrecision(2)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('LAT'), ogr.OFTReal) # Latitude
        tmpField.SetWidth(10)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('LONG'), ogr.OFTReal)   # Longitude
        tmpField.SetWidth(10)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('ALTITUDE'), ogr.OFTReal)  # Altitude
        tmpField.SetWidth(10)
        tmpField.SetPrecision(3)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('HEADING'), ogr.OFTReal)  # Heading
        tmpField.SetWidth(10)
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
            #~ outFeature.SetField(str('TIME'), float(t)) 
            outFeature.SetField(str('LAT'), float(lat)) 
            outFeature.SetField(str('LONG'), float(lng)) 
            outFeature.SetField(str('ALTITUDE'), float(alt)) 
            outFeature.SetField(str('HEADING'), float(heading))
            # 2.4 - On ajoute l'objet dans la couche de sortie
            outLayer.CreateFeature(outFeature)
            
        # 3 - Destroy the data sources to free resources
        outDataSource.Destroy()


##########################
        

if __name__ == '__main__':
    
    noval = -999900000

    record = 10
    azimuth_index = np.ones(record) + 50.
    range_index = np.ones(record)+50.
    classification = np.ones(record)
    pixel_area = np.ones(record)+10.
    ambiguity_altitude = np.ones(record)+3.
    latitude = np.arange(record)+10.
    longitude = np.arange(record)+50.
    height = np.ones(record)+50.
    crosstrack = np.ones(record)+20.
    range_tab = np.ones(record)+50.
    
    record = 200
    nadir_time = np.ones(record)+10000.
    nadir_latitude = np.arange(record)+10.
    nadir_longitude = np.arange(record)+50.
    nadir_altitude = np.ones(record)*890.
    nadir_heading = np.ones(record)*77.6
    nadir_x = np.ones(record)*20.
    nadir_y = np.ones(record)*30.
    nadir_z = np.ones(record)*40.
    nadir_near_range = np.ones(record)*10.0
    
    cycle_num = 1
    pass_num = 10
    tile_ref = "45N-L"
    nb_pix_range = 200
    nb_pix_azimuth = 20
    
    my_pixc = l2_hr_pixc(azimuth_index, range_index, classification, pixel_area, ambiguity_altitude, latitude, longitude, height, crosstrack, \
                         nadir_time, nadir_latitude, nadir_longitude, nadir_altitude, nadir_heading, nadir_x, nadir_y, nadir_z, nadir_near_range, \
                         cycle_num, pass_num, tile_ref, nb_pix_range, nb_pix_azimuth)
    
    my_pixc.write_main("D:\\Utilisateurs\\pottierc\\Documents\\workspace_eclipse\\swot_py\\res\\sisimp_real_pixc\\outputs\\", noval, True)
    my_pixc.write_sensor("D:\\Utilisateurs\\pottierc\\Documents\\workspace_eclipse\\swot_py\\res\\sisimp_real_pixc\\outputs\\", noval, True)
    my_pixc.write_main_asShp("D:\\Utilisateurs\\pottierc\\Documents\\workspace_eclipse\\swot_py\\res\\sisimp_real_pixc\\outputs\\")
    my_pixc.write_sensor_asShp("D:\\Utilisateurs\\pottierc\\Documents\\workspace_eclipse\\swot_py\\res\\sisimp_real_pixc\\outputs\\")
