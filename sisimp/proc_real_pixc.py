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

    def __init__(self, IN_azimuth_index, IN_range_index, IN_classification, IN_pixel_area, IN_latitude, IN_longitude, IN_height, IN_phase_noise_std,
                 IN_dh_dphi, IN_dlon_dphi, IN_dlat_dphi, IN_crosstrack,
                 IN_nadir_time, IN_nadir_latitude, IN_nadir_longitude, IN_nadir_altitude, IN_nadir_heading, IN_nadir_x, IN_nadir_y, IN_nadir_z, IN_nadir_vx, IN_nadir_vy, IN_nadir_vz, IN_nadir_near_range,
                 IN_mission_start_time, IN_cycle_duration, IN_cycle_num, IN_pass_num, IN_tile_ref, IN_nb_pix_range, IN_nb_pix_azimuth, IN_azimuth_spacing, IN_range_spacing, IN_near_range, IN_tile_coords, IN_interferogram):
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
        :param IN_phase_noise_std: phase noise standart deviation
        :type IN_phase_noise_std: 1D-array of float        
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
        :param IN_tile_coords: tile coordinates (inner_first, inner_last, outer_first, outer_last), inner_first=(lon, lat)
        :type IN_tile_coords: tuple of tuple of float
            
        + nb_water_pix(int) : number of water pixels, i.e. pixels in azimuth_index, ..., crosstrack vectors
        + nb_nadir_pix(int) : number of pixels on the nadir track, i.e. pixels in time, ..., near_range vectors
        + pattern(str): filename pattern
        """
        my_api.printInfo("[proc_real_pixc] == INIT ==")

        self.azimuth_index = IN_azimuth_index
        self.range_index = IN_range_index
        self.classification = IN_classification
        self.pixel_area = IN_pixel_area
        self.latitude = IN_latitude
        self.longitude = IN_longitude
        self.height = IN_height
        self.phase_noise_std = IN_phase_noise_std
        self.dh_dphi = IN_dh_dphi
        self.dlon_dphi = IN_dlon_dphi
        self.dlat_dphi = IN_dlat_dphi
        
        self.crosstrack = IN_crosstrack
        self.nb_water_pix = IN_azimuth_index.size
        self.water_frac = np.ones(self.nb_water_pix)

        # Modification to have sensor_s (sensor azimuth position for each pixel) to be compatible with HR simulator. It is a duplication of azimuth_index in the large scale simulator
        self.sensor_s = IN_azimuth_index
        self.nadir_time = IN_nadir_time

        if np.max(IN_azimuth_index) >= IN_nadir_time.size:
            exc = '[proc_realPixC] ERROR = Azimuth index max value %d over nb_nadir_pix %d' %(np.max(IN_azimuth_index), IN_nadir_time.size)
            exit(exc)

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
        self.interferogram = IN_interferogram
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

        (inner_first, inner_last, outer_first, outer_last) = IN_tile_coords
        self.inner_first = inner_first
        self.inner_last = inner_last
        self.outer_first = outer_first
        self.outer_last = outer_last


    #----------------------------------

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
        data.add_global_attribute("inner_first_longitude", self.inner_first[0])
        data.add_global_attribute("inner_first_latitude", self.inner_first[1])
        data.add_global_attribute("inner_last_longitude", self.inner_last[0])
        data.add_global_attribute("inner_last_latitude", self.inner_last[1])
        data.add_global_attribute("outer_first_longitude", self.outer_first[0])
        data.add_global_attribute("outer_first_latitude", self.outer_first[1])
        data.add_global_attribute("outer_last_longitude", self.outer_last[0])
        data.add_global_attribute("outer_last_latitude", self.outer_last[1])
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
        dic={'long_name':'rare interferogram azimuth index','units':'1','valid_min':'0','valid_max':'1','coordinates':'longitude latitude','comment':'Rare interferogram azimuth index'}
        data.add_variable_attributes('azimuth_index', dic, group=pixc)
        fill_vector_param(self.azimuth_index, 'azimuth_index', self.nb_water_pix, data, group=pixc)

        data.add_variable('range_index', np.int32, 'points', my_var.FV_NETCDF["int32"], compress, group=pixc)
        dic={'long_name':'rare interferogram range index','units':'1','valid_min':'0','valid_max':'999999','coordinates':'longitude latitude','comment':'Rare interferogram range index'}
        data.add_variable_attributes('range_index', dic, group=pixc)
        fill_vector_param(self.range_index, 'range_index', self.nb_water_pix, data, group=pixc)

        #------

        data.add_variable('interferogram', np.float32, ('points', 'depth'), my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'rare interferogram','units':'1','valid_min':'-999999','valid_max':'999999','coordinates':'longitude latitude','comment':'Complex unflattened rare interferogram'}
        data.add_variable_attributes('interferogram', dic, group=pixc)
        fill_vector_param(self.interferogram, 'interferogram', self.nb_water_pix, data, group=pixc)#np.ones([self.nb_water_pix, 2])/np.sqrt(2)

        data.add_variable('power_plus_y', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'power for plus_y channel','units':'1','valid_min':'0','valid_max':'999999','coordinates':'longitude latitude','comment':'Power for plus_y channel (arbitrary units that give sigma0 when nois substracted and normalized by the X factor)'}
        data.add_variable_attributes('power_plus_y', dic, group=pixc)
        fill_vector_param(np.ones(self.nb_water_pix), 'power_plus_y', self.nb_water_pix, data, group=pixc)
        
        data.add_variable('power_minus_y', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'power for minus_y channel','units':'1','valid_min':'0','valid_max':'999999','coordinates':'longitude latitude','comment':'Power for the minus_y channel (arbitrary units that give sigma0 when noise substracted and normalized by the X factor)'} 
        data.add_variable_attributes('power_minus_y', dic, group=pixc)
        fill_vector_param(np.ones(self.nb_water_pix), 'power_minus_y', self.nb_water_pix, data, group=pixc)   
        
        data.add_variable('coherent_power', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'coherent power combination of minus_y and plus_y channel','units':'1','valid_min':'0','valid_max':'999999','coordinates':'longitude latitude','comment':'Power computed by combining the plus_y and minus_y channels coherently by coaligning the phases (arbitrary units that give sigma0 when noise substracted and normalized by the X factor'}
        data.add_variable_attributes('coherent_power', dic, group=pixc)
        fill_vector_param(np.ones(self.nb_water_pix), 'coherent_power', self.nb_water_pix, data, group=pixc)
 
        data.add_variable('x_factor_plus_y', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'X factor for plus_y channel power','units':'1','valid_min':'0','valid_max':'999999','coordinates':'longitude latitude','comment':'X factor for the plus_y channel power in linear units (arbitrary units to normalize noise-substracted power to sigma0)'}
        data.add_variable_attributes('x_factor_plus_y', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'x_factor_plus_y', self.nb_water_pix, data, group=pixc)
        
        data.add_variable('x_factor_minus_y', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'X factor for minus_y channel power','units':'1','valid_min':'0','valid_max':'999999','coordinates':'longitude latitude','comment':'X factor for the minus_y channel power in linear units (arbitrary units to normalize noise-substracted power to sigma0'}
        data.add_variable_attributes('x_factor_minus_y', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'x_factor_minus_y', self.nb_water_pix, data, group=pixc)  
        
        #-----------

        data.add_variable('water_frac', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'water fraction','units':'1','valid_min':'-999999','valid_max':'999999','coordinates':'longitude latitude','comment':'Noisy estimate of the fraction of the pixel that is water'}
        data.add_variable_attributes('water_frac', dic, group=pixc)
        fill_vector_param(self.water_frac, 'water_frac', self.nb_water_pix, data, group=pixc)       
        
        data.add_variable('water_frac_uncert', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'water fraction uncertaincy','units':'1','valid_min':'0','valid_max':'999999','coordinates':'longitude latitude','comment':'Uncertaincy estimate of the water fraction estimate (width of noisy water frac estimate distribution)'}
        data.add_variable_attributes('water_frac_uncert', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'water_frac_uncert', self.nb_water_pix, data, group=pixc)              
        
        data.add_variable('classification', np.int8, 'points', my_var.FV_NETCDF["int8"], compress, group=pixc)
        dic={'long_name':'classification','flag_meanings':'land land_near_water water_near_land open_water land_near_dark_water dark_water_edge dark_water','flag_meanings':'1 2 3 4 22 23 24','valid_min':'1','valid_max':'24','coordinates':'longitude latitude','comment':'Flags indicating water detection results'}
        data.add_variable_attributes('classification', dic, group=pixc)
        fill_vector_param(self.classification, 'classification', self.nb_water_pix, data, group=pixc) 
        
        data.add_variable('false_detection_rate', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'false detection rate','units':'1','valid_min':'0','valid_max':'1','coordinates':'longitude latitude','comment':'Probability of falsely detecting water when there is none'}
        data.add_variable_attributes('false_detection_rate', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'false_detection_rate', self.nb_water_pix, data, group=pixc)
        
        data.add_variable('missed_detection_rate', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'missed detection rate','units':'1','valid_min':'0','valid_max':'1','coordinates':'longitude latitude','comment':'Probability of falsely detecting no water when there is water'}
        data.add_variable_attributes('missed_detection_rate', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'missed_detection_rate', self.nb_water_pix, data, group=pixc)
        
        data.add_variable('prior_water_prob', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'prior water probability','units':'1','valid_min':'0','valid_max':'1','coordinates':'longitude latitude','comment':'Prior probability of water occuring'}
        data.add_variable_attributes('prior_water_prob', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'prior_water_prob', self.nb_water_pix, data, group=pixc)
        
        data.add_variable('bright_land_flag', np.int8, 'points', my_var.FV_NETCDF["int8"], compress, group=pixc)
        dic={'long_name':'bright land flag','standard_name':'status_flag','flag_meanings':'not_bright_land bright_land','flag_values':'0 1','valid_min':'0','valid_max':'1','coordinates':'longitude latitude','comment':'Flag indicating areas that are not typically water but are expected to be bright (e.g., urban areas, ice). this flag can be used to exclude detected water pixels in downstream processing'}
        data.add_variable_attributes('bright_land_flag', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'bright_land_flag', self.nb_water_pix, data, group=pixc)          
        
        data.add_variable('layover_impact', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'layover impact','units':'m','valid_min':'-999999','valid_max':'999999','coordinates':'longitude latitude','comment':'Estimate of the height error caused by layover, which may not be reliable on a pixel by pixel basis, but may be useful to augment aggregated height uncertainties'}
        data.add_variable_attributes('layover_impact', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'layover_impact', self.nb_water_pix, data, group=pixc)
        
        #--------------------

        data.add_variable('eff_num_rare_looks', np.int8, 'points', my_var.FV_NETCDF["int8"], compress, group=pixc)
        dic={}
        data.add_variable_attributes('eff_num_rare_looks', dic, group=pixc)
        fill_vector_param(np.full(self.nb_water_pix, 7.), 'eff_num_rare_looks', self.nb_water_pix, data, group=pixc) 
        
        #---------------------

        data.add_variable('latitude', np.float64, 'points', my_var.FV_NETCDF["float64"], compress, group=pixc)
        data.add_variable_attribute('latitude', 'units', 'degrees_north', group=pixc)
        dic={'long_name':'latitude (positive N, negative S)','standard_name':'latitude','units':'degrees_north','valid_min':'-80','valid_max':'80','comment':'Geodetic latitude [-80,80] (degrees north of equator) of the pixel'}
        data.add_variable_attributes('latitude', dic, group=pixc)
        fill_vector_param(self.latitude, 'latitude', self.nb_water_pix, data, group=pixc)

        data.add_variable('longitude', np.float64, 'points', my_var.FV_NETCDF["float64"], compress, group=pixc)
        data.add_variable_attribute('longitude', 'units', 'degrees_east', group=pixc)
        dic={'long_name':'longitude (degrees East)','standard_name':'longitude','units':'degrees_north','valid_min':'-180','valid_max':'180','comment':'Longitude [-180,180] (east of the Greenwich meridian) of the pixel'}
        data.add_variable_attributes('longitude', dic, group=pixc)
        fill_vector_param(self.longitude, 'longitude', self.nb_water_pix, data, group=pixc)

        data.add_variable('height', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        data.add_variable_attribute('height', 'units', 'm', group=pixc)
        dic={'long_name':'height above reference ellipsoid','units':'m','valid_min':'-1500','valid_max':'15000','coordinates':'longitude latitude','comment':'Approximate cross-track location of the pixel'}
        data.add_variable_attributes('height', dic, group=pixc)
        fill_vector_param(self.height, 'height', self.nb_water_pix, data, group=pixc)
        
        #---------------------

        data.add_variable('cross_track', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'approximate cross-track location','units':'m','valid_min':'0','valid_max':'75000','coordinates':'longitude latitude','comment':'Approximate cross-track location of the pixel'}
        data.add_variable_attributes('cross_track', dic, group=pixc)
        fill_vector_param(self.crosstrack, 'cross_track', self.nb_water_pix, data, group=pixc)

        data.add_variable('pixel_area', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'pixel area','units':'m²','valid_min':'0','valid_max':'999999','coordinates':'longitude latitude','comment':'Pixel area'}
        data.add_variable_attributes('pixel_area', dic, group=pixc)
        fill_vector_param(self.pixel_area, 'pixel_area', self.nb_water_pix, data, group=pixc) 

        data.add_variable('inc', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'incidence angle','units':'degrees','valid_min':'0','valid_max':'999999','coordinates':'longitude latitude','comment':'incidence angle'}
        data.add_variable_attributes('inc', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'inc', self.nb_water_pix, data, group=pixc)

        #----------------------
        
        data.add_variable('phase_noise_std', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'phase noise standard deviation','units':'degrees','valid_min':'0','valid_max':'999999','coordinates':'longitude latitude','comment':'Estimate of the phase noise standard deviation'}
        data.add_variable_attributes('phase_noise_std', dic, group=pixc)
        fill_vector_param(self.phase_noise_std, 'phase_noise_std', self.nb_water_pix, data, group=pixc)

        data.add_variable('dlatitude_dphase', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'sensitivity of latitude estimate to interferogram phase','units':'degrees/radians','valid_min':'-999999','valid_max':'999999','coordinates':'longitude latitude','comment':'sensitivity of the latitude estimate to the interferogram'}
        data.add_variable_attributes('dlatitude_dphase', dic, group=pixc)
        fill_vector_param(self.dlat_dphi, 'dlatitude_dphase', self.nb_water_pix, data, group=pixc)

        data.add_variable('dlongitude_dphase', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'sensitivity of longitude estimate to interferogram phase','units':'degrees/radian','valid_min':'-999999','valid_max':'999999','coordinates':'longitude latitude','comment':'sensitivity of the longitude estimate to the interferogram'}
        data.add_variable_attributes('dlongitude_dphase', dic, group=pixc)
        fill_vector_param(self.dlon_dphi, 'dlongitude_dphase', self.nb_water_pix, data, group=pixc)  

        data.add_variable('dheight_dphase', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dicdheight_dphase={'long_name':'sensitivity of height estimate to interferogram phase','units':'m/radian','valid_min':'-999999','valid_max':'999999','coordinates':'longitude latitude','comment':'sensitivity of the height estimate to the interferogram phase'}
        data.add_variable_attributes('dheight_dphase', dicdheight_dphase, group=pixc)
        fill_vector_param(self.dh_dphi, 'dheight_dphase', self.nb_water_pix, data, group=pixc) 

        data.add_variable('dheight_droll', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'sensitivity of height estimates to spacecraft roll','units':'m/degrees','valid_min':'-999999','valid_max':'999999','coordinates':'longitude latitude','comment':'sensitivity of the height estimate to the spacecraft roll'}
        data.add_variable_attributes('dheight_droll', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dheight_droll', self.nb_water_pix, data, group=pixc)

        data.add_variable('dheight_dbaseline', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'sensitivity of height estimate to interferometric baseline','units':'m/m','valid_min':'-999999','valid_max':'999999','coordinates':'longitude latitude','comment':'sensitivity of the height estimate to the interferometric baseline'}
        data.add_variable_attributes('dheight_dbaseline', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dheight_dbaseline', self.nb_water_pix, data, group=pixc)  

        data.add_variable('dheight_drange', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'sensitivity of height estimate to range (delay)','units':'m/m','valid_min':'-999999','valid_max':'999999','coordinates':'longitude latitude','comment':'sensitivity of the height estimate to the range (delay)'}
        data.add_variable_attributes('dheight_drange', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'dheight_drange', self.nb_water_pix, data, group=pixc) 

        data.add_variable('darea_dheight', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'sensitivity of pixel area to reference height','units':'m²/m','valid_min':'-999999','valid_max':'999999','coordinates':'longitude latitude','comment':'sensitivity of the pixel area to the reference height'}
        data.add_variable_attributes('darea_dheight', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'darea_dheight', self.nb_water_pix, data, group=pixc)

        #--------------------------
        
        data.add_variable('illumination_time', np.float64, 'points', my_var.FV_NETCDF["float64"], compress, group=pixc)
        dic={'long_name':'time of illumination of each pixel (UTC)','standard_name':'time','calendar':'gregorian','leap_second':'YYYY-MM-DD hh:mm:ss','units':'seconds since 2000-01-01 00:00:00.000','comment':'Time of measurement in seconds in the UTC time scale since 1 jan 2000 00:00:00 UTC. [tai_utc_difference] is the difference between TAI and UTC reference time (seconds) for the first measurement of the data set. If a leap second occurs within the data set, the attribute leap_second is set to the UTC time at which the leap second occurs'}
        data.add_variable_attributes('illumination_time', dic, group=pixc)
        fill_vector_param(self.computeTime_UTC(self.illumination_time), 'illumination_time', self.nb_water_pix, data, group=pixc)
        data.add_variable('illumination_time_tai', np.float64, 'points', my_var.FV_NETCDF["float64"], compress, group=pixc)
        dic={'long_name':'time of illumination of each pixel (TAI)','standard_name':'time','calendar':'gregorian','units':'seconds since 2000-01-01 00:00:00.000','comment':'Time of measurement in seconds in the TAI time scale since 1 jan 2000 00:00:00.000 TAI. This time scale contains no leap seconds. The difference (in seconds) with time in UTC is given by the attribute [illumination_time:tai_utc_difference]'}
        data.add_variable_attributes('illumination_time_tai', dic, group=pixc)
        fill_vector_param(self.computeTime_TAI(self.illumination_time), 'illumination_time_tai', self.nb_water_pix, data, group=pixc)  # TODO: to improve
        
        data.add_variable('eff_num_medium_looks', np.int32, 'points', my_var.FV_NETCDF["int32"], compress, group=pixc)
        dic={'long_name':'effective number of medium looks','units':'1','valid_min':'0','valid_max':'999999','coordinates':'longitude latitude','comment':'Effective number of independent looks taken in forming the medium interferogram (after adaptative averaging)'}
        data.add_variable_attributes('eff_num_medium_looks', dic, group=pixc)
        fill_vector_param(np.full(self.nb_water_pix, 36), 'eff_num_medium_looks', self.nb_water_pix, data, group=pixc)

        data.add_variable('sig0', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'sigma0','units':'1','valid_min':'-999999','valid_max':'999999','coordinates':'longitude latitude','comment':'Normalized radar cross section, or backscatter brightness'}
        data.add_variable_attributes('sig0', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'sig0', self.nb_water_pix, data, group=pixc)
        data.add_variable('phase_unwrapping_region', np.int32, 'points', my_var.FV_NETCDF["int32"], compress, group=pixc)
        dic={'long_name':'phase unwrapping region index','units':'1','valid_min':'-1','valid_max':'99999999','coordinates':'longitude latitude','comment':'Phase unwrapping region index'}
        data.add_variable_attributes('phase_unwrapping_region', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'phase_unwrapping_region', self.nb_water_pix, data, group=pixc)
        
        data.add_variable('instrument_range_cor', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'instrument range correction','units':'m','valid_min':'-999999','valid_max':'999999','coordinates':'longitude latitude','comment':'Term that incorporates all calibration corrections applied to range before geolocation'}
        data.add_variable_attributes('instrument_range_cor', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'instrument_range_cor', self.nb_water_pix, data, group=pixc)
        data.add_variable('instrument_phase_cor', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'instrument phase correction','units':'radians','valid_min':'-999999','valid_max':'999999','coordinates':'longitude latitude','comment':'Term that incorporates all calibration corrections applied to attitude before geolocation'}
        data.add_variable_attributes('instrument_phase_cor', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'instrument_phase_cor', self.nb_water_pix, data, group=pixc)
        data.add_variable('instrument_baseline_cor', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'instrument baseline correction','units':'m','valid_min':'-999999','valid_max':'999999','coordinates':'longitude latitude','comment':'Term that incorporates all calibration corrections applied to attitude before geolocation'}
        data.add_variable_attributes('instrument_baseline_cor', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'instrument_baseline_cor', self.nb_water_pix, data, group=pixc)
        data.add_variable('instrument_attitude_cor', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'instrument attitude correction','units':'degrees','valid_min':'-999999','valid_max':'999999','coordinates':'longitude latitude','comment':'Term that incorporates all calibration corrections applied to attiude before geolocation'}
        data.add_variable_attributes('instrument_attitude_cor', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'instrument_attitude_cor', self.nb_water_pix, data, group=pixc)

        data.add_variable('model_dry_tropo_cor', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'dry troposphere vertical correction','source':'European Centre for Medium-Range Weather Forecasts','institution':'ECMWF','units':'m','valid_min':'-3','valid_max':'-1.5','coordinates':'longitude latitude','comment':'Equivalent vertical correction due to wet troposphere delay.The reported pixel height, latitude and longitude are computed after adding negative media corrections to uncorrected range along slant-range paths, accounting for the differential delay between two KaRIn antennas. The equivalent vertical correction is computed by applying obliquity factors to the slant-path correction. Adding the reported correction to the reported pixel height results in the uncorrected pixel height'}
        data.add_variable_attributes('model_dry_tropo_cor', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'model_dry_tropo_cor', self.nb_water_pix, data, group=pixc)
        data.add_variable('model_wet_tropo_cor', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'wet troposphere vertical correction','source':'European Centre for Medium-Range Weather Forecasts','institution':'ECMWF','units':'m','valid_min':'-1','valid_max':'0','coordinates':'longitude latitude','comment':'Equivalent vertical correction due to wet troposphere delay. The reported pixel height, latitude and longitude are computed after adding negative media corrections to uncorrected range along slant-range paths, accounting for the differential delay between the two KaRIn antennas. The equivalent vertical correction is computed by applying obliquity factors to the slant-path correction. Adding the reported correction to the reported pixel height results in the uncorrected pixel height'}
        data.add_variable_attributes('model_wet_tropo_cor', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'model_wet_tropo_cor', self.nb_water_pix, data, group=pixc)
        data.add_variable('iono_cor_gim_ka', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'ionosphere vertical correction','source':'Global Ionosphere Maps','institution':'JPL','units':'m','valid_min':'-0.5','valid_max':'0','coordinates':'longitude latitude','comment':'Equivalent vertical correction due to ionosphere delay. The reported pixel height, latitude and longitude are computed after adding negative media corrections to uncorrected range along slant-range paths, accounting for the differential delay between the two KaRIn antennas. The equivalent vertical correction is computed by applying obliquity factors to the slant-path correction. Adding the reported correction to the reported pixel height results in the uncorrected pixel height'}
        data.add_variable_attributes('iono_cor_gim_ka', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'iono_cor_gim_ka', self.nb_water_pix, data, group=pixc)   
        data.add_variable('height_cor_xover', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'height correction from KaRIn crossovers','units':'m','valid_min':'-10','valid_max':'10','coordinates':'longitude latitude','comment':'Height correction from KaRIn crossover calibration. The correction is applied before geolocation but reported as an equivalent height correction'}
        data.add_variable_attributes('height_cor_xover', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'height_cor_xover', self.nb_water_pix, data, group=pixc)
        # ~ data.add_variable('height_cor_xover', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        # ~ fill_vector_param(np.zeros(self.nb_water_pix), 'height_cor_xover', self.nb_water_pix, data, group=pixc)        

        data.add_variable('geoid', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'geoid height','standard_name':'geoid_height_above_reference_ellipsoid','source':'EGM2008(Pavlis et al., 2012','units':'m','valid_min':'-150','valid_max':'150','coordinates':'longitude latitude','comment':'Geoid height above the reference ellipsoid with a correction to refer the value to the mean tide system, i.e. includes the permanent tide (zero frequency)'}
        data.add_variable_attributes('geoid', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'geoid', self.nb_water_pix, data, group=pixc)
        data.add_variable('solid_earth_tide', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'solid Earth tide height','source':'Cartwright and Taylor (1971) and Cartwright and Edden (1973)','units':'m','valid_min':'-1','valid_max':'1','coordinates':'longitude latitude','comment':'Geoid height above the reference ellipsoid with a correction to refer the value to the main tide system, i.e. includes the permanent tide (zero frequency)'}
        data.add_variable_attributes('solid_earth_tide', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'solid_earth_tide', self.nb_water_pix, data, group=pixc)
        data.add_variable('load_tide_sol1', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'geocentric load tide height from model 1','source':'FES2014b(Carrere et al.,2016)','institution':'LEGOS/CNES','units':'m','valid_min':'-0.2','valid_max':'0.2','coordinates':'longitude latitude','comment':'Geocentric load tide height. The effect of the ocean tide loading of the earths crust. This value is reported for referenced for reference but is not applied to the reported height'}
        data.add_variable_attributes('load_tide_sol1', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'load_tide_sol1', self.nb_water_pix, data, group=pixc)
        data.add_variable('load_tide_sol2', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'geocentric load tide height from model 2','source':'GOT4.10c(Ray, 2013)','institution':'GSFC','units':'m','valid_min':'-0.2','valid_max':'0.2','coordinates':'longitude latitude','comment':'Geocentric load tide height. The effect of the ocean tide loading of the Earth’s crust. This value is reported for reference but is not applied to the reported height'}
        data.add_variable_attributes('load_tide_sol2', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'load_tide_sol2', self.nb_water_pix, data, group=pixc)
        data.add_variable('pole_tide', np.float32, 'points', my_var.FV_NETCDF["float32"], compress, group=pixc)
        dic={'long_name':'geocentric pole tide height','source':'Wahr(1985) and Desai et al.(2015)','units':'m','valid_min':'-0.2','valid_max':'0.2','coordinates':'longitude latitude','comment':'Geocentric pole tide height. The total of the contribution from the solid-Earth (body) pole tide height and the load pole tide height (i.e., the effect of the ocean pole tide loading of the Earth’s crust)'}
        data.add_variable_attributes('pole_tide', dic, group=pixc)
        fill_vector_param(np.zeros(self.nb_water_pix), 'pole_tide', self.nb_water_pix, data, group=pixc)
        data.add_variable('pixc_qual', np.int8, 'points', my_var.FV_NETCDF["int8"], compress, group=pixc)
        dic={'standard_name':'status_flag','flag_meanings':'good bad','flag_values':'0 1','valid_min':'0','valid_max':'1','coordinates':'longitude latitude','comment':'Quality flag dor pixel cloud data'}
        data.add_variable_attributes('pixc_qual', dic, group=pixc)
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
        dic={'long_name':'time in UTC','standard_name':'time','calendar':'gregorian','tai_utc_difference':'[Value of TAI - UTC at first record]','leap_second':'YYYY-MM-DD hh:mm:ss','units':'seconds since 2000-01-01 00:00:00.000','comment':'Time of measurement in seconds in the UTC time scale since 1 Jan 2000 00:00:00 UTC. [tai_utc_difference] is the difference between TAI and UTC reference time (seconds) for the first measurement of the data set. If a leap second occurs within the data set, the attribute leap_second is set to the UTC time at which the leap second occurs'}
        data.add_variable_attributes('time', dic, group=sensor)
        fill_vector_param(self.computeTime_UTC(self.nadir_time), 'time', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('time_tai', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'time in TAI','standard_name':'time','calendar':'gregorian','units':'second since 2000-01-01 00:00:00.000','comment':'Time of measurement in seconds in the UTC time scale since 1 Jan 2000 00:00:00 UTC. [tai_utc_difference] is the difference between TAI and UTC reference time (seconds) for the first measurement of the data set. If a leap second occurs within the data set, the attribute leap_second is set to the UTC time at which the leap second occurs'}
        data.add_variable_attributes('time_tai', dic, group=sensor)
        fill_vector_param(self.computeTime_TAI(self.nadir_time), 'time_tai', self.nb_nadir_pix, data, group=sensor)
        
        data.add_variable('latitude', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'latitude (positive N, negative S) of the spacecraft','standard_name':'latitude','units':'degrees_north','valid_min':'-80.0','valid_max':'80.0','comment':'Geodetic latitude of the KMSF origin with respect to the reference ellipsoid'}
        data.add_variable_attributes('latitude', dic, group=sensor)
        fill_vector_param(self.nadir_latitude, 'latitude', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('longitude', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'longitude (degrees east) of the spacecraft','standard_name':'longitude','units':'degrees_east','valid_min':'-180.0','valid_max':'180.0','comment':'Longitude of the KMSF origin with positive values indicating longitudes east of the greenwich meridian'}
        data.add_variable_attributes('longitude', dic, group=sensor)
        fill_vector_param(self.nadir_longitude, 'longitude', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('altitude', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'altitude of the spacecraft','units':'m','valid_min':'0.0','valid_max':'10000000.0','coordinates':'longitude latitude','comment':'Altitude above the reference ellipsoid of the KMSF origin'}
        data.add_variable_attributes('altitude', dic, group=sensor)
        fill_vector_param(self.nadir_altitude, 'altitude', self.nb_nadir_pix, data, group=sensor)
        
        data.add_variable('roll', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'roll of the spacecraft','units':'degrees','valid_min':'-180','valid_max':'180','coordinates':'longitude latitude','comment':'KMSF attitude roll angle; positive values move the +y antenna down'}
        data.add_variable_attributes('roll', dic, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'roll', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('pitch', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'pitch of the spacecraft','units':'degrees','valid_min':'-180','valid_max':'180','coordinates':'longitude latitude','comment':'KMSF attitude pitch angle; positive values move the KMSF +x axis up'}
        data.add_variable_attributes('pitch', dic, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'pitch', self.nb_nadir_pix, data, group=sensor)  
        data.add_variable('yaw', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'yaw of the spacecraft','units':'degrees','valid_min':'-180','valid_max':'180','coordinates':'longitude latitude','comment':'KMSF attitude yaw angle relative to the nadir track. The yaw angle is a right-handed rotation about the nadir (downward) direction. A yaw value of 0 deg indicates that the KMSF +x axis is aligned with the horizontal component of the Earth-relative velocity vector. A yaw value of 180 deg indicates that the spacecraft is in a yaw-flipped state, with the KMSF -x axis aligned with the horizontal component of the Earth-relative velocity vector'}
        data.add_variable_attributes('yaw', dic, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'yaw', self.nb_nadir_pix, data, group=sensor)   
        data.add_variable('velocity_heading', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'heading of the spacecraft Earth-relative velocity vector','units':'degrees','valid_min':'0','valid_max':'360','coordinates':'longitude latitude','comment':'KMSF attitude yaw angle relative to the nadir track. The yaw angle is a right-handed rotation about the nadir (downward) direction. A yaw value of 0 deg indicates that the KMSF +x axis is aligned with the horizontal component of the Earth-relative velocity vector. A yaw value of 180 deg indicates that the spacecraft is in a yaw-flipped state, with the KMSF -x axis aligned with the horizontal component of the Earth-relative velocity vector'}
        data.add_variable_attributes('velocity_heading', dic, group=sensor)
        fill_vector_param(self.nadir_heading, 'velocity_heading', self.nb_nadir_pix, data, group=sensor)
        
        data.add_variable('x', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'x coordinates of the KMSF in the ECEF frame','units':'m','valid_min':'-10000000.0','valid_max':'10000000.0','comment':'x coordinate of the KMSF origin in ECEF frame'}
        data.add_variable_attributes('x', dic, group=sensor)
        fill_vector_param(self.nadir_x, 'x', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('y', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'y coordinates of the KMSF in the ECEF frame','units':'m','valid_min':'-10000000.0','valid_max':'10000000.0','comment':'y coordinate of the KMSF origin in ECEF frame'}
        data.add_variable_attributes('y', dic, group=sensor)
        fill_vector_param(self.nadir_y, 'y', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('z', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'z coordinates of the KMSF in the ECEF frame','units':'m','valid_min':'-10000000.0','valid_max':'10000000.0','comment':'z coordinate of the KMSF origin in ECEF frame'}
        data.add_variable_attributes('z', dic, group=sensor)
        fill_vector_param(self.nadir_z, 'z', self.nb_nadir_pix, data, group=sensor)
    
        data.add_variable('vx', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'x component of the spacecraft velocity in the ECEF frame','units':'m/s','valid_min':'-10000.0','valid_max':'10000000.0','comment':'KMSF velocity component in x direction in the ECEF'}
        data.add_variable_attributes('vx', dic, group=sensor)
        fill_vector_param(self.nadir_vx, 'vx', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('vy', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'y component of the spacecraft velocity in the ECEF frame','units':'m/s','valid_min':'-10000.0','valid_max':'10000000.0','comment':'KMSF velocity component in y direction in the ECEF'}
        data.add_variable_attributes('vy', dic, group=sensor)
        fill_vector_param(self.nadir_vy, 'vy', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('vz', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'z component of the spacecraft velocity in the ECEF frame','units':'m/s','valid_min':'-10000.0','valid_max':'10000000.0','comment':'KMSF velocity component in z direction in the ECEF'}
        data.add_variable_attributes('vz', dic, group=sensor)
        fill_vector_param(self.nadir_vz, 'vz', self.nb_nadir_pix, data, group=sensor)
        
        data.add_variable('plus_y_antenna_x', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'x coordinate of plus_y antenna phase center in the ECEF frame','units':'m','valid_min':'-10000000.0','valid_max':'10000000.0','comment':'x coordinate of the plus_y antenna phase center in the ECEF frame'}
        data.add_variable_attributes('plus_y_antenna_x', dic, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'plus_y_antenna_x', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('plus_y_antenna_y', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'y coordinate of plus_y antenna phase center in the ECEF frame','units':'m','valid_min':'-10000000.0','valid_max':'10000000.0','comment':'x coordinate of the plus_y antenna phase center in the ECEF frame'}
        data.add_variable_attributes('plus_y_antenna_y', dic, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'plus_y_antenna_y', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('plus_y_antenna_z', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'z coordinate of plus_y antenna phase center in the ECEF frame','units':'m','valid_min':'-10000000.0','valid_max':'10000000.0','comment':'z coordinate of the plus_y antenna phase center in the ECEF frame'}
        data.add_variable_attributes('plus_y_antenna_z', dic, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'plus_y_antenna_z', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('minus_y_antenna_x', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'x coordinate of the minus_y antenna phase center in the ECEF frame','units':'m','valid_min':'-10000000.0','valid_max':'10000000.0','comment':'x coordinate of the minus_y antenna phase center in the ECEF frame'}
        data.add_variable_attributes('minus_y_antenna_x', dic, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'minus_y_antenna_x', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('minus_y_antenna_y', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'y coordinate of the minus_y antenna phase center in the ECEF frame','units':'m','valid_min':'-10000000.0','valid_max':'10000000.0','comment':'y coordinate of the minus_y antenna phase center in the ECEF frame'}
        data.add_variable_attributes('minus_y_antenna_y', dic, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'minus_y_antenna_y', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('minus_y_antenna_z', np.float64, 'num_tvps', my_var.FV_NETCDF["float64"], compress, group=sensor)
        dic={'long_name':'z coordinate of the minus_y antenna phase center in the ECEF frame','units':'m','valid_min':'-10000000.0','valid_max':'10000000.0','comment':'z coordinate of the minus_y antenna phase center in the ECEF frame'}
        data.add_variable_attributes('minus_y_antenna_z', dic, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'minus_y_antenna_z', self.nb_nadir_pix, data, group=sensor)
        data.add_variable('record_counter', np.int32, 'num_tvps', my_var.FV_NETCDF["int32"], compress, group=sensor)
        dic={'long_name':'record counter','units':'1','valid_min':'1','valid_max':'999999999','coordinates':'longitude latitude','comment':'Index of the TVP record used to align data samples across granules'}
        data.add_variable_attributes('record_counter', dic, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'record_counter', self.nb_nadir_pix, data, group=sensor) 
        data.add_variable('sc_event_flag', np.int8, 'num_tvps', my_var.FV_NETCDF["int8"], compress, group=sensor)
        dic={'standard_name':'status_flag','flag_meanings':'nominal not_nominal','flag_values':'0 1','valid_min':'0','valid_max':'1','coordinates':'longitude latitude','comment':'spacecraft event flag'}
        data.add_variable_attributes('sc_event_flag', dic, group=sensor)
        fill_vector_param(np.zeros(self.nb_nadir_pix), 'sc_event_flag', self.nb_nadir_pix, data, group=sensor) 
        data.add_variable('tvp_qual', np.int8, 'num_tvps', my_var.FV_NETCDF["int8"], compress, group=sensor)
        dic={'standard_name':'status_flag','flag_meanings':'good bad','flag_values':'0 1','valid_min':'0','valid_max':'1','coordinates':'longitude latitude','comment':'spacecraft event flag'}
        data.add_variable_attributes('tvp_qual', dic, group=sensor)
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
        dic={'long_name':'Noise estimate for the plus_y channel','units':'1','valid_min':'0.0','valid_max':'1e+20','comment':'Noise power for the plus_y channel expressed in linear units (arbitrary units that give noise-equivalent sigma0 when normalized by the X factor)'}
        data.add_variable_attributes('noise_plus_y', dic, group=noise)
        fill_vector_param(np.full(self.nb_nadir_pix, -116.845780895788), 'noise_plus_y', self.nb_nadir_pix, data, group=noise)
        data.add_variable('noise_minus_y', np.float32, 'num_lines', my_var.FV_NETCDF["float32"], compress, group=noise)
        dic={'long_name':'Noise estimate for the minus_y channel','units':'1','valid_min':'0.0','valid_max':'1e+20','comment':'Noise power for the minus_y channel expressed in linear units (arbitrary units that give noise-equivalent sigma0 when normalized by the X factor)'}
        data.add_variable_attributes('noise_minus_y', dic, group=noise)
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
        tmpField = ogr.FieldDefn(str('phi_std'), ogr.OFTReal)  # Phase noise standart deviation 
        tmpField.SetWidth(15)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('dlat_dph'), ogr.OFTReal)  # latitude error relatively to the phase
        tmpField.SetWidth(15)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('dlon_dph'), ogr.OFTReal)  # Longitude error relatively to the phase
        tmpField.SetWidth(15)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('dh_dphi'), ogr.OFTReal)  # Height error relatively to the phase
        tmpField.SetWidth(15)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        # 1.5 - On recupere la definition de la couche
        outLayerDefn = outLayer.GetLayerDefn()
        
        # 2 - On traite point par point
        for az_ind, range_index, classif, pixel_area, lat, lng, height, crosstrack, phase_noise_std, dlat_dphi, dlon_dphi, dh_dphi in zip(self.azimuth_index, self.range_index, self.classification, \
        self.pixel_area, self.latitude, self.longitude, self.height, self.crosstrack, self.phase_noise_std, self.dlat_dphi, self.dlon_dphi, self.dh_dphi):
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
            outFeature.SetField(str('phi_std'), float(phase_noise_std))
            outFeature.SetField(str('dlat_dph'), float(dlat_dphi))
            outFeature.SetField(str('dlon_dph'), float(dlon_dphi))
            outFeature.SetField(str('dh_dphi'), float(dh_dphi))
         
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

