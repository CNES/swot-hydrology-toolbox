# -*- coding: utf-8 -*-
"""
.. module proc_pixc.py
    :synopsis: Deal with L2_HR_PIXC (pixel cloud) files
    Created on 08/24/2017
    Totally modified on 04/14/2020 to use JPL XML description files

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

import product_netcdf as my_nc


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
        my_api.printInfo("[proc_pixc] == INIT ==")

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

    def write_pixc_file(self, IN_output_file):
        """
        Write the main file of real pixel cloud product (L2_HR_PIXC product, main file)

        :param IN_output_file: output full path
        :type IN_output_file: string
        """
        my_api.printInfo("[proc_pixc] == write_pixc_file : %s ==" % IN_output_file)
    
        # 1 - Init a PixcProduct object
        nc_writer = my_nc.PixcProduct()
        
        # 2 - Update global attributes
        tmp_metadata = {}
        tmp_metadata['cycle_number'] = self.cycle_num
        tmp_metadata['pass_number'] = np.int(self.pass_num)
        tmp_metadata['tile_number'] = int(self.tile_ref[0:-1])
        tmp_metadata['swath_side'] = self.tile_ref[-1]
        tmp_metadata['tile_name'] = "%03d_%03d%s" % (np.int(self.pass_num), int(self.tile_ref[0:-1]), self.tile_ref[-1])
        tmp_metadata['near_range'] = np.min(self.near_range)  # TODO: improve
        tmp_metadata['nominal_slant_range_spacing'] = self.range_spacing
        tmp_metadata['time_coverage_start'] = self.computeDate(self.nadir_time[0])
        tmp_metadata['time_coverage_end'] = self.computeDate(self.nadir_time[-1])
        tmp_metadata['geospatial_lon_min'] = min([self.inner_first[0], self.inner_last[0], self.outer_first[0], self.outer_last[0]])
        tmp_metadata['geospatial_lon_max'] = max([self.inner_first[0], self.inner_last[0], self.outer_first[0], self.outer_last[0]])
        tmp_metadata['geospatial_lat_min'] = min([self.inner_first[1], self.inner_last[1], self.outer_first[1], self.outer_last[1]])
        tmp_metadata['geospatial_lat_max'] = max([self.inner_first[1], self.inner_last[1], self.outer_first[1], self.outer_last[1]])
        tmp_metadata["inner_first_longitude"] = self.inner_first[0]
        tmp_metadata["inner_first_latitude"] = self.inner_first[1]
        tmp_metadata["inner_last_longitude"] = self.inner_last[0]
        tmp_metadata["inner_last_latitude"] = self.inner_last[1]
        tmp_metadata["outer_first_longitude"] = self.outer_first[0]
        tmp_metadata["outer_first_latitude"] = self.outer_first[1]
        tmp_metadata["outer_last_longitude"] = self.outer_last[0]
        tmp_metadata["outer_last_latitude"] = self.outer_last[1]
        nc_writer.set_metadata_val(tmp_metadata)

        # =======================
        # == Group pixel_cloud ==
        # =======================
        
        # Update some of group attributes   
        pixel_cloud_metadata = {}
        pixel_cloud_metadata['interferogram_size_azimuth'] = self.nb_pix_azimuth
        pixel_cloud_metadata['interferogram_size_range'] = self.nb_pix_range
        nc_writer.set_metadata_val(pixel_cloud_metadata, group="pixel_cloud")
  
        # Update group dimension
        nc_writer.set_dim_val("points", self.nb_water_pix)
        
        # Create dictionary with value of variables
        pixel_cloud_vars_val = {}
        pixel_cloud_vars_val['azimuth_index'] = self.azimuth_index
        pixel_cloud_vars_val['range_index'] = self.range_index
        #--------------------
        pixel_cloud_vars_val['interferogram'] = np.array(self.interferogram)
        pixel_cloud_vars_val['power_plus_y'] = np.ones(self.nb_water_pix)
        pixel_cloud_vars_val['power_minus_y'] = np.ones(self.nb_water_pix)
        pixel_cloud_vars_val['coherent_power'] = np.ones(self.nb_water_pix)
        pixel_cloud_vars_val['x_factor_plus_y'] = np.ones(self.nb_water_pix)
        pixel_cloud_vars_val['x_factor_minus_y'] = np.ones(self.nb_water_pix)
        #--------------------
        pixel_cloud_vars_val['water_frac'] = self.water_frac
        pixel_cloud_vars_val['water_frac_uncert'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['classification'] = self.classification
        pixel_cloud_vars_val['false_detection_rate'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['missed_detection_rate'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['prior_water_prob'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['bright_land_flag'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['layover_impact'] = np.zeros(self.nb_water_pix)
        #--------------------
        pixel_cloud_vars_val['eff_num_rare_looks'] = np.full(self.nb_water_pix, 7.)
        #--------------------
        pixel_cloud_vars_val['latitude'] = self.latitude
        pixel_cloud_vars_val['longitude'] = self.longitude
        pixel_cloud_vars_val['height'] = self.height
        #--------------------
        pixel_cloud_vars_val['cross_track'] = self.crosstrack
        pixel_cloud_vars_val['pixel_area'] = self.pixel_area
        pixel_cloud_vars_val['inc'] = np.zeros(self.nb_water_pix)
        #--------------------
        pixel_cloud_vars_val['phase_noise_std'] = self.phase_noise_std
        pixel_cloud_vars_val['dlatitude_dphase'] = self.dlat_dphi
        pixel_cloud_vars_val['dlongitude_dphase'] = self.dlon_dphi
        pixel_cloud_vars_val['dheight_dphase'] = self.dh_dphi
        pixel_cloud_vars_val['dheight_droll'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['dheight_dbaseline'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['dheight_drange'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['darea_dheight'] = np.zeros(self.nb_water_pix)
        #--------------------
        pixel_cloud_vars_val['illumination_time'] = self.computeTime_UTC(self.illumination_time)
        pixel_cloud_vars_val['illumination_time_tai'] = self.computeTime_TAI(self.illumination_time)
        pixel_cloud_vars_val['eff_num_medium_looks'] = np.full(self.nb_water_pix, 36)
        pixel_cloud_vars_val['sig0'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['phase_unwrapping_region'] = np.zeros(self.nb_water_pix)
        #--------------------
        pixel_cloud_vars_val['instrument_range_cor'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['instrument_phase_cor'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['instrument_baseline_cor'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['instrument_attitude_cor'] = np.zeros(self.nb_water_pix)
        #--------------------
        pixel_cloud_vars_val['model_dry_tropo_cor'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['model_wet_tropo_cor'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['iono_cor_gim_ka'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['height_cor_xover'] = np.zeros(self.nb_water_pix)
        #--------------------
        pixel_cloud_vars_val['geoid'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['solid_earth_tide'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['load_tide_sol1'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['load_tide_sol2'] = np.zeros(self.nb_water_pix)
        pixel_cloud_vars_val['pole_tide'] = np.zeros(self.nb_water_pix)
        #--------------------
        pixel_cloud_vars_val['pixc_qual'] = np.zeros(self.nb_water_pix)
        
        # ===============
        # == Group TVP ==
        # ===============
  
        # Update group dimension
        nc_writer.set_dim_val("num_tvps", self.nb_nadir_pix)
        
        # Create dictionary with value of variables
        tvp_vars_val = {}
        tvp_vars_val['time'] = self.computeTime_UTC(self.nadir_time)
        tvp_vars_val['time_tai'] = self.computeTime_TAI(self.nadir_time)
        #--------------------
        tvp_vars_val['latitude'] = self.nadir_latitude
        tvp_vars_val['longitude'] = self.nadir_longitude
        tvp_vars_val['altitude'] = self.nadir_altitude
        #--------------------
        tvp_vars_val['roll'] = np.zeros(self.nb_nadir_pix)
        tvp_vars_val['pitch'] = np.zeros(self.nb_nadir_pix)
        tvp_vars_val['yaw'] = np.zeros(self.nb_nadir_pix)
        tvp_vars_val['velocity_heading'] = self.nadir_heading
        #--------------------
        tvp_vars_val['x'] = self.nadir_x
        tvp_vars_val['y'] = self.nadir_y
        tvp_vars_val['z'] = self.nadir_z
        #--------------------
        tvp_vars_val['vx'] = self.nadir_vx
        tvp_vars_val['vy'] = self.nadir_vy
        tvp_vars_val['vz'] = self.nadir_vz
        #--------------------
        tvp_vars_val['plus_y_antenna_x'] = np.zeros(self.nb_nadir_pix)
        tvp_vars_val['plus_y_antenna_y'] = np.zeros(self.nb_nadir_pix)
        tvp_vars_val['plus_y_antenna_z'] = np.zeros(self.nb_nadir_pix)
        tvp_vars_val['minus_y_antenna_x'] = np.zeros(self.nb_nadir_pix)
        tvp_vars_val['minus_y_antenna_y'] = np.zeros(self.nb_nadir_pix)
        tvp_vars_val['minus_y_antenna_z'] = np.zeros(self.nb_nadir_pix)
        #--------------------
        tvp_vars_val['record_counter'] = np.zeros(self.nb_nadir_pix)
        tvp_vars_val['sc_event_flag'] = np.zeros(self.nb_nadir_pix)
        tvp_vars_val['tvp_qual'] = np.zeros(self.nb_nadir_pix)
                
        # =================
        # == Group Noise ==
        # =================
  
        # Update group dimension
        nc_writer.set_dim_val("num_lines", self.nb_nadir_pix)
        
        # Create dictionary with value of variables
        noise_vars_val = {}
        noise_vars_val['noise_plus_y'] = np.full(self.nb_nadir_pix, -116.845780895788)
        noise_vars_val['noise_minus_y'] = np.full(self.nb_nadir_pix, -116.845780895788)
                
        # =====================
        # == Write PIXC file ==
        # =====================

        # Group dictionaries of group variables in a parent dictionary
        all_vars_val_per_group = {}
        all_vars_val_per_group["pixel_cloud"] = pixel_cloud_vars_val
        all_vars_val_per_group["tvp"] = tvp_vars_val
        all_vars_val_per_group["noise"] = noise_vars_val
        
        # Write the PIXC product
        nc_writer.write_product(IN_output_file, in_vars_value=all_vars_val_per_group)
    
    #----------------------------------
 
    def write_annotation_file(self, IN_output_file, IN_pixc_file):
        """
        write the river-annotation.rdf file so that lake processor can run
        
        :param IN_output_file: output full path
        :type IN_output_file: string
        :param IN_pixc_file: PIXC full path
        :type IN_pixc_file: string
        """
        my_api.printInfo("[proc_pixc] == write_annotation_file : %s ==" % IN_output_file)
        
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
        my_api.printInfo("[proc_pixc] == write_pixc_asShp : %s ==" % IN_output_file) 
        
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
        my_api.printInfo("[proc_pixc] == write_tvp_asShp : %s ==" % IN_output_file) 
    
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

