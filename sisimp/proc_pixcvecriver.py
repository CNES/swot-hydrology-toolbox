# -*- coding: utf-8 -*-
"""
.. module proc_pixcvecriver.py
    :synopsis: Deal with L2_HR_PIXCVecRiver files (data file complementary to L2_HR_PIXC pixel cloud files)
    Created on 08/29/2017
    Totally modified on 04/15/2020 to use JPL XML description files
    
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


class l2_hr_pixc_vec_river(object):

    def __init__(self, IN_azimuth_index, IN_range_index,
                 IN_mission_start_time, IN_cycle_duration, IN_cycle_num, IN_pass_num, IN_tile_ref, IN_start_time, IN_stop_time, 
                 IN_nb_pix_range, IN_nb_pix_azimuth, IN_tile_coords):
        """
        Constructor of the pixel cloud vector attribute product

        :param IN_azimuth_index: azimuth indices
        :type IN_azimuth_index: 1D-array of int
        :param IN_range_index: range indices
        :type IN_range_index: 1D-array of int
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
            
        + nb_water_pix(int) : number of water pixels, i.e. pixels in azimuth_index, ..., crosstrack vectors
        + latitude_vectorproc(1D-array) : improved latitude values
        + longitude_vectorproc(1D-array) : improved longitude values
        + height_vectorproc(1D-array) : improved height values
        + river_tag(1D-array): reach ID associated to each pixel
        + pattern(str): filename pattern
        """
        my_api.printInfo("[l2_hr_pixc_vec_river] == INIT ==") 
        
        self.azimuth_index = IN_azimuth_index
        self.range_index = IN_range_index
        self.nb_water_pix = IN_azimuth_index.size
        self.latitude_vectorproc = None
        self.longitude_vectorproc = None
        self.height_vectorproc = None
        self.river_tag = None
        
        self.mission_start_time = IN_mission_start_time
        self.cycle_duration = IN_cycle_duration
        self.cycle_num = IN_cycle_num
        self.pass_num = IN_pass_num
        self.tile_ref = IN_tile_ref
        self.start_time = IN_start_time
        self.stop_time = IN_stop_time
        self.nb_pix_range = IN_nb_pix_range
        self.nb_pix_azimuth = IN_nb_pix_azimuth

        (inner_first, inner_last, outer_first, outer_last) = IN_tile_coords
        self.inner_first = inner_first
        self.inner_last = inner_last
        self.outer_first = outer_first
        self.outer_last = outer_last
        
        self.nb_pix_river = 0  # Number of PixC associated to river
        self.pixc_river_idx = None  # Indices of pixels in the L2_HR_PIXC product associated to a river

    #----------------------------------

    def write_file(self, IN_output_file):
        """
        Write the pixel cloud vector attribute for river product (L2_HR_PIXCVecRiver product)
        NB: if there is no river pixels, write an empty file (with points = 0)

        :param IN_output_file: output full path
        :type IN_output_file: string
        """ 
        my_api.printInfo("[l2_hr_pixc_vec_river] == write_file : %s ==" % IN_output_file)
    
        # 1 - Init a PixcVecRiverProduct object
        nc_writer = my_nc.PixcVecRiverProduct()
        
        # 2 - Get river pixels indices
        self.pixc_river_idx = np.where(self.river_tag != "")[0]  # Indices of river pixels in the L2_HR_PIXC product
        self.nb_pix_river = len(self.pixc_river_idx)  # Number of river pixels
        
        # 3 - Update global attributes
        tmp_metadata = {}
        tmp_metadata['cycle_number'] = int(self.cycle_num)
        tmp_metadata['pass_number'] = int(self.pass_num)
        tmp_metadata['tile_number'] = int(self.tile_ref[0:-1])
        tmp_metadata['swath_side'] = self.tile_ref[-1]
        tmp_metadata['tile_name'] = "%03d_%03d%s" % (np.int(self.pass_num), int(self.tile_ref[0:-1]), self.tile_ref[-1])
        tmp_metadata['interferogram_size_range'] = self.nb_pix_range 
        tmp_metadata['interferogram_size_azimuth'] = self.nb_pix_azimuth
        tmp_metadata['start_time'] = self.computeDate(self.start_time)
        tmp_metadata['stop_time'] = self.computeDate(self.stop_time)
        tmp_metadata["inner_first_longitude"] = self.inner_first[0]
        tmp_metadata["inner_first_latitude"] = self.inner_first[1]
        tmp_metadata["inner_last_longitude"] = self.inner_last[0]
        tmp_metadata["inner_last_latitude"] = self.inner_last[1]
        tmp_metadata["outer_first_longitude"] = self.outer_first[0]
        tmp_metadata["outer_first_latitude"] = self.outer_first[1]
        tmp_metadata["outer_last_longitude"] = self.outer_last[0]
        tmp_metadata["outer_last_latitude"] = self.outer_last[1]
        nc_writer.set_metadata_val(tmp_metadata)
        
        # 4 - Update dimension
        nc_writer.set_dim_val("points", self.nb_pix_river)
        
        # 5 - Create dictionary with value of variables
        pixcvecriver_vars_val = {}
        pixcvecriver_vars_val['azimuth_index'] = self.azimuth_index[self.pixc_river_idx]
        pixcvecriver_vars_val['range_index'] = self.range_index[self.pixc_river_idx]
        #--------------------
        pixcvecriver_vars_val['latitude_vectorproc'] = self.latitude_vectorproc[self.pixc_river_idx]
        pixcvecriver_vars_val['longitude_vectorproc'] = self.longitude_vectorproc[self.pixc_river_idx]
        pixcvecriver_vars_val['height_vectorproc'] = self.height_vectorproc[self.pixc_river_idx]
        #--------------------
        pixcvecriver_vars_val['reach_id'] = self.river_tag[self.pixc_river_idx]
        pixcvecriver_vars_val['node_id'] = self.river_tag[self.pixc_river_idx]
        #--------------------
        pixcvecriver_vars_val['ice_clim_f'] = np.zeros(self.nb_pix_river)
        pixcvecriver_vars_val['ice_dyn_f'] = np.zeros(self.nb_pix_river)
        #--------------------
        pixcvecriver_vars_val['pixc_index'] = self.pixc_river_idx
        pixcvecriver_vars_val['lake_flag'] = np.zeros(self.nb_pix_river)
        
        # 6 - Write the PIXC product
        nc_writer.write_product(IN_output_file, in_vars_value=pixcvecriver_vars_val)

    def write_file_asShp(self, IN_output_file):
        """
        Write the pixel cloud vector attribute for rivers file as a shapefile

        :param IN_output_file: output file full path
        :type IN_output_file: string
        """
        
        if self.nb_pix_river == 0:
            my_api.printInfo("[l2_hr_pixc_vec_river] [write_file_asShp] == write_file_asShp ==")
            my_api.printInfo("[l2_hr_pixc_vec_river] [write_file_asShp] NO river pixels => NO PIXCVecRiver shapefile generated")
            
        else:
            
            my_api.printInfo("[l2_hr_pixc_vec_river] [write_file_asShp] == write_file_asShp : %s ==" % IN_output_file)
            
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
            outLayer = outDataSource.CreateLayer(str(os.path.basename(IN_output_file).split(".")[0]), srs, geom_type=ogr.wkbPoint)
            # 1.4 - Creation des attributs
            outLayer.CreateField(ogr.FieldDefn(str('az_index'), ogr.OFTInteger))
            outLayer.CreateField(ogr.FieldDefn(str('r_index'), ogr.OFTInteger))
            tmpField = ogr.FieldDefn(str('lat2'), ogr.OFTReal)
            tmpField.SetWidth(15)
            tmpField.SetPrecision(6)
            outLayer.CreateField(tmpField)
            tmpField = ogr.FieldDefn(str('long2'), ogr.OFTReal)
            tmpField.SetWidth(15)
            tmpField.SetPrecision(6)
            outLayer.CreateField(tmpField)
            tmpField = ogr.FieldDefn(str('wse2'), ogr.OFTReal)
            tmpField.SetWidth(15)
            tmpField.SetPrecision(6)
            outLayer.CreateField(tmpField)
            outLayer.CreateField(ogr.FieldDefn(str('node_id'), ogr.OFTString))
            outLayerDefn = outLayer.GetLayerDefn()
            
            # 2 - On traite point par point
            for indp in self.pixc_river_idx:
                # 2.1 - On cree l'objet dans le format de la couche de sortie
                outFeature = ogr.Feature(outLayerDefn)
                # 2.2 - On lui assigne le point
                point = ogr.Geometry(ogr.wkbPoint)
                point.AddPoint(self.longitude_vectorproc[indp], self.latitude_vectorproc[indp])
                outFeature.SetGeometry(point)
                # 2.3 - On lui assigne les attributs
                outFeature.SetField(str('az_index'), int(self.azimuth_index[indp]))
                outFeature.SetField(str('r_index'), int(self.range_index[indp]))
                outFeature.SetField(str('lat2'), float(self.latitude_vectorproc[indp]))
                outFeature.SetField(str('long2'), float(self.longitude_vectorproc[indp]))
                outFeature.SetField(str('wse2'), float(self.height_vectorproc[indp]))
                outFeature.SetField(str('node_id'), str(self.river_tag[indp]))
                # 2.4 - On ajoute l'objet dans la couche de sortie
                outLayer.CreateFeature(outFeature)
            
            # 3 - Destroy the data sources to free resources
            outDataSource.Destroy()
    
    #----------------------------------
    
    def set_vectorproc(self, IN_latitude_vectorproc, IN_longitude_vectorproc, IN_height_vectorproc):
        """
        Set latitude / longitude / height vectorproc values
        
        :param IN_latitude_vectorproc: improved latitude values (must be same size and same order as self.azimuth_index et al.)
        :type IN_latitude_vectorproc: 1D-array of float
        :param IN_longitude_vectorproc: improved longitude values (must be same size and same order as self.azimuth_index et al.)
        :type IN_longitude_vectorproc: 1D-array of float
        :param IN_height_vectorproc: improved height values (must be same size and same order as self.azimuth_index et al.)
        :type IN_height_vectorproc: 1D-array of float
        """
        
        self.latitude_vectorproc = IN_latitude_vectorproc
        self.longitude_vectorproc = IN_longitude_vectorproc
        self.height_vectorproc = IN_height_vectorproc

    def set_river_lake_tag(self, IN_vFlag):
        """
        Set river tag value
        
        :param IN_vFlag: water body values (0=lake; 1=river) (must be same size and same order as self.azimuth_index et al.)
        :type IN_vFlag: 1D-array of int
        """
        
        # Init as a 1D-array of string
        self.river_tag = np.empty(len(IN_vFlag), dtype=object)
        self.river_tag[:] = ""
        
        # Indices of river pixels
        TMP_river_idx = np.where(IN_vFlag == 1)[0]
        
        # Set tag of river pixels to a specific prefix
        self.river_tag[TMP_river_idx] = "1_"
    
    #----------------------------------
        
    def computeDate(self, IN_sec_from_start):
        """
        Compute date
        
        :param IN_sec_from_start: number of seconds from mission start time
        :type IN_sec_from_start: int
        
        :return: out_date_time = date in UTC
        :rtype: out_date_time = string YYYY-MM-DD HH:MM:SS.SSSSSSZ
        """
        
        # Computation
        tmp_time_split = self.mission_start_time.split("-")
        date_in_sec = datetime(int(tmp_time_split[0]), int(tmp_time_split[1]), int(tmp_time_split[2])) + timedelta(seconds=IN_sec_from_start)
        
        # Format
        out_date_time = "%sZ" % datetime.strftime(date_in_sec, '%Y-%m-%d %H:%M:%S.%f')
        return out_date_time
        
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

