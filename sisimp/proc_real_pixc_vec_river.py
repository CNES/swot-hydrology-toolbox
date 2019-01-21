# -*- coding: utf8 -*-
"""
.. module proc_real_pixc_vec_river.py
    :synopsis: Deal with data file complementary to pixel cloud files (L2_HR_PIXCVecRiver)
    Created on 08/29/2017
    2018/11/30 (D. Desroches, V. Poughon, C. Pottier - CNES): change variables names wrt to new PIXCVecRiver naming convention
                                                                    + write only river pixels (if empty, produce NetCDF file with record = 0)
                                                                    + add pixc_index field = indices of river pixels within the L2_HR_PIXC product

.. module author: Claire POTTIER - CNES DSO/SI/TR

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import os
from osgeo import ogr, osr

import lib.my_api as my_api
import lib.my_netcdf_file as my_nc


def fill_vector_param(IN_variable, IN_variable_name, IN_ref_size, IN_data_param):
    """
    Fill variable data field
    
    :param IN_variable: input data
    :type IN_variable: array
    :param IN_variable_name: name for the variable in the NetCDF file
    :type IN_variable_name: string
    :param IN_ref_size: number of elements associated to the variable in the NetCDF file
    :type IN_ref_size: int
    :param IN_data_param: pointer to the NetCDF writer (=NcWrite(OUT_file))
    :type IN_data_param: -
    """
    
    if IN_variable is not None:
        xsize = len(IN_variable)
        if xsize != IN_ref_size:
            exc = '[proc_real_pixc_vec_river/fill_vector_param] ERROR = There is a problem with the size of ' + IN_variable_name
            exit(exc)
        else:
            IN_data_param.fill_variable(IN_variable_name, IN_variable)


#----------------------------------------


class l2_hr_pixc_vec_river(object):

    def __init__(self, IN_azimuth_index, IN_range_index, \
                 IN_mission_start_time, IN_cycle_duration, IN_cycle_num, IN_pass_num, IN_tile_ref, IN_nb_pix_range, IN_nb_pix_azimuth):
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
        self.nb_pix_range = IN_nb_pix_range
        self.nb_pix_azimuth = IN_nb_pix_azimuth
        
        self.nb_pix_river = 0  # Number of PixC associated to river
        self.pixc_river_idx = None  # Indices of pixels in the L2_HR_PIXC product associated to a river

    def write_file(self, IN_output_file, noval, compress=False):
        """
        Write the pixel cloud vector attribute for river product (L2_HR_PIXCVecRiver product)
        NB: if there is no river pixels, write an empty file (with record = 0)

        :param IN_output_file: output full path
        :type IN_output_file: string
        :param noval: No data value
        :type noval: float
        :param compress: parameter the define to compress or not the file
        :type compress: boolean
        """ 
        my_api.printInfo("[l2_hr_pixc_vec_river] == write_file : %s ==" % IN_output_file)
        
        # 0 - Init noval if necessary
        if noval is None:
            noval = -999000000.
    
        # 1 - Open NetCDF file in writing mode
        data = my_nc.myNcWriter(IN_output_file)
        
        # 2 - Get river pixels indices
        self.pixc_river_idx = np.where(self.river_tag != "")[0]  # Indices of river pixels in the L2_HR_PIXC product
        self.nb_pix_river = len(self.pixc_river_idx)  # Number of river pixels
        
        # 3 - Write file depending on the number of river pixels
        if self.nb_pix_river == 0:
            
            my_api.printInfo("NO river pixels => empty PIXCVecRiver file generated")
            data.add_dimension('record', 0)
            
        else:

            data.add_dimension('record', self.nb_pix_river)
    
            data.add_variable('pixc_index', np.int32, 'record', np.int(noval), compress)
            fill_vector_param(self.pixc_river_idx, 'pixc_index', self.nb_pix_river, data)
            data.add_variable('azimuth_index', np.int32, 'record', np.int(noval), compress)
            fill_vector_param(self.azimuth_index[self.pixc_river_idx], 'azimuth_index', self.nb_pix_river, data)
            data.add_variable('range_index', np.int32, 'record', np.int(noval), compress)
            fill_vector_param(self.range_index[self.pixc_river_idx], 'range_index', self.nb_pix_river, data)
        
            data.add_variable('latitude_vectorproc', np.double, 'record', noval, compress)
            data.add_variable_attribute('latitude_vectorproc', 'units', 'degrees_north')
            fill_vector_param(self.latitude_vectorproc[self.pixc_river_idx], 'latitude_vectorproc', self.nb_pix_river, data)
            data.add_variable('longitude_vectorproc', np.double, 'record', noval, compress)
            data.add_variable_attribute('longitude_vectorproc', 'units', 'degrees_east')
            fill_vector_param(self.longitude_vectorproc[self.pixc_river_idx], 'longitude_vectorproc', self.nb_pix_river, data)
            data.add_variable('height_vectorproc', np.float32, 'record', np.float(noval), compress)
            data.add_variable_attribute('height_vectorproc', 'units', 'm')
            fill_vector_param(self.height_vectorproc[self.pixc_river_idx], 'height_vectorproc', self.nb_pix_river, data)
        
            data.add_variable('river_tag', str, 'record', "", compress)
            fill_vector_param(self.river_tag[self.pixc_river_idx], 'river_tag', self.nb_pix_river, data)
        
        # Write global attributes even if empty
        data.add_global_attribute('description', 'L2_HR_PIXCVecRiver product obtained by CNES/LEGOS Large Scale Simulator')
        data.add_global_attribute('mission_start_time', self.mission_start_time)
        data.add_global_attribute('repeat_cycle_period', self.cycle_duration)
        data.add_global_attribute('cycle_number', self.cycle_num)
        data.add_global_attribute('pass_number', np.int(self.pass_num))
        data.add_global_attribute('tile_name', self.tile_ref)
        data.add_global_attribute('interferogram_size_range', self.nb_pix_range)    
        data.add_global_attribute('interferogram_size_azimuth', self.nb_pix_azimuth)   
        
        # 4 - Close output file
        data.close()

    def write_file_asShp(self, IN_output_file):
        """
        Write the pixel cloud vector attribute for rivers file as a shapefile

        :param IN_output_file: output file full path
        :type IN_output_file: string
        """
        
        if self.nb_pix_river == 0:
            my_api.printInfo("[l2_hr_pixc_vec_river] == write_file_asShp ==") 
            my_api.printInfo("NO river pixels => NO PIXCVecRiver shapefile generated")
            
        else:
            
            my_api.printInfo("[l2_hr_pixc_vec_river] == write_file_asShp : %s ==" % IN_output_file)
            
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
            outLayer.CreateField(ogr.FieldDefn(str('AZ_INDEX'), ogr.OFTInteger))
            outLayer.CreateField(ogr.FieldDefn(str('R_INDEX'), ogr.OFTInteger))
            tmpField = ogr.FieldDefn(str('LAT2'), ogr.OFTReal)
            tmpField.SetWidth(10)
            tmpField.SetPrecision(6)
            outLayer.CreateField(tmpField)
            tmpField = ogr.FieldDefn(str('LONG2'), ogr.OFTReal)
            tmpField.SetWidth(10)
            tmpField.SetPrecision(6)
            outLayer.CreateField(tmpField)
            tmpField = ogr.FieldDefn(str('HEIGHT2'), ogr.OFTReal)
            tmpField.SetWidth(10)
            tmpField.SetPrecision(6)
            outLayer.CreateField(tmpField)
            outLayer.CreateField(ogr.FieldDefn(str('TAG'), ogr.OFTString))
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
                outFeature.SetField(str('AZ_INDEX'), int(self.azimuth_index[indp]))
                outFeature.SetField(str('R_INDEX'), int(self.range_index[indp]))
                outFeature.SetField(str('LAT2'), float(self.latitude_vectorproc[indp]))
                outFeature.SetField(str('LONG2'), float(self.longitude_vectorproc[indp]))
                outFeature.SetField(str('HEIGHT2'), float(self.height_vectorproc[indp]))
                outFeature.SetField(str('TAG'), str(self.river_tag[indp]))
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
                

#######################################


if __name__ == '__main__':
    
    noval = -999900000

    record = 10
    azimuth_index = np.ones(record)+50.
    range_index = np.ones(record)+50.
    latitude = np.arange(record)+10.
    longitude = np.arange(record)+50.
    height = np.ones(record)+50.
    
    cycle_num = 1
    pass_num = 10
    tile_ref = "45N-L"
    nb_pix_range = 200
    nb_pix_azimuth = 20
    
    # Test avec vrai PixC
    pixc_file = os.path.join(os.getcwd(), "SWOT_L2_HR_PIXC_001_010_45N-L_main.nc")
    coord = my_nc.myNcReader(pixc_file)
    # Init PIXCVecRiver product
    my_pixc_vec = l2_hr_pixc_vec_river(coord.get_variable("azimuth_index"), \
                                       coord.get_variable("range_index"), \
                                       coord.get_global_attribute("cycle_number"), \
                                       coord.get_global_attribute("pass_number"), \
                                       coord.get_global_attribute("tile_ref"), \
                                       coord.get_global_attribute("nb_pix_range"), \
                                       coord.get_global_attribute("nb_pix_azimuth"))
    
