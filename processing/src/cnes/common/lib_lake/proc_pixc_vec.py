# -*- coding: utf8 -*-
"""
.. module:: proc_pixc_vec.py
    :synopsis: Deals with SWOT pixel cloud complementary file
    Created on 09/15/2017
    
..todo:: revoir la sélection des pixels rivières

.. moduleauthor: Claire POTTIER - CNES DSO/SI/TR

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from __future__ import absolute_import, division, print_function, unicode_literals 

import datetime
import numpy as np
import os
from osgeo import ogr, osr
import logging

import cnes.modules.geoloc.lib.geoloc as my_geoloc

import cnes.common.lib.my_api as my_api
import cnes.common.lib.my_netcdf_file as my_nc
import cnes.common.lib.my_variables as my_var
import cnes.common.lib_lake.locnes_variables as my_var2


class PixelCloudVec(object):

    def __init__(self, IN_productType, IN_pixc_vec):
        """
        Constructor: init with data retrieved from pixel cloud complementary file after river processing
        
        :param IN_productType: type of product among "SP"=LakeSP and "TILE"=LakeTile
        :type IN_productType: string
        :param IN_pixc_vec: full path of L2_HR_PIXCVec / pixel cloud complementary file (from PGE_RiverTile or PGE_LakeTile)
        :type IN_pixc_vec: string
        
        Variables of the object:
        - From L2_HR_PIXCVecRiver
            nb_pix_range / int: number of pixels in range dimension (= global attribute named nr_pixels in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
            nb_pix_azimuth / int: number of pixels in azimuth dimension (= global attribute named nr_lines in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
            cycle_num / int: cycle number (= global attribute named cycle_number in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
            pass_num / int: pass number (= global attribute named pass_number in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
            tile_ref / int: tile reference (= global attribute named tile_ref in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
            river_idx / 1D-array of int: indices of river pixels within PIXC arrays (= variable named pixc_index in L2_HR_PIXCVecRiver only)
            range_idx / 1D-array of int: range indices of water pixels (= variable named range_index in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
            azimuth_idx / 1D-array of int: azimuth indices of water pixels (= variable named azimuth_index in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
            longitude_vectorproc / 1D-array of float: improved longitude of water pixels (= variable named longitude_vectorproc in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
            latitude_vectorproc / 1D-array of float: improved latitude of water pixels (= variable named latitude_vectorproc in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
            height_vectorproc / 1D-array of float: improved height of water pixels (= variable named height_vectorproc in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
            tag / 1D-array of float: tag associated to river and lake databases (= variable named river_tag in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
        - From L2_HR_PIXC:
            continent / string: continent covered by the tile (if global var CONTINENT_FILE exists)
        - From processing:
            nb_water_pix / int: number of water pixels
            reject_idx / 1D-array of int: indices of pixels that are river only, ie not reservoirs or dams
            nb_river_pix / int: number of river pixels
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("[PixelCloudVec] == INIT ==")
        if IN_productType == "TILE":
            logger.info("[PixelCloudVec] L2_HR_PIXCVec file = %s" % IN_pixc_vec)
        elif IN_productType == "SP":
            logger.info("[PixelCloudVec] LakeTile_pixcvec file = %s" % IN_pixc_vec)
        else:
            logger.info("[PixelCloudVec] L2_HR_PIXCVecRiver file = %s" % IN_pixc_vec)
            
        # 0 - Init variables
        self.range_idx = None
        self.azimuth_idx = None
        self.nb_water_pix = 0
        self.longitude_vectorproc = None
        self.latitude_vectorproc = None
        self.height_vectorproc = None
        self.tag = None
        self.reservoirs_idx = None
        if IN_productType == "TILE":
            self.nb_river_pix = 0  # Number of river pixels (used for TILE processing)
            self.river_idx = None  # Indices of pixels processed by RiverTile (used in TILE processing)
            self.reject_idx = None  # Indices of river pixels (not reservoirs)
        
        # 1 - Retrieve needed information from complementary pixel cloud file
        pixc_vec_river = my_nc.myNcReader(IN_pixc_vec)
        # 1.1 - Number of records
        try:
            self.nb_water_pix = pixc_vec_river.getDimValue("points")  
        except:
            print("SWOT_hydrology_toolbox version, use record instead of points, to be corrected")
            self.nb_water_pix = pixc_vec_river.getDimValue("record")

        # 1.2 - Global attributes
        self.cycle_num = pixc_vec_river.getAttValue("cycle_number")
        self.pass_num = pixc_vec_river.getAttValue("pass_number")
        self.tile_ref = pixc_vec_river.getAttValue("tile_name")
        
        shape = pixc_vec_river.getAttValue("interferogram_shape")
        self.nb_pix_azimuth = int(shape.split(",")[0])
        self.nb_pix_range = int((shape.split(",")[1]).split("(")[0])
        
        # Retrieve others if there are river pixels
        if self.nb_water_pix != 0:
                
            # 1.3 - Range index
            self.range_idx = pixc_vec_river.getVarValue("range_index")
            # 1.4 - Longitude
            self.azimuth_idx = pixc_vec_river.getVarValue("azimuth_index")
            # 1.5 - Longitude
            self.longitude_vectorproc = pixc_vec_river.getVarValue("longitude_vectorproc")
            # 1.6 - Latitude
            self.latitude_vectorproc = pixc_vec_river.getVarValue("latitude_vectorproc")
            # 1.7 - Height
            self.height_vectorproc = pixc_vec_river.getVarValue("height_vectorproc")
            # 1.8 - Reference to a priori
            if IN_productType == "TILE":
                TMP_tag = [str(i) for i in pixc_vec_river.getVarValue("river_tag")]
                self.tag = np.array(TMP_tag, dtype=object)
            elif IN_productType == "SP":
                TMP_tag = [str(i) for i in pixc_vec_river.getVarValue("river_lake_other_tag")]
                self.tag = np.array(TMP_tag, dtype=object)
            else:
                self.tag = np.empty(self.nb_water_pix, dtype=object)
                self.tag[:] = ""
            
            if IN_productType == "TILE":
                
                # Number of pixels of PixC already processed by PGE_RiverTile
                self.nb_river_pix = self.nb_water_pix
                
                # Indices of pixels of PixC already processed by PGE_RiverTile
                self.river_idx = pixc_vec_river.getVarValue("pixc_index") 
                
                # Indices of pixels of PixC not to remove from LakeTile processing = river only pixels (ie reservoirs kept)
                try:
                    TMP_lake_flag = pixc_vec_river.getVarValue("lake_flag")  # Copy of lakeflag variable from river DB: river (lakeflag=0), lake/reservoir (lakeflag=1), tidal river (lakeflag=2), or canal (lakeflag=3)
                except: 
                    TMP_lake_flag = np.zeros(self.nb_river_pix)  # If lake_flag not in PIXCVecRiver product, consider all river objects as river pixels (lakeflag=0)
                self.reject_idx = self.river_idx[np.where(TMP_lake_flag != 1)[0]]  # lakeFlag == 1 for lakes/reservoirs
                 
        # 2 - Close file
        pixc_vec_river.close()
    
    def reshape(self, IN_objPixc):
        """
        Reshape PIXCVecRiver arrays to new arrays of size of related PIXC arrays
        
        :param IN_objPixc: pixel cloud associated to current PIXCVecRiver object
        :type IN_objPixc: proc_pixc.PixelCloud
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("[PixelCloudVec] == reshape ==")
        
        # 1 - Azimuth and range indices arrays, number of water pixels
        self.range_idx = IN_objPixc.origin_range_idx
        self.azimuth_idx = IN_objPixc.origin_azimuth_idx
        self.nb_water_pix = len(self.azimuth_idx)
        
        # 2 - Init the other arrays
        TMP_longitude_vectorproc = np.zeros(self.nb_water_pix)
        TMP_latitude_vectorproc = np.zeros(self.nb_water_pix)
        TMP_height_vectorproc = np.zeros(self.nb_water_pix)
        TMP_tag = np.empty(self.nb_water_pix, dtype=object)
        TMP_tag[:] = ""  # For isolated pixels, not processed by neither RiverTile nor LakeTile, tag=""
        
        # 3 - Include river pixels info if there is
        if self.nb_river_pix == 0:
            logger.info("[PixelCloudVec] No pixel associated to river")
        else:
            
            logger.info("[PixelCloudVec] %d pixels associated to rivers" % self.nb_river_pix)
            TMP_longitude_vectorproc[self.river_idx] = self.longitude_vectorproc
            TMP_latitude_vectorproc[self.river_idx] = self.latitude_vectorproc
            TMP_height_vectorproc[self.river_idx] = self.height_vectorproc
            TMP_tag[self.river_idx] = self.tag
            
        # 4 - Save arrays
        self.longitude_vectorproc = TMP_longitude_vectorproc
        self.latitude_vectorproc = TMP_latitude_vectorproc
        self.height_vectorproc = TMP_height_vectorproc
        self.tag = TMP_tag
        
        # 5 - Save continent var
        self.continent = IN_objPixc.continent
    
    def setContinent(self, IN_continent):
        """
        Setter for continent variable
        
        :param IN_continent: continent covered by the tile (if global var CONTINENT_FILE exists)
        :type IN_continent: string
        """
        self.continent = IN_continent
    
    # ----------------------------------------
    
    def getTileInfo(self):
        """
        Getter of cycle_num, pass_num and tile_ref
        """
        return self.cycle_num, self.pass_num, self.tile_ref
        
    # ----------------------------------------
    
    def computeRiverPix(self):
        """
        Compute indices of pixels related to river
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("[PixelCloudVec] == computeRiverPix ==")
        
        # 1 - Get river pixels, ie for which tag != 0
        self.river_idx = np.where(self.tag != "")[0]
        
        # 2 - Deduce river pixels number
        self.nb_river_pix = self.river_idx.size
        if self.nb_river_pix == 0:
            logger.info("[PixelCloudVec] No pixel associated to river")
        else:
            logger.info("[PixelCloudVec] %d pixels associated to rivers" % self.river_idx.size)
        
    # ----------------------------------------
    
    def write_file(self, IN_filename, noval, compress=False):
        """
        Write the pixel cloud vector attribute product (L2_HR_PIXCVec product)

        :param IN_filename: full path of the output file
        :type IN_filename: string
        :param noval: No data value
        :type noval: float
        :param compress: parameter the define to compress or not the file
        :type compress: boolean
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("[PixelCloudVec] == write_file ==")
        logger.info("[PixelCloudVec] Output L2_HR_PIXCVec NetCDF file = %s" % IN_filename)
    
        # Open file in writing mode
        data = my_nc.myNcWriter(IN_filename)
        
        if noval is None:
            noval = -999000000.

        data.add_dimension('record', self.nb_water_pix)
    
        data.add_variable('azimuth_index', np.int32, 'record', IN_fill_value=np.int(noval), IN_compress=compress)
        data.fill_variable('azimuth_index', self.azimuth_idx)
        data.add_variable('range_index', np.int32, 'record', IN_fill_value=np.int(noval), IN_compress=compress)
        data.fill_variable('range_index', self.range_idx)
        
        data.add_variable('latitude_vectorproc', np.double, 'record', IN_fill_value=noval, IN_compress=compress)
        data.add_variable_attribute('latitude_vectorproc', 'units', 'degrees_north')
        data.fill_variable('latitude_vectorproc', self.latitude_vectorproc)
        data.add_variable('longitude_vectorproc', np.double, 'record', IN_fill_value=noval, IN_compress=compress)
        data.add_variable_attribute('longitude_vectorproc', 'units', 'degrees_east')
        data.fill_variable('longitude_vectorproc', self.longitude_vectorproc)
        data.add_variable('height_vectorproc', np.float64, 'record', IN_fill_value=np.float(noval), IN_compress=compress)
        data.add_variable_attribute('height_vectorproc', 'units', 'm')
        data.fill_variable('height_vectorproc', self.height_vectorproc)
        
        #~ data.add_variable('river_lake_other_tag', str, 'record', IN_fill_value="", IN_compress=compress)
        #~ print("self.tag=", self.tag)
        #~ data.fill_variable('river_lake_other_tag', self.tag)
        
        data.add_global_attribute('producer', my_var2.PRODUCER)
        data.add_global_attribute('creation_date', str(datetime.datetime.now()))
        data.add_global_attribute('cycle_number', self.cycle_num)
        data.add_global_attribute('pass_number', self.pass_num)
        data.add_global_attribute('tile_ref', self.tile_ref)
        if my_var2.CONTINENT_FILE is not None:
            data.add_global_attribute('continent', self.continent)
        data.add_global_attribute('nr_pixels', self.nb_pix_range)
        data.add_global_attribute('nr_lines', self.nb_pix_azimuth)
        data.add_global_attribute('lake_db', my_var2.LAKE_DB)
        if my_var2.CONTINENT_FILE is not None:
            data.add_global_attribute('continent_file', my_var2.CONTINENT_FILE)
        data.add_global_attribute('flag_water', my_var2.FLAG_WATER)
        data.add_global_attribute('flag_dark', my_var2.FLAG_DARK)
        data.add_global_attribute('flag_layover', my_var2.FLAG_LAYOVER)
        data.add_global_attribute('min_size', my_var2.MIN_SIZE)
        TMP_geoloc = 0
        if my_var2.IMP_GEOLOC:
            TMP_geoloc = 1
        data.add_global_attribute('imp_geoloc', TMP_geoloc)
        data.add_global_attribute('hull_method', my_var2.HULL_METHOD)
        data.add_global_attribute('std_height_max', my_var2.STD_HEIGHT_MAX)
        data.add_global_attribute('biglake_model', my_var2.BIGLAKE_MODEL)
        data.add_global_attribute('biglake_min_size', my_var2.BIGLAKE_MIN_SIZE)
        data.add_global_attribute('biglake_grid_spacing', my_var2.BIGLAKE_GRID_SPACING)
        data.add_global_attribute('biglake_grid_res', my_var2.BIGLAKE_GRID_RES)
        
        # Close file
        data.close()

    def write_file_asShp(self, IN_filename, IN_classif=None):
        """
        Write the pixel cloud vector product as a shapefile

        :param IN_filename: full path of the output file
        :type IN_filename: string
        :param IN_classif: classification flag for each pixel (retrieved from the PIXC data)
        :type IN_classif: 1D array of int
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("[PixelCloudVec] == write_file_asShp ==")
        logger.info("[PixelCloudVec] Output L2_HR_PIXCVec shapefile = %s" % IN_filename)
        
        # 1 - Initialisation du fichier de sortie
        # 1.1 - Driver
        shpDriver = ogr.GetDriverByName(str("ESRI Shapefile"))
        # 1.2 - Creation du fichier
        if os.path.exists(IN_filename):
            shpDriver.DeleteDataSource(IN_filename)
        outDataSource = shpDriver.CreateDataSource(IN_filename)
        # 1.3 - Creation de la couche
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)  # WGS84
        outLayer = outDataSource.CreateLayer(str(str(os.path.basename(IN_filename)).replace('.shp', '')), srs, geom_type=ogr.wkbPoint)
        # 1.4 - Creation des attributs
        outLayer.CreateField(ogr.FieldDefn(str('AZ_INDEX'), ogr.OFTInteger))
        outLayer.CreateField(ogr.FieldDefn(str('R_INDEX'), ogr.OFTInteger))
        if IN_classif is not None:
            outLayer.CreateField(ogr.FieldDefn(str('CLASSIF'), ogr.OFTInteger))
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
        tmpField.SetPrecision(4)
        outLayer.CreateField(tmpField)
        outLayer.CreateField(ogr.FieldDefn(str('TAG'), ogr.OFTString))
        outLayerDefn = outLayer.GetLayerDefn()
        
        # 2 - On traite point par point
        for indp in range(self.nb_water_pix):
            # 2.1 - On cree l'objet dans le format de la couche de sortie
            outFeature = ogr.Feature(outLayerDefn)
            # 2.2 - On lui assigne le point
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(self.longitude_vectorproc[indp], self.latitude_vectorproc[indp])
            outFeature.SetGeometry(point)
            # 2.3 - On lui assigne les attributs
            outFeature.SetField(str('AZ_INDEX'), int(self.azimuth_idx[indp]))
            outFeature.SetField(str('R_INDEX'), int(self.range_idx[indp]))
            if IN_classif is not None:
                outFeature.SetField(str('CLASSIF'), int(IN_classif[indp]))
            outFeature.SetField(str('LAT2'), float(self.latitude_vectorproc[indp]))
            outFeature.SetField(str('LONG2'), float(self.longitude_vectorproc[indp]))
            outFeature.SetField(str('HEIGHT2'), float(self.height_vectorproc[indp]))
            outFeature.SetField(str('TAG'), str(self.tag[indp]))
            # 2.4 - On ajoute l'objet dans la couche de sortie
            outLayer.CreateFeature(outFeature)
            
        # 3 - Destroy the data sources to free resources
        outDataSource.Destroy()
        
        
#######################################
    
def computeImpGeoloc(IN_productType, IN_objPixc, IN_indices, IN_height):
    """
    Refines geolocation for IN_objPixc pixels corresponding to indices IN_idx (in IN_objPix)
        
    :param IN_productType: type of product among "SP"=LakeSP and "TILE"=LakeTile
    :type IN_productType: string
    :param IN_objPixc: pixel cloud from which to compute improved geolocation
    :type IN_objPixc: proc_pixc.PixelCloud or proc_pixc_sp.PixelCloudSP object
    :param IN_indices: indices of pixels related to the same object
    :type IN_indices: 1D-array of int
    :param IN_height: new height of each point used for improved geoloc
    :type IN_height: 1D-array of float
        
    :return OUT_lon_corr: improved longitude of pixels of IN_indices
    :rtype: 1D-array of float
    :return OUT_lat_corr: improved latitude of pixels of IN_indices
    :rtype: 1D-array of float
    :return OUT_lat_corr: improved latitude of pixels of IN_indices
    :rtype: 1D-array of float
    """
    logger = logging.getLogger("proc_pixc_vec")
    logger.debug("[PixelCloudVec] == computeImpGeoloc ==")
    nb_pix = IN_indices.size
    logger.debug("> %d PixC to deal with" % nb_pix)
        
    # 1 - Prepare data for Damien's algo
    # 1.1 - Convert geodetic coordinates (lat, lon, height) to cartesian coordinates (x, y, z)
    x, y, z = my_geoloc.convert_llh2ecef(IN_objPixc.latitude[IN_indices], IN_objPixc.longitude[IN_indices], IN_objPixc.height[IN_indices], my_var.GEN_RAD_EARTH_EQ, my_var.GEN_RAD_EARTH_POLE)
    # 1.2 - Get position of associated along-track pixels (in cartesian coordinates)
    nadir_x = IN_objPixc.nadir_x[IN_indices]
    nadir_y = IN_objPixc.nadir_y[IN_indices]
    nadir_z = IN_objPixc.nadir_z[IN_indices]
    # 1.3 - Get velocity of associated along-track pixels (in cartesian coordinates)
    nadir_vx = IN_objPixc.nadir_vx[IN_indices]
    nadir_vy = IN_objPixc.nadir_vy[IN_indices]
    nadir_vz = IN_objPixc.nadir_vz[IN_indices]
    # 1.4 - Get distance from satellite to target point
    ri = IN_objPixc.near_range + IN_objPixc.range_idx[IN_indices] * my_var.GEN_RANGE_SPACING
    
    # 2 - Use Damien's algo
    # 2.1 - Init output vectors
    OUT_lat_corr = np.zeros(nb_pix)  # Improved latitudes
    OUT_lon_corr = np.zeros(nb_pix)  # Improved longitudes
    OUT_height_corr = np.zeros(nb_pix)  # Improved heights
    # 2.2 - Improve geolocation
    p_final, p_final_llh, h_mu, (iter_grad,nfev_minimize_scalar) = my_geoloc.pointcloud_height_geoloc_vect(np.transpose(np.array([x, y, z])), IN_objPixc.height[IN_indices],
                                                                                                           np.transpose(np.array([nadir_x, nadir_y, nadir_z])),
                                                                                                           np.transpose(np.array([nadir_vx, nadir_vy, nadir_vz])),
                                                                                                           ri, IN_height, 
                                                                                                           recompute_Doppler=True, recompute_R=True, verbose = False, 
                                                                                                           max_iter_grad=1, height_goal = 1.e-3, safe_flag=True)
    # 2.3 - Save output variables
    OUT_lat_corr, OUT_lon_corr, OUT_height_corr = np.transpose(p_final_llh)
    
    return OUT_lon_corr, OUT_lat_corr, OUT_height_corr
