# -*- coding: utf8 -*-
"""
.. module:: proc_pixc_vec.py
    :synopsis: Deals with SWOT pixel cloud complementary file
    Created on 09/15/2017

.. moduleauthor: Claire POTTIER - CNES DSO/SI/TR

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""

import numpy as np
import os
from osgeo import ogr, osr
import logging

import cnes.modules.geoloc.lib.geoloc as my_geoloc

import cnes.common.lib.my_netcdf_file as my_nc
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_variables as my_var
import cnes.common.lib_lake.locnes_products_netcdf as nc_file


class PixelCloudVec(object):

    def __init__(self, in_productType, in_pixcvec_file=None):
        """
        Constructor: init variables and set them with data retrieved from pixel cloud complementary file if asked
        
        :param in_productType: type of product among "SP"=LakeSP and "TILE"=LakeTile
        :type in_productType: string
        :param in_pixcvec_file: full path of pixel cloud complementary file 
                                    (L2_HR_PIXCVecRiver file if from PGE_RiverTile 
                                    or LakeTile_piexcvec if from PGE_LakeTile)
        :type in_pixcvec_file: string
        
        Variables of the object:
            
        - From L2_HR_PIXCVecRiver:
            river_index / 1D-array of int: indices of river pixels within PIXC arrays (= variable named pixc_index in L2_HR_PIXCVecRiver only)
        
        - From both L2_HR_PIXCVecRiver and LakeTile_pixcvec:
            range_index / 1D-array of int: range indices of water pixels (= variable named range_index in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
            azimuth_index / 1D-array of int: azimuth indices of water pixels (= variable named azimuth_index in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
            longitude_vectorproc / 1D-array of float: improved longitude of water pixels (= variable named longitude_vectorproc in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
            latitude_vectorproc / 1D-array of float: improved latitude of water pixels (= variable named latitude_vectorproc in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
            height_vectorproc / 1D-array of float: improved height of water pixels (= variable named height_vectorproc in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
            river_reach_tag / 1D-array of float: tag associated to river reach database (= variable named reach_index in L2_HR_PIXCVecRiver and river_reach_tag in LakeTile_pixcvec)
            river_node_tag / 1D-array of float: tag associated to river node database (= variable named node_index in L2_HR_PIXCVecRiver and river_node_tag in LakeTile_pixcvec)
            pixcvec_metadata / dict: dictionary of PIXCVec file metadata
        
        - From L2_HR_PIXC:
            continent / string: continent covered by the tile (if global var CONTINENT_FILE exists)
        
        - From processing:
            nb_water_pix / int: number of water pixels
            reject_index / 1D-array of int: indices of pixels that are river only, ie not reservoirs or dams
            nb_river_pix / int: number of river pixels
            greenwich_idx / 1D-array of int: indices of pixels related to lakes crossing Greenwich meridian
            lake_tag / 1D-array of str: tag associated to lake database (= variable named lake_tag in LakeTile_pixcvec)
            other_tag / 1D-array of str: tag associated to unknown object (= variable named other_tag in LakeTile_pixcvec)
            prior_ice_flag / 1D-array of int: TBD
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # Init type
        if in_productType in ("TILE", "SP"):
            self.product_type = in_productType
        else:
            logger.debug("Product type %s unknown, set to TILE" % in_productType)
            self.product_type = "TILE"
            
        # Init dimension
        self.nb_water_pix = 0
        
        # Init PIXCVec variables
        self.azimuth_index = None
        self.range_index = None
        self.longitude_vectorproc = None
        self.latitude_vectorproc = None
        self.height_vectorproc = None
        self.river_reach_tag = None
        self.river_node_tag = None
        self.lake_tag = None
        self.other_tag = None
        self.prior_ice_flag = None
        
        # Fill variables if filename available
        if in_pixcvec_file is not None:
            self.set_from_pixcvec_file(in_pixcvec_file)
            
        # Init dictionary of PIXCVec metadata
        self.pixcvec_metadata = {}
        self.pixcvec_metadata["cycle_number"] = -9999
        self.pixcvec_metadata["pass_number"] = -9999
        self.pixcvec_metadata["tile_number"] = -9999
        self.pixcvec_metadata["swath_side"] = ""
        self.pixcvec_metadata["tile_name"] = ""
        self.pixcvec_metadata["start_time"] = ""
        self.pixcvec_metadata["stop_time"] = ""
        self.pixcvec_metadata["inner_first_latitude"] = -9999.0
        self.pixcvec_metadata["inner_first_longitude"] = -9999.0
        self.pixcvec_metadata["inner_last_latitude"] = -9999.0
        self.pixcvec_metadata["inner_last_longitude"] = -9999.0
        self.pixcvec_metadata["outer_first_latitude"] = -9999.0
        self.pixcvec_metadata["outer_first_longitude"] = -9999.0
        self.pixcvec_metadata["outer_last_latitude"] = -9999.0
        self.pixcvec_metadata["outer_last_longitude"] = -9999.0
        self.pixcvec_metadata["continent"] = ""
        self.pixcvec_metadata["ellipsoid_semi_major_axis"] = ""
        self.pixcvec_metadata["ellipsoid_flattening"] = ""
        self.pixcvec_metadata["xref_static_river_db_file"] = ""
        
        # Variables specific to processing
        self.continent = None
        # Specific to LakeTile processing
        if self.product_type == "TILE":
            self.nb_river_pix = 0  # Number of river pixels (used for TILE processing)
            self.river_index = None  # Indices of pixels processed by RiverTile (used in TILE processing)
            self.reject_index = None  # Indices of river pixels (not reservoirs)
            self.greenwich_idx = []  # Indices of pixels related to lakes crossing Greenwich meridian
            
    def set_from_pixcvec_file(self, in_pixcvec_file):
        """
        Set variables from PIXCVec file
        
        :param in_pixcvec_file: full path of pixel cloud complementary file 
                                    (L2_HR_PIXCVecRiver file if from PGE_RiverTile 
                                    or LakeTile_piexcvec if from PGE_LakeTile)
        :type in_pixcvec_file: string
        """
        logger = logging.getLogger(self.__class__.__name__)
        if self.product_type == "TILE":
            logger.info("L2_HR_PIXCVec file = %s", in_pixcvec_file)
        elif self.product_type == "SP":
            logger.info("LakeTile_pixcvec file = %s", in_pixcvec_file)
        else:
            logger.debug("Product type %s unknown, set to TILE" % self.product_type)
            self.product_type = "TILE"
        
        # 1 - Open file in reading mode
        pixcvec_reader = my_nc.myNcReader(in_pixcvec_file)
        
        # 2 - Retrieve the number of records
        try:
            self.nb_water_pix = pixcvec_reader.getDimValue("points")
        except:
            logger.info("Warning, using record instead of points for nb_water_pixels")
            self.nb_water_pix = pixcvec_reader.getDimValue("record")
            
        # 3 - Retrieve needed global attributes
        pixcvec_keys = pixcvec_reader.getListAtt()
        for key, value in self.pixcvec_metadata.items():
            if key in pixcvec_keys:
                self.pixcvec_metadata[key] = pixcvec_reader.getAttValue(key)
        
        # 4 - Retrieve variables if there are river pixels
        if self.nb_water_pix != 0:
                
            # 4.1 - Range index
            self.range_index = pixcvec_reader.getVarValue("range_index")
            # 4.2 - Azimuth index
            self.azimuth_index = pixcvec_reader.getVarValue("azimuth_index")
            # 4.3 - Longitude
            self.longitude_vectorproc = pixcvec_reader.getVarValue("longitude_vectorproc")
            # 4.4 - Latitude
            self.latitude_vectorproc = pixcvec_reader.getVarValue("latitude_vectorproc")
            # 4.5 - Height
            self.height_vectorproc = pixcvec_reader.getVarValue("height_vectorproc")
            # 4.6 - References to a priori
            if self.product_type == "TILE":
                self.river_reach_tag = pixcvec_reader.getVarValue("reach_index")
                self.river_node_tag = pixcvec_reader.getVarValue("node_index")
            elif self.product_type == "SP":
                self.river_reach_tag = pixcvec_reader.getVarValue("river_reach_tag")
                self.river_node_tag = pixcvec_reader.getVarValue("river_node_tag")
                tmp_tag = [str(i) for i in pixcvec_reader.getVarValue("lake_tag")]
                self.lake_tag = np.array(tmp_tag, dtype=object)
                tmp_tag = [str(i) for i in pixcvec_reader.getVarValue("other_tag")]
                self.other_tag = np.array(tmp_tag, dtype=object)
            # 4.7 - Prior ice flag
            try:
                self.prior_ice_flag = pixcvec_reader.getVarValue("prior_ice_flag")
            except:
                logger.debug("prior_ice_flag variable not in PIXCVec file => set to 1D-array of _FillValue")
                self.prior_ice_flag = np.zeros(self.nb_water_pix, dtype=np.uint8) + my_var.FV_NETCDF['uint8']
            
            # 5 - Reject river pixels from the list of pixel to process
            if self.product_type == "TILE":
                self.computeRiverPix(pixcvec_reader)
                 
        # 6 - Close file
        pixcvec_reader.close()
    
    def computeRiverPix(self, in_pixcvec_reader):
        """
        Compute indices of pixels already processed by RiverTile (except reservoirs)
        
        :param in_pixcvec_reader: reader of L2_HR_PIXCVecRiver file
        :type in_pixcvec_reader: my_netcdf_file.myNcReader
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
                
        # Indices of pixels of PixC already processed by PGE_RiverTile
        self.river_index = in_pixcvec_reader.getVarValue("pixc_index") 
        
        # Number of pixels of PixC already processed by PGE_RiverTile
        self.nb_river_pix = self.river_index.size
                
        # Indices of pixels of PixC not to remove from LakeTile processing = river only pixels (ie reservoirs kept)
        try:
            tmp_lake_flag = in_pixcvec_reader.getVarValue("lake_flag")  # Copy of lakeflag variable from river DB: river (lakeflag=0), lake/reservoir (lakeflag=1), tidal river (lakeflag=2), or canal (lakeflag=3)
        except: 
            logger.debug("lake_flag variable not in PIXCVecRiver product => consider all river objects as river pixels (lakeflag=0)")
            tmp_lake_flag = np.zeros(self.nb_river_pix)  # If lake_flag not in PIXCVecRiver product, consider all river objects as river pixels (lakeflag=0)
        self.reject_index = self.river_index[np.where(tmp_lake_flag != 1)[0]]  # lakeFlag == 1 for lakes/reservoirs
        
        # 2 - Deduce river pixels number
        if self.nb_river_pix == 0:
            logger.info("No pixel associated to river")
        else:
            logger.info("%d pixels associated to rivers", self.nb_river_pix)
            logger.info("%d pixels associated to reservoirs", self.nb_river_pix-self.reject_index.size)
        
    # ----------------------------------------
    
    def reshape(self, in_objPixc):
        """
        Reshape PIXCVecRiver arrays to new arrays of size of related PIXC arrays
        
        :param in_objPixc: pixel cloud associated to current PIXCVecRiver object
        :type in_objPixc: proc_pixc.PixelCloud
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Azimuth and range indices arrays, number of water pixels
        self.range_index = in_objPixc.origin_range_index
        self.azimuth_index = in_objPixc.origin_azimuth_index
        self.nb_water_pix = len(self.azimuth_index)
        
        # 2 - Init the other arrays
        tmp_longitude_vectorproc = np.zeros(self.nb_water_pix, dtype=np.float64) + my_var.FV_NETCDF['float64']
        tmp_latitude_vectorproc = np.zeros(self.nb_water_pix, dtype=np.float64) + my_var.FV_NETCDF['float64']
        tmp_height_vectorproc = np.zeros(self.nb_water_pix, dtype=np.float32) + my_var.FV_NETCDF['float32']
        tmp_river_reach_tag = np.zeros(self.nb_water_pix, dtype=np.int32) + my_var.FV_NETCDF['int32']
        tmp_river_node_tag = np.zeros(self.nb_water_pix, dtype=np.int32) + my_var.FV_NETCDF['int32']
        self.lake_tag = np.empty(self.nb_water_pix, dtype=object)
        self.lake_tag[:] = ""
        self.other_tag = np.empty(self.nb_water_pix, dtype=object)
        self.other_tag[:] = ""
        tmp_prior_ice_flag = np.zeros(self.nb_water_pix, dtype=np.uint8) + my_var.FV_NETCDF['uint8']
        
        # 3 - Include river pixels info if there is
        if self.nb_river_pix == 0:
            logger.info("No pixel associated to river")
            
        else:
            
            logger.info("%d pixels associated to rivers", self.nb_river_pix)
            tmp_longitude_vectorproc[self.river_index] = self.longitude_vectorproc
            tmp_latitude_vectorproc[self.river_index] = self.latitude_vectorproc
            tmp_height_vectorproc[self.river_index] = self.height_vectorproc
            tmp_river_reach_tag[self.river_index] = self.river_reach_tag
            tmp_river_node_tag[self.river_index] = self.river_node_tag
            tmp_prior_ice_flag[self.river_index] = self.prior_ice_flag
            
        # 4 - Save arrays
        self.longitude_vectorproc = tmp_longitude_vectorproc
        self.latitude_vectorproc = tmp_latitude_vectorproc
        self.height_vectorproc = tmp_height_vectorproc
        self.river_reach_tag = tmp_river_reach_tag
        self.river_node_tag = tmp_river_node_tag
        self.prior_ice_flag = tmp_prior_ice_flag
        
        # 5 - Save continent var
        self.continent = in_objPixc.continent
        self.pixcvec_metadata["continent"] = self.continent
    
    # ----------------------------------------
    
    def getTileInfo(self):
        """
        Getter of cycle_num, pass_num and tile_ref
        """
        return self.cycle_num, self.pass_num, self.tile_ref
        
    # ----------------------------------------
    
    def write_file(self, in_filename, in_proc_metadata):
        """
        Write the pixel cloud vector attribute product (L2_HR_PIXCVec product)

        :param in_filename: full path of the output file
        :type in_filename: string
        :param in_proc_metadata: processing metadata
        :type in_proc_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Output L2_HR_PIXCVec NetCDF file = %s", in_filename)
        
        # 1 - Init product
        pixcvec_file = nc_file.LakeTilePixcvec_product(in_pixcvecriver_metadata=self.pixcvec_metadata, 
                                                       in_proc_metadata=in_proc_metadata)
        
        # 2 - Form dictionary with variables to write
        vars_to_write = {}
        vars_to_write["azimuth_index"] = self.azimuth_index
        vars_to_write["range_index"] = self.range_index
        vars_to_write["latitude_vectorproc"] = self.latitude_vectorproc
        vars_to_write["longitude_vectorproc"] = self.longitude_vectorproc
        vars_to_write["height_vectorproc"] = self.height_vectorproc
        vars_to_write["river_reach_tag"] = self.river_reach_tag
        vars_to_write["river_node_tag"] = self.river_node_tag
        vars_to_write["lake_tag"] = self.lake_tag
        vars_to_write["other_tag"] = self.other_tag
        vars_to_write["prior_ice_flag"] = self.prior_ice_flag
        
        # 3 - Write file
        pixcvec_file.write_product(in_filename, self.nb_water_pix, vars_to_write)

    def write_file_asShp(self, in_filename, in_objPixc):
        """
        Write the pixel cloud vector product as a shapefile

        :param in_filename: full path of the output file
        :type in_filename: string
        :param in_objPixc: PixelCloud object associated to this PixelCloudVec object
        :type in_objPixc: proc_pixc.PixelCloud
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Output L2_HR_PIXCVec shapefile = %s", in_filename)
        
        # 1 - Init output file
        # 1.1 - Driver
        shpDriver = ogr.GetDriverByName(str("ESRI Shapefile"))
        # 1.2 - Create file
        if os.path.exists(in_filename):
            shpDriver.DeleteDataSource(in_filename)
        outDataSource = shpDriver.CreateDataSource(in_filename)
        # 1.3 - Create layer
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)  # WGS84
        outLayer = outDataSource.CreateLayer(str(str(os.path.basename(in_filename)).replace('.shp', '')), srs, geom_type=ogr.wkbPoint)
        # 1.4 - Create attributes
        outLayer.CreateField(ogr.FieldDefn(str('az_index'), ogr.OFTInteger))
        outLayer.CreateField(ogr.FieldDefn(str('r_index'), ogr.OFTInteger))
        outLayer.CreateField(ogr.FieldDefn(str('classif'), ogr.OFTInteger))
        tmpField = ogr.FieldDefn(str('height2'), ogr.OFTReal)
        tmpField.SetWidth(12)
        tmpField.SetPrecision(4)
        outLayer.CreateField(tmpField)
        outLayer.CreateField(ogr.FieldDefn(str('reach_tag'), ogr.OFTInteger))
        outLayer.CreateField(ogr.FieldDefn(str('node_tag'), ogr.OFTInteger))
        outLayer.CreateField(ogr.FieldDefn(str('lake_tag'), ogr.OFTString))
        outLayer.CreateField(ogr.FieldDefn(str('other_tag'), ogr.OFTString))
        outLayer.CreateField(ogr.FieldDefn(str('ice_flag'), ogr.OFTInteger))
        outLayerDefn = outLayer.GetLayerDefn()
        
        # 2 - Retrieve indices of pixels having finite longitude and latitude
        tmp_idx = np.where(self.longitude_vectorproc < my_var.FV_DOUBLE)[0]
        tmp_idx2 = np.where(self.latitude_vectorproc[tmp_idx] < my_var.FV_DOUBLE)[0]
        indices_to_write = tmp_idx[tmp_idx2]
        
        # 3 - Conversions
        # 3.1 - Longitudes for objects crossing Greenwich meridian
        tmp_longitude = self.longitude_vectorproc
        if len(self.greenwich_idx) > 0:
            tmp_longitude = my_tools.convert_to_m180_180(self.longitude_vectorproc[self.greenwich_idx])
        # 3.2 - Fill values
        tmp_height2 = my_tools.convert_fillvalue(self.height_vectorproc)
        tmp_reach_tag = my_tools.convert_fillvalue(self.river_reach_tag)
        tmp_node_tag = my_tools.convert_fillvalue(self.river_node_tag)
        tmp_lake_tag = self.lake_tag.astype('U')
        tmp_other_tag = self.other_tag.astype('U')
        tmp_ice_flag = my_tools.convert_fillvalue(self.prior_ice_flag)
        
        # 4 - Process each point with improved geolocation
        for indp in indices_to_write:
            # 4.1 - Create feature
            outFeature = ogr.Feature(outLayerDefn)
            # 4.2 - Set point with improved geoloc as geometry
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(tmp_longitude[indp], self.latitude_vectorproc[indp])
            outFeature.SetGeometry(point)
            # 4.3 - Set attributes
            outFeature.SetField(str('az_index'), int(self.azimuth_index[indp]))
            outFeature.SetField(str('r_index'), int(self.range_index[indp]))
            outFeature.SetField(str('classif'), int(in_objPixc.origin_classif[indp]))
            outFeature.SetField(str('height2'), float(tmp_height2[indp]))
            outFeature.SetField(str('reach_tag'), int(tmp_reach_tag[indp]))
            outFeature.SetField(str('node_tag'), int(tmp_node_tag[indp]))
            outFeature.SetField(str('lake_tag'), str(tmp_lake_tag[indp]))
            outFeature.SetField(str('other_tag'), str(tmp_other_tag[indp]))
            outFeature.SetField(str('ice_flag'), int(tmp_ice_flag[indp]))
            # 4.4 - On ajoute l'objet dans la couche de sortie
            outLayer.CreateFeature(outFeature)
            
        # 5 - Destroy the data sources to free resources
        outDataSource.Destroy()
        
        
#######################################
    
def computeImpGeoloc(in_productType, in_objPixc, in_indices, in_height):
    """
    Refines geolocation for in_objPixc pixels corresponding to indices in_indices (in in_objPix)
        
    :param in_productType: type of product among "SP"=LakeSP and "TILE"=LakeTile
    :type in_productType: string
    :param in_objPixc: pixel cloud from which to compute improved geolocation
    :type in_objPixc: proc_pixc.PixelCloud or proc_pixc_sp.PixelCloudSP object
    :param in_indices: indices of pixels related to the same object
    :type in_indices: 1D-array of int
    :param in_height: new height of each point used for improved geoloc
    :type in_height: 1D-array of float
        
    :return out_lon_corr: improved longitude of pixels of in_indices
    :rtype: 1D-array of float
    :return out_lat_corr: improved latitude of pixels of in_indices
    :rtype: 1D-array of float
    :return out_lat_corr: improved latitude of pixels of in_indices
    :rtype: 1D-array of float
    """
    logger = logging.getLogger("proc_pixc_vec")
    nb_pix = in_indices.size
    logger.debug("> %d PixC to deal with", nb_pix)
        
    # 1 - Prepare data for Damien's algo
    # 1.1 - Convert geodetic coordinates (lat, lon, height) to cartesian coordinates (x, y, z)
    x, y, z = my_geoloc.convert_llh2ecef(in_objPixc.latitude[in_indices], in_objPixc.longitude[in_indices], in_objPixc.height[in_indices], my_var.GEN_RAD_EARTH_EQ, my_var.GEN_RAD_EARTH_POLE)
    # 1.2 - Get position of associated along-track pixels (in cartesian coordinates)
    nadir_x = in_objPixc.nadir_x[in_indices]
    nadir_y = in_objPixc.nadir_y[in_indices]
    nadir_z = in_objPixc.nadir_z[in_indices]
    # 1.3 - Get velocity of associated along-track pixels (in cartesian coordinates)
    nadir_vx = in_objPixc.nadir_vx[in_indices]
    nadir_vy = in_objPixc.nadir_vy[in_indices]
    nadir_vz = in_objPixc.nadir_vz[in_indices]
    # 1.4 - Get distance from satellite to target point
    ri = in_objPixc.near_range + in_objPixc.range_index[in_indices] * my_var.GEN_RANGE_SPACING
    
    # 2 - Use Damien's algo
    # 2.1 - Init output vectors
    out_lat_corr = np.zeros(nb_pix)  # Improved latitudes
    out_lon_corr = np.zeros(nb_pix)  # Improved longitudes
    out_height_corr = np.zeros(nb_pix)  # Improved heights
    # 2.2 - Improve geolocation
    p_final, p_final_llh, h_mu, (iter_grad,nfev_minimize_scalar) = my_geoloc.pointcloud_height_geoloc_vect(np.transpose(np.array([x, y, z])), in_objPixc.height[in_indices],
                                                                                                           np.transpose(np.array([nadir_x, nadir_y, nadir_z])),
                                                                                                           np.transpose(np.array([nadir_vx, nadir_vy, nadir_vz])),
                                                                                                           ri, in_height, 
                                                                                                           recompute_Doppler=True, recompute_R=True, verbose=False, 
                                                                                                           max_iter_grad=1, height_goal=1.e-3, safe_flag=True)
    # 2.3 - Save output variables
    out_lat_corr, out_lon_corr, out_height_corr = np.transpose(p_final_llh)
    
    # 2.4 - Return output between 0 and 360 degrees (pixel cloud format)
    return my_tools.convert_to_0_360(out_lon_corr), out_lat_corr, out_height_corr
