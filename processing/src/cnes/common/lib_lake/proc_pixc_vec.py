# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2019/05/17:version initiale.
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: proc_pixc_vec.py
    :synopsis: Deals with SWOT pixel cloud complementary file
     Created on 2017/09/15

.. moduleauthor: Claire POTTIER - CNES DSO/SI/TR

..
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
import cnes.common.service_config_file as service_config_file
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib_lake.locnes_products_netcdf as nc_file
import cnes.common.lib.my_variables as my_var


class PixelCloudVec(object):
    """
        class PixelCloudVec
    """
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
                - river_idx / 1D-array of int: indices of river pixels within PIXC arrays (= variable named pixc_index in L2_HR_PIXCVecRiver only)
            - From both L2_HR_PIXCVecRiver and LakeTile_pixcvec:
                - range_idx / 1D-array of int: range indices of water pixels (= variable named range_index in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
                - azimuth_idx / 1D-array of int: azimuth indices of water pixels (= variable named azimuth_index in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
                - longitude_vectorproc / 1D-array of float: improved longitude of water pixels (= variable named longitude_vectorproc in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
                - latitude_vectorproc / 1D-array of float: improved latitude of water pixels (= variable named latitude_vectorproc in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
                - height_vectorproc / 1D-array of float: improved height of water pixels (= variable named height_vectorproc in L2_HR_PIXCVecRiver and LakeTile_pixcvec)
                - node_id / 1D-array of float: identifier associated to river node database (= variable named node_index in L2_HR_PIXCVecRiver and node_id in LakeTile_pixcvec)
                - pixcvec_metadata / dict: dictionary of PIXCVec file metadata
            - From L2_HR_PIXC:
                - continent / string: continent covered by the tile (if global var CONTINENT_FILE exists)
            - From processing:
                - nb_water_pix / int: number of water pixels
                - reject_index / 1D-array of int: indices of pixels that are river only, ie not reservoirs or dams
                - nb_river_pix / int: number of river pixels
                - lakedb_id / 1D-array of str: identifier from the lake database (= variable named lakedb_id in LakeTile_pixcvec)
                - lakeobs_id / 1D-array of str: identifier associated to unknown object (= variable named lakeobs_id in LakeTile_pixcvec)
                - ice_clim_flag / 1D-array of int: climatological ice flag
                - ice_dyn_flag / 1D-array of int: dynamical ice flag
        """
        # Get instance of service config file
        self.cfg = service_config_file.get_instance()
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
        self.node_id = None
        self.lakedb_id = None
        self.lakeobs_id = None
        self.ice_clim_flag = None
        self.ice_dyn_flag = None
        self.pixcvec_metadata = {}
        # Fill variables if filename available
        if in_pixcvec_file is not None:
            self.set_from_pixcvec_file(in_pixcvec_file)
            
        # Init dictionary of PIXCVec metadata
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
            self.longitude_vectorproc = my_tools.convert_to_m180_180(pixcvec_reader.getVarValue("longitude_vectorproc"))
            # 4.4 - Latitude
            self.latitude_vectorproc = pixcvec_reader.getVarValue("latitude_vectorproc")
            # 4.5 - Height
            try:
                self.height_vectorproc = pixcvec_reader.getVarValue("wse_vectorproc")  # new variable name
            except:
                logger.debug("wse_vectorproc variable not in PIXCVec file => height_vectorproc used instead")
                self.height_vectorproc = pixcvec_reader.getVarValue("height_vectorproc")  # old variable name
            # 4.6 - References entities
            # Node identifier
            try:
                tmp_id = [str(i) for i in pixcvec_reader.getVarValue("node_id")]  # new variable name
            except:
                logger.debug("node_id variable not in PIXCVec file => node_index used instead")
                tmp_id = [str(i) for i in pixcvec_reader.getVarValue("node_index")]  # old variable name
            self.node_id = np.array(tmp_id, dtype=object) 
            # Specific in LakeTile product for LakeSP product    
            if self.product_type == "SP":
                tmp_id = [str(i) for i in pixcvec_reader.getVarValue("lakedb_id")]
                self.lakedb_id = np.array(tmp_id, dtype=object)
                tmp_id = [str(i) for i in pixcvec_reader.getVarValue("lakeobs_id")]
                self.lakeobs_id = np.array(tmp_id, dtype=object)
            # 4.7 - Ice flags
            # Climato
            try:
                self.ice_clim_flag = pixcvec_reader.getVarValue("ice_clim_flag")
            except:
                logger.debug("ice_clim_flag variable not in PIXCVec file => set to 1D-array of _FillValue")
                self.ice_clim_flag = np.zeros(self.nb_water_pix, dtype=np.uint8) + my_var.FV_NETCDF['uint8']
            # Dynamical
            try:
                self.ice_dyn_flag = pixcvec_reader.getVarValue("ice_dyn_flag")
            except:
                logger.debug("ice_dyn_flag variable not in PIXCVec file => set to 1D-array of _FillValue")
                self.ice_dyn_flag = np.zeros(self.nb_water_pix, dtype=np.uint8) + my_var.FV_NETCDF['uint8']
            
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
        tmp_node_id = np.zeros(self.nb_water_pix, dtype=object)
        tmp_node_id[:] = ""
        self.lakedb_id = np.empty(self.nb_water_pix, dtype=object)
        self.lakedb_id[:] = ""
        self.lakeobs_id = np.empty(self.nb_water_pix, dtype=object)
        self.lakeobs_id[:] = ""
        tmp_ice_clim_flag = np.zeros(self.nb_water_pix, dtype=np.uint8) + my_var.FV_NETCDF['uint8']
        tmp_ice_dyn_flag = np.zeros(self.nb_water_pix, dtype=np.uint8) + my_var.FV_NETCDF['uint8']
        
        # 3 - Include river pixels info if there is
        if self.nb_river_pix == 0:
            logger.info("No pixel associated to river")
            
        else:
            
            logger.info("%d pixels associated to rivers", self.nb_river_pix)
            tmp_longitude_vectorproc[self.river_index] = self.longitude_vectorproc
            tmp_latitude_vectorproc[self.river_index] = self.latitude_vectorproc
            tmp_height_vectorproc[self.river_index] = self.height_vectorproc
            tmp_node_id[self.river_index] = self.node_id
            tmp_ice_clim_flag[self.river_index] = self.ice_clim_flag
            tmp_ice_dyn_flag[self.river_index] = self.ice_dyn_flag
            
        # 4 - Save arrays
        self.longitude_vectorproc = tmp_longitude_vectorproc
        self.latitude_vectorproc = tmp_latitude_vectorproc
        self.height_vectorproc = tmp_height_vectorproc
        self.node_id = tmp_node_id
        self.ice_clim_flag = tmp_ice_clim_flag
        self.ice_dyn_flag = tmp_ice_dyn_flag
        
        # 5 - Save continent var
        self.continent = in_objPixc.continent
        self.pixcvec_metadata["continent"] = self.continent
    
    # ----------------------------------------
    
    def getTileInfo(self):
        """
        Getter of cycle_num, pass_num and tile_ref
        """
        return self.pixcvec_metadata["cycle_number"], self.pixcvec_metadata["pass_number"], self.pixcvec_metadata["tile_number"]
#        return self.cycle_num, self.pass_num, self.tile_ref
        
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
        vars_to_write["wse_vectorproc"] = self.height_vectorproc
        vars_to_write["node_id"] = self.node_id
        vars_to_write["lakedb_id"] = self.lakedb_id
        vars_to_write["lakeobs_id"] = self.lakeobs_id
        vars_to_write["ice_climatological_flag"] = self.ice_clim_flag
        vars_to_write["ice_dynamical_flag"] = self.ice_dyn_flag
        
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
        if self.product_type == "TILE" :
            outLayer.CreateField(ogr.FieldDefn(str('classif'), ogr.OFTInteger))
        tmpField = ogr.FieldDefn(str('height2'), ogr.OFTReal)
        tmpField.SetWidth(12)
        tmpField.SetPrecision(4)
        outLayer.CreateField(tmpField)
        outLayer.CreateField(ogr.FieldDefn(str('node_id'), ogr.OFTInteger))
        outLayer.CreateField(ogr.FieldDefn(str('lakedb_id'), ogr.OFTString))
        outLayer.CreateField(ogr.FieldDefn(str('lakeobs_id'), ogr.OFTString))
        outLayer.CreateField(ogr.FieldDefn(str('ice_clim_f'), ogr.OFTInteger))
        outLayer.CreateField(ogr.FieldDefn(str('ice_dyn_f'), ogr.OFTInteger))
        outLayerDefn = outLayer.GetLayerDefn()
        
        # 2 - Retrieve indices of pixels having finite longitude and latitude
        tmp_idx = np.where(self.longitude_vectorproc < my_var.FV_DOUBLE)[0]
        tmp_idx2 = np.where(self.latitude_vectorproc[tmp_idx] < my_var.FV_DOUBLE)[0]
        indices_to_write = tmp_idx[tmp_idx2]

        # 3 - Conversions of fill values
        tmp_height2 = my_tools.convert_fillvalue(self.height_vectorproc)
        tmp_node_id = self.node_id.astype('U')
        tmp_lakedb_id = self.lakedb_id.astype('U')
        tmp_lakeobs_id = self.lakeobs_id.astype('U')
        tmp_ice_clim_flag = my_tools.convert_fillvalue(self.ice_clim_flag)
        tmp_ice_dyn_flag = my_tools.convert_fillvalue(self.ice_dyn_flag)
        
        for indp in indices_to_write:
            # 4.1 - Create feature
            outFeature = ogr.Feature(outLayerDefn)
            # 4.2 - Set point with improved geoloc as geometry
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(self.longitude_vectorproc[indp], self.latitude_vectorproc[indp])
            outFeature.SetGeometry(point)
            # 4.3 - Set attributes
            outFeature.SetField(str('az_index'), int(self.azimuth_index[indp]))
            outFeature.SetField(str('r_index'), int(self.range_index[indp]))
            if self.product_type == "TILE":
                outFeature.SetField(str('classif'), int(in_objPixc.origin_classif[indp]))
            outFeature.SetField(str('height2'), float(tmp_height2[indp]))
            outFeature.SetField(str('node_id'), str(tmp_node_id[indp]))
            outFeature.SetField(str('lakedb_id'), str(tmp_lakedb_id[indp]))
            outFeature.SetField(str('lakeobs_id'), str(tmp_lakeobs_id[indp]))
            outFeature.SetField(str('ice_clim_f'), int(tmp_ice_clim_flag[indp]))
            outFeature.SetField(str('ice_dyn_f'), int(tmp_ice_dyn_flag[indp]))
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
    x, y, z = my_geoloc.convert_llh2ecef(in_objPixc.latitude[in_indices], in_objPixc.longitude[in_indices],
                                         in_objPixc.height[in_indices], my_var.GEN_RAD_EARTH_EQ,
                                         my_var.GEN_RAD_EARTH_POLE)
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
                                                                                                           recompute_doppler=True, recompute_range=True, verbose=False,
                                                                                                           max_iter_grad=1, height_goal=1.e-3)
    # 2.3 - Save output variables
    out_lat_corr, out_lon_corr, out_height_corr = np.transpose(p_final_llh)
    
    # 2.4 - Return output (pixel cloud format)
    return out_lon_corr, out_lat_corr, out_height_corr
