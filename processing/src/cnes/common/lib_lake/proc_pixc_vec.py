# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2019/05/17:version initiale.
# VERSION:2.0.0:DM:#91:2020/07/03:Poursuite industrialisation
# VERSION:3.0.0:DM:#91:2021/03/12:Poursuite industrialisation
# VERSION:3.1.0:DM:#91:2021/05/21:Poursuite industrialisation
# VERSION:4.0.0:DM:#91:2022/05/05:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: proc_pixc_vec.py
    :synopsis: Deals with SWOT pixel cloud complementary file, i.e. PIXCVec and PIXCVecRiver products
     Created on 2017/09/15

.. moduleauthor: Claire POTTIER - CNES DSO/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""

import logging
import netCDF4
import numpy as np
import os
from osgeo import ogr, osr

import cnes.common.service_config_file as service_config_file

import cnes.modules.geoloc.lib.geoloc as my_geoloc

import cnes.common.lib.my_netcdf_file as my_nc
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_variables as my_var
import cnes.common.lib_lake.locnes_products_netcdf as nc_file


class PixelCloudVec(object):
    """
    class PixelCloudVec
    Manage PIXCVec and PIXCVecRiver products
    """
    
    def __init__(self, in_product_type):
        """
        Constructor: init PIXCVec(River) object
        
        :param in_product_type: type of product among "SP"=LakeSP and "TILE"=LakeTile
        :type in_product_type: string

        Variables of the object:
            - From L2_HR_PIXCVecRiver:
                - river_index / 1D-array of int: indices of river pixels within PIXC arrays (= variable named pixc_index in L2_HR_PIXCVecRiver only)
            - From both L2_HR_PIXCVecRiver and LakeTile_PIXCVec:
                - range_index / 1D-array of int: range indices of water pixels (= variable named range_index in L2_HR_PIXCVecRiver 
                  and LakeTile_PIXCVec)
                - azimuth_index / 1D-array of int: azimuth indices of water pixels (= variable named azimuth_index in L2_HR_PIXCVecRiver 
                  and LakeTile_PIXCVec)
                - longitude_vectorproc / 1D-array of float: improved longitude of water pixels (= variable named longitude_vectorproc 
                  in L2_HR_PIXCVecRiver and LakeTile_PIXCVec)
                - latitude_vectorproc / 1D-array of float: improved latitude of water pixels (= variable named latitude_vectorproc 
                  in L2_HR_PIXCVecRiver and LakeTile_PIXCVec)
                - height_vectorproc / 1D-array of float: improved height of water pixels (= variable named height_vectorproc 
                  in L2_HR_PIXCVecRiver and LakeTile_PIXCVec)
                - reach_id / 1D-array of int or char(11): identifier associated to river reach database (= variable named reach_id 
                  in L2_HR_PIXCVecRiver and LakeTile_PIXCVec)
                - node_id / 1D-array of int or char(14): identifier associated to river node database (= variable named node_id 
                  in L2_HR_PIXCVecRiver and LakeTile_PIXCVec)
                - ice_clim_f / 1D-array of int: climatological ice flag
                - ice_dyn_f / 1D-array of int: dynamical ice flag
                - pixcvec_metadata / dict: dictionary of PIXCVec file metadata
            - From L2_HR_PIXC:
                - continent / string: continent covered by the tile (if global var CONTINENT_FILE exists)
                - nadir_time / 1D-array of float: UTC time of all pixels 
            - From processing:
                - nb_water_pix / int: number of water pixels
                - lake_id / 1D-array of char(10): identifier from the lake database (= variable named lake_id in LakeTile_PIXCVec)
                - obs_id / 1D-array of char(13): identifier associated to unknown object (= variable named obs_id in LakeTile_PIXCVec)
                - reject_index / 1D-array of int: indices of pixels that are river only, ie not connected lakes or dams
                - nb_river_pix / int: number of river pixels
        """
        # Get instance of service config file
        self.cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # Init type
        if in_product_type in ("TILE", "SP"):
            self.product_type = in_product_type
        else:
            message = "Product type %s unknown; should be TILE or SP"
            logger.error(message, exc_info=True)
            raise
            
        # Init dimension
        self.nb_water_pix = 0
        
        # Init PIXCVec variables
        self.azimuth_index = None
        self.range_index = None
        self.longitude_vectorproc = None
        self.latitude_vectorproc = None
        self.height_vectorproc = None
        self.reach_id = None
        self.node_id = None
        self.lake_id = None
        self.obs_id = None
        self.ice_clim_f = None
        self.ice_dyn_f = None
            
        # Init dictionary of PIXCVec metadata
        self.pixcvec_metadata = {}
        self.pixcvec_metadata["source"] = ""
        self.pixcvec_metadata["cycle_number"] = -9999
        self.pixcvec_metadata["pass_number"] = -9999
        self.pixcvec_metadata["tile_number"] = -9999
        self.pixcvec_metadata["swath_side"] = ""
        self.pixcvec_metadata["tile_name"] = ""
        self.pixcvec_metadata["continent_id"] = ""
        self.pixcvec_metadata["continent_code"] = ""
        self.pixcvec_metadata["time_granule_start"] = ""
        self.pixcvec_metadata["time_granule_end"] = ""
        self.pixcvec_metadata["time_coverage_start"] = ""
        self.pixcvec_metadata["time_coverage_end"] = ""
        self.pixcvec_metadata["geospatial_lon_min"] = -9999.0
        self.pixcvec_metadata["geospatial_lon_max"] = -9999.0
        self.pixcvec_metadata["geospatial_lat_min"] = -9999.0
        self.pixcvec_metadata["geospatial_lat_max"] = -9999.0
        self.pixcvec_metadata["inner_first_latitude"] = -9999.0
        self.pixcvec_metadata["inner_first_longitude"] = -9999.0
        self.pixcvec_metadata["inner_last_latitude"] = -9999.0
        self.pixcvec_metadata["inner_last_longitude"] = -9999.0
        self.pixcvec_metadata["outer_first_latitude"] = -9999.0
        self.pixcvec_metadata["outer_first_longitude"] = -9999.0
        self.pixcvec_metadata["outer_last_latitude"] = -9999.0
        self.pixcvec_metadata["outer_last_longitude"] = -9999.0
        self.pixcvec_metadata["xref_l2_hr_pixc_file"] = ""
        self.pixcvec_metadata["xref_l2_hr_pixcvecriver_file"] = ""
        self.pixcvec_metadata["xref_prior_river_db_file"] = ""
        self.pixcvec_metadata["xref_prior_lake_db_file"] = ""
        self.pixcvec_metadata["xref_reforbittrack_files"] = ""
        self.pixcvec_metadata["xref_param_l2_hr_laketile_file"] = ""
        self.pixcvec_metadata["ellipsoid_semi_major_axis"] = -9999.0
        self.pixcvec_metadata["ellipsoid_flattening"] = -9999.0
        
        # Variables specific to processing
        self.continent_id = None
        # Specific to LakeTile processing
        if self.product_type == "TILE":
            self.nb_river_pix = 0  # Number of river pixels (used for TILE processing)
            self.river_index = None  # Indices of pixels processed by RiverTile (used in TILE processing)
            self.reject_index = None  # Indices of river pixels (not connected lakes)
            
    def set_from_pixcvec_file(self, in_pixcvec_file):
        """
        Set variables from PIXCVec file
        
        :param in_pixcvec_file: full path of pixel cloud complementary file 
                                    (L2_HR_PIXCVecRiver file if from PGE_RiverTile 
                                    or LakeTile_PIXCVec if from PGE_LakeTile)
        :type in_pixcvec_file: string
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        if self.product_type == "TILE":
            logger.debug("L2_HR_PIXCVecRiver file = %s", in_pixcvec_file)
        elif self.product_type == "SP":
            logger.debug("LakeTile_PIXCVec file = %s", in_pixcvec_file)
        else:
            message = "Product type %s unknown; should be TILE or SP"
            logger.error(message, exc_info=True)
            raise
        
        # 1 - Open file in reading mode
        pixcvec_reader = my_nc.MyNcReader(in_pixcvec_file)
        
        # 2 - Retrieve the number of records
        self.nb_water_pix = pixcvec_reader.get_dim_value("points")
            
        # 3 - Retrieve needed global attributes
        pixcvec_keys = pixcvec_reader.get_list_att()
        for key in self.pixcvec_metadata.keys():
            if key in pixcvec_keys:
                self.pixcvec_metadata[key] = pixcvec_reader.get_att_value(key)
        
        # 4 - Retrieve variables if there are river pixels
        if self.nb_water_pix == 0:
            logger.debug("NO pixel associated to rivers")
            
        else:
                
            # 4.1 - Range index
            self.range_index = pixcvec_reader.get_var_value("range_index")
            # 4.2 - Azimuth index
            self.azimuth_index = pixcvec_reader.get_var_value("azimuth_index")
            # 4.3 - Longitude
            self.longitude_vectorproc = my_tools.convert_to_m180_180(pixcvec_reader.get_var_value("longitude_vectorproc"))
            # 4.4 - Latitude
            self.latitude_vectorproc = pixcvec_reader.get_var_value("latitude_vectorproc")
            # 4.5 - Height
            self.height_vectorproc = pixcvec_reader.get_var_value("height_vectorproc")
            # 4.6 - References entities
            # Reach identifier
            if self.product_type == "TILE":
                tmp_id = pixcvec_reader.get_var_value("reach_id")
            elif self.product_type == "SP":
                tmp_id = netCDF4.chartostring(pixcvec_reader.get_var_value("reach_id"))
            self.reach_id = np.char.asarray(tmp_id, itemsize=11)
            # Node identifier
            if self.product_type == "TILE":
                tmp_id = pixcvec_reader.get_var_value("node_id")
            elif self.product_type == "SP":
                tmp_id = netCDF4.chartostring(pixcvec_reader.get_var_value("node_id"))
            self.node_id = np.char.asarray(tmp_id, itemsize=14) 
            # Specific in LakeTile product for LakeSP product    
            if self.product_type == "SP":
                tmp_id = netCDF4.chartostring(pixcvec_reader.get_var_value("lake_id"))
                self.lake_id = np.char.asarray(tmp_id, itemsize=10)
                tmp_id = netCDF4.chartostring(pixcvec_reader.get_var_value("obs_id"))
                self.obs_id = np.char.asarray(tmp_id, itemsize=13)
            # 4.7 - Ice flags
            # Climato
            self.ice_clim_f = pixcvec_reader.get_var_value("ice_clim_f")
            # Dynamical
            self.ice_dyn_f = pixcvec_reader.get_var_value("ice_dyn_f")
            
            # 5 - Compute indices of pixels to remove from lake processing (= river pixels - connected lakes)
            if self.product_type == "TILE":
                self.compute_pix_to_reject(pixcvec_reader)
                 
        # 6 - Close file
        pixcvec_reader.close()
    
    def compute_pix_to_reject(self, in_pixcvec_reader):
        """
        Compute indices of pixels to remove from LakeTile processing, i.e.:
        - pixels already processed by RiverTile with correct output values (except connected lakes with reach_id finishing by Type digit = 3)
        
        :param in_pixcvec_reader: reader of L2_HR_PIXCVecRiver file
        :type in_pixcvec_reader: my_netcdf_file.MyNcReader
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # =0 for valid river pixels; =1 for bad river pixels
        # River pixels with classif == 1 will be processed again by lake processor
        classif = np.zeros(self.nb_water_pix)
        # 1.1 - Remove pixels with range_index and azimuth_index = _FillValue
        vars_to_look_for_fillvalue = {}
        vars_to_look_for_fillvalue["range_index"] = self.range_index
        vars_to_look_for_fillvalue["azimuth_index"] = self.azimuth_index
        vars_to_look_for_fillvalue["pixc_index"] = in_pixcvec_reader.get_var_value("pixc_index")
        for var_name, var_value in vars_to_look_for_fillvalue.items():
            fv_idx = np.argwhere(var_value > 2e9)
            nb_fv = len(fv_idx)
            if nb_fv != 0:
                logger.debug("%d pixels have _FillValue in %s variable => will be rejected" % (nb_fv, var_name))
                classif[fv_idx] = 1
        # 1.2 - Remove pixels with longitude, latitude, and height = numpy.nan
        vars_to_look_for_nans = {}
        vars_to_look_for_nans["longitude"] = self.longitude_vectorproc
        vars_to_look_for_nans["latitude"] = self.latitude_vectorproc
        vars_to_look_for_nans["height"] = self.height_vectorproc
        for var_name, var_value in vars_to_look_for_nans.items():
            nan_idx = np.argwhere(np.isnan(var_value))
            nb_nan = len(nan_idx)
            if nb_nan != 0:
                logger.debug("%d pixels have NaN in %s variable => will be rejected" % (nb_nan, var_name))
                classif[nan_idx] = 1
                
        # 2 - Indices of pixels of PixC already processed by PGE_RiverTile
        tmp_river_index = in_pixcvec_reader.get_var_value("pixc_index")
        ind_good_pixcvec = np.where(classif == 0)[0]
        self.river_index = tmp_river_index[ind_good_pixcvec]
        # Number of pixels of PixC already processed by PGE_RiverTile
        self.nb_river_pix = self.river_index.size
        
        if self.nb_river_pix == 0:
            logger.debug("No pixel associated to river")
            
        else:
                
            # Indices of PIXC to remove from LakeTile processing = river only pixels (ie connected lakes kept)
            ind_type_not_3 = []
            for ind in ind_good_pixcvec:
                # OLD = Check for FillValue in index list
                # OLD = if (not self.reach_id[ind][-1:] == b'3' and self.river_index[ind] < my_var.FV_NETCDF['int']):
                if not self.reach_id[ind][-1:] == b'3':
                    ind_type_not_3.append(ind)
            self.reject_index = self.river_index[ind_type_not_3]  # reach_id not ending with 3 (connected lakes)
            
            logger.debug("%d pixels associated to rivers", self.nb_river_pix)
            logger.debug("%d pixels associated to connected lakes", self.nb_river_pix-self.reject_index.size)
        
    # ----------------------------------------
    
    def reshape(self, in_obj_pixc):
        """
        Reshape PIXCVecRiver arrays to new arrays of size of related PIXC arrays
        
        :param in_obj_pixc: pixel cloud associated to current PIXCVecRiver object
        :type in_obj_pixc: proc_pixc.PixelCloud
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # 1 - Init azimuth and range indices arrays, and number of water pixels
        self.range_index = in_obj_pixc.origin_range_index
        self.azimuth_index = in_obj_pixc.origin_azimuth_index
        self.nb_water_pix = len(self.azimuth_index)
        
        # 2 - Init the other arrays
        tmp_longitude_vectorproc = np.zeros(self.nb_water_pix, dtype=np.float64) + my_var.FV_NETCDF['float64']
        tmp_latitude_vectorproc = np.zeros(self.nb_water_pix, dtype=np.float64) + my_var.FV_NETCDF['float64']
        tmp_height_vectorproc = np.zeros(self.nb_water_pix, dtype=np.float32) + my_var.FV_NETCDF['float32']
        tmp_reach_id = np.chararray([self.nb_water_pix], itemsize=11)
        tmp_reach_id[:] = my_var.FV_NETCDF['str']
        tmp_node_id = np.chararray([self.nb_water_pix], itemsize=14)
        tmp_node_id[:] = my_var.FV_NETCDF['str']
        self.lake_id = np.chararray([self.nb_water_pix], itemsize=10)
        self.lake_id[:] = my_var.FV_NETCDF['str']
        self.obs_id = np.chararray([self.nb_water_pix], itemsize=13)
        self.obs_id[:] = my_var.FV_NETCDF['str']
        tmp_ice_clim_f = np.zeros(self.nb_water_pix, dtype=np.uint8) + my_var.FV_NETCDF['uint8']
        tmp_ice_dyn_f = np.zeros(self.nb_water_pix, dtype=np.uint8) + my_var.FV_NETCDF['uint8']
        
        # 3 - Include river pixels info if there is
        if self.nb_river_pix == 0:
            logger.debug("No pixel associated to river")
            
        else:
            """
            OLD
            # Manage _FillValue in PIXCVecRiver
            # Get valid index, ie not a _FillValue
            valid_index = []
            for ind in range(len(self.river_index)):
                if (self.river_index[ind] < my_var.FV_NETCDF['int']):
                    valid_index.append(ind)
            logger.debug("%d valid pixels associated to rivers", len(self.river_index[valid_index]))
            # Copy only valid value
            tmp_longitude_vectorproc[self.river_index[valid_index]] = self.longitude_vectorproc[valid_index]
            tmp_latitude_vectorproc[self.river_index[valid_index]] = self.latitude_vectorproc[valid_index]
            tmp_height_vectorproc[self.river_index[valid_index]] = self.height_vectorproc[valid_index]
            tmp_reach_id[self.river_index[valid_index]] = self.reach_id[valid_index]
            tmp_node_id[self.river_index[valid_index]] = self.node_id[valid_index]
            tmp_ice_clim_f[self.river_index[valid_index]] = self.ice_clim_f[valid_index]
            tmp_ice_dyn_f[self.river_index[valid_index]] = self.ice_dyn_f[valid_index]
            """
            logger.debug("%d pixels associated to rivers", self.nb_river_pix)
            tmp_longitude_vectorproc[self.river_index] = self.longitude_vectorproc
            tmp_latitude_vectorproc[self.river_index] = self.latitude_vectorproc
            tmp_height_vectorproc[self.river_index] = self.height_vectorproc
            tmp_reach_id[self.river_index] = self.reach_id
            tmp_node_id[self.river_index] = self.node_id
            tmp_ice_clim_f[self.river_index] = self.ice_clim_f
            tmp_ice_dyn_f[self.river_index] = self.ice_dyn_f
            tmp_ice_dyn_f[self.river_index] = self.ice_dyn_f
            
        # 4 - Save arrays
        self.longitude_vectorproc = tmp_longitude_vectorproc
        self.latitude_vectorproc = tmp_latitude_vectorproc
        self.height_vectorproc = tmp_height_vectorproc
        self.reach_id = tmp_reach_id
        self.node_id = tmp_node_id
        self.ice_clim_f = tmp_ice_clim_f
        self.ice_dyn_f = tmp_ice_dyn_f
    
    # ----------------------------------------
    
    def get_tile_info(self):
        """
        Getter of cycle_num, pass_num and tile_ref
        """
        return self.pixcvec_metadata["cycle_number"], self.pixcvec_metadata["pass_number"], self.pixcvec_metadata["tile_number"]
        
    # ----------------------------------------
    
    def update_metadata(self, in_obj_pixc):
        """
        Update PIXCVec global metadata with global attributes of PIXC file
        
        :param in_obj_pixc: PIXC object
        :type in_obj_pixc: lake_tile.proc_pixc
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        if self.product_type != "TILE":
            message = "Product type should be TILE"
            logger.error(message, exc_info=True)
            raise
        
        # 1 - Copy attributes from PIXC 
        pixc_keys = in_obj_pixc.pixc_metadata.keys()
        for key in self.pixcvec_metadata.keys():
            if key in pixc_keys:
                self.pixcvec_metadata[key] = in_obj_pixc.pixc_metadata[key]
                
        # 2 - Update time_coverage
        # 2.1 - Retrieve indices of pixels having finite longitude and latitude
        tmp_idx = np.where(self.longitude_vectorproc < my_var.FV_DOUBLE)[0]
        tmp_idx2 = np.where(self.latitude_vectorproc[tmp_idx] < my_var.FV_DOUBLE)[0]
        self.indices_improved_pixels = tmp_idx[tmp_idx2]
        # 2.2 - Estimate time_coverage_start and time_coverage_end
        if len(self.indices_improved_pixels) > 0:
            nadir_time_improved_pixels = in_obj_pixc.nadir_time_all[self.indices_improved_pixels]
            self.pixcvec_metadata["time_coverage_start"] = my_tools.convert_utc_to_str(min(nadir_time_improved_pixels), in_format=2)
            self.pixcvec_metadata["time_coverage_end"] = my_tools.convert_utc_to_str(max(nadir_time_improved_pixels), in_format=2)
        else:
            self.pixcvec_metadata["time_coverage_start"] = ""
            self.pixcvec_metadata["time_coverage_end"] = ""
        
    # ----------------------------------------
    
    def write_file(self, in_filename, in_proc_metadata, in_compress=True):
        """
        Write the pixel cloud vector attribute product (L2_HR_PIXCVec product)
        Update with PIXC metadata if needed.

        :param in_filename: full path of the output file
        :type in_filename: string
        :param in_proc_metadata: processing metadata
        :type in_proc_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        if self.product_type == "TILE":
            logger.debug("Writing output LakeTile_PIXCVec NetCDF file = %s", in_filename)
        else:
            logger.debug("Writing output L2_HR_PIXCVec NetCDF file = %s", in_filename)
            
        # 0 - Compute geospatial_lon|lat_min|max if not done
        if self.pixcvec_metadata["geospatial_lon_min"] == -9999.0:
            self.pixcvec_metadata["geospatial_lon_min"] = min(self.pixcvec_metadata["inner_first_longitude"],
                                                              self.pixcvec_metadata["inner_last_longitude"],
                                                              self.pixcvec_metadata["outer_first_longitude"],
                                                              self.pixcvec_metadata["outer_last_longitude"])
            self.pixcvec_metadata["geospatial_lon_max"] = max(self.pixcvec_metadata["inner_first_longitude"],
                                                              self.pixcvec_metadata["inner_last_longitude"],
                                                              self.pixcvec_metadata["outer_first_longitude"],
                                                              self.pixcvec_metadata["outer_last_longitude"])
            self.pixcvec_metadata["geospatial_lat_min"] = min(self.pixcvec_metadata["inner_first_latitude"],
                                                              self.pixcvec_metadata["inner_last_latitude"],
                                                              self.pixcvec_metadata["outer_first_latitude"],
                                                              self.pixcvec_metadata["outer_last_latitude"])
            self.pixcvec_metadata["geospatial_lat_max"] = max(self.pixcvec_metadata["inner_first_latitude"],
                                                              self.pixcvec_metadata["inner_last_latitude"],
                                                              self.pixcvec_metadata["outer_first_latitude"],
                                                              self.pixcvec_metadata["outer_last_latitude"])
        
        # 1 - Init product
        if self.product_type == "TILE":
            pixcvec_file = nc_file.LakeTilePixcvecProduct(in_pixcvecriver_metadata=self.pixcvec_metadata, 
                                                          in_proc_metadata=in_proc_metadata)
        else:
            pixcvec_file = nc_file.PixcvecProduct(in_laketile_pixcvec_metadata=self.pixcvec_metadata, 
                                                  in_proc_metadata=in_proc_metadata)
        
        # 2 - Form dictionary with variables to write
        vars_to_write = {}
        vars_to_write["azimuth_index"] = self.azimuth_index
        vars_to_write["range_index"] = self.range_index
        vars_to_write["latitude_vectorproc"] = self.latitude_vectorproc
        vars_to_write["longitude_vectorproc"] = self.longitude_vectorproc
        vars_to_write["height_vectorproc"] = self.height_vectorproc
        vars_to_write["reach_id"] = self.reach_id
        vars_to_write["node_id"] = self.node_id
        if self.lake_id is not None:
            vars_to_write["lake_id"] = self.lake_id
        else:
            self.lake_id = np.chararray([self.nb_water_pix], itemsize=10)
            self.lake_id[:] = my_var.FV_NETCDF['str']
        if self.obs_id is not None:
            vars_to_write["obs_id"] = self.obs_id
        else:
            self.obs_id = np.chararray([self.nb_water_pix], itemsize=10)
            self.obs_id[:] = my_var.FV_NETCDF['str']
        vars_to_write["ice_clim_f"] = self.ice_clim_f
        vars_to_write["ice_dyn_f"] = self.ice_dyn_f
        
        # 3 - Write file
        pixcvec_file.write_product(in_filename, self.nb_water_pix, vars_to_write, in_compress=in_compress)

    def write_file_as_shp(self, in_filename, in_obj_pixc):
        """
        Write the pixel cloud vector product as a shapefile

        :param in_filename: full path of the output file
        :type in_filename: string
        :param in_obj_pixc: PixelCloud object associated to this PixelCloudVec object
        :type in_obj_pixc: proc_pixc.PixelCloud
        """
        logger = logging.getLogger(self.__class__.__name__)
        if self.product_type == "TILE":
            logger.debug("Writing output LakeTile_PIXCVec shapefile = %s", in_filename)
        else:
            logger.debug("Writing output L2_HR_PIXCVec shapefile = %s", in_filename)
        
        # 1 - Init output file
        # 1.1 - Driver
        shp_driver = ogr.GetDriverByName(str("ESRI Shapefile"))
        # 1.2 - Create file
        if os.path.exists(in_filename):
            shp_driver.DeleteDataSource(in_filename)
        out_data_source = shp_driver.CreateDataSource(in_filename)
        # 1.3 - Create layer
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)  # WGS84
        out_layer = out_data_source.CreateLayer(str(str(os.path.basename(in_filename)).replace('.shp', '')), srs, geom_type=ogr.wkbPoint)
        # 1.4 - Create attributes
        out_layer.CreateField(ogr.FieldDefn(str('az_index'), ogr.OFTInteger))
        out_layer.CreateField(ogr.FieldDefn(str('r_index'), ogr.OFTInteger))
        if self.product_type == "TILE":
            out_layer.CreateField(ogr.FieldDefn(str('classif'), ogr.OFTInteger))
        tmp_field = ogr.FieldDefn(str('height2'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(4)
        out_layer.CreateField(tmp_field)
        out_layer.CreateField(ogr.FieldDefn(str('reach_id'), ogr.OFTString))
        out_layer.CreateField(ogr.FieldDefn(str('node_id'), ogr.OFTString))
        out_layer.CreateField(ogr.FieldDefn(str('lake_id'), ogr.OFTString))
        out_layer.CreateField(ogr.FieldDefn(str('obs_id'), ogr.OFTString))
        out_layer.CreateField(ogr.FieldDefn(str('ice_clim_f'), ogr.OFTInteger))
        out_layer.CreateField(ogr.FieldDefn(str('ice_dyn_f'), ogr.OFTInteger))
        out_layer_defn = out_layer.GetLayerDefn()
        
        if self.nb_water_pix != 0:
    
            # 2 - Conversions of fill values
            tmp_height2 = my_tools.convert_fillvalue(self.height_vectorproc)
            tmp_reach_id = self.reach_id.astype('U')
            tmp_node_id = self.node_id.astype('U')
            tmp_lake_id = self.lake_id.astype('U')
            tmp_obs_id = self.obs_id.astype('U')
            tmp_ice_clim_f = my_tools.convert_fillvalue(self.ice_clim_f)
            tmp_ice_dyn_f = my_tools.convert_fillvalue(self.ice_dyn_f)
            # 2.1 - Retrieve indices of pixels having finite longitude and latitude
            if self.product_type == "SP":
                tmp_idx = np.where(self.longitude_vectorproc < my_var.FV_DOUBLE)[0]
                tmp_idx2 = np.where(self.latitude_vectorproc[tmp_idx] < my_var.FV_DOUBLE)[0]
                self.indices_improved_pixels = tmp_idx[tmp_idx2]
            
            for indp in self.indices_improved_pixels:
                # 3.1 - Create feature
                out_feature = ogr.Feature(out_layer_defn)
                # 3.2 - Set point with improved geoloc as geometry
                point = ogr.Geometry(ogr.wkbPoint)
                point.AddPoint(self.longitude_vectorproc[indp], self.latitude_vectorproc[indp])
                out_feature.SetGeometry(point)
                # 3.3 - Set attributes
                out_feature.SetField(str('az_index'), int(self.azimuth_index[indp]))
                out_feature.SetField(str('r_index'), int(self.range_index[indp]))
                if self.product_type == "TILE":
                    out_feature.SetField(str('classif'), int(in_obj_pixc.origin_classif[indp]))
                out_feature.SetField(str('height2'), float(tmp_height2[indp]))
                out_feature.SetField(str('reach_id'), str(tmp_reach_id[indp]))
                out_feature.SetField(str('node_id'), str(tmp_node_id[indp]))
                out_feature.SetField(str('lake_id'), str(tmp_lake_id[indp]))
                out_feature.SetField(str('obs_id'), str(tmp_obs_id[indp]))
                out_feature.SetField(str('ice_clim_f'), int(tmp_ice_clim_f[indp]))
                out_feature.SetField(str('ice_dyn_f'), int(tmp_ice_dyn_f[indp]))
                # 3.4 - On ajoute l'objet dans la couche de sortie
                out_layer.CreateFeature(out_feature)

        # 4 - Destroy the data sources to free resources
        out_data_source.Destroy()


#######################################


def compute_imp_geoloc(in_product_type, in_obj_pixc, in_indices, in_height):
    """
    Refines geolocation for in_obj_pixc pixels corresponding to indices in_indices (in in_objPix)
        
    :param in_product_type: type of product among "SP"=LakeSP and "TILE"=LakeTile
    :type in_product_type: string
    :param in_obj_pixc: pixel cloud from which to compute improved geolocation
    :type in_obj_pixc: proc_pixc.PixelCloud or proc_pixc_sp.PixelCloudSP object
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
    :return out_p_final: position of pixels in geocentric coordinates
    :rtype: numpy 2D-array of float; size (3=x/y/z, nb_points)
    """
    logger = logging.getLogger("proc_pixc_vec")
    nb_pix = in_indices.size
    logger.debug("> %d PixC to deal with", nb_pix)
        
    # 1 - Prepare data for Damien's algo
    # 1.1 - Convert geodetic coordinates (lat, lon, height) to cartesian coordinates (x, y, z)
    x, y, z = my_geoloc.convert_llh2ecef(in_obj_pixc.latitude[in_indices], in_obj_pixc.longitude[in_indices],
                                         in_obj_pixc.height[in_indices], my_var.GEN_RAD_EARTH_EQ,
                                         my_var.GEN_RAD_EARTH_POLE)
    # 1.2 - Get position of associated along-track pixels (in cartesian coordinates)
    nadir_x = in_obj_pixc.nadir_x[in_indices]
    nadir_y = in_obj_pixc.nadir_y[in_indices]
    nadir_z = in_obj_pixc.nadir_z[in_indices]
    # 1.3 - Get velocity of associated along-track pixels (in cartesian coordinates)
    nadir_vx = in_obj_pixc.nadir_vx[in_indices]
    nadir_vy = in_obj_pixc.nadir_vy[in_indices]
    nadir_vz = in_obj_pixc.nadir_vz[in_indices]
    # 1.4 - Get distance from satellite to target point
    if in_product_type == "TILE":
        ri = in_obj_pixc.near_range + in_obj_pixc.range_index[in_indices] * my_var.GEN_RANGE_SPACING
    else:
        ri = in_obj_pixc.get_near_range(in_indices) + in_obj_pixc.range_index[in_indices] * my_var.GEN_RANGE_SPACING

    # 2 - Use Damien's algo
    # 2.1 - Init output vectors
    out_lat_corr = np.zeros(nb_pix)  # Improved latitudes
    out_lon_corr = np.zeros(nb_pix)  # Improved longitudes
    out_height_corr = np.zeros(nb_pix)  # Improved heights
    # 2.2 - Improve geolocation
    out_p_final, p_final_llh, h_mu, (iter_grad, nfev_minimize_scalar) = \
                      my_geoloc.pointcloud_height_geoloc_vect(np.transpose(np.array([x, y, z])), in_obj_pixc.height[in_indices],
                                                              np.transpose(np.array([nadir_x, nadir_y, nadir_z])),
                                                              np.transpose(np.array([nadir_vx, nadir_vy, nadir_vz])),
                                                              ri, in_height, recompute_doppler=True, recompute_range=True, 
                                                              verbose=False, max_iter_grad=1, height_goal=1.e-3)
    # 2.3 - Save output variables
    out_lat_corr, out_lon_corr, out_height_corr = np.transpose(p_final_llh)
    
    # 2.4 - Return output (pixel cloud format)
    return out_lon_corr, out_lat_corr, out_height_corr, out_p_final
