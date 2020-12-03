#!/usr/bin/env python
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
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: pge_lake_tile.py
    :synopsis: Process PGE_L2_HR_LakeTile, i.e. generate L2_HR_LakeTile
    product from one tile of L2_HR_PIXC product and associated
    L2_HR_PIXCVecRiver product

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
import configparser
import datetime
import logging
import os
import sys

import cnes.common.service_config_file as service_config_file
import cnes.common.service_error as service_error
import cnes.common.service_logger as service_logger

import cnes.common.lib.my_timer as my_timer
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib_lake.lake_db as lake_db
import cnes.common.lib_lake.locnes_filenames as locnes_filenames
import cnes.common.lib_lake.proc_lake as proc_lake
import cnes.common.lib_lake.proc_pixc_vec as proc_pixc_vec

import cnes.sas.lake_tile.sas_lake_tile as sas_lake_tile
import cnes.sas.lake_tile.proc_pixc as proc_pixc


#######################################


class PGELakeTile():
    """
    Class PGELakeTile
    Pseudo PGE class to launch LakeTile SAS processing
    """

    def __init__(self, cmd_file):
        """
        Constructor of PGELakeTile

        :param cmdFile: command file full path
        :type cmdFile: string
        """

        # 0 - Init timer
        self.timer = my_timer.Timer()
        self.timer.start()
        
        # 1 - Load command file
        self.cmd_file = cmd_file
        my_tools.test_file(cmd_file, in_extent=".cfg")  # Test existance and extension
        my_params = self._read_cmd_file()  # Read parameters
        self.pixc_file = my_params["pixc_file"]
        self.pixc_vec_river_file = my_params["pixc_vec_river_file"]
        self.output_dir = my_params["output_dir"]

        # 2 - Load parameter file
        # 2.1 - Read value from command file
        file_config = my_params["param_file"]
        # 2.2 - Test existence
        if not os.path.exists(file_config):
            raise service_error.DirFileError(file_config)
        # 2.3 - Load parameters
        self.cfg = service_config_file.ServiceConfigFile(file_config)

        # 3 - Put command parameter inside cfg
        self._put_cmd_value(my_params)
        # Retrieve options which are optionnal
        self.flag_prod_shp = self.cfg.getboolean("OPTIONS", "Produce shp")
        self.flag_inc_file_counter = self.cfg.getboolean("OPTIONS", "Increment file counter")
        
        # 4 - Init data object
        self.obj_pixc_vec = None
        self.obj_pixc = None
        self.obj_lake_db = None
        self.obj_lake = None
        self.lake_tile_filenames = None

        # 5 - Initiate logging service
        service_logger.ServiceLogger()
        logger = logging.getLogger(self.__class__.__name__)

        # 6 - Print info
        logger.sigmsg("======================================")
        logger.sigmsg("===== LakeTileProcessing = BEGIN =====")
        logger.sigmsg("======================================")
        logger.info("> INPUT PIXC file = " + str(self.cfg.get("PATHS", "PIXC file")))
        logger.info("> INPUT PIXCVecRiver file = " + str(self.cfg.get("PATHS", "PIXCVecRiver file")))
        logger.info("> OUTPUT directory = " + str(self.cfg.get("PATHS", "Output directory")))
        logger.info("======================================")
        message = "> Command file: " + str(self.cmd_file)
        logger.info(message)
        message = "> " + str(self.cfg)
        logger.info(message)
        logger.info("======================================")
        
        # 7 - Test input parameters
        logger.info(">> Check config parameters type and value")
        self._check_config_parameters()

        # 8 - Form processing metadata dictionary
        self.proc_metadata = {}
        if self.cfg.get("DATABASES", "LAKE_DB") is not None:
            self.proc_metadata["xref_static_lake_db_file"] = self.cfg.get("DATABASES", "LAKE_DB")
        else:
            self.proc_metadata["xref_static_lake_db_file"] = ""
        self.proc_metadata["xref_input_l2_hr_pixc_file"] = self.pixc_file
        self.proc_metadata["xref_input_l2_hr_pixc_vec_river_file"] = self.pixc_vec_river_file
        self.proc_metadata["xref_l2_hr_lake_tile_config_parameter_file"] = file_config
        
        logger.info("")
        
    # -------------------------------------------

    def _read_cmd_file(self):
        """
        Read the parameter file and store parameters in a dictionary

        :return: dictionary containing parameters
        :rtype: dict
        """
        
        # 0 - Init output dictionary
        out_params = {}
        # Default values
        out_params["flag_prod_shp"] = "False"
        out_params["flag_inc_file_counter"] = "True"

        # 1 - Read parameter file
        config = configparser.ConfigParser()
        config.read(self.cmd_file)

        # 2 - Retrieve PATHS
        try:
            out_params["param_file"] = os.path.expandvars(config.get("PATHS", "param_file"))
        except:
            out_params["param_file"] = os.path.join(sys.path[0], "lake_tile_param.cfg")
            #out_params["param_file"] = os.path.join(os.getcwd(), "lake_tile_param.cfg")
        out_params["pixc_file"] = config.get("PATHS", "PIXC file")
        out_params["pixc_vec_river_file"] = config.get("PATHS", "PIXCVecRiver file")
        out_params["output_dir"] = config.get("PATHS", "Output directory")

        # 3 - Retrieve DATABASES
        out_params["LAKE_DB"] = None
        out_params["LAKE_DB_ID"] = None
        # TODO : None in filename if no lake DATABASES or basins file
        if "DATABASES" in config.sections():
            list_db = config.options("DATABASES")
            # Prior lake database
            if "lake_db" in list_db:
                out_params["LAKE_DB"] = config.get("DATABASES", "LAKE_DB")
            if "lake_db_id" in list_db:
                out_params["LAKE_DB_ID"] = config.get("DATABASES", "LAKE_DB_ID")

        # 4 - Retrieve OPTIONS
        if "OPTIONS" in config.sections():
            list_options = config.options("OPTIONS")
            # Flag to also produce LakeTile_edge and LakeTile_pixcvec as shapefiles (=True); else=False (default)
            if "produce shp" in list_options:
                out_params["flag_prod_shp"] = config.get("OPTIONS", "Produce shp")
            # Flag to increment the file counter in the output filenames (=True, default); else=False
            if "increment file counter" in list_options:
                out_params["flag_inc_file_counter"] = config.get("OPTIONS", "Increment file counter")

        # 5 - Retrieve LOGGING
        error_file = config.get("LOGGING", "errorFile")
        out_params["errorFile"] = os.path.splitext(error_file)[0] + "_" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + os.path.splitext(error_file)[1]
        log_file = config.get("LOGGING", "logFile")
        out_params["logFile"] = os.path.splitext(log_file)[0] + "_" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + os.path.splitext(log_file)[1]
        out_params["logfilelevel"] = config.get("LOGGING", "logfilelevel")
        out_params["logConsole"] = config.get("LOGGING", "logConsole")
        out_params["logconsolelevel"] = config.get("LOGGING", "logconsolelevel")
        
        # 6 - Retrieve FILE_INFORMATION
        # CRID
        out_params["CRID"] = config.get("FILE_INFORMATION", "CRID")
        # Producer
        out_params["PRODUCER"] = config.get("FILE_INFORMATION", "PRODUCER")
        # Method of production of the original data
        out_params["SOURCE"] = config.get("FILE_INFORMATION", "SOURCE")
        # Software version
        out_params["SOFTWARE_VERSION"] = config.get("FILE_INFORMATION", "SOFTWARE_VERSION")
        # Contact
        out_params["CONTACT"] = config.get("FILE_INFORMATION", "CONTACT")

        return out_params

    def _put_cmd_value(self, param_list):
        """
        Add command parameters to the class.
        All parameters added here come from the command file.
        This is usefull to centralize all processing parameters.
        
        :param param_list: dictionary with command parameters
        :type param_list: dict
        """

        try:
            
            # Add LOGGING section and parameters
            section = "LOGGING"
            self.cfg.add_section(section)
            self.cfg.set(section, "errorFile", param_list["errorFile"])
            self.cfg.set(section, "logFile", param_list["logFile"])
            self.cfg.set(section, "logfilelevel", param_list["logfilelevel"])
            self.cfg.set(section, "logConsole", param_list["logConsole"])
            self.cfg.set(section, "logconsolelevel", param_list["logconsolelevel"])
            
            # Add PATHS section and parameters
            section = "PATHS"
            self.cfg.add_section(section)
            self.cfg.set(section, "PIXC file", param_list["pixc_file"])
            self.cfg.set(section, "PIXCVecRiver file", param_list["pixc_vec_river_file"])
            self.cfg.set(section, "Output directory", param_list["output_dir"])
            
            # Add DATABASES section and parameters
            section = "DATABASES"
            self.cfg.add_section(section)
            self.cfg.set(section, "LAKE_DB", param_list["LAKE_DB"])
            self.cfg.set(section, "LAKE_DB_ID", param_list["LAKE_DB_ID"])

            # Add OPTIONS section and parameters
            section = "OPTIONS"
            self.cfg.add_section(section)
            self.cfg.set(section, "Produce shp", param_list["flag_prod_shp"])
            self.cfg.set(section, "Increment file counter", param_list["flag_inc_file_counter"])
            
            # Add section FILE_INFORMATION
            section = "FILE_INFORMATION"
            self.cfg.add_section(section)
            self.cfg.set(section, "CRID", param_list["CRID"])
            self.cfg.set(section, "PRODUCER", param_list["PRODUCER"])
            self.cfg.set(section, "SOURCE", param_list["SOURCE"])
            self.cfg.set(section, "SOFTWARE_VERSION", param_list["SOFTWARE_VERSION"])
            self.cfg.set(section, "CONTACT", param_list["CONTACT"])

        except Exception:
            print("Something wrong happened in _put_cmd_value !")
            raise

    def _check_config_parameters(self):
        """
        Check parameters coherence for LakeTile parameter file
        
        :return: True if OK
        :rtype: boolean
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        try:

            # 1 - Config parameters from command file
            
            # 1.1 - PATH section
            # PIXC file
            self.cfg.test_var_config_file('PATHS', 'PIXC file', str, logger=logger)
            my_tools.test_file(self.cfg.get('PATHS', 'PIXC file'))
            logger.debug('OK - PIXC file = ' + str(self.cfg.get('PATHS', 'PIXC file')))
            # PIXCVecRiver file
            self.cfg.test_var_config_file('PATHS', 'PIXCVecRiver file', str, logger=logger)
            my_tools.test_file(self.cfg.get('PATHS', 'PIXCVecRiver file'))
            logger.debug('OK - PIXCVecRiver file = ' + str(self.cfg.get('PATHS', 'PIXCVecRiver file')))
            # Output directory
            self.cfg.test_var_config_file('PATHS', 'Output directory', str, logger=logger)
            my_tools.test_dir(self.cfg.get('PATHS', 'Output directory'))
            logger.debug('OK - Output directory = ' + str(self.cfg.get('PATHS', 'Output directory')))

            # 1.2 - DATABASES section
            # Lake database full path
            if self.cfg.get('DATABASES', 'LAKE_DB') is None:
                logger.debug('WARNING - LAKE_DB not filled => LakeTile product not linked to a lake database')
            elif self.cfg.get('DATABASES', 'LAKE_DB').endswith(".shp"):
                self.cfg.test_var_config_file('DATABASES', 'LAKE_DB', str, logger=logger)
                my_tools.test_file(self.cfg.get('DATABASES', 'LAKE_DB'))
                logger.debug('OK - LAKE_DB = ' + str(self.cfg.get('DATABASES', 'LAKE_DB')))
                # Lake identifier attribute name in the database
                if self.cfg.get('DATABASES', 'LAKE_DB_ID'):
                    self.cfg.test_var_config_file('DATABASES', 'LAKE_DB_ID', str, logger=logger)
                    logger.debug('OK - LAKE_DB_ID = ' + str(self.cfg.get('DATABASES', 'LAKE_DB_ID')))
                else :
                    logger.warning('WARNING - LAKE_DB file given but the lake_id fieldname is missing')
            elif self.cfg.get('DATABASES', 'LAKE_DB').endswith(".sqlite"):
                self.cfg.test_var_config_file('DATABASES', 'LAKE_DB', str, logger=logger)
                my_tools.test_file(self.cfg.get('DATABASES', 'LAKE_DB'))
                logger.debug('OK - LAKE_DB = ' + str(self.cfg.get('DATABASES', 'LAKE_DB')))
            else :
                logger.debug('WARNING - Unknown LAKE_DB file format for file: %s => LakeTile product not linked to a lake database' % (self.cfg.get('DATABASES', 'LAKE_DB')))
                
            # 1.3 - OPTIONS section
            # Shapefile production
            self.cfg.test_var_config_file('OPTIONS', 'Produce shp', bool, logger=logger)
            logger.debug('OK - Produce shp = ' + str(self.cfg.get('OPTIONS', 'Produce shp')))
            # Increment output file counter
            self.cfg.test_var_config_file('OPTIONS', 'Increment file counter', bool, logger=logger)
            logger.debug('OK - Increment file counter = ' + str(self.cfg.get('OPTIONS', 'Increment file counter')))

            # 2 - Config parameters from parameter file

            # 2.1 - CONFIG_PARAMS section
            
            # Water flag = 3=water near land edge  4=interior water
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'FLAG_WATER', str, val_default="3;4", logger=logger)
            logger.debug('OK - FLAG_WATER = ' + str(self.cfg.get('CONFIG_PARAMS', 'FLAG_WATER')))
            # Dark water flag = 23=darkwater near land  24=interior dark water
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'FLAG_DARK', str, val_default="23;24", logger=logger)
            logger.debug('OK - FLAG_DARK = ' + str(self.cfg.get('CONFIG_PARAMS', 'FLAG_DARK')))
            
            # Min size for a lake to generate a lake product (=polygon + attributes) for it
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'MIN_SIZE', float, val_default=0.01, logger=logger)
            logger.debug('OK - MIN_SIZE = ' + str(self.cfg.get('CONFIG_PARAMS', 'MIN_SIZE')))
            # Maximal standard deviation of height inside a lake (-1 = do not compute lake height segmentation)
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'STD_HEIGHT_MAX', float, val_default=-1.0, logger=logger)
            logger.debug('OK - STD_HEIGHT_MAX = ' + str(self.cfg.get('CONFIG_PARAMS', 'STD_HEIGHT_MAX')))
            
            # To improve PixC golocation (=True) or not (=False)
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'IMP_GEOLOC', bool, val_default=True, logger=logger)
            logger.debug('OK - IMP_GEOLOC = ' + str(self.cfg.get('CONFIG_PARAMS', 'IMP_GEOLOC')))
            # Method to compute lake boundary or polygon hull
            # 0 = convex hull 
            # 1.0 = concave hull computed in ground geometry, based on Delaunay triangulation - using CGAL library
            # 1.1 = concave hull computed in ground geometry, based on Delaunay triangulation
            # 1.2 = concave hull computed in ground geometry, based on Delaunay triangulation - with alpha parameter varying across-track
            # 2 = edge computed in radar geometry, then converted in ground geometry (default)
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'HULL_METHOD', float, valeurs=[0, 1.0, 1.1, 1.2, 2], val_default=2, logger=logger)
            logger.debug('OK - HULL_METHOD = ' + str(self.cfg.get('CONFIG_PARAMS', 'HULL_METHOD')))
            # max number of pixel for hull computation 1            
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'NB_PIX_MAX_DELAUNEY', int, val_default=100000, logger=logger)
            logger.debug('OK - NB_PIX_MAX_DELAUNEY = ' + str(self.cfg.get('CONFIG_PARAMS', 'NB_PIX_MAX_DELAUNEY')))
            # max number of contour points for hull computation 2
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'NB_PIX_MAX_CONTOUR', int, val_default=8000, logger=logger)
            logger.debug('OK - NB_PIX_MAX_CONTOUR = ' + str(self.cfg.get('CONFIG_PARAMS', 'NB_PIX_MAX_CONTOUR')))

            # Big lakes parameters for improved geoloc
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'BIGLAKE_MODEL', str, valeurs=["polynomial", "no"], val_default="polynomial", logger=logger)
            logger.debug('OK - BIGLAKE_MODEL = ' + str(self.cfg.get('CONFIG_PARAMS', 'BIGLAKE_MODEL')))
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'BIGLAKE_MIN_SIZE', float, val_default=50.0, logger=logger)
            logger.debug('OK - BIGLAKE_MIN_SIZE = ' + str(self.cfg.get('CONFIG_PARAMS', 'BIGLAKE_MIN_SIZE')))
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'BIGLAKE_GRID_SPACING', float, val_default=4000, logger=logger)
            logger.debug('OK - BIGLAKE_GRID_SPACING = ' + str(self.cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_SPACING')))
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'BIGLAKE_GRID_RES', float, val_default=8000, logger=logger)
            logger.debug('OK - BIGLAKE_GRID_RES = ' + str(self.cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_RES')))
            
            # 2.2 - ID section
            # Nb digits for counter of lakes in a tile or pass
            self.cfg.test_var_config_file('ID', 'NB_DIGITS', str, val_default=6, logger=logger)
            logger.debug('OK - NB_DIGITS = ' + str(self.cfg.get('ID', 'NB_DIGITS')))
            
            # 2.3 - FILE_INFORMATION section
            # Composite Release IDentifier for LakeTile processing
            self.cfg.test_var_config_file('FILE_INFORMATION', 'CRID', str, logger=logger)
            logger.debug('OK - CRID = ' + str(self.cfg.get('FILE_INFORMATION', 'CRID')))
            # Producer
            self.cfg.test_var_config_file('FILE_INFORMATION', 'PRODUCER', str, logger=logger)
            logger.debug('OK - PRODUCER = ' + str(self.cfg.get('FILE_INFORMATION', 'PRODUCER')))
            # Method of production of the original data
            self.cfg.test_var_config_file('FILE_INFORMATION', 'SOURCE', str, logger=logger)
            logger.debug('OK - SOURCE = ' + str(self.cfg.get('FILE_INFORMATION', 'SOURCE')))
            # Software version
            self.cfg.test_var_config_file('FILE_INFORMATION', 'SOFTWARE_VERSION', str, logger=logger)
            logger.debug('OK - SOFTWARE_VERSION = ' + str(self.cfg.get('FILE_INFORMATION', 'SOFTWARE_VERSION')))
            # Contact
            self.cfg.test_var_config_file('FILE_INFORMATION', 'CONTACT', str, logger=logger)
            logger.debug('OK - CONTACT = ' + str(self.cfg.get('FILE_INFORMATION', 'CONTACT')))

        # Error managed
        except service_error.ConfigFileError:
            message = "Error in the configuration file " + self.cfg.path_conf
            logger.error(message, exc_info=True)
            raise
        # Warning error not managed !
        except Exception:
            logger.error("Something wrong happened during configuration file check!", exc_info=True)
            raise

        return True
        
    # -------------------------------------------

    def _read_input_data(self):
        """
        This method prepares and init input data classes
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Init PIXCVec object by retrieving data from the pixel cloud complementary file after river processing
        logger.info("> 1.1 - Init pixel cloud complementary file...")
        self.obj_pixc_vec = proc_pixc_vec.PixelCloudVec("TILE")
        self.obj_pixc_vec.set_from_pixcvec_file(self.pixc_vec_river_file)
        logger.info("")

        # 2 - Init PIXC object by retrieving needed data from the pixel cloud
        logger.info("> 1.2 - Retrieving needed data from the pixel cloud...")
        self.obj_pixc = proc_pixc.PixelCloud()
        self.obj_pixc.set_from_pixc_file(self.pixc_file, self.obj_pixc_vec.reject_index)
        logger.info("")

    def _read_lake_db(self):
        """
        This method prepares Lake DB class, and init with Lake DB data
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        lake_db_file = self.cfg.get("DATABASES", "LAKE_DB")
        
        if (lake_db_file == "") or (lake_db_file is None):  # No PLD
            logger.warning("NO database specified -> NO link of SWOT obs with a priori lake")
            self.obj_lake_db = lake_db.LakeDb()
            
        else:  # Init PLD object wrt to file type
            
            type_db = lake_db_file.split('.')[-1]  # File type
            
            if type_db == "shp":  # Shapefile format
                self.obj_lake_db = lake_db.LakeDbShp(lake_db_file, self.obj_pixc.tile_poly)
            elif type_db == "sqlite":  # SQLite format
                self.obj_lake_db = lake_db.LakeDbSqlite(lake_db_file, self.obj_pixc.tile_poly)
            elif os.path.isdir(lake_db_file):  # Directory containing Sqlite files
                self.obj_lake_db = lake_db.LakeDbDirectory(lake_db_file, self.obj_pixc.tile_poly)
            else:
                message = "Prior lake database format (%s) is unknown: must be .shp or .sqlite => set to None" % type_db
                logger.warning(message)
                self.obj_lake_db = lake_db.LakeDb()
        
    # -------------------------------------------

    def _prepare_output_data(self):
        """
        This method creates output data classes
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Retrieve orbit info from PIXC filename and compute output filenames
        logger.info("> 3.1 - Retrieving tile infos from PIXC filename...")
        self.lake_tile_filenames = locnes_filenames.LakeTileFilenames(self.pixc_file, self.pixc_vec_river_file, self.output_dir,
                                                                      flag_inc=self.flag_inc_file_counter)
        
        # 2 - Initialize LakeTile product object
        logger.info("> 3.2 - Initialize LakeTile product object...")
        self.obj_lake = proc_lake.LakeTileProduct(self.obj_pixc,
                                                  self.obj_pixc_vec,
                                                  self.obj_lake_db,
                                                  os.path.basename(self.lake_tile_filenames.lake_tile_shp_file_obs).split(".")[0])

    def _write_output_data(self):
        """
        This method write output data
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Write LakeTile shapefiles
        logger.info("> 8.1 - Writing LakeTile memory layer to shapefiles...")
        # 1.1 - Obs-oriented shapefile
        logger.info(". Obs-oriented shapefile = %s" % self.lake_tile_filenames.lake_tile_shp_file_obs)
        self.obj_lake.write_obs_file(self.lake_tile_filenames.lake_tile_shp_file_obs, self.proc_metadata)
        # 1.2 - PLD-oriented shapefile
        logger.info(". PLD-oriented shapefile = %s" % self.lake_tile_filenames.lake_tile_shp_file_prior)
        self.obj_lake.write_prior_file(self.lake_tile_filenames.lake_tile_shp_file_prior, self.proc_metadata)
        # 1.3 - Shapefile of unassigned water features
        logger.info(". Shapefile of unassigned water features = %s" % self.lake_tile_filenames.lake_tile_shp_file_unknown)
        self.obj_lake.write_unknown_file(self.lake_tile_filenames.lake_tile_shp_file_unknown, self.proc_metadata)
        # 1.4 - Close memory layer
        self.obj_lake.free_memory()
        logger.info("")
        
        # 2 - Write PIXCVec for objects entirely inside tile
        logger.info("> 8.2 - Writing LakeTile_PIXCVec file...")
        self.obj_pixc_vec.write_file(self.lake_tile_filenames.lake_tile_pixcvec_file, 
                                     self.proc_metadata, 
                                     in_pixc_metadata=self.obj_pixc.pixc_metadata)
        if self.flag_prod_shp:
            if self.obj_pixc_vec.nb_water_pix == 0:
                logger.info("NO water pixel => NO LakeTile_PIXCVec shapefile produced")
            else:
                self.obj_pixc_vec.write_file_as_shp(self.lake_tile_filenames.lake_tile_pixcvec_file.replace(".nc", ".shp"), self.obj_pixc)
        logger.info("")
        
        # 3 - Write intermediate NetCDF file with indices of pixels (and associated label) related to objects at the top/bottom edges of the tile
        logger.info("> 8.3 - Writing LakeTile_Edge file...")
        self.obj_pixc.write_edge_file(self.lake_tile_filenames.lake_tile_edge_file, self.proc_metadata)
        if self.flag_prod_shp:
            if self.obj_pixc.nb_edge_pix == 0:
                logger.info("NO edge pixel => NO LakeTile_edge shapefile produced")
            else:
                self.obj_pixc.write_edge_file_as_shp(self.lake_tile_filenames.lake_tile_edge_file.replace(".nc", ".shp"))
        logger.info("")
        
    # -------------------------------------------

    def start(self):
        """
        Main pseudoPGE method
        Start computation
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        try:
            
            logger.info("***************************************************")
            logger.info("***** Read input data and prepare output data *****")
            logger.info("***************************************************")

            # 1 - Read input data
            logger.info("> 1 - Init and format input objects...")
            self._read_input_data()
            logger.info("")
            
            # 2 - Read lake db
            logger.info("> 2 - Init and read lake DB object...")
            self._read_lake_db()
            logger.info("")

            # 3 - Prepare output data
            logger.info("> 3 - Prepare output objects...")
            self._prepare_output_data()
            logger.info("")
            
            logger.info("****************************")
            logger.info("***** Run SAS LakeTile *****")
            logger.info("****************************")
            
            # 4 - Initialization
            logger.info("> 4 - Initialization of SASLakeTile class")
            my_lake_tile = sas_lake_tile.SASLakeTile(self.obj_pixc, self.obj_pixc_vec, self.obj_lake_db, self.obj_lake)
            logger.info(self.timer.info(0)) 
            logger.info("")         
            
            # 5 - Run pre-processing
            logger.info("> 5 - Run SASpre-processing")
            my_lake_tile.run_preprocessing()
            logger.info(self.timer.info(0))
            logger.info("")
            
            # 6 - Run processing
            logger.info("> 6 - Run SASprocessing")
            my_lake_tile.run_processing()
            logger.info(self.timer.info(0))
            logger.info("")
    
            # 7 - Run post-processing
            logger.info("> 7 - Run SASpost-processing")
            my_lake_tile.run_postprocessing()
            logger.info(self.timer.info(0))
            logger.info("")
            
            logger.info("*****************************")
            logger.info("***** Write output data *****")
            logger.info("*****************************")
            
            # 8 - Write output data
            logger.info("> 8 - Write output data")
            self._write_output_data()
            
        except service_error.SwotError:
            raise
            
        except Exception:
            message = "Fatal error catch in PGE LakeTile"
            logger.error(message, exc_info=True)
            raise
        
    def stop(self):
        """
        pseudoPGE method stop
        Close all services
        """
        logger = logging.getLogger(self.__class__.__name__)
            
        logger.info("******************************")
        logger.info("***** Close all services *****")
        logger.info("******************************")
        logger.info("")
        
        # 1 - Close lake database
        logger.info("> 1 - Closing lake database...")
        self.obj_lake_db.close_db()

        logger.info("")
        logger.info(self.timer.stop())
        logger.sigmsg("====================================")
        logger.sigmsg("===== LakeTileProcessing = END =====")
        logger.sigmsg("====================================")

        # 2 - Stop logger
        instance_logger = service_logger.get_instance()
        if instance_logger is not None:
            instance_logger.flush_log()
            instance_logger.close()
            
        # 3 - Clean configuration file
        service_config_file.clear_config()
    

#######################################


if __name__ == '__main__':

    # 0 - Parse inline parameters
    PARSER = argparse.ArgumentParser(description="Compute SWOT LakeTile product\
            from one tile of PIXC product and its associated PIXCVecRiver product. ")
    PARSER.add_argument("command_file", help="command file (*.cfg)")
    ARGS = PARSER.parse_args()
    
    PGE = None
    
    try:
        # 1 - Instantiate PGE
        PGE = PGELakeTile(ARGS.command_file)

        # 2 - Start PGE LakeTile
        PGE.start()
        
    finally:
        if PGE is not None:
            # 3 - Stop PGE LakeTile
            PGE.stop()
