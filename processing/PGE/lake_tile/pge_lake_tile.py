#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
.. module:: pge_lake_tile.py
    :synopsis: Process PGE_L2_HR_LakeTile, i.e. generate L2_HR_LakeTile
    product from one tile of L2_HR_PIXC product and associated
    L2_HR_PIXCVec product
    Created on 02/27/2017

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

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

import cnes.common.lib.my_timer as my_timer
import cnes.common.lib.my_tools as my_tools
import cnes.sas.lake_tile.sas_lake_tile as sas_lake_tile
import cnes.common.lib_lake.locnes_filenames as locnes_filenames
import cnes.common.lib_lake.proc_pixc_vec as proc_pixc_vec
import cnes.sas.lake_tile.proc_pixc as proc_pixc
import cnes.common.lib_lake.lake_db as lake_db
import cnes.common.lib_lake.proc_lake as proc_lake
import cnes.common.lib.my_shp_file as my_shp

import cnes.common.service_config_file as service_config_file
import cnes.common.service_error as service_error
import cnes.common.service_logger as service_logger


#######################################


class PGELakeTile():
    """
    Class PGELakeTile
    Pseudo PGE class to launch SAS LakeTile computation
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
        my_tools.testFile(cmd_file, IN_extent=".cfg")  # Test existance and extension
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
        
        # 4 - Init data object
        self.obj_pixc_vec = None
        self.obj_pixc = None
        self.obj_lake_db = None
        self.obj_lake = None
        self.lake_tile_filenames = None

        self.flag_prod_shp = self.cfg.getboolean("OPTIONS", "Produce shp")

        # 4 - Initiate logging service
        service_logger.ServiceLogger()
        logger = logging.getLogger(self.__class__.__name__)

        # 5 - Print info
        logger.info("======================================")
        logger.info("===== lakeTileProcessing = BEGIN =====")
        logger.info("======================================")
        message = "> Command file: " + str(self.cmd_file)
        logger.info(message)
        message = "> " + str(self.cfg)
        logger.info(message)
        logger.info("")
        
        # 6 - Test input parameters
        logger.info(">> Test input parameters")
        self._check_config_parameters()
        logger.info("")

        # 7 - Update global variables
        # TODO replace all global variables by call to service_config_file
        # my_var.tmpGetConfigFromServiceConfigFile()
        
        # 8 - Form processing metadata dictionary
        self.proc_metadata = {}
        self.proc_metadata["xref_static_lake_db_file"] = self.cfg.get("DATABASES", "LAKE_DB")
        self.proc_metadata["xref_input_l2_hr_pixc_file"] = self.pixc_file
        self.proc_metadata["xref_input_l2_hr_pixc_vec_river_file"] = self.pixc_vec_river_file
        self.proc_metadata["xref_l2_hr_lake_tile_param_file"] = file_config
        
        logger.info("")
        logger.info("")

    def _test_input_data(self):
        """
        This method test mandatory input files
        """
        logger = logging.getLogger(self.__class__.__name__)
        # 1.1 - PIXC file
        message = "> 1.1  INPUT PIXC file = %s" % self.pixc_file
        logger.info(message)
        my_tools.testFile(self.pixc_file, IN_extent=".nc")
        # 1.2 - PIXCVecRiver file
        message = "> 1.2  INPUT PIXCVecRiver file = %s" % self.pixc_vec_river_file
        logger.info(message)
        my_tools.testFile(self.pixc_vec_river_file, IN_extent=".nc")

    def _test_output_directory(self):
        """
        This method test output directory
        """
        logger = logging.getLogger(self.__class__.__name__)
        message = "[lakeTileProcessing]   OUTPUT DIR = %s" % self.output_dir
        logger.info(message)
        my_tools.testDir(self.output_dir)

    def _read_input_data(self):
        """
        This method read input data and prepare data class
        """
        logger = logging.getLogger(self.__class__.__name__)
        # 2.1 - Init PIXCVec product by retrieving data from the pixel cloud complementary file after river processing
        logger.info("> 2.1 - Init pixel cloud complementary file...")
        #self.obj_pixc_vec = proc_pixc_vec.PixelCloudVec("TILE", self.pixc_vec_river_file)
        self.obj_pixc_vec = proc_pixc_vec.PixelCloudVec("TILE")
        self.obj_pixc_vec.set_from_pixcvec_file(self.pixc_vec_river_file)
        logger.info("")

        # 2.2 - Init Pixc product by retrieving needed data from the pixel cloud
        logger.info("> 2.2 - Retrieving needed data from the pixel cloud...")
        #self.obj_pixc = proc_pixc.PixelCloud(self.pixc_file, self.obj_pixc_vec.reject_idx)
        self.obj_pixc = proc_pixc.PixelCloud()
        self.obj_pixc.set_from_pixc_file(self.pixc_file, self.obj_pixc_vec.reject_index)
        logger.info("")

    def _read_lake_db(self):
        """
        This method create output data class
        """
        logger = logging.getLogger(self.__class__.__name__)
        # 3 - Retrieve lake Db layer
        lake_db_file = self.cfg.get("DATABASES", "LAKE_DB")
        if (lake_db_file == "") or (lake_db_file is None):
            logger.warning("NO database specified -> NO link of SWOT obs with a priori lake")
        else:
            if os.path.exists(lake_db_file):
                type_db = lake_db_file.split('.')[-1]  # Type of database
                if type_db == "shp":  # Shapefile format
                    self.obj_lake_db = lake_db.LakeDbShp(lake_db_file, self.obj_pixc.tile_poly)
                elif type_db == "sqlite":  # SQLite format
                    self.obj_lake_db = lake_db.LakeDbSqlite(lake_db_file, self.obj_pixc.tile_poly)
                else:
                    message = "Lake a priori database format (%s) is unknown: must be .shp or .sqlite" % type_db
                    raise service_error.ProcessingError(message, logger)
            else:
                message = "ERROR = %s doesn't exist" % lake_db_file
                raise service_error.ProcessingError(message, logger)

    def _prepare_output_data(self):
        """
        This method create output data class
        """
        # 4 - Retrieve orbit info from PIXC filename and compute output filenames
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("> 4.1 - Retrieving tile infos from PIXC filename...")
        self.lake_tile_filenames = locnes_filenames.lakeTileFilenames(self.pixc_file, self.pixc_vec_river_file, self.output_dir)
        # 4 - Initialize lake product
        logger.info("> 4.2 - Initialize lake product object...")
        self.obj_lake = proc_lake.LakeProduct("TILE",
                                              self.obj_pixc,
                                              self.obj_pixc_vec,
                                              self.obj_lake_db,
                                              os.path.basename(self.lake_tile_filenames.lake_tile_shp_file).split(".")[0],
                                              in_id_prefix=self.lake_tile_filenames.lake_id_prefix)


    def _write_output_data(self, in_proc_metadata):
        """
        This method write output data
        """
        logger = logging.getLogger(self.__class__.__name__)
        # 1 - Write LakeTile shapefile
        logger.info("9.1 - Writing LakeTile memory layer to shapefile...")
        my_shp.write_mem_layer_as_shp(self.obj_lake.shp_mem_layer.layer, self.lake_tile_filenames.lake_tile_shp_file)
        self.obj_lake.shp_mem_layer.free()  # Close memory layer
        logger.info("")
        # Write XML metadatafile for shapefile
        #self.obj_lake.writeMetadataFile("%s.xml" % self.lake_tile_filenames.lake_tile_shp_file)
        # Write XML metadatafile for shapefile
        self.obj_lake.shp_mem_layer.update_and_write_metadata("%s.xml" % self.lake_tile_filenames.lake_tile_shp_file, 
                                                              in_pixc_metadata=self.obj_pixc.pixc_metadata,
                                                              in_proc_metadata=in_proc_metadata)
        # 2 - Write PIXCVec for objects entirely inside tile
        logger.info("9.2 - Writing LakeTile_pixcvec file...")
        self.obj_pixc_vec.write_file(self.lake_tile_filenames.lake_tile_pixcvec_file, in_proc_metadata)
        if self.flag_prod_shp and (self.obj_pixc_vec.nb_water_pix != 0):
            self.obj_pixc_vec.write_file_asShp(self.lake_tile_filenames.lake_tile_pixcvec_file.replace(".nc", ".shp"), self.obj_pixc)
        logger.info("")

        # 3 - Write intermediate NetCDF file with indices of pixels (and associated label) related to objects at the top/bottom edges of the tile
        logger.info("9.3 - Writing LakeTile_edge file...")
        self.obj_pixc.write_edge_file(self.lake_tile_filenames.lake_tile_edge_file, in_proc_metadata)
        if self.flag_prod_shp and (self.obj_pixc.nb_edge_pix != 0):
            self.obj_pixc.write_edge_file_asShp(self.lake_tile_filenames.lake_tile_edge_file.replace(".nc", ".shp"))
        logger.info("")


    def start(self):
        """
        Main pseudoPGE method
        Start computation
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        try:
            # 1 - Test existance and file format of input paths
            logger.info("> 1 - Testing existence of input paths...")
            self._test_input_data()
            self._test_output_directory()

            # 2 - Read input data
            logger.info("> 2 - Init and format intput objects...")
            self._read_input_data()
            
            # 3 - Read lake db
            logger.info("> 3 - Init and read lake db object...")
            self._read_lake_db()

            # 4 - Prepare output data
            logger.info("> 4 - Prepare output object...")
            self._prepare_output_data()

            # 5 - Initialization
            logger.info("> 5 - Initialization of SASLakeTile class")
            my_lake_tile = sas_lake_tile.SASLakeTile(self.obj_pixc, self.obj_pixc_vec, self.obj_lake_db, self.obj_lake)
            logger.info(self.timer.info(0))          
            
            # 6 - Run pre-processing
            logger.info("> 6 - Run SASpre-processing")
            my_lake_tile.run_preprocessing()
            logger.info(self.timer.info(0))
            
            # 7 - Run processing
            logger.info("> 7 - Run SASprocessing")
            my_lake_tile.run_processing()
            logger.info(self.timer.info(0))
            logger.info("")
    
            # 8 - Run post-processing
            logger.info("> 8 - Run SASpost-processing")
            my_lake_tile.run_postprocessing()
            logger.info(self.timer.info(0))
            logger.info("")
            
            # 9 - Write output data
            logger.info("> 9 - Write output data")
            self._write_output_data(self.proc_metadata)
            
        except service_error.SwotError:
            raise
            
        except Exception:
            message = "Fatal error catch in PGE lake_tile"
            logger.error(message, exc_info=True)
            raise
        
    def stop(self):
        """
        pseudoPGE method stop
        Close all services
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Close lake database
        if self.obj_lake_db is not None:
            logger.info("Closing lake database...")
            self.obj_lake_db.close_db()

        logger.info("")
        logger.info(self.timer.stop())
        logger.info("====================================")
        logger.info("===== lakeTileProcessing = END =====")
        logger.info("====================================")

        # 2 - Stop logger
        instance_logger = service_logger.get_instance()
        if instance_logger is not None:
            instance_logger.flush_log()
            instance_logger.close()
            
        # 3 - Clean configuration file
        service_config_file.clear_config()
        
    # -------------------------------------------

    def _read_cmd_file(self):
        """
            Read the parameter file in input and store parameters in a dictionary

            :param filename: parameter file full path
            :type filename: string

            :return: dictionary containing parameters
            :rtype: dict
        """
        
        # 0 - Init output dictionary
        out_params = {}
        # Default values
        out_params["CONTINENT_FILE"] = None
        out_params["Produce_shp"] = False

        # 1 - Read parameter file
        config = configparser.ConfigParser()
        config.read(self.cmd_file)

        # 2 - Retrieve PATHS
        try:
            out_params["param_file"] = os.path.expandvars(config.get("PATHS", "param_file"))
        except:
            out_params["param_file"] = "lake_tile_param.cfg"
        out_params["pixc_file"] = config.get("PATHS", "PIXC file")
        out_params["pixc_vec_river_file"] = config.get("PATHS", "PIXCVecRiver file")
        out_params["output_dir"] = config.get("PATHS", "Output directory")

        # 3 - Retrieve DATABASES
        out_params["LAKE_DB"] = None
        out_params["LAKE_DB_ID"] = None
        out_params["CONTINENT_FILE"] = None
        # TODO : None in filename if no lake DATABASES or continent file
        if "DATABASES" in config.sections():
            list_db = config.options("DATABASES")
            # Lake a priori database
            if "lake_db" in list_db:
                out_params["LAKE_DB"] = config.get("DATABASES", "LAKE_DB")
            if "lake_db_id" in list_db:
                out_params["LAKE_DB_ID"] = config.get("DATABASES", "LAKE_DB_ID")
            # Continent file
            if "continent_file" in list_db:
                out_params["CONTINENT_FILE"] = config.get("DATABASES", "CONTINENT_FILE")

        # 4 - Retrieve OPTIONS
        if "OPTIONS" in config.sections():
            list_options = config.options("OPTIONS")
            # Flag to also produce LakeTile_edge and LakeTile_pixcvec as shapefiles (=True); else=False (default)
            if "produce shp" in list_options:
                out_params["Produce_shp"] = config.get("OPTIONS", "Produce shp")

        # 5 - Retrieve LOGGING
        log_file = config.get("LOGGING", "logFile")
        out_params["logFile"] = os.path.splitext(log_file)[0] + "_" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + os.path.splitext(log_file)[1]
        out_params["logfilelevel"] = config.get("LOGGING", "logfilelevel")
        out_params["logConsole"] = config.get("LOGGING", "logConsole")
        out_params["logconsolelevel"] = config.get("LOGGING", "logconsolelevel")
        
        # 6 - Retrieve FILE informations
        out_params["PRODUCER"] = config.get("FILE_INFORMATION", "PRODUCER")
        out_params["LAKE_TILE_CRID"] = config.get("CRID", "LAKE_TILE_CRID")
        # PIXC
        out_params["PIXC_PREFIX"] = config.get("FILE_INFORMATION", "PIXC_PREFIX")
        out_params["PIXC_PATTERN_PRINT"] = config.get("FILE_INFORMATION", "PIXC_PATTERN_PRINT")
        out_params["PIXC_PATTERN_IND"] = config.get("FILE_INFORMATION", "PIXC_PATTERN_IND")
        # PIXCVEC_RIVER
        out_params["PIXCVEC_RIVER_PREFIX"] = config.get("FILE_INFORMATION", "PIXCVEC_RIVER_PREFIX")
        out_params["PIXCVEC_RIVER_PATTERN_PRINT"] = config.get("FILE_INFORMATION", "PIXCVEC_RIVER_PATTERN_PRINT")
        out_params["PIXCVEC_RIVER_PATTERN_IND"] = config.get("FILE_INFORMATION", "PIXCVEC_RIVER_PATTERN_IND")
        # LAKE_TILE
        out_params["LAKE_TILE_PREFIX"] = config.get("FILE_INFORMATION", "LAKE_TILE_PREFIX")
        # ATTENTION % pas supporté dans config parser
        #out_params["LAKE_TILE_PATTERN"] = config.get("FILE_INFORMATION", "LAKE_TILE_PATTERN")
        out_params["LAKE_TILE_PATTERN"] = out_params["LAKE_TILE_PREFIX"] + "%03d_%03d_%s_%s_%s_%s_%02d%s"  # LakeTile filename with %03d=cycle number %03d=pass number %s=tile ref %s=swath %s=begin date %s=end date %s=CRID %s=counter %s=suffix
        # ATTENTION % pas supporté dans config parser
        #out_params["LAKE_TILE_PATTERN_PRINT"] = config.get("FILE_INFORMATION", "LAKE_TILE_PATTERN_PRINT")
        out_params["LAKE_TILE_PATTERN_PRINT"] = out_params["LAKE_TILE_PREFIX"] + "%s<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>"
        out_params["LAKE_TILE_PATTERN_IND"] = config.get("FILE_INFORMATION", "LAKE_TILE_PATTERN_IND")
        out_params["LAKE_TILE_SHP_SUFFIX"] = config.get("FILE_INFORMATION", "LAKE_TILE_SHP_SUFFIX")
        out_params["LAKE_TILE_SHP_META_SUFFIX"] = config.get("FILE_INFORMATION", "LAKE_TILE_SHP_META_SUFFIX")
        out_params["LAKE_TILE_EDGE_SUFFIX"] = config.get("FILE_INFORMATION", "LAKE_TILE_EDGE_SUFFIX")
        out_params["LAKE_TILE_PIXCVEC_SUFFIX"] = config.get("FILE_INFORMATION", "LAKE_TILE_PIXCVEC_SUFFIX")
      
        # PIXCVEC
        out_params["PIXCVEC_PREFIX"] = config.get("FILE_INFORMATION", "PIXCVEC_PREFIX")
        out_params["PIXCVEC_SUFFIX"] = config.get("FILE_INFORMATION", "PIXCVEC_SUFFIX")
        # ATTENTION % pas supporté dans config parser
        #out_params["PIXCVEC_PATTERN"] = config.get("FILE_INFORMATION", "PIXCVEC_PATTERN")
        out_params["PIXCVEC_PATTERN"] = out_params["PIXCVEC_PREFIX"] + "%03d_%03d_%s_%s_%s_%s_%02d" + out_params["PIXCVEC_SUFFIX"]  # PIXCVec filename with %03d=cycle number %03d=pass number %s=tile ref %s=begin date %s=end date %s=CRID %s=counter 
        # LAKE_SP
        out_params["LAKE_SP_CRID"] = config.get("CRID", "LAKE_SP_CRID")
        out_params["LAKE_SP_PREFIX"] = config.get("FILE_INFORMATION", "LAKE_SP_PREFIX")
        # ATTENTION % pas supporté dans config parser
        #out_params["LAKE_SP_PATTERN"] = config.get("FILE_INFORMATION", "LAKE_SP_PATTERN")
        out_params["LAKE_SP_PATTERN"]  = out_params["LAKE_SP_PREFIX"] + "%03d_%03d_%s_%s_%s_%s_%02d.shp"  # LakeSP filename with %03d=cycle number %03d=pass number %s=continent %s=begin date %s=end date %s=CRID %s=counter
        out_params["LAKE_SP_PATTERN_NO_CONT"] = out_params["LAKE_SP_PREFIX"] + "%03d_%03d_%s_%s_%s_%02d.shp"  # LakeSP filename without continent info with %03d=cycle number %03d=pass number %s=begin date %s=end date %s=CRID %s=counter

        return out_params

    def _put_cmd_value(self, param_list):
        """
            Add command parameters to the class.
            All parameters added here come from the command file.
            This is usefull to centralize all processing parameters.
        
            :param param_list: dictionary with command parameters
        """

        try:
            
            # Add LOGGING section and parameters
            section = "LOGGING"
            self.cfg.add_section(section)
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
            print(param_list["LAKE_DB"])
            self.cfg.set(section, "LAKE_DB", param_list["LAKE_DB"])
            self.cfg.set(section, "LAKE_DB_ID", param_list["LAKE_DB_ID"])
            self.cfg.set(section, "CONTINENT_FILE", param_list["CONTINENT_FILE"])
            
            # Add OPTIONS section and parameters
            section = "OPTIONS"
            self.cfg.add_section(section)
            self.cfg.set(section, "Produce shp", param_list["Produce_shp"])

            # add section CRID
            section = "CRID"
            self.cfg.add_section(section)
            self.cfg.set(section, "LAKE_TILE_CRID", param_list["LAKE_TILE_CRID"])
            self.cfg.set(section, "LAKE_SP_CRID", param_list["LAKE_SP_CRID"])
            
            # add section FILE_INFORMATION
            section = "FILE_INFORMATION"
            if section not in self.cfg.sections():
                self.cfg.add_section(section)
                self.cfg.set(section, "PRODUCER", param_list["PRODUCER"])
                self.cfg.set(section, "PIXC_PREFIX", param_list["PIXC_PREFIX"])
                self.cfg.set(section, "PIXC_PATTERN_PRINT", param_list["PIXC_PATTERN_PRINT"])
                self.cfg.set(section, "PIXC_PATTERN_IND", param_list["PIXC_PATTERN_IND"])
                self.cfg.set(section, "PIXCVEC_RIVER_PREFIX", param_list["PIXCVEC_RIVER_PREFIX"])
                self.cfg.set(section, "PIXCVEC_RIVER_PATTERN_PRINT", param_list["PIXCVEC_RIVER_PATTERN_PRINT"])
                self.cfg.set(section, "PIXCVEC_RIVER_PATTERN_IND", param_list["PIXCVEC_RIVER_PATTERN_IND"])
                self.cfg.set(section, "LAKE_TILE_PREFIX", param_list["LAKE_TILE_PREFIX"])
                self.cfg.set(section, "LAKE_TILE_PATTERN", param_list["LAKE_TILE_PATTERN"])
                self.cfg.set(section, "LAKE_TILE_PATTERN_PRINT", param_list["LAKE_TILE_PATTERN_PRINT"])
                self.cfg.set(section, "LAKE_TILE_PATTERN_IND", param_list["LAKE_TILE_PATTERN_IND"])
                self.cfg.set(section, "LAKE_TILE_SHP_SUFFIX", param_list["LAKE_TILE_SHP_SUFFIX"])
                self.cfg.set(section, "LAKE_TILE_SHP_META_SUFFIX", param_list["LAKE_TILE_SHP_META_SUFFIX"])
                self.cfg.set(section, "LAKE_TILE_EDGE_SUFFIX", param_list["LAKE_TILE_EDGE_SUFFIX"])
                self.cfg.set(section, "LAKE_TILE_PIXCVEC_SUFFIX", param_list["LAKE_TILE_PIXCVEC_SUFFIX"])
                self.cfg.set(section, "PIXCVEC_PREFIX", param_list["PIXCVEC_PREFIX"])
                self.cfg.set(section, "PIXCVEC_SUFFIX", param_list["PIXCVEC_SUFFIX"])
                self.cfg.set(section, "PIXCVEC_PATTERN", param_list["PIXCVEC_PATTERN"])
                self.cfg.set(section, "LAKE_SP_PREFIX", param_list["LAKE_SP_PREFIX"])
                self.cfg.set(section, "LAKE_SP_PATTERN", param_list["LAKE_SP_PATTERN"])
                self.cfg.set(section, "LAKE_TILE_SHP_SUFFIX", param_list["LAKE_TILE_SHP_SUFFIX"])
            
        except Exception:
            print("Something wrong happened in _put_cmd_value !")
            raise

    def _check_config_parameters(self):
        """
        Check parameters coherence for LakeTile parameter file
        
        :return: True if OK
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        try:

            # 1 - Config parameters from command file
            
            # 1.1 - PATH section
            # PIXC file
            self.cfg.test_var_config_file('PATHS', 'PIXC file', str)
            my_tools.testFile(self.cfg.get('PATHS', 'PIXC file'))
            logger.debug('PIXC file = ' + str(self.cfg.get('PATHS', 'PIXC file')))
            # PIXCVecRiver file
            self.cfg.test_var_config_file('PATHS', 'PIXCVecRiver file', str)
            my_tools.testFile(self.cfg.get('PATHS', 'PIXCVecRiver file'))
            logger.debug('PIXCvecRiver file = ' + str(self.cfg.get('PATHS', 'PIXCVecRiver file')))
            # Output directory
            self.cfg.test_var_config_file('PATHS', 'Output directory', str)
            my_tools.testDir(self.cfg.get('PATHS', 'Output directory'))
            logger.debug('Output directory = ' + str(self.cfg.get('PATHS', 'Output directory')))

            # 1.2 - DATABASES section
            # Lake database full path
            self.cfg.test_var_config_file('DATABASES', 'LAKE_DB', str)
            my_tools.testFile(self.cfg.get('DATABASES', 'LAKE_DB'))
            logger.debug('LAKE_DB = ' + str(self.cfg.get('DATABASES', 'LAKE_DB')))
            # Lake identifier attribute name in the database
            self.cfg.test_var_config_file('DATABASES', 'LAKE_DB_ID', str)
            logger.debug('LAKE_DB_ID = ' + str(self.cfg.get('DATABASES', 'LAKE_DB_ID')))
            # Continent file if want LakeSP product split per continent
            if self.cfg.get('DATABASES', 'CONTINENT_FILE') is None:
                logger.debug('CONTINENT_FILE not filled => LakeTile product not linked to a continent')
            else:
                self.cfg.test_var_config_file('DATABASES', 'CONTINENT_FILE', str)
                my_tools.testFile(self.cfg.get('DATABASES', 'CONTINENT_FILE'))
                logger.debug('CONTINENT_FILE = ' + str(self.cfg.get('DATABASES', 'CONTINENT_FILE')))

            # 1.3 - OPTIONS section
            # Shapefile production
            self.cfg.test_var_config_file('OPTIONS', 'Produce shp', bool)
            logger.debug('Produce shp = ' + str(self.cfg.get('OPTIONS', 'Produce shp')))

            # 2 - Config parameters from parameter file

            # 2.1 - CONFIG_PARAMS section
            
            # Water flag = 3=water near land edge  4=interior water
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'FLAG_WATER', str, val_default="3;4", logger=logger)
            logger.debug('FLAG_WATER = ' + str(self.cfg.get('CONFIG_PARAMS', 'FLAG_WATER')))
            # Dark water flag = 23=darkwater near land  24=interior dark water
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'FLAG_DARK', str, val_default="23;24", logger=logger)
            logger.debug('FLAG_DARK = ' + str(self.cfg.get('CONFIG_PARAMS', 'FLAG_DARK')))
            
            # Min size for a lake to generate a lake product (=polygon + attributes) for it
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'MIN_SIZE', float, val_default=1.0, logger=logger)
            logger.debug('MIN_SIZE = ' + str(self.cfg.get('CONFIG_PARAMS', 'MIN_SIZE')))
            # Maximal standard deviation of height inside a lake
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'STD_HEIGHT_MAX', float, val_default=10)
            logger.debug('STD_HEIGHT_MAX = ' + str(self.cfg.get('CONFIG_PARAMS', 'STD_HEIGHT_MAX')))
            
            # To improve PixC golocation (=True) or not (=False)
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'IMP_GEOLOC', bool, val_default=True, logger=logger)
            logger.debug('IMP_GEOLOC = ' + str(self.cfg.get('CONFIG_PARAMS', 'IMP_GEOLOC')))
            # Method to compute lake boundary or polygon hull
            # 0 = convex hull 
            # 1.0 = concave hull computed in ground geometry, based on Delaunay triangulation - using CGAL library (default) 
            # 1.1 = concave hull computed in ground geometry, based on Delaunay triangulation - with alpha parameter varying across-track
            # 1.1 = concave hull computed in ground geometry, based on Delaunay triangulation - without alpha parameter varying across-track
            # 2 = edge computed in radar geometry, then converted in ground geometry
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'HULL_METHOD', float, valeurs=[0, 1.0, 1.1, 1.2, 2], val_default=1.0, logger=logger)
            logger.debug('HULL_METHOD = ' + str(self.cfg.get('CONFIG_PARAMS', 'HULL_METHOD')))
            # max number of pixel for hull computation 1            
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'NB_PIX_MAX_DELAUNEY', int, val_default=100000, logger=logger)
            logger.debug('NB_PIX_MAX_DELAUNEY = ' + str(self.cfg.get('CONFIG_PARAMS', 'NB_PIX_MAX_DELAUNEY')))
            # max number of contour points for hull computation 2
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'NB_PIX_MAX_CONTOUR', int, val_default=8000, logger=logger)
            logger.debug('NB_PIX_MAX_CONTOUR = ' + str(self.cfg.get('CONFIG_PARAMS', 'NB_PIX_MAX_CONTOUR')))

            # Big lakes parameters for improved geoloc
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'BIGLAKE_MODEL', str, valeurs=["polynomial", "no"], val_default="polynomial", logger=logger)
            logger.debug('BIGLAKE_MODEL = ' + str(self.cfg.get('CONFIG_PARAMS', 'BIGLAKE_MODEL')))
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'BIGLAKE_MIN_SIZE', float, val_default=50000000.0, logger=logger)
            logger.debug('BIGLAKE_MIN_SIZE = ' + str(self.cfg.get('CONFIG_PARAMS', 'BIGLAKE_MIN_SIZE')))
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'BIGLAKE_GRID_SPACING', float, val_default=4000, logger=logger)
            logger.debug('BIGLAKE_GRID_SPACING = ' + str(self.cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_SPACING')))
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'BIGLAKE_GRID_RES', float, val_default=8000, logger=logger)
            logger.debug('BIGLAKE_GRID_RES = ' + str(self.cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_RES')))
            
            # 2.2 - ID section
            # Nb digits for counter of lakes in a tile or pass
            self.cfg.test_var_config_file('ID', 'NB_DIGITS', str, val_default=4, logger=logger)
            logger.debug('NB_DIGITS = ' + str(self.cfg.get('ID', 'NB_DIGITS')))
            
            # 2.3 - FILE_INFORMATION section
            
            # Product generator
            self.cfg.test_var_config_file('FILE_INFORMATION', 'PRODUCER', str)
            logger.debug('PRODUCER = ' + str(self.cfg.get('FILE_INFORMATION', 'PRODUCER')))
            # Composite Release IDentifier for LakeTile processing
            self.cfg.test_var_config_file('CRID', 'LAKE_TILE_CRID', str)
            logger.debug('LAKE_TILE_CRID = ' + str(self.cfg.get('CRID', 'LAKE_TILE_CRID')))
            # Composite Release IDentifier for LakeSP processing
            self.cfg.test_var_config_file('CRID', 'LAKE_SP_CRID', str)
            logger.debug('LAKE_SP_CRID = ' + str(self.cfg.get('CRID', 'LAKE_SP_CRID')))
            
            # PIXC product
            self.cfg.test_var_config_file('FILE_INFORMATION', 'PIXC_PREFIX', str)
            logger.debug('PIXC_PREFIX = ' + str(self.cfg.get('FILE_INFORMATION', 'PIXC_PREFIX')))
            self.cfg.test_var_config_file('FILE_INFORMATION', 'PIXC_PATTERN_PRINT', str)
            logger.debug('PIXC_PATTERN_PRINT = ' + str(self.cfg.get('FILE_INFORMATION', 'PIXC_PATTERN_PRINT')))
            self.cfg.test_var_config_file('FILE_INFORMATION', 'PIXC_PATTERN_IND', str)
            logger.debug('PIXC_PATTERN_IND = ' + str(self.cfg.get('FILE_INFORMATION', 'PIXC_PATTERN_IND')))
            
            # PIXCVecRiver product
            self.cfg.test_var_config_file('FILE_INFORMATION', 'PIXCVEC_RIVER_PREFIX', str)
            logger.debug('PIXCVEC_RIVER_PREFIX = ' + str(self.cfg.get('FILE_INFORMATION', 'PIXCVEC_RIVER_PREFIX')))
            self.cfg.test_var_config_file('FILE_INFORMATION', 'PIXCVEC_RIVER_PATTERN_PRINT', str)
            logger.debug('PIXCVEC_RIVER_PATTERN_PRINT = ' + str(self.cfg.get('FILE_INFORMATION', 'PIXCVEC_RIVER_PATTERN_PRINT')))
            self.cfg.test_var_config_file('FILE_INFORMATION', 'PIXCVEC_RIVER_PATTERN_IND', str)
            logger.debug('PIXCVEC_RIVER_PATTERN_IND = ' + str(self.cfg.get('FILE_INFORMATION', 'PIXCVEC_RIVER_PATTERN_IND')))
            
            # LakeTile product
            self.cfg.test_var_config_file('FILE_INFORMATION', 'LAKE_TILE_PREFIX', str)
            logger.debug('LAKE_TILE_PREFIX = ' + str(self.cfg.get('FILE_INFORMATION', 'LAKE_TILE_PREFIX')))
            self.cfg.test_var_config_file('FILE_INFORMATION', 'LAKE_TILE_PATTERN', str)
            logger.debug('LAKE_TILE_PATTERN = ' + str(self.cfg.get('FILE_INFORMATION', 'LAKE_TILE_PATTERN')))
            self.cfg.test_var_config_file('FILE_INFORMATION', 'LAKE_TILE_PATTERN_PRINT', str)
            logger.debug('LAKE_TILE_PATTERN_PRINT = ' + str(self.cfg.get('FILE_INFORMATION', 'LAKE_TILE_PATTERN_PRINT')))
            self.cfg.test_var_config_file('FILE_INFORMATION', 'LAKE_TILE_PATTERN_IND', str)
            logger.debug('LAKE_TILE_PATTERN_IND = ' + str(self.cfg.get('FILE_INFORMATION', 'LAKE_TILE_PATTERN_IND')))
            self.cfg.test_var_config_file('FILE_INFORMATION', 'LAKE_TILE_SHP_SUFFIX', str)
            logger.debug('LAKE_TILE_SHP_SUFFIX = ' + str(self.cfg.get('FILE_INFORMATION', 'LAKE_TILE_SHP_SUFFIX')))
            self.cfg.test_var_config_file('FILE_INFORMATION', 'LAKE_TILE_SHP_META_SUFFIX', str)
            logger.debug('LAKE_TILE_SHP_META_SUFFIX = ' + str(self.cfg.get('FILE_INFORMATION', 'LAKE_TILE_SHP_META_SUFFIX')))
            self.cfg.test_var_config_file('FILE_INFORMATION', 'LAKE_TILE_EDGE_SUFFIX', str)
            logger.debug('LAKE_TILE_EDGE_SUFFIX = ' + str(self.cfg.get('FILE_INFORMATION', 'LAKE_TILE_EDGE_SUFFIX')))
            self.cfg.test_var_config_file('FILE_INFORMATION', 'LAKE_TILE_PIXCVEC_SUFFIX', str)
            logger.debug('LAKE_TILE_PIXCVEC_SUFFIX = ' + str(self.cfg.get('FILE_INFORMATION', 'LAKE_TILE_PIXCVEC_SUFFIX')))
            
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

        # 2 - Start PGE Lake Tile
        PGE.start()
    finally:
        if PGE is not None:
            # 3 - Stop PGE Lake Tile
            PGE.stop()
