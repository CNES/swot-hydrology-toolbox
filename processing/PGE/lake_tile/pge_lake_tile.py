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
import cnes.common.lib_lake.locnes_variables as my_var
import cnes.sas.lake_tile.sas_lake_tile as sas_lake_tile

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
        print(my_params)
        # 2.2 - Test existence
        if not os.path.exists(file_config):
            raise service_error.DirFileError(file_config)
        # 2.3 - Load parameters
        self.cfg = service_config_file.ServiceConfigFile(file_config)

        # 3 - Put command parameter inside cfg
        self._put_cmd_value(my_params)
        
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
        my_var.tmpGetConfigFromServiceConfigFile()
        
        logger.info("")
        logger.info("")

    def start(self):
        """
        Main pseudoPGE method
        Start computation
        """
        logger = logging.getLogger(self.__class__.__name__)

        # 1 - Initialization
        shp_option = self.cfg.getboolean("OPTIONS", "Produce shp")
        my_lake_tile = sas_lake_tile.SASLakeTile(self.pixc_file, self.pixc_vec_river_file, self.output_dir, IN_shp_option=shp_option)
        logger.info(self.timer.info(0))
        logger.info("")

        # 3 - Run pre-processing
        my_lake_tile.run_preprocessing()
        logger.info(self.timer.info(0))
        
        # 4 - Run processing
        my_lake_tile.run_processing()
        logger.info(self.timer.info(0))
        logger.info("")

        # 5 - Run post-processing
        my_lake_tile.run_postprocessing()
        logger.info(self.timer.info(0))
        logger.info("")
        
    def stop(self):
        """
        pseudoPGE method stop
        Close all services
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        logger.info("")
        logger.info(self.timer.stop())
        logger.info("====================================")
        logger.info("===== lakeTileProcessing = END =====")
        logger.info("====================================")

        # Stop logger
        instance_logger = service_logger.get_instance()
        if instance_logger is not None:
            instance_logger.flush_log()
            instance_logger.close()
            
        # Clean configuration file
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
            self.cfg.set(section, "LAKE_DB", param_list["LAKE_DB"])
            self.cfg.set(section, "LAKE_DB_ID", param_list["LAKE_DB_ID"])
            self.cfg.set(section, "CONTINENT_FILE", param_list["CONTINENT_FILE"])
            
            # Add OPTIONS section and parameters
            section = "OPTIONS"
            self.cfg.add_section(section)
            self.cfg.set(section, "Produce shp", param_list["Produce_shp"])

        except Exception:
            print("Something wrong happened in ServiceConfigFile !")
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
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'FLAG_WATER', str)
            logger.debug('FLAG_WATER = ' + str(self.cfg.get('CONFIG_PARAMS', 'FLAG_WATER')))
            # Dark water flag = 23=darkwater near land  24=interior dark water
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'FLAG_DARK', str)
            logger.debug('FLAG_DARK = ' + str(self.cfg.get('CONFIG_PARAMS', 'FLAG_DARK')))
            # Layover flag
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'FLAG_LAYOVER', str)
            logger.debug('FLAG_LAYOVER = ' + str(self.cfg.get('CONFIG_PARAMS', 'FLAG_LAYOVER')))
            
            # Min size for a lake to generate a lake product (=polygon + attributes) for it
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'MIN_SIZE', float)
            logger.debug('MIN_SIZE = ' + str(self.cfg.get('CONFIG_PARAMS', 'MIN_SIZE')))
            # Maximal standard deviation of height inside a lake
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'STD_HEIGHT_MAX', float)
            logger.debug('STD_HEIGHT_MAX = ' + str(self.cfg.get('CONFIG_PARAMS', 'STD_HEIGHT_MAX')))
            
            # To improve PixC golocation (=True) or not (=False)
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'IMP_GEOLOC', bool)
            logger.debug('IMP_GEOLOC = ' + str(self.cfg.get('CONFIG_PARAMS', 'IMP_GEOLOC')))
            # Method to compute lake boundary or polygon hull
            # 0=convex hull 1=concav hull (1.0=with alpha param (default) 1.1=without) 2=concav hull radar vectorisation
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'HULL_METHOD', float, valeurs=[1, 1.1, 2])
            logger.debug('HULL_METHOD = ' + str(self.cfg.get('CONFIG_PARAMS', 'HULL_METHOD')))
            
            # Big lakes parameters for improved geoloc
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'BIGLAKE_MODEL', str, valeurs=["polynomial"])
            logger.debug('BIGLAKE_MODEL = ' + str(self.cfg.get('CONFIG_PARAMS', 'BIGLAKE_MODEL')))
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'BIGLAKE_MIN_SIZE', float)
            logger.debug('BIGLAKE_MIN_SIZE = ' + str(self.cfg.get('CONFIG_PARAMS', 'BIGLAKE_MIN_SIZE')))
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'BIGLAKE_GRID_SPACING', float)
            logger.debug('BIGLAKE_GRID_SPACING = ' + str(self.cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_SPACING')))
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'BIGLAKE_GRID_RES', float)
            logger.debug('BIGLAKE_GRID_RES = ' + str(self.cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_RES')))
            
            # 2.2 - ID section
            # Nb digits for counter of lakes in a tile or pass
            self.cfg.test_var_config_file('ID', 'NB_DIGITS', str)
            logger.debug('NB_DIGITS = ' + str(self.cfg.get('ID', 'NB_DIGITS')))
            
            # 2.3 - FILENAMES_PATTERN section
            
            # Product generator
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'PRODUCER', str)
            logger.debug('PRODUCER = ' + str(self.cfg.get('FILENAMES_PATTERN', 'PRODUCER')))
            # Composite Release IDentifier for LakeTile processing
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'LAKE_TILE_CRID', str)
            logger.debug('LAKE_TILE_CRID = ' + str(self.cfg.get('FILENAMES_PATTERN', 'LAKE_TILE_CRID')))
            # Composite Release IDentifier for LakeSP processing
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'LAKE_SP_CRID', str)
            logger.debug('LAKE_SP_CRID = ' + str(self.cfg.get('FILENAMES_PATTERN', 'LAKE_SP_CRID')))
            
            # PIXC product
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'PIXC_PREFIX', str)
            logger.debug('PIXC_PREFIX = ' + str(self.cfg.get('FILENAMES_PATTERN', 'PIXC_PREFIX')))
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'PIXC_PATTERN_PRINT', str)
            logger.debug('PIXC_PATTERN_PRINT = ' + str(self.cfg.get('FILENAMES_PATTERN', 'PIXC_PATTERN_PRINT')))
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'PIXC_PATTERN_IND', str)
            logger.debug('PIXC_PATTERN_IND = ' + str(self.cfg.get('FILENAMES_PATTERN', 'PIXC_PATTERN_IND')))
            
            # PIXCVecRiver product
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'PIXCVEC_RIVER_PREFIX', str)
            logger.debug('PIXCVEC_RIVER_PREFIX = ' + str(self.cfg.get('FILENAMES_PATTERN', 'PIXCVEC_RIVER_PREFIX')))
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'PIXCVEC_RIVER_PATTERN_PRINT', str)
            logger.debug('PIXCVEC_RIVER_PATTERN_PRINT = ' + str(self.cfg.get('FILENAMES_PATTERN', 'PIXCVEC_RIVER_PATTERN_PRINT')))
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'PIXCVEC_RIVER_PATTERN_IND', str)
            logger.debug('PIXCVEC_RIVER_PATTERN_IND = ' + str(self.cfg.get('FILENAMES_PATTERN', 'PIXCVEC_RIVER_PATTERN_IND')))
            
            # LakeTile product
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'LAKE_TILE_PREFIX', str)
            logger.debug('LAKE_TILE_PREFIX = ' + str(self.cfg.get('FILENAMES_PATTERN', 'LAKE_TILE_PREFIX')))
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'LAKE_TILE_PATTERN', str)
            logger.debug('LAKE_TILE_PATTERN = ' + str(self.cfg.get('FILENAMES_PATTERN', 'LAKE_TILE_PATTERN')))
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'LAKE_TILE_PATTERN_PRINT', str)
            logger.debug('LAKE_TILE_PATTERN_PRINT = ' + str(self.cfg.get('FILENAMES_PATTERN', 'LAKE_TILE_PATTERN_PRINT')))
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'LAKE_TILE_PATTERN_IND', str)
            logger.debug('LAKE_TILE_PATTERN_IND = ' + str(self.cfg.get('FILENAMES_PATTERN', 'LAKE_TILE_PATTERN_IND')))
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'LAKE_TILE_SHP_SUFFIX', str)
            logger.debug('LAKE_TILE_SHP_SUFFIX = ' + str(self.cfg.get('FILENAMES_PATTERN', 'LAKE_TILE_SHP_SUFFIX')))
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'LAKE_TILE_SHP_META_SUFFIX', str)
            logger.debug('LAKE_TILE_SHP_META_SUFFIX = ' + str(self.cfg.get('FILENAMES_PATTERN', 'LAKE_TILE_SHP_META_SUFFIX')))
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'LAKE_TILE_EDGE_SUFFIX', str)
            logger.debug('LAKE_TILE_EDGE_SUFFIX = ' + str(self.cfg.get('FILENAMES_PATTERN', 'LAKE_TILE_EDGE_SUFFIX')))
            self.cfg.test_var_config_file('FILENAMES_PATTERN', 'LAKE_TILE_PIXCVEC_SUFFIX', str)
            logger.debug('LAKE_TILE_PIXCVEC_SUFFIX = ' + str(self.cfg.get('FILENAMES_PATTERN', 'LAKE_TILE_PIXCVEC_SUFFIX')))
            
        # Error managed
        except service_error.ConfigFileError:
            logger.error("Error in the configuration file ", self.cfg.path_conf)
            raise
        # Warning error not managed !
        except Exception:
            logger.error("Something wrong happened during configuration file check!")
            raise

        return True
    

#######################################


if __name__ == '__main__':

    # 0 - Parse inline parameters
    PARSER = argparse.ArgumentParser(description="Compute SWOT LakeTile product\
            from one tile of PIXC product and its associated PIXCVecRiver product. ")
    PARSER.add_argument("command_file", help="command file (*.cfg)")
    ARGS = PARSER.parse_args()

    # 1 - Instantiate PGE
    PGE = PGELakeTile(ARGS.command_file)

    # 2 - Start PGE Lake Tile
    PGE.start()

    # 3 - Stop PGE Lake Tile
    PGE.stop()
