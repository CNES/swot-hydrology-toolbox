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
# VERSION:3.0.0:DM:#91:2021/03/12:Poursuite industrialisation
# VERSION:3.1.0:DM:#91:2021/05/21:Poursuite industrialisation
# VERSION:3.2.0:DM:#91:2021/10/27:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: pge_lake_sp.py
    :synopsis: Process PGE_L2_HR_LakeSP, i.e. generate L2_HR_LakeSP
    product from all tiles of L2_HR_LakeTile products

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR
                  Cécile Cazals - CS

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
from lxml import etree as ET
import multiprocessing
import numpy as np
import os

import cnes.common.service_config_file as service_config_file
import cnes.common.service_error as service_error
import cnes.common.service_logger as service_logger

import cnes.common.lib.my_timer as my_timer
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_variables as my_var
import cnes.common.lib_lake.lake_db as lake_db
import cnes.common.lib_lake.locnes_filenames as locnes_filenames
import cnes.common.lib_lake.locnes_products_shapefile as locnes_products_shapefile
import cnes.common.lib_lake.proc_lake as proc_lake

import cnes.sas.lake_sp.proc_pixc_sp as proc_pixc_sp
import cnes.sas.lake_sp.proc_pixc_vec_sp as proc_pixc_vec_sp


#######################################


class PGELakeSP():
    """
    Class PGELakeSP
    Pseudo PGE class to launch LakeSP computation
    """

    def __init__(self, cmd_file):
        """
        Constructor of PGELakeSP

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

        self.laketile_dir = my_params["laketile_dir"]
        self.output_dir = my_params["output_dir"]

        self.cycle_num = my_params["cycle_num"]
        self.pass_num = my_params["pass_num"]
        self.continent_id = my_params["continent_id"]

        # 2.1 - Init parameter file
        self.cfg = service_config_file.ServiceConfigFile(None)
        # 2.2 - Put command parameter inside cfg
        self._put_cmd_value(my_params)
        # 2.3 - Retrieve parameters which are optionnal
        self.flag_prod_shp = self.cfg.getboolean("OPTIONS", "Produce shp")
        self.flag_inc_file_counter = self.cfg.getboolean("OPTIONS", "Increment file counter")
        self.flag_write_full_path = self.cfg.getboolean("OPTIONS", "Write full path")
        self.flag_del_tmp_shp = self.cfg.getboolean("OPTIONS", "Delete temporary shp")

        # 3 - Set log filenames
        self.logfile = self.cfg.get('LOGGING', 'logFile')
        self.logfile_begin = self.logfile.replace(".log", "_BEGIN.log")
        self.logfile_end = self.logfile.replace(".log", "_END.log")
        self.logfile_swath = dict()
        self.cfg.set("LOGGING", "logFile", self.logfile_begin)

        # 4 - Init other class properties
        self.list_swath_sides = ["R", "L"]  # Swath sides to process (R=right and L=left)
        self.list_laketile_obs = list()  # List of LakeTile_Obs files
        self.list_laketile_prior = list()  # List of LakeTile_Prior files
        self.list_laketile_unknown = list()  # List of LakeTile_Unassigned files
        self.list_laketile_edge = dict()  # List of LakeTile_Edge files, organized per swath (R=right or L=left)
        self.list_laketile_pixcvec = dict()  # List of LakeTile_PIXCVec files, organized per swath (R=right or L=left)
        for cur_swath_side in self.list_swath_sides:
            self.list_laketile_edge[cur_swath_side] = list()
            self.list_laketile_pixcvec[cur_swath_side] = list()
            self.logfile_swath[cur_swath_side] = self.logfile.replace(".log", "_%s.log" % cur_swath_side)
        self.sum_task = 2  # Tasks may be separated per swath side
        self.nb_proc = int(my_params["nb_proc"])  # Number of cpu used 
        self.list_tmp_lakesp_shp = dict()
        self.list_tmp_lakesp_shp["obs"] = list()  # List of temporary LakeSP_Obs shapefiles (to be later combined)
        self.list_tmp_lakesp_shp["prior"] = list()  # List of temporary LakeSP_Prior shapefiles (to be later combined)
        self.list_tmp_lakesp_shp["unknown"] = list()  # List of temporary LakeSP_Unassigned shapefiles (to be later combined)

        # 5 - Initiate logging service
        service_logger.ServiceLogger()
        logger = logging.getLogger(self.__class__.__name__)

        # 6 - Print info
        logger.sigmsg("====================================")
        logger.sigmsg("===== LakeSPProcessing = BEGIN =====")
        logger.sigmsg("====================================")
        logger.info("> INPUT LakeTile directory = " + str(self.cfg.get("PATHS", "LakeTile directory")))
        logger.info("> OUTPUT directory = " + str(self.cfg.get("PATHS", "Output directory")))
        logger.info("======================================")
        message = "> Command file: " + str(self.cmd_file)
        logger.info(message)
        logger.info("=============service_logger=========================")
        
        # 7 - Test input parameters
        logger.info(">> Check config parameters type and value")
        self._check_config_parameters()

        # 8 - Form processing metadata dictionary
        self.proc_metadata = {}
        self.proc_metadata["history"] = "%sZ: Creation" % my_tools.swot_timeformat(datetime.datetime.utcnow(), in_format=3)
        self.proc_metadata["cycle_number"] = "%03d" % self.cycle_num
        self.proc_metadata["pass_number"] = "%03d" % self.pass_num
        self.proc_metadata["continent_id"] = self.continent_id
        self.proc_metadata["continent_code"] = lake_db.compute_continent_code(self.continent_id)
        if self.cfg.get("DATABASES", "LAKE_DB") is not None:
            self.proc_metadata["xref_prior_lake_db_file"] = self.cfg.get("DATABASES", "LAKE_DB")
        else:
            self.proc_metadata["xref_prior_lake_db_file"] = ""

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
        out_params["flag_write_full_path"] = "False"
        out_params["nb_proc"] = 1
        out_params["flag_del_tmp_shp"] = "True"

        # 1 - Read parameter file
        config = configparser.ConfigParser()
        config.read(self.cmd_file)

        # 2 - Retrieve PATHS
        out_params["laketile_dir"] = config.get("PATHS", "LakeTile directory")
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

        # 4 - Retrieve cycle and pass number
        out_params["cycle_num"] = config.getint("TILES_INFOS", "Cycle number")
        out_params["pass_num"] = config.getint("TILES_INFOS", "Pass number")
        list_tile_infos = config.options("TILES_INFOS")
        if "continent identifier" in list_tile_infos:
            out_params["continent_id"] = config.get("TILES_INFOS", "Continent identifier")
        else:
            out_params["continent_id"] = ""

        # 5 - Retrieve OPTIONS
        if "OPTIONS" in config.sections():
            list_options = config.options("OPTIONS")
            # Flag to also produce LakeTile_Edge and LakeTile_PIXCVec as shapefiles (=True); else=False (default)
            if "produce shp" in list_options:
                out_params["flag_prod_shp"] = config.get("OPTIONS", "Produce shp")
            # Flag to increment the file counter in the output filenames (=True, default); else=False
            if "increment file counter" in list_options:
                out_params["flag_inc_file_counter"] = config.get("OPTIONS", "Increment file counter")
            # To write full path in global attributes (=True); to write only basename=False
            if "write full path" in list_options:
                out_params["flag_write_full_path"] = config.get("OPTIONS", "Write full path")
            # Number of processors to use (default=1)
            if "nb_proc" in list_options:
                out_params["nb_proc"] = config.getint("OPTIONS", "Nb_proc")
                if out_params["nb_proc"] > 2:
                    out_params["nb_proc"] = 2
            # Flag to delete temporary swath LakeSP shapefiles (=True, default) or not (=False)
            if "delete temporary shp" in list_options:
                out_params["flag_del_tmp_shp"] = config.get("OPTIONS", "Delete temporary shp")

        # 6 - Retrieve LOGGING
        error_file = config.get("LOGGING", "errorFile")
        out_params["errorFile"] = os.path.splitext(error_file)[0] + "_" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + os.path.splitext(error_file)[1]
        log_file = config.get("LOGGING", "logFile")
        out_params["logFile"] = os.path.splitext(log_file)[0] + "_" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + os.path.splitext(log_file)[1]
        out_params["logfilelevel"] = config.get("LOGGING", "logfilelevel")
        out_params["logConsole"] = config.get("LOGGING", "logConsole")
        out_params["logconsolelevel"] = config.get("LOGGING", "logconsolelevel")
        
        # 7 - Retrieve FILE_INFORMATION
        # Name of producing agency
        out_params["INSTITUTION"] = config.get("FILE_INFORMATION", "INSTITUTION")
        # Product version
        out_params["PRODUCT_VERSION"] = config.get("FILE_INFORMATION", "PRODUCT_VERSION")
        # Composite Release IDentifier for LakeTile processing
        out_params["CRID_LAKETILE"] = config.get("FILE_INFORMATION", "CRID_LAKETILE")
        # Composite Release IDentifier for LakeSP processing
        out_params["CRID_LAKESP"] = config.get("FILE_INFORMATION", "CRID_LAKESP")
        # Version identifier of the product generation executable (PGE)
        out_params["PGE_VERSION"] = config.get("FILE_INFORMATION", "PGE_VERSION")
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
            self.cfg.set(section, "LakeTile directory", param_list["laketile_dir"])
            self.cfg.set(section, "Output directory", param_list["output_dir"])

            # Add DATABASES section and parameters
            section = "DATABASES"
            self.cfg.add_section(section)
            self.cfg.set(section, "LAKE_DB", param_list["LAKE_DB"])
            self.cfg.set(section, "LAKE_DB_ID", param_list["LAKE_DB_ID"])

            # Add TILES INFOS section and parameters
            section = "TILES_INFOS"
            self.cfg.add_section(section)
            self.cfg.set(section, "Cycle number", param_list["cycle_num"])
            self.cfg.set(section, "Pass number", param_list["pass_num"])
            self.cfg.set(section, "Continent identifier", param_list["continent_id"])

            # Add OPTIONS section and parameters
            section = "OPTIONS"
            self.cfg.add_section(section)
            self.cfg.set(section, "Produce shp", param_list["flag_prod_shp"])
            self.cfg.set(section, "Increment file counter", param_list["flag_inc_file_counter"])
            self.cfg.set(section, "Write full path", param_list["flag_write_full_path"])
            self.cfg.set(section, "Nb_proc", param_list["nb_proc"])
            self.cfg.set(section, "Delete temporary shp", param_list["flag_del_tmp_shp"])
            
            # Add section FILE_INFORMATION
            section = "FILE_INFORMATION"
            self.cfg.add_section(section)
            self.cfg.set(section, "INSTITUTION", param_list["INSTITUTION"])
            self.cfg.set(section, "PRODUCT_VERSION", param_list["PRODUCT_VERSION"])
            self.cfg.set(section, "CRID_LAKETILE", param_list["CRID_LAKETILE"])
            self.cfg.set(section, "CRID_LAKESP", param_list["CRID_LAKESP"])
            self.cfg.set(section, "PGE_VERSION", param_list["PGE_VERSION"])
            self.cfg.set(section, "CONTACT", param_list["CONTACT"])

        except Exception:
            print("Something wrong happened in _put_cmd_value !")
            raise

    def _check_config_parameters(self):
        """
        Check parameters coherence for LakeSP parameter file
        
        :return: True if OK
        :rtype: boolean
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        try:

            # Config parameters from command file
            
            # 1 - PATH section
            # LakeTile files
            self.cfg.test_var_config_file('PATHS', 'LakeTile directory', str)
            my_tools.test_dir(self.cfg.get('PATHS', 'LakeTile directory'))
            logger.debug('LakeTile directory = ' + str(self.cfg.get('PATHS', 'LakeTile directory')))
            # Output directory
            self.cfg.test_var_config_file('PATHS', 'Output directory', str)
            my_tools.test_dir(self.cfg.get('PATHS', 'Output directory'))
            logger.debug('Output directory = ' + str(self.cfg.get('PATHS', 'Output directory')))

            # 2 - DATABASES section
            # Lake database full path
            if self.cfg.get('DATABASES', 'LAKE_DB') is None:
                logger.debug('WARNING - LAKE_DB not filled => LakeTile product not linked to a lake database and continent')
            elif self.cfg.get('DATABASES', 'LAKE_DB').endswith(".shp"):
                self.cfg.test_var_config_file('DATABASES', 'LAKE_DB', str, logger=logger)
                my_tools.test_file(self.cfg.get('DATABASES', 'LAKE_DB'))
                logger.debug('OK - LAKE_DB = ' + str(self.cfg.get('DATABASES', 'LAKE_DB')))
                # Lake identifier attribute name in the database
                if self.cfg.get('DATABASES', 'LAKE_DB_ID'):
                    self.cfg.test_var_config_file('DATABASES', 'LAKE_DB_ID', str, logger=logger)
                    logger.debug('OK - LAKE_DB_ID = ' + str(self.cfg.get('DATABASES', 'LAKE_DB_ID')))
                else:
                    logger.warning('WARNING - LAKE_DB file given but the lake_id fieldname is missing')
            elif self.cfg.get('DATABASES', 'LAKE_DB').endswith(".sqlite"):
                self.cfg.test_var_config_file('DATABASES', 'LAKE_DB', str, logger=logger)
                my_tools.test_file(self.cfg.get('DATABASES', 'LAKE_DB'))
                logger.debug('OK - LAKE_DB = ' + str(self.cfg.get('DATABASES', 'LAKE_DB')))
            elif os.path.isdir(self.cfg.get('DATABASES', 'LAKE_DB')):
                self.cfg.test_var_config_file('DATABASES', 'LAKE_DB', str, logger=logger)
                my_tools.test_dir(self.cfg.get('DATABASES', 'LAKE_DB'))
                logger.debug('OK - LAKE_DB directory = ' + str(self.cfg.get('DATABASES', 'LAKE_DB')))
            else:
                logger.debug('WARNING - Unknown LAKE_DB file format for file: %s => LakeSP product not linked to a lake database' % self.cfg.get('DATABASES', 'LAKE_DB'))
              
            # 3 - TILES_INFO section
            # Cycle number
            self.cfg.test_var_config_file('TILES_INFOS', 'Cycle number', int)
            logger.debug('OK - Cycle number = ' + str(self.cfg.get('TILES_INFOS', 'Cycle number')))
            # Pass number
            self.cfg.test_var_config_file('TILES_INFOS', 'Pass number', int)
            logger.debug('OK - Pass number = ' + str(self.cfg.get('TILES_INFOS', 'Pass number')))
            # Continent identifier
            self.cfg.test_var_config_file('TILES_INFOS', 'Continent identifier', str)
            logger.debug('OK - Continent identifier = ' + str(self.cfg.get('TILES_INFOS', 'Continent identifier')))
            
            # 4 - OPTIONS section
            # Shapefile production
            self.cfg.test_var_config_file('OPTIONS', 'Produce shp', bool, logger=logger)
            logger.debug('OK - Produce shp = ' + str(self.cfg.get('OPTIONS', 'Produce shp')))
            # Increment output file counter
            self.cfg.test_var_config_file('OPTIONS', 'Increment file counter', bool, logger=logger)
            logger.debug('OK - Increment file counter = ' + str(self.cfg.get('OPTIONS', 'Increment file counter')))
            # Write full path in global attributes
            self.cfg.test_var_config_file('OPTIONS', 'Write full path', bool, logger=logger)
            logger.debug('OK - Write full path = ' + str(self.cfg.get('OPTIONS', 'Write full path')))
            # Number of processors to use 
            self.cfg.test_var_config_file('OPTIONS', 'Nb_proc', int, logger=logger)
            logger.debug('OK - Nb_proc = ' + str(self.cfg.get('OPTIONS', 'Nb_proc')))
            # To delete temporary swath LakeSP shapefiles (=True, default) or not (=False)
            self.cfg.test_var_config_file('OPTIONS', 'Delete temporary shp', bool, val_default=True, logger=logger)
            logger.debug('OK - Delete temporary shp = ' + str(self.cfg.get('OPTIONS', 'Delete temporary shp')))
            
            # 5 - FILE_INFORMATION section
            # Name of producing agency
            self.cfg.test_var_config_file('FILE_INFORMATION', 'INSTITUTION', str, logger=logger)
            logger.debug('OK - INSTITUTION = ' + str(self.cfg.get('FILE_INFORMATION', 'INSTITUTION')))
            # Product version
            self.cfg.test_var_config_file('FILE_INFORMATION', 'PRODUCT_VERSION', str, logger=logger)
            logger.debug('OK - PRODUCT_VERSION = ' + str(self.cfg.get('FILE_INFORMATION', 'PRODUCT_VERSION')))
            # Composite Release IDentifier for LakeTile processing
            self.cfg.test_var_config_file('FILE_INFORMATION', 'CRID_LAKETILE', str, logger=logger)
            logger.debug('OK - CRID_LAKETILE = ' + str(self.cfg.get('FILE_INFORMATION', 'CRID_LAKETILE')))
            # Composite Release IDentifier for LakeSP processing
            self.cfg.test_var_config_file('FILE_INFORMATION', 'CRID_LAKESP', str, logger=logger)
            logger.debug('OK - CRID_LAKESP = ' + str(self.cfg.get('FILE_INFORMATION', 'CRID_LAKESP')))
            # Version identifier of the product generation executable (PGE)
            self.cfg.test_var_config_file('FILE_INFORMATION', 'PGE_VERSION', str, logger=logger)
            logger.debug('OK - PGE_VERSION = ' + str(self.cfg.get('FILE_INFORMATION', 'PGE_VERSION')))
            # Contact
            self.cfg.test_var_config_file('FILE_INFORMATION', 'CONTACT', str, logger=logger)
            logger.debug('OK - CONTACT = ' + str(self.cfg.get('FILE_INFORMATION', 'CONTACT')))

        # Error managed
        except service_error.ConfigFileError:
            message = "Error in the command file " + self.cfg.path_conf
            logger.error(message, exc_info=True)
            raise
            
        # Warning error not managed !
        except Exception:
            logger.error("Something wrong happened during command file check!", exc_info=True)
            raise

        return True
        
    # -------------------------------------------

    def _select_input_data(self):
        """
        This method selects appropriate LakeTile shapefiles, and classify them depending on the swath related to them
        """
        logger = logging.getLogger(self.__class__.__name__)

        # 0 - Compute lake DB filename if directory is given
        lake_db_file = self.cfg.get("DATABASES", "LAKE_DB")
        if (lake_db_file is not None) and os.path.isdir(lake_db_file):
            obj_lakedb_filenames = locnes_filenames.PldFilenames("SP", 
                                                                 lake_db_file, 
                                                                 cycle_num=self.cycle_num, 
                                                                 pass_num=self.pass_num)
            lake_db_file = obj_lakedb_filenames.filename

        # 1 - List all files in LakeTile_shp directory
        tmp_list_laketile_files = os.listdir(self.laketile_dir)

        # 2 - Compute LakeTile_Obs prefix regarding cycle and pass numbers
        cond_prefix = locnes_filenames.LAKE_TILE_PREFIX["obs"]  # Generic LakeTile_Obs prefix 
        cond_prefix += "%03d" % self.cycle_num  # Add cycle number to prefix
        cond_prefix += "_%03d" % self.pass_num  # Add pass number to prefix
        logger.debug("Prefix to select LakeTile_Obs shapefiles = %s" % cond_prefix)

        # 3 - For each selected LakeTile_Obs shapefile
        # -> verify if they correspond to specified cycle and pass numbers
        # -> get related LakeTile_Prior, LakeTile_Unassigned, LakeTile_Edge and LakeTile_PIXCVec files
        #    (NB: if at least one of these files is missing, then the corresponding tile is not processed)
        
        # Init variables
        proc_metadata_keys = self.list_swath_sides + ["all"]
        cpt_input_products = dict()  # Input files counter
        list_laketile_root = dict()  # List of LakeTile root name
        list_pixc_files = dict()  # List of PIXC files
        list_basin_code = dict()  # List of basin_code
        for key in proc_metadata_keys:
            cpt_input_products[key] = 0
            list_laketile_root[key] = list()
            list_pixc_files[key] = list()
            list_basin_code[key] = set()
        continent_code = str(lake_db.compute_continent_code(self.continent_id))
        
        for cur_file in tmp_list_laketile_files:  # Foreach file in LakeTile directory

            # 3.1 - Test if file meets the condition of LakeTile_Obs file
            if cur_file.startswith(cond_prefix) and cur_file.endswith(locnes_filenames.LAKE_TILE_SHP_SUFFIX):
                logger.debug("Working with current LakeTile_Obs file: %s" % cur_file)

                # Construct LakeTile files filenames and verify existance
                cur_laketile_obs_path = os.path.join(self.laketile_dir, cur_file)  # LakeTile_Obs filename
                flag_laketile_product_complete = True  # Init flag LakeTile complete, ie (_Obs, _Unassigned, _Edge, _PIXCVec) files exist
                # LakeTile_Prior
                cur_laketile_prior_path = os.path.join(self.laketile_dir, 
                                                       cur_file.replace(locnes_filenames.LAKE_TILE_PREFIX["obs"],
                                                                        locnes_filenames.LAKE_TILE_PREFIX["prior"]))
                if not os.path.exists(cur_laketile_prior_path):
                    logger.debug("-> Associated LakeTile_Prior file is missing: %s" % cur_laketile_prior_path)
                    flag_laketile_product_complete = False
                # LakeTile_Unassigned
                cur_laketile_unknown_path = os.path.join(self.laketile_dir, 
                                                         cur_file.replace(locnes_filenames.LAKE_TILE_PREFIX["obs"],
                                                                          locnes_filenames.LAKE_TILE_PREFIX["unknown"]))
                if not os.path.exists(cur_laketile_unknown_path):
                    logger.debug("-> Associated LakeTile_Unassigned file is missing: %s" % cur_laketile_unknown_path)
                    flag_laketile_product_complete = False
                # LakeTile_Edge
                cur_laketile_root = cur_file.replace(locnes_filenames.LAKE_TILE_PREFIX["obs"],
                                                     locnes_filenames.LAKE_TILE_PREFIX_BASE).replace(locnes_filenames.LAKE_TILE_SHP_SUFFIX, "")
                cur_laketile_edge_path = os.path.join(self.laketile_dir,
                                                      cur_file.replace(locnes_filenames.LAKE_TILE_PREFIX["obs"],
                                                                       locnes_filenames.LAKE_TILE_PREFIX["edge"]).replace(locnes_filenames.LAKE_TILE_SHP_SUFFIX, locnes_filenames.LAKE_TILE_EDGE_SUFFIX))  # LakeTile_Edge filename
                if not os.path.exists(cur_laketile_edge_path):
                    logger.debug("-> Associated LakeTile_Edge file is missing: %s" % cur_laketile_edge_path)
                    flag_laketile_product_complete = False
                # LakeTile_PIXCVec
                cur_laketile_pixcvec_path = os.path.join(self.laketile_dir,
                                                         cur_file.replace(locnes_filenames.LAKE_TILE_PREFIX["obs"],
                                                                          locnes_filenames.LAKE_TILE_PREFIX["pixcvec"]).replace(locnes_filenames.LAKE_TILE_SHP_SUFFIX, locnes_filenames.LAKE_TILE_PIXCVEC_SUFFIX))  # LakeTile_PIXCVec filename
                if not os.path.exists(cur_laketile_pixcvec_path):
                    logger.debug("-> Associated LakeTile_PIXCVec file is missing: %s" % cur_laketile_pixcvec_path)
                    flag_laketile_product_complete = False
                
                # 3.2 - If the LakeTile product is complete, ie (_Obs, _Prior, _Unassigned, _Edge, _PIXCVec) files exist
                if flag_laketile_product_complete:
                    logger.debug("-> LakeTile product is COMPLETE")
                    
                    # Get continent, PIXC and PIXCVecRiver information from LakeTile_shp metadata file (ie .shp.xml)
                    metadata = ET.parse(cur_laketile_obs_path + ".xml")
                    
                    # Test and overwrite configuration parameters with XML file content
                    try:
                        locnes_products_shapefile.set_config_from_xml(metadata)
                        locnes_products_shapefile.check_lake_db(metadata, lake_db_file)
                        locnes_products_shapefile.check_lake_db_id(metadata)
                    except:
                        message = "ERROR: problem with LakeTile_Obs XML file: " + cur_laketile_obs_path + ".xml"
                        raise service_error.ProcessingError(message, logger)

                    # Retrieve continent_id
                    cur_continent_id_txt = metadata.xpath("//swot_product/global_metadata/continent_id")[0].text
                    if not cur_continent_id_txt:
                        cur_continent_id_list = [""]
                    else:
                        cur_continent_id_list = cur_continent_id_txt.split(my_var.SEP_ATT)
                        
                    # Retrieve basin_code 
                    cur_basin_code_txt = metadata.xpath("//swot_product/global_metadata/basin_code")[0].text
                    if not cur_basin_code_txt:
                        cur_basin_code_list = [""]
                    else:
                        cur_basin_code_list = cur_basin_code_txt.split(my_var.SEP_ATT)
                        
                    # Retrieve LakeTile processing input files (only PIXC)
                    cur_pixc_path = metadata.xpath("//swot_product/global_metadata/xref_l2_hr_pixc_file")[0].text
                                        
                    # 3.3 - Add current information to lists related to files to process
                    cur_swath_side = locnes_filenames.get_info_from_filename(cur_file, "LakeTile")["tile_ref"][-1]
                    if (self.continent_id == "") or (self.continent_id in cur_continent_id_list):

                        # LakeTile_Obs full path
                        self.list_laketile_obs.append(cur_laketile_obs_path)
                        # LakeTile_Prior full path
                        self.list_laketile_prior.append(cur_laketile_prior_path)
                        # LakeTile_Unassigned full path
                        self.list_laketile_unknown.append(cur_laketile_unknown_path)
                        # LakeTile_Edge full path
                        self.list_laketile_edge[cur_swath_side].append(cur_laketile_edge_path)
                        # LakeTile_PIXCVec full path
                        self.list_laketile_pixcvec[cur_swath_side].append(cur_laketile_pixcvec_path)
                        # LakeTile root names
                        if self.flag_write_full_path:
                            cur_laketile_root = os.path.join(self.laketile_dir, cur_laketile_root)
                        list_laketile_root[cur_swath_side].append(cur_laketile_root)
                        list_laketile_root["all"].append(cur_laketile_root)

                        # PIXC full path
                        list_pixc_files[cur_swath_side].append(cur_pixc_path)
                        list_pixc_files["all"].append(cur_pixc_path)
                        
                        # Update basin_code set 
                        for cur_basin_code in cur_basin_code_list:
                            if str(cur_basin_code).startswith(continent_code):
                                list_basin_code[cur_swath_side].add(cur_basin_code)
                                list_basin_code["all"].add(cur_basin_code)
                                
                        # Increment counters
                        cpt_input_products[cur_swath_side] += 1
                        cpt_input_products["all"] += 1

                    # 4 - Populate processing metadata

                    # Parameter file from LakeTile processing
                    self.proc_metadata["source"] = metadata.xpath("//swot_product/global_metadata/source")[0].text
                    self.proc_metadata["xref_param_l2_hr_laketile_file"] = metadata.xpath("//swot_product/global_metadata/xref_param_l2_hr_laketile_file")[0].text

        # Parameters specific to swath and global
        self.proc_metadata_specific = dict()
        for key in proc_metadata_keys:
            self.proc_metadata_specific[key] = dict()  # Init
            
            # basin_code list
            self.proc_metadata_specific[key]["basin_code"] = my_var.SEP_ATT.join(sorted(list_basin_code[key]))

            # PIXC filenames
            self.proc_metadata_specific[key]["xref_l2_hr_pixc_files"] = my_var.SEP_FILES.join(sorted(list_pixc_files[key]))
            # LakeTile filenames
            self.proc_metadata_specific[key]["xref_l2_hr_laketile_files"] = my_var.SEP_FILES.join(sorted(list_laketile_root[key]))
        
        # Print some numbers
        for cur_swath_side in self.list_swath_sides:
            logger.debug("Number of input LakeTile products for swath %s = %d" % (cur_swath_side, cpt_input_products[cur_swath_side]))
            if cpt_input_products[cur_swath_side] > 0:
                # Get all tile ref and sort them
                list_tiles = list()
                for cur_file in self.list_laketile_edge[cur_swath_side]:
                    list_tiles.append(locnes_filenames.get_info_from_filename(cur_file, "LakeTile")["tile_ref"])
                list_tiles_sorted = sorted(list_tiles)
                # Prints
                logger.debug("  First tile = %s" % list_tiles_sorted[0])
                logger.debug("  Last tile = %s" % list_tiles_sorted[-1])
                # Compute list of all tiles
                for cur_num in list(range(int(list_tiles_sorted[0][:-1]), int(list_tiles_sorted[-1][:-1])+1)):
                    cur_tile = "%03d%s" % (cur_num, cur_swath_side)
                    if cur_tile in list_tiles_sorted:
                        list_tiles_sorted.remove(cur_tile)
                # Print
                if len(list_tiles_sorted) == 0:
                    logger.debug("  NO missing tiles")
                else:
                    logger.debug("  Missing tiles = %s" % ";".join(list_tiles_sorted))
                    
        if len(self.list_laketile_obs) > 0:

            # 5 - Retrieve tile info from LakeTile filename and compute output filenames
            logger.debug("Retrieving tile infos from LakeTile filename...")
            self.lake_sp_filenames = locnes_filenames.LakeSPFilenames(self.list_laketile_obs,
                                                                      self.continent_id,
                                                                      self.output_dir,
                                                                      flag_inc=self.flag_inc_file_counter)

    def _read_lake_db(self, in_continent_id, in_poly, in_swath):
        """
        This method prepares Lake DB class, and init with Lake DB data
        
        :param in_continent_id: 1-digit continent identifier on which to focus
        :type in_continent_id: int
        :param in_poly: polygon area for selection of PLD lakes
        :type in_poly: OGRPolygon
        :param in_swath: "R" or "L" for Right ou Left swath
        :type in_swath: String

        :return: out_obj_lake_db = subset of the PLD 
        :rtype: cnes.common.lib_lake.lake_db.<LakeDb|LakeDbShp|LakeDbSqlite>
        
        Possible combinations:
            - in_continent_id=not None, others=None => LakeSP processing with operational PLD
            - in_tile_poly=not None, others=None => LakeSP processing with user lake DB
        """
        logger = logging.getLogger(self.__class__.__name__)
        timer = my_timer.Timer()
        timer.start()

        lake_db_file = self.cfg.get("DATABASES", "LAKE_DB")
        
        if (lake_db_file == "") or (lake_db_file is None):  # No PLD
            logger.warning("NO database specified -> NO link of SWOT obs with a priori lake")
            out_obj_lake_db = lake_db.LakeDb()
                
        else:  # Init PLD object wrt to file type
            
            type_db = lake_db_file.split('.')[-1]  # File type
            continent_code = str(lake_db.compute_continent_code(in_continent_id))
            
            if type_db == "shp":  # User lake DB file in shapefile format
                logger.debug("Use of personal user lake database in shapefile format")
                out_obj_lake_db = lake_db.LakeDbShp(lake_db_file, 
                                                    in_poly=in_poly)
                
            elif type_db == "sqlite":  # User or operational lake DB file in SQLite format
                if os.path.basename(lake_db_file).startswith(locnes_filenames.PLD_PREFIX_BASE):  # Operationnal PLD filename
                    logger.debug("Use of operational Prior Lake Database (specific filename given)")
                    out_obj_lake_db = lake_db.LakeDbSqlite(lake_db_file, 
                                                           in_continent_id=continent_code,
                                                           in_swath=in_swath)
                else:
                    logger.debug("Use of personal user lake database in SQLite format")
                    out_obj_lake_db = lake_db.LakeDbSqlite(lake_db_file, 
                                                           in_poly=in_poly)
                    
            elif os.path.isdir(lake_db_file):  # Directory of operational Prior Lake Database (PLD)
                logger.debug("Use of operational Prior Lake Database (upper directory given)")
                # Compute PLD filename
                lakedb_filenames = locnes_filenames.PldFilenames("SP", 
                                                                 lake_db_file, 
                                                                 cycle_num=self.cycle_num,
                                                                 pass_num=self.pass_num)
                # Init PLD object
                if os.path.exists(lakedb_filenames.filename):
                    out_obj_lake_db = lake_db.LakeDbSqlite(lakedb_filenames.filename, 
                                                           in_continent_id=continent_code,
                                                           in_swath=in_swath)
                    self.proc_metadata["xref_prior_lake_db_file"] = lakedb_filenames.filename
                else:
                    logger.warning("NO related operational Prior Lake Database file found -> NO database specified")
                    self.obj_lake_db = lake_db.LakeDb()
                    
            else:
                message = "Prior Lake Database format (%s) is unknown: must be .shp or .sqlite => set to None" % type_db
                logger.warning(message)
                out_obj_lake_db = lake_db.LakeDb()

        logger.debug(timer.stop("PLD opening time: "))

        return out_obj_lake_db
        
    # -------------------------------------------

    def _write_output_data(self):
        """
        This method writes output data
        """
        logger = logging.getLogger(self.__class__.__name__)
        timer = my_timer.Timer()
        timer.start()

        # 1 - Write LakeSP shapefiles
        logger.info("> Combining LakeSP temporary swath shapefiles to unique LakeSP shapefiles...")
        # 1.1 - Obs-oriented shapefile
        logger.info(". Obs-oriented shapefile = %s" % self.lake_sp_filenames.lake_sp_file_obs)
        proc_lake.write_lakesp_file(self.list_tmp_lakesp_shp["obs"], 
                                    self.list_laketile_obs,
                                    self.lake_sp_filenames.lake_sp_file_obs, 
                                    {**self.proc_metadata, **self.proc_metadata_specific["all"]},
                                    continent_id=self.continent_id,
                                    flag_del_tmp_shp=self.flag_del_tmp_shp)
        # 1.2 - Shapefile of prior water features
        logger.info(". PLD-oriented shapefile = %s" % self.lake_sp_filenames.lake_sp_file_prior)
        proc_lake.write_lakesp_file(self.list_tmp_lakesp_shp["prior"], 
                                    self.list_laketile_prior,
                                    self.lake_sp_filenames.lake_sp_file_prior,
                                    {**self.proc_metadata, **self.proc_metadata_specific["all"]},
                                    continent_id=self.continent_id,
                                    flag_del_tmp_shp=self.flag_del_tmp_shp,
                                    flag_add_nogeom_feat=True)
        # 1.3 - Shapefile of unassigned water features
        logger.info(". Shapefile of unassigned water features = %s" % self.lake_sp_filenames.lake_sp_file_unknown)
        proc_lake.write_lakesp_file(self.list_tmp_lakesp_shp["unknown"], 
                                    self.list_laketile_unknown,
                                    self.lake_sp_filenames.lake_sp_file_unknown,
                                    {**self.proc_metadata, **self.proc_metadata_specific["all"]},
                                    flag_del_tmp_shp=self.flag_del_tmp_shp,
                                    continent_id=self.continent_id)
        logger.debug(timer.stop("LakeSP file writing time: "))
        logger.info("")

    # -------------------------------------------
    
    def main_process(self, in_list_laketile_edge, in_list_laketile_pixcvec):
        """
        LakeSP main process
        
        :param in_list_laketile_edge: dictionnary of lists of LakeTile_Edge files to process; key = swath side (R or L)
        :type in_list_laketile_edge: dict
        :param in_list_laketile_pixcvec: dictionnary of lists of LakeTile_PIXCVec files to process; key = swath side (R or L)
        :type in_list_laketile_pixcvec: dict
        
        :return: out_list_tmp_lakesp_shp = list of temporary output LakeSP shapefiles (key = obs | prior | unknown); 
                    there is a one-to-one correspondance with the swath sides above
        :rtype: out_list_tmp_lakesp_shp = dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 0.1 - Init list of temporary LakeSP shapefiles
        out_list_tmp_lakesp_shp = dict()
        out_list_tmp_lakesp_shp["obs"] = list()
        out_list_tmp_lakesp_shp["prior"] = list()
        out_list_tmp_lakesp_shp["unknown"] = list()

        for cur_swath in in_list_laketile_edge.keys():
            timer = my_timer.Timer()
            timer.start()

            logger.info("=============================")
            logger.info("== Processing swath side %s ==" % cur_swath)
            logger.info("=============================")
            logger.info("")
            
            # 1 - Init and populate PixCEdgeSwath object
            
            # 1.1 - Init PixCEdgeSwath object
            logger.info("> 1.1 - Init PixCEdgeSwath object")
            obj_pixc_sp = proc_pixc_sp.PixCEdgeSwath(self.cycle_num,
                                                     self.pass_num,
                                                     self.continent_id,
                                                     cur_swath)
            logger.info("")
            
            # 1.2 - Populate PixCEdgeSwath object with LakeTile_Edge file information
            logger.info("> 1.2 - Fill PixCEdgeSwath object with LakeTile_Edge file information")
            obj_pixc_sp.set_from_laketile_edge_file(in_list_laketile_edge[cur_swath])
            logger.debug("=> %d edge pixels to manage" % obj_pixc_sp.nb_pixels)
            logger.info("")
            
            # 1.3 - Init PIXCVec_SP
            logger.info("> 1.3 - Init PIXCVec_SP object")
            obj_pixcvec_sp = proc_pixc_vec_sp.PixCVecSwath(in_list_laketile_pixcvec[cur_swath],
                                                           obj_pixc_sp,
                                                           self.continent_id)
            logger.info("")
            
            # Processing only if there are edge pixels to manage
            if obj_pixc_sp.nb_pixels > 0:
            
                # 2 - Read lake DB
                logger.info("> 2 - Init and read lake DB object...")
                obj_lakedb = self._read_lake_db(self.continent_id, obj_pixc_sp.tile_poly, cur_swath)
                logger.info("")
                
                # 3 - Initialize LakeSP object
                logger.info("> 3 - Initialize LakeSP object...")
                obj_lake_sp = proc_lake.LakeProduct("SP",
                                                    obj_pixc_sp,
                                                    obj_pixcvec_sp,
                                                    obj_lakedb,
                                                    os.path.basename(self.lake_sp_filenames.lake_sp_file_obs).split(".")[0])
                logger.info("")
                
                # 4 - Process swath
                
                # 4.1 - Gather edge pixels in separate entities for all the tiles of this swath
                logger.info("4.1 - Gathering edge pixels in separate entities for all the tiles...")
                obj_pixc_sp.swath_global_relabeling()
                logger.info("")

                # 4.2 - Compute LakeSP product for this swath
                logger.info("4.2 - Computing LakeSP features at along-track edges...")
                obj_lake_sp.compute_lake_features(np.unique(obj_pixc_sp.labels))
                logger.info("")
                
                # 5 - Close lake database
                logger.info("> 5 - Closing lake database...")
                obj_lakedb.close_db()
                
            else:
                
                # 3 - Initialize LakeSP object
                logger.info("> 3 - Initialize EMPTY LakeSP object...")
                obj_lake_sp = proc_lake.LakeProduct("SP",
                                                    obj_pixc_sp,
                                                    None,
                                                    None,
                                                    os.path.basename(self.lake_sp_filenames.lake_sp_file_obs).split(".")[0])
                logger.info("")
            
            # 6 - Write temporary output shapefiles
            logger.info("> 6 - Writing LakeSP temporary shapefiles...")
            
            # 6.1 - Obs-oriented shapefile
            # 6.1.1 - Compute filename
            tmp_filename_obs = os.path.join(self.output_dir, "%s_%s.shp" % (os.path.splitext(os.path.basename(self.lake_sp_filenames.lake_sp_file_obs))[0], cur_swath))
            out_list_tmp_lakesp_shp["obs"].append(tmp_filename_obs)
            # 6.1.2 - Write shapefile
            logger.info(". Obs-oriented shapefile = %s" % tmp_filename_obs)
            obj_lake_sp.write_obs_file(tmp_filename_obs, {**self.proc_metadata, **self.proc_metadata_specific[cur_swath]})
            
            # 6.2 - Shapefile of prior water features
            # 6.2.1 - Compute filename
            tmp_filename_prior = os.path.join(self.output_dir, "%s_%s.shp" % (os.path.splitext(os.path.basename(self.lake_sp_filenames.lake_sp_file_prior))[0], cur_swath))
            out_list_tmp_lakesp_shp["prior"].append(tmp_filename_prior)
            # 6.2.2 - Write shapefile
            logger.info(". PLD-oriented shapefile = %s" % tmp_filename_prior)
            obj_lake_sp.write_prior_file(tmp_filename_prior, {**self.proc_metadata, **self.proc_metadata_specific[cur_swath]})
            
            # 6.3 - Shapefile of unassigned water features
            # 6.3.1 - Compute filename
            tmp_filename_unknown = os.path.join(self.output_dir, "%s_%s.shp" % (os.path.splitext(os.path.basename(self.lake_sp_filenames.lake_sp_file_unknown))[0], cur_swath))
            out_list_tmp_lakesp_shp["unknown"].append(tmp_filename_unknown)
            # 6.3.2 - Write shapefile
            logger.info(". Shapefile of unassigned water features = %s" % tmp_filename_unknown)
            obj_lake_sp.write_unknown_file(tmp_filename_unknown, {**self.proc_metadata, **self.proc_metadata_specific[cur_swath]})
            logger.info("")
            
            # 7 - Update PIXCVec files
            logger.info("> 7 - Updating L2_HR_PIXCVec files...")
            obj_pixcvec_sp.update_pixcvec(self.output_dir,
                                          in_write_to_shp=self.flag_prod_shp,
                                          in_proc_metadata=self.proc_metadata)
            logger.info("")
                
            # 8 - Free memory
            del obj_pixc_sp
            del obj_pixcvec_sp
            obj_lake_sp.free_memory()
            del obj_lake_sp

            logger.debug(timer.stop("%s swath computation time: " %cur_swath))
            logger.info("")

        return out_list_tmp_lakesp_shp

    def start(self):
        """
        Main pseudoPGE method
        Start computation
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        ##################################
        # Definition of process function #
        ##################################
        def process_function(in_list_laketile_edge, 
                             in_list_laketile_pixcvec, 
                             queue):
            """
            Process function
        
            :param in_list_laketile_edge: dictionnary of lists of LakeTile_Edge files to process; key = swath side (R or L)
            :type in_list_laketile_edge: dict
            :param in_list_laketile_pixcvec: dictionnary of lists of LakeTile_PIXCVec files to process; key = swath side (R or L)
            :type in_list_laketile_pixcvec: dict
            :param queue: process queue
            :type queue: multiprocessing.Queue
            """
            logger = logging.getLogger(self.__class__.__name__)

            # 1 - Open new log file for each swath if a single swath is processed
            if len(in_list_laketile_edge) == 1:  # if a single swath is processed
                swath = next(iter(in_list_laketile_edge))
                self.cfg.set("LOGGING", "logFile", self.logfile_swath[swath])
                service_logger.ServiceLogger()
                logger = logging.getLogger(self.__class__.__name__)
                logger.info("=============================")
                logger.info("Principal log file in written to file %s" % self.logfile)
                logger.info("=============================")
                logger.info("Log file for swath %s will be written to log file %s " % (swath, self.logfile_swath[swath]))
                logger.info("=============================")

            # 2 - Run main process
            try:
                list_tmp_lakesp_shp = self.main_process(in_list_laketile_edge, in_list_laketile_pixcvec)
                queue.put(list_tmp_lakesp_shp)
                logger.info("STOP process")
                
            except:
                proc_pid = multiprocessing.current_process().pid
                proc_name = multiprocessing.current_process().name
                logger.error("Process " + proc_name + " with pid " + str(proc_pid) + " failed", exc_info=True)
                queue.put(None)
                raise

            # 3 - Close log file for each  swath if a single swath is processed
            if len(in_list_laketile_edge) == 1:
                swath = next(iter(in_list_laketile_edge))
                logger.debug("Closing temporary log file %s" % self.logfile_swath[swath])
                instance_logger = service_logger.get_instance()
                if instance_logger is not None:
                    instance_logger.flush_log()
                    instance_logger.close()

        #########################################
        # End of definition of process function #
        #########################################
        
        try:
            
            logger.info("*****************************")
            logger.info("***** Select input data *****")
            logger.info("*****************************")
            logger.info("")
            
            # 1 - Select input data
            logger.info("> Step 1 - Select and classify input data per swath side")
            self._select_input_data()
            logger.info("")
            
            if len(self.list_laketile_obs) == 0:
                logger.warning("!!! NO selected LakeTile product in input => NO LakeSP product in output !!!")
                logger.info("")
                
            else:

                if self.nb_proc == 1:
                    # 2 - Run main process directly
                    logger.info("> Step 2 - Main processing using ONLY 1 processor")
                    logger.info("")
                    self.list_tmp_lakesp_shp = self.main_process(self.list_laketile_edge, self.list_laketile_pixcvec)
                    
                else:
                    # 2 - Run main process after spreading over available processors
                    logger.info("> Step 2 - Main processing using %d processors" % self.nb_proc)
                    nb_proc = int(self.nb_proc)
                    nb_task_per_proc = int(self.sum_task/nb_proc)
                    message = "sum_nb_task: " + str(self.sum_task)
                    logger.debug(message)
                    message = "nb_task_per_proc: " + str(nb_task_per_proc)
                    logger.debug(message)
                    logger.info("")
    
                    logger.info("****************************************")
                    logger.info("***** Prepare and launch processes *****")
                    logger.info("****************************************")
                    logger.info("")
         
                    list_process = []
                    list_queue = []
        
                    # 2.1 - Prepare processes
                    logger.info("> 2.1 - Prepare processes")
                    for proc_i, cur_swath_side in zip(range(nb_proc), self.list_swath_sides):
                        # 2.1.1 - Select list of LakeTile files to process, defined by their related swath_side
                        list_laketile_edge_to_process = dict()
                        list_laketile_edge_to_process[cur_swath_side] = self.list_laketile_edge[cur_swath_side]
                        list_laketile_pixcvec_to_process = dict()
                        list_laketile_pixcvec_to_process[cur_swath_side] = self.list_laketile_pixcvec[cur_swath_side]
                        logger.debug("swath side for process nb " + str(proc_i) + ": " + str(cur_swath_side))
                        # 2.1.2 - Create Queue
                        queue = multiprocessing.Queue()
                        list_queue.append(queue)
                        # 2.1.3 - Create Process
                        nom = "sub_process_" + str(proc_i)
                        p = multiprocessing.Process(name=nom, target=process_function, args=(list_laketile_edge_to_process, list_laketile_pixcvec_to_process, queue))
                        list_process.append(p)
    
                        logger.debug("Log for swath %s will be written to log %s" % (cur_swath_side, self.logfile_swath[cur_swath_side]))
    
                    logger.info("")
        
                    # 2.2 - Launch processes
                    logger.info("> 2.2 - Launch processes")
                    proc_i = 0
    
                    # 2.3.0 - Write log to file before running multiprocessing
                    logger.debug("Closing temporary log %s to write a log for each swath" % self.logfile_begin)
                    instance_logger = service_logger.get_instance()
                    if instance_logger is not None:
                        instance_logger.flush_log()
                        instance_logger.close()
    
                    # 2.3.1 - Run subproc
                    for p in list_process:
                        logger.debug("Start process " + str(proc_i))
                        proc_i = proc_i + 1
                        p.start()
    
                    # 2.3.2 -  Reopen principal log file
                    self.cfg.set("LOGGING", "logFile", self.logfile_end)
                    service_logger.ServiceLogger()
                    logger = logging.getLogger(self.__class__.__name__)
                    logger.debug("Open log file %s" % self.logfile_end)
                    logger.info("")
        
                    # 2.4 - Get results from processes
                    logger.info("> 2.3 - Get results")
                    for q in list_queue:
                        list_tmp_lakesp_shp = q.get()
                        if list_tmp_lakesp_shp is not None:
                            for key in self.list_tmp_lakesp_shp.keys():
                                if key in list_tmp_lakesp_shp.keys():
                                    self.list_tmp_lakesp_shp[key] = [*self.list_tmp_lakesp_shp[key], *list_tmp_lakesp_shp[key]]
                        else:
                            message = "Fatal Error in one subProcess ! Stop PGE !"
                            #logger.error(message, exc_info=True)
                            raise Exception(message)
                    logger.info("")
        
                    logger.debug("Wait until process finish")
                    # Wait until process finish
                    for p in list_process:
                        p.join()
                    
                logger.info("")
                logger.info("*****************************")
                logger.info("***** Write output data *****")
                logger.info("*****************************")
                logger.info("")
                
                # 3 - Write output data
                logger.info("> Step 3 - Merge temporary LakeSP_Obs|_Prior|_Unassigned shapefiles into single output LakeSP_Obs|_Prior|_Unassigned shapefile")
                self._write_output_data()
            
        except service_error.SwotError:
            raise
            
        except Exception:
            message = "Fatal error catch in PGE LakeSP"
            logger.error(message, exc_info=True)
            raise
        
    def stop(self):
        """
        pseudoPGE method stop
        Close all services
        """
        logger = logging.getLogger(self.__class__.__name__)

        logger.info("")
        logger.info("******************************")
        logger.info("***** Close all services *****")
        logger.info("******************************")
        logger.info("")

        logger.info("")
        logger.info(self.timer.stop())
        logger.sigmsg("==================================")
        logger.sigmsg("===== LakeSPProcessing = END =====")
        logger.sigmsg("==================================")

        # 2 - Stop logger
        instance_logger = service_logger.get_instance()
        if instance_logger is not None:
            instance_logger.flush_log()
            instance_logger.close()

        # 3 - Clean configuration file
        service_config_file.clear_config()

        # 4 - Concatenate log files
        filenames = [self.logfile_begin] + list(self.logfile_swath.values()) + [self.logfile_end]
        with open(self.logfile, 'w') as out_logfile:
            for fname in filenames:
                if os.path.exists(fname):
                    print("Merge logfile %s into %s" %(fname, self.logfile))
                    with open(fname) as infile:
                        for line in infile:
                            out_logfile.write(line)
                    if self.flag_del_tmp_shp:
                        os.remove(fname)
    

#######################################


if __name__ == '__main__':

    # 0 - Parse inline parameters
    PARSER = argparse.ArgumentParser(description="Compute SWOT LakeSP product\
            from all LakeTile products related to the same cycle and pass. ")
    PARSER.add_argument("command_file", help="command file (*.cfg)")
    ARGS = PARSER.parse_args()
    
    PGE = None
    
    try:
        # 1 - Instantiate PGE
        PGE = PGELakeSP(ARGS.command_file)

        # 2 - Start PGE LakeSP
        PGE.start()
        
    finally:
        if PGE is not None:
            # 3 - Stop PGE LakeSP
            PGE.stop()
