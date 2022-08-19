#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:3.2.0:DM:#91:2021/10/27:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: pge_lake_tile.py
    :synopsis: Process PGE_L2_HR_LakeAvg, i.e. generate L2_HR_LakeAvg
    product from one cycle of L2_HR_LakeSP products over one basin 

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
from lxml import etree as ET
import multiprocessing
import os
import sys

import cnes.common.service_config_file as service_config_file
import cnes.common.service_error as service_error
import cnes.common.service_logger as service_logger

import cnes.common.lib.my_timer as my_timer
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_variables as my_var
import cnes.common.lib_lake.lake_db as lake_db
import cnes.common.lib_lake.locnes_filenames as locnes_filenames

import cnes.sas.lake_avg.proc_lakeavg as proc_lakeavg
import cnes.sas.lake_avg.proc_lakesp as proc_lakesp


#######################################


class PGELakeAvg():
    """
    Class PGELakeAvg
    Pseudo PGE class to launch LakeAvg SAS processing
    """

    def __init__(self, cmd_file):
        """
        Constructor of PGELakeAvg

        :param cmd_file: command file full path
        :type cmd_file: string
        """

        # 0 - Init timer
        self.timer = my_timer.Timer()
        self.timer.start()
        
        # 1 - Load command file
        self.cmd_file = cmd_file
        my_tools.test_file(cmd_file, in_extent=".cfg")  # Test existance and extension
        my_params = self._read_cmd_file()  # Read parameters

        self.lakesp_dir = my_params["lakesp_dir"]
        self.output_dir = my_params["output_dir"]

        self.cycle_num = my_params["cycle_num"]
        self.basin_code = my_params["basin_code"]
        self.continent_id = lake_db.compute_continent_id_from_basin_code(self.basin_code[0])

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
        self.flag_inc_file_counter = self.cfg.getboolean("OPTIONS", "Increment file counter")
        self.flag_write_full_path = self.cfg.getboolean("OPTIONS", "Write full path")
        
        # 4 - Set log filenames
        self.logfile = self.cfg.get('LOGGING', 'logFile')
        self.logfile_begin = self.logfile.replace(".log", "_BEGIN.log")
        self.logfile_end = self.logfile.replace(".log", "_END.log")
        self.cfg.set("LOGGING", "logFile", self.logfile_begin)

        # 5 - Init other class properties
        self.list_lakesp_prior = dict()  # List of _Prior files, organized per basin_id
        self.sum_task = 0  # Number of total LakeSP_Prior files to compute
        self.nb_proc = int(my_params["nb_proc"])  # Number of cpu used 
        self.list_tmp_lakeavg_shp = list()  # List of temporary LakeAvg shapefiles (to be later combined)
        self.lake_avg_filename = locnes_filenames.compute_lakeavg_filename(self.output_dir,
                                                                           self.cycle_num, 
                                                                           self.continent_id, 
                                                                           self.basin_code)

        # 6 - Initiate logging service
        service_logger.ServiceLogger()
        logger = logging.getLogger(self.__class__.__name__)

        # 7 - Print info
        logger.sigmsg("=====================================")
        logger.sigmsg("===== LakeAvgProcessing = BEGIN =====")
        logger.sigmsg("=====================================")
        logger.info("> INPUT LakeSP_Prior directory = " + str(self.cfg.get("PATHS", "LakeSP_Prior directory")))
        logger.info("> OUTPUT directory = " + str(self.cfg.get("PATHS", "Output directory")))
        logger.info("> OUTPUT filename = " + str(self.lake_avg_filename))
        logger.info("======================================")
        message = "> Command file: " + str(self.cmd_file)
        logger.info(message)
        message = "> " + str(self.cfg)
        logger.info(message)
        logger.info("======================================")
        
        # 8 - Test input parameters
        logger.info(">> Check config parameters type and value")
        self._check_config_parameters()

        # 9 - Form processing metadata dictionary
        self.proc_metadata = {}
        self.proc_metadata["history"] = "%sZ: Creation" % my_tools.swot_timeformat(datetime.datetime.utcnow(), in_format=3)
        self.proc_metadata["cycle_number"] = "%03d" % self.cycle_num
        self.proc_metadata["continent_id"] = self.continent_id
        self.proc_metadata["continent_code"] = lake_db.compute_continent_code(self.continent_id)
        self.proc_metadata["basin_code"] = self.basin_code
        self.proc_metadata["time_granule_start"] = "None"
        self.proc_metadata["time_granule_end"] = "None"
        geospatial_keys = ["lon_min", "lon_max", "lat_min", "lat_max"]
        for key in geospatial_keys:
            self.proc_metadata["geospatial_%s" % key] = "None"
        if self.cfg.get("DATABASES", "LAKE_DB") is not None:
            self.proc_metadata["xref_prior_lake_db_file"] = self.cfg.get("DATABASES", "LAKE_DB")
        else:
            self.proc_metadata["xref_prior_lake_db_file"] = ""
        self.proc_metadata["xref_param_l2_hr_lakeavg_file"] = file_config
        
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
        out_params["flag_inc_file_counter"] = "True"
        out_params["flag_write_full_path"] = "False"

        # 1 - Read parameter file
        config = configparser.ConfigParser()
        config.read(self.cmd_file)

        # 2 - Retrieve PATHS
        try:
            out_params["param_file"] = os.path.expandvars(config.get("PATHS", "param_file"))
        except:
            out_params["param_file"] = os.path.join(sys.path[0], "lake_avg_param.cfg")
            #out_params["param_file"] = os.path.join(os.getcwd(), "lake_avg_param.cfg")
        out_params["lakesp_dir"] = config.get("PATHS", "LakeSP_Prior directory")
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

        # 4 - Retrieve cycle number and basin code
        out_params["cycle_num"] = config.getint("PASS_INFOS", "Cycle number")
        out_params["basin_code"] = config.get("PASS_INFOS", "Basin code")

        # 5 - Retrieve OPTIONS
        if "OPTIONS" in config.sections():
            list_options = config.options("OPTIONS")
            # Flag to increment the file counter in the output filename (=True, default); else=False
            if "increment file counter" in list_options:
                out_params["flag_inc_file_counter"] = config.get("OPTIONS", "Increment file counter")
            # To write full path in global attributes (=True); to write only basename=False (default)
            if "write full path" in list_options:
                out_params["flag_write_full_path"] = config.get("OPTIONS", "Write full path")
            # Number of processors to use (default=1)
            if "nb_proc" in list_options:
                out_params["nb_proc"] = config.get("OPTIONS", "Nb_proc")
            else:
                out_params["nb_proc"] = 1

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
        # Composite Release IDentifier for LakeSP processing
        out_params["CRID_LAKESP"] = config.get("FILE_INFORMATION", "CRID_LAKESP")
        # Composite Release IDentifier for LakeAvg processing
        out_params["CRID_LAKEAVG"] = config.get("FILE_INFORMATION", "CRID_LAKEAVG")
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
            self.cfg.set(section, "LakeSP_Prior directory", param_list["lakesp_dir"])
            self.cfg.set(section, "Output directory", param_list["output_dir"])
            
            # Add DATABASES section and parameters
            section = "DATABASES"
            self.cfg.add_section(section)
            self.cfg.set(section, "LAKE_DB", param_list["LAKE_DB"])
            self.cfg.set(section, "LAKE_DB_ID", param_list["LAKE_DB_ID"])

            # Add TILES INFOS section and parameters
            section = "PASS_INFOS"
            self.cfg.add_section(section)
            self.cfg.set(section, "Cycle number", param_list["cycle_num"])
            self.cfg.set(section, "Basin code", param_list["basin_code"])

            # Add OPTIONS section and parameters
            section = "OPTIONS"
            self.cfg.add_section(section)
            self.cfg.set(section, "Increment file counter", param_list["flag_inc_file_counter"])
            self.cfg.set(section, "Write full path", param_list["flag_write_full_path"])
            self.cfg.set(section, "Nb_proc", param_list["nb_proc"])
            
            # Add section FILE_INFORMATION
            section = "FILE_INFORMATION"
            self.cfg.add_section(section)
            self.cfg.set(section, "INSTITUTION", param_list["INSTITUTION"])
            self.cfg.set(section, "PRODUCT_VERSION", param_list["PRODUCT_VERSION"])
            self.cfg.set(section, "CRID_LAKESP", param_list["CRID_LAKESP"])
            self.cfg.set(section, "CRID_LAKEAVG", param_list["CRID_LAKEAVG"])
            self.cfg.set(section, "PGE_VERSION", param_list["PGE_VERSION"])
            self.cfg.set(section, "CONTACT", param_list["CONTACT"])

        except Exception:
            print("Something wrong happened in _put_cmd_value !")
            raise

    def _check_config_parameters(self):
        """
        Check parameters coherence for LakeAvg parameter file
        
        :return: True if OK
        :rtype: boolean
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        try:

            # 1 - Config parameters from command file
            
            # 1.1 - PATH section
            # LakeSP_Prior directory
            self.cfg.test_var_config_file('PATHS', 'LakeSP_Prior directory', str, logger=logger)
            my_tools.test_dir(self.cfg.get('PATHS', 'LakeSP_Prior directory'))
            logger.debug('OK - LakeSP_Prior directory = ' + str(self.cfg.get('PATHS', 'LakeSP_Prior directory')))
            # Output directory
            self.cfg.test_var_config_file('PATHS', 'Output directory', str, logger=logger)
            my_tools.test_dir(self.cfg.get('PATHS', 'Output directory'))
            logger.debug('OK - Output directory = ' + str(self.cfg.get('PATHS', 'Output directory')))

            # 1.2 - DATABASES section
            # Lake database full path
            if self.cfg.get('DATABASES', 'LAKE_DB') is None:
                logger.debug('WARNING - LAKE_DB not filled => LakeAvg product not linked to a lake database')
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
                logger.debug('WARNING - Unknown LAKE_DB file format for file: %s => LakeAvg product not linked to a lake database' % self.cfg.get('DATABASES', 'LAKE_DB'))
                
            # 1.3 - TILES_INFO section
            # Cycle number
            self.cfg.test_var_config_file('PASS_INFOS', 'Cycle number', int)
            logger.debug('OK - Cycle number = ' + str(self.cfg.get('PASS_INFOS', 'Cycle number')))
            # Basin code
            self.cfg.test_var_config_file('PASS_INFOS', 'Basin code', int)
            logger.debug('OK - Basin code = ' + str(self.cfg.get('PASS_INFOS', 'Basin code')))
            
            # 1.4 - OPTIONS section
            # Increment output file counter
            self.cfg.test_var_config_file('OPTIONS', 'Increment file counter', bool, logger=logger)
            logger.debug('OK - Increment file counter = ' + str(self.cfg.get('OPTIONS', 'Increment file counter')))
            # Write full path in global attributes
            self.cfg.test_var_config_file('OPTIONS', 'Write full path', bool, logger=logger)
            logger.debug('OK - Write full path = ' + str(self.cfg.get('OPTIONS', 'Write full path')))
            # Number of processors to use 
            self.cfg.test_var_config_file('OPTIONS', 'Nb_proc', int, logger=logger)
            logger.debug('OK - Nb_proc = ' + str(self.cfg.get('OPTIONS', 'Nb_proc')))
            
            # 1.5 - FILE_INFORMATION section
            # Name of producing agency
            self.cfg.test_var_config_file('FILE_INFORMATION', 'INSTITUTION', str, logger=logger)
            logger.debug('OK - INSTITUTION = ' + str(self.cfg.get('FILE_INFORMATION', 'INSTITUTION')))
            # Product version
            self.cfg.test_var_config_file('FILE_INFORMATION', 'PRODUCT_VERSION', str, logger=logger)
            logger.debug('OK - PRODUCT_VERSION = ' + str(self.cfg.get('FILE_INFORMATION', 'PRODUCT_VERSION')))
            # Composite Release IDentifier for LakeSP processing
            self.cfg.test_var_config_file('FILE_INFORMATION', 'CRID_LAKESP', str, logger=logger)
            logger.debug('OK - CRID_LAKESP = ' + str(self.cfg.get('FILE_INFORMATION', 'CRID_LAKESP')))
            # Composite Release IDentifier for LakeAvg processing
            self.cfg.test_var_config_file('FILE_INFORMATION', 'CRID_LAKEAVG', str, logger=logger)
            logger.debug('OK - CRID_LAKEAVG = ' + str(self.cfg.get('FILE_INFORMATION', 'CRID_LAKEAVG')))
            # Version identifier of the product generation executable (PGE)
            self.cfg.test_var_config_file('FILE_INFORMATION', 'PGE_VERSION', str, logger=logger)
            logger.debug('OK - PGE_VERSION = ' + str(self.cfg.get('FILE_INFORMATION', 'PGE_VERSION')))
            # Contact
            self.cfg.test_var_config_file('FILE_INFORMATION', 'CONTACT', str, logger=logger)
            logger.debug('OK - CONTACT = ' + str(self.cfg.get('FILE_INFORMATION', 'CONTACT')))

            # 2 - Config parameters from parameter file

            # 2.1 - CONFIG_PARAMS section
            
            # To add not-observed PLD lakes as empty features to the LakeAvg product (=True, default) or not (=False)
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'ADD_ALL', bool, val_default=True, logger=logger)
            logger.debug('OK - ADD_ALL = ' + str(self.cfg.get('CONFIG_PARAMS', 'ADD_ALL')))
            
            # To delete temporary sub-basin LakeAvg shapefiles (=True, default) or not (=False)
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'DEL_TMP_SHP', bool, val_default=True, logger=logger)
            logger.debug('OK - DEL_TMP_SHP = ' + str(self.cfg.get('CONFIG_PARAMS', 'DEL_TMP_SHP')))

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

    def _select_input_data(self):
        """
        This method selects appropriate LakeSP_Prior shapefiles, and classify them depending on the basin_code related to them
        """
        logger = logging.getLogger(self.__class__.__name__)

        # 1 - List all files in LakeSP_Prior directory
        tmp_list_lakesp_files = os.listdir(self.lakesp_dir)

        # 2 - Compute LakeSP_Prior prefix regarding cycle number and basin_code
        cond_prefix = locnes_filenames.LAKE_SP_PREFIX["prior"]  # Generic LakeSP_Prior prefix 
        cond_prefix += "%03d" % self.cycle_num # Add cycle number to prefix
        logger.debug("Prefix to select LakeSP_Prior shapefiles = %s" % cond_prefix)
        
        # 3 - Compute continent pattern condition wrt specified continent_id
        cond_continent = "_%s_" % self.continent_id

        # 4 - Select LakeSP_Prior shapefiles corresponding to specified cycle number and basin code 
        
        # Init variables
        cpt_input_files = 0  # Input files counter
        self.proc_metadata["source"] = None
        self.proc_metadata["time_coverage_start"] = None
        self.proc_metadata["time_coverage_end"] = None
        list_files_ok = list()
        
        for cur_file in tmp_list_lakesp_files:  # Foreach file in LakeSP directory
            
            flag_file_ok = False

            # 4.1 - Test if file meets the condition of LakeSP_Prior file
            if cur_file.startswith(cond_prefix) and cur_file.endswith(locnes_filenames.LAKE_SP_SHP_SUFFIX) and (cond_continent in cur_file):
                logger.debug("Working with current LakeSP_Prior shapefile: %s" % cur_file)

                # 4.2 - Get basin_code list from LakeSP_Prior metadata file (ie .shp.xml)
                metadata = ET.parse(os.path.join(self.lakesp_dir, cur_file + ".xml"))
                # Retrieve basin_code 
                cur_basin_code_txt = metadata.xpath("//swot_product/global_metadata/basin_code")[0].text
                if not cur_basin_code_txt:
                    cur_basin_code_list = [""]
                else:
                    cur_basin_code_list = cur_basin_code_txt.split(";")
                    
                # 4.3 - Add LakeSP_Prior filename in list of input files for each selected basin_code
                for cur_basin_code in cur_basin_code_list:
                                        
                    if cur_basin_code.startswith(self.basin_code):
                        
                        if cur_basin_code not in self.list_lakesp_prior.keys():
                            logger.debug("basin_code=%s UNKNOWN => add to dictionary" % cur_basin_code)
                            self.list_lakesp_prior[cur_basin_code] = set()
                            
                        logger.debug("basin_code=%s -> Add %s shapefile" % (cur_basin_code, cur_file))
                        self.list_lakesp_prior[cur_basin_code].add(os.path.join(self.lakesp_dir, cur_file))
                        
                        # Retrieve source global attribute
                        cur_source_txt = metadata.xpath("//swot_product/global_metadata/source")[0].text
                        if self.proc_metadata["source"] is None:
                            self.proc_metadata["source"] = cur_source_txt
                        else:
                            if cur_source_txt != self.proc_metadata["source"]:
                                message = "ERROR bad source. In XML file: " + str(cur_source_txt) + " WHEREAS previously set at: " \
                                            + str(self.proc_metadata["source"])
                                raise service_error.ProcessingError(message, logger)
                        
                        # Compare used database for input LakeSP_Prior file
                        cur_lakedb_txt = metadata.xpath("//swot_product/global_metadata/xref_prior_lake_db_file")[0].text
                        if cur_lakedb_txt != self.proc_metadata["xref_prior_lake_db_file"]:
                            message = "ERROR bad xref_prior_lake_db_file. In XML file: " + str(cur_lakedb_txt) + " WHEREAS in LakeAvg command file: " \
                                            + str(self.proc_metadata["xref_prior_lake_db_file"])
                            logger.warning(message)
                            #raise service_error.ProcessingError(message, logger)
                            
                        # Retrieve time_coverage_start and _end, and update values for current LakeAvg product
                        # Time start
                        cur_time_start = metadata.xpath("//swot_product/global_metadata/time_coverage_start")[0].text
                        if self.proc_metadata["time_coverage_start"] is None:
                            self.proc_metadata["time_coverage_start"] = cur_time_start
                        else:
                            if cur_time_start and self.proc_metadata["time_coverage_start"] > cur_time_start:
                                self.proc_metadata["time_coverage_start"] = cur_time_start
                        # Time stop
                        cur_time_stop = metadata.xpath("//swot_product/global_metadata/time_coverage_end")[0].text
                        if self.proc_metadata["time_coverage_end"] is None:
                            self.proc_metadata["time_coverage_end"] = cur_time_stop
                        else:
                            if cur_time_stop and self.proc_metadata["time_coverage_end"] < cur_time_stop:
                                self.proc_metadata["time_coverage_end"] = cur_time_stop
                        
                        flag_file_ok = True
                        
                    else:
                        logger.debug("basin_code=%s NOT SELECTED" % cur_basin_code)
                        
                if flag_file_ok:
                    list_files_ok.append(cur_file)
                    cpt_input_files += 1
                    
        # 5 - Sort list of files per basin_code
        for cur_basin_code in self.list_lakesp_prior.keys():
            self.list_lakesp_prior[cur_basin_code] = sorted(self.list_lakesp_prior[cur_basin_code])
            self.sum_task += len(self.list_lakesp_prior[cur_basin_code])
        
        # 6 - Set global attribute for list of selected LakeSP_Prior files
        list_files_ok_sorted = sorted(list_files_ok)
        self.proc_metadata["xref_l2_hr_lakesp_files"] = my_var.SEP_FILES.join(list_files_ok_sorted)
        
        # Print some numbers
        logger.debug("Number of selected input files = %d" % cpt_input_files)
        logger.debug("Number of level-3 basins to process = %d" % len(self.list_lakesp_prior.keys()))

    def _read_lake_db(self, in_basin_code):
        """
        This method prepares Lake DB class, and init with Lake DB data
        
        :param in_basin_code: basin code to select lakes from DB having lake_id starting with in_basin_code
        :type in_basin_code: string
        
        :return: out_obj_lake_db = subset of the PLD 
        :rtype: cnes.common.lib_lake.lake_db.<LakeDb|LakeDbShp|LakeDbSqlite>
        """
        logger = logging.getLogger(self.__class__.__name__)
        timer = my_timer.Timer()
        timer.start()

        lake_db_file = self.cfg.get("DATABASES", "LAKE_DB")
        
        if (lake_db_file == "") or (lake_db_file is None):  # No PLD
            logger.warning("NO database specified")
            out_obj_lake_db = lake_db.LakeDb()
            
        else:  # Init PLD object wrt to file type
            
            type_db = lake_db_file.split('.')[-1]  # File type
            
            if type_db == "shp":  # User lake DB file in shapefile format
                logger.debug("Use of personal user lake database in shapefile format")
                out_obj_lake_db = lake_db.LakeDbShp(lake_db_file)
                
            elif type_db == "sqlite": # User or operational lake DB file in SQLite format
                if os.path.basename(lake_db_file).startswith(locnes_filenames.PLD_PREFIX_BASE):  # Operationnal PLD filename
                    logger.debug("Use of operational Prior Lake Database (specific filename given)")
                    out_obj_lake_db = lake_db.LakeDbSqlite(lake_db_file, 
                                                           in_basin_id=in_basin_code)
                else:
                    logger.debug("Use of personal user lake database in SQLite format")
                    out_obj_lake_db = lake_db.LakeDbSqlite(lake_db_file)
                    
            elif os.path.isdir(lake_db_file):  # Directory containing SQLite PLD files # Directory of operational Prior Lake Database (PLD)
                logger.debug("Use of operational Prior Lake Database (upper directory given)")
                # Compute PLD filename
                lakedb_filenames = locnes_filenames.PldFilenames("AVG", 
                                                                 lake_db_file, 
                                                                 cycle_num=self.cycle_num,
                                                                 basin_code=in_basin_code)
                # Init PLD object
                if os.path.exists(lakedb_filenames.filename):
                    out_obj_lake_db = lake_db.LakeDbSqlite(lakedb_filenames.filename, 
                                                           in_basin_id=in_basin_code)
                    self.proc_metadata["xref_prior_lake_db_file"] = lakedb_filenames.filename
                else:
                    logger.warning("NO related operational Prior Lake Database file found -> NO database specified")
                    self.obj_lake_db = lake_db.LakeDb()
                    
            else:
                message = "Prior Lake Database format (%s) is unknown: must be .shp or .sqlite" % type_db
                logger.error(message, exc_info=True)
                raise
                
            if out_obj_lake_db.basin_coords is not None:
                for key, value in out_obj_lake_db.basin_coords.items():
                    if self.proc_metadata["geospatial_%s" % key] == "None":
                        self.proc_metadata["geospatial_%s" % key] = value
                    elif ("min" in key) and (value < self.proc_metadata["geospatial_%s" % key]):
                            self.proc_metadata["geospatial_%s" % key] = value
                    elif ("max" in key) and (value > self.proc_metadata["geospatial_%s" % key]):
                            self.proc_metadata["geospatial_%s" % key] = value

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

        # 1 - Write LakeAvg shapefile
        logger.info("> Combining %d LakeAvg temporary shapefiles to output LakeAvg shapefile = %s" % (len(self.list_tmp_lakeavg_shp), self.lake_avg_filename))
        proc_lakeavg.write_file(self.list_tmp_lakeavg_shp, self.lake_avg_filename, self.proc_metadata)
        logger.debug(timer.stop("LakeSP file writing time: "))
        logger.info("")
        
    # -------------------------------------------
    
    def main_process(self, in_list_lakesp_prior):
        """
        LakeAvg main process
        
        :param in_list_lakesp_prior: dictionnary of LakeSP_Prior files to process; key = basin codes
        :type in_list_lakesp_prior: dict
        
        :return: out_list_tmp_lakeavg_shp = list of temporary output LakeAvg shapefiles; 
                    there is a one-to-one correspondance with the basin codes above
        :rtype: out_list_tmp_lakeavg_shp = list
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 0.1 - Init list of temporary LakeAvg shapefiles
        out_list_tmp_lakeavg_shp = list()
        
        for cur_basin_code, cur_list_lakesp_files in in_list_lakesp_prior.items():            
            timer = my_timer.Timer()
            timer.start()

            logger.info("==========================")
            logger.info("== Processing basin %03d ==" % int(cur_basin_code))
            logger.info("==========================")
            logger.info("")
        
            # 0.2 - Compute temporary filename
            tmp_filename = os.path.join(self.output_dir, "%s_%s.shp" % (os.path.splitext(os.path.basename(self.lake_avg_filename))[0], cur_basin_code))
            out_list_tmp_lakeavg_shp.append(tmp_filename)
        
            # 1 - Read lake DB and filter over cur_basin_code
            logger.info("> 1 - Init and read lake DB object over basin...")
            obj_lakedb = self._read_lake_db(in_basin_code=cur_basin_code)
            if obj_lakedb.lake_layer is not None:
                logger.debug("%d PLD lakes over basin" % obj_lakedb.lake_layer.GetFeatureCount())
            logger.info("")
        
            # 2 - Read input LakeSP_Prior shapefiles
            logger.info("> 2 - Read input files...")
            obj_lakesp = proc_lakesp.LakeSPProduct(cur_basin_code)
            obj_lakesp.set_from_lakesp_files(self.list_lakesp_prior[cur_basin_code])
            logger.debug("%d PLD lakes observed during the cycle on this basin" % len(obj_lakesp.lakesp_archive.keys()))
            logger.info("")
    
            # 3 - Prepare output data
            logger.info("> 3 - Initialize LakeAvg object......")
            obj_lakeavg = proc_lakeavg.LakeAvgProduct(obj_lakedb, obj_lakesp,
                                                  os.path.splitext(os.path.basename(tmp_filename))[0])
            logger.info("")
        
            # 4 - Compute LakeAvg features
            logger.info("> 4 - Compute LakeAvg features...")
            obj_lakeavg.compute_lake_features()
            logger.info("")
        
            # 5 - Write temporary output file
            logger.info("> 5 - Write temporary LakeAvg file...")
            obj_lakeavg.write_tmp_file(tmp_filename)
            logger.info("")
        
            # 6 - Close lake database
            logger.info("> 6 - Closing lake database...")
            obj_lakedb.close_db()
        
            # 8 - Free memory
            del obj_lakesp
            obj_lakeavg.free_memory()
            del obj_lakeavg

            logger.debug(timer.stop("Sub basin %s computation time: " %cur_basin_code))
            logger.info("")
            
        return out_list_tmp_lakeavg_shp

    def start(self):
        """
        Main pseudoPGE method
        Start computation
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        ##################################
        # Definition of process function #
        ##################################
        def process_function(in_list_lakesp_prior, queue):
            """
            Process function

            :param in_list_lakesp_prior: dictionnary of LakeSP_Prior files to process; key = basin codes
            :type in_list_lakesp_prior: dict
            :param queue: process queue
            :type queue: multiprocessing.Queue
            """
            proc_pid = multiprocessing.current_process().pid

            # 1 - Open new log file for each subproc
            cur_logfile_subproc = self.logfile.replace(".log", "_subproc%d.log" %proc_pid)

            self.cfg.set("LOGGING", "logFile", cur_logfile_subproc)
            service_logger.ServiceLogger()
            logger = logging.getLogger(self.__class__.__name__)
            logger.debug("=============================")
            logger.debug("Beginning of log is written to file %s" % self.logfile_begin)
            logger.debug("=============================")
            logger.debug("Log for basin codes %s will be written to log %s" % (";".join(in_list_lakesp_prior.keys()), cur_logfile_subproc))
            logger.debug("=============================")
            logger.debug("")

            # 2 - Run main process
            try:
                list_tmp_lakeavg_shp = self.main_process(in_list_lakesp_prior)
                queue.put(list_tmp_lakeavg_shp)
                logger.debug("STOP process")
                
            except:
                proc_pid = multiprocessing.current_process().pid
                proc_name = multiprocessing.current_process().name
                logger.error("Process " + proc_name + " with pid " + str(proc_pid) + " failed", exc_info=True)
                queue.put(None)
                raise

            # 3 - Close log file for each subproc
            logger.debug("Closing temporary log file %s" % cur_logfile_subproc)
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
            logger.info("> Step 1 - Select and classify input data per Level-3 basin identifier")
            self._select_input_data()
            logger.info("")

            if self.nb_proc == 1:
                # 2 - Run main process directly
                logger.info("> Step 2 - Main processing using ONLY 1 processor")
                logger.info("")
                self.list_tmp_lakeavg_shp = self.main_process(self.list_lakesp_prior)
                
            else:
                # 2 - Run main process after spreading over available processors
                logger.info("> Step 2 - Main processing using %d processors" % self.nb_proc)
                logger.info("")
                
                logger.info("*********************************************************")
                logger.info("***** Spread list of basin_code to process per proc *****")
                logger.info("*********************************************************")
                logger.info("")
                
                # 2.1 - Spread basin_code list per processor
                logger.info("> 2.1 - Spread basin_code list per processor")
            
                # 2.1.1 - Init values
                nb_proc = int(self.nb_proc)
                nb_task_per_proc = int(self.sum_task/nb_proc)
                message = "sum_nb_task: " + str(self.sum_task)
                logger.debug(message)
                message = "nb_task_per_proc: " + str(nb_task_per_proc)
                logger.debug(message)
    
                # 2.1.2 - Calculates the repartition of basin_id in function of number of LakeSP_Prior files to read
                # Generate lists of basin_id with the corresponding number of LakeSP_prior files to process
                list_of_basin_id = []
                list_of_nb_file = []
                for cur_basin_id in self.list_lakesp_prior.keys():
                    list_of_basin_id.append(cur_basin_id)
                    list_of_nb_file.append(len(self.list_lakesp_prior[cur_basin_id]))
                # Sort list in function of number of LakeSP_Prior files
                indice_list_of_nb_file = sorted(range(len(list_of_nb_file)), key=list_of_nb_file.__getitem__, reverse=True)
                list_of_basin_id_sorted = []
                for i in indice_list_of_nb_file:
                    list_of_basin_id_sorted.append(list_of_basin_id[i])
                # Generate list of basin_id in function of number of LakeSP_Prior files
                tmp_list_process_bassin_id = []
                list_process_bassin_id = []
                tmp_nb_task_per_proc = 0
                proc_i = 0
                for cur_basin_id in list_of_basin_id_sorted:
                    tmp_list_process_bassin_id.append(cur_basin_id)
                    # Calculate the number of task for this proc
                    tmp_nb_task_per_proc = tmp_nb_task_per_proc + len(self.list_lakesp_prior[cur_basin_id])
                    # If the number of task per proc is greater than the mean value, stop to append
                    if (tmp_nb_task_per_proc > nb_task_per_proc):
                        list_process_bassin_id.append(tmp_list_process_bassin_id)
                        proc_i = proc_i + 1
                        tmp_nb_task_per_proc = 0
                        tmp_list_process_bassin_id = []
                # Manage the last proc
                if (tmp_nb_task_per_proc > 0):
                    list_process_bassin_id.append(tmp_list_process_bassin_id)
                    proc_i = proc_i + 1
    
                # 2.1.3 - Update of number of processors in function of input data
                nb_proc = proc_i
                # If number of processor change log info
                if (nb_proc != self.nb_proc):
                    message = "The number of proc is reduced: from " + str(self.nb_proc) + " to " + str(nb_proc)
                    logger.debug(message)
                else:
                    message = "nb_proc: " + str(nb_proc)
                    logger.debug(message)
                message = "list_process_bassin_id: " + str(list_process_bassin_id)
                logger.debug(message)
    
                # 2.1.4 - Creation of list of LakeSP_Prior files per proc
                select_process = []
                for proc_i in range(nb_proc):
                    tmp_select_process = dict()
                    for cur_basin_id in list_process_bassin_id[proc_i]:
                        tmp_select_process[cur_basin_id] = self.list_lakesp_prior[cur_basin_id]
                    select_process.append(tmp_select_process)
                    
                logger.info("")

                logger.info("****************************************")
                logger.info("***** Prepare and launch processes *****")
                logger.info("****************************************")
                logger.info("")
     
                list_process = []
                list_queue = []
    
                # 2.2 - Prepare processes
                logger.info("> 2.2 - Prepare processes")
                for proc_i in range(nb_proc):
                    # 2.2.1 - Select list of LakeSP_Prior files to process, defined by their related basin_code
                    cur_list_lakesp_prior = select_process[proc_i]
                    logger.debug("basin_id for process nb " + str(proc_i) + ": ")
                    logger.debug(cur_list_lakesp_prior)
                    # 2.2.2 - Create Queue
                    queue = multiprocessing.Queue()
                    list_queue.append(queue)
                    # 2.2.3 - Create Process
                    nom = "sub_process_" + str(proc_i)
                    p = multiprocessing.Process(name=nom, target=process_function, args=(cur_list_lakesp_prior, queue))
                    list_process.append(p)
                logger.info("")
    
                # 2.3 - Launch processes
                logger.info("> 2.3 - Launch processes")
                proc_i = 0

                # 2.3.0 - Write log to file before running multiprocessing
                logger.debug("Closing log %s to write a log for each subproc" % self.logfile_begin)
                instance_logger = service_logger.get_instance()
                if instance_logger is not None:
                    instance_logger.flush_log()
                    instance_logger.close()

                # 2.3.1 - Run subproc
                for p in list_process:
                    logger.debug("Start process " + str(proc_i))
                    proc_i = proc_i + 1
                    p.start()
                logger.info("")

                # 2.3.2 -  Reopen principal log file
                self.cfg.set("LOGGING", "logFile", self.logfile_end)
                service_logger.ServiceLogger()
                logger = logging.getLogger(self.__class__.__name__)
                logger.debug("Open log file %s" % self.logfile_end)
                logger.info("")

                # 2.4 - Get results from processes
                logger.info("> 2.4 - Get results")
                for q in list_queue:
                    list_tmp_lakeavg_shp = q.get()
                    if list_tmp_lakeavg_shp is not None:
                        for tmp_lakeavg_shp in list_tmp_lakeavg_shp:
                            # Add temporary LakeAvg filenames to list
                            self.list_tmp_lakeavg_shp.append(tmp_lakeavg_shp)
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
            logger.info("> Step 3 - Merge temporary LakeAvg shapefiles into a single output LakeAvg shapefile")
            self._write_output_data()
            
        except service_error.SwotError:
            raise
            
        except Exception:
            message = "Fatal error catch in PGE LakeAvg"
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
        logger.sigmsg("===================================")
        logger.sigmsg("===== LakeAvgProcessing = END =====")
        logger.sigmsg("===================================")

        # 2 - Stop logger
        instance_logger = service_logger.get_instance()
        if instance_logger is not None:
            instance_logger.flush_log()
            instance_logger.close()
            
        # 3 - Clean configuration file
        service_config_file.clear_config()

        # 4 - Concatenate log files
        # 4.1 - Retrieve logfiles from subproc
        files = os.listdir(os.path.dirname(self.logfile))
        logfile_name = os.path.basename(self.logfile)
        logfile_subproc = [os.path.join(os.path.dirname(self.logfile), file) for file in files if file.startswith(logfile_name.replace(".log", "_subproc"))]

        # 4.2 - Concatenate logfile
        filenames = [self.logfile_begin] + logfile_subproc + [self.logfile_end]
        with open(self.logfile, 'w') as out_logfile:
            for fname in filenames:
                if os.path.exists(fname):
                    print("Merge logfile %s into %s" %(fname, self.logfile))
                    with open(fname) as infile:
                        for line in infile:
                            out_logfile.write(line)

        # 4.3 - Delete temporary files
        if self.cfg.getboolean("CONFIG_PARAMS", "DEL_TMP_SHP"):
            for file in filenames:
                if os.path.exists(file):
                    os.remove(file)
    

#######################################


if __name__ == '__main__':

    # 0 - Parse inline parameters
    PARSER = argparse.ArgumentParser(description="Compute SWOT LakeAvg product\
            from one cycle of L2_HR_LakeSP products over one basin. ")
    PARSER.add_argument("command_file", help="command file (*.cfg)")
    ARGS = PARSER.parse_args()
    
    PGE = None
    
    try:
        # 1 - Instantiate PGE
        PGE = PGELakeAvg(ARGS.command_file)

        # 2 - Start PGE LakeAvg
        PGE.start()
        
    finally:
        if PGE is not None:
            # 3 - Stop PGE LakeAvg
            PGE.stop()
    
