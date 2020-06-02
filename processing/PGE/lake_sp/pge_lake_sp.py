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
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: pge_lake_sp.py
    :synopsis: Process PGE_L2_HR_LakeSP, i.e. generate L2_HR_LakeSP
    product from all tile of L2_HR_LAKE_TILE product
    Created on 27/09/2017

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR
                  Cécile Cazals - CS

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
import os
import sys

import cnes.common.service_config_file as service_config_file
import cnes.common.service_error as service_error
import cnes.common.service_logger as service_logger

import cnes.common.lib.my_shp_file as my_shp
import cnes.common.lib.my_timer as my_timer
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib_lake.lake_db as lake_db
import cnes.common.lib_lake.locnes_filenames as locnes_filenames
import cnes.common.lib_lake.locnes_products_shapefile as locnes_products_shapefile
import cnes.common.lib_lake.proc_lake as proc_lake

import cnes.sas.lake_sp.sas_lake_sp as sas_lake_sp
import cnes.sas.lake_sp.proc_pixc_sp as proc_pixc_sp
import cnes.sas.lake_sp.proc_pixc_vec_sp as proc_pixc_vec_sp


#######################################


class PGELakeSP():
    """
    Class PGELakeSP
    Pseudo PGE class to launch SAS LakeSP computation
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

        self.laketile_shp_dir = my_params["laketile_shp_directory"]
        self.laketile_edge_dir = my_params["laketile_edge_directory"]
        self.laketile_pixcvec_dir = my_params["laketile_pixcvec_directory"]
        self.output_dir = my_params["output_dir"]

        self.cycle_num = my_params["cycle_num"]
        self.pass_num = my_params["pass_num"]
        self.continent_list = set()

        # 2 - Load parameter file
        # 2.1 - Read value from command file
        self.file_config = my_params["param_file"]

        # 2.2 - Test existence
        if not os.path.exists(self.file_config):
            raise service_error.DirFileError(self.file_config)
        # 2.3 - Load parameters
        self.cfg = service_config_file.ServiceConfigFile(self.file_config)

        # 3 - Put command parameter inside cfg
        self._put_cmd_value(my_params)
        
        # 4 - Init data dict. Each key / value will be continent_code : object
        self.obj_pixc_sp = {}
        self.obj_pixcvec_sp = {}
        self.obj_lake_db = {}
        self.obj_lake_r = {}
        self.obj_lake_l = {}
        self.lake_sp_filenames = {}
        self.laketile_shp_path_dict = {}
        self.flag_prod_shp = self.cfg.getboolean("OPTIONS", "Produce shp")

        # 5 - Initiate logging service
        service_logger.ServiceLogger()
        logger = logging.getLogger(self.__class__.__name__)

        # 6 - Print info
        logger.sigmsg("====================================")
        logger.sigmsg("===== lakeSPProcessing = BEGIN =====")
        logger.sigmsg("====================================")
        message = "> Command file: " + str(self.cmd_file)
        logger.info(message)
        message = "> " + str(self.cfg)
        logger.info(message)
        logger.info("")
        
        # 7 - Test input parameters
        logger.info(">> Test input parameters")
        self._check_config_parameters()
        logger.info("")

        # 8 - Form processing metadata dictionary
        self.proc_metadata = {}

        logger.info("")
        logger.info("")

    def _test_input_directories(self):
        """
        This method tests mandatory input directories
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1.1 - LakeTile_shp directory
        message = "> 1.1 - INPUT LakeTile_shp directory = %s" % self.laketile_shp_dir
        logger.info(message)
        my_tools.test_dir(self.laketile_shp_dir)


        # 1.2 - LakeTile_edge directory
        message = "> 1.2 - INPUT LakeTile_edge directory = %s" % self.laketile_edge_dir
        logger.info(message)
        my_tools.test_dir(self.laketile_edge_dir)

        # 1.3 - LakeTile_pixcvec directory
        message = "> 1.3 - INPUT LakeTile_pixcvec directory = %s" % self.laketile_pixcvec_dir
        logger.info(message)
        my_tools.test_dir(self.laketile_pixcvec_dir)

    def _test_output_directory(self):
        """
        This method tests output directory
        """
        logger = logging.getLogger(self.__class__.__name__)
        message = "  OUTPUT DIR = %s" % self.output_dir
        logger.info(message)
        my_tools.test_dir(self.output_dir)

    def _read_input_data(self):
        """
        This method retrieves needed input data, reads some info in it and prepares associated data classes
        """
        logger = logging.getLogger(self.__class__.__name__)

        # 1 - List all files in LakeTile_shp directory
        tmp_list_laketile_shp = os.listdir(self.laketile_shp_dir)

        # 2 - Compute LakeTile prefix regarding cycle and pass numbers
        cond_prefix = locnes_filenames.LAKE_TILE_PREFIX  # Generic LakeTile 
        cond_prefix += "%03d" % self.cycle_num # Add cycle number to prefix
        cond_prefix += "_%03d" % self.pass_num # Add pass number to prefix

        # 3 - For each listed LakeTile shapefile, get related LakeTile_edge and _pixcvec files if they exist
        
        # Init variables
        laketile_edge_path_dict = {}  # List of _edge files
        laketile_pixcvec_path_dict = {}  # List of _pixcvec files
        pixc_path_dict = {}  # List of PIXC files
        pixcvec_river_path_dict = {}  # List of associated files
        
        for cur_laketile_shp_file in tmp_list_laketile_shp:  # Foreach file in LakeTile shapefile directory

            # 3.1 - Test if file meets the condition of LakeTile_shp file
            if cur_laketile_shp_file.startswith(cond_prefix) and cur_laketile_shp_file.endswith(locnes_filenames.LAKE_TILE_SHP_SUFFIX):

                # Construct LakeTile files filenames
                cur_laketile_shp_path = os.path.join(self.laketile_shp_dir, cur_laketile_shp_file)  # LakeTile_shp filename
                cur_laketile_edge_path = os.path.join(self.laketile_edge_dir,
                                                      cur_laketile_shp_file.replace(locnes_filenames.LAKE_TILE_SHP_SUFFIX,
                                                                       locnes_filenames.LAKE_TILE_EDGE_SUFFIX))  # LakeTile_edge filename
                cur_laketile_pixcvec_path = os.path.join(self.laketile_pixcvec_dir,
                                                         cur_laketile_shp_file.replace(locnes_filenames.LAKE_TILE_SHP_SUFFIX,
                                                                          locnes_filenames.LAKE_TILE_PIXCVEC_SUFFIX))  # LakeTile_pixcvec filename

                # 3.2 - Verify if the LakeTile triplet (_shp, _edge, _pixcvec) exists
                if os.path.exists(cur_laketile_edge_path) and os.path.exists(cur_laketile_pixcvec_path):
                    
                    # Get continent, PIXC and PIXCVecRiver information from LakeTile_shp metadata file (ie .shp.xml)
                    metadata = ET.parse(cur_laketile_shp_path + ".xml")
                    
                    # Test and overwrite configuration parameters with XML file content
                    try:
                        locnes_products_shapefile.set_config_from_xml(metadata)
                    except:
                        message = "ERROR: problem with LakeTile XML file: " + cur_laketile_shp_path + ".xml"
                        raise service_error.ProcessingError(message, logger)

                    # Retrieve continent 
                    cur_continent_txt = metadata.xpath("//swot_product/global_attributes/continent")[0].text
                    if not cur_continent_txt:
                        cur_continent_list = [""]
                    else :
                        cur_continent_list = cur_continent_txt.split(";")
                        
                    # Retrieve LakeTile processing input files (PIXC and PIXCVecRiver)
                    cur_pixc_path = metadata.xpath("//swot_product/global_attributes/xref_input_l2_hr_pixc_file")[0].text
                    cur_pixcvec_river_path = metadata.xpath("//swot_product/global_attributes/xref_input_l2_hr_pixc_vec_river_file")[0].text
                                        
                    # 3.3 - Add current information to lists related to files to process
                    for cur_continent in cur_continent_list:
                        
                        # Current continent
                        self.continent_list.add(cur_continent)

                        # LakeTile_shp full path 
                        self.laketile_shp_path_dict.setdefault(cur_continent, []).append(cur_laketile_shp_path)
                        # LakeTile_edge full path
                        laketile_edge_path_dict.setdefault(cur_continent, []).append(cur_laketile_edge_path)
                        # LakeTile_pixcvec full path
                        laketile_pixcvec_path_dict.setdefault(cur_continent, []).append(cur_laketile_pixcvec_path)

                        # PIXC full path
                        pixc_path_dict.setdefault(cur_continent, []).append(cur_pixc_path)
                        # PIXCVecRiver full path
                        pixcvec_river_path_dict.setdefault(cur_continent, []).append(cur_pixcvec_river_path)

        # 4 - Fill processing metadata with filenames
        for cur_continent in self.continent_list:
            
            # Init
            self.proc_metadata[cur_continent] = {}
            
            # PLD filename
            if not self.cfg.get("DATABASES", "LAKE_DB"):
                self.proc_metadata[cur_continent]["xref_static_lake_db_file"] = ""
            else :
                self.proc_metadata[cur_continent]["xref_static_lake_db_file"] = self.cfg.get("DATABASES", "LAKE_DB")
            # PIXC filename
            self.proc_metadata[cur_continent]["xref_input_l2_hr_pixc_file"] = ", ".join(pixc_path_dict[cur_continent])
            # PIXCVecRiver filename
            self.proc_metadata[cur_continent]["xref_input_l2_hr_pixc_vec_river_file"] = ", ".join(pixcvec_river_path_dict[cur_continent])
            # LakeSP processor parameter file
            self.proc_metadata[cur_continent]["xref_l2_hr_lake_sp_param_file"] = self.file_config

        # 5 - Init PIXC_SP and PIXCVec_SP objects by retrieving data from the LakeTile files for each continent
        logger.info("> 2 - Init PIXC_edge_SP and PIXCVec_SP for each continent")
        for cur_continent in self.continent_list:

            logger.info(". Init objects for continent %s" % cur_continent)

            # Init PIXC_edge_SP object
            self.obj_pixc_sp[cur_continent] = proc_pixc_sp.PixCEdge(laketile_edge_path_dict[cur_continent], 
                                                                    self.cycle_num,
                                                                    self.pass_num, 
                                                                    cur_continent)
            
            # Fill PIXC_edge_SP object with LakeTile_edge file information
            self.obj_pixc_sp[cur_continent].set_pixc_edge_from_laketile_edge_file()

            # Init PIXCVec_SP
            self.obj_pixcvec_sp[cur_continent] = proc_pixc_vec_sp.PixCVecSP(laketile_pixcvec_path_dict[cur_continent], 
                                                                            self.obj_pixc_sp[cur_continent],
                                                                            self.output_dir, 
                                                                            cur_continent)

    def _read_lake_db(self):
        """
        This method retrieves needed PLD info
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Retrieve lake Db layer
        lake_db_file = self.cfg.get("DATABASES", "LAKE_DB")
        
        if (lake_db_file == "") or (lake_db_file is None):  # No PLD
            logger.warning("NO database specified -> NO link of SWOT obs with a priori lake")
            for cur_continent in self.continent_list:
                self.obj_lake_db[cur_continent] = lake_db.LakeDb()
                
        else:  # Init PLD object wrt to file type
            
            type_db = lake_db_file.split('.')[-1]  # Type of file
            
            if type_db == "shp":  # Shapefile format
                for cur_continent in self.continent_list:
                    self.obj_lake_db[cur_continent] = lake_db.LakeDbShp(lake_db_file, self.obj_pixc_sp[cur_continent].tile_poly)
            elif type_db == "sqlite":  # SQLite format
                for cur_continent in self.continent_list:
                    self.obj_lake_db[cur_continent] = lake_db.LakeDbSqlite(lake_db_file, self.obj_pixc_sp[cur_continent].tile_poly)
            else:
                self.obj_lake_db = lake_db.LakeDb()
                logger.warning("Prior Lake Database format (%s) is unknown: must be .shp or .sqlite => NO database used" % type_db)

    def _prepare_output_data(self):
        """
        This method creates output data class
        """
        
        # 4 - Retrieve orbit info from Laketile filename and compute output filenames
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("> 4 - Preparing output objects for each continent...")

        for cur_continent in self.continent_list:
            logger.info("== Continent %s ==" % cur_continent)

            logger.info("Retrieving tile infos from LakeTile filename")
            self.lake_sp_filenames[cur_continent] = locnes_filenames.LakeSPFilenames(self.laketile_shp_path_dict[cur_continent], 
                                                                                     cur_continent, 
                                                                                     self.output_dir)

            logger.info("Init lake product object for swath R")
            self.obj_lake_r[cur_continent] = proc_lake.LakeProduct("SP",
                                                                   self.obj_pixc_sp[cur_continent].pixc_edge_r,
                                                                   self.obj_pixcvec_sp[cur_continent].pixcvec_r,
                                                                   self.obj_lake_db[cur_continent],
                                                                   os.path.basename(self.lake_sp_filenames[cur_continent].lake_sp_file).split(".")[0] + "_R")

            logger.info("Init lake product object for swath L")
            self.obj_lake_l[cur_continent] = proc_lake.LakeProduct("SP",
                                                                   self.obj_pixc_sp[cur_continent].pixc_edge_l,
                                                                   self.obj_pixcvec_sp[cur_continent].pixcvec_l,
                                                                   self.obj_lake_db[cur_continent],
                                                                   os.path.basename(self.lake_sp_filenames[cur_continent].lake_sp_file).split(".")[0] + "_L")

    def _write_output_data(self, in_proc_metadata_list):
        """
        This method write output data
        """
        logger = logging.getLogger(self.__class__.__name__)

        for cur_continent in self.continent_list:
            
            # 1 - Merge shapefiles to get LakeSP product
            logger.info("4 - Merging shapefiles to get LakeSP product %s..." % os.path.basename(self.lake_sp_filenames[cur_continent].lake_sp_file))
            # 1.1 - Merging Right and Left SP layers
            logger.debug("> Merging right and left SP layers...")
            dataSource_sp1, layer_sp = my_shp.merge_2_layers(self.obj_lake_r[cur_continent].shp_mem_layer.layer, self.obj_lake_l[cur_continent].shp_mem_layer.layer, cur_continent)
            # 1.2 - Merge SP layer with shapefiles retrieved from PGE_LakeTile
            logger.debug("> Merging SP layer with shapefiles retrieved from PGE_LakeTile...")
            dataSource_sp2, layer_sp = my_shp.merge_mem_layer_with_shp(self.laketile_shp_path_dict[cur_continent], layer_sp, cur_continent)
            # 1.3 - Write LakeSP shapefile product
            logger.debug("> Writing L2_HR_LakeSP shapefile = %s" % os.path.basename(self.lake_sp_filenames[cur_continent].lake_sp_file))
            my_shp.write_mem_layer_as_shp(layer_sp, self.lake_sp_filenames[cur_continent].lake_sp_file)
            # 1.4 - Write XML metadatafile for LakeSP shapefile product
            if self.obj_pixc_sp[cur_continent].pixc_edge_r.pass_num != 0:
                self.obj_lake_r[cur_continent].shp_mem_layer.update_and_write_metadata("%s.xml" % self.lake_sp_filenames[cur_continent].lake_sp_file,
                                                                                       in_pixc_metadata=self.obj_pixc_sp[cur_continent].pixc_metadata,
                                                                                       in_proc_metadata=in_proc_metadata_list[cur_continent])
            elif self.obj_pixc_sp[cur_continent].pixc_edge_l.pass_num != 0:
                self.obj_lake_l[cur_continent].shp_mem_layer.update_and_write_metadata("%s.xml" % self.lake_sp_filenames[cur_continent].lake_sp_file,
                                                                                       in_pixc_metadata=self.obj_pixc_sp[cur_continent].pixc_metadata,
                                                                                       in_proc_metadata=in_proc_metadata_list[cur_continent])

            # 2 - Close dataSources
            dataSource_sp1.Destroy()
            dataSource_sp2.Destroy()
            self.obj_lake_r[cur_continent].shp_mem_layer.free()
            self.obj_lake_l[cur_continent].shp_mem_layer.free()
            logger.info("")

            # 3 - Update PIXCVec files
            logger.info("5 - Updating L2_HR_PIXCVec files...")
            self.obj_pixcvec_sp[cur_continent].pixcvec_r.update_pixc_vec(self.flag_prod_shp)
            self.obj_pixcvec_sp[cur_continent].pixcvec_l.update_pixc_vec(self.flag_prod_shp)
            logger.info("")

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
            self._test_input_directories()
            self._test_output_directory()

            # 2 - read input data
            logger.info("> 2 - Init and format intput objects...")
            self._read_input_data()

            # 3 - read lake db
            logger.info("> 3 - Read lake DB and init object for each continent...")
            self._read_lake_db()

            # 4 - prepare output data
            logger.info("> 4 - Preparing output objects for each continent...")
            self._prepare_output_data()

            for cur_continent in self.continent_list :
                logger.info("=============================")
                logger.info("== Processing continent %s ==" % cur_continent)
                logger.info("=============================")

                # 5 - Initialization
                logger.info("> 5 - Initialization of SASLakeSP class")
                my_lake_sp = sas_lake_sp.SASLakeSP(self.obj_pixc_sp[cur_continent], self.obj_pixcvec_sp[cur_continent],
                                                   self.obj_lake_db[cur_continent], self.obj_lake_l[cur_continent], self.obj_lake_r[cur_continent])
                logger.info(self.timer.info(0))

                # 6 - Run pre-processing
                logger.info("> 6 - Run SASpre-processing")
                my_lake_sp.run_preprocessing()
                logger.info(self.timer.info(0))

                # 7 - Run processing
                logger.info("> 7 - Run SASprocessing")
                my_lake_sp.run_processing()
                logger.info(self.timer.info(0))
                logger.info("")
    
                # 8 - Run post-processing
                logger.info("> 8 - Run SASpost-processing")
                my_lake_sp.run_postprocessing()
                logger.info(self.timer.info(0))
                logger.info("")
            
            logger.info("============E N D============")

            # 9 - Write output data
            logger.info("> 9 - Write output data")
            self._write_output_data(self.proc_metadata)
            
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

        for cur_continent in self.continent_list:
            # 1 - Close lake database
            logger.info("Closing lake database...")
            self.obj_lake_db[cur_continent].close_db()

        logger.info("")
        logger.info(self.timer.stop())
        logger.sigmsg("==================================")
        logger.sigmsg("===== lakeSPProcessing = END =====")
        logger.sigmsg("==================================")

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
        out_params["Produce_shp"] = False

        # 1 - Read parameter file
        config = configparser.ConfigParser()
        config.read(self.cmd_file)

        # 2.1 - Retrieve PATHS
        try:
            out_params["param_file"] = os.path.expandvars(config.get("PATHS", "param_file"))
        except:
            out_params["param_file"] = os.path.join(sys.path[0], "lake_sp_param.cfg")
        out_params["laketile_shp_directory"] = config.get("PATHS", "LakeTile shp directory")
        out_params["laketile_edge_directory"] = config.get("PATHS", "LakeTile edge directory")
        out_params["laketile_pixcvec_directory"] = config.get("PATHS", "LakeTile pixcvec directory")
        out_params["output_dir"] = config.get("PATHS", "Output directory")

        # 2.2 - Retrieve cycle and pass number
        out_params["cycle_num"] = config.getint("TILES_INFOS", "Cycle number")
        out_params["pass_num"] = config.getint("TILES_INFOS", "Pass number")

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
                out_params["Produce_shp"] = config.get("OPTIONS", "Produce shp")

        # 5 - Retrieve LOGGING
        log_file = config.get("LOGGING", "logFile")
        out_params["logFile"] = os.path.splitext(log_file)[0] + "_" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + os.path.splitext(log_file)[1]
        out_params["logfilelevel"] = config.get("LOGGING", "logfilelevel")
        out_params["logConsole"] = config.get("LOGGING", "logConsole")
        out_params["logconsolelevel"] = config.get("LOGGING", "logconsolelevel")
        
        # 6 - Retrieve FILE informations
        out_params["PRODUCER"] = config.get("FILE_INFORMATION", "PRODUCER")
        # LAKE_TILE
        out_params["LAKE_TILE_CRID"] = config.get("CRID", "LAKE_TILE_CRID")
        # LAKE_SP
        out_params["LAKE_SP_CRID"] = config.get("CRID", "LAKE_SP_CRID")

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
            if "LogFile" in param_list["logFile"]:
                log_error_file = param_list["logFile"].replace("LogFile", "LogError")
            else:
                log_error_file = os.path.join(param_list["output_dir"], "LogError_" 
                                                  + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + ".log")
            self.cfg.set(section, "errorFile", log_error_file)

            
            # Add PATHS section and parameters
            section = "PATHS"
            self.cfg.add_section(section)
            self.cfg.set(section, "LakeTile shp directory", param_list["laketile_shp_directory"])
            self.cfg.set(section, "LakeTile edge directory", param_list["laketile_edge_directory"])
            self.cfg.set(section, "LakeTile pixcvec directory", param_list["laketile_pixcvec_directory"])
            self.cfg.set(section, "Output directory", param_list["output_dir"])

            # Add TILES INFOS section and parameters
            section = "TILES_INFOS"
            self.cfg.add_section(section)
            self.cfg.set(section, "Cycle number", param_list["cycle_num"])
            self.cfg.set(section, "Pass number", param_list["pass_num"])

            # Add DATABASES section and parameters
            section = "DATABASES"
            self.cfg.add_section(section)
            self.cfg.set(section, "LAKE_DB", param_list["LAKE_DB"])
            self.cfg.set(section, "LAKE_DB_ID", param_list["LAKE_DB_ID"])

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
            # Laketile shp files
            self.cfg.test_var_config_file('PATHS', 'LakeTile shp directory', str)
            my_tools.test_dir(self.cfg.get('PATHS', 'LakeTile shp directory'))
            logger.debug('LakeTile shp directory = ' + str(self.cfg.get('PATHS', 'LakeTile shp directory')))
            # Laketile edge files
            self.cfg.test_var_config_file('PATHS', 'LakeTile edge directory', str)
            my_tools.test_dir(self.cfg.get('PATHS', 'LakeTile edge directory'))
            logger.debug('LakeTile edge directory = ' + str(self.cfg.get('PATHS', 'LakeTile edge directory')))
            # Laketile pixcvec files
            self.cfg.test_var_config_file('PATHS', 'LakeTile pixcvec directory', str)
            my_tools.test_dir(self.cfg.get('PATHS', 'LakeTile pixcvec directory'))
            logger.debug('LakeTile pixcvec directory = ' + str(self.cfg.get('PATHS', 'LakeTile pixcvec directory')))

            # Output directory
            self.cfg.test_var_config_file('PATHS', 'Output directory', str)
            my_tools.test_dir(self.cfg.get('PATHS', 'Output directory'))
            logger.debug('Output directory = ' + str(self.cfg.get('PATHS', 'Output directory')))

            # 1.2 - TILES_INFO section
            # Cycle number
            self.cfg.test_var_config_file('TILES_INFOS', 'Cycle number', int)
            logger.debug('Cycle number = ' + str(self.cfg.get('TILES_INFOS', 'Cycle number')))
            # Pass number
            self.cfg.test_var_config_file('TILES_INFOS', 'Pass number', int)
            logger.debug('Pass number = ' + str(self.cfg.get('TILES_INFOS', 'Pass number')))

            # 1.3 - DATABASES section
            # Lake database full path
            if self.cfg.get('DATABASES', 'LAKE_DB') is None:
                logger.debug('LAKE_DB not filled => LakeTile product not linked to a lake database and continent')
            elif self.cfg.get('DATABASES', 'LAKE_DB').endswith(".shp"):
                self.cfg.test_var_config_file('DATABASES', 'LAKE_DB', str)
                my_tools.test_file(self.cfg.get('DATABASES', 'LAKE_DB'))
                logger.debug('LAKE_DB = ' + str(self.cfg.get('DATABASES', 'LAKE_DB')))
                # Lake identifier attribute name in the database
                if self.cfg.get('DATABASES', 'LAKE_DB_ID'):
                    self.cfg.test_var_config_file('DATABASES', 'LAKE_DB_ID', str)
                    logger.debug('LAKE_DB_ID = ' + str(self.cfg.get('DATABASES', 'LAKE_DB_ID')))
                else :
                    logger.warning('LAKE_DB file given but the lake_id fieldname is missing')
            elif self.cfg.get('DATABASES', 'LAKE_DB').endswith(".sqlite"):
                self.cfg.test_var_config_file('DATABASES', 'LAKE_DB', str)
                my_tools.test_file(self.cfg.get('DATABASES', 'LAKE_DB'))
                logger.debug('LAKE_DB = ' + str(self.cfg.get('DATABASES', 'LAKE_DB')))
            else :
                logger.debug('Unknown LAKE_DB file format for file : %s' % (self.cfg.get('DATABASES', 'LAKE_DB')))

            # 1.4 - OPTIONS section
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
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'MIN_SIZE', float, val_default=0.01, logger=logger)
            logger.debug('MIN_SIZE = ' + str(self.cfg.get('CONFIG_PARAMS', 'MIN_SIZE')))
            # Maximal standard deviation of height inside a lake
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'STD_HEIGHT_MAX', float, val_default=10)
            logger.debug('STD_HEIGHT_MAX = ' + str(self.cfg.get('CONFIG_PARAMS', 'STD_HEIGHT_MAX')))

            # To improve PixC golocation (=True) or not (=False)
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'IMP_GEOLOC', bool, val_default=True, logger=logger)
            logger.debug('IMP_GEOLOC = ' + str(self.cfg.get('CONFIG_PARAMS', 'IMP_GEOLOC')))
            # Method to compute lake boundary or polygon hull
            # 0=convex hull 1=concav hull (1.0=with alpha param (default) 1.1=without) 2=concav hull radar vectorisation 3=alpha shape with CGAL
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'HULL_METHOD', float, valeurs=[0, 1, 1.1, 2, 3], val_default=2.0, logger=logger)

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
            self.cfg.test_var_config_file('CONFIG_PARAMS', 'BIGLAKE_MIN_SIZE', float, val_default=50, logger=logger)
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
    PARSER = argparse.ArgumentParser(description="Compute SWOT SP product\
            from all tiles of LAKE TILE. ")
    PARSER.add_argument("command_file", help="command file (*.cfg)")
    ARGS = PARSER.parse_args()
    PGE = None
    try:
        # 1 - Instantiate PGE
        PGE = PGELakeSP(ARGS.command_file)

        # 2 - Start PGE Lake Tile
        PGE.start()
    finally:
        if PGE is not None:
            # 3 - Stop PGE Lake Tile
            PGE.stop()
