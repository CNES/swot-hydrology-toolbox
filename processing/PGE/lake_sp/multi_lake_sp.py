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
.. module:: multi_lake_sp.py
    :synopsis Process PGE_L2_HR_LakeSP, i.e. generate L2_HR_LakeSP
    product from all tile of L2_HR_LAKE_TILE product
    Created on 08/23/2018

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR
                  Cécile Cazals - C-S
This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
import configparser as cfg
import datetime
import os
import sys

import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_timer as my_timer
import cnes.common.lib_lake.locnes_filenames as locnes_filenames
import pge_lake_sp


class MultiLakeSP(object):
    """
    Class MultiLakeTile
    Main class to run many SAS LakeTile processes
    """
    
    def __init__(self, in_params):
        """
        Constructor: initialize variables
        
        :param in_params: input parameters to run the processor
        :type in_params: dictionary
        """
        print("[multiLakeSPProcessing] == INIT ==")

        # Directories
        self.param_file = in_params["param_file"]  # Parameter file
        self.laketile_shp_dir = in_params["laketile_shp_dir"]  # PIXC files directory
        self.laketile_edge_dir = in_params["laketile_edge_dir"]  # PIXCVecRiver files directory
        self.laketile_pixcvec_dir = in_params["laketile_pixcvec_dir"]  # PIXCVecRiver files directory
        self.output_dir = in_params["output_dir"]  # Output directory

        # BDLac and basins file
        self.lake_db = in_params["LAKE_DB"]
        self.lake_db_id = in_params["LAKE_DB_ID"]

        # Tiles info
        self.cycle_num = in_params["cycle_num"]  # Cycle number
        self.pass_num = in_params["pass_num"]  # Pass number
        self.cycle_pass_set = set()

        # Flag to produce LakeTile_edge and LakeTile_pixcvec shapefiles
        self.flag_prod_shp = in_params["flag_prod_shp"]
        
        # Log level
        self.log_file = in_params["logFile"]
        self.log_file_level = in_params["logfilelevel"]
        self.log_console = in_params["logConsole"]
        self.log_console_level = in_params["logconsolelevel"]

        # List of input files
        self.list_laketile_shp_files = []
        self.list_laketile_edge_files = []
        self.list_laketile_pixcvec_files = []
        self.nb_input = 0  # Nb of input files

    def run_preprocessing(self):
        """
        Retrieve the list of input files, i.e. L2_HR_PIXC tile file main and associated L2_HR_PIXCVec
        """

        print("")
        print("")
        print("[multiLakeSPProcessing] PRE-PROCESSING...")
        print("")

        # 1 - Test existence of directories
        print("[multiLakeSPProcessing] > 1 - Testing existence of working directories...")
        # 1.1 - LakeTile shp directory
        print("[multiLakeSPProcessing]   INPUT DIR for LakeTile shp directory = %s" % self.laketile_shp_dir)
        my_tools.test_dir(self.laketile_shp_dir)
        # 1.2 - LakeTile egde directory
        print("[multiLakeSPProcessing]   INPUT DIR for LakeTile edge directory = %s" % self.laketile_edge_dir)
        my_tools.test_dir(self.laketile_edge_dir)
        # 1.3 - LakeTile pixcvec directory
        print("[multiLakeSPProcessing]   INPUT DIR for LakeTile pixcvec directory = %s" % self.laketile_pixcvec_dir)
        my_tools.test_dir(self.laketile_pixcvec_dir)
        # 1.4 - Output directory
        print("[multiLakeSPProcessing]   OUTPUT DIR = %s" % self.output_dir)
        my_tools.test_dir(self.output_dir)
        print("")

        # 2 - Get input files
        print("[multiLakeSPProcessing] > 2 - Retrieving input files ...")

        # 2.1 - List all files in self.laketile_shp_dir
        tmp_list = os.listdir(self.laketile_shp_dir)

        # 2.2 - Compute file prefix regarding cycle / pass / tile conditions
        lake_tile_prefix = "SWOT_L2_HR_LakeTile_" # from my_var.LAKE_TILE_PREFIX
        lake_tile_shp_suffix = ".shp" # from my_var.LAKE_TILE_SHP_SUFFIX
        lake_tile_edge_suffix = "_edge.nc" # from my_var.LAKE_TILE_EDGE_SUFFIX
        lake_tile_pixcvec_suffix = "_pixcvec.nc" # from my_var.LAKE_TILE_PIXCVEC_SUFFIX
        cond_prefix = lake_tile_prefix  # Deal with all laketile files in laketile directories
        if (self.cycle_num is None) or (self.cycle_num == "-1"):
            print("[multiLakeSPProcessing]   All LakeTile files in the input directories")

        else:  # Deal with LakeTile files with cycle = self.cycle_num
            cond_prefix += "%03d" % self.cycle_num
            if (self.pass_num is None) or (self.pass_num == "-1"):
                print("[multiLakeSPProcessing]   Lake tile files with cycle=%03d" % self.cycle_num)
            else:  # Deal with LakeTile files with cycle = self.cycle_num and pass = self.pass_num
                print("[multiLakeSPProcessing]   Lake tile files with cycle=%03d and pass=%03d" % (self.cycle_num, self.pass_num))
                cond_prefix += "_%03d" % self.pass_num

        # 2.3 - For each listed laketile shp file, get related laketile edge and pixcvec files if they exist
        for cur_item in tmp_list:

            # Test if file meets the condition of laketile shp file
            if cur_item.startswith(cond_prefix) and cur_item.endswith(lake_tile_shp_suffix):

                information = locnes_filenames.get_info_from_filename(cur_item, "LakeTile")
                cycle_num = information["cycle"]
                pass_num = information["pass"]

                self.cycle_pass_set.add((cycle_num, pass_num))

        print("[multiLakeSPProcessing]   --> %d cycle and pass to deal with" % len(self.cycle_pass_set))
        for (cycle_num, pass_num) in self.cycle_pass_set:
            print("[multiLakeSPProcessing]  cycle %s pass %s " % (cycle_num, pass_num))
        print("")

    def run_processing(self):
        """
        Process SAS_L2_HR_LakeTile for each input tile
        """
        print("")
        print("")
        print("[multiLakeSPProcessing] PROCESSING...")
        print("")
        print("")

        timer_proc = my_timer.Timer()
        timer_proc.start()

        for (cycle_num, pass_num) in self.cycle_pass_set:  # Deal with all selected files

            print("****************************************************************************************************************")
            print("***** Dealing with cycle %s and pass %s *****" % (cycle_num, pass_num))
            print("****************************************************************************************************************")
            print("")
            print("")

            # 1 - Initialization
            cmd_file = self.create_cmd_file(cycle_num, pass_num)

            my_lake_tile = pge_lake_sp.PGELakeSP(cmd_file)

            # 2 - Run
            my_lake_tile.start()

            # 3 - Stop
            my_lake_tile.stop()

            print("")
            print("")
            print("")

        print(timer_proc.stop())
        print("****************************************************************************************************************")
        print("")
        print("")

    def create_cmd_file(self, cycle_num, pass_num):
        """
        Create command file for PGE_L2_HR_LakeSP for each input cycle and pass
        
        :param cycle_num: cycle number
        :type cycle_num: int
        :param pass_num: pass number
        :type pass_num: int
        
        :return: out_cmd_file = command file full path
        :rtype: string
        """

        # 1 - Init log filename
        if self.log_file is None:
            log_file = os.path.join(self.output_dir, "LogFile_" + str(cycle_num) + "_" +
                                    str(pass_num) + ".log")
            print("Log file : " + log_file)
        else:
            log_file = os.path.splitext(self.log_file)[0] + "_" + str(cycle_num) + "_" + str(pass_num) + os.path.splitext(self.log_file)[1]

        # 2 - Init command filename
        cmd_filename = "lake_sp_command_" + str(cycle_num) + "_" + str(pass_num) + ".cfg"
        out_cmd_file = os.path.join(self.output_dir, cmd_filename)
        
        # 3 - Write command variables in command file
        writer_command_file = open(out_cmd_file, "w")  # Open file in writing mode
        
        # 3.1 - Fill PATHS section
        writer_command_file.write("[PATHS]\n")
        if self.param_file is not None:
            writer_command_file.write("param_file = %s\n" % self.param_file)
        else:
            # needed by jenkins script
            writer_command_file.write("param_file = " + os.path.join(sys.path[0], "lake_sp_param.cfg") + "\n")

        writer_command_file.write("LakeTile shp directory = " + self.laketile_shp_dir + "\n")
        writer_command_file.write("LakeTile edge directory = " + self.laketile_edge_dir + "\n")
        writer_command_file.write("LakeTile pixcvec directory = " + self.laketile_pixcvec_dir + "\n")
        writer_command_file.write("Output directory = " + self.output_dir + "\n")
        writer_command_file.write("\n")

        writer_command_file.write("[TILES_INFOS]\n")
        writer_command_file.write("Cycle number = " + str(cycle_num) +"\n")
        writer_command_file.write("Pass number = " + str(pass_num) +"\n")
        writer_command_file.write("\n")

        # 3.2 - Fill DATABASES section
        writer_command_file.write("[DATABASES]\n")
        writer_command_file.write("# Prior lake database\n")
        if self.lake_db is not None:
            writer_command_file.write("LAKE_DB = " + self.lake_db + "\n")
        writer_command_file.write("# Lake identifier attribute name in the database\n")
        if self.lake_db_id is not None:
            writer_command_file.write("LAKE_DB_ID = " + self.lake_db_id + "\n")
        writer_command_file.write("\n")
        
        # 3.3 - Fill OPTIONS section
        writer_command_file.write("[OPTIONS]\n")
        writer_command_file.write("# To also produce LakeTile_edge and LakeTile_pixcvec as shapefiles (=True); else=False (default)\n")
        writer_command_file.write("Produce shp = " + str(self.flag_prod_shp) + "\n")
        writer_command_file.write("\n")
        
        # 3.4 - Fill LOGGING section
        writer_command_file.write("[LOGGING]\n")
        writer_command_file.write("# Log filename\n")
        writer_command_file.write("logFile = " + log_file + "\n")
        writer_command_file.write("# Log level put inside the file\n")
        writer_command_file.write("logfilelevel = " + self.log_file_level + "\n")
        writer_command_file.write("# Is log console output ?\n")
        writer_command_file.write("logConsole = " + str(self.log_console) + "\n")
        writer_command_file.write("# Log level print in console\n")
        writer_command_file.write("logconsolelevel = " + self.log_console_level + "\n")
        writer_command_file.write("\n")
        
        # 3.5 - Fill CRID section
        writer_command_file.write("# CRID information\n")
        writer_command_file.write("[CRID]\n")
        writer_command_file.write("# Composite Release IDentifier for LakeTile processing\n")        
        writer_command_file.write("LAKE_TILE_CRID = Dx0000 \n")
        writer_command_file.write("# Composite Release IDentifier for LakeSP processing\n")
        writer_command_file.write("LAKE_SP_CRID = Dx0000\n")
        writer_command_file.write("\n")

        # 3.6 - Fill FILE_INFORMATION section
        writer_command_file.write("# File informations\n")
        writer_command_file.write("[FILE_INFORMATION]\n")
        writer_command_file.write("# Product generator\n")
        writer_command_file.write("PRODUCER = CNES\n")
        writer_command_file.write("#######################################\n")        
        writer_command_file.close()  # Close command file
        return out_cmd_file

#######################################
        

def read_command_file(in_filename):
    """
    Read the command file in input and store parameters in a dictionary
    
    :param in_filename: parameter file full path
    :type in_filename: string
    
    :return: out_params = dictionary containing parameters
    :rtype: dict
    """
    print("[multiLakeSPProcessing] == read_command_file = %s ==" % in_filename)

    # 0 - Init output dictionary
    out_params = {}
    # Default values
    out_params["param_file"] = None
    out_params["LAKE_DB"] = None
    out_params["LAKE_DB_ID"] = None
    out_params["cycle_num"] = None
    out_params["pass_num"] = None
    out_params["flag_prod_shp"] = False
    out_params["logFile"] = None
    out_params["logfilelevel"] = "DEBUG"
    out_params["logConsole"] = True
    out_params["logconsolelevel"] = "DEBUG"

    # 1 - Read parameter file
    config = cfg.ConfigParser()
    config.read(in_filename)

    # 2 - Retrieve PATHS
    list_paths = config.options("PATHS")
    if "param_file" in list_paths:
        out_params["param_file"] = config.get("PATHS", "param_file")
    out_params["laketile_shp_dir"] = config.get("PATHS", "LakeTile shp directory")
    out_params["laketile_edge_dir"] = config.get("PATHS", "LakeTile edge directory")
    out_params["laketile_pixcvec_dir"] = config.get("PATHS", "LakeTile pixcvec directory")
    out_params["output_dir"] = config.get("PATHS", "Output directory")

    # 3 - Retrieve DATABASES
    if "DATABASES" in config.sections():
        list_db = config.options("DATABASES")
        # Prior lake database
        if "lake_db" in list_db:
            out_params["LAKE_DB"] = config.get("DATABASES", "LAKE_DB")
        if "lake_db_id" in list_db:
            out_params["LAKE_DB_ID"] = config.get("DATABASES", "LAKE_DB_ID")

    # 4 - Retrieve TILES_INFOS
    if "TILES_INFOS" in config.sections():
        list_options = config.options("TILES_INFOS")
        # Cycle number
        if "cycle number" in list_options:
            out_params["cycle_num"] = config.getint("TILES_INFOS", "Cycle number")
        # Pass number
        if "pass number" in list_options:
            out_params["pass_num"] = config.getint("TILES_INFOS", "Pass number")

    # 5 - Retrieve OPTIONS
    if "OPTIONS" in config.sections():
        list_options = config.options("OPTIONS")
        # Flag to also produce PIXCVEC nc files as shapefiles (=True); else=False (default)
        if "produce shp" in list_options:
            out_params["flag_prod_shp"] = config.getboolean("OPTIONS", "Produce shp")

    # 5 - Retrieve optionnal CONFIG_OVERWRITE
    # Still useful?
    #my_var.overwriteConfig_from_cfg(config)
    # 6 - Retrieve LOGGING
    if "LOGGING" in config.sections():
        out_params["logFile"] = config.get("LOGGING", "logFile")
        if out_params["logFile"] is not None:
            out_params["logFile"] = out_params["logFile"].replace("<date>", datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
        out_params["logfilelevel"] = config.get("LOGGING", "logfilelevel")
        out_params["logConsole"] = config.get("LOGGING", "logConsole")
        out_params["logconsolelevel"] = config.get("LOGGING", "logconsolelevel")

    return out_params


#######################################


if __name__ == '__main__':

    # 0 - Parse inline parameters
    parser = argparse.ArgumentParser(description="Compute multiple SWOT LakeSP products from multiple tiles of LAKE TILE product.")
    parser.add_argument("command_file", help="command file (*.cfg)")
    args = parser.parse_args()

    print("===== multiLakeSPProcessing = BEGIN =====")
    print("")
    timer = my_timer.Timer()
    timer.start()

    # 1 - Read command file
    print("WORKING VARIABLES")
    print()
    my_tools.test_file(args.command_file, in_extent=".cfg")  # Test existance and extension
    my_params = read_command_file(args.command_file)  # Read variables in command file

    # 2 - Initialization
    multi_lake_sp = MultiLakeSP(my_params)
    print(timer.info(0))

    # 3 - Run pre-processing
    multi_lake_sp.run_preprocessing()
    print(timer.info(0))

    # 4 - Run processing
    multi_lake_sp.run_processing()
    print(timer.info(0))

    print("")
    print(timer.stop())
    print("===== multiLakeTileProcessing = END =====")
