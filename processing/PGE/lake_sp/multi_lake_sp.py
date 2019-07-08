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
import cnes.common.lib_lake.locnes_variables as my_var
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

        # BDLac and continent file
        self.lake_db = in_params["LAKE_DB"]
        self.influence_lake_db = in_params["INFLUENCE_LAKE_DB"]
        self.lake_db_id = in_params["LAKE_DB_ID"]
        self.continent_file = in_params["CONTINENT_FILE"]

        # Tiles info
        self.cycle_num = in_params["cycle_num"]  # Cycle number
        self.pass_num = in_params["pass_num"]  # Pass number

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
        print("[multiLakeSPProcessing]   INPUT DIR for LakeTile shp files = %s" % self.laketile_shp_dir)
        my_tools.testDir(self.laketile_shp_dir)
        # 1.2 - LakeTile egde directory
        print("[multiLakeSPProcessing]   INPUT DIR for LakeTile edge files = %s" % self.laketile_edge_dir)
        my_tools.testDir(self.laketile_edge_dir)
        # 1.3 - LakeTile pixcvec directory
        print("[multiLakeSPProcessing]   INPUT DIR for LakeTile pixcvec files = %s" % self.laketile_pixcvec_dir)
        my_tools.testDir(self.laketile_pixcvec_dir)
        # 1.4 - Output directory
        print("[multiLakeSPProcessing]   OUTPUT DIR = %s" % self.output_dir)
        my_tools.testDir(self.output_dir)
        print("")

        # 2 - Get input files
        print("[multiLakeSPProcessing] > 2 - Retrieving input files ...")

        # 2.1 - List all files in self.laketile_shp_dir
        tmp_list = os.listdir(self.laketile_shp_dir)

        # 2.2 - Compute file prefix regarding cycle / pass / tile conditions
        cond_prefix = my_var.LAKE_TILE_PREFIX  # Deal with all laketile files in laketile directories
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
        self.laketile_files_list = []
        self.cycle_pass_list = []
        self.nb_tiles_list= []
        for cur_item in tmp_list:

            # Test if file meets the condition of laketile shp file
            if cur_item.startswith(cond_prefix) and cur_item.endswith(my_var.LAKE_TILE_SHP_SUFFIX):
                # find edge and pixcvec associated files
                cur_laketile_shp_file = os.path.join(self.laketile_shp_dir, cur_item)
                cur_laketile_edge_file = os.path.join(self.laketile_edge_dir, cur_item.replace(my_var.LAKE_TILE_SHP_SUFFIX, my_var.LAKE_TILE_EDGE_SUFFIX))
                cur_laketile_pixcvec_file = os.path.join(self.laketile_pixcvec_dir, cur_item.replace(my_var.LAKE_TILE_SHP_SUFFIX, my_var.LAKE_TILE_PIXCVEC_SUFFIX))

                # If associated edge and pixcvec associated files exist, add tuple of filenames
                if os.path.exists(cur_laketile_edge_file) and os.path.exists(cur_laketile_edge_file):
                    information = locnes_filenames.getInfoFromFilename(cur_item, "LakeTile")
                    cycle_num = information["cycle"]
                    pass_num = information["pass"]

                    if (cycle_num, pass_num) not in self.cycle_pass_list : # if new pair of cycle and pass, create new list of tiles
                        self.cycle_pass_list.append((cycle_num, pass_num))
                        self.laketile_files_list.append(([cur_laketile_shp_file], [cur_laketile_edge_file], [cur_laketile_pixcvec_file]))
                        self.nb_tiles_list.append(0)
                    else : # if pair of (cycle and path) exists, add tile to the list
                        cycle_pass_idx = self.cycle_pass_list.index((cycle_num, pass_num))
                        self.laketile_files_list[cycle_pass_idx][0].append(cur_laketile_shp_file)
                        self.laketile_files_list[cycle_pass_idx][1].append(cur_laketile_edge_file)
                        self.laketile_files_list[cycle_pass_idx][2].append(cur_laketile_pixcvec_file)
                        self.nb_tiles_list[cycle_pass_idx] += 1

        print("[multiLakeSPProcessing]   --> %d cycle and pass to deal with" % len(self.cycle_pass_list))
        for i, (cycle_num, pass_num) in enumerate(self.cycle_pass_list):
            print("[multiLakeSPProcessing]  cycle %s pass %s contains %d tiles" % (cycle_num, pass_num, self.nb_tiles_list[i]))
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

        for i, (cycle_num, pass_num) in enumerate(self.cycle_pass_list):  # Deal with all selected files

            print("****************************************************************************************************************")
            print("***** Dealing with cycle %s and pass %s *****" % (cycle_num, pass_num))
            print("****************************************************************************************************************")
            print("")
            print("")

            # 1 - Initialization
            cmd_file = self.create_cmd_file(i)

            my_lake_tile = pge_lake_sp.PGELakeSP(cmd_file)

            # 2 - Run
            my_lake_tile.start()

            # 3 - Stop
            my_lake_tile.stop()

            print("")
            print(timer.stop())
            print("")
            print("")

        print("****************************************************************************************************************")
        print("")
        print("")

    def create_cmd_file(self, indf):
        """
        Create command file for PGE_L2_HR_LakeSP for each input cycle and pass
        
        :param indf: index of current input cycle and pass
        :type indf: int
        
        :return: out_cmd_file = command file full path
        :rtype: string
        """

        # 1 - Get LakeTile files full path
        (laketile_shp_files, laketile_edge_files, laketile_pixcvec_files) = self.laketile_files_list[indf]

        # 2 - Retrieve SP infos
        (cycle_num, pass_num) = self.cycle_pass_list[indf]

        # 3 - Init log filename
        if self.log_file is None:
            log_file = os.path.join(self.output_dir, "LogFile_" + str(cycle_num) + "_" +
                                    str(pass_num) + ".log")
            print("Log file : " + log_file)
        else:
            log_file = os.path.splitext(self.log_file)[0] + "_" + str(cycle_num) + "_" + str(pass_num) + os.path.splitext(self.log_file)[1]

        # 4 - Init command filename
        cmd_filename = "lake_sp_command_" + str(cycle_num) + "_" + str(pass_num) + ".cfg"
        out_cmd_file = os.path.join(self.output_dir, cmd_filename)
        
        # 5 - Write command variables in command file
        writer_command_file = open(out_cmd_file, "w")  # Open file in writing mode
        
        # 5.1 - Fill PATHS section
        writer_command_file.write("[PATHS]\n")
        if self.param_file is not None:
            writer_command_file.write("param_file = %s\n" % self.param_file)
        else:
            # needed by jenkins script
            writer_command_file.write("param_file = " + os.path.join(sys.path[0], "lake_sp_param.cfg") + "\n")

        writer_command_file.write("LakeTile shp files = " + ";".join(laketile_shp_files) + "\n")
        writer_command_file.write("LakeTile edge files = " + ";".join(laketile_edge_files) + "\n")
        writer_command_file.write("LakeTile pixcvec files = " + ";".join(laketile_pixcvec_files) + "\n")
        writer_command_file.write("Output directory = " + self.output_dir + "\n")
        writer_command_file.write("\n")

        writer_command_file.write("[TILES_INFOS]\n")
        writer_command_file.write("Cycle number = " + cycle_num +"\n")
        writer_command_file.write("Pass number = " + pass_num +"\n")
        writer_command_file.write("\n")


        # 5.2 - Fill DATABASES section
        writer_command_file.write("[DATABASES]\n")
        writer_command_file.write("# Lake a priori database\n")
        if self.lake_db is not None:
            writer_command_file.write("LAKE_DB = " + self.lake_db + "\n")
        writer_command_file.write("# Lake identifier attribute name in the database\n")
        if self.influence_lake_db is not None:
            writer_command_file.write("INFLUENCE_LAKE_DB = " + self.influence_lake_db + "\n")
        if self.lake_db_id is not None:
            writer_command_file.write("LAKE_DB_ID = " + self.lake_db_id + "\n")
        writer_command_file.write("# Shapefile with polygons of continents\n")
        if self.continent_file is not None:
            writer_command_file.write("CONTINENT_FILE = " + self.continent_file + "\n")
        else:
            # needed by script for Jenkins
            writer_command_file.write("CONTINENT_FILE = /work/ALT/swot/swotpub/BD/major_basins/FAO/major_hydrobasins.shp\n")
        writer_command_file.write("\n")
        
        # 5.3 - Fill OPTIONS section
        writer_command_file.write("[OPTIONS]\n")
        writer_command_file.write("# To also produce LakeTile_edge and LakeTile_pixcvec as shapefiles (=True); else=False (default)\n")
        writer_command_file.write("Produce shp = " + str(self.flag_prod_shp) + "\n")
        writer_command_file.write("\n")
        
        # 5.4 - Fill LOGGING section
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
        
        # 5.5 - Fill CRID section
        writer_command_file.write("# CRID information\n")
        writer_command_file.write("[CRID]\n")
        writer_command_file.write("# Composite Release IDentifier for LakeTile processing\n")        
        writer_command_file.write("LAKE_TILE_CRID = Dx0000 \n")
        writer_command_file.write("# Composite Release IDentifier for LakeSP processing\n")
        writer_command_file.write("LAKE_SP_CRID = Dx0000 # Product generator\n")
        writer_command_file.write("\n")

        # 5.6 - Fill FILE_INFORMATION section
        writer_command_file.write("# File informations\n")
        writer_command_file.write("[FILE_INFORMATION]\n")
        writer_command_file.write("# Product generator\n")
        writer_command_file.write("PRODUCER = CNES\n")
        writer_command_file.write("\n")
        writer_command_file.write("PIXC_PREFIX = SWOT_L2_HR_PIXC_\n")
        writer_command_file.write('PIXC_PATTERN_PRINT = PIXC_PREFIX + "<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>.nc"\n')
        writer_command_file.write('# Indices when PIXC_PATTERN.split("_"); None if value not in filename\n')
        writer_command_file.write("\n")
        writer_command_file.write('PIXC_PATTERN_IND = {"cycle": 4, "pass": 5, "tile_ref": 6, "start_date": 7, "stop_date": 8, "crid": 9, "counter": 10} \n')
        writer_command_file.write('PIXCVEC_RIVER_PREFIX = SWOT_L2_HR_PIXCVecRiver_\n')
        writer_command_file.write('PIXCVEC_RIVER_PATTERN_PRINT = PIXCVEC_RIVER_PREFIX + "<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>.nc"\n')
        writer_command_file.write('# Indices when PIXCVEC_RIVER_PATTERN.split("_"); None if value not in filename\n')
        writer_command_file.write('PIXCVEC_RIVER_PATTERN_IND = {"cycle": 4, "pass": 5, "tile_ref": 6, "start_date": 7, "stop_date": 8, "crid": 9, "counter": 10}\n')
        writer_command_file.write('LAKE_TILE_PREFIX = SWOT_L2_HR_LakeTile_\n')
        writer_command_file.write("# LakeTile filename with %03d=cycle number %03d=pass number %s=tile ref %s=swath %s=begin date %s=end date %s=CRID %s=counter %s=suffix \n")
        writer_command_file.write("\n")
        writer_command_file.write('LAKE_TILE_PATTERN = LAKE_TILE_PREFIX + "%03d_%03d_%s_%s_%s_%s_%02d%s"\n')
        writer_command_file.write('LAKE_TILE_PATTERN_PRINT = LAKE_TILE_PREFIX + "%s<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>"\n')
        writer_command_file.write('# Indices when LAKE_TILE_*_PATTERN.split("_"); None if value not in filename \n')
        writer_command_file.write("\n")
        writer_command_file.write('LAKE_TILE_PATTERN_IND = {"cycle": 4, "pass": 5, "tile_ref": 6, "start_date": 7, "stop_date": 8, "crid": 9, "counter": 10}\n')
        writer_command_file.write('LAKE_TILE_SHP_SUFFIX = .shp\n')
        writer_command_file.write('LAKE_TILE_SHP_META_SUFFIX = .shp.xml\n')
        writer_command_file.write('LAKE_TILE_EDGE_SUFFIX = _edge.nc\n')
        writer_command_file.write('LAKE_TILE_PIXCVEC_SUFFIX = _pixcvec.nc\n')
        writer_command_file.write("\n")
        writer_command_file.write("PIXCVEC_PREFIX = SWOT_L2_HR_PIXCVec_\n")
        writer_command_file.write("PIXCVEC_SUFFIX = .nc\n")
        writer_command_file.write("\n")
        writer_command_file.write("# PIXCVec filename with %03d=cycle number %03d=pass number %s=tile ref %s=begin date %s=end date %s=CRID %s=counter  \n")
        writer_command_file.write('PIXCVEC_PATTERN = PIXCVEC_PREFIX + "%03d_%03d_%s_%s_%s_%s_%02d" + PIXCVEC_SUFFIX \n')
        writer_command_file.write("\n")
        writer_command_file.write("LAKE_SP_PREFIX = SWOT_L2_HR_LakeSP_\n")
        writer_command_file.write("# LakeSP filename with %03d=cycle number %03d=pass number %s=continent %s=begin date %s=end date %s=CRID %s=counter \n")
        writer_command_file.write('LAKE_SP_PATTERN = LAKE_SP_PREFIX + "%03d_%03d_%s_%s_%s_%s_%02d.shp"\n')
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
    out_params["INFLUENCE_LAKE_DB"] = None
    out_params["LAKE_DB_ID"] = None
    out_params["CONTINENT_FILE"] = None
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
        # Lake a priori database
        if "lake_db" in list_db:
            out_params["LAKE_DB"] = config.get("DATABASES", "LAKE_DB")
        if "influence_lake_db" in list_db:
            out_params["INFLUENCE_LAKE_DB"] = config.get("DATABASES", "INFLUENCE_LAKE_DB")
        if "lake_db_id" in list_db:
            out_params["LAKE_DB_ID"] = config.get("DATABASES", "LAKE_DB_ID")
        # Continent file
        if "continent_file" in list_db:
            out_params["CONTINENT_FILE"] = config.get("DATABASES", "CONTINENT_FILE")

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
    my_tools.testFile(args.command_file, IN_extent=".cfg")  # Test existance and extension
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
