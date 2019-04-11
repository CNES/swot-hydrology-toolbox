#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
.. module:: multi_lake_tile.py
    :synopsis: Process PGE_L2_HR_LakeTile (i.e. generate L2_HR_LakeTile product from one tile of L2_HR_PIXC product and associated L2_HR_PIXCVec product) for multiple tiles
    Created on 08/23/2018

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

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
import pge_lake_tile


class MultiLakeTile(object):
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
        print("[multiLakeTileProcessing] == INIT ==")

        # Directories
        self.param_file = in_params["param_file"]  # Parameter file
        self.pixc_dir = in_params["pixc_dir"]  # PIXC files directory
        self.pixc_vec_river_dir = in_params["pixc_vec_river_dir"]  # PIXCVecRiver files directory
        self.output_dir = in_params["output_dir"]  # Output directory

        # BDLac and continent file
        self.lake_db = in_params["LAKE_DB"]
        self.lake_db_id = in_params["LAKE_DB_ID"]
        self.continent_file = in_params["CONTINENT_FILE"]

        # Tiles info
        self.cycle_num = in_params["cycle_num"]  # Cycle number
        self.pass_num = in_params["pass_num"]  # Pass number
        self.tile_ref = in_params["tile_ref"]  # Tile id

        # Flag to produce LakeTile_edge and LakeTile_pixcvec shapefiles
        self.flag_prod_shp = in_params["flag_prod_shp"]
        
        # Log level
        self.log_file = in_params["logFile"]
        self.log_file_level = in_params["logfilelevel"]
        self.log_console = in_params["logConsole"]
        self.log_console_level = in_params["logconsolelevel"]

        # List of input files
        self.list_pixc = []  # PIXC files
        self.list_pixc_vec_river = []  # PIXCVecRiver
        self.nb_input = 0  # Nb of input files

    def run_preprocessing(self):
        """
        Retrieve the list of input files, i.e. L2_HR_PIXC tile file main and associated L2_HR_PIXCVec
        """

        print("")
        print("")
        print("[multiLakeTileProcessing] PRE-PROCESSING...")
        print("")

        # 1 - Test existence of directories
        print("[multiLakeTileProcessing] > 1 - Testing existence of working directories...")
        # 1.1 - PIXC directory
        print("[multiLakeTileProcessing]   INPUT DIR for PIXC files = %s" % self.pixc_dir)
        my_tools.testDir(self.pixc_dir)
        # 1.3 - PIXCVecRiver directory
        print("[multiLakeTileProcessing]   INPUT DIR for PIXCVecRiver files = %s" % self.pixc_vec_river_dir)
        my_tools.testDir(self.pixc_vec_river_dir)
        # 1.4 - Output directory
        print("[multiLakeTileProcessing]   OUTPUT DIR = %s" % self.output_dir)
        my_tools.testDir(self.output_dir)
        print("")

        # 2 - Get input files
        print("[multiLakeTileProcessing] > 2 - Retrieving input files...")

        # 2.1 - Compute file prefix regarding cycle / pass / tile conditions
        cond_prefix = my_var.PIXC_PREFIX  # Deal with all PixC files in self.pixc_dir
        if (self.cycle_num is None) or (self.cycle_num == "-1"):
            print("[multiLakeTileProcessing]   All PixC files in the input directory")
        else:  # Deal with PixC files with cycle = self.cycle_num
            cond_prefix += "%03d" % self.cycle_num
            if (self.pass_num is None) or (self.pass_num == "-1"):
                print("[multiLakeTileProcessing]   PixC files with cycle=%03d" % self.cycle_num)
            else:  # Deal with PixC files with cycle = self.cycle_num and pass = self.pass_num
                cond_prefix += "_%03d" % self.pass_num
                if self.tile_ref is not None:  # Deal with PixC files with cycle = self.cycle_num, pass = self.pass_num and tile id = self.tile_ref
                    print("[multiLakeTileProcessing]   PixC files with cycle=%03d , pass=%03d , tile=%s" % (self.cycle_num, self.pass_num, self.tile_ref))
                    cond_prefix += "_%s" % self.tile_ref
                else:
                    print("[multiLakeTileProcessing]   PixC files with cycle=%03d and pass=%03d" % (self.cycle_num, self.pass_num))

        # 2.2 - List all files in self.pixc_dir
        tmp_list = os.listdir(self.pixc_dir)

        # 2.3 - For each listed file, get related PIXCVecRiver files if they exist
        cur_pixc_vec_river = None
        for cur_item in tmp_list:

            # Test if file meets the condition
            if cur_item.startswith(cond_prefix):  # Test if it's a wanted PIXC file

                # Associated PIXCVecRiver file name
                cur_pixc_vec_river = cur_item.replace(my_var.PIXC_PREFIX, my_var.PIXCVEC_RIVER_PREFIX)

                # If associated PIXCVecRiver file exists, add pair of filenames
                if os.path.exists(os.path.join(self.pixc_vec_river_dir, cur_pixc_vec_river)):
                    self.list_pixc.append(cur_item)
                    self.list_pixc_vec_river.append(cur_pixc_vec_river)
                    self.nb_input += 1

        print("[multiLakeTileProcessing]   --> %d tile(s) to deal with" % self.nb_input)
        print("")

    def run_processing(self):
        """
        Process SAS_L2_HR_LakeTile for each input tile
        """
        print("")
        print("")
        print("[multiLakeTileProcessing] PROCESSING...")
        print("")
        print("")

        timer_proc = my_timer.Timer()
        timer_proc.start()

        if self.nb_input != 0:

            for indf in range(self.nb_input):  # Deal with all selected files

                print("****************************************************************************************************************")
                print("***** Dealing with tile %d / %d = %s *****" % (indf+1, self.nb_input, self.list_pixc[indf]))
                print("****************************************************************************************************************")
                print("")
                print("")

                # 1 - Initialization
                cmd_file = self.create_cmd_file(indf)
                print(cmd_file)
                my_lake_tile = pge_lake_tile.PGELakeTile(cmd_file)

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
        Create command file for PGE_L2_HR_LakeTile for each input tile
        
        :param indf: index of current input PIXC file within the list of PIXC files
        :type indf: int
        
        :return: out_cmd_file = command file full path
        :rtype: string
        """

        # 1 - Get PIXC file full path
        pixc_file = os.path.join(self.pixc_dir, self.list_pixc[indf])
        
        # 2 - Retrieve tile info from the PIXC filename
        information = locnes_filenames.getInfoFromFilename(pixc_file, "PIXC")
        cycle_num = information["cycle"]
        pass_num = information["pass"]
        tile_ref = information["tile_ref"]

        # 3 - Init log filename (without date because it is computed in pge_lake_tile.py)
        if self.log_file is None:
            log_file = os.path.join(self.output_dir, "LogFile_" + str(cycle_num) + "_" +
                                    str(pass_num) + "_" + str(tile_ref) + ".log")
            print("Log file : " + log_file)
        else:
            log_file = os.path.splitext(self.log_file)[0] + "_" + str(cycle_num) + "_" \
                                        + str(pass_num) + "_" + str(tile_ref) + os.path.splitext(self.log_file)[1]
        
        # 4 - Init command filename
        cmd_filename = "lake_tile_command_" + str(cycle_num) + "_" + str(pass_num) + "_" + str(tile_ref) + ".cfg"
        out_cmd_file = os.path.join(self.output_dir, cmd_filename)
        
        # 5 - Write command variables in command file
        writer_command_file = open(out_cmd_file, "w")  # Open file in writing mode
        
        # 5.1 - Fill PATHS section
        writer_command_file.write("[PATHS]\n")
        if self.param_file is not None:
            writer_command_file.write("param_file = %s\n" % self.param_file)
        else:
            # needed by jenkins script
            writer_command_file.write("param_file = " + os.path.join(sys.path[0], "lake_tile_param.cfg") + "\n")

        writer_command_file.write("PIXC file = " + pixc_file + "\n")
        writer_command_file.write("PIXCVecRiver file = " + os.path.join(self.pixc_vec_river_dir, self.list_pixc_vec_river[indf]) + "\n")
        writer_command_file.write("Output directory = " + self.output_dir + "\n")
        writer_command_file.write("\n")
        
        # 5.2 - Fill DATABASES section
        writer_command_file.write("[DATABASES]\n")
        writer_command_file.write("# Lake a priori database\n")
        if self.lake_db is not None:
            writer_command_file.write("LAKE_DB = " + self.lake_db + "\n")
        writer_command_file.write("# Lake identifier attribute name in the database\n")
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
        
        # 5.5 - Fill CRID section
        writer_command_file.write("# CRID information\n")
        writer_command_file.write("[CRID]\n")
        writer_command_file.write("# Composite Release IDentifier for LakeTile processing\n")        
        writer_command_file.write("LAKE_TILE_CRID = Dx0000 \n")
        writer_command_file.write("# Composite Release IDentifier for LakeSP processing\n")
        writer_command_file.write("LAKE_SP_CRID = Dx0000 # Product generator\n")

        # 5.6 - Fill FILE_INFORMATION section
        writer_command_file.write("# File informations\n")
        writer_command_file.write("[FILE_INFORMATION]\n")
        writer_command_file.write("# Product generator\n")
        writer_command_file.write("PRODUCER = CNES\n")
        writer_command_file.write("PIXC_PREFIX = SWOT_L2_HR_PIXC_\n")
        writer_command_file.write('PIXC_PATTERN_PRINT = PIXC_PREFIX + "<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>.nc"\n')
        writer_command_file.write('# Indices when PIXC_PATTERN.split("_"); None if value not in filename\n')
        writer_command_file.write('PIXC_PATTERN_IND = {"cycle": 4, "pass": 5, "tile_ref": 6, "start_date": 7, "stop_date": 8, "crid": 9, "counter": 10} \n')
        writer_command_file.write('PIXCVEC_RIVER_PREFIX = SWOT_L2_HR_PIXCVecRiver_\n')
        writer_command_file.write('PIXCVEC_RIVER_PATTERN_PRINT = PIXCVEC_RIVER_PREFIX + "<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>.nc"\n')
        writer_command_file.write('# Indices when PIXCVEC_RIVER_PATTERN.split("_"); None if value not in filename\n')
        writer_command_file.write('PIXCVEC_RIVER_PATTERN_IND = {"cycle": 4, "pass": 5, "tile_ref": 6, "start_date": 7, "stop_date": 8, "crid": 9, "counter": 10}\n')
        writer_command_file.write('LAKE_TILE_PREFIX = SWOT_L2_HR_LakeTile_\n')
        writer_command_file.write("# LakeTile filename with %03d=cycle number %03d=pass number %s=tile ref %s=swath %s=begin date %s=end date %s=CRID %s=counter %s=suffix \n")
        writer_command_file.write('LAKE_TILE_PATTERN = LAKE_TILE_PREFIX + "%03d_%03d_%s_%s_%s_%s_%02d%s"\n')
        writer_command_file.write('LAKE_TILE_PATTERN_PRINT = LAKE_TILE_PREFIX + "%s<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>"\n')
        writer_command_file.write('# Indices when LAKE_TILE_*_PATTERN.split("_"); None if value not in filename \n')
        writer_command_file.write('LAKE_TILE_PATTERN_IND = {"cycle": 4, "pass": 5, "tile_ref": 6, "start_date": 7, "stop_date": 8, "crid": 9, "counter": 10}\n')
        writer_command_file.write('LAKE_TILE_SHP_SUFFIX = .shp\n')
        writer_command_file.write('LAKE_TILE_SHP_META_SUFFIX = .shp.xml\n')
        writer_command_file.write('LAKE_TILE_EDGE_SUFFIX = _edge.nc\n')
        writer_command_file.write('LAKE_TILE_PIXCVEC_SUFFIX = _pixcvec.nc\n')
        writer_command_file.write("PIXCVEC_PREFIX = SWOT_L2_HR_PIXCVec_\n")
        writer_command_file.write("PIXCVEC_SUFFIX = .nc\n")
        writer_command_file.write("# PIXCVec filename with %03d=cycle number %03d=pass number %s=tile ref %s=begin date %s=end date %s=CRID %s=counter  \n")
        writer_command_file.write('PIXCVEC_PATTERN = PIXCVEC_PREFIX + "%03d_%03d_%s_%s_%s_%s_%02d" + PIXCVEC_SUFFIX \n')
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
    print("[multiLakeTileProcessing] == read_command_file = %s ==" % in_filename)

    # 0 - Init output dictionary
    out_params = {}
    # Default values
    out_params["param_file"] = None
    out_params["LAKE_DB"] = None
    out_params["LAKE_DB_ID"] = None
    out_params["CONTINENT_FILE"] = None
    out_params["cycle_num"] = None
    out_params["pass_num"] = None
    out_params["tile_ref"] = None
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
    out_params["pixc_dir"] = config.get("PATHS", "PIXC directory")
    out_params["pixc_vec_river_dir"] = config.get("PATHS", "PIXCVecRiver directory")
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

    # 4 - Retrieve TILES_INFOS
    if "TILES_INFOS" in config.sections():
        list_options = config.options("TILES_INFOS")
        # Cycle number
        if "cycle number" in list_options:
            out_params["cycle_num"] = config.getint("TILES_INFOS", "Cycle number")
        # Pass number
        if "pass number" in list_options:
            out_params["pass_num"] = config.getint("TILES_INFOS", "Pass number")
        # Tile ref
        if "tile ref" in list_options:
            out_params["tile_ref"] = config.get("TILES_INFOS", "Tile ref")

    # 5 - Retrieve OPTIONS
    if "OPTIONS" in config.sections():
        list_options = config.options("OPTIONS")
        # Flag to also produce LakeTile_edge and LakeTile_pixcvec as shapefiles (=True); else=False (default)
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
    parser = argparse.ArgumentParser(description="Compute SWOT LakeTile products from multiple tiles of PIXC products and their associated PIXCVecRiver products.")
    parser.add_argument("command_file", help="command file (*.cfg)")
    args = parser.parse_args()

    print("===== multiLakeTileProcessing = BEGIN =====")
    print("")
    timer = my_timer.Timer()
    timer.start()

    # 1 - Read command file
    print("WORKING VARIABLES")
    print()
    my_tools.testFile(args.command_file, IN_extent=".cfg")  # Test existance and extension
    my_params = read_command_file(args.command_file)  # Read variables in command file

    # 2 - Initialization
    multi_lake_tile = MultiLakeTile(my_params)
    print(timer.info(0))

    # 3 - Run pre-processing
    multi_lake_tile.run_preprocessing()
    print(timer.info(0))

    # 4 - Run processing
    multi_lake_tile.run_processing()
    print(timer.info(0))

    print("")
    print(timer.stop())
    print("===== multiLakeTileProcessing = END =====")
