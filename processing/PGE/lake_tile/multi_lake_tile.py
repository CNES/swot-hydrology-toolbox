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
import multiprocessing as mp

import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_timer as my_timer
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
        self.pixcvec_river_dir = in_params["pixcvec_river_dir"]  # PIXCVecRiver files directory
        self.output_dir = in_params["output_dir"]  # Output directory

        # BDLac and basins file
        self.lake_db = in_params["LAKE_DB"]
        self.lake_db_id = in_params["LAKE_DB_ID"]

        # Tiles info
        self.cycle_num = in_params["cycle_num"]  # Cycle number
        self.pass_num = in_params["pass_num"]  # Pass number
        self.tile_ref = in_params["tile_ref"]  # Tile id

        # Flag to produce LakeTile_Edge and LakeTile_PIXCVec shapefiles
        self.flag_prod_shp = in_params["flag_prod_shp"]
        # Flag to increment output file counter
        self.flag_inc_file_counter = in_params["flag_inc_file_counter"]
        # Flag to write full path in global attributes
        self.flag_write_full_path = in_params["flag_write_full_path"]
        
        # Log level
        self.error_file = in_params["errorFile"]
        self.log_file = in_params["logFile"]
        self.log_file_level = in_params["logfilelevel"]
        self.log_console = in_params["logConsole"]
        self.log_console_level = in_params["logconsolelevel"]
        
        # File information
        self.institution = in_params["institution"]
        self.product_version = in_params["product_version"]
        self.crid = in_params["crid"]
        self.pge_version = in_params["pge_version"]
        self.contact = in_params["contact"]

        # List of input files
        self.list_pixc = []  # PIXC files
        self.list_pixcvec_river = []  # PIXCVecRiver
        self.nb_input = 0  # Nb of input files
        self.cmd_file_path_list = [] # list of command files

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
        my_tools.test_dir(self.pixc_dir)
        # 1.3 - PIXCVecRiver directory
        print("[multiLakeTileProcessing]   INPUT DIR for PIXCVecRiver files = %s" % self.pixcvec_river_dir)
        my_tools.test_dir(self.pixcvec_river_dir)
        # 1.4 - Output directory
        print("[multiLakeTileProcessing]   OUTPUT DIR = %s" % self.output_dir)
        my_tools.test_dir(self.output_dir)
        print("")

        # 2 - Get input files
        print("[multiLakeTileProcessing] > 2 - Retrieving input files...")

        # 2.1 - Compute file prefix regarding cycle / pass / tile conditions
        pixc_prefix = locnes_filenames.PIXC_PREFIX
        pixc_suffix = locnes_filenames.PIXC_SUFFIX
        pixcvec_river_prefix = locnes_filenames.PIXCVEC_RIVER_PREFIX
        cond_prefix = pixc_prefix  # Deal with all PixC files in self.pixc_dir 

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
        cur_pixcvec_river = None
        for cur_item in tmp_list:

            # Test if file meets the condition
            if cur_item.startswith(cond_prefix) and cur_item.endswith(pixc_suffix):  # Test if it's a wanted PIXC file

                # Associated PIXCVecRiver file name
                cur_pixcvec_river = cur_item.replace(pixc_prefix, pixcvec_river_prefix)

                # If associated PIXCVecRiver file exists, add pair of filenames
                if os.path.exists(os.path.join(self.pixcvec_river_dir, cur_pixcvec_river)):
                    self.list_pixc.append(cur_item)
                    self.list_pixcvec_river.append(cur_pixcvec_river)
                    self.nb_input += 1

        print("[multiLakeTileProcessing]   --> %d tile(s) to deal with" % self.nb_input)
        print("")

        for indf in range(self.nb_input):  # Deal with all selected files
            print("[multiLakeTileProcessing]  Writing command files %d / %d" % (indf+1, self.nb_input))
            cmd_file_path = self.create_cmd_file(indf)
            self.cmd_file_path_list.append(cmd_file_path)
            print("")

    def run_processing(self, cmd_file = None):
        """
        Process SAS_L2_HR_LakeTile for each input tile
        """

        print("")
        print("")
        print("[multiLakeTileProcessing] PROCESSING...")
        print("")
        print("")
        
        timer_proc = my_timer.Timer()
        timer_proc.start()  # Init timer
        
        if not cmd_file:
            for indf, cmd_file in enumerate(self.cmd_file_path_list) :
                print("")
                print("")
                print("***********************************************")
                print("***** Dealing with command file %d / %d *****" % (indf+1, len(self.cmd_file_path_list)))
                print("***********************************************")
                print("")
                print("")

                try:
                    # 1 - Init
                    my_lake_tile = pge_lake_tile.PGELakeTile(cmd_file)

                    # 2 - Run
                    my_lake_tile.start()

                    # 3 - Stop
                    my_lake_tile.stop()
                except:
                    print("WARNING : cmd_file %s failed" % cmd_file)

        else:
            print("")
            print("")
            print("***********************************************")
            print("***** Dealing with command file %s *****" % cmd_file)
            print("***********************************************")
            print("")
            print("")

            try:
                # 1 - Init
                my_lake_tile = pge_lake_tile.PGELakeTile(cmd_file)

                # 2 - Run
                my_lake_tile.start()

                # 3 - Stop
                my_lake_tile.stop()
            except:
                print("WARNING : cmd_file %s failed" % cmd_file)

        print("")
        print("")
        print(timer_proc.stop())  # Print tile process duration
        print("")
        print("")

    def run_multiprocessing(self):
        n_cores = int(mp.cpu_count())
        pool = mp.Pool(n_cores)
        with pool:
            print('Running map')
            print("[multiLakeTileProcessing] PROCESSING...")
            tmp_result = pool.map(self.run_processing, self.cmd_file_path_list)
            pool.close()
            pool.join()

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
        information = locnes_filenames.get_info_from_filename(pixc_file, "PIXC")
        cycle_num = information["cycle"]
        pass_num = information["pass"]
        tile_ref = information["tile_ref"]

        # 3.1 - Init error log filename (without date because it is computed in pge_lake_tile.py)
        if self.error_file is None:
            error_file = os.path.join(self.output_dir, "ErrorFile_" + str(cycle_num) + "_" +
                                    str(pass_num) + "_" + str(tile_ref) + ".log")
        else:
            error_file = os.path.splitext(self.error_file)[0] + "_" + str(cycle_num) + "_" \
                                        + str(pass_num) + "_" + str(tile_ref) + os.path.splitext(self.error_file)[1]
        print("Error log file: " + error_file)

        # 3.2 - Init log filename (without date because it is computed in pge_lake_tile.py)
        if self.log_file is None:
            log_file = os.path.join(self.output_dir, "LogFile_" + str(cycle_num) + "_" +
                                    str(pass_num) + "_" + str(tile_ref) + ".log")
        else:
            log_file = os.path.splitext(self.log_file)[0] + "_" + str(cycle_num) + "_" \
                                        + str(pass_num) + "_" + str(tile_ref) + os.path.splitext(self.log_file)[1]
        print("Log file: " + log_file)
        
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
        writer_command_file.write("PIXCVecRiver file = " + os.path.join(self.pixcvec_river_dir, self.list_pixcvec_river[indf]) + "\n")
        writer_command_file.write("Output directory = " + self.output_dir + "\n\n")
        
        # 5.2 - Fill DATABASES section
        writer_command_file.write("[DATABASES]\n")
        if self.lake_db is not None:
            writer_command_file.write("LAKE_DB = " + self.lake_db + "\n")
        if self.lake_db_id is not None:
            writer_command_file.write("LAKE_DB_ID = " + self.lake_db_id + "\n")
        writer_command_file.write("\n")
        
        # 5.3 - Fill OPTIONS section
        writer_command_file.write("[OPTIONS]\n")
        writer_command_file.write("Produce shp = " + str(self.flag_prod_shp) + "\n")
        writer_command_file.write("Increment file counter = " + str(self.flag_inc_file_counter) + "\n")
        writer_command_file.write("Write full path = " + str(self.flag_write_full_path) + "\n")
        writer_command_file.write("\n")
        
        # 5.4 - Fill LOGGING section
        writer_command_file.write("[LOGGING]\n")
        writer_command_file.write("errorFile = " + error_file + "\n")
        writer_command_file.write("logFile = " + log_file + "\n")
        writer_command_file.write("logfilelevel = " + self.log_file_level + "\n")
        writer_command_file.write("logConsole = " + str(self.log_console) + "\n")
        writer_command_file.write("logconsolelevel = " + self.log_console_level + "\n")
        writer_command_file.write("\n")

        # 5.6 - Fill FILE_INFORMATION section
        writer_command_file.write("[FILE_INFORMATION]\n")   
        writer_command_file.write("INSTITUTION = " + self.institution + "\n")
        writer_command_file.write("PRODUCT_VERSION = " + self.product_version + "\n")
        writer_command_file.write("CRID = " + self.crid + "\n")
        writer_command_file.write("PGE_VERSION = " + self.pge_version + "\n")
        writer_command_file.write("CONTACT = " + self.contact + "\n")
        writer_command_file.write("\n")
        
        # 5.7 - Close command file
        writer_command_file.close() 
        
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
    out_params["cycle_num"] = None
    out_params["pass_num"] = None
    out_params["tile_ref"] = None
    out_params["flag_prod_shp"] = False
    out_params["flag_inc_file_counter"] = True
    out_params["flag_write_full_path"] = False
    out_params["errorFile"] = None
    out_params["logFile"] = None
    out_params["logfilelevel"] = "DEBUG"
    out_params["logConsole"] = True
    out_params["logconsolelevel"] = "DEBUG"
    out_params["institution"] = "CNES"
    out_params["product_version"] = "0.0"
    out_params["crid"] = "Dx0000"
    out_params["pge_version"] = "0.0"
    out_params["contact"] = "test@cnes.fr"

    # 1 - Read parameter file
    config = cfg.ConfigParser()
    config.read(in_filename)

    # 2 - Retrieve PATHS
    list_paths = config.options("PATHS")
    if "param_file" in list_paths:
        out_params["param_file"] = config.get("PATHS", "param_file")
    out_params["pixc_dir"] = config.get("PATHS", "PIXC directory")
    out_params["pixcvec_river_dir"] = config.get("PATHS", "PIXCVecRiver directory")
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
        # Tile ref
        if "tile ref" in list_options:
            out_params["tile_ref"] = config.get("TILES_INFOS", "Tile ref")

    # 5 - Retrieve OPTIONS
    if "OPTIONS" in config.sections():
        list_options = config.options("OPTIONS")
        # Flag to also produce LakeTile_Edge and LakeTile_PIXCVec as shapefiles (=True); else=False (default)
        if "produce shp" in list_options:
            out_params["flag_prod_shp"] = config.getboolean("OPTIONS", "Produce shp")
        # Flag to increment the file counter in the output filenames (=True, default); else=False
        if "increment file counter" in list_options:
            out_params["flag_inc_file_counter"] = config.get("OPTIONS", "Increment file counter")
        # To write full path in global attributes (=True); to write only basename=False
        if "write full path" in list_options:
            out_params["flag_write_full_path"] = config.get("OPTIONS", "Write full path")

    # 6 - Retrieve LOGGING
    if "LOGGING" in config.sections():
        list_options = config.options("LOGGING")
        if "errorFile" in list_options:
            out_params["errorFile"] = config.get("LOGGING", "errorFile")
            if out_params["errorFile"] is not None:
                out_params["errorFile"] = out_params["errorFile"].replace("<date>", datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
        if "logFile" in list_options:
            out_params["logFile"] = config.get("LOGGING", "logFile")
            if out_params["logFile"] is not None:
                out_params["logFile"] = out_params["logFile"].replace("<date>", datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
        if "logfilelevel" in list_options:
            out_params["logfilelevel"] = config.get("LOGGING", "logfilelevel")
        if "logConsole" in list_options:
            out_params["logConsole"] = config.get("LOGGING", "logConsole")
        if "logconsolelevel" in list_options:
            out_params["logconsolelevel"] = config.get("LOGGING", "logconsolelevel")
            
    # 7 - Retrieve FILE_INFORMATION
    if "FILE_INFORMATION" in config.sections():
        list_options = config.options("FILE_INFORMATION")
        if "institution" in list_options:
            out_params["institution"] = config.get("FILE_INFORMATION", "INSTITUTION")
        if "product_version" in list_options:
            out_params["product_version"] = config.get("FILE_INFORMATION", "PRODUCT_VERSION")
        if "crid" in list_options:
            out_params["crid"] = config.get("FILE_INFORMATION", "CRID")
        if "pge_version" in list_options:
            out_params["pge_version"] = config.get("FILE_INFORMATION", "PGE_VERSION")
        if "contact" in list_options:
            out_params["contact"] = config.get("FILE_INFORMATION", "CONTACT")

    return out_params


if __name__ == '__main__':

    # 0 - Parse inline parameters
    parser = argparse.ArgumentParser(description="Compute SWOT LakeTile products from multiple tiles of PIXC products and their associated PIXCVecRiver products.")
    parser.add_argument("command_file", help="command file (*.cfg)")
    parser.add_argument("-mp", "--multiproc", help="if true, tiles will be computed in parallel", nargs='?', type=bool,
                        default=False, const=True)
    args = parser.parse_args()

    print("===== multiLakeTileProcessing = BEGIN =====")
    print("")
    timer = my_timer.Timer()
    timer.start()

    # 1 - Read command file
    print("WORKING VARIABLES")
    print()
    my_tools.test_file(args.command_file, in_extent=".cfg")  # Test existance and extension
    my_params = read_command_file(args.command_file)  # Read variables in command file

    # 2 - Initialization
    multi_lake_tile = MultiLakeTile(my_params)
    print(timer.info(0))

    # 3 - Run pre-processing
    multi_lake_tile.run_preprocessing()
    print(timer.info(0))

    # 4 - Run processing
    if args.multiproc:
        multi_lake_tile.run_multiprocessing()
    else:
        multi_lake_tile.run_processing()

    print(timer.info(0))

    print("")
    print(timer.stop())
    print("===== multiLakeTileProcessing = END =====")
