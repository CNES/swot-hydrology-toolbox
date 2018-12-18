#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
.. module:: multi_lake_sp.py
    :synopsis: Process PGE_L2_HR_LakeSP (i.e. generate L2_HR_LakeSP product from the tiles of L2_HR_LakeTile products related to 1 single pass) for multiple passes
    Created on 08/12/2018

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

import cnes.common.lib.my_api as my_api
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_timer as my_timer
import cnes.common.lib_lake.locnes_filenames as my_names
import cnes.common.lib_lake.locnes_variables as my_var
import pge_lake_sp


class Processing(object):

    def __init__(self, IN_laketile_dir, IN_output_dir, IN_cycle_num, IN_pass_num, IN_shp_option=False):
        """
        Constructor: initialize variables
        
        :param IN_laketile_dir: full path of directory containing LakeTile products
        :type IN_laketile_dir: string
        :param IN_output_dir: output directory full path
        :type IN_output_dir: string
        :param IN_cycle_num: cycle number
        :type IN_cycle_num: int
        :param IN_pass_num: pass number
        :type IN_pass_num: int
        :param IN_shp_option: to also produce PIXCVec as shapefile (=True); else=False (default)
        :type IN_shp_option: boolean
        """
        my_api.printInfo("[multiLakeSpProcessing] == INIT ==")

        # Input paths
        self.lake_tile_dir = IN_laketile_dir  # Full path of directory containing LakeTile products
        self.output_dir = IN_output_dir  # Output directory
        
        # Pass infos
        self.cycle_num = IN_cycle_num  # Cycle number
        self.pass_num = IN_pass_num  # Pass number

        # Flag to produce PIXCVec shapefile
        self.flag_prod_shp = False
        if IN_shp_option:
            self.flag_prod_shp = True
            
        # Init variables
        self.list_cycle = []  # List of cycle numbers to deal with
        self.list_pass = []  # List of pass numbers to deal with

    def run_preprocessing(self):
        """
        Retrieve the list of input files, i.e. L2_HR_LakeTile products for wanted passes = shapefile + PIXC_edge + PIXCVec files
        """

        my_api.printInfo("")
        my_api.printInfo("")
        my_api.printInfo("[multiLakeSPProcessing] PRE-PROCESSING...")
        my_api.printInfo("")

        # 1 - Test existence of directories
        my_api.printInfo("[multiLakeSPProcessing] > 1 - Testing existence of working directories...")
        # 1.1 - LakeTile directory
        my_api.printInfo("[multiLakeSPProcessing]   INPUT LakeTile DIR = %s" % self.lake_tile_dir)
        my_tools.testDir(self.lake_tile_dir)
        # 1.2 - Output directory
        my_api.printInfo("[multiLakeSPProcessing]   OUTPUT DIR = %s" % self.output_dir)
        my_tools.testDir(self.output_dir)
        my_api.printInfo("")

        # 2 - Get input files
        my_api.printInfo("[multiLakeSPProcessing] > 2 - Retrieving input files...")

        # 2.1 - Compute file prefix regarding cycle / pass / tile conditions
        cond_prefix = my_var.LAKE_TILE_PREFIX  # Deal with all LakeTile products in self.lake_tile_dir
        if (self.cycle_num is None) or (self.cycle_num == "-1"):
            my_api.printInfo("[multiLakeSPProcessing]   All LakeTile files in the input directory")
        else:  # Deal with LakeTile files with cycle = self.cycle_num
            cond_prefix += "%03d" % self.cycle_num
            if (self.pass_num is None) or (self.pass_num == "-1"):
                my_api.printInfo("[multiLakeSPProcessing]   LakeTile files with cycle=%03d" % self.cycle_num)
            else:  # Deal with LakeTile files with cycle = self.cycle_num and pass = self.pass_num
                cond_prefix += "_%03d" % self.pass_num
                my_api.printInfo("[multiLakeSPProcessing]   LakeTile files with cycle=%03d and pass=%03d" % (self.cycle_num, self.pass_num))

        # 2.2 - List all files in input directory
        TMP_list = os.listdir(self.lake_tile_dir)
            
        # 2.3 - For each listed file, get (cycle, pass) pair
        for curFile in TMP_list:

            # Test if file meets the condition
            if curFile.startswith(cond_prefix) and curFile.endswith(my_var.LAKE_TILE_SHP_META_SUFFIX):  # Test if it's a wanted LakeTile_shp file

                TMP_infos = my_names.getInfoFromFilename(curFile, "LakeTile")
                if int(TMP_infos["pass"]) in self.list_pass:
                    TMP_ind = [indice for indice, valeur in enumerate(self.list_pass) if valeur==int(TMP_infos["pass"])]  # All occurrences of current pass number in the list of passes
                    TMP_subset_cycle = [self.list_cycle[ind] for ind in TMP_ind]  # Subset of cycle numbers related to current pass
                    if int(TMP_infos["cycle"]) not in TMP_subset_cycle:  # If current cycle not in subset = (cycle, pass) pair not listed
                        self.list_cycle.append(int(TMP_infos["cycle"]))
                        self.list_pass.append(int(TMP_infos["pass"]))
                else:
                    self.list_cycle.append(int(TMP_infos["cycle"]))
                    self.list_pass.append(int(TMP_infos["pass"]))

        my_api.printInfo("[multiLakeSPProcessing]   --> %d (cycle, pass) pair(s) to deal with" % len(self.list_cycle))
        my_api.printInfo("")

    def run_processing(self):
        """
        Process SAS_L2_HR_LakeSP for each (cycle, pass) pair
        """
        my_api.printInfo("")
        my_api.printInfo("")
        my_api.printInfo("[multiLakeSPProcessing] PROCESSING...")
        my_api.printInfo("")
        my_api.printInfo("")

        timer_proc = my_timer.Timer()
        timer_proc.start()

        if len(self.list_pass) != 0:

            for cur_cycle, cur_pass in zip(self.list_cycle, self.list_pass):  # Deal with all selected pairs

                my_api.printInfo("*******************************************************************")
                my_api.printInfo("*****           Dealing with cycle %03d and pass %03d           *****" % (cur_cycle, cur_pass))
                if my_api.GEN_ENV != 2:
                    print()
                    print("> Dealing with cycle %03d and pass %03d" % (cur_cycle, cur_pass))
                my_api.printInfo("*******************************************************************")
                my_api.printInfo("")
                my_api.printInfo("")
                
                # 1 - Initialization
                myLakeSP = pge_lake_sp.Processing(self.lake_tile_dir, self.output_dir, 
                                                  cur_cycle, cur_pass, 
                                                  IN_shp_option=self.flag_prod_shp)
                my_api.printInfo(timer.info(0))
                
                # 2 - Run pre-processing
                myLakeSP.run_preprocessing()
                my_api.printInfo(timer.info(0))
                
                # 3 - Run processing
                myLakeSP.run_processing()
                my_api.printInfo(timer.info(0))
                
                # 4 - Run post-processing
                myLakeSP.run_postprocessing()
                my_api.printInfo(timer.info(0))
                
                my_api.printInfo("")
                my_api.printInfo(timer.stop())
                my_api.printInfo("")
                my_api.printInfo("")
                
            
            my_api.printInfo("*******************************************************************")
            my_api.printInfo("")
            my_api.printInfo("")


#######################################
            
            
def readParamFile(IN_filename):
    """
    Read the parameter file in input and store parameters in a dictionary
    
    :param IN_filename: parameter file full path
    :type IN_filename: string
    
    :return: dictionary containing parameters
    :rtype: dict
    """
    print("[multiLakeSPProcessing] == readParamFile = %s ==" % IN_filename)
    
    # 0 - Init output dictionary
    OUT_params = {}
    # Default values
    OUT_params["cycle_num"] = None
    OUT_params["pass_num"] = None
    OUT_params["flag_prod_shp"] = False
    
    # 1 - Read parameter file
    config = cfg.ConfigParser()
    config.read(IN_filename)
    
    # 2 - Retrieve PATHS
    OUT_params["laketile_dir"] = config.get("PATHS", "LakeTile directory")
    OUT_params["output_dir"] = config.get("PATHS", "Output directory")
    
    #• 3 - Retrieve TILES_INFOS
    if "TILES_INFOS" in config.sections():
        list_infos = config.options("TILES_INFOS")
        if "cycle number" in list_infos:
            OUT_params["cycle_num"] = config.getint("TILES_INFOS", "Cycle number")
        if "pass number" in list_infos:
            OUT_params["pass_num"] = config.getint("TILES_INFOS", "Pass number")
    
    # 4 - Retrieve OPTIONS
    if "OPTIONS" in config.sections():
        list_options = config.options("OPTIONS")
        # Flag to also produce PIXCVec file as shapefile (=True); else=False (default)
        if "produce shp" in list_options:
            OUT_params["flag_prod_shp"] = config.getboolean("OPTIONS", "Produce shp")
    
    return OUT_params


#######################################


if __name__ == '__main__':
    
    # 0 - Parse inline parameters
    parser = argparse.ArgumentParser(description="Compute multiple SWOT LakeSP products from LakeTile products corresponding to one or more specific (cycle, pass). \
                                     If indir_or_param_file is a parameter file (*.cfg), all input parameters are only read in the parameter file.")
    parser.add_argument("indir_or_param_file", help="LakeTile directory or parameter file (*.cfg)")
    parser.add_argument("output_dir", help="output directory", nargs='?')
    parser.add_argument("-cycle", help="cycle number", type=int)
    parser.add_argument("-pass", help="pass number", type=int)
    parser.add_argument("-shp", help="convert output NetCDF file as shapefile", action="store_true")
    parser.add_argument("-l", "--logfile", help="write prints to a logfile", action="store_true")  # To print logs on screen (=False, default) or in a logfile (=True)
    parser.add_argument("-v", "--verbose", help="verbose level", choices=["DEBUG", "INFO"], default="INFO")  # Verbose level
    args = parser.parse_args()
    
    print("===== multiLakeSPProcessing = BEGIN =====")
    print("")
    timer = my_timer.Timer()
    timer.start()
    
    # 1 - Working variables
    
    # 1.1 - Init variables
    paramFile = None  # Parameter file full path
    laketile_dir = None  # Input directory, containing LakeTile products
    output_dir = None  # Output directory
    cycle_num = None  # Cycle number
    pass_num = None  # Pass number
    shp_option = False  # To also produce PIXCVec file as shapefile (=True); else=False (default)
    
    # 1.2 - Read values according to indir_or_param_file values
    if os.path.isdir(args.indir_or_param_file):  # Read inline parameters
        location = "inline command"
        laketile_dir = args.indir_or_param_file 
        output_dir = args.output_dir
        cycle_num = args.cycle_num
        pass_num = args.pass_num
        shp_option = args.shp
        
    else:
        
        file_base, file_extent = os.path.splitext(args.indir_or_param_file)
        
        if file_extent == ".cfg":  # Read parameter file
            location = "parameter file"
            my_tools.testFile(args.indir_or_param_file, IN_extent=".cfg")  # Test existance and extension
            my_params = readParamFile(args.indir_or_param_file)  # Read parameters
            laketile_dir = my_params["laketile_dir"]
            output_dir = my_params["output_dir"]
            cycle_num = my_params["cycle_num"]
            pass_num = my_params["pass_num"]
            shp_option = my_params["flag_prod_shp"]
        
        else:
            print("[ERROR]")
            print("Run by multi_lake_sp.py param_file.cfg [-l] [-v VERBOSE]")
            print("OR multi_lake_sp.py lake_tile_dir output_dir [cycle_num [pass_num]] [-shp] [-l] [-v VERBOSE]")
            sys.exit("indir_or_param_file is %s, not .cfg" % file_extent)
    
    # 1.3 - Test input params have been filled
    # 1.3.1 - LakeTile directory
    if laketile_dir is None:
        my_api.exitWithError("LakeTile directory is missing in %s" % location)
    # 1.3.2 - Output directory
    if output_dir is None:
        my_api.exitWithError("Output directory is missing in %s" % location)
    my_tools.testDir(output_dir)  # Test existence of output directory here and not in pre-proc because used in 1.5
    
    # 1.4 - Init environment for verbose level
    verbose_level = my_api.setVerbose(args.verbose)
    print("> Verbose level = %s" % verbose_level)
        
    # 1.5 - Init environment for log
    if args.logfile:
        logFile = os.path.join(output_dir, "multi_lake_sp_" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + ".log")
        my_api.initLogger(logFile, verbose_level)
        print("> Log file = %s" % logFile)
    else:
        print("> No log file ; print info on screen")
        print()
        print()

    # 2 - Initialization
    myLakeSP = Processing(laketile_dir, output_dir, cycle_num, pass_num, IN_shp_option=shp_option)
    my_api.printInfo(timer.info(0))

    # 3 - Run pre-processing
    myLakeSP.run_preprocessing()
    my_api.printInfo(timer.info(0))

    # 4 - Run processing
    myLakeSP.run_processing()
    my_api.printInfo(timer.info(0))

    my_api.printInfo("")
    my_api.printInfo(timer.stop())
    
    # Close logger
    if args.logfile:
        my_api.closeLogger()

    print("")
    print((timer.stop()))
    print("===== multiLakeSPProcessing = END =====")
