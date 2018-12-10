#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
.. module:: multi_lake_tile.py
    :synopsis: Process PGE_L2_HR_LakeTile (i.e. generate L2_HR_LakeTile product from one tile of L2_HR_PIXC product and associated L2_HR_PIXCVec product) for multiple tiles 
    08/23/2018 - Creation

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

Copyright (c) 2018 CNES. All rights reserved.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import os
import datetime
import argparse
import configparser as cfg

import cnes.common.lib.my_api as my_api
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_timer as my_timer
import cnes.common.lib_lake.locnes_variables as my_var
import pge_lake_tile


class Processing(object):

    def __init__(self, IN_params):
        """
        Constructor: initialize variables
        
        :param IN_params: input parameters to run the processor
        :type IN_params: dictionary
        """
        my_api.printInfo("[multiLakeTileProcessing] == INIT ==")

        # Directories
        self.pixc_dir = IN_params["pixc_dir"]  # PIXC files directory
        self.pixc_vec_river_dir = IN_params["pixc_vec_river_dir"]  # PIXCVecRiver files directory
        self.output_dir = IN_params["output_dir"]  # Output directory
        
        # Tiles info
        self.cycle_num = IN_params["cycle_num"]  # Cycle number
        self.pass_num = IN_params["pass_num"]  # Pass number
        self.tile_ref = IN_params["tile_ref"]  # Tile id

        # List of input files
        self.list_pixc = []  # PIXC files
        self.list_pixc_vec_river = []  # PIXCVecRiver
        self.nb_input = 0  # Nb of input files

        # Flag to produce LakeTile_edge and LakeTile_pixcvec shapefiles
        self.flag_prod_shp = IN_params["flag_prod_shp"]

    def run_preprocessing(self):

        """
        Retrieve the list of input files, i.e. L2_HR_PIXC tile file main and associated L2_HR_PIXCVec
        """

        my_api.printInfo("")
        my_api.printInfo("")
        my_api.printInfo("[multiLakeTileProcessing] PRE-PROCESSING...")
        my_api.printInfo("")

        # 1 - Test existence of directories
        my_api.printInfo("[multiLakeTileProcessing] > 1 - Testing existence of working directories...")
        # 1.1 - PIXC directory
        my_api.printInfo("[multiLakeTileProcessing]   INPUT DIR for PIXC files = %s" % self.pixc_dir)
        my_tools.testDir(self.pixc_dir)
        # 1.3 - PIXCVecRiver directory
        my_api.printInfo("[multiLakeTileProcessing]   INPUT DIR for PIXCVecRiver files = %s" % self.pixc_vec_river_dir)
        my_tools.testDir(self.pixc_vec_river_dir)
        # 1.4 - Output directory
        my_api.printInfo("[multiLakeTileProcessing]   OUTPUT DIR = %s" % self.output_dir)
        my_tools.testDir(self.output_dir)
        my_api.printInfo("")

        # 2 - Get input files
        my_api.printInfo("[multiLakeTileProcessing] > 2 - Retrieving input files...")

        # 2.1 - Compute file prefix regarding cycle / pass / tile conditions
        cond_prefix = my_var.PIXC_PREFIX  # Deal with all PixC files in self.pixc_dir
        if (self.cycle_num is None) or (self.cycle_num == "-1"):
            my_api.printInfo("[multiLakeTileProcessing]   All PixC files in the input directory")
        else:  # Deal with PixC files with cycle = self.cycle_num
            cond_prefix += "%03d" % self.cycle_num
            if (self.pass_num is None) or (self.pass_num == "-1"):
                my_api.printInfo("[multiLakeTileProcessing]   PixC files with cycle=%03d" % self.cycle_num)
            else:  # Deal with PixC files with cycle = self.cycle_num and pass = self.pass_num
                cond_prefix += "_%03d" % self.pass_num
                if self.tile_ref is not None:  # Deal with PixC files with cycle = self.cycle_num, pass = self.pass_num and tile id = self.tile_ref
                    my_api.printInfo("[multiLakeTileProcessing]   PixC files with cycle=%03d , pass=%03d , tile=%s" % (self.cycle_num, self.pass_num, self.tile_ref))
                    cond_prefix += "_%s" % self.tile_ref
                else:
                    my_api.printInfo("[multiLakeTileProcessing]   PixC files with cycle=%03d and pass=%03d" % (self.cycle_num, self.pass_num))

        # 2.2 - List all files in self.pixc_dir
        TMP_list = os.listdir(self.pixc_dir)
            
        # 2.3 - For each listed file, get related PIXCVecRiver files if they exist
        cur_pixc_vec_river = None
        for curItem in TMP_list:

            # Test if file meets the condition
            if curItem.startswith(cond_prefix):  # Test if it's a wanted PIXC file

                # Associated PIXCVecRiver file name
                cur_pixc_vec_river = curItem.replace(my_var.PIXC_PREFIX, my_var.PIXCVEC_RIVER_PREFIX)
                
                # If associated PIXCVecRiver file exists, add pair of filenames
                if os.path.exists(os.path.join(self.pixc_vec_river_dir, cur_pixc_vec_river)):  
                    self.list_pixc.append(curItem)
                    self.list_pixc_vec_river.append(cur_pixc_vec_river)
                    self.nb_input += 1

        my_api.printInfo("[multiLakeTileProcessing]   --> %d tile(s) to deal with" % self.nb_input)
        my_api.printInfo("")

    def run_processing(self):
        """
        Process SAS_L2_HR_LakeTile for each input tile
        """
        my_api.printInfo("")
        my_api.printInfo("")
        my_api.printInfo("[multiLakeTileProcessing] PROCESSING...")
        my_api.printInfo("")
        my_api.printInfo("")

        timer_proc = my_timer.Timer()
        timer_proc.start()

        if self.nb_input != 0:

            for indf in range(self.nb_input):  # Deal with all selected files

                my_api.printInfo("****************************************************************************************************************")
                my_api.printInfo("***** Dealing with tile %d / %d = %s *****" % (indf+1, self.nb_input, self.list_pixc[indf]))
                if my_api.GEN_ENV != 2:
                    print()
                    print("> Dealing with tile %d / %d = %s" % (indf+1, self.nb_input, self.list_pixc[indf]))
                my_api.printInfo("****************************************************************************************************************")
                my_api.printInfo("")
                my_api.printInfo("")
                
                # 1 - Initialization
                myLakeTile = pge_lake_tile.Processing(os.path.join(self.pixc_dir, self.list_pixc[indf]),
                                                      os.path.join(self.pixc_vec_river_dir, self.list_pixc_vec_river[indf]), 
                                                      self.output_dir, 
                                                      IN_shp_option=self.flag_prod_shp)
                my_api.printInfo(timer.info(0))
                
                # 2 - Run pre-processing
                myLakeTile.run_preprocessing()
                my_api.printInfo(timer.info(0))
                
                # 3 - Run processing
                myLakeTile.run_processing()
                my_api.printInfo(timer.info(0))
                
                # 4 - Run post-processing
                myLakeTile.run_postprocessing()
                my_api.printInfo(timer.info(0))
                
                my_api.printInfo("")
                my_api.printInfo(timer.stop())
                my_api.printInfo("")
                my_api.printInfo("")
                
            my_api.printInfo("****************************************************************************************************************")
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
    print("[multiLakeTileProcessing] == readParamFile = %s ==" % IN_filename)
    
    # 0 - Init output dictionary
    OUT_params = {}
    # Default values
    OUT_params["cycle_num"] = None
    OUT_params["pass_num"] = None
    OUT_params["tile_ref"] = None
    OUT_params["flag_prod_shp"] = False
    
    # 1 - Read parameter file
    config = cfg.ConfigParser()
    config.read(IN_filename)
    
    # 2 - Retrieve PATHS
    OUT_params["pixc_dir"] = config.get("PATHS", "PIXC directory")
    OUT_params["pixc_vec_river_dir"] = config.get("PATHS", "PIXCVecRiver directory")
    OUT_params["output_dir"] = config.get("PATHS", "Output directory")
    
    # 3 - Retrieve tiles infos
    if "TILES_INFOS" in config.sections():
        list_options = config.options("TILES_INFOS")
        # Cycle number
        if "cycle number" in list_options:
            OUT_params["cycle_num"] = config.getint("TILES_INFOS", "Cycle number")
        # Pass number
        if "pass number" in list_options:
            OUT_params["pass_num"] = config.getint("TILES_INFOS", "Pass number")
        # Tile ref
        if "tile ref" in list_options:
            OUT_params["tile_ref"] = config.get("TILES_INFOS", "Tile ref")
    
    # 4 - Retrieve OPTIONS
    if "OPTIONS" in config.sections():
        list_options = config.options("OPTIONS")
        # Flag to also produce LakeTile_edge and LakeTile_pixcvec as shapefiles (=True); else=False (default)
        if "produce shp" in list_options:
            OUT_params["flag_prod_shp"] = config.getboolean("OPTIONS", "Produce shp")
            
    # 5 - Retrieve optionnal CONFIG_OVERWRITE
    my_var.overwriteConfig_from_cfg(config)
    
    print()
            
    return OUT_params


#######################################


if __name__ == '__main__':
    
    # 0 - Parse inline parameters
    parser = argparse.ArgumentParser(description="Compute SWOT LakeTile products from multiple tiles of PIXC products and their associated PIXCVecRiver products.")
    parser.add_argument("param_file", help="parameter file (*.cfg)")
    parser.add_argument("-l", "--logfile", help="write prints to a logfile", action="store_true")
    parser.add_argument("-v", "--verbose", help="verbose level", choices=["DEBUG", "INFO"], default="INFO")
    args = parser.parse_args()

    print("===== multiLakeTileProcessing = BEGIN =====")
    print("")
    timer = my_timer.Timer()
    timer.start()

    # 1 - Working variables
    print("WORKING VARIABLES")
    print()
    logFile = None  # Log file full path
    
    # 1.1 - Init variables
    log_option = False  # To print logs on screen (=False, default) or in a logfile (=True)
    verbose_option = "INFO"  # Verbose level
    
    # 1.2 - Read parameter file
    my_tools.testFile(args.param_file, IN_extent=".cfg")  # Test existance and extention
    my_params = readParamFile(args.param_file)  # Read parameters
    
    # 1.3 - Init environment for verbose level
    verbose_level = my_api.setVerbose(args.verbose)
    print("> Verbose level = %s" % verbose_level)
        
    # 1.4 - Init environment for log
    if args.logfile:
        logFile = os.path.join(my_params["output_dir"], "multi_lake_tile_" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + ".log")
        my_api.initLogger(logFile, verbose_level)
        print("> Log file = %s" % logFile)
    else:
        # Init log system for compatibilies with new log system
        logFormatter = logging.Formatter(fmt='%(asctime)s.%(msecs)03d     %(levelname)s:%(name)s::%(funcName)s: %(message)s',datefmt='%Y-%m-%dT%H:%M:%S')
        rootLogger = logging.getLogger()
        rootLogger.setLevel("DEBUG")
        consoleHandler = logging.StreamHandler()
#        consoleHandler.setFormatter(logFormatter)
        rootLogger.addHandler(consoleHandler)

        print("> No log file ; print info on screen")
        print()
        print()

    # 2 - Initialization
    myMultiLakeTile = Processing(my_params)
    my_api.printInfo(timer.info(0))

    # 2 - Run pre-processing
    myMultiLakeTile.run_preprocessing()
    my_api.printInfo(timer.info(0))

    # 3 - Run processing
    myMultiLakeTile.run_processing()
    my_api.printInfo(timer.info(0))

    my_api.printInfo("")
    my_api.printInfo(timer.stop())
    
    # Close logger
    if args.logfile:
        my_api.closeLogger()

    print("")
    print(timer.stop())
    print("===== multiLakeTileProcessing = END =====")
