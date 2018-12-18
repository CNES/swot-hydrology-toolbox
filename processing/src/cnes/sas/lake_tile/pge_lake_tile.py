#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
.. module:: pge_lake_tile.py
    :synopsis: Process PGE_L2_HR_LakeTile, i.e. generate L2_HR_LakeTile product from one tile of L2_HR_PIXC product and associated L2_HR_PIXCVec product
    Created on 02/27/2017

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
import logging

import cnes.sas.lake_tile.proc_pixc as proc_pixc
import cnes.common.lib.my_api as my_api
import cnes.common.lib.my_shp_file as my_shp
import cnes.common.lib.my_timer as my_timer
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib_lake.lake_db as lake_db
import cnes.common.lib_lake.locnes_filenames as my_names
import cnes.common.lib_lake.locnes_variables as my_var
import cnes.common.lib_lake.proc_lake as proc_lake
import cnes.common.lib_lake.proc_pixc_vec as proc_pixc_vec


class Processing(object):

    def __init__(self, IN_pixc_file, IN_pixc_vec_river_file, IN_output_dir, IN_shp_option=False):
        """
        Constructor: initialize variables
        
        :param IN_pixc_file: PIXC file full path
        :type IN_pixc_file: string
        :param IN_pixc_vec_river_file: PIXCVecRiver file full path
        :type IN_pixc_vec_river_file: string
        :param IN_output_dir: output directory full path
        :type IN_output_dir: string
        :param IN_shp_option: to also produce LakeTile_edge and LakeTile_pixcvec as shapefiles (=True); else=False (default)
        :type IN_shp_option: boolean
        """
        my_api.printInfo("[lakeTileProcessing] == INIT ==")

        # Input paths
        self.pixc_file = IN_pixc_file  # PixC file
        self.pixc_vec_river_file = IN_pixc_vec_river_file  # Associated PIXCVecRiver file
        self.output_dir = IN_output_dir  # Output directory

        # Flag to produce LakeTile_edge and LakeTile_pixcvec shapefiles
        self.flag_prod_shp = False
        if IN_shp_option:
            self.flag_prod_shp = True
        
        # LakeTile filenames
        self.lake_tile_filenames = None
        
        # Objects
        self.objLakeDb = None  # Lake DB object
        self.objPixc = None  # PIXC object
        self.objPixcVec = None  # PIXCVecRiver object
        self.objLake = None  # LakeTile object

    def run_preprocessing(self):

        """
        Process PGE_LakeTile IN, i.e. test input paths, retrieve orbit infos, open lake database and init objects
        """

        my_api.printInfo("")
        my_api.printInfo("")
        my_api.printInfo("[lakeTileProcessing] PRE-PROCESSING...")
        my_api.printInfo("")
        
        # 1 - Test existance and file format of input paths
        my_api.printInfo("[lakeTileProcessing] > 1 - Testing existence of input paths...")
        # 1.1 - PIXC file
        my_api.printInfo("[lakeTileProcessing]   INPUT PIXC file = %s" % self.pixc_file)
        my_tools.testFile(self.pixc_file, IN_extent=".nc")
        # 1.2 - PIXCVecRiver file
        my_api.printInfo("[lakeTileProcessing]   INPUT PIXCVecRiver file = %s" % self.pixc_vec_river_file)
        my_tools.testFile(self.pixc_vec_river_file, IN_extent=".nc")
        # 1.3 - Output directory
        my_api.printInfo("[lakeTileProcessing]   OUTPUT DIR = %s" % self.output_dir)
        my_tools.testDir(self.output_dir)
        my_api.printInfo("") 
        
        # 2 - Retrieve orbit info from PIXC filename and compute output filenames
        my_api.printInfo("[lakeTileProcessing] > 2 - Retrieving tile infos from PIXC filename...")
        self.lake_tile_filenames = my_names.lakeTileFilenames(self.pixc_file, self.pixc_vec_river_file, self.output_dir)
        my_api.printInfo("")  
        
        # 3 - Objects initialisation
        my_api.printInfo("[lakeTileProcessing] > 3 - Init and format intput objects...")
        my_api.printInfo("")
        
        # 3.1 - Init PIXCVec product by retrieving data from the pixel cloud complementary file after river processing
        my_api.printInfo("[lakeTileProcessing] > 3a - Init pixel cloud complementary file...")
        self.objPixcVec = proc_pixc_vec.PixelCloudVec("TILE", self.pixc_vec_river_file)
        my_api.printInfo("")
        
        # 3.2 - Retrieve needed data from the pixel cloud
        my_api.printInfo("[lakeTileProcessing] > 3b - Retrieving needed data from the pixel cloud...")
        self.objPixc = proc_pixc.PixelCloud(self.pixc_file, self.objPixcVec.reject_idx)
        my_api.printInfo("")
        
        # 3.3 - Reshape PIXCVec arrays
        my_api.printInfo("[lakeTileProcessing] > 3c - Reshape PIXCVecRiver arrays...")
        self.objPixcVec.reshape(self.objPixc)
        my_api.printInfo("")   

        # 4 - Retrieve lake Db layer
        my_api.printInfo("[lakeTileProcessing] > 4 - Retrieving lake database layer...")
        if my_var.LAKE_DB == "":
            my_api.printInfo("[lakeTileProcessing] NO database specified -> NO link of SWOT obs with a priori lake")
        else:
            if os.path.exists(my_var.LAKE_DB):
                type_db = my_var.LAKE_DB.split('.')[-1]  # Type of database
                if type_db == "shp":  # Shapefile format
                    self.objLakeDb = lake_db.LakeDb_shp(my_var.LAKE_DB, self.objPixc.tile_poly)
                elif type_db == "sqlite":  # SQLite format
                    self.objLakeDb = lake_db.LakeDb_sqlite(my_var.LAKE_DB, self.objPixc.tile_poly)
                else:
                    my_api.exitWithError("[lakeTileProcessing] Lake a priori database format (%s) is unknown: must be .shp or .sqlite" % type_db)
            else:
                my_api.exitWithError("[lakeTileProcessing]   ERROR = %s doesn't exist" % my_var.LAKE_DB)
        my_api.printInfo("")
        
        # 4 - Initialize lake product
        my_api.printInfo("[lakeTileProcessing] > 5 - Init lake product...")
        self.objLake = proc_lake.LakeProduct("TILE",
                                             self.objPixc,
                                             self.objPixcVec,
                                             self.objLakeDb,
                                             os.path.basename(self.lake_tile_filenames.lake_tile_shp_file).split(".")[0],
                                             IN_id_prefix=self.lake_tile_filenames.lake_id_prefix)
        my_api.printInfo("")  

    def run_processing(self):
        """
        Process SAS_L2_HR_LakeTile
        """
        my_api.printInfo("")
        my_api.printInfo("")
        my_api.printInfo("[lakeTileProcessing] PROCESSING...")
        my_api.printInfo("")

        timer_proc = my_timer.Timer()
        timer_proc.start()

        # Processing only if PixC pixels are selected
        if self.objPixc.nb_selected != 0:

            # 2 - F2-F3-F3b = Identify all separate entities in the water mask
            my_api.printInfo("[lakeTileProcessing] 1 - Identifying all separate entities in the water mask...")
            self.objPixc.computeSeparateEntities()
            my_api.printInfo("[lakeTileProcessing] " + timer_proc.info(0))
            my_api.printInfo("")

            # 3 - F4 = Retrieve pixels corresponding to lakes and new objects entirely inside the tile
            my_api.printInfo("[lakeTileProcessing] 2 - Getting pixels corresponding to lakes and new objects entirely inside the tile...")
            self.objPixc.computeObjInsideTile()
            my_api.printInfo("[lakeTileProcessing] " + timer_proc.info(0))
            my_api.printInfo("")

            # 4 - F6 = Fill lake product
            my_api.printInfo("[lakeTileProcessing] 3 - Filling LakeTile product...")
            self.objLake.computeLakeProducts(self.objPixc.labels_inside)
            my_api.printInfo("[lakeTileProcessing] " + timer_proc.info(0))
            my_api.printInfo("")

        else:
            my_api.printInfo("[lakeTileProcessing] NO selected PixC => empty lake tile product generated")
            my_api.printInfo("")

    def run_postprocessing(self):
        """
        Process PGE_L2_HR_LakeTile OUT, i.e. convert output data in L2_HR_LakeTile product and close files
        """
        my_api.printInfo("")
        my_api.printInfo("")
        my_api.printInfo("[lakeTileProcessing] POST-PROCESSING...")
        my_api.printInfo("")

        # 1 - Write LakeTile shapefile
        my_api.printInfo("[lakeTileProcessing] 1 - Writing LakeTile memory layer to shapefile...")
        my_shp.write_mem_layer_as_shp(self.objLake.layer, self.lake_tile_filenames.lake_tile_shp_file)
        self.objLake.dataSource.Destroy()  # Close memory layer
        my_api.printInfo("")
        # Write XML metadatafile for shapefile
        self.objLake.writeMetadataFile("%s.xml" % self.lake_tile_filenames.lake_tile_shp_file)

        # 2 - Write PIXCVec for objects entirely inside tile
        my_api.printInfo("[lakeTileProcessing] 2 - Writing LakeTile_pixcvec file...")
        self.objPixcVec.write_file(self.lake_tile_filenames.lake_tile_pixcvec_file, None, True)
        if self.flag_prod_shp and (self.objPixcVec.nb_water_pix != 0):
            self.objPixcVec.write_file_asShp(self.lake_tile_filenames.lake_tile_pixcvec_file.replace(".nc", ".shp"), IN_classif=self.objPixc.origin_classif)
        my_api.printInfo("")

        # 3 - Write intermediate NetCDF file with indices of pixels (and associated label) related to objects at the top/bottom edges of the tile
        my_api.printInfo("[lakeTileProcessing] 3 - Writing LakeTile_edge file...")
        self.objPixc.write_edge_file(self.lake_tile_filenames.lake_tile_edge_file, None, True)
        if self.flag_prod_shp and (self.objPixc.nb_edge_pix != 0):
            self.objPixc.write_edge_file_asShp(self.lake_tile_filenames.lake_tile_edge_file.replace(".nc", ".shp"))
        my_api.printInfo("")

        # 4 - Close lake database
        if self.objLakeDb is not None:
            my_api.printInfo("[lakeTileProcessing] 4 - Closing lake database...")
            self.objLakeDb.close_db()
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
    print("[lakeTileProcessing] == readParamFile = %s ==" % IN_filename)
    
    # 0 - Init output dictionary
    OUT_params = {}
    # Default values
    OUT_params["flag_prod_shp"] = False
    
    # 1 - Read parameter file
    config = cfg.ConfigParser()
    config.read(IN_filename)
    
    # 2 - Retrieve PATHS
    OUT_params["pixc_file"] = config.get("PATHS", "PIXC file")
    OUT_params["pixc_vec_river_file"] = config.get("PATHS", "PIXCVecRiver file")
    OUT_params["output_dir"] = config.get("PATHS", "Output directory")
    
    # 3 - Retrieve OPTIONS
    if "OPTIONS" in config.sections():
        list_options = config.options("OPTIONS")
        # Flag to also produce LakeTile_edge and LakeTile_pixcvec as shapefiles (=True); else=False (default)
        if "produce shp" in list_options:
            OUT_params["flag_prod_shp"] = config.getboolean("OPTIONS", "Produce shp")
            
    # 4 - Retrieve optionnal CONFIG_OVERWRITE
    my_var.overwriteConfig_from_cfg(config)
    
    print()
    
    return OUT_params


#######################################


if __name__ == '__main__':
    
    # 0 - Parse inline parameters
    parser = argparse.ArgumentParser(description="Compute SWOT LakeTile product from a PIXC product and its associated PIXCVecRiver product. \
                                     If pixc_or_param_file is a parameter file (*.cfg), all input parameters are only read in the parameter file.")
    parser.add_argument("pixc_or_param_file", help="PIXC file (*.nc) or parameter file (*.cfg)")
    parser.add_argument("pixc_vec_river_file", help="associated PIXCVecRiver file (*.nc)", nargs='?')
    parser.add_argument("output_dir", help="output directory", nargs='?')
    parser.add_argument("-shp", help="convert output NetCDF file as shapefile", action="store_true")
    parser.add_argument("-l", "--logfile", help="write prints to a logfile", action="store_true")  # To print logs on screen (=False, default) or in a logfile (=True)
    parser.add_argument("-v", "--verbose", help="verbose level", choices=["DEBUG", "INFO"], default="INFO")  # Verbose level
    args = parser.parse_args()

    print("===== lakeTileProcessing = BEGIN =====")
    print("")
    timer = my_timer.Timer()
    timer.start()

    # 1 - Working variables
    print("WORKING VARIABLES")
    print()
    
    # 1.1 - Init variables
    paramFile = None  # Parameter file full path
    pixc_file = None  # PIXC file full path
    pixc_vec_river_file = None  # PIXCVecRiver file full path
    output_dir = None  # Output directory
    shp_option = False  # To also produce LakeTile_edge and LakeTile_pixcvec as shapefiles (=True); else=False (default)
    
    # 1.2 - Read values according to pixc_or_param_file extent
    file_base, file_extent = os.path.splitext(args.pixc_or_param_file)
    if file_extent == ".cfg":  # Read parameter file
        location = "parameter file"
        my_tools.testFile(args.pixc_or_param_file, IN_extent=".cfg")  # Test existance and extention
        my_params = readParamFile(args.pixc_or_param_file)  # Read parameters
        pixc_file = my_params["pixc_file"]
        pixc_vec_river_file = my_params["pixc_vec_river_file"]
        output_dir = my_params["output_dir"]
        shp_option = my_params["flag_prod_shp"]
        
    elif file_extent == ".nc":  # Read inline parameters
        location = "inline command"
        pixc_file = args.pixc_or_param_file 
        pixc_vec_river_file = args.pixc_vec_river_file
        output_dir = args.output_dir
        shp_option = args.shp
        
    else:
        print("[ERROR]")
        print("Run by pge_lake_tile.py param_file.cfg [-l] [-v VERBOSE]")
        print("OR pge_lake_tile.py pixc_file.nc pixc_vec_river_file.nc output_dir [-shp] [-l] [-v VERBOSE]")
        sys.exit("pixc_or_param_file is %s" % file_extent)
    
    # 1.3 - Test input paths have been filled
    # 1.3.1 - PIXC file
    if pixc_file is None:
        my_api.exitWithError("PIXC file is missing in %s" % location)
    # 1.3.2 - PIXCVecRiver file
    if pixc_vec_river_file is None:
        my_api.exitWithError("PIXCVecRiver file is missing in %s" % location)
    # 1.3.3 - Output directory
    if output_dir is None:
        my_api.exitWithError("Output directory is missing in %s" % location)
    my_tools.testDir(output_dir)  # Test existence of output directory here and not in pre-proc because used in 1.5
    
    # 1.4 - Init environment for verbose level
    verbose_level = my_api.setVerbose(args.verbose)
    print("> Verbose level = %s" % verbose_level)
        
    # 1.5 - Init environment for log
    if args.logfile:
        logFile = os.path.join(output_dir, "pge_lake_tile_" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + ".log")
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
    
    try:
        # 2 - Initialization
        myLakeTile = Processing(pixc_file, pixc_vec_river_file, output_dir, IN_shp_option=shp_option)
        my_api.printInfo(timer.info(0))

        # 3 - Run pre-processing
        myLakeTile.run_preprocessing()
        my_api.printInfo(timer.info(0))

        # 4 - Run processing
        myLakeTile.run_processing()
        my_api.printInfo(timer.info(0))

        # 5 - Run post-processing
        myLakeTile.run_postprocessing()
        my_api.printInfo(timer.info(0))

        my_api.printInfo("")
        my_api.printInfo(timer.stop())

    finally:
        # Close logger
        if args.logfile:
            my_api.closeLogger()
        else:
            consoleHandler.close()
            rootLogger.removeHandler(consoleHandler)

    print("")
    print(timer.stop())
    print("===== lakeTileProcessing = END =====")
