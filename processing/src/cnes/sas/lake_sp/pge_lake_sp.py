# -*- coding: utf8 -*-
"""
.. module:: pge_lake_sp.py
    :synopsis: Process PGE_L2_HR_LakeSP, i.e. generate L2_HR_LakeSP shp and update L2_HR_PIXCVec NetCDF products from files produced by PGE_L2_HR_LakeTile
    Created on 27/09/2017

.. moduleauthor:: Cécile Cazals - CS

Copyright (c) 2017 CNES. All rights reserved.
"""
from __future__ import absolute_import, division, print_function, unicode_literals 

import argparse
import configparser as cfg
import datetime
from lxml import etree as ET
import numpy as np
import os
import sys

import cnes.sas.lake_sp.proc_pixc_sp as proc_pixc_sp
import cnes.sas.lake_sp.proc_pixc_vec_sp as proc_pixc_vec_sp
import cnes.common.lib.my_api as my_api
import cnes.common.lib.my_shp_file as my_shp
import cnes.common.lib.my_timer as my_timer
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib_lake.lake_db as lake_db
import cnes.common.lib_lake.locnes_filenames as my_names
import cnes.common.lib_lake.locnes_variables as my_var
import cnes.common.lib_lake.proc_lake as proc_lake


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
        my_api.printInfo("[lakeSPProcessing] == INIT ==")

        # Input paths
        self.lake_tile_dir = IN_laketile_dir  # Full path of directory containing LakeTile products
        self.output_dir = IN_output_dir  # Output directory
        
        # Pass infos
        self.cycle_num = IN_cycle_num  # Cycle number
        self.pass_num = IN_pass_num  # Pass number
        self.ascending = True  # Orientation: True=ascending False=descending
        self.asc_dict = {True: "Ascending", False: "Descending"}  # Associated dictionary

        # Flag to produce PIXCVec shapefile
        self.flag_prod_shp = False
        if IN_shp_option:
            self.flag_prod_shp = True

        # Number of tiles to process
        self.nb_input_tiles = 0
        # List of continent(s) processed for the pass
        self.list_continent = []
        # List of all LakeTile shapefiles to merge, organized by continent
        self.lake_tile_shp_file_path_list = {}
        # List of files to process
        self.lake_tile_pixcvec_path_list_R = {}  # List of _pixcvec files for right swath, organized by continent
        self.lake_tile_edge_path_list_R = {}  # List of _edge files for right swath, organized by continent
        self.lake_tile_pixcvec_path_list_L = {}  # List of _pixcvec files for left swath, organized by continent
        self.lake_tile_edge_path_list_L = {}  # List of _edge files for left swath, organized by continent
        
        # LakeSP filenames
        self.lake_sp_filenames = None
        
        # Objects
        self.objLakeDb = None  # Lake DB object
        self.objPixc = None  # PIXC object
        self.objPixcVec = None  # PIXCVecRiver object
        self.objLake = None  # LakeTile object
        
        # SP objects
        # Objects related to edge PixC for right (R) and left (L) swaths
        self.objPixc_SP_L = None
        self.objPixc_SP_R = None
        # Objects related to PIXCVec for right (R) and left (L) swaths
        self.objPixc_Vec_SP_L = None
        self.objPixc_Vec_SP_R = None
        # Objects related to LakeSP for right (R) and left (L) swaths
        self.objLake_SP_L = None
        self.objLake_SP_R = None

    def run_preprocessing(self):
        """
        Retrieve the list of input files, i.e. L2_HR_LakeTile products for wanted pass = shapefile + PIXC_edge + PIXCVec files
        """
        my_api.printInfo("")
        my_api.printInfo("")
        my_api.printInfo("[lakeSPProcessing] PRE-PROCESSING...")
        my_api.printInfo("")

        # 1 - Test existence of directories
        my_api.printInfo("[lakeSPProcessing] > 1 - Testing existence of working directories ...")
        # 1.1 - LakeTile directory
        my_api.printInfo("[lakeSPProcessing]   INPUT LakeTile DIR = %s" % self.lake_tile_dir)
        my_tools.testDir(self.lake_tile_dir)
        # 1.2 - Output directory
        my_api.printInfo("[lakeSPProcessing]   OUTPUT DIR = %s" % self.output_dir)
        my_tools.testDir(self.output_dir)
        my_api.printInfo("")   
        
        # 2 - Get list of input files
        my_api.printInfo("[lakeSPProcessing] > 2 - Retrieving input files ...")
        
        # 2.1 - Get ascending or descending orientation
        self.ascending = (self.pass_num%2 == 0)
        # TODO: replace when new orbits
        # NB: Ascending if pass_num is odd, descending if pass_num is pair.
        # self.ascending = ((self.pass_num%2 - 1) == 0)

        # 2.2 - Compute file prefix regarding cycle and pass conditions
        cond_prefix = my_var.LAKE_TILE_PREFIX
        # Add cycle number condition
        cond_prefix += "%03d" % self.cycle_num
        # Add orbit number condition
        cond_prefix += "_%03d" % self.pass_num
        
        my_api.printInfo("[lakeSPProcessing]   LakeTile files with cycle=%03d and orbit=%03d (%s)" % (self.cycle_num, self.pass_num, self.asc_dict[self.ascending]))

        # 2.3 - List all files in LakeTile directory
        lake_tile_list = os.listdir(self.lake_tile_dir)
            
        # 2.4 - For each listed file, if it's a metadata file related to LakeTile_shp, get related _edge and _pixcvec files if they exist
        
        # Init
        tile_ref_list_R = {}  # List of tile reference for right swath, organized by continent
        tile_ref_list_L = {}  # List of tile reference for left swath, organized by continent
        
        # Loop on LakeTile_pixcvec files
        flag_first = True  # Flag if it's the first tile to deal with
        for curFile in lake_tile_list:
            
            # Test if file meets the condition
            if curFile.startswith(cond_prefix) and curFile.endswith(my_var.LAKE_TILE_SHP_META_SUFFIX):  # Test if it's a wanted LakeTile_shp file (NB: .shp.xml used instead of .shp because _edge and _pixcvec may also have associated .shp file)

                # Init add_tuple flag to True; set to False if 1 file of LakeTile product is missing
                add_tuple = True

                # Shapefile
                cur_shp = curFile.replace(my_var.LAKE_TILE_SHP_META_SUFFIX, my_var.LAKE_TILE_SHP_SUFFIX)
                if not os.path.exists(os.path.join(self.lake_tile_dir, cur_shp)):  # Test if associated _shp file exists
                    add_tuple = False

                # Shapefile
                cur_pixcvec = curFile.replace(my_var.LAKE_TILE_SHP_META_SUFFIX, my_var.LAKE_TILE_PIXCVEC_SUFFIX)
                if not os.path.exists(os.path.join(self.lake_tile_dir, cur_shp)):  # Test if associated _shp file exists
                    add_tuple = False
                    
                # Edge file
                cur_edge = curFile.replace(my_var.LAKE_TILE_SHP_META_SUFFIX, my_var.LAKE_TILE_EDGE_SUFFIX)
                if not os.path.exists(os.path.join(self.lake_tile_dir, cur_edge)):  # Test if associated _edge file exists
                    add_tuple = False
                    
                # Add tuple if exists
                if add_tuple:
                    
                    self.nb_input_tiles += 1
                    
                    # Get metadata
                    metadata = ET.parse(os.path.join(self.lake_tile_dir, curFile))
                    try:
                        cur_continent = metadata.xpath("//LakeTile_shp/tile_info/continent")[0].text
                    except:
                        cur_continent = "WORLD"
                    
                    if flag_first:
                        flag_first = False
                        # Init list of continents for the pass
                        self.list_continent = [cur_continent]
                        # Init lists for continent
                        self.lake_tile_shp_file_path_list[cur_continent] = []
                        self.lake_tile_pixcvec_path_list_R[cur_continent] = []
                        self.lake_tile_edge_path_list_R[cur_continent] = []
                        tile_ref_list_R[cur_continent] = []
                        self.lake_tile_pixcvec_path_list_L[cur_continent] = []
                        self.lake_tile_edge_path_list_L[cur_continent] = []
                        tile_ref_list_L[cur_continent] = []
                        # Overwrite metadata if 1st file processed
                        print()
                        print("WORKING VARIABLES retrieved from LakeTile processing")
                        my_var.overwriteConfig_from_xml(metadata)
                        print()
                        
                    else:
                        # Test if new continent
                        if not cur_continent in self.list_continent:
                            # Add new continent to the list
                            self.list_continent.append(cur_continent)
                            # Init lists for new continent
                            self.lake_tile_shp_file_path_list[cur_continent] = []
                            self.lake_tile_pixcvec_path_list_R[cur_continent] = []
                            self.lake_tile_edge_path_list_R[cur_continent] = []
                            tile_ref_list_R[cur_continent] = []
                            self.lake_tile_pixcvec_path_list_L[cur_continent] = []
                            self.lake_tile_edge_path_list_L[cur_continent] = []
                            tile_ref_list_L[cur_continent] = []
                        # Metadata should be the same as the others
                        my_var.compareConfig_to_xml(metadata)
                    
                    # Add LakeTile_shp to list
                    self.lake_tile_shp_file_path_list[cur_continent].append(os.path.join(self.lake_tile_dir, cur_shp))
                    
                    # Get latitude from filename
                    # TODO: change when tile numbering is fixed
                    TMP_infos = my_names.getInfoFromFilename(curFile, "LakeTile")
                    TMP_tile = TMP_infos["tile_ref"].split("-")[0]
                    TMP_lat = int(TMP_tile[:-1])
                    if TMP_tile.endswith("S"):
                        TMP_lat = -TMP_lat
                        
                    # In Right swath list
                    if "-R" in curFile:
                        self.lake_tile_pixcvec_path_list_R[cur_continent].append(os.path.join(self.lake_tile_dir, cur_pixcvec))
                        self.lake_tile_edge_path_list_R[cur_continent].append(os.path.join(self.lake_tile_dir, cur_edge))
                        tile_ref_list_R[cur_continent].append(TMP_lat)
                        
                    # In Left swath list
                    elif "-L" in curFile:
                        self.lake_tile_pixcvec_path_list_L[cur_continent].append(os.path.join(self.lake_tile_dir, cur_pixcvec))
                        self.lake_tile_edge_path_list_L[cur_continent].append(os.path.join(self.lake_tile_dir, cur_edge))
                        tile_ref_list_L[cur_continent].append(TMP_lat)
                        
        # 2.5 - Test list of continents
        if ("WORLD" in self.list_continent) and (len(self.list_continent) > 1):
            my_api.exitWithError("[lakeSPProcessing] Mix of continent and no continent split; look at tiles process")
        
        # 2.6 - Sort files from south to north, continent per continent
        for curContinent in self.list_continent:
            sorted_idx_R = np.argsort(tile_ref_list_R[curContinent])
            self.lake_tile_pixcvec_path_list_R[curContinent] = [self.lake_tile_pixcvec_path_list_R[curContinent][ind] for ind in sorted_idx_R]
            self.lake_tile_edge_path_list_R[curContinent] = [self.lake_tile_edge_path_list_R[curContinent][ind] for ind in sorted_idx_R]
            sorted_idx_L = np.argsort(tile_ref_list_L[curContinent])
            self.lake_tile_pixcvec_path_list_L[curContinent] = [self.lake_tile_pixcvec_path_list_L[curContinent][ind] for ind in sorted_idx_L]
            self.lake_tile_edge_path_list_L[curContinent] = [self.lake_tile_edge_path_list_L[curContinent][ind] for ind in sorted_idx_L]
            
        # 2.7 - Print list of files, per continent
        for curContinent in self.list_continent:
            my_api.printInfo("[lakeSPProcessing]   > Continent %s --> %d tile(s) to deal with" % (curContinent, len(self.lake_tile_shp_file_path_list[curContinent])))
            for curFile in self.lake_tile_shp_file_path_list[curContinent]:
                my_api.printInfo("[lakeSPProcessing]   %s" % os.path.basename(curFile))
            my_api.printInfo("")
        my_api.printInfo("[lakeSPProcessing]   --> %d tile(s) to deal with, over %d continent(s)" % (self.nb_input_tiles, len(self.list_continent)))
        my_api.printInfo("") 

        # 3 - Retrieve lake Db layer
        my_api.printInfo("[lakeTileProcessing] > 3 - Retrieving lake database layer...")
        if my_var.LAKE_DB == "":
            my_api.printInfo("[lakeTileProcessing] NO database specified -> NO link of SWOT obs with a priori lake")
        else:
            if os.path.exists(my_var.LAKE_DB):
                type_db = my_var.LAKE_DB.split('.')[-1]  # Type of database
                if type_db == "shp":  # Shapefile format
                    self.objLakeDb = lake_db.LakeDb_shp(my_var.LAKE_DB)
                elif type_db == "sqlite":  # SGLite format
                    self.objLakeDb = lake_db.LakeDb_sqlite(my_var.LAKE_DB)
                else:
                    my_api.exitWithError("[lakeTileProcessing] Lake a priori database format (%s) is unknown: must be .shp or .sqlite" % type_db)
            else:
                my_api.exitWithError("[lakeTileProcessing]   ERROR = %s doesn't exist" % my_var.LAKE_DB)
        my_api.printInfo("")

    def run_processing(self):
        """
        Process SAS_L2_HR_LakeTile
        """
        my_api.printInfo("")
        my_api.printInfo("")
        my_api.printInfo("[lakeSPProcessing] PROCESSING tiles of cycle %03d and pass %03d (%s)..." % (self.cycle_num, self.pass_num, self.asc_dict[self.ascending]))
        my_api.printInfo("")
        
        for curContinent in self.list_continent:
        
            my_api.printInfo("")
            my_api.printInfo("==========================================")
            my_api.printInfo("[lakeSPProcessing] Processing continent %s" % curContinent)
            my_api.printInfo("==========================================")
            my_api.printInfo("")
            
            # 1 - Compute output filenames
            my_api.printInfo("[lakeTileProcessing] 1 - Computing LakeSP filenames...")
            self.lake_sp_filenames = my_names.lakeSPFilenames(self.lake_tile_shp_file_path_list[curContinent], curContinent, self.output_dir)
            my_api.printInfo("")
            
            # 2 - Objects initialisation
            my_api.printInfo("[lakeSPProcessing] 2 - Init objects...")
            my_api.printInfo("")
            
            # 2.1 - Right swath
            my_api.printInfo("[lakeSPProcessing] 2a - Right swath")
            # PIXC_edge
            self.objPixc_SP_R = proc_pixc_sp.PixC_Edge_SP(self.ascending, self.lake_tile_edge_path_list_R[curContinent])
            # PIXCVec
            self.objPixc_Vec_SP_R = proc_pixc_vec_sp.PixC_Vec_SP(self.lake_tile_pixcvec_path_list_R[curContinent], self.objPixc_SP_R, self.output_dir, curContinent)
            # LakeSP
            self.objLake_SP_R = proc_lake.LakeProduct("SP", 
                                                      self.objPixc_SP_R, 
                                                      self.objPixc_Vec_SP_R, 
                                                      self.objLakeDb, 
                                                      "tmp_R", 
                                                      IN_id_prefix=self.lake_sp_filenames.lake_id_prefix)
            my_api.printInfo("") 

            # 2.2 - Left swath
            my_api.printInfo("[lakeSPProcessing] 2b - Left swath")
            # PIXC_edge
            self.objPixc_SP_L = proc_pixc_sp.PixC_Edge_SP(self.ascending, self.lake_tile_edge_path_list_L[curContinent])
            # PIXCVec
            self.objPixc_Vec_SP_L = proc_pixc_vec_sp.PixC_Vec_SP(self.lake_tile_pixcvec_path_list_L[curContinent], self.objPixc_SP_L, self.output_dir, curContinent)
            # LakeSP
            self.objLake_SP_L = proc_lake.LakeProduct("SP", 
                                                      self.objPixc_SP_L, 
                                                      self.objPixc_Vec_SP_L, 
                                                      self.objLakeDb, 
                                                      "tmp_L", 
                                                      IN_id_prefix=self.lake_sp_filenames.lake_id_prefix)
            my_api.printInfo("")
            
            # 3 - Compute lake products
            my_api.printInfo("[lakeSPProcessing] 3 - Computing lake products...")
            
            # 3.1 - Processing swath R if there are pixels to process
            my_api.printInfo("")
            my_api.printInfo(">>> RIGHT SWATH <<<")
            
            if self.objPixc_SP_R.nb_pixels > 0:
                
                my_api.printInfo("")
                
                # 3.1.1 - Gather pixels by entities for all the tiles of this swath
                self.objPixc_SP_R.edgeGlobalRelabeling() 
                
                # 3.1.2 - Compute lake products
                self.objLake_SP_R.computeLakeProducts(np.unique(self.objPixc_SP_R.labels))
                
            else:
                
                my_api.printInfo("No pixel to process")
            
            # 3.2 - Processing swath L if there are pixels to process
            my_api.printInfo("")
            my_api.printInfo(">>> LEFT SWATH <<<")
            
            if self.objPixc_SP_L.nb_pixels > 0:
                
                my_api.printInfo("")
                
                # 3.2.1 - Gather pixels by entities for all the tiles of this swath
                self.objPixc_SP_L.edgeGlobalRelabeling() 
                
                # 3.2.2 - Compute lake products
                self.objLake_SP_L.computeLakeProducts(np.unique(self.objPixc_SP_L.labels))
                
            else:
                
                my_api.printInfo("No pixel to process")
                
            my_api.printInfo("")
            
            # 4 - Merge shapefiles to get LakeSP product
            my_api.printInfo("[lakeSPProcessing] 4 - Merging shapefiles to get LakeSP product %s..." % os.path.basename(self.lake_sp_filenames.lake_sp_file))
            # 4.1 - Merging Right and Left SP layers
            my_api.printDebug("[lakeSPProcessing] > Merging right and left SP layers...")
            dataSource_sp1, layer_sp = my_shp.merge2Layers(self.objLake_SP_R.layer, self.objLake_SP_L.layer)
            # 4.2 - Merge SP layer with shapefiles retrieved from PGE_LakeTile
            my_api.printDebug("[lakeSPProcessing] > Merging SP layer with shapefiles retrieved from PGE_LakeTile...")
            dataSource_sp2, layer_sp = my_shp.mergeMemLayerWithShp(self.lake_tile_shp_file_path_list[curContinent], layer_sp)
            # 4.3 - Write LakeSP shapefile product
            my_api.printDebug("[lakeSPProcessing] > Writing L2_HR_LakeSP shapefile = %s" % os.path.basename(self.lake_sp_filenames.lake_sp_file))
            my_shp.writeMemLayer_asShp(layer_sp, self.lake_sp_filenames.lake_sp_file)
            # 4.4 - Write XML metadatafile for shapefile
            if self.objPixc_SP_R.pass_num != 0:
                self.objLake_SP_R.writeMetadataFile("%s.xml" % self.lake_sp_filenames.lake_sp_file)
            elif self.objPixc_SP_L.pass_num != 0:
                self.objLake_SP_L.writeMetadataFile("%s.xml" % self.lake_sp_filenames.lake_sp_file)
            # 4.5 - Close dataSources
            dataSource_sp1.Destroy()
            dataSource_sp2.Destroy()
            self.objLake_SP_R.dataSource.Destroy()
            self.objLake_SP_L.dataSource.Destroy()
            my_api.printInfo("")

            # 5 - Update PIXCVec files
            my_api.printInfo("[lakeSPProcessing] 5 - Updating L2_HR_PIXCVec files...")
            self.objPixc_Vec_SP_R.updatePixcVec(self.flag_prod_shp)
            self.objPixc_Vec_SP_L.updatePixcVec(self.flag_prod_shp)
            my_api.printInfo("")
            
        my_api.printInfo("")

    def run_postprocessing(self):
        """
        Process PGE_L2_HR_LakeTile, i.e. convert output data in L2_HR_LakeTile product
        """
        my_api.printInfo("")
        my_api.printInfo("")
        my_api.printInfo("[lakeSPProcessing] POST-PROCESSING...")
        my_api.printInfo("")

        # 1 - Close lake database
        if self.objLakeDb is not None:
            my_api.printInfo("[lakeTileProcessing] 1 - Closing lake database...")
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
    print("[lakeSPProcessing] == readParamFile = %s ==" % IN_filename)
    
    # 0 - Init output dictionary
    OUT_params = {}
    # Default values
    OUT_params["flag_prod_shp"] = False
    
    # 1 - Read parameter file
    config = cfg.ConfigParser()
    config.read(IN_filename)
    
    # 2 - Retrieve PATHS
    OUT_params["laketile_dir"] = config.get("PATHS", "LakeTile directory")
    OUT_params["output_dir"] = config.get("PATHS", "Output directory")
    
    #• 3 - Retrieve TILES_INFOS
    OUT_params["cycle_num"] = config.getint("TILES_INFOS", "Cycle number")
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
    parser = argparse.ArgumentParser(description="Compute SWOT LakeSP product from LakeTile products corresponding to a specific (cycle, pass). \
                                     If indir_or_param_file is a parameter file (*.cfg), all input parameters are only read in the parameter file.")
    parser.add_argument("indir_or_param_file", help="LakeTile directory or parameter file (*.cfg)")
    parser.add_argument("output_dir", help="output directory", nargs='?')
    parser.add_argument("cycle_num", help="cycle number", type=int, nargs='?')
    parser.add_argument("pass_num", help="pass number", type=int, nargs='?')
    parser.add_argument("-shp", help="convert output NetCDF file as shapefile", action="store_true")
    parser.add_argument("-l", "--logfile", help="write prints to a logfile", action="store_true")  # To print logs on screen (=False, default) or in a logfile (=True)
    parser.add_argument("-v", "--verbose", help="verbose level", choices=["DEBUG", "INFO"], default="INFO")  # Verbose level
    args = parser.parse_args()
    
    print("===== lakeSPProcessing = BEGIN =====")
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
            print("Run by pge_lake_sp.py param_file.cfg [-l] [-v VERBOSE]")
            print("OR pge_lake_sp.py lake_tile_dir output_dir cycle_num pass_num [-shp] [-l] [-v VERBOSE]")
            sys.exit("indir_or_param_file is %s, not .cfg" % file_extent)
    
    # 1.3 - Test input params have been filled
    # 1.3.1 - LakeTile directory
    if laketile_dir is None:
        my_api.exitWithError("LakeTile directory is missing in %s" % location)
    # 1.3.2 - Output directory
    if output_dir is None:
        my_api.exitWithError("Output directory is missing in %s" % location)
    my_tools.testDir(output_dir)  # Test existence of output directory here and not in pre-proc because used in 1.5
    # 1.3.3 - Cycle number
    if cycle_num is None:
        my_api.exitWithError("Cycle number is missing in %s" % location)
    # 1.3.4 - Pass number
    if pass_num is None:
        my_api.exitWithError("Pass number is missing in %s" % location)
    
    # 1.4 - Init environment for verbose level
    verbose_level = my_api.setVerbose(args.verbose)
    print("> Verbose level = %s" % verbose_level)
        
    # 1.5 - Init environment for log
    if args.logfile:
        logFile = os.path.join(output_dir, "pge_lake_sp_" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + ".log")
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

    # 5 - Run post-processing
    myLakeSP.run_postprocessing()
    my_api.printInfo(timer.info(0))

    my_api.printInfo("")
    my_api.printInfo(timer.stop())
    
    # Close logger
    if args.logfile:
        my_api.closeLogger()

    print("")
    print((timer.stop()))
    print("===== lakeSPProcessing = END =====")
