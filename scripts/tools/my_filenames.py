# -*- coding: utf-8 -*-
"""
.. module:: my_filenames.py
    :synopsis: Deal with filenames

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2019 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from __future__ import absolute_import, division, print_function, unicode_literals

import os


# Filenames pattern
PRODUCER = "CNES"  # Product generator
RIVER_TILE_CRID = "Dx0000"  # Composite Release IDentifier for RiverTile processing
LAKE_TILE_CRID = "Dx0000"  # Composite Release IDentifier for LakeTile processing
LAKE_SP_CRID = "Dx0000"  # Composite Release IDentifier for LakeSP processing

PIXC_PREFIX = "SWOT_L2_HR_PIXC_"
PIXC_PATTERN_PRINT = PIXC_PREFIX + "<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>.nc"
PIXC_PATTERN_IND = {"cycle": 4, "pass": 5, "tile_ref": 6, "start_date": 7, "stop_date": 8, "crid": 9, "counter": 10}  # Indices when PIXC_PATTERN.split("_"); None if value not in filename

PIXCVEC_RIVER_PREFIX = "SWOT_L2_HR_PIXCVecRiver_"
PIXCVEC_RIVER_SUFFIX = ".nc"
PIXCVEC_RIVER_PATTERN = PIXCVEC_RIVER_PREFIX + "%03d_%03d_%s_%s_%s_%s_%02d" + PIXCVEC_RIVER_SUFFIX  # PIXCVecRiver filename with %03d=cycle number %03d=pass number %s=tile ref %s=swath %s=begin date %s=end date %s=CRID %s=counter
PIXCVEC_RIVER_PATTERN_PRINT = PIXCVEC_RIVER_PREFIX + "<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>.nc"
PIXCVEC_RIVER_PATTERN_IND = {"cycle": 4, "pass": 5, "tile_ref": 6, "start_date": 7, "stop_date": 8, "crid": 9, "counter": 10}  # Indices when PIXCVEC_RIVER_PATTERN.split("_"); None if value not in filename

RIVER_TILE_PREFIX = "SWOT_L2_HR_RiverTile_"
RIVER_TILE_PATTERN = RIVER_TILE_PREFIX + "%03d_%03d_%s_%s_%s_%s_%02d%s"  # RiverTile filename with %03d=cycle number %03d=pass number %s=tile ref %s=swath %s=begin date %s=end date %s=CRID %s=counter %s=suffix
RIVER_TILE_PATTERN_PRINT = RIVER_TILE_PREFIX + "%s<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>"
RIVER_TILE_PATTERN_IND = {"cycle": 4, "pass": 5, "tile_ref": 6, "start_date": 7, "stop_date": 8, "crid": 9, "counter": 10}  # Indices when RIVER_TILE_*_PATTERN.split("_"); None if value not in filename
RIVER_TILE_SUFFIX = ".nc"
RIVER_TILE_NODES_SUFFIX = "_nodes"
RIVER_TILE_REACHES_SUFFIX = "_reaches"

PATTERN_FILE_ANNOT = "river-annotation_%03d_%03d_%s.rdf"  # RiverTile annotation filename with %03d=cycle number %03d=pass number %s=tile ref

LAKE_TILE_PREFIX = "SWOT_L2_HR_LakeTile_"
LAKE_TILE_PATTERN = LAKE_TILE_PREFIX + "%03d_%03d_%s_%s_%s_%s_%02d%s"  # LakeTile filename with %03d=cycle number %03d=pass number %s=tile ref %s=swath %s=begin date %s=end date %s=CRID %s=counter %s=suffix
LAKE_TILE_PATTERN_PRINT = LAKE_TILE_PREFIX + "%s<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>"
LAKE_TILE_PATTERN_IND = {"cycle": 4, "pass": 5, "tile_ref": 6, "start_date": 7, "stop_date": 8, "crid": 9, "counter": 10}  # Indices when LAKE_TILE_*_PATTERN.split("_"); None if value not in filename
LAKE_TILE_SHP_SUFFIX = ".shp"
LAKE_TILE_SHP_META_SUFFIX = ".shp.xml"
LAKE_TILE_EDGE_SUFFIX = "_edge.nc"
LAKE_TILE_PIXCVEC_SUFFIX = "_pixcvec.nc"

LAKE_SP_PREFIX = "SWOT_L2_HR_LakeSP_"
LAKE_SP_PATTERN = LAKE_SP_PREFIX + "%03d_%03d_%s_%s_%s_%s_%02d.shp"  # LakeSP filename with %03d=cycle number %03d=pass number %s=continent %s=begin date %s=end date %s=CRID %s=counter
LAKE_SP_PATTERN_NO_CONT = LAKE_SP_PREFIX + "%03d_%03d_%s_%s_%s_%02d.shp"  # LakeSP filename without continent info with %03d=cycle number %03d=pass number %s=begin date %s=end date %s=CRID %s=counter

PIXCVEC_PREFIX = "SWOT_L2_HR_PIXCVec_"
PIXCVEC_SUFFIX = ".nc"
PIXCVEC_PATTERN = PIXCVEC_PREFIX + "%03d_%03d_%s_%s_%s_%s_%02d" + PIXCVEC_SUFFIX  # PIXCVec filename with %03d=cycle number %03d=pass number %s=tile ref %s=begin date %s=end date %s=CRID %s=counter 


#######################################


def getInfoFromFilename(in_filename, in_type):
    """
    Retrieve orbit info from in_filename

    :param in_filename: input full path
    :type in_filename: string
    :param in_type: type of product =PIXC =PIXCVecRiver =LakeTile
    :type in_type: string

    :return: dictionnary containing values found in in_filename
    :rtype: dict
    """
    
    # 0 - Init variables
    # 0.1 - Output dictionary
    out_dict = {}
    # 0.2 - Get pattern variables depending on in_type
    if in_type == "PIXC":
        if not os.path.basename(in_filename).startswith(PIXC_PREFIX):
            print("Filename %s doesn't match with PIXC pattern %s",
                  os.path.basename(in_filename), PIXC_PATTERN_PRINT)
            out_dict = {}
        pattern = PIXC_PATTERN_IND
    elif in_type == "PIXCVecRiver":
        if not os.path.basename(in_filename).startswith(PIXCVEC_RIVER_PREFIX):
            print("Filename %s doesn't match with PIXCVecRiver pattern %s",
                  os.path.basename(in_filename), PIXCVEC_RIVER_PATTERN_PRINT)
            out_dict = {}
        pattern = PIXCVEC_RIVER_PATTERN_IND
    elif in_type == "LakeTile":
        if not os.path.basename(in_filename).startswith(LAKE_TILE_PREFIX):
            message = "Filename %s doesn't match with LakeTile pattern %s" % (os.path.basename(in_filename), LAKE_TILE_PATTERN_PRINT)
            raise ValueError(message)
        pattern = LAKE_TILE_PATTERN_IND
    else:
        message = "Type %s unknown ; should be PIXC, PIXCVecRiver or LakeTile" % in_type
        raise ValueError(message)

    # 1 - Split basename
    basename_split = os.path.splitext(os.path.basename(in_filename))[0].split("_")

    # 2 - Get values in filename
    # Add error check
    try:
        for key, val in pattern.items():  # Loop on keys
            tmp_val = None  # Init output value
            if val is not None:
                tmp_val = basename_split[val]  # Read value in filename if not None
            out_dict[key] = tmp_val  # Store value in output dictionary
    except:
        message = "Filename %s doesn't match with LakeTile pattern %s" % (os.path.basename(in_filename), LAKE_TILE_PATTERN_PRINT)
        raise ValueError(message)

    # 3 - Check if cycle and pass fields could be convert into int
    try:
        tmp_val = int(out_dict["cycle"])
        tmp_val = int(out_dict["pass"])
    except:
        message = "Filename %s doesn't match with LakeTile pattern %s" % (os.path.basename(in_filename), LAKE_TILE_PATTERN_PRINT)
        raise ValueError(message)

    return out_dict


#######################################


class riverTileFilenames(object):

    def __init__(self, IN_pixc_file=None):
        """
        Constructor of basenames of the output files of l2pixc_to_rivertile.py
        
        :param IN_pixc_file: PIXC filename
        :IN_pixc_file: string
        """
        
        # Init tile variables
        self.cycle_num = None
        self.pass_num = None
        self.tile_ref = None
        self.start_date = None
        self.stop_date = None
        
        # Init filenames
        if IN_pixc_file is not None:
            self.updateWithPixcFilename(IN_pixc_file)
        else:
            self.rivertile_file = None
            self.rivertile_nodes_file = None
            self.rivertile_reaches_file = None
            self.pixc_vec_river_file = None
            self.annot_file = None
    
    #----------------------------------
        
    def updateWithPixcFilename(self, IN_pixc_file):
        """
        Update output filenames from PIXC filename
        
        :param IN_pixc_file: PIXC filename
        :IN_pixc_file: string
        """
        
        # Set self attributes
        self.set_vars_from_filename(IN_pixc_file)
        
        # Set filenames
        self.computeRiverTileFilename()
        self.computeRiverTileNodesFilename()
        self.computeRiverTileReachesFilename()
        self.computePixcVecRiverFilename()
        self.computeAnnotFileFilename()
            
    def set_vars_from_filename(self, IN_pixc_filename):
        """
        Update self variables from PIXC filename
        
        :param IN_pixc_file: PIXC filename
        :IN_pixc_file: string
        """
        
        print()
        
        # 1 - Retrieve info from PIXC filename
        tmp_dict = getInfoFromFilename(IN_pixc_filename, "PIXC")
        
        # 2 - Set self variables
        # 2.1 - Cycle number
        cycle_num = tmp_dict["cycle"]
        if cycle_num is None:
            self.cycle_num = 999
            print("WARNING: cycle number has not been found in PIXC filename %s -> set to default value = %03d",
                  os.path.basename(self.pixc_file), self.cycle_num)
        else:
            self.cycle_num = int(cycle_num)
            print("Cycle number = %03d" % self.cycle_num)
        # 2.2 - Pass number
        pass_num = int(tmp_dict["pass"])
        if pass_num is None:
            self.pass_num = 999
            print("WARNING: pass number has not been found in PIXC filename %s -> set to default value = %03d",
                  os.path.basename(self.pixc_file), self.pass_num)
        else:
            self.pass_num = int(pass_num)
            print("Pass number = %03d" % self.pass_num)
        # 2.3 - Tile ref
        self.tile_ref = tmp_dict["tile_ref"]
        if self.tile_ref is None:
            self.tile_ref = "ttts"
            print("WARNING: tile ref has not been found in PIXC filename %s -> set to default value = %s",
                  os.path.basename(self.pixc_file), self.tile_ref)
        else:
            print("Tile ref = %s" % self.tile_ref)
        # 2.4 - Start date
        self.start_date = tmp_dict["start_date"]
        if self.start_date is None:
            self.start_date = "yyyyMMddThhmmss"
            print("WARNING: start date has not been found in PIXC filename %s -> set to default value = %s",
                  os.path.basename(self.pixc_file), self.start_date)
        else:
            print("Start date = %s" % self.start_date)
        # 2.5 - Stop date
        self.stop_date = tmp_dict["stop_date"]
        if self.stop_date is None:
            self.stop_date = "yyyyMMddThhmmss"
            print("WARNING: stop date has not been found in PIXC filename %s -> set to default value = %s",
                  os.path.basename(self.pixc_file), self.stop_date)
        else:
            print("Stop date = %s" % self.stop_date)
            
        print()
    
    #----------------------------------
    
    def computeRiverTileFilename(self):
        """
        Compute RiverTile NetCDF basename
        """
        self.rivertile_file = RIVER_TILE_PATTERN % (self.cycle_num, self.pass_num, self.tile_ref, self.start_date, self.stop_date, RIVER_TILE_CRID, 1, RIVER_TILE_SUFFIX)
    
    def computeRiverTileNodesFilename(self):
        """
        Compute RiverTile nodes basename
        """
        self.rivertile_nodes_file = RIVER_TILE_PATTERN % (self.cycle_num, self.pass_num, self.tile_ref, self.start_date, self.stop_date, RIVER_TILE_CRID, 1, RIVER_TILE_NODES_SUFFIX)
        
    def computeRiverTileReachesFilename(self):
        """
        Compute RiverTile reaches basename
        """
        self.rivertile_reaches_file = RIVER_TILE_PATTERN % (self.cycle_num, self.pass_num, self.tile_ref, self.start_date, self.stop_date, RIVER_TILE_CRID, 1, RIVER_TILE_REACHES_SUFFIX)
        
    def computePixcVecRiverFilename(self):
        """
        Compute PIXCVecRiver basename
        """
        self.pixc_vec_river_file = PIXCVEC_RIVER_PATTERN % (self.cycle_num, self.pass_num, self.tile_ref, self.start_date, self.stop_date, RIVER_TILE_CRID, 1)
        
    def computeAnnotFileFilename(self):
        """
        Compute file annotation full path
        """
        self.annot_file = PATTERN_FILE_ANNOT % (self.cycle_num, self.pass_num, self.tile_ref)