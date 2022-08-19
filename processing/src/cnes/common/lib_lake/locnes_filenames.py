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
# VERSION:4.0.0:DM:#91:2022/05/05:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: locnes_filenames.py
    :synopsis: Deal with filenames used by LOCNES. There is one specific class per processor.

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import os

import cnes.common.service_error as service_error
import cnes.common.service_config_file as service_config_file

import cnes.common.lib.my_tools as my_tools
import cnes.common.lib_lake.lake_db as lake_db


#####################
# Filenames pattern #
#####################

# 1 - PIXC
PIXC_PREFIX = "SWOT_L2_HR_PIXC_"
PIXC_PATTERN_PRINT = PIXC_PREFIX + "<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>.nc"
# Indices when PIXC_PATTERN.split("_"); None if value not in filename
PIXC_PATTERN_IND = {"cycle": 4, "pass": 5, "tile_ref": 6, "start_date": 7, "stop_date": 8, "crid": 9, "counter": 10}
PIXC_SUFFIX = ".nc"

# 2 - PIXCVecRiver 
PIXCVEC_RIVER_PREFIX = "SWOT_L2_HR_PIXCVecRiver_"
PIXCVEC_RIVER_PATTERN_PRINT = PIXCVEC_RIVER_PREFIX \
                              + "<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>.nc"
# Indices when PIXCVEC_RIVER_PATTERN.split("_"); None if value not in filename
PIXCVEC_RIVER_PATTERN_IND = {"cycle": 4, "pass": 5, "tile_ref": 6, "start_date": 7, "stop_date": 8, "crid": 9, "counter": 10}

# 3 - LakeTile
LAKE_TILE_PREFIX_BASE = "SWOT_L2_HR_LakeTile_"
# Distinguish each file
LAKE_TILE_PREFIX = {}
LAKE_TILE_PREFIX["obs"] = LAKE_TILE_PREFIX_BASE + "Obs_"
LAKE_TILE_PREFIX["prior"] = LAKE_TILE_PREFIX_BASE + "Prior_"
LAKE_TILE_PREFIX["unknown"] = LAKE_TILE_PREFIX_BASE + "Unassigned_"
LAKE_TILE_PREFIX["edge"] = LAKE_TILE_PREFIX_BASE + "Edge_"
LAKE_TILE_PREFIX["pixcvec"] = LAKE_TILE_PREFIX_BASE + "PIXCVec_"
# LakeTile filename with %03d=cycle number %03d=pass number %s=tile ref %s=swath %s=begin date %s=end date %s=CRID %s=counter %s=suffix
LAKE_TILE_PATTERN_BASE = LAKE_TILE_PREFIX_BASE + "%03d_%03d_%s_%s_%s_%s_%02d%s"
LAKE_TILE_PATTERN = {}
LAKE_TILE_PATTERN_PRINT = {}
for key, value in LAKE_TILE_PREFIX.items():
    LAKE_TILE_PATTERN[key] = value + "%03d_%03d_%s_%s_%s_%s_%02d%s"
    LAKE_TILE_PATTERN_PRINT[key] = value + "%s<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>"
# Indices when LAKE_TILE_*_PATTERN.split("_"); None if value not in filename
LAKE_TILE_PATTERN_IND = {"cycle": 5, "pass": 6, "tile_ref": 7, "start_date": 8, "stop_date": 9, "crid": 10, "counter": 11}
LAKE_TILE_SHP_SUFFIX = ".shp"
LAKE_TILE_SHP_META_SUFFIX = ".shp.xml"
LAKE_TILE_EDGE_SUFFIX = ".nc"
LAKE_TILE_PIXCVEC_SUFFIX = ".nc"

# 4 - LakeSP
LAKE_SP_PREFIX_BASE = "SWOT_L2_HR_LakeSP_"
# Distinguish each file
LAKE_SP_PREFIX = {}
LAKE_SP_PREFIX["obs"] = LAKE_SP_PREFIX_BASE + "Obs_"
LAKE_SP_PREFIX["prior"] = LAKE_SP_PREFIX_BASE + "Prior_"
LAKE_SP_PREFIX["unknown"] = LAKE_SP_PREFIX_BASE + "Unassigned_"
LAKE_SP_SHP_SUFFIX = ".shp"
# LakeSP filename with %03d=cycle number %03d=pass number %s=continent %s=begin date %s=end date %s=CRID %s=counter
LAKE_TILE_PATTERN_BASE = LAKE_TILE_PREFIX_BASE + "%03d_%03d_%s_%s_%s_%s_%02d" + LAKE_SP_SHP_SUFFIX
LAKE_SP_PATTERN = {}
LAKE_SP_PATTERN_PRINT = {}
LAKE_SP_PATTERN_NO_CONT = {}
for key, value in LAKE_SP_PREFIX.items():
    # LakeSP filename with %03d=cycle number %03d=pass number %s=continent %s=begin date %s=end date %s=CRID %s=counter
    LAKE_SP_PATTERN[key] = value + "%03d_%03d_%s_%s_%s_%s_%02d" + LAKE_SP_SHP_SUFFIX
    LAKE_SP_PATTERN_PRINT[key] = value + "%s<CycleID>_<PassID>_<ContinentID>_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>"
    # LakeSP filename without continent info with %03d=cycle number %03d=pass number %s=begin date %s=end date %s=CRID %s=counter
    LAKE_SP_PATTERN_NO_CONT[key] = value + "%03d_%03d_%s_%s_%s_%02d" + LAKE_SP_SHP_SUFFIX
# Indices when LAKE_TILE_*_PATTERN.split("_"); None if value not in filename
LAKE_SP_PATTERN_IND = {"cycle": 5, "pass": 6, "continent": 7, "start_date": 8, "stop_date": 9, "crid": 10, "counter": 11}

# 5 - PIXCVec
PIXCVEC_PREFIX = "SWOT_L2_HR_PIXCVec_"
PIXCVEC_SUFFIX = ".nc"
# PIXCVec filename with %03d=cycle number %03d=pass number %s=tile ref %s=begin date %s=end date %s=CRID %s=counter
PIXCVEC_PATTERN = PIXCVEC_PREFIX + "%03d_%03d_%s_%s_%s_%s_%02d" + PIXCVEC_SUFFIX

# 6 - LakeAvg
LAKE_AVG_PREFIX = "SWOT_L2_HR_LakeAvg_"
LAKE_AVG_SUFFIX = ".shp"
# LakeAvg filename with %03d=cycle number %s=continent %s=basin idenfier %s=begin date %s=end date %s=CRID %s=counter
LAKE_AVG_PATTERN = LAKE_AVG_PREFIX + "%03d_%s_%s_%s_%s_%s_%02d" + LAKE_AVG_SUFFIX

# 7 - Operationnal Prior Lake Database (PLD)
PLD_PREFIX_BASE = "SWOT_LakeDatabase_"
PLD_SUFFIX = ".sqlite"
# Distinguish each file
PLD_PREFIX = {}
PLD_PREFIX["calval"] = PLD_PREFIX_BASE + "Cal_"
PLD_PREFIX["science"] = PLD_PREFIX_BASE + "Nom_"
PLD_PREFIX["basin"] = PLD_PREFIX_BASE
# PLD filenames
PLD_PATTERN_START = {}
for key, value in PLD_PREFIX.items():
    if key == "basin":
        # PLD filename with %s=2-letter continent identifier %02d=2-digit basin code %s=begin date %s=end date %file creation date 
        # %03d=version number (Vvv)
        PLD_PATTERN_START[key] = value + "%s_%s"
    else:
        # PLD filename with %03d=pass number %s=begin date %s=end date %file creation date %03d=version number (Vvv)
        PLD_PATTERN_START[key] = value + "%03d"


#######################################


def get_info_from_filename(in_filename, in_type):
    """
    Retrieve orbit info from in_filename

    :param in_filename: input full path
    :type in_filename: string
    :param in_type: type of product =PIXC =PIXCVecRiver =LakeTile
    :type in_type: string

    :return: dictionnary containing values found in in_filename
    :rtype: dict
    """
    logger = logging.getLogger("locnes_filename")
    
    # 0 - Init variables
    # 0.1 - Output dictionary
    out_dict = {}
    
    # 1 - Get pattern variables depending on in_type
    if in_type == "PIXC":
        if not os.path.basename(in_filename).startswith(PIXC_PREFIX):
            logger.debug("Filename %s doesn't match with PIXC pattern %s",
                        os.path.basename(in_filename), PIXC_PATTERN_PRINT)
            pattern = None
        else:
            pattern = PIXC_PATTERN_IND
    elif in_type == "PIXCVecRiver":
        if not os.path.basename(in_filename).startswith(PIXCVEC_RIVER_PREFIX):
            logger.debug("Filename %s doesn't match with PIXCVecRiver pattern %s",
                        os.path.basename(in_filename), PIXCVEC_RIVER_PATTERN_PRINT)
            pattern = None
        else:
            pattern = PIXCVEC_RIVER_PATTERN_IND
    elif in_type == "LakeTile":
        flag_ok = False
        for prefix in LAKE_TILE_PREFIX.values():
            if os.path.basename(in_filename).startswith(prefix):
                flag_ok = True
        if not flag_ok:
            message = "Filename %s doesn't match with any LakeTile pattern %s" % (os.path.basename(in_filename), \
                                                                                  ", ".join(LAKE_TILE_PATTERN_PRINT.values()))
            raise service_error.SASLakeTileError(message, logger)
        pattern = LAKE_TILE_PATTERN_IND
    elif in_type == "LakeSP":
        flag_ok = False
        for prefix in LAKE_SP_PREFIX.values():
            if os.path.basename(in_filename).startswith(prefix):
                flag_ok = True
        if not flag_ok:
            message = "Filename %s doesn't match with any LakeSP pattern %s" % (os.path.basename(in_filename), \
                                                                                  ", ".join(LAKE_SP_PATTERN_PRINT.values()))
            raise service_error.SASLakeTileError(message, logger)
        pattern = LAKE_SP_PATTERN_IND
    else:
        message = "Type %s unknown ; should be PIXC, PIXCVecRiver or LakeTile" % in_type
        raise service_error.SASLakeTileError(message, logger)
        
    if pattern is not None:

        # 1 - Split basename
        basename_split = os.path.splitext(os.path.basename(in_filename))[0].split("_")

        # 2 - Get values in filename
        try:
            for key, val in pattern.items():  # Loop on keys
                tmp_val = None  # Init output value
                if val is not None:
                    tmp_val = basename_split[val]  # Read value in filename if not None
                out_dict[key] = tmp_val  # Store value in output dictionary
        except:
            message = "Filename %s doesn't match with %s pattern %s" % (os.path.basename(in_filename), in_type, pattern)
            raise service_error.ProcessingError(message, logger)

        # 3 - Check if cycle and pass fields can be convert into int
        # Cycle value
        try:
            tmp_val = int(out_dict["cycle"])
        except:
            message = "Cycle number cannot be converted to integer: %s" % out_dict["cycle"]
            raise service_error.ProcessingError(message, logger)
        # Pass value
        try:
            tmp_val = int(out_dict["pass"])
        except:
            message = "Pass number cannot be converted to integer: %s" % out_dict["pass"]
            raise service_error.ProcessingError(message, logger)

    return out_dict


#######################################


class LakeTileFilenames(object):
    """
    class LakeTileFilenames
    Manage LakeTile product filenames
    """
    
    def __init__(self, in_pixc_file, in_pixcvec_river_file, 
                 in_out_dir, flag_inc=True):
        """
        Constructor of LakeTile filenames

        :param in_pixc_file: PIXC file full path
        :type in_pixc_file: string
        :param in_pixcvec_river_file: PIXCVecRiver file full path
        :type in_pixcvec_river_file: string
        :param in_out_dir: output directory
        :type in_out_dir: string
        :param flag_inc: flag to increment output file counter if True (default); else=False
        :type flag_inc: boolean

        Variables of the object:
            - pixcvec_river_file / string: PIXCVecRiver full path
            - out_dir / string: output directory
            - flag_inc_file_counter / boolean: True to increment output file counter
            - product_counter / int: product counter
            - crid_laketile / string: CRID
            - cycle_num / int: cycle number
            - pass_num / int: pass number
            - tile_ref / str: tile reference, i.e. <tile_number>[R/L]
            - start_date / str: begin of tile time coverage, as AAAAMMDDThhmmss
            - stop_date / str: end of tile time coverage, as AAAAMMDDThhmmss
            - lake_tile_shp_file_obs / str: LakeTile_shp full path of file dedicated to obs-oriented product
            - lake_tile_shp_file_prior / str: LakeTile_shp full path of file dedicated to PLD-oriented product
            - lake_tile_shp_file_unknown / str: LakeTile_shp full path of file of unassigned water features
            - lake_tile_edge_file / str: LakeTile_Edge full path
            - lake_tile_pixcvec_file / str: LakeTile_PIXCVec full path
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        cfg = service_config_file.get_instance()  # Get config file

        # 1 - Init variables
        pixc_file = in_pixc_file  # PIXC file full path
        self.pixcvec_river_file = in_pixcvec_river_file  # PIXCVecRiver full path
        self.out_dir = in_out_dir  # Output directory
        self.flag_inc_file_counter = flag_inc  # Flag to increment file counter
        self.product_counter = 1  # Product counter
        # Get CRID from configuration file
        try:
            self.crid_laketile = cfg.get("FILE_INFORMATION", "CRID")
        except: 
            self.crid_laketile = cfg.get("FILE_INFORMATION", "CRID_LAKETILE")

        # 2 - Retrieve info from PIXC filename
        tmp_dict = get_info_from_filename(pixc_file, "PIXC")
        # 2.1 - Cycle number
        cycle_num = tmp_dict["cycle"]
        if cycle_num is None:
            self.cycle_num = 999
            logger.debug("WARNING: cycle number has not been found in PIXC filename %s -> set to default value = %03d",
                        os.path.basename(pixc_file), self.cycle_num)
        else:
            self.cycle_num = int(cycle_num)
            logger.debug("Cycle number = %03d", self.cycle_num)
        # 2.2 - Pass number
        pass_num = tmp_dict["pass"]
        if pass_num is None:
            self.pass_num = 999
            logger.debug("WARNING: pass number has not been found in PIXC filename %s -> set to default value = %03d",
                        os.path.basename(pixc_file), self.pass_num)
        else:
            self.pass_num = int(pass_num)
            logger.debug("Pass number = %03d", self.pass_num)
        # 2.3 - Tile ref
        self.tile_ref = tmp_dict["tile_ref"]
        if self.tile_ref is None:
            self.tile_ref = "ttts"
            logger.debug("WARNING: tile ref has not been found in PIXC filename %s -> set to default value = %s",
                        os.path.basename(pixc_file), self.tile_ref)
        else:
            logger.debug("Tile ref = %s", self.tile_ref)
        # 2.4 - Start date
        self.start_date = tmp_dict["start_date"]
        if self.start_date is None:
            self.start_date = "yyyyMMddThhmmss"
            logger.debug("WARNING: start date has not been found in PIXC filename %s -> set to default value = %s",
                        os.path.basename(pixc_file), self.start_date)
        else:
            logger.debug("Start date = %s", self.start_date)
        # 2.5 - Stop date
        self.stop_date = tmp_dict["stop_date"]
        if self.stop_date is None:
            self.stop_date = "yyyyMMddThhmmss"
            logger.debug("WARNING: stop date has not been found in PIXC filename %s -> set to default value = %s",
                        os.path.basename(pixc_file), self.stop_date)
        else:
            logger.debug("Stop date = %s", self.stop_date)

        # 3 - Test compatibility of PIXCVecRiver filename with PIXC filename
        tmp_ok = self.test_pixcvec_river_filename()
        if tmp_ok:
            logger.debug("PIXCVecRiver basename %s is compatible with PIXC basename %s",
                        os.path.basename(self.pixcvec_river_file), os.path.basename(pixc_file))
        else:
            logger.debug("WARNING: PIXCVecRiver basename %s IS NOT compatible with PIXC basename %s (cf. above)",
                        os.path.basename(self.pixcvec_river_file), os.path.basename(pixc_file))

        # 4 - Compute output filenames
        self.compute_lake_tile_filename_shp()  # LakeTile_shp filenames
        self.compute_lake_tile_filename_edge()  # LakeTile_Edge filename
        self.compute_lake_tile_filename_pixcvec()  # LakeTile_PIXCVec filename

    #----------------------------------

    def compute_lake_tile_filename_shp(self):
        """
        Compute LakeTile_shp full path of:
        - shapefile dedicated to obs-oriented product
        - shapefile dedicated to PLD-oriented product
        - shapefile of unassigned water features
        """
        
        # 1 - Full path of obs-oriented product
        filename = LAKE_TILE_PATTERN["obs"] % (self.cycle_num, self.pass_num, self.tile_ref, self.start_date, self.stop_date, \
                                             self.crid_laketile, self.product_counter, LAKE_TILE_SHP_SUFFIX)
        self.lake_tile_shp_file_obs = os.path.join(self.out_dir, filename)

        # Test existence and modify self.product_counter if wanted and needed
        if self.flag_inc_file_counter:
            while os.path.exists(self.lake_tile_shp_file_obs):
                self.product_counter += 1
                filename = LAKE_TILE_PATTERN["obs"] % (self.cycle_num, self.pass_num, self.tile_ref, self.start_date, self.stop_date, \
                                                        self.crid_laketile, self.product_counter, LAKE_TILE_SHP_SUFFIX)
                self.lake_tile_shp_file_obs = os.path.join(self.out_dir, filename)
            
        # 2 - Full path of shapefile of PLD-oriented product
        self.lake_tile_shp_file_prior = self.lake_tile_shp_file_obs.replace(LAKE_TILE_PREFIX["obs"], LAKE_TILE_PREFIX["prior"])
        
        # 3 - Full path of shapefile of unassigned water features
        self.lake_tile_shp_file_unknown = self.lake_tile_shp_file_obs.replace(LAKE_TILE_PREFIX["obs"], LAKE_TILE_PREFIX["unknown"])

    def compute_lake_tile_filename_edge(self):
        """
        Compute LakeTile_Edge full path
        """
        
        filename = LAKE_TILE_PATTERN["edge"] % (self.cycle_num, self.pass_num, self.tile_ref, self.start_date, self.stop_date, \
                                            self.crid_laketile, self.product_counter, LAKE_TILE_EDGE_SUFFIX)
        self.lake_tile_edge_file = os.path.join(self.out_dir, filename)

    def compute_lake_tile_filename_pixcvec(self):
        """
        Compute LakeTile_PIXCVec full path
        """
        
        filename = LAKE_TILE_PATTERN["pixcvec"] % (self.cycle_num, self.pass_num, self.tile_ref, self.start_date, self.stop_date, \
                                                self.crid_laketile, self.product_counter, LAKE_TILE_PIXCVEC_SUFFIX)
        self.lake_tile_pixcvec_file = os.path.join(self.out_dir, filename)

    #----------------------------------

    def test_pixcvec_river_filename(self):
        """
        Test if PIXCVecRiver filename is coherent with PIXC file, wrt (cycle, pass, tile) triplet

        :return: True if it's the case, or False if not
        :rtype: boolean
        """
        logger = logging.getLogger(self.__class__.__name__)

        # 0 - Init output boolean
        out_ok = True

        # 1 - Retrieve info from PIXCVecRiver filename
        tmp_dict = get_info_from_filename(self.pixcvec_river_file, "PIXCVecRiver")
        if len(tmp_dict) == 0:
            out_ok = False
        else:
            pixcvec_cycle_num = int(tmp_dict["cycle"])  # Cycle number
            pixcvec_pass_num = int(tmp_dict["pass"])  # Pass number
            pixcvec_tile_ref = tmp_dict["tile_ref"]  # Tile ref

            # 2 - Test each value
            # 2.1 - Cycle number
            if pixcvec_cycle_num != self.cycle_num:
                out_ok = False
                logger.debug("WARNING: cycle number is not the same in PIXC (=%03d) and PIXCVecRiver (=%03d) filenames", \
                            self.cycle_num, pixcvec_cycle_num)
            # 2.2 - Pass number
            if pixcvec_pass_num != self.pass_num:
                out_ok = False
                logger.debug("WARNING: pass number is not the same in PIXC (=%03d) and PIXCVecRiver (=%03d) filenames", \
                            self.pass_num, pixcvec_pass_num)
            # 2.3 - Tile ref
            if pixcvec_tile_ref != self.tile_ref:
                out_ok = False
                logger.debug("WARNING: tile ref not the same in PIXC (=%s) and PIXCVecRiver (=%s) filenames", self.tile_ref, \
                            pixcvec_tile_ref)

        return out_ok


#######################################


class LakeSPFilenames(object):
    """
    class LakeSPFilenames
    Manage LakeSP product filenames
    """
    
    def __init__(self, in_lake_tile_file_list, in_continent_id, 
                 in_out_dir, flag_inc=True):
        """
        Constructor of LakeSP filenames

        :param in_lake_tile_file_list: list of LakeTile_shp files full path
        :type in_lake_tile_file_list: list of string
        :param in_continent_id: 2-letter continent identifier
        :type in_continent_id: string
        :param in_out_dir: output directory
        :type in_out_dir: string
        :param flag_inc: flag to increment output file counter if True (default); else=False
        :type flag_inc: boolean

        Variables of the object:
            - lake_tile_file_list / string: list of LakeTile_shp files full path
            - out_dir / string: output directory
            - flag_inc_file_counter / boolean: True to increment output file counter
            - product_counter / int: product counter
            - crid_lakesp / string: CRID
            - cycle_num / int: cycle number
            - pass_num / int: pass number
            - continent / str: continent identifier
            - start_date / str: begin of continent/pass time coverage, as AAAAMMDDThhmmss
            - stop_date / str: end of continent/pass time coverage, as AAAAMMDDThhmmss
            - lake_sp_file_obs / str: full path of file dedicated to obs-oriented product
            - lake_sp_file_prior / str: full path of file dedicated to PLD-oriented product
            - lake_sp_file_unknown / str: full path of file of unassigned water features
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        # Get config file
        cfg = service_config_file.get_instance()

        # 1 - Init variables
        self.lake_tile_file_list = in_lake_tile_file_list  # List of LakeTile_shp files full path
        self.out_dir = in_out_dir  # Output directory
        self.flag_inc_file_counter = flag_inc  # Flag to increment file counter
        self.product_counter = 1  # Product counter
        # Get CRID from configuration file
        self.crid_lakesp = cfg.get("FILE_INFORMATION", "CRID_LAKESP")

        # 2 - Retrieve info from LakeTile filename
        tmp_dict = get_info_from_filename(self.lake_tile_file_list[0], "LakeTile")
        # 2.1 - Cycle number
        cycle_num = tmp_dict["cycle"]
        if cycle_num is None:
            self.cycle_num = 999
            logger.debug("WARNING: cycle number has not been found in LakeTile filename %s -> set to default value = %03d",
                        os.path.basename(self.lake_tile_file_list[0]), self.cycle_num)
        else:
            self.cycle_num = int(cycle_num)
            logger.debug("Cycle number = %03d", self.cycle_num)
        # 2.2 - Pass number
        pass_num = int(tmp_dict["pass"])
        if pass_num is None:
            self.pass_num = 999
            logger.debug("WARNING: pass number has not been found in LakeTile filename %s -> set to default value = %03d",
                        os.path.basename(self.lake_tile_file_list[0]), self.pass_num)
        else:
            self.pass_num = int(pass_num)
            logger.debug("Pass number = %03d", self.pass_num)
        # 2.3 - Continent identifier
        if in_continent_id == "":
            self.continent_id = ""
            logger.debug("WARNING: no continental split")
        else:
            self.continent_id = in_continent_id
            logger.debug("Continent identifier = %s", self.continent_id)

        # 3 - Retrieve start and stop dates from LakeTile_shp filenames
        self.compute_start_stop_dates()

        # 4 - Compute output filenames
        self.compute_lake_sp_filenames()  # LakeSP filenames

    #----------------------------------

    def compute_start_stop_dates(self):
        """
        Compute start and stop dates from a list of LakeTile_shp filenames
        """

        # Init with 1st filename
        tmp_dict = get_info_from_filename(self.lake_tile_file_list[0], "LakeTile")
        list_start = []
        list_start.append(tmp_dict["start_date"])
        list_stop = []
        list_stop.append(tmp_dict["stop_date"])

        # Process other files
        for cur_file in self.lake_tile_file_list[1:]:
            tmp_dict = get_info_from_filename(cur_file, "LakeTile")
            list_start.append(tmp_dict["start_date"])
            list_stop.append(tmp_dict["stop_date"])

        # Sort dates
        list_start.sort()
        list_stop.sort()

        # Retrieve first and last dates
        self.start_date = list_start[0]
        self.stop_date = list_stop[-1]

    #----------------------------------

    def compute_lake_sp_filenames(self):
        """
        Compute LakeSP full path of:
        - shapefile dedicated to obs-oriented product
        - shapefile dedicated to PLD-oriented product
        - shapefile of unassigned water features
        """

        # 1 - Full path of obs-oriented product
        if self.continent_id is "":
            filename = LAKE_SP_PATTERN_NO_CONT["obs"] % (self.cycle_num, self.pass_num, self.start_date, self.stop_date, self.crid_lakesp, \
                                                          self.product_counter)
        else:
            filename = LAKE_SP_PATTERN["obs"] % (self.cycle_num, self.pass_num, self.continent_id, self.start_date, self.stop_date, \
                                                 self.crid_lakesp, self.product_counter)
        self.lake_sp_file_obs = os.path.join(self.out_dir, filename)

        # Test existence and modify self.product_counter if wanted and needed
        if self.flag_inc_file_counter:
            while os.path.exists(self.lake_sp_file_obs):
                self.product_counter += 1
                if self.continent_id is "":
                    filename = LAKE_SP_PATTERN_NO_CONT["obs"] % (self.cycle_num, self.pass_num, self.start_date, self.stop_date, self.crid_lakesp, \
                                                                  self.product_counter)
                else:
                    filename = LAKE_SP_PATTERN["obs"] % (self.cycle_num, self.pass_num, self.continent_id, self.start_date, self.stop_date, \
                                                         self.crid_lakesp, self.product_counter)
                self.lake_sp_file_obs = os.path.join(self.out_dir, filename)
            
        # 2 - Full path of shapefile of PLD-oriented product
        self.lake_sp_file_prior = self.lake_sp_file_obs.replace(LAKE_SP_PREFIX["obs"], LAKE_SP_PREFIX["prior"])
            
        # 3 - Full path of shapefile of unassigned water features
        self.lake_sp_file_unknown = self.lake_sp_file_obs.replace(LAKE_SP_PREFIX["obs"], LAKE_SP_PREFIX["unknown"])


#######################################


def compute_pixcvec_filename(in_laketile_pixcvec_filename, in_output_dir):
    """
    Compute L2_HR_PIXCVec filename from L2_HR_LakeTile_PIXCVec filename

    :param in_laketile_pixcvec_filename: L2_HR_LakeTile_PIXCVec filename to convert
    :type in_laketile_pixcvec_filename: string
    :param in_output_dir: output directory
    :type in_output_dir: string
    
    :return: out_filename = L2_HR_PIXCVec full path
    :rtype: out_filename = string
    """
    # Get config file
    cfg = service_config_file.get_instance()
    
    # Init variables
    product_counter = 1
    # Build filename pattern
    pixcvec_pattern = PIXCVEC_PATTERN
    pixcvec_prefix = PIXCVEC_PREFIX
    pixcvec_suffix = PIXCVEC_SUFFIX
    pixcvec_pattern = pixcvec_pattern.replace("PIXCVEC_PREFIX + ", pixcvec_prefix).replace('"', '').replace(" + PIXCVEC_SUFFIX", pixcvec_suffix)
    # Get infos from config file
    crid_lakesp = cfg.get("FILE_INFORMATION", "CRID_LAKESP")
    flag_inc_file_counter = cfg.getboolean("OPTIONS", "Increment file counter")
    # Get infos from input LakeTile_PIXCVec filename
    tmp_dict = get_info_from_filename(in_laketile_pixcvec_filename, "LakeTile")

    # Compute associated PIXCVec filename
    filename = pixcvec_pattern % (int(tmp_dict["cycle"]), int(tmp_dict["pass"]), tmp_dict["tile_ref"], tmp_dict["start_date"], \
                                      tmp_dict["stop_date"], crid_lakesp, product_counter)
    out_filename = os.path.join(in_output_dir, filename)

    # Test existence and modify self.product_counter if wanted and needed
    if flag_inc_file_counter:
        while os.path.exists(out_filename):
            product_counter += 1
            filename = pixcvec_pattern % (int(tmp_dict["cycle"]), int(tmp_dict["pass"]), tmp_dict["tile_ref"], tmp_dict["start_date"], \
                                              tmp_dict["stop_date"], crid_lakesp, product_counter)
            out_filename = os.path.join(in_output_dir, filename)

    return out_filename


#######################################
    

def compute_lakeavg_filename(in_output_dir, in_cycle_number, in_continent_id, in_basin_code):
    """
    Compute L2_HR_LakeAvg filename from cycle number and basin code

    :param in_output_dir: output directory
    :type in_output_dir: string
    :param in_cycle_number: cycle number
    :type in_cycle_number: int
    :param in_continent_id: 2-letter continent identifier
    :type in_continent_id: string
    :param in_basin_code: 2-digit basin code CB
    :type in_basin_code: string
    
    :return: out_filename = L2_HR_LakeAvg full path
    :rtype: out_filename = string
    """
    # Get config file
    cfg = service_config_file.get_instance()
    
    # Init variables
    product_counter = 1
    # Build filename pattern
    lakeavg_pattern = LAKE_AVG_PATTERN
    lakeavg_prefix = LAKE_AVG_PREFIX
    lakeavg_suffix = LAKE_AVG_SUFFIX
    lakeavg_pattern = lakeavg_pattern.replace("LAKE_AVG_PREFIX + ", lakeavg_prefix).replace('"', '').replace(" + LAKE_AVG_SUFFIX", lakeavg_suffix)
    # Get infos from config file
    crid = cfg.get("FILE_INFORMATION", "CRID_LAKEAVG")
    flag_inc_file_counter = cfg.getboolean("OPTIONS", "Increment file counter")

    # Compute LakeAvg filename
    start_datetime = "YYYYMMDDThhmmss"
    stop_datetime = "YYYYMMDDThhmmss"
    filename = lakeavg_pattern % (in_cycle_number, in_continent_id, in_basin_code, \
                                  start_datetime, stop_datetime, \
                                  crid, product_counter)
    out_filename = os.path.join(in_output_dir, filename)

    # Test existence and modify product_counter if wanted and needed
    if flag_inc_file_counter:
        while os.path.exists(out_filename):
            product_counter += 1
            filename = lakeavg_pattern % (in_cycle_number, in_continent_id, in_basin_code, \
                                          start_datetime, stop_datetime, \
                                          crid, product_counter)
            out_filename = os.path.join(in_output_dir, filename)

    return out_filename


#######################################
    

class PldFilenames(object):
    """
    class PldFilenames
    Manage operational Prior Lake Database filenames
    """
    
    def __init__(self, in_proc_type, in_root_dir, cycle_num=None, pass_num=None, basin_code=None):
        """
        Constructor of PLD filenames

        :param in_proc_type: processor type, among "TILE"=LakeTile, "SP"=LakeSP, and "AVG"=LakeAvg
        :type in_proc_type: string
        :param in_root_dir: full path of the PLD root directory
        :type in_root_dir: string
        :param cycle_num: cycle number
        :type cycle_num: int
        :param pass_num: (used for in_proc_type="SP" only) pass number
        :type pass_num: int
        :param basin_code: (used for in_proc_type="AVG" only) 2-digit basin code CB
        :type basin_code: str

        Variables of the object:
            - root_dir / str: full path of the PLD root directory
            - orbit_phase / str: orbit phase, among "calval" and "science", depending on cycle number
            - pass_num / int: pass number
            - basin_code / int: 2-digit basin code CB
            - filename / str: full path of the selected PLD file
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # Init variables
        self.root_dir = in_root_dir
        self.filename = None
        # Orbit phase
        self.orbit_phase = "science"
        if cycle_num is not None:
            self.orbit_phase = my_tools.swot_orbit_phase(cycle_num)

        # Filename depends on processing type
        if (in_proc_type == "TILE") or (in_proc_type == "SP"):
            self.pass_num = pass_num
            self.compute_pld_filename_pass()
            
        elif in_proc_type == "AVG":
            self.basin_code = basin_code
            self.compute_pld_filename_basin()
            
        else:
            message = "Processor type %s should be TILE, SP, or AVG" % in_proc_type
            raise service_error.ProcessingError(message, logger)

    #----------------------------------

    def compute_pld_filename_pass(self):
        """
        Compute PLD full path to be used for LakeTile or LakeSP processing
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Filename pattern
        filename_pattern = PLD_PATTERN_START[self.orbit_phase] % self.pass_num
        
        # 2 - Directory
        if self.orbit_phase == "science":
            pld_dir = os.path.join(self.root_dir, "Nom")
        else:
            pld_dir = os.path.join(self.root_dir, "Cal")
        
        # 3 - Full path
        # 3.1 - List of PLD files in directory
        tmp_list = os.listdir(pld_dir)
        # 3.2 - Retrieve PLD file corresponding to filename pattern
        for cur_file in tmp_list:
            # Test if file meets the condition
            if cur_file.startswith(filename_pattern) and cur_file.endswith(PLD_SUFFIX):
                self.filename = os.path.join(pld_dir, cur_file)
                break
            
        if self.filename is None:
            logger.warning("PLD file starting with %s has not been found" % filename_pattern)

    def compute_pld_filename_basin(self):
        """
        Compute PLD full path to be used for LakeAvg processing
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Filename pattern
        continent_id = lake_db.compute_continent_id_from_basin_code(self.basin_code)
        filename_pattern = PLD_PATTERN_START["basin"] % (continent_id, self.basin_code[0:2])
        
        # 2 - Directory
        pld_dir = os.path.join(self.root_dir, "Basin")
        
        # 3.1 - List of PLD files in directory
        tmp_list = os.listdir(pld_dir)
        # 3.2 - Retrieve PLD file corresponding to filename pattern
        for cur_file in tmp_list:
            # Test if file meets the condition
            if cur_file.startswith(filename_pattern) and cur_file.endswith(PLD_SUFFIX):
                self.filename = os.path.join(pld_dir, cur_file)
                break
            
        if self.filename is None:
            logger.warning("PLD file starting with %s has not been found" % filename_pattern)
        
