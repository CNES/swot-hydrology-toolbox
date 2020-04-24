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
LAKE_TILE_PREFIX = "SWOT_L2_HR_LakeTile_"
# LakeTile filename with %03d=cycle number %03d=pass number %s=tile ref %s=swath %s=begin date %s=end date %s=CRID %s=counter %s=suffix
LAKE_TILE_PATTERN = LAKE_TILE_PREFIX + "%03d_%03d_%s_%s_%s_%s_%02d%s"
LAKE_TILE_PATTERN_PRINT = LAKE_TILE_PREFIX + "%s<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>"
# Indices when LAKE_TILE_*_PATTERN.split("_"); None if value not in filename
LAKE_TILE_PATTERN_IND = {"cycle": 4, "pass": 5, "tile_ref": 6, "start_date": 7, "stop_date": 8, "crid": 9, "counter": 10}
LAKE_TILE_SHP_SUFFIX = ".shp"
LAKE_TILE_SHP_META_SUFFIX = ".shp.xml"
LAKE_TILE_EDGE_SUFFIX = "_edge.nc"
LAKE_TILE_PIXCVEC_SUFFIX = "_pixcvec.nc"

# 4 - LakeSP
LAKE_SP_PREFIX = "SWOT_L2_HR_LakeSP_"
# LakeSP filename with %03d=cycle number %03d=pass number %s=continent %s=begin date %s=end date %s=CRID %s=counter
LAKE_SP_PATTERN = LAKE_SP_PREFIX + "%03d_%03d_%s_%s_%s_%s_%02d.shp"
# LakeSP filename without continent info with %03d=cycle number %03d=pass number %s=begin date %s=end date %s=CRID %s=counter
LAKE_SP_PATTERN_NO_CONT = LAKE_SP_PREFIX + "%03d_%03d_%s_%s_%s_%02d.shp"

# 5 - PIXCVec
PIXCVEC_PREFIX = "SWOT_L2_HR_PIXCVec_"
PIXCVEC_SUFFIX = ".nc"
# PIXCVec filename with %03d=cycle number %03d=pass number %s=tile ref %s=begin date %s=end date %s=CRID %s=counter
PIXCVEC_PATTERN = PIXCVEC_PREFIX + "%03d_%03d_%s_%s_%s_%s_%02d" + PIXCVEC_SUFFIX

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
            logger.info("Filename %s doesn't match with PIXC pattern %s",
                        os.path.basename(in_filename), PIXC_PATTERN_PRINT)
            pattern = None
        else:
            pattern = PIXC_PATTERN_IND
    elif in_type == "PIXCVecRiver":
        if not os.path.basename(in_filename).startswith(PIXCVEC_RIVER_PREFIX):
            logger.info("Filename %s doesn't match with PIXCVecRiver pattern %s",
                        os.path.basename(in_filename), PIXCVEC_RIVER_PATTERN_PRINT)
            pattern = None
        else:
            pattern = PIXCVEC_RIVER_PATTERN_IND
    elif in_type == "LakeTile":
        if not os.path.basename(in_filename).startswith(LAKE_TILE_PREFIX):
            message = "Filename %s doesn't match with LakeTile pattern %s" % (os.path.basename(in_filename), LAKE_TILE_PATTERN_PRINT)
            raise service_error.SASLakeTileError(message, logger)
        pattern = LAKE_TILE_PATTERN_IND
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
    
    def __init__(self, in_pixc_file, in_pixc_vec_river_file, in_out_dir):
        """
        Constructor of LakeTile filenames

        :param IN_pixc_file: PIXC file full path
        :type IN_pixc_file: string
        :param IN_pixc_vec_river_file: PIXCVecRiver file full path
        :type IN_pixc_vec_river_file: string
        :param IN_out_dir: output directory
        :type IN_out_dir: string
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        cfg = service_config_file.get_instance()  # Get config file

        # 1 - Init variables
        pixc_file = in_pixc_file  # PIXC file full path
        self.pixc_vec_river_file = in_pixc_vec_river_file  # PIXCVecRiver full path
        self.out_dir = in_out_dir  # Output directory
        self.product_counter = 1  # Product counter
        # Get CRID from configuration file
        self.lake_tile_crid = cfg.get("CRID", "LAKE_TILE_CRID")

        # 2 - Retrieve info from PIXC filename
        tmp_dict = get_info_from_filename(pixc_file, "PIXC")
        # 2.1 - Cycle number
        cycle_num = tmp_dict["cycle"]
        if cycle_num is None:
            self.cycle_num = 999
            logger.info("WARNING: cycle number has not been found in PIXC filename %s -> set to default value = %03d",
                        os.path.basename(pixc_file), self.cycle_num)
        else:
            self.cycle_num = int(cycle_num)
            logger.info("Cycle number = %03d", self.cycle_num)
        # 2.2 - Pass number
        pass_num = tmp_dict["pass"]
        if pass_num is None:
            self.pass_num = 999
            logger.info("WARNING: pass number has not been found in PIXC filename %s -> set to default value = %03d",
                        os.path.basename(pixc_file), self.pass_num)
        else:
            self.pass_num = int(pass_num)
            logger.info("Pass number = %03d", self.pass_num)
        # 2.3 - Tile ref
        self.tile_ref = tmp_dict["tile_ref"]
        if self.tile_ref is None:
            self.tile_ref = "ttts"
            logger.info("WARNING: tile ref has not been found in PIXC filename %s -> set to default value = %s",
                        os.path.basename(pixc_file), self.tile_ref)
        else:
            logger.info("Tile ref = %s", self.tile_ref)
        # 2.4 - Start date
        self.start_date = tmp_dict["start_date"]
        if self.start_date is None:
            self.start_date = "yyyyMMddThhmmss"
            logger.info("WARNING: start date has not been found in PIXC filename %s -> set to default value = %s",
                        os.path.basename(pixc_file), self.start_date)
        else:
            logger.info("Start date = %s", self.start_date)
        # 2.5 - Stop date
        self.stop_date = tmp_dict["stop_date"]
        if self.stop_date is None:
            self.stop_date = "yyyyMMddThhmmss"
            logger.info("WARNING: stop date has not been found in PIXC filename %s -> set to default value = %s",
                        os.path.basename(pixc_file), self.stop_date)
        else:
            logger.info("Stop date = %s", self.stop_date)

        # 3 - Test compatibility of PIXCVecRiver filename with PIXC filename
        tmp_ok = self.test_pixc_vec_river_filename()
        if tmp_ok:
            logger.info("PIXCVecRiver basename %s is compatible with PIXC basename %s",
                        os.path.basename(self.pixc_vec_river_file), os.path.basename(pixc_file))
        else:
            logger.info("WARNING: PIXCVecRiver basename %s IS NOT compatible with PIXC basename %s (cf. above)",
                        os.path.basename(self.pixc_vec_river_file), os.path.basename(pixc_file))

        # 4 - Compute output filenames
        self.compute_lake_tile_filename_shp()  # LakeTile_shp filename
        self.compute_lake_tile_filename_edge()  # LakeTile_edge filename
        self.compute_lake_tile_filename_pixcvec()  # LakeTile_pixcvec filename

    #----------------------------------

    def compute_lake_tile_filename_shp(self):
        """
        Compute LakeTile_shp full path
        """
        
        filename = LAKE_TILE_PATTERN % (self.cycle_num, self.pass_num, self.tile_ref, self.start_date, self.stop_date, \
                                             self.lake_tile_crid, self.product_counter, LAKE_TILE_SHP_SUFFIX)
        self.lake_tile_shp_file = os.path.join(self.out_dir, filename)

        # Test existence and modify self.product_counter if needed
        while os.path.exists(self.lake_tile_shp_file):
            self.product_counter += 1
            filename = LAKE_TILE_PATTERN % (self.cycle_num, self.pass_num, self.tile_ref, self.start_date, self.stop_date, \
                                                 self.lake_tile_crid, self.product_counter, LAKE_TILE_SHP_SUFFIX)
            self.lake_tile_shp_file = os.path.join(self.out_dir, filename)

    def compute_lake_tile_filename_edge(self):
        """
        Compute LakeTile_edge full path
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        filename = LAKE_TILE_PATTERN % (self.cycle_num, self.pass_num, self.tile_ref, self.start_date, self.stop_date, self.lake_tile_crid, \
                                             self.product_counter, LAKE_TILE_EDGE_SUFFIX)
        self.lake_tile_edge_file = os.path.join(self.out_dir, filename)

        # Filename must not exist
        if os.path.exists(self.lake_tile_edge_file):
            message = "ERROR = %s already exists" % self.lake_tile_edge_file
            raise service_error.SASLakeTileError(message, logger)

    def compute_lake_tile_filename_pixcvec(self):
        """
        Compute LakeTile_pixcvec full path
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        filename = LAKE_TILE_PATTERN % (self.cycle_num, self.pass_num, self.tile_ref, self.start_date, self.stop_date, self.lake_tile_crid, \
                                             self.product_counter, LAKE_TILE_PIXCVEC_SUFFIX)
        self.lake_tile_pixcvec_file = os.path.join(self.out_dir, filename)

        # Filename must not exist
        if os.path.exists(self.lake_tile_pixcvec_file):
            message = "ERROR = %s already exists" % self.lake_tile_pixcvec_file
            raise service_error.SASLakeTileError(message, logger)

    #----------------------------------

    def test_pixc_vec_river_filename(self):
        """
        Test if PIXCVecRiver filename is coherent with PIXC file, wrt (cycle, pass, tile) triplet

        :return: True if it's the case, or False if not
        :rtype: boolean
        """
        logger = logging.getLogger(self.__class__.__name__)

        # 0 - Init output boolean
        out_ok = True

        # 1 - Retrieve info from PIXCVecRiver filename
        tmp_dict = get_info_from_filename(self.pixc_vec_river_file, "PIXCVecRiver")
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
                logger.info("WARNING: cycle number is not the same in PIXC (=%03d) and PIXCVecRiver (=%03d) filenames", \
                            self.cycle_num, pixcvec_cycle_num)
            # 2.2 - Pass number
            if pixcvec_pass_num != self.pass_num:
                out_ok = False
                logger.info("WARNING: pass number is not the same in PIXC (=%03d) and PIXCVecRiver (=%03d) filenames", \
                            self.pass_num, pixcvec_pass_num)
            # 2.3 - Tile ref
            if pixcvec_tile_ref != self.tile_ref:
                out_ok = False
                logger.info("WARNING: tile ref not the same in PIXC (=%s) and PIXCVecRiver (=%s) filenames", self.tile_ref, \
                            pixcvec_tile_ref)

        return out_ok


#######################################


class LakeSPFilenames(object):
    """
        class LakeSPFilenames
    """
    def __init__(self, in_lake_tile_file_list, in_continent, in_out_dir):
        """
        Constructor of LakeSP filenames

        :param in_lake_tile_file_list: list of LakeTile_shp files full path
        :type in_lake_tile_file_list: list of string
        :param in_continent: continent concerning the LakeSP
        :type in_continent: string
        :param in_out_dir: output directory
        :type in_out_dir: string
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        # Get config file
        cfg = service_config_file.get_instance()

        # 1 - Init variables
        self.lake_tile_file_list = in_lake_tile_file_list  # List of LakeTile_shp files full path
        self.out_dir = in_out_dir  # Output directory
        self.product_counter = 1  # Product counter
        # Get CRID from configuration file
        self.lake_sp_crid = cfg.get("CRID", "LAKE_SP_CRID")

        # 2 - Retrieve info from LakeTile filename
        tmp_dict = get_info_from_filename(self.lake_tile_file_list[0], "LakeTile")
        # 2.1 - Cycle number
        cycle_num = tmp_dict["cycle"]
        if cycle_num is None:
            self.cycle_num = 999
            logger.info("WARNING: cycle number has not been found in LakeTile filename %s -> set to default value = %03d",
                        os.path.basename(self.lake_tile_file_list[0]), self.cycle_num)
        else:
            self.cycle_num = int(cycle_num)
            logger.info("Cycle number = %03d", self.cycle_num)
        # 2.2 - Pass number
        pass_num = int(tmp_dict["pass"])
        if pass_num is None:
            self.pass_num = 999
            logger.info("WARNING: pass number has not been found in LakeTile filename %s -> set to default value = %03d",
                        os.path.basename(self.lake_tile_file_list[0]), self.pass_num)
        else:
            self.pass_num = int(pass_num)
            logger.info("Pass number = %03d", self.pass_num)
        # 2.3 - Continent
        if in_continent == "":
            self.continent = ""
            logger.info("WARNING: no continental split")
        else:
            self.continent = in_continent
            logger.info("Continent = %s", self.continent)

        # 3 - Retrieve start and stop dates from LakeTile_shp filenames
        self.compute_start_stop_dates()

        # 4 - Compute output filenames
        self.compute_lake_sp_filename()  # LakeSP filename

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

    def compute_lake_sp_filename(self):
        """
        Compute LakeSP full path
        """

        if self.continent is "":
            filename = LAKE_SP_PATTERN_NO_CONT % (self.cycle_num, self.pass_num, self.start_date, self.stop_date, self.lake_sp_crid, self.product_counter)
        else:
            filename = LAKE_SP_PATTERN % (self.cycle_num, self.pass_num, self.continent, self.start_date, self.stop_date, self.lake_sp_crid, self.product_counter)
        self.lake_sp_file = os.path.join(self.out_dir, filename)

        # Test existence and modify self.product_counter if needed
        while os.path.exists(self.lake_sp_file):
            self.product_counter += 1
            if self.continent is "":
                filename = LAKE_SP_PATTERN_NO_CONT % (self.cycle_num, self.pass_num, self.start_date, self.stop_date, self.lake_sp_crid, self.product_counter)
            else:
                filename = LAKE_SP_PATTERN % (self.cycle_num, self.pass_num, self.continent, self.start_date, self.stop_date, self.lake_sp_crid, self.product_counter)
            self.lake_sp_file = os.path.join(self.out_dir, filename)


#######################################


def compute_pixcvec_filename(in_laketile_pixcvec_filename, in_output_dir):
    """
    Compute L2_HR_PIXCVec filename from L2_HR_LakeTile_pixcvec filename

    :param in_laketile_pixcvec_filename: L2_HR_LakeTile_pixcvec filename to convert
    :type IN_laketile_filename: string
    :param in_output_dir: output directory
    :type in_output_dir: string
    """
    # Get config file
    cfg = service_config_file.get_instance()
    # Init variables
    product_counter = 1
    # get parameters
    pixcvec_pattern = PIXCVEC_PATTERN
    pixcvec_prefix = PIXCVEC_PREFIX
    pixcvec_suffix = PIXCVEC_SUFFIX
    pixcvec_pattern = pixcvec_pattern.replace("PIXCVEC_PREFIX + ", pixcvec_prefix).replace('"', '').replace(" + PIXCVEC_SUFFIX", pixcvec_suffix)
    lake_sp_crid = cfg.get("CRID", "LAKE_SP_CRID")
    # Get infos from input LakeTile_pixcvec filename
    tmp_dict = get_info_from_filename(in_laketile_pixcvec_filename, "LakeTile")

    # Compute associated PIXCVec filename
    filename = pixcvec_pattern % (int(tmp_dict["cycle"]), int(tmp_dict["pass"]), tmp_dict["tile_ref"], tmp_dict["start_date"], \
                                      tmp_dict["stop_date"], lake_sp_crid, product_counter)
    out_filename = os.path.join(in_output_dir, filename)

    while os.path.exists(out_filename):
        product_counter += 1
        filename = pixcvec_pattern % (int(tmp_dict["cycle"]), int(tmp_dict["pass"]), tmp_dict["tile_ref"], tmp_dict["start_date"], \
                                          tmp_dict["stop_date"], lake_sp_crid, product_counter)
        out_filename = os.path.join(in_output_dir, filename)

    return out_filename
