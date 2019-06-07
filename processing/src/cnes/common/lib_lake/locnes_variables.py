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
.. module:: locnes_variables.py
    :synopsis: Gather global variables used by LOCNES
     Created on 2018/08/21

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from __future__ import absolute_import, division, print_function, unicode_literals 
import logging

import cnes.common.service_config_file as service_config_file
import cnes.common.service_error as service_error


# Lake a priori database
LAKE_DB = "/work/ALT/swot/swotpub/BD/BD_lakes/20181002_EU/apriori_db_lakes_EU.shp"
INFLUENCE_LAKE_DB = "/work/ALT/swot/swotpub/BD/BD_lakes/20181002_EU/influence_area_map.shp"

# Lake identifier attribute name in the database
LAKE_DB_ID = "lake_id"

# Shapefile with polygons of continents
CONTINENT_FILE = "/work/ALT/swot/swotpub/BD/major_basins/FAO/major_hydrobasins.shp"

# Flags
FLAG_WATER = "3;4"  # Water flag  3=water near land edge  4=interior water
FLAG_DARK = "23;24"  # Dark water flag  23=darkwater near land  24=interior dark water

# Min size for a lake to generate a lake product (=polygon + attributes) for it
MIN_SIZE = 10000.0  # In m2

# To improve PixC golocation (=True) or not (=False)
IMP_GEOLOC = True

# Method to compute lake boundary or polygon hull
# 0 = convex hull 
# 1.0 = concave hull computed in ground geometry, based on Delaunay triangulation - using CGAL library 
# 1.1 = concave hull computed in ground geometry, based on Delaunay triangulation - with alpha parameter varying across-track
# 2 = edge computed in radar geometry, then converted in ground geometry (default)
HULL_METHOD = 2
# If HULL_METHOD=1.1: max number of pixels for Delaunay triangulation
NB_PIX_MAX_DELAUNEY = 1e5
# If HULL_METHOD=2: max number of contour points
NB_PIX_MAX_CONTOUR = 8000  

# Maximal standard deviation of height inside a lake  (-1.0 = do not compute lake height segmentation)
STD_HEIGHT_MAX = -1.0

# Big lakes parameters for improved geoloc; used only if imp geoloc=1
BIGLAKE_MODEL = "polynomial"  # =polynomial or =grid
BIGLAKE_MIN_SIZE = 50000000  # In m^2; if None, disable biglake model
BIGLAKE_GRID_SPACING = 4000  # Grid spacing for lake height smoothing; in m
BIGLAKE_GRID_RES = 8000  # Grid resolution for lake height smoothing; in m


# Filenames pattern
PRODUCER = "CNES"  # Product generator
LAKE_TILE_CRID = "Dx0000"  # Composite Release IDentifier for LakeTile processing
LAKE_SP_CRID = "Dx0000"  # Composite Release IDentifier for LakeSP processing

PIXC_PREFIX = "SWOT_L2_HR_PIXC_"
PIXC_PATTERN_PRINT = PIXC_PREFIX + "<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>.nc"
PIXC_PATTERN_IND = {"cycle": 4, "pass": 5, "tile_ref": 6, "start_date": 7, "stop_date": 8, "crid": 9, "counter": 10}  # Indices when PIXC_PATTERN.split("_"); None if value not in filename

PIXCVEC_RIVER_PREFIX = "SWOT_L2_HR_PIXCVecRiver_"
PIXCVEC_RIVER_PATTERN_PRINT = PIXCVEC_RIVER_PREFIX + "<CycleID>_<PassID>_<TileID>[L/R]_<RangeBeginDateTime>_<RangeEndingDateTime>_<CRID>_<ProductCounter>.nc"
PIXCVEC_RIVER_PATTERN_IND = {"cycle": 4, "pass": 5, "tile_ref": 6, "start_date": 7, "stop_date": 8, "crid": 9, "counter": 10}  # Indices when PIXCVEC_RIVER_PATTERN.split("_"); None if value not in filename

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

# Nb digits for counter of lakes in a tile or pass
NB_DIGITS = 4


# ----------------------------------------
# Parameter overwrite functions
# ----------------------------------------

def tmpGetConfigFromServiceConfigFile():
    """
    Set global variables from serviceConfgiFile
    This function is temporary. It will be delete in the future
    when serviceConfigFile will be used by lake_tile
    """
    cfg = service_config_file.get_instance()
    logger = logging.getLogger("locnes_variables")


    # Lake database
    global LAKE_DB
    lake_db_file = cfg.get("DATABASES", "LAKE_DB")
    LAKE_DB = lake_db_file
    logger.info("> LAKE_DB = %s" % LAKE_DB)

    # Influence lake database
    global INFLUENCE_LAKE_DB
    influence_lake_db_file = cfg.get("DATABASES", "INFLUENCE_LAKE_DB")
    INFLUENCE_LAKE_DB = influence_lake_db_file
    logger.info("> INFLUENCE_LAKE_DB = %s" % INFLUENCE_LAKE_DB)

    # Lake identifier attribute name in the database
    global LAKE_DB_ID
    lake_db_id = cfg.get("DATABASES", "LAKE_DB_ID")
    LAKE_DB_ID = lake_db_id
    logger.info("> LAKE_DB_ID = %s" % LAKE_DB_ID)

    # Shapefile with polygons of continents
    global CONTINENT_FILE
    continent_file = cfg.get("DATABASES", "CONTINENT_FILE")
    CONTINENT_FILE = continent_file
    logger.info("> CONTINENT_FILE = %s" % CONTINENT_FILE)

    # Water flags
    global FLAG_WATER
    FLAG_WATER = cfg.get("CONFIG_PARAMS", "FLAG_WATER")
    logger.info("> FLAG_WATER = %s" % FLAG_WATER)

    # Dark water flags
    global FLAG_DARK
    FLAG_DARK = cfg.get("CONFIG_PARAMS", "FLAG_DARK")
    logger.info("> FLAG_DARK = %s" % FLAG_DARK)

    # Layover flags
    global FLAG_LAYOVER
    FLAG_LAYOVER = cfg.get("CONFIG_PARAMS", "FLAG_LAYOVER")
    logger.info("> FLAG_LAYOVER = %s" % FLAG_LAYOVER)

    # Hull method
    global HULL_METHOD
    HULL_METHOD = cfg.getfloat("CONFIG_PARAMS", "HULL_METHOD")
    logger.info("> HULL_METHOD = %s" % HULL_METHOD)

    # Maximal standard deviation of height inside a lake
    global STD_HEIGHT_MAX
    STD_HEIGHT_MAX = cfg.getfloat("CONFIG_PARAMS", "STD_HEIGHT_MAX")
    logger.info("> STD_HEIGHT_MAX = %s" % STD_HEIGHT_MAX)
    
    # Model to deal with big lake processing
    global BIGLAKE_MODEL
    BIGLAKE_MODEL = cfg.get("CONFIG_PARAMS", "BIGLAKE_MODEL")
    logger.info("> BIGLAKE_MODEL = %s" % BIGLAKE_MODEL)

    # Min size for lake to be considered as big
    global BIGLAKE_MIN_SIZE
    BIGLAKE_MIN_SIZE = cfg.getfloat("CONFIG_PARAMS", "BIGLAKE_MIN_SIZE")
    logger.info("> BIGLAKE_MIN_SIZE = %s" % BIGLAKE_MIN_SIZE)

    # Grid spacing for lake height smoothing
    global BIGLAKE_GRID_SPACING
    BIGLAKE_GRID_SPACING = cfg.getint("CONFIG_PARAMS", "BIGLAKE_GRID_SPACING")
    logger.info("> BIGLAKE_GRID_SPACING = %s" % BIGLAKE_GRID_SPACING)

    # Grid resolution for lake height smoothing
    global BIGLAKE_GRID_RES
    BIGLAKE_GRID_RES = cfg.getint("CONFIG_PARAMS", "BIGLAKE_GRID_RES")
    logger.info("> BIGLAKE_GRID_RES = %s" % BIGLAKE_GRID_RES)


def overwriteConfig_from_cfg(IN_config):
    """
    Set global variables if overwritten in a parameter file
    
    :param IN_config: configuration parameters
    :type IN_config: ConfigParser.reader
    """
    
    if "CONFIG_PARAMS" in IN_config.sections():
        
        list_over = IN_config.options("CONFIG_PARAMS")
        
        # Lake database
        if "lake_db" in list_over:
            global LAKE_DB
            lake_db_file = IN_config.get("CONFIG_PARAMS", "LAKE_DB")
            import cnes.common.lib.my_tools as my_tools
            my_tools.testFile(lake_db_file)  # Test existence of file
            LAKE_DB = lake_db_file
            print("> LAKE_DB = %s" % LAKE_DB)
        else:
            print("> Default value for LAKE_DB = %s" % LAKE_DB)

        # Influence lake database
        if "influence_lake_db" in list_over:
            global INFLUENCE_LAKE_DB
            influence_lake_db_file = IN_config.get("CONFIG_PARAMS", "INFLUENCE_LAKE_DB")
            import cnes.common.lib.my_tools as my_tools
            my_tools.testFile(influence_lake_db_file)  # Test existence of file
            INFLUENCE_LAKE_DB = influence_lake_db_file
            print("> INFLUENCE_LAKE_DB = %s" % INFLUENCE_LAKE_DB)
        else:
            print("> Default value for INFLUENCE_LAKE_DB = %s" % INFLUENCE_LAKE_DB)

        # Lake identifier attribute name in the database
        if "lake_db_id" in list_over:
            global LAKE_DB_ID
            lake_db_id = IN_config.get("CONFIG_PARAMS", "LAKE_DB_ID")
            LAKE_DB_ID = lake_db_id
            print("> LAKE_DB_ID = %s" % LAKE_DB_ID)
        else:
            print("> Default value for LAKE_DB_ID = %s" % LAKE_DB_ID)
        
        # Shapefile with polygons of continents
        if "continent_file" in list_over:
            global CONTINENT_FILE
            continent_file = IN_config.get("CONFIG_PARAMS", "CONTINENT_FILE")
            import cnes.common.lib.my_tools as my_tools
            my_tools.testFile(continent_file)  # Test existence of file
            CONTINENT_FILE = continent_file
            print("> CONTINENT_FILE = %s" % CONTINENT_FILE)
        else:
            try:
                print("> Default value for CONTINENT_FILE = %s" % CONTINENT_FILE)
            except:
                CONTINENT_FILE = None
                print("> No value for CONTINENT_FILE => LAKE_TILE product not linked to a continent")
            
        # Water flags
        if "flag_water" in list_over:
            global FLAG_WATER
            FLAG_WATER = IN_config.get("CONFIG_PARAMS", "FLAG_WATER")
            print("> FLAG_WATER = %s" % FLAG_WATER)
        else:
            print("> Default value for FLAG_WATER = %s" % FLAG_WATER)
            
        # Dark water flags
        if "flag_dark" in list_over:
            global FLAG_DARK
            FLAG_DARK = IN_config.get("CONFIG_PARAMS", "FLAG_DARK")
            print("> FLAG_DARK = %s" % FLAG_DARK)
        else:
            print("> Default value for FLAG_DARK = %s" % FLAG_DARK)
            
        # Hull method
        if "hull_method" in list_over:
            global HULL_METHOD
            HULL_METHOD = IN_config.getfloat("CONFIG_PARAMS", "HULL_METHOD")
            print("> HULL_METHOD = %s" % HULL_METHOD)
        else:
            print("> Default value for HULL_METHOD = %s" % HULL_METHOD)

        # Std height max
        if "std_height_max" in list_over:
            global STD_HEIGHT_MAX
            STD_HEIGHT_MAX = IN_config.getfloat("CONFIG_PARAMS", "STD_HEIGHT_MAX")
            print("> STD_HEIGHT_MAX = %s" % STD_HEIGHT_MAX)
        else:
            print("> Default value for STD_HEIGHT_MAX = %s" % STD_HEIGHT_MAX)
                        
        # Model to deal with big lake processing
        if "biglake_model" in list_over:
            global BIGLAKE_MODEL
            BIGLAKE_MODEL = IN_config.get("CONFIG_PARAMS", "BIGLAKE_MODEL")
            print("> BIGLAKE_MODEL = %s" % BIGLAKE_MODEL)
        else:
            print("> Default value for BIGLAKE_MODEL = %s" % BIGLAKE_MODEL)
        
        # Min size for lake to be considered as big
        if "biglake_min_size" in list_over:
            global BIGLAKE_MIN_SIZE
            BIGLAKE_MIN_SIZE = IN_config.getfloat("CONFIG_PARAMS", "BIGLAKE_MIN_SIZE")
            print("> BIGLAKE_MIN_SIZE = %s" % BIGLAKE_MIN_SIZE)
        else:
            print("> Default value for BIGLAKE_MIN_SIZE = %s" % BIGLAKE_MIN_SIZE)
        
        # Grid spacing for lake height smoothing
        if "biglake_grid_spacing" in list_over:
            global BIGLAKE_GRID_SPACING
            BIGLAKE_GRID_SPACING = IN_config.getint("CONFIG_PARAMS", "BIGLAKE_GRID_SPACING")
            print("> BIGLAKE_GRID_SPACING = %s" % BIGLAKE_GRID_SPACING)
        else:
            print("> Default value for BIGLAKE_GRID_SPACING = %s" % BIGLAKE_GRID_SPACING)
        
        # Grid resolution for lake height smoothing
        if "biglake_grid_res" in list_over:
            global BIGLAKE_GRID_RES
            BIGLAKE_GRID_RES = IN_config.getint("CONFIG_PARAMS", "BIGLAKE_GRID_RES")
            print("> BIGLAKE_GRID_RES = %s" % BIGLAKE_GRID_RES)
        else:
            print("> Default value for BIGLAKE_GRID_RES = %s" % BIGLAKE_GRID_RES)


def overwriteConfig_from_xml(IN_xml_tree):
    """
    Set global variables from an XML tree
    
    :param IN_xml_tree: XML tree (typically from .shp.xml file)
    :type IN_xml_tree: etree.parse
    """
    
    # Lake database
    global LAKE_DB
    lake_db_file = IN_xml_tree.xpath("//LakeTile_shp/config_params/lake_db")[0].text
    import cnes.common.lib.my_tools as my_tools
    my_tools.testFile(lake_db_file)  # Test existence of file
    LAKE_DB = lake_db_file
    print("> LAKE_DB = %s" % LAKE_DB)

    # Influence lake database
    global INFLUENCE_LAKE_DB
    influence_lake_db_file = IN_xml_tree.xpath("//LakeTile_shp/config_params/influence_lake_db")[0].text
    import cnes.common.lib.my_tools as my_tools
    my_tools.testFile(influence_lake_db_file)  # Test existence of file
    INFLUENCE_LAKE_DB = influence_lake_db_file
    print("> INFLUENCE_LAKE_DB = %s" % INFLUENCE_LAKE_DB)
        
    # Lake identifier attribute name in the database
    global LAKE_DB_ID
    lake_db_id = IN_xml_tree.xpath("//LakeTile_shp/config_params/lake_db_id")[0].text
    LAKE_DB_ID = lake_db_id
    print("> LAKE_DB_ID = %s" % LAKE_DB_ID)
        
    # Shapefile with polygons of continents
    global CONTINENT_FILE
    try:
        continent_file = IN_xml_tree.xpath("//LakeTile_shp/config_params/continent_file")[0].text
        import cnes.common.lib.my_tools as my_tools
        my_tools.testFile(continent_file)  # Test existence of file
        CONTINENT_FILE = continent_file
        print("> CONTINENT_FILE = %s" % CONTINENT_FILE)
    except:
        CONTINENT_FILE = None
        print("> No value for CONTINENT_FILE => LAKE_SP product not linked to a continent")
            
    # Water flags
    global FLAG_WATER
    FLAG_WATER = IN_xml_tree.xpath("//LakeTile_shp/config_params/flag_water")[0].text
    print("> FLAG_WATER = %s" % FLAG_WATER)
            
    # Dark water flags
    global FLAG_DARK
    FLAG_DARK = IN_xml_tree.xpath("//LakeTile_shp/config_params/flag_dark")[0].text
    print("> FLAG_DARK = %s" % FLAG_DARK)
            
    # Min size for lake product computation
    global MIN_SIZE
    MIN_SIZE = float(IN_xml_tree.xpath("//LakeTile_shp/config_params/min_size")[0].text)
    print("> MIN_SIZE = %s" % MIN_SIZE)
    
    # Improve geolocation or not 
    global IMP_GEOLOC
    IMP_GEOLOC = bool(IN_xml_tree.xpath("//LakeTile_shp/config_params/imp_geoloc")[0].text)
    print("> IMP_GEOLOC = %s" % IMP_GEOLOC)
            
    # Hull method
    global HULL_METHOD
    HULL_METHOD = float(IN_xml_tree.xpath("//LakeTile_shp/config_params/hull_method")[0].text)
    print("> HULL_METHOD = %s" % HULL_METHOD)
            
    # Maximal standard deviation of height inside a lake
    global STD_HEIGHT_MAX
    STD_HEIGHT_MAX = float(IN_xml_tree.xpath("//LakeTile_shp/config_params/std_height_max")[0].text)
    print("> STD_HEIGHT_MAX = %s" % STD_HEIGHT_MAX)
            
    # Model to deal with big lake processing
    global BIGLAKE_MODEL
    BIGLAKE_MODEL = IN_xml_tree.xpath("//LakeTile_shp/config_params/biglake_model")[0].text
    print("> BIGLAKE_MODEL = %s" % BIGLAKE_MODEL)
        
    # Min size for lake to be considered as big
    global BIGLAKE_MIN_SIZE
    BIGLAKE_MIN_SIZE = float(IN_xml_tree.xpath("//LakeTile_shp/config_params/biglake_min_size")[0].text)
    print("> BIGLAKE_MIN_SIZE = %s" % BIGLAKE_MIN_SIZE)
        
    # Grid spacing for lake height smoothing
    global BIGLAKE_GRID_SPACING
    BIGLAKE_GRID_SPACING = int(IN_xml_tree.xpath("//LakeTile_shp/config_params/biglake_grid_spacing")[0].text)
    print("> BIGLAKE_GRID_SPACING = %s" % BIGLAKE_GRID_SPACING)
        
    # Grid resolution for lake height smoothing
    global BIGLAKE_GRID_RES
    BIGLAKE_GRID_RES = int(IN_xml_tree.xpath("//LakeTile_shp/config_params/biglake_grid_res")[0].text)
    print("> BIGLAKE_GRID_RES = %s" % BIGLAKE_GRID_RES)


def compareConfig_to_xml(IN_xml_tree):
    """
    Compare global variables to the input XML tree
    
    :param IN_xml_tree: XML tree (typically from .shp.xml file)
    :type IN_xml_tree: etree.parse
    """
    logger = logging.getLogger("locnes_variables")
    
    # Lake database
    TMP_lake_db = IN_xml_tree.xpath("//LakeTile_shp/config_params/lake_db")[0].text
    if TMP_lake_db != LAKE_DB:
        message = "At least 2 different values of LAKE_DB for one processing: %s vs %s" % (LAKE_DB, TMP_lake_db)
        raise service_error.ProcessingError(message, logger)

    # Influence lake database
    TMP_influence_lake_db = IN_xml_tree.xpath("//LakeTile_shp/config_params/influence_lake_db")[0].text
    if TMP_influence_lake_db != INFLUENCE_LAKE_DB:
        message = "At least 2 different values of INFLUENCE_LAKE_DB for one processing: %s vs %s" % (INFLUENCE_LAKE_DB, TMP_influence_lake_db)
        raise service_error.ProcessingError(message, logger)

    # Lake identifier attribute name in the database
    TMP_lake_db_id = IN_xml_tree.xpath("//LakeTile_shp/config_params/lake_db_id")[0].text
    if TMP_lake_db_id != LAKE_DB_ID:
        message = "At least 2 different values of LAKE_DB_ID for one processing: %s vs %s" % (LAKE_DB_ID, TMP_lake_db_id)
        raise service_error.ProcessingError(message, logger)
        
    # Shapefile with polygons of continents
    try:
        TMP_continent_file = IN_xml_tree.xpath("//LakeTile_shp/config_params/continent_file")[0].text
    except: 
        TMP_continent_file = None
    if TMP_continent_file != CONTINENT_FILE:
        message = "At least 2 different values of CONTINENT_FILE for one processing: %s vs %s" % (CONTINENT_FILE, TMP_continent_file)
        raise service_error.ProcessingError(message, logger)
            
    # Water flags
    TMP_flag_water = IN_xml_tree.xpath("//LakeTile_shp/config_params/flag_water")[0].text
    if TMP_flag_water != FLAG_WATER:
        message = "At least 2 different values of FLAG_WATER for one processing: %s vs %s" % (FLAG_WATER, TMP_flag_water)
        raise service_error.ProcessingError(message, logger)
            
    # Dark water flags
    TMP_flag_dark = IN_xml_tree.xpath("//LakeTile_shp/config_params/flag_dark")[0].text
    if TMP_flag_dark != FLAG_DARK:
        message = "At least 2 different values of FLAG_DARK for one processing: %s vs %s" % (FLAG_DARK, TMP_flag_dark)
        raise service_error.ProcessingError(message, logger)
            
    # Min size for lake product computation
    TMP_min_size = float(IN_xml_tree.xpath("//LakeTile_shp/config_params/min_size")[0].text)
    if TMP_min_size != MIN_SIZE:
        message = "At least 2 different values of MIN_SIZE for one processing: %s vs %s" % (MIN_SIZE, TMP_min_size)
        raise service_error.ProcessingError(message, logger)
    
    # Improve geolocation or not 
    TMP_imp_geoloc = bool(IN_xml_tree.xpath("//LakeTile_shp/config_params/imp_geoloc")[0].text)
    if TMP_imp_geoloc != IMP_GEOLOC:
        message = "At least 2 different values of IMP_GEOLOC for one processing: %s vs %s" % (IMP_GEOLOC, TMP_imp_geoloc)
        raise service_error.ProcessingError(message, logger)
            
    # Hull method
    TMP_hull_method = float(IN_xml_tree.xpath("//LakeTile_shp/config_params/hull_method")[0].text)
    if TMP_hull_method != HULL_METHOD:
        message = "At least 2 different values of HULL_METHOD for one processing: %s vs %s" % (HULL_METHOD, TMP_hull_method)
        raise service_error.ProcessingError(message, logger)
            
    # Maximal standard deviation of height inside a lake
    TMP_std_height_max = float(IN_xml_tree.xpath("//LakeTile_shp/config_params/std_height_max")[0].text)
    if TMP_std_height_max != STD_HEIGHT_MAX:
        message = "At least 2 different values of STD_HEIGHT_MAX for one processing: %s vs %s" % (STD_HEIGHT_MAX, TMP_std_height_max)
        raise service_error.ProcessingError(message, logger)
            
    # Model to deal with big lake processing
    TMP_biglake_model = IN_xml_tree.xpath("//LakeTile_shp/config_params/biglake_model")[0].text
    if TMP_biglake_model != BIGLAKE_MODEL:
        message = "At least 2 different values of BIGLAKE_MODEL for one processing: %s vs %s" % (BIGLAKE_MODEL, TMP_biglake_model)
        raise service_error.ProcessingError(message, logger)
        
    # Min size for lake to be considered as big
    TMP_biglake_min_size = float(IN_xml_tree.xpath("//LakeTile_shp/config_params/biglake_min_size")[0].text)
    if TMP_biglake_min_size != BIGLAKE_MIN_SIZE:
        message = "At least 2 different values of BIGLAKE_MIN_SIZE for one processing: %s vs %s" % (BIGLAKE_MIN_SIZE, TMP_biglake_min_size)
        raise service_error.ProcessingError(message, logger)
        
    # Grid spacing for lake height smoothing
    TMP_biglake_grid_spacing = int(IN_xml_tree.xpath("//LakeTile_shp/config_params/biglake_grid_spacing")[0].text)
    if TMP_biglake_grid_spacing != BIGLAKE_GRID_SPACING:
        message = "At least 2 different values of BIGLAKE_GRID_SPACING for one processing: %s vs %s" % (BIGLAKE_GRID_SPACING, TMP_biglake_grid_spacing)
        raise service_error.ProcessingError(message, logger)
        
    # Grid resolution for lake height smoothing
    TMP_biglake_grid_res = int(IN_xml_tree.xpath("//LakeTile_shp/config_params/biglake_grid_res")[0].text)
    if TMP_biglake_grid_res != BIGLAKE_GRID_RES:
        message = "At least 2 different values of BIGLAKE_GRID_RES for one processing: %s vs %s" % (BIGLAKE_GRID_RES, TMP_biglake_grid_res)
        raise service_error.ProcessingError(message, logger)
