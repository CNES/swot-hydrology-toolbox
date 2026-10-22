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
# VERSION:4.0.0:DM:#91:2022/05/05:Poursuite industrialisation
# VERSION:5.0.0:DM:#91:2022/12/08:Poursuite industrialisation
# VERSION:5.0.1:DM:#91:2023/02/24:Poursuite industrialisation
# VERSION:5.2.1:DM:#91:2023/07/24:Poursuite industrialisation
# VERSION:5.3.0:DM:#91:2023/10/04:Poursuite industrialisation
# VERSION:6.1.0:DM:#91:2024/09/20:Poursuite industrialisation
# VERSION:6.2.0:DM:#91:2024/11/29:Poursuite industrialisation
# VERSION:6.2.2:DM:#822:2025/06/30:Nouveau patch Lake_TILE_V6.2.1
# VERSION:6.2.3:FA:#839:2025/11/10:[LAKE_TILE] L2_HR_LAKETILE fails on C032P429T305R with message ValueError: list.remove(x): x not in list
# VERSION:6.2.6:DM:#900:2026/07/06: [PGE] [LAKES] Ajout d'une interface niveau moyen des lacs
# VERSION:6.2.6:DM:#900:2026/07/06:SSALTO-565 [PGE] [LAKES] Ajout d'une interface coeff NN pour calcul ice_dyn_f
# VERSION:6.2.7:DM:#906:2026/07/28:SSALTO-872 [PGE][LAKES] Gestion des erreurs/warning GDAL
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: my_variables.py
    :synopsis: Gather generic variables
     Created on 2018/03/08

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from osgeo import ogr


######################
# GENERIC PARAMETERS #
######################

# Earth parameters
GEN_RAD_EARTH_EQ = 6378137.0  # Radius of the Earth model (WGS84 ellipsoid) at the equator
GEN_RAD_EARTH_POLE = 6356752.31425  # Radius of the Earth model at the pole
GEN_APPROX_RAD_EARTH = (2*GEN_RAD_EARTH_EQ + GEN_RAD_EARTH_POLE)/3  # Radius (in meters) of the sphere equivalent to ellipsoid

# SWOT parameters
GEN_RANGE_SPACING = 0.75  # Range spacing of SWOT

# Compression parameter
GEN_NETCDF_COMPLEVEL = 4 # Level of compression use for netcdf file

###############
# _FillValues #
###############

# FillValues for NetCDF files 
# double based on netcdf.h value
FV_DOUBLE = 9.9692099683868690e+36
# float based on value in netcdf.h file
FV_FLOAT = 9.96921e+36
# int64 based on value in netcdf.h file
FV_INT64 = 9223372036854775807
# int inspired by value in netcdf.h file (positive instead of negative)
FV_INT = 2147483647
# uint based on value in netcdf.h file
FV_UINT = 4294967295
# short inspired by value in netcdf.h file (positive instead of negative)
FV_SHORT = 32767
# ushort based on value in netcdf.h file
FV_USHORT = 65535
# byte inspired by value in netcdf.h file (positive instead of negative)
FV_BYTE = 127
# ubyte based on value in netcdf.h file
FV_UBYTE = 255
# no fillvalue for char
FV_CHAR = ""
# no fillvalue for string
FV_STRING = ""
FV_NETCDF = {'int': FV_INT,
             'integer': FV_INT, 
             'int4': FV_BYTE, 
             'byte': FV_BYTE,
             'bytes': FV_BYTE,
             'int8': FV_BYTE,
             'short': FV_SHORT,
             'int16': FV_SHORT,
             'int32': FV_INT,
             'int64': FV_INT64,
             'uint8': FV_UBYTE,
             'uint16': FV_USHORT,
             'uint32': FV_UINT,
             'float': FV_FLOAT,
             'float32': FV_FLOAT,
             'double': FV_DOUBLE,
             'float64': FV_DOUBLE,
             'char': FV_CHAR,
             'str': FV_STRING,
             'str32': FV_STRING,
             'object': FV_STRING}

# FillValues for Shapefile
FV_REAL = -999999999999
FV_INT9_SHP = -99999999
FV_INT_SHP = -999
FV_STRING_SHP = "no_data"
FV_SHP = {'int': FV_INT_SHP,  
          'int4': FV_INT_SHP, 
          'int8': FV_INT_SHP, 
          'int9': FV_INT9_SHP,
          'int16': FV_INT_SHP,
          'int32': FV_INT_SHP,
          'integer': FV_INT_SHP,
          'uint8': FV_INT_SHP,
          'uint16': FV_INT_SHP,
          'uint32': FV_INT_SHP,
          'float': FV_REAL,
          'float32': FV_REAL,
          'double': FV_REAL,
          'float64': FV_REAL,
          'real': FV_REAL,
          'str': FV_STRING_SHP,
          'text': FV_STRING_SHP,
          'str32': FV_STRING_SHP,
          'string': FV_STRING_SHP, 
          'object': FV_STRING_SHP}


########################################
# Conversion metadata type to OGR type #
########################################

FORMAT_OGR = {'int': ogr.OFTInteger,
              'int4': ogr.OFTInteger,
              'int9': ogr.OFTInteger,
              'integer': ogr.OFTInteger,
              'float': ogr.OFTReal,
              'real': ogr.OFTReal,
              'text': ogr.OFTString,
              'string': ogr.OFTString}

FORMAT_OGR_STR = {'integer': "ogr.OFTInteger",
                  'int': "ogr.OFTInteger",
                  'int4': "ogr.OFTInteger",
                  'int9': "ogr.OFTInteger",
                  'float': "ogr.OFTReal",
                  'real': "ogr.OFTReal",
                  'text': "ogr.OFTString",
                  'string': "ogr.OFTString"}


########################################
# Conversion metadata type to SQL type #
########################################

FORMAT_SQL_STR = {'integer': "INTEGER",
                  'int': "INTEGER",
                  'int4': "INTEGER",
                  'int9': "INTEGER",
                  'float': "DOUBLE",
                  'real': "DOUBLE",
                  'text': "TEXT",
                  'string': "TEXT"}


#################################
# Prior Lake Database structure #
# Operational format in SQLite  #
#################################

# Table names
PLD_TABLE_LAKE = "lake"
PLD_TABLE_BASIN = "basin"
PLD_TABLE_LAKE_INFL = "lake_influence"
PLD_TABLE_LAKE_CATCH = "lake_catchment"
PLD_TABLE_PASS_TILE = "pass_tile"

# Fields names
# lake table
PLD_FIELD_LAKE = dict()  # Key = general key
PLD_FIELD_LAKE["id"] = "lake_id"
PLD_FIELD_LAKE["names"] = "names"
PLD_FIELD_LAKE["res_id"] = "res_id"
PLD_FIELD_LAKE["reach"] = "reach_id_list"
PLD_FIELD_LAKE["lon"] = "lon"
PLD_FIELD_LAKE["lat"] = "lat"
PLD_FIELD_LAKE["ice_clim"] = "ice_clim_flag"
PLD_FIELD_LAKE["ice_start"] = "ice_start"
PLD_FIELD_LAKE["ice_end"] = "ice_end"
PLD_FIELD_LAKE["max_wse"] = "ref_wse"
PLD_FIELD_LAKE["max_wse_u"] = "ref_wse_u"
PLD_FIELD_LAKE["max_area"] = "ref_area"
PLD_FIELD_LAKE["max_area_u"] = "ref_area_u"
PLD_FIELD_LAKE["ref_date"] = "date_t0"
PLD_FIELD_LAKE["ref_ds"] = "ds_t0"
PLD_FIELD_LAKE["storage"] = "storage"
PLD_FIELD_LAKE["overlap_avg_l"] = "overlap_avg_l"
PLD_FIELD_LAKE["nb_tiles_l"] = "nb_tiles_l"
PLD_FIELD_LAKE["overlap_avg_r"] = "overlap_avg_r"
PLD_FIELD_LAKE["nb_tiles_r"] = "nb_tiles_r"
PLD_FIELD_LAKE["is_mls"] = "is_mls"
PLD_FIELD_LAKE_TYPE = dict()  # Type of previous fields; key = PLD fieldname
PLD_FIELD_LAKE_TYPE["lake_id"] = "text"
PLD_FIELD_LAKE_TYPE["names"] = "text"
PLD_FIELD_LAKE_TYPE["res_id"] = "int9"
PLD_FIELD_LAKE_TYPE["reach_id_list"] = "text"
PLD_FIELD_LAKE_TYPE["lon"] = "float"
PLD_FIELD_LAKE_TYPE["lat"] = "float"
PLD_FIELD_LAKE_TYPE["ice_clim_flag"] = "int4"
PLD_FIELD_LAKE_TYPE["ice_start"] = "text"
PLD_FIELD_LAKE_TYPE["ice_end"] = "text"
PLD_FIELD_LAKE_TYPE["ref_wse"] = "float"
PLD_FIELD_LAKE_TYPE["ref_wse_u"] = "float"
PLD_FIELD_LAKE_TYPE["ref_area"] = "float"
PLD_FIELD_LAKE_TYPE["ref_area_u"] = "float"
PLD_FIELD_LAKE_TYPE["date_t0"] = "text"
PLD_FIELD_LAKE_TYPE["ds_t0"] = "float"
PLD_FIELD_LAKE_TYPE["storage"] = "float"
PLD_FIELD_LAKE_TYPE["overlap_avg_l"] = "float"
PLD_FIELD_LAKE_TYPE["nb_tiles_l"] = "int"
PLD_FIELD_LAKE_TYPE["overlap_avg_r"] = "float"
PLD_FIELD_LAKE_TYPE["nb_tiles_r"] = "int"
PLD_FIELD_LAKE_TYPE["is_mls"] = "int"
# basin table
PLD_FIELD_BASIN = dict()  # Key = general key
PLD_FIELD_BASIN["id"] = "basin_id"
PLD_FIELD_BASIN["lon_min"] = "lon_min"
PLD_FIELD_BASIN["lon_max"] = "lon_max"
PLD_FIELD_BASIN["lat_min"] = "lat_min"
PLD_FIELD_BASIN["lat_max"] = "lat_max"
# pass_tile table
PLD_FIELD_TILE = dict()  # Key = general key
PLD_FIELD_TILE["id"] = "ref_tile"
PLD_FIELD_TILE["cont"] = "continent"
PLD_FIELD_TILE["az_0_line"] = "az_0_line"
PLD_FIELD_TILE["az_max_line"] = "az_max_line"

# Mapping between prior and PLD fieldnames
MAPPING_PRIOR_ATT_WRT_PLD = dict()
MAPPING_PRIOR_ATT_WRT_PLD["lake_id"] = PLD_FIELD_LAKE["id"]
MAPPING_PRIOR_ATT_WRT_PLD["lake_name"] = PLD_FIELD_LAKE["names"]
MAPPING_PRIOR_ATT_WRT_PLD["p_res_id"] = PLD_FIELD_LAKE["res_id"]
MAPPING_PRIOR_ATT_WRT_PLD["reach_id"] = PLD_FIELD_LAKE["reach"]
MAPPING_PRIOR_ATT_WRT_PLD["p_lon"] = PLD_FIELD_LAKE["lon"]
MAPPING_PRIOR_ATT_WRT_PLD["p_lat"] = PLD_FIELD_LAKE["lat"]
MAPPING_PRIOR_ATT_WRT_PLD["ice_clim_f"] = PLD_FIELD_LAKE["ice_clim"]
MAPPING_PRIOR_ATT_WRT_PLD["p_ref_wse"] = PLD_FIELD_LAKE["max_wse"]
MAPPING_PRIOR_ATT_WRT_PLD["p_ref_area"] = PLD_FIELD_LAKE["max_area"]
MAPPING_PRIOR_ATT_WRT_PLD["p_date_t0"] = PLD_FIELD_LAKE["ref_date"]
MAPPING_PRIOR_ATT_WRT_PLD["p_ds_t0"] = PLD_FIELD_LAKE["ref_ds"]
MAPPING_PRIOR_ATT_WRT_PLD["p_storage"] = PLD_FIELD_LAKE["storage"]

# Prior fields to keep in _Obs shapefile
PLD_FIELD_TO_KEEP_IN_OBS = ["lake_name", "p_res_id", "ice_clim_f", "reach_id"]


#############################
# PIXC classification flags #
#############################

# Water flags
CLASSIF_LAND_EDGE = [2]
CLASSIF_WATER_EDGE = [3, 6]
CLASSIF_INTERIOR_WATER = [4, 7]
CLASSIF_OPEN_WATER = [4]

# Dark water flags
CLASSIF_DARK = [5]

# Low-coherence flags
CLASSIF_LOW = [6, 7]


########
# MISC #
########

SEP_FILES = ", "  # Separator for list of files
SEP_ATT = ";"     # Separator for list of values of an attribute
BASIN_4_OCEAN = "000"  # 3-digits identifier for coastal zones
BASIN_DEFAULT = "010"  # 3-digits identifier for default basin
