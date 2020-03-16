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


###############
# _FillValues #
###############

# FillValues for NetCDF files 
FV_DOUBLE = 9.9692099683868690e+36
FV_FLOAT = 9.96921e+36
FV_INT = 2147483647
FV_UINT = 4294967295
FV_SHORT = 32767
FV_USHORT = 65535
FV_BYTE = 127
FV_UBYTE = 255
FV_CHAR = ""
FV_STRING = ""
FV_NETCDF = {'int': FV_INT,
              'int8': FV_BYTE,
              'int16': FV_SHORT,
              'int32': FV_INT,
              'int64': FV_INT,
              'uint8': FV_UBYTE,
              'uint16': FV_USHORT,
              'uint32': FV_UINT,
              'float': FV_FLOAT,
              'float32': FV_FLOAT,
              'double': FV_DOUBLE,
              'float64': FV_DOUBLE,
              'str': FV_STRING,
              'str32': FV_STRING,
              'object': FV_STRING}

# FillValues for Shapefile
FV_REAL = -999999999999
FV_INT9_SHP = -99999999
FV_INT_SHP = -999
FV_STRING_SHP = "no_data"
FV_SHP = {'int4': FV_INT_SHP, 
          'int8': FV_INT_SHP, 
          'int9': FV_INT9_SHP,
          'int16': FV_INT_SHP,
          'int32': FV_INT_SHP,
          'uint8': FV_INT_SHP,
          'uint16': FV_INT_SHP,
          'uint32': FV_INT_SHP,
          'float': FV_REAL,
          'float32': FV_REAL,
          'double': FV_REAL,
          'float64': FV_REAL,
          'str': FV_STRING_SHP,
          'text': FV_STRING_SHP,
          'str32': FV_STRING_SHP,
          'object': FV_STRING_SHP}


########################################
# Conversion metadata type to OGR type #
########################################

FORMAT_OGR = {'int4': ogr.OFTInteger,
              'int9': ogr.OFTInteger,
              'float': ogr.OFTReal,
              'text': ogr.OFTString}


#################################
# Prior Lake Database structure #
# Operational format in SQLite  #
#################################

# Table names
PLD_TABLE_LAKE = "lake"
PLD_TABLE_BASIN = "basin"
PLD_TABLE_LAKE_INFL = "lake_influence"

# Fields names
PLD_FIELD_LAKE_ID = "lakedb_id"
PLD_FIELD_BASIN_ID = "basin_id"
PLD_FIELD_REF_AREA = "ref_area"
PLD_FIELD_REF_HEIGHT = "ref_height"
PLD_FIELD_NAMES = "names"
PLD_FIELD_GRAND_ID = "grand_id"
