# -*- coding: utf-8 -*-
"""
.. module my_variables.py
    :synopsis: Gather global variables used in SISIMP

.. module author: CNES DSO/SI/TR

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from math import pi

# Earth parameters
GEN_RAD_EARTH_EQ = 6378137.0  # Radius of the Earth model (WGS84 ellipsoid) at the equator
GEN_RAD_EARTH_POLE = 6356752.31425  # Radius of the Earth model at the pole
GEN_APPROX_RAD_EARTH = (2*GEN_RAD_EARTH_EQ + GEN_RAD_EARTH_POLE)/3  # Radius (in meters) of the sphere equivalent to ellipsoid

# Instrument parameters
SWATH_WIDTH = 120000.000000  # Swath width (m)
NR_CROSS_TRACK = 5000.000  # NR cross track (m)
SENSOR_WAVELENGTH = 0.008385803  # Sensor wavelength (m)
NB_PIX_RANGE = 3117  # Number of pixels in range 3117 (10-60km) or 3500 (extended swath 5-65km)
RANGE_SAMPLING = 0.75  # Range spacing of SWOT
BASELINE = 10.

# Noise parameters
GEOLOCATION_IMPROVEMENT = 'no'
NOISE_MULTIPLIER_FACTOR = 0.5  # Noise multiplier factor 1/sqrt(Nl) where Nl=4(nb_multilook)
HEIGHT_BIAS_STD = 0.1  # Height bias std (m)

NOISE_FILE_PATH = "$SWOT_HYDROLOGY_TOOLBOX/sisimp/data/height_noise_presum2.txt"
NOISE_FILE_PATH_FOR_LAND = "$SWOT_HYDROLOGY_TOOLBOX/sisimp/data/height_noise_presum2_land.txt"
NOISE_FILE_PATH_FOR_DW = "$SWOT_HYDROLOGY_TOOLBOX/sisimp/data/height_noise_presum2_land.txt"
TILE_DATABASE_PATH = "$SWOT_HYDROLOGY_TOOLBOX/sisimp/data/tiles_full.txt.zip"
GEOID_PATH = "$SWOT_HYDROLOGY_TOOLBOX/sisimp/data/egm2008-5.pgm"

# Orbit parameters
ORBIT_JITTER = 1000  # Orbit jitter (m)
MULTIPLE_ORBIT = 'yes'

# Height model
HEIGHT_MODEL = None   # polynomial / gaussian (default)
HEIGHT_MODEL_MIN_AREA = 100  # Optionnal argument to add complex 2D height model

# Constant height model (always applied, set HEIGHT_MODEL_A to 0 to desactivate)
HEIGHT_MODEL_A = 10  # Height model A (m)
HEIGHT_MODEL_t0 = 47076  # Height model t0
HEIGHT_MODEL_PERIOD = 365.25  # Height model period (days)

# Polynomial parameters for polynomial model for big lake
COEFF_X2 = 1.e-9
COEFF_Y2 = 1.e-9
COEFF_X = 5.e-5
COEFF_Y = 5.e-5
COEFF_XY = 1.e-9
COEFF_CST = 0.

# Dark water
FACT_ECHELLE_DW = 2.
DW_PERCENT = 10.
DARKWATER_FLAG = 24
SCALE_FACTOR_NON_DETECTED_DW = 0.5
DW_DETECTED_PERCENT = 90
DW_DETECTED_NOISE_FACTOR = 5
DW_CORRELATION_LENGTH = 50

# Gaussian parameter for gaussian model
FACT_ECHELLE = 2.
HEIGHT_MODEL_STDV = 0.1

# Water flag
WATER_FLAG = 4  # Water flag
WATER_LAND_FLAG = 3
LAND_FLAG = 1
LAND_WATER_FLAG = 2
LAND_DETECTED_NOISE_FACTOR = 10.

# Degrees / radians convertors
RAD2DEG = 180. / pi
DEG2RAD = pi / 180.

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
FV_NETCDF = {'int8': FV_BYTE,
              'int16': FV_SHORT,
              'int32': FV_INT,
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

# Filenames pattern
PATTERN_FOOTPRINT = "footprint_%03d_%03d.shp"  # Footprint filename with %03d=cycle number %03d=pass number
PATTERN_PIXC = "SWOT_L2_HR_PIXC_%03d_%03d_%s_%s_%s_Dx0000_01"  # Pixel cloud filename with %03d=cycle number %03d=pass number %s=tile ref %s=begin date %s=end date
PATTERN_FILE_ANNOT = "pixc_annotation_%03d_%03d_%s.rdf"  # PixC annotation filename with %03d=cycle number %03d=pass number %s=tile ref
PATTERN_PIXC_VEC_RIVER = "SWOT_L2_HR_PIXCVecRiver_%03d_%03d_%s_%s_%s_Dx0000_01"  # PIXCVecRiver filename with %03d=cycle number %03d=pass number %s=tile ref %s=begin date %s=end date
