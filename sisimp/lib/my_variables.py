# -*- coding: utf8 -*-
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
NR_CROSS_TRACK = 10000.000  # NR cross track (m)
SENSOR_WAVELENGTH = 0.008385803  # Sensor wavelength (m)
NB_PIX_RANGE = 3117  # Number of pixels in range 3117 (10-60km) or 3500 (extended swath 5-65km)
RANGE_SAMPLING = 0.75  # Range spacing of SWOT

# Noise parameters
GEOLOCATION_IMPROVEMENT = 'no'
NOISE_MULTIPLIER_FACTOR = 0.5  # Noise multiplier factor 1/sqrt(Nl) where Nl=4(nb_multilook)
HEIGHT_BIAS_STD = 0.1  # Height bias std (m)

# Orbit parameters
ORBIT_JITTER = 1000  # Orbit jitter (m)
MULTIPLE_ORBIT = 'yes'

# Height model
HEIGHT_MODEL = 'gaussian'   # polynomial / gaussian (default)
HEIGHT_MODEL_MIN_AREA = 100  # Optionnal argument to add complex 2D height model

# Constant height model (always applied, set HEIGHT_MODEL_A to 0 to desactivate)
HEIGHT_MODEL_A = 10  # Height model A (m)
HEIGHT_MODEL_t0 = 47076  # Height model t0
HEIGHT_MODEL_PERIOD = 365.25  # Height model period (days)

# Polynomial parameters for polynomial model for big lake
COEFF_X2 = 0.
COEFF_Y2 = 0.
COEFF_X = 5.e-5
COEFF_Y = 5.e-5
COEFF_XY = 0.
COEFF_CST = 0.

# Gaussian parameter for gaussian model
FACT_ECHELLE = 2.
HEIGHT_MODEL_STDV = 0.1

# Dark water
FACT_ECHELLE_DW = 2.0  # Float to parameterize the correlation length of simulated dark water areas
DW_PERCENT = 10  # Probability of total dark water simulated (int between 0 and 100)
DARKWATER_FLAG = 24  # Classification flag for detected dark water (usually 24)

## Error models

# Water flag
WATER_FLAG = 4  # Water flag

# Degrees / radians convertors
RAD2DEG = 180. / pi
DEG2RAD = pi / 180.

# Filenames pattern
PATTERN_FOOTPRINT = "footprint_%03d_%03d.shp"  # Footprint filename with %03d=cycle number %03d=pass number
PATTERN_PIXC = "SWOT_L2_HR_PIXC_%03d_%03d_%s_%s_%s_Dx0000_01"  # Pixel cloud filename with %03d=cycle number %03d=pass number %s=tile ref %s=begin date %s=end date
PATTERN_FILE_ANNOT = "pixc_annotation_%03d_%03d_%s.rdf"  # PixC annotation filename with %03d=cycle number %03d=pass number %s=tile ref
PATTERN_PIXC_VEC_RIVER = "SWOT_L2_HR_PIXCVecRiver_%03d_%03d_%s_%s_%s_Dx0000_01"  # PIXCVecRiver filename with %03d=cycle number %03d=pass number %s=tile ref %s=begin date %s=end date
