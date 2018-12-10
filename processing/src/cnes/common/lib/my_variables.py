# -*- coding: utf8 -*-
"""
.. module:: my_variables.py
    :synopsis: Gather generic variables
    Created on 03/08/2018

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

Copyright (c) 2018 CNES. All rights reserved.
"""

# Earth parameters
GEN_RAD_EARTH_EQ = 6378137.0  # Radius of the Earth model (WGS84 ellipsoid) at the equator
GEN_RAD_EARTH_POLE = 6356752.31425  # Radius of the Earth model at the pole
GEN_APPROX_RAD_EARTH = (2*GEN_RAD_EARTH_EQ + GEN_RAD_EARTH_POLE)/3  # Radius (in meters) of the sphere equivalent to ellipsoid

# SWOT parameters
GEN_RANGE_SPACING = 0.75  # Range spacing of SWOT
