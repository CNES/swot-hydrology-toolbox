#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
.. module:: storage_change.py
    :synopsis: STOCC module = deal with storage change computation wrt to a reference
    Created on 12/19/2019

.. moduleauthor:: Manon QUELLEC (LEGOS) & Claire POTTIER (CNES DSO/SI/TR)

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import math


def STOCC_linear(in_height, in_area, in_ref_height, in_ref_area, in_ref_flood_dem=None):
    """
    Compute linear storage change between a state of a water body given by its (height, area)
    and a reference given by its (height, area), eventually flooded DEM (ie several levels of area/height/polygon).
    
    :param in_height: height of water body to compare
    :type in_height: float 
    :param in_area: area of water body to compare (in ha)
    :type in_area: float
    :param in_ref_height: reference height used to compute storage change
    :type in_ref_height: float
    :param in_ref_area: reference area used to compute storage change (in ha)
    :type in_ref_area: float
    :param in_ref_flood_dem: 
    :type in_ref_flood_dem: TBD (dict with height/area/polygon?)
    
    :return out_stoc_val: linear storage change value
    :rtype out_stoc_val: float
    :return out_stoc_u: linear storage change uncertainty
    :rtype out_stoc_u: float
    """
    
    if (in_height is None) or (in_area is None) or (in_ref_height is None) or (in_ref_area is None):
        return None, None
    
    # Volume variation between both surfaces
    out_stoc_val = (in_height - in_ref_height)/2. * (in_area + in_ref_area)
    
    # Associated uncertainty
    out_stoc_u = -9999.
        
    # Return
    return out_stoc_val, out_stoc_u


def STOCC_quadratic(in_height, in_area, in_ref_height, in_ref_area, in_ref_flood_dem=None):
    """
    Compute quadratic storage change between a state of a water body given by its (height, area)
    and a reference given by its (height, area), eventually flooded DEM (ie several levels of area/height/polygon).
    
    :param in_height: height of water body to compare
    :type in_height: float 
    :param in_area: area of water body to compare (in ha)
    :type in_area: float
    :param in_ref_height: reference height used to compute storage change
    :type in_ref_height: float
    :param in_ref_area: reference area used to compute storage change (in ha)
    :type in_ref_area: float
    :param in_ref_flood_dem: 
    :type in_ref_flood_dem: TBD (dict with height/area/polygon?)
    
    :return out_stoc_val: quadratic storage change value
    :rtype out_stoc_val: float
    :return out_stoc_u: quadratic storage change uncertainty
    :rtype out_stoc_u: float
    """
    
    if (in_height is None) or (in_area is None) or (in_ref_height is None) or (in_ref_area is None):
        return None, None
    
    # Volume variation between both surfaces
    out_stoc_val = (in_height - in_ref_height)/3. * (in_area + in_ref_area + math.sqrt(in_area * in_ref_area))
    
    # Associated uncertainty
    out_stoc_u = -9999.
        
    # Return
    return out_stoc_val, out_stoc_u