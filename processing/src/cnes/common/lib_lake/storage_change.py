#!/usr/bin/env python
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
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: storage_change.py
    :synopsis: STOCC module = deal with storage change computation wrt to a reference
     Created on 2018/12/19

.. moduleauthor:: Manon QUELLEC (LEGOS) & Claire POTTIER (CNES DSO/SI/TR)

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import math

import cnes.common.lib.my_variables as my_var


def stocc_linear(in_height, in_area, in_ref_height, in_ref_area, in_ref_flood_dem=None):
    """
    Compute linear storage change between a state of a water body given by its (height, area)
    and a reference given by its (height, area), eventually flooded DEM (ie several levels of area/height/polygon).
    
    :param in_height: height of water body to compare (in m)
    :type in_height: float 
    :param in_area: area of water body to compare (in km2)
    :type in_area: float
    :param in_ref_height: reference height used to compute storage change (in m)
    :type in_ref_height: float
    :param in_ref_area: reference area used to compute storage change (in km2)
    :type in_ref_area: float
    :param in_ref_flood_dem: 
    :type in_ref_flood_dem: TBD (dict with height/area/polygon?)
    
    :return out_stoc_val: linear storage change value (in km3)
    :rtype out_stoc_val: float
    :return out_stoc_u: linear storage change uncertainty
    :rtype out_stoc_u: float
    """
    
    if (in_height is None) or (in_area is None) or (in_ref_height is None) or (in_ref_area is None):
        retour = None, None
        
    else:
        # Volume variation between both surfaces
        out_stoc_val = (in_height - in_ref_height)/2. * (in_area + in_ref_area)
        out_stoc_val /= 10**3  # Conversion in km3 (height is in m instead of km)
        
        # Associated uncertainty
        out_stoc_u = my_var.FV_REAL
            
        # Return
        retour = out_stoc_val, out_stoc_u
    
    return retour


def stocc_quadratic(in_height, in_area, in_ref_height, in_ref_area, in_ref_flood_dem=None):
    """
    Compute quadratic storage change between a state of a water body given by its (height, area)
    and a reference given by its (height, area), eventually flooded DEM (ie several levels of area/height/polygon).
    
    :param in_height: height of water body to compare (in m)
    :type in_height: float 
    :param in_area: area of water body to compare (in km2)
    :type in_area: float
    :param in_ref_height: reference height used to compute storage change (in m)
    :type in_ref_height: float
    :param in_ref_area: reference area used to compute storage change (in km2)
    :type in_ref_area: float
    :param in_ref_flood_dem: 
    :type in_ref_flood_dem: TBD (dict with height/area/polygon?)
    
    :return out_stoc_val: quadratic storage change value (in km3)
    :rtype out_stoc_val: float
    :return out_stoc_u: quadratic storage change uncertainty
    :rtype out_stoc_u: float
    """
    
    if (in_height is None) or (in_area is None) or (in_ref_height is None) or (in_ref_area is None):
        retour = None, None
        
    else:
        # Volume variation between both surfaces
        out_stoc_val = (in_height - in_ref_height)/3. * (in_area + in_ref_area + math.sqrt(in_area * in_ref_area))
        out_stoc_val /= 10**3  # Conversion in km3 (height is in m instead of km)
        
        # Associated uncertainty
        out_stoc_u = my_var.FV_REAL
            
        # Return
        retour = out_stoc_val, out_stoc_u
        
    return retour
