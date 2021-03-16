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


def stocc_linear_basic(in_list_obs, in_ref_area, in_ref_area_u, in_ref_wse, in_ref_wse_u):
    """
    Compute linear storage change between a state of a water body given by its (water surface elevation, area)
    and a reference given by its (wse, area).
    
    :param in_list_obs: list of observed features intersecting the prior lake
                            - in_list_obs[obs_id]["area"] = area of observed lake obs_id
                            - in_list_obs[obs_id]["area_u"] = uncertainty over area of observed lake obs_id
                            - in_list_obs[obs_id]["wse"] = water surface elevation of observed lake obs_id
                            - in_list_obs[obs_id]["wse_u"] = uncertainty over water surface elevation of observed lake obs_id
                            - in_list_obs[obs_id]["alpha"] = proprotionnal coefficient related to area of observed lake wrt
                                                                all observed lakes linked to the prior lake
    :type in_list_obs: dict
    :param in_ref_area: reference area used for storage change computation (in km2)
    :type in_ref_area: float
    :param in_ref_area_u: uncertainty over reference area used for storage change computation (in km2)
    :type in_ref_area_u: float
    :param in_ref_wse: reference water surface elevation used for storage change computation (in m)
    :type in_ref_wse: float
    :param in_ref_wse_u: uncertainty over reference water surface elevation used for storage change computation (in m)
    :type in_ref_wse_u: float
    
    :return out_stoc_val: linear storage change value (in km3)
    :rtype out_stoc_val: float
    :return out_stoc_u: linear storage change uncertainty (in km3)
    :rtype out_stoc_u: float
    """
    
    if in_ref_area_u is None:
        in_ref_area_u = 0.0
    if in_ref_wse_u is None:
        in_ref_wse_u = 0.0
    
    if (in_list_obs is None) or (in_ref_area is None) or (in_ref_wse is None):
        out_stoc_val = None
        out_stoc_u = None
        
    elif len(in_list_obs) == 1:
        
        obs_id = list(in_list_obs.keys())[0]
        
        # Volume variation between both surfaces
        out_stoc_val = (in_list_obs[obs_id]["wse"] - in_ref_wse)/2. * (in_ref_area + in_list_obs[obs_id]["area"])
        out_stoc_val /= 10**3  # Convert value in km3 (wse is in m instead of km)
            
        # Associated uncertainty
        dV_dhi = (in_ref_area + in_list_obs[obs_id]["area"])/2.
        dV_dhref = -dV_dhi
        dV_dAi = (in_list_obs[obs_id]["wse"] - in_ref_wse)/2.
        dV_dAref = dV_dAi
        out_stoc_u = math.sqrt( (dV_dhi*in_list_obs[obs_id]["wse_u"])**2 + (dV_dhref*in_ref_wse_u)**2 + \
                               (dV_dAi*in_list_obs[obs_id]["area_u"])**2 + (dV_dAref*in_ref_area_u)**2 )
        out_stoc_u /= 10**3  # Convert value in km3 (wse and wse_u are in m instead of km)
        
    else:
        
        out_stoc_val = 0
        tmp_stoc_u = 0
        for obs_id in in_list_obs.keys():
            
            # Volume variation between both surfaces
            out_stoc_val += (in_list_obs[obs_id]["wse"] - in_ref_wse)/2. * (in_list_obs[obs_id]["alpha"]*in_ref_area + in_list_obs[obs_id]["area"])
            
            # Associated uncertainty
            dV_dhi = (in_list_obs[obs_id]["alpha"]*in_ref_area + in_list_obs[obs_id]["area"])/2.
            dV_dhref = -dV_dhi
            dV_dAi = (in_list_obs[obs_id]["wse"] - in_ref_wse)/2.
            dV_dAref = in_list_obs[obs_id]["alpha"]*dV_dAi
            tmp_stoc_u += (dV_dhi*in_list_obs[obs_id]["wse_u"])**2 + (dV_dhref*in_ref_wse_u)**2 + \
                            (dV_dAi*in_list_obs[obs_id]["area_u"])**2 + (dV_dAref*in_ref_area_u)**2
            
        # Convert value in km3 (wse is in m instead of km)
        out_stoc_val /= 10**3
        
        # Compute uncertainty and convert in km3 (wse and wse_u are in m instead of km)
        out_stoc_u = math.sqrt(tmp_stoc_u) / 10**3
    
    return out_stoc_val, out_stoc_u


def stocc_quadratic_basic(in_list_obs, in_ref_area, in_ref_area_u, in_ref_wse, in_ref_wse_u):
    """
    Compute quadratic storage change between a state of a water body given by its (water surface elevation, area)
    and a reference given by its (wse, area).
    
    :param in_list_obs: list of observed features intersecting the prior lake
                            - in_list_obs[obs_id]["area"] = area of observed lake obs_id
                            - in_list_obs[obs_id]["area_u"] = uncertainty over area of observed lake obs_id
                            - in_list_obs[obs_id]["wse"] = water surface elevation of observed lake obs_id
                            - in_list_obs[obs_id]["wse_u"] = uncertainty over water surface elevation of observed lake obs_id
                            - in_list_obs[obs_id]["alpha"] = proprotionnal coefficient related to area of observed lake wrt
                                                                all observed lakes linked to the prior lake
    :type in_list_obs: dict
    :param in_ref_area: reference area used for storage change computation (in km2)
    :type in_ref_area: float
    :param in_ref_area_u: uncertainty over reference area used for storage change computation (in km2)
    :type in_ref_area_u: float
    :param in_ref_wse: reference water surface elevation used for storage change computation (in m)
    :type in_ref_wse: float
    :param in_ref_wse_u: uncertainty over reference water surface elevation used for storage change computation (in m)
    :type in_ref_wse_u: float
    
    :return out_stoc_val: quadratic storage change value (in km3)
    :rtype out_stoc_val: float
    :return out_stoc_u: quadratic storage change uncertainty (in km3)
    :rtype out_stoc_u: float
    """
    
    if in_ref_area_u is None:
        in_ref_area_u = 0.0
    if in_ref_wse_u is None:
        in_ref_wse_u = 0.0
    
    if (in_list_obs is None) or (in_ref_wse is None) or (in_ref_area is None):
        out_stoc_val = None
        out_stoc_u = None
        
    elif len(in_list_obs) == 1:
        
        obs_id = list(in_list_obs.keys())[0]
        
        # Volume variation between both surfaces
        out_stoc_val = (in_list_obs[obs_id]["wse"] - in_ref_wse)/3. * (in_ref_area + in_list_obs[obs_id]["area"] + math.sqrt(in_ref_area * in_list_obs[obs_id]["area"]))
        out_stoc_val /= 10**3  # Convert value in km3 (wse is in m instead of km)
            
        # Associated uncertainty
        dV_dhi = (in_ref_area + in_list_obs[obs_id]["area"] + math.sqrt(in_ref_area * in_list_obs[obs_id]["area"]))/3.
        dV_dhref = -dV_dhi
        dV_dAi = (in_list_obs[obs_id]["wse"] - in_ref_wse)/3. * ( 1. + math.sqrt(in_ref_area/in_list_obs[obs_id]["area"])/2. )
        dV_dAref = (in_list_obs[obs_id]["wse"] - in_ref_wse)/3. * ( 1. + math.sqrt(in_list_obs[obs_id]["area"]/in_ref_area)/2. )
        out_stoc_u = math.sqrt( (dV_dhi*in_list_obs[obs_id]["wse_u"])**2 + (dV_dhref*in_ref_wse_u)**2 + \
                               (dV_dAi*in_list_obs[obs_id]["area_u"])**2 + (dV_dAref*in_ref_area_u)**2 )
        out_stoc_u /= 10**3  # Convert value in km3 (wse and wse_u are in m instead of km)
        
    else:
        
        out_stoc_val = 0
        tmp_stoc_u = 0
        for obs_id in in_list_obs.keys():
            
            # Volume variation between both surfaces
            out_stoc_val += (in_list_obs[obs_id]["wse"] - in_ref_wse)/3. * (in_list_obs[obs_id]["alpha"]*in_ref_area + in_list_obs[obs_id]["area"]
                                + math.sqrt(in_list_obs[obs_id]["alpha"]*in_ref_area * in_list_obs[obs_id]["area"]))
            
            # Associated uncertainty
            dV_dhi = (in_list_obs[obs_id]["alpha"]*in_ref_area + in_list_obs[obs_id]["area"] + math.sqrt(in_list_obs[obs_id]["alpha"]*in_ref_area * in_list_obs[obs_id]["area"]))/3.
            dV_dhref = -dV_dhi
            dV_dAi = (in_list_obs[obs_id]["wse"] - in_ref_wse)/3. * ( 1. + math.sqrt(in_list_obs[obs_id]["alpha"]*in_ref_area/in_list_obs[obs_id]["area"])/2. )
            dV_dAref = (in_list_obs[obs_id]["wse"] - in_ref_wse)/3. * ( in_list_obs[obs_id]["alpha"] + math.sqrt(in_list_obs[obs_id]["area"]/(in_list_obs[obs_id]["alpha"]*in_ref_area))/2. )
            tmp_stoc_u += (dV_dhi*in_list_obs[obs_id]["wse_u"])**2 + (dV_dhref*in_ref_wse_u)**2 + \
                            (dV_dAi*in_list_obs[obs_id]["area_u"])**2 + (dV_dAref*in_ref_area_u)**2
            
        # Convert value in km3 (wse is in m instead of km)
        out_stoc_val /= 10**3
        
        # Compute uncertainty and convert in km3 (wse and wse_u are in m instead of km)
        out_stoc_u = math.sqrt(tmp_stoc_u) / 10**3
    
    return out_stoc_val, out_stoc_u
