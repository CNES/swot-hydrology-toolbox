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
.. module:: my_basins.py
    :synopsis: Deal with river basins shapefile.
     Created on 2018/11/09

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""

import logging
import cnes.common.service_config_file as service_config_file
from osgeo import ogr

def link_poly_to_continent(in_poly):
    """
    Link a polygon to a list of continent(s) by considering intersection of both
    
    :param in_poly: polygon to link to a continent
    :type in_poly: ogr.Polygon
    
    :return: list of continent(s) associated to polygon
    :rtype: list of string
    """
    # Get instance of service config file
    cfg = service_config_file.get_instance()
    logger = logging.getLogger("my_bassin")
    
    # Retrieve continent file
    continent_file = cfg.get("DATABASES", "CONTINENT_FILE")
    
    # Case when no continent file
    if continent_file is None:
        logger.debug("> No continent file")
        retour = None
        
    else:
        logger.debug("> Continent file = %s" % continent_file)
        
        # 1 - Open continent shapefile in read-only mode
        shp_driver = ogr.GetDriverByName(str("ESRI Shapefile"))
        data_source = shp_driver.Open(continent_file, 0)
        continent_layer = data_source.GetLayer()
        
        # 2 - Compute intersection
        continent_layer.SetSpatialFilter(in_poly)
        
        # 3 - Get continent name
        out_continent = []
        for item in continent_layer:
            out_continent.append(compute_continent_from_basin_id(str(item.GetField("MAJ_BAS"))))
        
        # 4 - Close continent file
        data_source.Destroy()
    
        # 5 - Return continent
        if len(out_continent) == 0:
            retour = "OCEAN"
        else:
            retour = out_continent[0]

    return retour


#######################################
    
    
def compute_continent_from_basin_id(in_basin_id):
    """
    Compute continent from basin ID (FAO nomenclature)
    
    :param in_basin_id: basin identifier
    :type in_basin_id: string
    
    :return: 2 letters related to the continent
    :rtype: string
    """
    
    retour = ""
    
    if in_basin_id.startswith("1"):
        retour = "NA"
    elif in_basin_id.startswith("2"):
        retour = "CA"
    elif in_basin_id.startswith("3"):
        retour = "SA"
    elif in_basin_id.startswith("4"):
        retour = "EU"
    elif in_basin_id.startswith("5"):
        retour = "EA"
    elif in_basin_id.startswith("6"):
        retour = "WA"
    elif in_basin_id.startswith("7"):
        retour = "AF"
    elif in_basin_id.startswith("8"):
        retour = "OC"
    elif in_basin_id.startswith("9"):
        retour = "AN"
    else:
        retour = "xx"
        
    return retour
