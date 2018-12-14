# -*- coding: utf-8 -*-
"""
.. module:: my_basins.py
    :synopsis: Deal with river basins shapefile.
    Created on 09/11/2018

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

Copyright (c) 2018 CNES. All rights reserved.
"""

from osgeo import ogr

import cnes.common.lib.my_api as my_api
import cnes.common.lib_lake.locnes_variables as my_var

    
def link_poly_to_continent(in_poly):
    """
    Link a polygon to a list of continent(s) by considering intersection of both
    
    :param in_poly: polygon to link to a continent
    :type in_poly: ogr.Polygon
    
    :return: list of continent(s) associated to polygon
    :rtype: list of string
    """
    my_api.printDebug("[my_tools] == link_poly_to_continent ==")
    
    # Case when no continent file
    if my_var.CONTINENT_FILE is None:
        my_api.printDebug("> No continent file")
        retour = None
    else:
        my_api.printDebug("> Continent file = %s" % my_var.CONTINENT_FILE)
        
        # 1 - Open continent shapefile in read-only mode
        shp_driver = ogr.GetDriverByName(str("ESRI Shapefile"))
        data_source = shp_driver.Open(my_var.CONTINENT_FILE, 0)
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