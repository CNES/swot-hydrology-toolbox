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

    
def linkPolyToContinent(IN_poly):
    """
    Link a polygon to a list of continent(s) by considering intersection of both
    
    :param IN_poly: polygon to link to a continent
    :type IN_poly: ogr.Polygon
    
    :return: list of continent(s) associated to polygon
    :rtype: list of string
    """
    my_api.printDebug("[my_tools] == linkPolyToContinent ==")
    
    # Case when no continent file
    if my_var.CONTINENT_FILE is None:
        my_api.printDebug("> No continent file")
        return None
    
    my_api.printDebug("> Continent file = %s" % my_var.CONTINENT_FILE)
    
    # 1 - Open continent shapefile in read-only mode
    shpDriver = ogr.GetDriverByName(str("ESRI Shapefile"))
    dataSource = shpDriver.Open(my_var.CONTINENT_FILE, 0)
    continent_layer = dataSource.GetLayer()
    
    # 2 - Compute intersection
    continent_layer.SetSpatialFilter(IN_poly)
    
    # 3 - Get continent name
    OUT_continent = []
    for item in continent_layer:
        OUT_continent.append(computeContinentFromBasinId(str(item.GetField("MAJ_BAS"))))
    
    # 4 - Close continent file
    dataSource.Destroy()
    
    # 5 - Return continent
    return OUT_continent[0]


#######################################
    
    
def computeContinentFromBasinId(IN_basin_id):
    """
    Compute continent from basin ID (FAO nomenclature)
    
    :param IN_basin_id: basin identifier
    :type IN_basin_id: string
    
    :return: 2 letters related to the continent
    :rtype: string
    """
    
    if IN_basin_id.startswith("1"):
        return "NA"
    elif IN_basin_id.startswith("2"):
        return "CA"
    elif IN_basin_id.startswith("3"):
        return "SA"
    elif IN_basin_id.startswith("4"):
        return "EU"
    elif IN_basin_id.startswith("5"):
        return "EA"
    elif IN_basin_id.startswith("6"):
        return "WA"
    elif IN_basin_id.startswith("7"):
        return "AF"
    elif IN_basin_id.startswith("8"):
        return "OC"
    elif IN_basin_id.startswith("9"):
        return "AN"
    else:
        return "xx"
