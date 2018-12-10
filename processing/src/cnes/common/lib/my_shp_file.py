# -*- coding: utf8 -*-
"""
.. module:: my_shp_file.py
    :synopsis: Deal with shapefiles and memory layers
    Created on 09/14/2018

.. moduleauthor: Claire POTTIER (CNES DSO/SI/TR)

Copyright (c) 2018 CNES. All rights reserved.
"""
from __future__ import absolute_import, division, print_function, unicode_literals 

import os
from osgeo import ogr, osr

import cnes.common.lib.my_api as my_api


def mergeMemLayerWithShp(IN_listShp, IN_layer):
    """
    This function merges shapefiles listed in IN_listShp with the layer IN_layer (typically LakeTile shp with LakeSP memory layer).
    All polygons in all shapefiles are copied in the current layer. All fields are copied.

    :param IN_listShp: list of shapefiles full path
    :type IN_listShp: list of string
    :param IN_layer: layer in which to merge all objects
    :type IN_layer: OGRlayer
    
    :return OUT_dataSource: data source of output layer
    :rtype OUT_dataSource: OGRdataSource
    :return OUT_layer: output layer
    :rtype OUT_layer: OGRlayer
    """
    my_api.printDebug("[LakeProduct] == mergeMemLayerWithShp ==")

    # 1 - Create memory data source
    memDriver = ogr.GetDriverByName(str('MEMORY'))  # Memory driver
    OUT_dataSource = memDriver.CreateDataSource('memData')

    # 2 - Copy input layer to memory layer
    OUT_layer = OUT_dataSource.CopyLayer(IN_layer, str('tmp'))

    # 3 - For each input shapefile
    for curShp in IN_listShp:
        my_api.printDebug("[LakeProduct] > Adding %s" % os.path.basename(curShp))

        # 3.1 - Open the shapefile and get layer
        shpDriver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Shapefile driver
        cur_dataSource = shpDriver.Open(curShp, 0)  # Open in reading mode
        cur_layer = cur_dataSource.GetLayer()  # Get the layer

        # 3.2 - Merge layer into output layer
        OUT_dataSource, OUT_layer = merge2Layers(OUT_layer, cur_layer)

        # 3.3 - Close layer
        cur_dataSource.Destroy()

    return OUT_dataSource, OUT_layer


def merge2Layers(IN_layer1, IN_layer2):
    """
    Merge 2 memory layers
    
    :param IN_layer1: first layer
    :type IN_layer1: OGRlayer
    :param IN_layer2: second layer
    :type IN_layer2: OGRlayer
    
    :return OUT_dataSource: data source of output layer
    :rtype OUT_dataSource: OGRdataSource
    :return OUT_layer: output layer
    :rtype OUT_layer: OGRlayer
    """
    my_api.printDebug("[LakeProduct] == mergeLayerRL ==")

    # 1 - Get layer definitions
    layerDefn1 = IN_layer1.GetLayerDefn()
    layerDefn2 = IN_layer2.GetLayerDefn()

    # 2 - Check if layer definition of IN_layer1 and IN_layer2 are equal
    fields_name_1 = [layerDefn1.GetFieldDefn(i).GetName() for i in range(layerDefn1.GetFieldCount())]
    fields_name_2 = [layerDefn2.GetFieldDefn(i).GetName() for i in range(layerDefn2.GetFieldCount())]
    if fields_name_1 != fields_name_2:
        my_api.exitWithError("[LakeProduct] ERROR = fields of layer %s and layer %s are not identical" % (layerDefn1.GetName(), layerDefn2.GetName()))

    # 3 - Create output memory driver and data source
    memDriver = ogr.GetDriverByName(str('MEMORY'))  # Memory driver
    OUT_dataSource = memDriver.CreateDataSource('memData')

    # 4 - Projection
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)  # WGS84

    # 5 - Create output layer
    OUT_layer = OUT_dataSource.CreateLayer(str('tmp'), srs, geom_type=ogr.wkbMultiPolygon)

    # 6 - Init output layer fields
    for i in range(layerDefn1.GetFieldCount()):
        idField = layerDefn1.GetFieldDefn(i)
        OUT_layer.CreateField(idField)

    # 7 - Add IN_layer1 and IN_layer2 to the output layer
    IN_layer1.Update(IN_layer2, OUT_layer)

    return OUT_dataSource, OUT_layer


def writeMemLayer_asShp(IN_mem_layer, IN_shp_filename):
    """
    Write memory layer IN_mem_layer into a shapefile
    
    :param IN_mem_layer: memory layer
    :type IN_mem_layer: ogr.Layer
    :param IN_shp_filename: shapefile full path
    :type IN_shp_filename: string
    """
    
    shpDriver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Driver for shapefiles
    
    # 1 - Delete output file if already exists
    if os.path.exists(IN_shp_filename):
        shpDriver.DeleteDataSource(IN_shp_filename)
        
    # 2 - Create output file
    dataSource = shpDriver.CreateDataSource(IN_shp_filename)
    
    # 3 - Copy memory layer to output 
    layer_name = os.path.basename(IN_shp_filename).split(".")[0]
    dataSource.CopyLayer(IN_mem_layer, str(layer_name))
    
    # 4 - Close output file
    dataSource.Destroy()
