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
# VERSION:3.0.0:DM:#91:2021/03/12:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: my_shp_file.py
    :synopsis: Deal with shapefiles and memory layers
     Created on 2018/09/14

.. moduleauthor: Claire POTTIER (CNES DSO/SI/TR)

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import absolute_import, division, print_function, unicode_literals 

import logging
import os
from osgeo import ogr, osr

import cnes.common.service_error as service_error

import cnes.common.lib_lake.lake_db as lake_db


def write_mem_layer_as_shp(in_mem_layer, in_shp_filename):
    """
    Write memory layer in_mem_layer into a shapefile
    
    :param in_mem_layer: memory layer
    :type in_mem_layer: ogr.Layer
    :param in_shp_filename: shapefile full path
    :type in_shp_filename: string
    """
    logger = logging.getLogger("my_shp_file")
    logger.debug("== write_mem_layer_as_shp ==")
    
    shp_driver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Driver for shapefiles
    
    # 1 - Delete output file if already exists
    if os.path.exists(in_shp_filename):
        logger.warning("Output shapefile %s already exists => delete file" % in_shp_filename)
        shp_driver.DeleteDataSource(in_shp_filename)
        
    # 2 - Create output file
    data_source = shp_driver.CreateDataSource(in_shp_filename)
    
    # 3 - Copy memory layer to output 
    layer_name = os.path.basename(in_shp_filename).split(".")[0]
    data_source.CopyLayer(in_mem_layer, str(layer_name))
    
    # 4 - Close output file
    data_source.Destroy()


#######################################


def merge_shp(in_list_shp, flag_add_nogeom_feat=False):
    """
    This function merges shapefiles listed in in_list_shp (typically sub-basin LakeAvg shp into basin-level LakeAvg).
    All polygons in all shapefiles are copied in a memory layer. All fields are copied.

    :param in_list_shp: list of shapefiles full path
    :type in_list_shp: list of string
    :param flag_add_nogeom_feat: =True to add features with no geometry in the output layer, =False otherwise (default)
    :type flag_add_nogeom_feat: boolean
    
    :return out_data_source: data source of output layer
    :rtype out_data_source: OGRdata_source
    :return out_layer: output memory layer
    :rtype out_layer: OGRlayer
    """
    logger = logging.getLogger("my_shp_file")
    logger.debug("== merge_shp ==")
    
    nb_shp = len(in_list_shp)
    logger.debug("%d shapefiles to merge into a single one" % nb_shp)

    # 1 - Create memory data source
    mem_driver = ogr.GetDriverByName(str('MEMORY'))  # Memory driver
    out_data_source = mem_driver.CreateDataSource('memData')
    
    # 2 - Copy first shapefile to memory layer
    # 2.1 - Open input datasource
    shp_driver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Driver for shapefiles
    data_source_1 = shp_driver.Open(in_list_shp[0], 0)  # Open in reading mode
    layer_1 = data_source_1.GetLayer()  # Get the layer
    # 2.2 - Copy layer to memory
    logger.debug("> Copy layer %s (%d features) to memory" % (os.path.basename(in_list_shp[0]), layer_1.GetFeatureCount()))
    out_layer = out_data_source.CopyLayer(layer_1, str('tmp'))

    # 3 - For each other input shapefile
    if nb_shp > 1:
        for cur_shp in in_list_shp[1:]:
    
            # 3.1 - Open the shapefile and get layer
            shp_driver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Shapefile driver
            cur_data_source = shp_driver.Open(cur_shp, 0)  # Open in reading mode
            cur_layer = cur_data_source.GetLayer()  # Get the layer
            logger.debug("> Adding %s (%d features) to memory layer" % (os.path.basename(cur_shp), cur_layer.GetFeatureCount()))
    
            # 3.2 - Merge layer into output layer
            out_data_source, out_layer = merge_2_layers(out_layer, cur_layer, flag_add_nogeom_feat=flag_add_nogeom_feat)
    
            # 3.3 - Close layer
            cur_data_source.Destroy()

    return out_data_source, out_layer


def merge_mem_layer_with_shp(in_list_shp, in_layer, in_continent_id=None, flag_add_nogeom_feat=False):
    """
    This function merges shapefiles listed in in_list_shp with the layer in_layer (typically LakeTile shp with LakeSP memory layer).
    All polygons in all shapefiles are copied in the output memory layer. All fields are copied.

    :param in_list_shp: list of shapefiles full path
    :type in_list_shp: list of string
    :param in_layer: layer in which to merge all objects
    :type in_layer: OGRlayer
    :param in_continent_id: 2-letter identifier of the processed continent
    :type in_continent_id: string
    :param flag_add_nogeom_feat: =True to add features with no geometry in the output layer, =False otherwise (default)
    :type flag_add_nogeom_feat: boolean
    
    :return out_data_source: data source of output layer
    :rtype out_data_source: OGRdata_source
    :return out_layer: output memory layer
    :rtype out_layer: OGRlayer
    """
    logger = logging.getLogger("my_shp_file")
    logger.debug("== merge_mem_layer_with_shp ==")

    # 1 - Create memory data source
    mem_driver = ogr.GetDriverByName(str('MEMORY'))  # Memory driver
    out_data_source = mem_driver.CreateDataSource('memData')

    # 2 - Copy input layer to memory layer
    out_layer = out_data_source.CopyLayer(in_layer, str('tmp'))

    # 3 - For each input shapefile
    for cur_shp in in_list_shp:
        logger.debug("> Adding %s" % os.path.basename(cur_shp))

        # 3.1 - Open the shapefile and get layer
        shp_driver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Shapefile driver
        cur_data_source = shp_driver.Open(cur_shp, 0)  # Open in reading mode
        cur_layer = cur_data_source.GetLayer()  # Get the layer

        # 3.2 - Merge layer into output layer
        out_data_source, out_layer = merge_2_layers(out_layer, cur_layer, in_continent_id=in_continent_id, 
                                                    flag_add_nogeom_feat=flag_add_nogeom_feat)

        # 3.3 - Close layer
        cur_data_source.Destroy()

    return out_data_source, out_layer


def merge_2_layers(in_layer1, in_layer2, in_continent_id=None, flag_add_nogeom_feat=False):
    """
    Merge 2 memory layers
    
    :param in_layer1: first layer
    :type in_layer1: OGRlayer
    :param in_layer2: second layer
    :type in_layer2: OGRlayer
    :param in_continent_id: 2-letter identifier of the processed continent
    :type in_continent_id: string
    :param flag_add_nogeom_feat: =True to add features with no geometry in the output layer, =False otherwise (default)
    :type flag_add_nogeom_feat: boolean
    
    :return out_data_source: data source of output layer
    :rtype out_data_source: OGRdata_source
    :return out_layer: output layer
    :rtype out_layer: OGRlayer
    """
    logger = logging.getLogger("my_shp_file")
    logger.debug("== merge_2_layers ==")
    
    # 1 - Get layer definitions
    layer_defn1 = in_layer1.GetLayerDefn()
    layer_defn2 = in_layer2.GetLayerDefn()

    # 2 - Check if layer definition of in_layer1 and in_layer2 are equal
    fields_name_1 = [layer_defn1.GetFieldDefn(i).GetName() for i in range(layer_defn1.GetFieldCount())]
    fields_name_2 = [layer_defn2.GetFieldDefn(i).GetName() for i in range(layer_defn2.GetFieldCount())]

    if fields_name_1 != fields_name_2:
        message = "Fields of layer %s and layer %s are not identical" % (layer_defn1.GetName(), layer_defn2.GetName())
        logger.error(fields_name_1)
        logger.error(fields_name_2)
        raise service_error.ProcessingError(message, logger)

    # 3 - Retrieve only features corresponding to the specified continent code
    if in_continent_id:
        logger.debug("Filter layer by continent %s " % str(in_continent_id))
        # Find attribute filter depending on the list of fields
        continent_pfaf_id = lake_db.compute_continent_code(in_continent_id)
        if "lake_id" in fields_name_1:
            cond = "lake_id like '" + continent_pfaf_id + "%' or obs_id like '" + continent_pfaf_id + "%'"  # Case of _Obs and _Prior layers
            # Run filter
            in_layer1.SetAttributeFilter(cond)
            in_layer2.SetAttributeFilter(cond)
        elif "obs_id" in fields_name_1:
            cond = "obs_id like '" + continent_pfaf_id + "%'"  # Case of _Unassigned layer
            # Run filter
            in_layer1.SetAttributeFilter(cond)
            in_layer2.SetAttributeFilter(cond)

    nb_feature1 = in_layer1.GetFeatureCount()
    nb_feature2 = in_layer2.GetFeatureCount()
    logger.debug("Merge %d features of layer 1 to %d features of layer 2" %(nb_feature1, nb_feature2))

    # 4 - Create output memory driver and data source
    mem_driver = ogr.GetDriverByName(str('MEMORY'))  # Memory driver
    out_data_source = mem_driver.CreateDataSource('memData')

    # 5 - Projection
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)  # WGS84

    # 6 - Create output layer
    out_layer = out_data_source.CreateLayer(str('tmp'), srs, geom_type=ogr.wkbMultiPolygon)

    # 7 - Init output layer fields
    for ind in range(layer_defn1.GetFieldCount()):
        id_field = layer_defn1.GetFieldDefn(ind)
        out_layer.CreateField(id_field)

    # 8 - Add in_layer1 and in_layer2 to the output layer
    # NB: the Update function doesn't keep feature without geometry (may
    # occur with _Prior layer) => these must be added afterward
    in_layer1.Update(in_layer2, out_layer)
    
    # 9 - Add not observed features (due to use of Update function, cf. above)
    if flag_add_nogeom_feat:
        for cur_layer in [in_layer1, in_layer2]:
            # 9.1 - Retrieve features without a geometry
            if "obs_id" in fields_name_1:
                cur_layer.SetAttributeFilter("obs_id = 'no_data'")
            elif "pass_kept" in fields_name_1:
                cur_layer.SetAttributeFilter("pass_kept = 'no_data'")
            # 9.2 - Add them to output layer
            for cur_feature in cur_layer:
                out_layer.CreateFeature(cur_feature)

    in_layer1.SetAttributeFilter(None)
    in_layer2.SetAttributeFilter(None)

    return out_data_source, out_layer
