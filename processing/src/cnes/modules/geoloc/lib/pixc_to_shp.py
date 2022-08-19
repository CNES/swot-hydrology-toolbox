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
# VERSION:3.0.0:DM:#91:2021/03/12:Poursuite industrialisation
# VERSION:3.1.0:DM:#91:2021/05/21:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
'''
.. module:: pixc_to_shp.py

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


import sys, os
from collections import OrderedDict
import argparse
import math
import numpy as np
from osgeo import ogr
import numpy.ma as ma

from collections import OrderedDict
import netCDF4 as nc
import cnes.common.lib.my_tools as my_tools

def get_valid_variable(var_name, var_data):
    try:
        var_data = ma.filled(var_data, fill_value=-9999.0).astype(np.float)
    except:
        if var_name == "node_id":
            tmp_id = nc.chartostring(var_data)
            var_data = np.char.asarray(tmp_id, itemsize=14)
            var_data_array = np.ones(var_data.shape, dtype=np.float) * -9999.0
            var_data_array[var_data != ''] = np.asarray(var_data[var_data != ''], dtype=np.float)
            var_data = var_data_array
        elif var_name == "reach_id":
            tmp_id = nc.chartostring(var_data)
            var_data = np.char.asarray(tmp_id, itemsize=11)
            var_data_array = np.ones(var_data.shape, dtype=np.float) * -9999.0
            var_data_array[var_data != ''] = np.asarray(var_data[var_data != ''], dtype=np.float)
            var_data = var_data_array
        elif var_name == "lake_id":
            tmp_id = nc.chartostring(var_data)
            var_data = np.char.asarray(tmp_id, itemsize=10)
            var_data_array = np.ones(var_data.shape, dtype=np.float) * -9999.0
            var_data_array[var_data != ''] = np.asarray(var_data[var_data != ''], dtype=np.float)
            var_data = var_data_array
        elif var_name == "obs_id":
            tmp_id = nc.chartostring(var_data)
            var_data = np.char.asarray(tmp_id, itemsize=13)
            var_data_array = np.ones(var_data.shape, dtype=np.float) * -9999.0
            var_data_array[var_data != ''] = np.asarray(var_data[var_data != ''], dtype=np.float)
            var_data = var_data_array

    return var_data

def create_shp(output_name, var_names):
    shp_driver = ogr.GetDriverByName("ESRI Shapefile")
    shp_ds = shp_driver.CreateDataSource(output_name)

    # Set spatial projection
    srs = ogr.osr.SpatialReference()
    srs.ImportFromEPSG(4326)

    # Set the lake_layer name
    layer_name = str((os.path.splitext(os.path.basename(output_name))[0]).lower())

    # Creating output lake_layer
    shp_layer = shp_ds.CreateLayer(str(layer_name), srs, geom_type=ogr.wkbPoint)

    # Create fields
    for i, var_name in enumerate(var_names):
        print('Creating output fields: %s %s' % (var_name, ogr.OFTReal))
        if len(var_name) > 10:
            print('         output field name reduced to: %s' % (var_name[:10]))
            var_name = var_name[:10]
        field_defn = ogr.FieldDefn(var_name, ogr.OFTReal)
        field_defn.SetWidth(20)
        field_defn.SetPrecision(4)
        shp_layer.CreateField(field_defn)

    return shp_ds, shp_layer

def pixc_to_shp(input_name, output_name, lat_name, lon_name, var_names, group_name=None, progress=False, wateronly=False):

    pixc = nc.Dataset(input_name, "r")
    variables = {}
    if group_name is None:
        latitude = get_valid_variable(lat_name, pixc.variables[lat_name][:])
        longitude = get_valid_variable(lon_name, my_tools.convert_to_m180_180(pixc.variables[lon_name][:]))
        for var_name in var_names :
            var_data = pixc.variables[var_name][:]
            var_data = get_valid_variable(var_name, var_data)
            variables[var_name] = var_data
    else:
        latitude = get_valid_variable(lat_name, pixc.groups[group_name].variables[lat_name][:])
        longitude = get_valid_variable(lon_name, my_tools.convert_to_m180_180(pixc.groups[group_name].variables[lon_name][:]))
        for var_name in var_names :
            variables[var_name] = get_valid_variable(var_name, pixc.groups[group_name].variables[var_name][:])

    if wateronly:
        idx = np.where(np.logical_and(np.logical_and(variables["classification"] != 1, variables["classification"] != 2), variables["classification"] != 22))
        latitude = latitude[idx]
        longitude = longitude[idx]
        for var_name in var_names:
            variables[var_name] = variables[var_name][idx]

    nb_points = latitude.size

    shp_ds, shp_layer = create_shp(output_name, var_names)

    print("Writting oupt shapefile %s with attributes %s " %(output_name, " ".join(var_names)))

    for i in range(nb_points):
        feature = ogr.Feature(shp_layer.GetLayerDefn())

        point = ogr.CreateGeometryFromWkt("POINT(%f %f)" % (longitude[i], latitude[i]))
        feature.SetGeometry(point)

        for var_name in var_names:
            var_value = variables[var_name][i]
            feature.SetField(var_name[:10], var_value)

        shp_layer.CreateFeature(feature)
        if progress and i % 10000 == 0:
            print("Writing shp points: {:.1f}% done ({} / {})".format(100*(i+1)/nb_points, i+1, nb_points))
    print()
    shp_ds.Destroy()
    pixc.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input netcdf file', type=str)
    parser.add_argument('output', help='output shp file', type=str)
    parser.add_argument("lat", help="Name of the latitude variable", type=str)
    parser.add_argument("lon", help="Name of the longitude variable", type=str)
    parser.add_argument('-v','--variables', nargs='+', help='List of variables', required=True)
    parser.add_argument('-g', '--group', help="Optional netcdf group", required=False, default=None)
    args = parser.parse_args()

    pixc_to_shp(args.input, args.output, args.lat, args.lon, args.variables, group_name=args.group, progress=True)
