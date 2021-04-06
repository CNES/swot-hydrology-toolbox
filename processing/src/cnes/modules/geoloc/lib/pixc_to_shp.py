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
'''
.. module:: pixc_to_shp.py

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


import sys
from collections import OrderedDict
import argparse
import math
import numpy as np

import fiona
import fiona.crs
import shapely.geometry as geometry
import netCDF4 as nc
import cnes.common.lib.my_tools as my_tools

def pixc_to_shp(input_name, output_name, lat_name, lon_name, var_names, group_name=None, progress=False):
    pixc = nc.Dataset(input_name, "r")

    if group_name is None:
        latitude = pixc.variables[lat_name][:]
        longitude = my_tools.convert_to_m180_180(pixc.variables[lon_name][:])
        variables = [pixc.variables[var_name][:] for var_name in var_names]
    else:
        latitude = pixc.groups[group_name].variables[lat_name][:]
        longitude = my_tools.convert_to_m180_180(pixc.groups[group_name].variables[lon_name][:])
        variables = [pixc.groups[group_name].variables[var_name][:] for var_name in var_names]

    nb_points = latitude.size

    driver = "ESRI Shapefile"
    crs = fiona.crs.from_epsg(4326) # WGS84

    schema = {'properties': OrderedDict([(lon_name, 'float:24.15'), (lat_name, 'float:24.15')] + [(var_name, 'float:24.15') for var_name in var_names]), 'geometry': 'Point'}

    sys.stdout.write("Writing shp points")
    with fiona.open(output_name,'w', driver=driver, crs=crs, schema=schema) as c:
        for i in range(nb_points):
            point = geometry.Point(longitude[i], latitude[i])

            prop = {lon_name: float(point.coords.xy[0][0]),
                    lat_name: float(point.coords.xy[1][0])}

            for var_name, var_values in zip(var_names, variables):
                if np.ma.is_masked(var_values[i]) or math.isnan(var_values[i]) :
                    prop[var_name] = float(-9999.)
                else:
                    prop[var_name] = float(var_values[i])

            c.write({'geometry': geometry.mapping(point), 'properties': prop})

            if progress and i % 100 == 0:
                sys.stdout.write("\rWriting shp points: {:.1f}% done ({} / {})".format(100*(i+1)/nb_points, i+1, nb_points))
    sys.stdout.write("\n")

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
