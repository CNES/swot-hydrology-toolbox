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
# FIN-HISTORIQUE
# ======================================================

"""
extract occurence ROI from pekel with margins

depends on Gdal et netCDF Gdal driver

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""

import argparse
from subprocess import run, PIPE
import shlex
import json
import numpy as np
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract occurrence from gdem dem")
    parser.add_argument("--gdem", type=str, help="input DEM file .nc file", required=True)
    parser.add_argument("--input", type=str, help="input pekel vrt", required=True)
    parser.add_argument("--output", type=str, help="output pekel file (default: occurrence.tif)", default="occurrence.tif")
    parser.add_argument("--margin_lon", type=float, help="margin in longitude (default: 2.0)", default=2.0)
    parser.add_argument("--margin_lat", type=float, help="margin in latitude (default: 2.0)", default=2.0)

    args = parser.parse_args()

    input_gdem_dem = args.gdem
    input_pekel = args.input
    output_pekel = args.output

    # Retrieve informations from dem
    cmd = "gdalinfo -json NETCDF:\"{}\":elevation".format(input_gdem_dem)
    info = run(shlex.split(cmd),stdout=PIPE)
    if info.returncode != 0:
        sys.stderr.write("Error during retrieving gdem informations\n")
        sys.exit(1)
    info_str = info.stdout.decode("utf-8")
    info_json = json.loads(info_str)
    upper_left = np.array(info_json['cornerCoordinates']['upperLeft'])
    lower_right = np.array(info_json['cornerCoordinates']['lowerRight'])
 
    # Define margins
    margin_lon = args.margin_lon
    margin_lat = args.margin_lat
    extended_upper_left = upper_left + np.array([-margin_lon, +margin_lat])
    extended_lower_right = lower_right + np.array([+margin_lon, -margin_lat])

    # Extract Pekel ROI
    info = run(["gdal_translate",
                         input_pekel,
                         output_pekel,
                         "-projwin",
                         str(extended_upper_left[0]),
                         str(extended_upper_left[1]),
                         str(extended_lower_right[0]),
                         str(extended_lower_right[1]),
                        "-ot","Float32"])
    if info.returncode != 0:
        sys.stderr.write("Error during Pekel extraction\n")
        sys.exit(1)

