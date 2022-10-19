# -*- coding: utf8 -*-
"""

.. module author: Claire POTTIER - CNES DSO/SI/TR
.. modified by : Cécile CAZALS - CS
.. modified by : Manon DELHOUME - CS Group

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from __future__ import absolute_import, division, print_function, unicode_literals

import Utils.my_netcdf_file as my_nc
import Utils.product_netcdf as my_prod_nc
import numpy as np

import argparse

def shift_pixc(in_pixc_path, out_pixc_path, degree_to_shift, tile_name):
    """
    This class is dedicated to produce modified pixc (with new longitudes values) from real pixc.
    :param in_pixc_path: real pixc path
    :param out_pixc_path: output pixc path
    :param degree_to_shift: distance (in degree) to shift the pixc longitudes.
    Negative value if you want to move the pixels to the West.
    :param tile_name: name of the new destination tile (where the pixels will be moved)
    :return: -
    """

    # 1 - Open pixel cloud file in reading mode
    pixc_reader = my_nc.MyNcReader(in_pixc_path)
    pixc_group = pixc_reader.content.groups['pixel_cloud']
    sensor_group = pixc_reader.content.groups['tvp']
    noise_group = pixc_reader.content.groups['noise']

    # 2 - Init a PixcProduct object
    nc_writer = my_prod_nc.PixcProduct()

    # 3 - Extract swath informations from tile_name
    swath_side = tile_name[-1]
    pass_number = tile_name.split("_")[0]
    tile_number = tile_name.split("_")[-1][0:3]

    # 4 - Modify some pixc attributes related to longitudes
    metadata = {}
    for key in (pixc_reader.get_list_att()):
        if key == "geospatial_lon_min" or key == "inner_first_longitude" or key == "geospatial_lon_max" \
                or key == "inner_last_longitude" or key == "outer_last_longitude" or key == "outer_first_longitude" :
            metadata[key] = float(pixc_reader.get_att_value(key)) + degree_to_shift
        elif key == "pass_number":
            metadata[key] = "%s" % pass_number
        elif key == "tile_number":
            metadata[key] = "%s" % tile_number
        elif key == "swath_side":
            metadata[key] = "%s" % swath_side
        elif key == "tile_name":
            metadata[key] = "%s" % tile_name
        else:
            metadata[key] = pixc_reader.get_att_value(key)
    nc_writer.set_metadata_val(metadata)

    # Copy of the invariant attributes from init pixc
    metadata = {}
    for key in (pixc_reader.get_list_att( in_group=pixc_group)):
        metadata[key] = pixc_reader.get_att_value(key,  in_group=pixc_group)
    nc_writer.set_metadata_val(metadata, group="pixel_cloud")

    metadata = {}
    for key in (pixc_reader.get_list_att( in_group=sensor_group)):
        metadata[key] = pixc_reader.get_att_value(key,  in_group=sensor_group)
    nc_writer.set_metadata_val(metadata, group="tvp")

    # 5 - Update group dimension (points, num_tvps, num_lines)
    # And moves the pixels beyond -180° to 180°
    nc_writer.set_dim_val("points", 0)
    pixel_cloud_vars_val = {}
    for key in pixc_reader.get_list_var(in_group=pixc_group):
        pixel_cloud_vars_val[key] = pixc_reader.get_var_value(key,  in_group=pixc_group)
        if key == "longitude":
            pixel_cloud_vars_val[key] += degree_to_shift
            pixel_cloud_vars_val[key][pixel_cloud_vars_val[key]<=-180] += 360

    nc_writer.set_dim_val("num_tvps", pixc_reader.get_dim_value( "num_tvps", in_group=sensor_group))
    tvp_vars_val = {}
    for key in pixc_reader.get_list_var(in_group=sensor_group):
        tvp_vars_val[key] = pixc_reader.get_var_value(key, in_group=sensor_group)
        if key == "longitude":
            tvp_vars_val[key] += degree_to_shift

    nc_writer.set_dim_val("num_lines", pixc_reader.get_dim_value("num_lines", in_group=noise_group))
    noise_vars_val = {}
    for key in pixc_reader.get_list_var(in_group=noise_group):
        noise_vars_val[key] = pixc_reader.get_var_value(key,  in_group=noise_group)
        if key == "longitude":
            noise_vars_val[key] += degree_to_shift

    # 5 - Group dictionaries of group variables in a parent dictionary
    all_vars_val_per_group = {}
    all_vars_val_per_group["pixel_cloud"] = pixel_cloud_vars_val
    all_vars_val_per_group["tvp"] = tvp_vars_val
    all_vars_val_per_group["noise"] = noise_vars_val

    # 6 - Write the pixc product
    nc_writer.write_product(out_pixc_path, in_vars_value=all_vars_val_per_group)



def shift_pixc_vec(in_pixc_vec_river_path, out_pixc_vec_river_path, degree_to_shift, tile_name):
    """
    This class is dedicated to produce modified pixcvec river (with new longitudes values) from real pixcvec river.
    :param in_pixc_vec_river_path: real pixcvec river path
    :param out_pixc_vec_river_path: output pixcvec river path
    :param degree_to_shift: distance (in degree) to shift the pixcvec longitudes.
    Negative value if you want to move the pixels to the West.
    :param tile_name: name of the new destination tile (where the pixels will be moved)
    :return: -
    """

    # 1 - Open pixcvec river file in reading mode
    pixc_reader = my_nc.MyNcReader(in_pixc_vec_river_path)

    # 2 - Init a PixcVecRiverProduct object
    nc_writer = my_prod_nc.PixcVecRiverProduct()

    # 3 -Extract swath informations from tile_name
    swath_side = tile_name[-1]
    pass_number = tile_name.split("_")[0]
    tile_number = tile_name.split("_")[-1][0:3]

    # 4 - Modify some pixcvec attributes related to longitudes
    metadata = {}
    for key in (pixc_reader.get_list_att()):
        if key == "geospatial_lon_min" or key == "inner_first_longitude" or key == "geospatial_lon_max" \
                or key == "inner_last_longitude" or key == "outer_last_longitude" or key == "outer_first_longitude" :
            metadata[key] = float(pixc_reader.get_att_value(key)) + degree_to_shift
        elif key == "pass_number":
            metadata[key] = "%s" % pass_number
        elif key == "tile_number":
            metadata[key] = "%s" % tile_number
        elif key == "swath_side":
            metadata[key] = "%s" % swath_side
        elif key == "tile_name":
            metadata[key] = "%s" % tile_name
        else:
            metadata[key] = pixc_reader.get_att_value(key)
    nc_writer.set_metadata_val(metadata)


    # 5 - Update group dimension
    # And moves the pixels beyond -180° to 180°
    nc_writer.set_dim_val("points", 0)
    pixel_cloud_vars_val = {}
    for key in pixc_reader.get_list_var():
        pixel_cloud_vars_val[key] = pixc_reader.get_var_value(key)
        if key == "longitude":
            pixel_cloud_vars_val[key] += degree_to_shift
            pixel_cloud_vars_val[key][pixel_cloud_vars_val[key]<=-180] += 360


    # 6 - Group dictionaries of group variables in a parent dictionary
    all_vars_val_per_group = {}
    all_vars_val_per_group["points"] = pixel_cloud_vars_val

    # 7 - Write the pixcvec river product
    nc_writer.write_product(out_pixc_vec_river_path, in_vars_value=all_vars_val_per_group)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-pixc_in', dest='pixc_in', help='input netcdf file', type=str, default=False)
    parser.add_argument('-pixc_out', dest='pixc_out', help='output netcdf file', type=str, default=False)
    parser.add_argument('-pixc_vec_river_in', dest='pixc_vec_river_in', help='input netcdf file', type=str, default=False)
    parser.add_argument('-pixc_vec_river_out', dest='pixc_vec_river_out', help='output netcdf file', type=str, default=False)
    parser.add_argument('-degree_to_shift', help='value in degree to shift the longitude', type=str, default=-4.1)
    parser.add_argument('-tile_name', help='name of the new destination tile', type=str, default="573_274L")
    args = parser.parse_args()

    if args.pixc_in and args.pixc_out and args.pixc_vec_river_in and args.pixc_vec_river_out:
        empty_pixc(args.pixc_in, args.pixc_out, args.degree_to_shift, args.tile_name)
        empty_pixc_vec(args.pixc_vec_river_in, args.pixc_vec_river_out, args.degree_to_shift, args.tile_name)
    else :
        print("Unknown option")