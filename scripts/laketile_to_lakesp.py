#!/usr/bin/env python
'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''

import argparse
import configparser as cfg
import glob
import shutil
from netCDF4 import Dataset
import os, sys
from os.path import join, abspath
import subprocess

import tools.my_rdf as my_rdf


def make_input_symlinks(links_dir, laketile_shp_file, laketile_edge_file, laketile_pixcvec_file):
    """ Makes symlinks to pixc with the right name for locnes input """

    def swot_symlink(target, link_name):

        if not os.path.isdir(links_dir):
            os.makedirs(links_dir)
        outpath = join(links_dir, link_name)
        if os.path.islink(outpath): # Overwrite if existing
            print("Overwritting existing {}".format(outpath))
            os.remove(outpath)
        try:
            os.symlink(target, outpath)
        except:
            shutil.copyfile(target, outpath)
            print("symlink impossible => copying file [%s] to [%s]" % (target, outpath))
            print()



    for ext in [".dbf", ".prj", ".shp", ".shp.xml", ".shx"]:
        laketile_shp_file_tmp = laketile_shp_file.replace(".shp", ext)
        laketile_shp_name = os.path.basename(laketile_shp_file_tmp)
        swot_symlink(laketile_shp_file_tmp, laketile_shp_name)


    laketile_edge_name = os.path.basename(laketile_edge_file)
    swot_symlink(laketile_edge_file, laketile_edge_name)

    laketile_pixcvec_name = os.path.basename(laketile_pixcvec_file)
    swot_symlink(laketile_pixcvec_file, laketile_pixcvec_name)

def set_swot_cnes_env():
    print("== Run Locnes from swotCNES ==")
    try:
        lakesrc = os.environ['SWOT_CNES']
        pge_lake_sp_src = os.path.join(lakesrc, "PGE" + os.path.sep + "lake_sp")
        sys.path.insert(0, pge_lake_sp_src)
        if os.environ['SWOT_HYDROLOGY_TOOLBOX']:
            for path in sys.path:
                if path.startswith(os.environ['SWOT_HYDROLOGY_TOOLBOX']):
                    sys.path.remove(path)
        import multi_lake_sp as multi_lake_sp

    except:
        print("swotCNES environment not found : please set SWOT_CNES variable")
        multi_lake_sp = None
        exit()
    return multi_lake_sp

def set_swot_hydro_env():
    print("== Run Locnes from swot-hydrology-toolbox ==")
    try:
        try:
            tbx_path = os.environ['SWOT_HYDROLOGY_TOOLBOX']
        except:
            tbx_path = os.getcwd().replace(os.sep + "scripts", "")
        pge_lake_sp_src = os.path.join(tbx_path, "processing" + os.path.sep + "PGE" + os.path.sep + "lake_sp")
        sys.path.insert(0, pge_lake_sp_src)
        if os.environ['SWOT_CNES']:
            swot_cnes_path = str(os.environ['SWOT_CNES'])
            tmp_sys_path = sys.path.copy()
            for path in tmp_sys_path:
                if str(path).startswith(swot_cnes_path):
                    sys.path.remove(path)

        import multi_lake_sp as multi_lake_sp
    except:
        print("swot-hydrology-toolbox environment not found : please set SWOT_HYDROLOGY_TOOLBOX variable")
        multi_lake_sp = None
        exit()
    return multi_lake_sp

def call_multi_lake_sp(parameter_lakesp, lakesp_dir, laketile_dir, multi_lake_sp):

    config = cfg.ConfigParser()
    config.read(parameter_lakesp)


    # Fill missing values in pge_lake_tile rdf file
    config.set('PATHS', "LakeTile shp directory", laketile_dir)
    config.set('PATHS', "LakeTile edge directory", laketile_dir)
    config.set('PATHS', "LakeTile pixcvec directory", laketile_dir)
    config.set('PATHS', "Output directory", lakesp_dir)
    config.set('LOGGING', "logFile", os.path.join(lakesp_dir, "LogFile.log"))


    # Write the parameter file
    lakesp_cfg = join(lakesp_dir, "lakesp.cfg")

    with open(lakesp_cfg, 'w') as cfg_file:
        config.write(cfg_file)


    print("== Run LakeSP ==")

    my_params = multi_lake_sp.read_command_file(lakesp_cfg)  # Read variables in command file
    # try:
    # 2 - Initialization
    multi_lake_sp_proc = multi_lake_sp.MultiLakeSP(my_params)

    # 3 - Run pre-processing
    multi_lake_sp_proc.run_preprocessing()

    # 4 - Run processing
    multi_lake_sp_proc.run_processing()


    print()


#######################################


def main():
    """When run as a script"""
    parser = argparse.ArgumentParser()
    parser.add_argument('output_dir', help="output directory", type=str)
    parser.add_argument('parameter_multi_lakesp', help="parameter file multi lake sp", type=str)
    parser.add_argument('laketile_annotation_file_list', help="laketile annotation file list (output from rivertile_to_laketile.py) ", type=str, nargs='+')
    parser.add_argument('-swotCNES', action='store_true', dest='swotCNES', default=None, help='only for CNES developer users, switch to swotCNES processing env')

    args = vars(parser.parse_args())

    if args['swotCNES']:
        multi_lake_sp = set_swot_hydro_env()
    else:
        multi_lake_sp = set_swot_hydro_env()

    print("===== laketile_to_lakesp = BEGIN =====")
    print("")

    # Print on console
    if len(args['laketile_annotation_file_list']) == 0:
        print("> NO laketile annotation file to deal with")
    else:
        print("> %d laketile annotation file(s) to deal with" % len(args['laketile_annotation_file_list']))
    print()

    # Create symlinks to input with the right name convention
    links_dir = join(args['output_dir'], "inputs")
    # Read config from config repository by default, or from the script arguments
    parameter_lakesp = args["parameter_multi_lakesp"]

    print()
    for laketile_annotation in args['laketile_annotation_file_list']:

        print(">>>>> Dealing with lake annotation file %s <<<<<" % laketile_annotation)
        print()

        # Load annotation file
        ann_rdf = my_rdf.myRdfReader(os.path.abspath(laketile_annotation))
        ann_cfg = ann_rdf.parameters
        print(ann_cfg)
        # Prepare args for pge_lake_tile.py
        laketile_shp_file = abspath(ann_cfg['laketile_shp file'])
        print("Laketile_shp file = %s" % laketile_shp_file)
        laketile_edge_file = abspath(ann_cfg['laketile_edge file'])
        print("Laketile_shp file = %s" % laketile_edge_file)
        laketile_pixcvec_file = abspath(ann_cfg['laketile_pixcvec file'])
        print("Laketile_pixcvec file = %s" % laketile_pixcvec_file)

        print("logFile = %s" % args['output_dir'])
        make_input_symlinks(links_dir, laketile_shp_file, laketile_edge_file, laketile_pixcvec_file)


    print()
    print()

    call_multi_lake_sp(parameter_lakesp, args['output_dir'], links_dir, multi_lake_sp)

    print("===== laketile_to_lakesp = END =====")


#######################################


if __name__ == "__main__":
    main()
