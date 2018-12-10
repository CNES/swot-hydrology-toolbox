#!/usr/bin/env python

import os
import glob
from os.path import join, abspath
import sys
import subprocess
from netCDF4 import Dataset
import argparse
import numpy as np
import cnes.modules.geoloc.lib.pixc_to_shp
import my_rdf
import csv
import configparser as cfg

def make_tail_proc_path(annotation_file, output_dir, suffix):
    "Clone the tail of proc directory structure"

    dir_parts = os.path.dirname(abspath(annotation_file)).split(os.path.sep)
    tail_dir = join(output_dir, dir_parts[-2], dir_parts[-1])
    river_dir = join(tail_dir, suffix)

    if not os.path.isdir(river_dir):
        os.makedirs(river_dir)

    return abspath(tail_dir), abspath(river_dir)
    

def make_input_symlinks(links_dir, pixc_file, pixcvec_file, cycle_number, pass_number, tile_number, mission_start_time):
    """ Makes symlinks to pixc with the right name for locnes input """

    def swot_symlink(target, name):
        stop_date_time = mission_start_time
        crid = "Dx0000"
        product_counter = "01"
        outname = name.format(cycle_number, pass_number, tile_number, mission_start_time, stop_date_time, crid, product_counter)
        outpath = join(links_dir, outname)
        if os.path.islink(outpath): # Overwrite if existing
            print("Overwritting existing {}".format(outpath))
            os.remove(outpath)
        os.symlink(target, outpath)
        return outpath

    pixcname = swot_symlink(pixc_file, "SWOT_L2_HR_PIXC_{}_{}_{}_{}_{}_{}_{}.nc")
    pixcvecname = swot_symlink(pixcvec_file, "SWOT_L2_HR_PIXCVecRiver_{}_{}_{}_{}_{}_{}_{}.nc")
    return pixcname, pixcvecname
    

def call_pge_lake_tile(parameter_laketile, lake_dir, pixc_file, pixcvec_file, cycle_number, pass_number, tile_number, mission_start_time, force_disable_improved_geolocation=False):


    config = cfg.ConfigParser()
    config.read(parameter_laketile)
    
    # Create symlinks to input with the right name convention
    links_dir = join(lake_dir, "inputs")
    if not os.path.isdir(links_dir):
        os.makedirs(links_dir)

    pixcname, pixcvecname = make_input_symlinks(links_dir, pixc_file, pixcvec_file, cycle_number, pass_number, tile_number, mission_start_time)

    # Fill missing values in pge_lake_tile rdf file
    config.set('PATHS', "PIXC file", pixcname)
    config.set('PATHS', "PIXCVecRiver file", pixcvecname)
    config.set('PATHS', "Output directory", lake_dir)

    if force_disable_improved_geolocation:
        config.set('OPTIONS', 'Improve geolocation', 0)
        
    # Write the parameter file
    laketile_cfg = join(lake_dir, "laketile.cfg")

    with open(laketile_cfg, 'w') as cfg_file:
        config.write(cfg_file)

    # Call pge_lake_tile
    pge_lake_tile = join(os.environ['SWOT_HYDROLOGY_TOOLBOX'], 'processing', 'src', 'cnes', 'sas', 'lake_tile', 'pge_lake_tile.py')
    subprocess.check_call([pge_lake_tile, laketile_cfg, '-shp'])

def main():
    """When run as a script"""
    parser = argparse.ArgumentParser()
    parser.add_argument('river_annotation_file', type=str)
    parser.add_argument('output_dir', type=str)

    parser.add_argument('--parameter_laketile', default= None, type=str)
    parser.add_argument("--nogdem",
        help="If true, don't call riverobs with gdem",
        nargs='?', type=bool, default=False, const=True)
    parser.add_argument(
        '-f', '--force', action='store_true', dest='force', default=False,
        help='Force overwrite existing outputs; default is to quit')
    args = vars(parser.parse_args())

    # unit or multi river annotation file test
    if os.path.isfile(args['river_annotation_file']):
        river_files = glob.glob(args['river_annotation_file'])
        print(river_files)
    else:
        river_files = glob.glob(os.path.join(args['river_annotation_file'], "*.rdf"))

    print(river_files, args['river_annotation_file'])
    for river_annotation in river_files:
    
	    # Load annotation file
	    print(river_annotation)
	    ann_rdf = my_rdf.myRdfReader(os.path.abspath(river_annotation))
	    ann_cfg = ann_rdf.parameters
	
	    # Clone the tail of proc directory structure
	    tail_dir, lake_dir_pixc = make_tail_proc_path(river_annotation, args['output_dir'], 'pixc')
	
	    if os.path.exists(lake_dir_pixc) and not args["force"]:
	        print("Output lake directory exists. Stopping (use -f to force).")
	        return
	
	    # Prepare args for pge_lake_tile.py
	    pixc_file    = abspath(ann_cfg['pixc file'])
	    pixcvec_file = abspath(ann_cfg['pixcvec file'])
	
	    # Read pixc file attributes and format with leading zeros
	    with Dataset(ann_cfg["pixc file"], "r") as pixc_dataset:
	                
	        cycle_number = "{:03d}".format(pixc_dataset.groups['pixel_cloud'].getncattr("cycle_number"))
	        pass_number = "{:03d}".format(pixc_dataset.groups['pixel_cloud'].getncattr("pass_number"))
	        tile_number = pixc_dataset.groups['pixel_cloud'].getncattr("tile_ref")
	        try:
	            mission_start_time = pixc_dataset.groups['pixel_cloud'].getncattr("mission_start_time")
	        except : 
	            mission_start_time = pixc_dataset.groups['pixel_cloud'].getncattr("start_time")
	            
	    # Read config from config repository by default, or from the script arguments
	    parameter_laketile = args["parameter_laketile"]
	
	    print(cycle_number, pass_number, tile_number, mission_start_time)
	    call_pge_lake_tile(parameter_laketile, lake_dir_pixc,
	                       ann_cfg["pixc file"], ann_cfg["pixcvec file"],
	                       cycle_number, pass_number, tile_number, mission_start_time)

if __name__ == "__main__":
    main()

