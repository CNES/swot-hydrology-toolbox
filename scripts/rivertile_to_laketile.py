#!/usr/bin/env python
'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''

import argparse
import configparser as cfg
import glob
from netCDF4 import Dataset
import os
from os.path import join, abspath
import subprocess

import tools.my_rdf as my_rdf


def make_input_symlinks(links_dir, pixc_file, pixcvec_file, cycle_number, pass_number, tile_number, start_time, stop_time):
    """ Makes symlinks to pixc with the right name for locnes input """

    def swot_symlink(target, name):
        crid = "Dx0000"
        product_counter = "01"
        outname = name.format(cycle_number, pass_number, tile_number, start_time, stop_time, crid, product_counter)
        if not os.path.isdir(links_dir):
            os.makedirs(links_dir)
        outpath = join(links_dir, outname)
        if os.path.islink(outpath): # Overwrite if existing
            print("Overwritting existing {}".format(outpath))
            os.remove(outpath)
        try:
            os.symlink(target, outpath)
            flag_rename = False
        except:
            outpath = os.path.join(os.path.dirname(target), outname)
            os.rename(target, outpath)
            flag_rename = True
            print("symlink impossible => change filename [%s] to [%s]" % (target, outpath))
            print()
        return outpath, flag_rename

    flag_rename_pixc = False
    if not os.path.basename(pixc_file).startswith("SWOT_L2_HR_PIXC"):
        pixcname, flag_rename_pixc = swot_symlink(pixc_file, "SWOT_L2_HR_PIXC_{}_{}_{}_{}_{}_{}_{}.nc")
    else:
        pixcname = pixc_file
    
    flag_rename_pixcvec = False
    if not os.path.basename(pixcvec_file).startswith("SWOT_L2_HR_PIXCVecRiver"):
        pixcvecname, flag_rename_pixcvec = swot_symlink(pixcvec_file, "SWOT_L2_HR_PIXCVecRiver_{}_{}_{}_{}_{}_{}_{}.nc")
    else:
        pixcvecname = pixcvec_file
    
    return pixcname, flag_rename_pixc, pixcvecname, flag_rename_pixcvec


def execute(cmd):
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='') # process line here

    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, p.args)
    

def call_pge_lake_tile(parameter_laketile, lake_dir, pixc_file, pixcvec_file, cycle_number, pass_number, tile_number, start_time, stop_time, env=None):

    config = cfg.ConfigParser()
    config.read(parameter_laketile)
    
    # Create symlinks to input with the right name convention
    links_dir = join(lake_dir, "inputs")

    pixcname, flag_rename_pixc, pixcvecname, flag_rename_pixcvec = make_input_symlinks(links_dir, pixc_file, pixcvec_file, cycle_number, pass_number, tile_number, start_time, stop_time)

    # Fill missing values in pge_lake_tile rdf file
    config.set('PATHS', "PIXC file", pixcname)
    config.set('PATHS', "PIXCVecRiver file", pixcvecname)
    config.set('PATHS', "Output directory", lake_dir)
    logFile_path = config.get('LOGGING', "logFile").replace("REPLACE_ME/", lake_dir+os.path.sep)
    config.set('LOGGING', "logFile", logFile_path)
        
    # Write the parameter file
    laketile_cfg = join(lake_dir, "laketile.cfg")

    with open(laketile_cfg, 'w') as cfg_file:
        config.write(cfg_file)

    # Call pge_lake_tile
    # Path to top of toolbox
    try:
        tbx_path = os.environ['SWOT_HYDROLOGY_TOOLBOX']
    except:
        tbx_path = os.getcwd().replace(os.sep+"scripts", "")
        
    # Build LOCNES/lake_tile main lib
    if (env == 'swotCNES'):
        print("Switch to swotCNES env instead of swot-hydrology-toolbox processing")
        pge_lake_tile = '/work/ALT/swot/swotdev/desrochesd/swot-sds/swotCNES/PGE/lake_tile/pge_lake_tile.py'
    else:
        pge_lake_tile = join(tbx_path, 'processing', 'PGE', 'lake_tile', 'pge_lake_tile.py')

        
    # Build command
    cmd = "{} {}".format(pge_lake_tile, laketile_cfg)
    print ("> executing:", cmd) 
    print()
    #subprocess.check_call(cmd, shell=True)
    execute(cmd)
    print()
    print("== Execution OK")
    print()
        
    # Rename to old input filenames if files had been renamed
    if flag_rename_pixc:
        print("Get back to old PixC filename [%s] to [%s]" % (pixcname, pixc_file))
        print()
        os.rename(pixcname, pixc_file)
    if flag_rename_pixcvec:
        print("Get back to old PIXCVecRiver filename [%s] to [%s]" % (pixcvecname, pixcvec_file))
        print()
        os.rename(pixcvecname, pixcvec_file)


#######################################


def main():
    """When run as a script"""
    parser = argparse.ArgumentParser()
    parser.add_argument('river_annotation_file', help="river annotation file (output from l2pixc_to_rivertile) or directory", type=str)
    parser.add_argument('output_dir', help="output directory", type=str)
    parser.add_argument('parameter_laketile', help="parameter file", type=str)
    parser.add_argument("--nogdem", help="if true, don't call riverobs with gdem", nargs='?', type=bool, default=False, const=True)
    parser.add_argument(
        '-swotCNES', action='store_true', dest='swotCNES', default=None,
        help='only for CNES developer users, switch to swotCNES processing env')
            
    args = vars(parser.parse_args())
                  
    if args['swotCNES']:    
        env = 'swotCNES'
    else:
        env = None
        
    print("===== rivertile_to_laketile = BEGIN =====")
    print("")

    # unit or multi river annotation file test
    if os.path.isfile(args['river_annotation_file']):
        river_files = glob.glob(args['river_annotation_file'])
    else:
        river_files = glob.glob(os.path.join(args['river_annotation_file'], "*.rdf"))
    # Print on console
    if len(river_files) == 0:
        print("> NO river annotation file to deal with")
    else:
        print("> %d river annotation file(s) to deal with" % len(river_files))
    print()
    
    for river_annotation in river_files:
        
        print(">>>>> Dealing with river annotation file %s <<<<<" % river_annotation)
        print()
    
        # Load annotation file
        ann_rdf = my_rdf.myRdfReader(os.path.abspath(river_annotation))
        ann_cfg = ann_rdf.parameters
	
        # Prepare args for pge_lake_tile.py
        pixc_file    = abspath(ann_cfg['pixc file'])
        print(". PixC file = %s" % pixc_file)
        pixcvec_file = abspath(ann_cfg['pixcvec file'])
        print(". PIXCVecRiver file = %s" % pixcvec_file)
	
        # Read pixc file attributes and format with leading zeros
        with Dataset(ann_cfg["pixc file"], "r") as pixc_dataset:
	                
            cycle_number = "{:03d}".format(pixc_dataset.getncattr("cycle_number"))
            pass_number = "{:03d}".format(pixc_dataset.getncattr("pass_number"))
            tile_number = pixc_dataset.getncattr("tile_name").split("_")[1]
            try:
                start_time = pixc_dataset.getncattr("start_time")
            except: 
                start_time = "YYYYMMDDThhmmss"
            try:
                stop_time = pixc_dataset.getncattr("stop_time")
            except: 
                stop_time = "YYYYMMDDThhmmss"
            
            print(". Cycle number = %s" % cycle_number)
            print(". Orbit number = %s" % pass_number)
            print(". Tile number = %s" % tile_number)
            print(". Start time = %s" % start_time)
            print(". Stop time = %s" % stop_time)
            
        print()
        	            
        # Read config from config repository by default, or from the script arguments
        parameter_laketile = args["parameter_laketile"]
	
        call_pge_lake_tile(parameter_laketile, args['output_dir'],
                           ann_cfg["pixc file"], ann_cfg["pixcvec file"],
                           cycle_number, pass_number, tile_number, start_time, stop_time, env=env)
        
        print()
        print()

    print("===== rivertile_to_laketile = END =====")


#######################################


if __name__ == "__main__":
    main()
