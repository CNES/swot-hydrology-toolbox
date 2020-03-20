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
import numpy as np
import os, sys
from os.path import join, abspath
import subprocess

import tools.my_rdf as my_rdf
import tools.my_filenames as my_names
#Import of PGE
try:
    tbx_path = os.environ['SWOT_HYDROLOGY_TOOLBOX']
except:
    tbx_path = os.getcwd().replace(os.sep + "scripts", "")

try :
    lakesrc = os.environ['SWOT_CNES']
    sys.path.insert(0, lakesrc)
    from PGE.lake_tile import pge_lake_tile as pge_lake_tile
except:
    sys.path.insert(0, tbx_path)
    from processing.PGE.lake_tile import pge_lake_tile as pge_lake_tile


def write_annotation_file(ann_file, laketile_shp_file, laketile_edge_file, laketile_pixcvec_file):
    """
    write the river-annotation.rdf file so that lake processor can run
    """
    f = open(ann_file,'w')
    f.write("laketile_shp file = %s\n"%laketile_shp_file)
    f.write("laketile_edge file = %s\n"%laketile_edge_file)
    f.write("laketile_pixcvec file = %s\n"%laketile_pixcvec_file)
    f.close()

def gdem_pixc_to_pixc( pixc_in, pixc_out):
    print("Copying file form %s to %s " %(pixc_in, pixc_out))
    if os.path.exists(pixc_out):
        os.remove(pixc_out)

    with Dataset(pixc_in) as src, Dataset(pixc_out, "w") as dst:
        src_pixc_group = src.groups['pixel_cloud']
        src_sensor_group = src.groups['tvp']

        dst_pixc_group = dst.createGroup('pixel_cloud')
        dst_sensor_group = dst.createGroup('tvp')

        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)

        copy_group(src_pixc_group, dst_pixc_group)
        copy_group(src_sensor_group, dst_sensor_group)


def copy_group(src, dst):
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
    # copy all file data except for the excluded

    for name, variable in src.variables.items():
        data = np.array(src[name][:])
        if name == 'classification':
            data[np.where(data == 0)] = 1
            data[np.where(data == 1)] = 4

        x = dst.createVariable(name, variable.datatype, variable.dimensions)
        # copy variable attributes all at once via dictionary
        dst[name].setncatts(src[name].__dict__)
        dst[name][:] = data


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
    config.set('LOGGING', "logFile", os.path.join(lake_dir, "LogFile.log"))


    # Write the parameter file
    laketile_cfg = join(lake_dir, "laketile.cfg")

    with open(laketile_cfg, 'w') as cfg_file:
        config.write(cfg_file)

    if env == "swotCNES":
        print("== Run LakeTile from swotCNES ==")
        PGE = None
        try:
            # 1 - Instantiate PGE
            # print(os.environ)
            print(os.path.abspath(inspect.getfile(pge_lake_tile)))

            PGE = pge_lake_tile.PGELakeTile(laketile_cfg)

            # 2 - Start PGE Lake Tile
            PGE.start()
        finally:
            if PGE is not None:
                # 3 - Stop PGE Lake Tile
                PGE.stop()
    else :
        print("== Run LakeTile from swot-hydrology-toolbox ==")
        PGE = None
        try:
            # 1 - Instantiate PGE
            PGE = pge_lake_tile.PGELakeTile(laketile_cfg)

            # 2 - Start PGE Lake Tile
            PGE.start()
        finally:
            if PGE is not None:
                # 3 - Stop PGE Lake Tile
                PGE.stop()

    print("== Run LakeTile OK ==")

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
    return pixcname


#######################################


def main():
    """When run as a script"""
    parser = argparse.ArgumentParser()
    parser.add_argument('river_annotation_file', help="river annotation file (output from l2pixc_to_rivertile) or directory", type=str)
    parser.add_argument('output_dir', help="output directory", type=str)
    parser.add_argument('parameter_laketile', help="parameter file", type=str)
    parser.add_argument("-gdem", action='store_true', dest='gdem', default=None, help='if true, compute tile with pixc_gdem')
    parser.add_argument(
        '-swotCNES', action='store_true', dest='swotCNES', default=None, help='only for CNES developer users, switch to swotCNES processing env')
            
    args = vars(parser.parse_args())
                  
    if args['swotCNES']:
        print("Running swotCNES environnment")
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

    if not os.path.isdir(join(args['output_dir'], "inputs")):
        os.makedirs(join(args['output_dir'], "inputs"))


    for river_annotation in river_files:
        
        print(">>>>> Dealing with river annotation file %s <<<<<" % river_annotation)
        print()
    
        # Load annotation file
        ann_rdf = my_rdf.myRdfReader(os.path.abspath(river_annotation))
        ann_cfg = ann_rdf.parameters
	
        # Prepare args for pge_lake_tile.py
        if args['gdem']:
            links_dir = join(args['output_dir'], "inputs")
            pixc_file = abspath(ann_cfg["gdem pixc file"])
            pixcvec_file = abspath(ann_cfg['pixcvec file from gdem'])
            pixc_name = os.path.basename(pixc_file)
            gdem_pixc_to_pixc(pixc_file, os.path.join(links_dir, pixc_name))
            pixc_file = os.path.join(abspath(links_dir), pixc_name)
        else :
            pixc_file    = abspath(ann_cfg['pixc file'])
            pixcvec_file = abspath(ann_cfg['pixcvec file'])

        print(". PixC file = %s" % pixc_file)
        print(". PIXCVecRiver file = %s" % pixcvec_file)
	
        # Read pixc file attributes and format with leading zeros
        with Dataset(pixc_file, "r") as pixc_dataset:
	                
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

        pixc_name = call_pge_lake_tile(parameter_laketile, args['output_dir'],
                           pixc_file, pixcvec_file,
                           cycle_number, pass_number, tile_number, start_time, stop_time, env=env)
        
        print()
        print()
        # Write annotation file

        # Create symlinks to input with the right name convention
        links_dir = join(args['output_dir'], "inputs")

        # File path
        pixc_file = os.path.join(links_dir, pixc_name)
        lake_filenames = my_names.lakeTileFilenames(IN_pixc_file=pixc_file)
        laketile_shp_file = os.path.join(args['output_dir'], lake_filenames.laketile_shp_file)
        laketile_pixcvec_file = os.path.join(args['output_dir'], lake_filenames.laketile_pixcvec_file)
        laketile_edge_file = os.path.join(args['output_dir'], lake_filenames.laketile_edge_file)
        lake_annot_file = os.path.join(args['output_dir'], lake_filenames.lake_annot_file)

        write_annotation_file(lake_annot_file, laketile_shp_file, laketile_edge_file, laketile_pixcvec_file)


    print("===== rivertile_to_laketile = END =====")


#######################################


if __name__ == "__main__":
    main()
