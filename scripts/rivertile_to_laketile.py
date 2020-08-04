#!/usr/bin/env python
'''
.. module:: rivertile_to_laketile.py
    :synopsis: Run SWOT PGE_L2_HR_LakeTile over one or more tile(s) of PIXC, i.e. generate SWOT L2_HR_LakeTile
                product from one or more tile(s) of L2_HR_PIXC and L2_HR_PIXCVecRiver products
                If specified, run afterwards SWOT PGE_L2_HR_LakeSP the previous output, i.e. generate SWOT L2_HR_LakeSP
                product(s) from one or more tile(s) of L2_HR_LakeTile product
    
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''

import argparse
import configparser as cfg
import glob
from netCDF4 import Dataset
import os, sys
from os.path import join, abspath

import tools.my_rdf as my_rdf


def set_swot_cnes_env():
    print("== Run LOCNES from swotCNES ==")
    
    try:
        
        # 1 - Modify PYTHONPATH
        # 1.1 - Add src dir to PYTHONPATH
        lakesrc = os.environ['SWOT_CNES']
        sys.path.append(lakesrc)
        # 1.2 - Add PGE/lake_sp dir to PYTHONPATH
        pge_lake_sp_src = os.path.join(lakesrc, "PGE" + os.path.sep + "lake_sp")
        sys.path.append(pge_lake_sp_src)
        # 1.3 - Remove tbx directories from PYTHONPATH
        if os.environ['SWOT_HYDROLOGY_TOOLBOX']:
            for path in sys.path:
                if path.startswith(os.environ['SWOT_HYDROLOGY_TOOLBOX']):
                    sys.path.remove(path)
        
        # 2 - Import packages          
        # 2.1 - Import pge_lake_tile package
        from PGE.lake_tile import pge_lake_tile as pge_lake_tile
        # 2.2 - Import multi_lake_sp package
        import multi_lake_sp as multi_lake_sp

    except:
        print("swotCNES environment not found: please, set SWOT_CNES variable")
        pge_lake_tile = None
        multi_lake_sp = None
        exit()
        
    return pge_lake_tile, multi_lake_sp


def set_swot_hydro_env():
    print("== Run LOCNES from swot-hydrology-toolbox ==")
    
    try:
        
        # 1 - Add tbx to PYTHONPATH
        # 1.1 - Add src dir to PYTHONPATH
        try:
            tbx_path = os.environ['SWOT_HYDROLOGY_TOOLBOX']
        except:
            tbx_path = os.getcwd().replace(os.sep + "scripts", "")
        sys.path.append(tbx_path)
        # 1.2 - Add PGE/lake_sp dir to PYTHONPATH
        pge_lake_sp_src = os.path.join(tbx_path, "processing" + os.path.sep + "PGE" + os.path.sep + "lake_sp")
        sys.path.append(pge_lake_sp_src)
        
        # 2 - Import packages          
        # 2.1 - Import pge_lake_tile package
        from processing.PGE.lake_tile import pge_lake_tile as pge_lake_tile
        # 2.2 - Import multi_lake_sp package
        import multi_lake_sp as multi_lake_sp
        
    except:
        print("swot-hydrology-toolbox environment not found: please, set SWOT_HYDROLOGY_TOOLBOX variable")
        pge_lake_tile = None
        multi_lake_sp = None
        exit()
        
    return pge_lake_tile, multi_lake_sp


#######################################


def make_input_symlinks(output_dir, pixc_file, pixcvec_file, cycle_number, pass_number, tile_number, start_time, stop_time, gdem=None):
    """ Makes symlinks to pixc with the right name for locnes input """

    def swot_symlink(target, name):

        # Create symlinks to input with the right name convention
        links_dir = join(output_dir, "inputs")
        if not os.path.isdir(links_dir):
            os.makedirs(links_dir)
            
        crid = "Dx0000"
        product_counter = "01"
        outname = name.format(cycle_number, pass_number, tile_number, start_time, stop_time, crid, product_counter)
        outpath = abspath(join(links_dir, outname))

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


def rename_old_input_files(pixc_file, pixcvec_file, pixcname, flag_rename_pixc, pixcvecname, flag_rename_pixcvec):
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


def call_pge_lake_tile(cycle_number, pass_number, tile_ref, 
                       parameter_laketile, laketile_dir, pixc_file, pixcvec_file, pge_lake_tile):

    config = cfg.ConfigParser()
    config.read(parameter_laketile)

    # Fill missing values in pge_lake_tile config file
    config.set('PATHS', "PIXC file", pixc_file)
    config.set('PATHS', "PIXCVecRiver file", pixcvec_file)
    config.set('PATHS', "Output directory", laketile_dir)
    config.set('LOGGING', "errorFile", os.path.join(laketile_dir, "ErrorFile_%s_%s_%s.log" % (cycle_number, pass_number, tile_ref)))
    config.set('LOGGING', "logFile", os.path.join(laketile_dir, "LogFile_%s_%s_%s.log" % (cycle_number, pass_number, tile_ref)))

    # Write the parameter file
    laketile_cfg = join(laketile_dir, "parameter_laketile_%s_%s_%s.cfg" % (cycle_number, pass_number, tile_ref))

    with open(laketile_cfg, 'w') as cfg_file:
        config.write(cfg_file)

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


def call_pge_lake_sp(parameter_laketile, laketile_dir, lakesp_dir, multi_lake_sp):

    config = cfg.ConfigParser()
    config.read(parameter_laketile)

    # 0.1 - Fill missing values in pge_lake_sp config file
    config.remove_option("PATHS", "PIXC file")
    config.remove_option("PATHS", "PIXCVecRiver file")
    config.set('PATHS', "LakeTile directory", laketile_dir)
    config.set('PATHS', "Output directory", lakesp_dir)
    config.set('LOGGING', "errorFile", os.path.join(lakesp_dir, "ErrorFile.log"))
    config.set('LOGGING', "logFile", os.path.join(lakesp_dir, "LogFile.log"))
    config.set("FILE_INFORMATION", "CRID_LAKETILE", config.get("FILE_INFORMATION", "CRID"))
    config.set("FILE_INFORMATION", "CRID_LAKESP", config.get("FILE_INFORMATION", "CRID"))
    config.remove_option("FILE_INFORMATION", "CRID")

    # 0.2 - Write the parameter file
    lakesp_cfg = join(lakesp_dir, "parameter_lakesp.cfg")

    with open(lakesp_cfg, 'w') as cfg_file:
        config.write(cfg_file)

    # 0.3 - Read variables in command file
    my_params = multi_lake_sp.read_command_file(lakesp_cfg)
    
    # 1 - Initialization
    multi_lake_sp_proc = multi_lake_sp.MultiLakeSP(my_params)

    # 2 - Run pre-processing
    multi_lake_sp_proc.run_preprocessing()

    # 3 - Run processing
    multi_lake_sp_proc.run_processing()

    print("== Run LakeSP OK ==")
    print()


#######################################


def main():
    """When run as a script"""

    # 0 - Parse inline parameters
    parser = argparse.ArgumentParser(description="Run SWOT LakeTile processing. \
                                     If output_dir_lakesp exists, run SWOT LakeSP processing over the output of previous LakeTile processing.")
    parser.add_argument('river_annotation_file', help="river annotation file (output from l2pixc_to_rivertile) or directory", type=str)
    parser.add_argument('output_dir', help="output directory", type=str)
    parser.add_argument('parameter_laketile', help="param file", type=str)
    parser.add_argument('-output_dir_lakesp', help="output directory for LakeSP processing", type=str)
    parser.add_argument("-gdem", action='store_true', dest='gdem', default=None, help='if true, compute tile with pixc_gdem')
    parser.add_argument("-trueassign", action='store_true', dest='trueassign', default=None, help='if true, compute tile with pixc_trueassign')
    parser.add_argument(
        '-swotCNES', action='store_true', dest='swotCNES', default=None, help='only for CNES developer users, switch to swotCNES processing env; SWOT_CNES environment variable must be set to related directory')
            
    args = vars(parser.parse_args())

    print("===== rivertile_to_laketile = BEGIN =====")
    print("")

    # 1 - Set environment
    if args['swotCNES']:
        pge_lake_tile, multi_lake_sp = set_swot_cnes_env()
    else:
        pge_lake_tile, multi_lake_sp = set_swot_hydro_env()

    # 2 - Run LakeTile processing
    
    # 2.0 - Test if uniq or multi river annotation file(s)
    if os.path.isfile(args['river_annotation_file']):
        river_files = glob.glob(args['river_annotation_file'])
    else:
        river_files = glob.glob(os.path.join(args['river_annotation_file'], "river-annotation*.rdf"))
    # Print on console
    if len(river_files) == 0:
        print("> NO river annotation file to deal with")
    else:
        print("> %d river annotation file(s) to deal with" % len(river_files))
    print()

    for river_annotation in river_files:
        
        print(">>>>> Dealing with river annotation file %s <<<<<" % river_annotation)
        print()
    
        # 2.1 - Load annotation file
        ann_rdf = my_rdf.myRdfReader(os.path.abspath(river_annotation))
        ann_cfg = ann_rdf.parameters

        # 2.2 - Set pixc sources
        if args['gdem']:
            pixc_path = abspath(ann_cfg["gdem pixc file"])
            pixcvec_path = abspath(ann_cfg['pixcvec file from gdem'])
        elif args['trueassign']:
            pixc_path = abspath(ann_cfg["pixc true assign file"])
            pixcvec_path = abspath(ann_cfg['pixcvec true assign file'])
        else :
            pixc_path    = abspath(ann_cfg['pixc file'])
            pixcvec_path = abspath(ann_cfg['pixcvec file'])
	
        # 2.3 - Read pixc file attributes and format with leading zeros
        with Dataset(pixc_path, "r") as pixc_dataset:
	                
            cycle_number = "{:03d}".format(pixc_dataset.getncattr("cycle_number"))
            pass_number = "{:03d}".format(pixc_dataset.getncattr("pass_number"))
            tile_ref = pixc_dataset.getncattr("tile_name").split("_")[1]
            try:
                start_time = pixc_dataset.getncattr("time_coverage_start")
            except: 
                start_time = "YYYYMMDDThhmmss"
            try:
                stop_time = pixc_dataset.getncattr("time_coverage_end")
            except: 
                stop_time = "YYYYMMDDThhmmss"
            
            print(". Cycle number = %s" % cycle_number)
            print(". Orbit number = %s" % pass_number)
            print(". Tile ref = %s" % tile_ref)
            print(". Start time = %s" % start_time)
            print(". Stop time = %s" % stop_time)
            print()

        print("Input PixC file : %s" % pixc_path)
        print("Input PIXCVecRiver file : %s" % pixcvec_path)

        # 2.4 - Select the correct input files 
        sym_pixc_path, flag_rename_pixc, sym_pixcvec_path, flag_rename_pixcvec = make_input_symlinks(args['output_dir'], pixc_path,
                                                                                           pixcvec_path, cycle_number,
                                                                                           pass_number, tile_ref,
                                                                                           start_time, stop_time, args['gdem'])
        print("Renamed PixC file : %s" % sym_pixc_path)
        print("Renamed PIXCVecRiver file : %s" % sym_pixcvec_path)

        # 2.5 - Run LakeTile processing for current tile
        call_pge_lake_tile(cycle_number, pass_number, tile_ref,
                           args["parameter_laketile"], args['output_dir'], sym_pixc_path, sym_pixcvec_path, pge_lake_tile)

        # 2.6 - Rename old input files if necessary
        rename_old_input_files(pixc_path, pixcvec_path, sym_pixc_path, flag_rename_pixc, sym_pixcvec_path, flag_rename_pixcvec)

        print()
        print(">>>>> end-of-tile <<<<<")
        print()
        print()
        
    # 3 - Run LakeSP processing
    if args['output_dir_lakesp']:
        call_pge_lake_sp(args["parameter_laketile"], args['output_dir'], args['output_dir_lakesp'], multi_lake_sp)

    print("===== rivertile_to_laketile = END =====")


#######################################


if __name__ == "__main__":
    main()
