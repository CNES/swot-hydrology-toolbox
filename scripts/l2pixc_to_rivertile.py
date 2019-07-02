#!/usr/bin/env python
'''
This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

'''
import argparse
import glob
import netCDF4 as nc
import numpy as np
import os
import subprocess
import sys

import my_rdf
import cnes.modules.geoloc.lib.pixc_to_shp

#~ from cnes.sas.lib import base_classes


def write_annotation_file(ann_file, 
                          pixc_file,
                          pixcvec_file,
                          pixcvec_file_gdem = None):
    """
    write the river-annotation.rdf file so that lake processor can run
    """
    f = open(ann_file,'w')
    f.write("pixc file = %s\n"%pixc_file)
    f.write("pixcvec file = %s\n"%pixcvec_file)
    if pixcvec_file_gdem is not None:
        f.write("pixcvec file from gdem = %s\n"%pixcvec_file_gdem)
    f.close()


#######################################
    

def execute(cmd):
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='') # process line here

    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, p.args)


def main():
    """When run as a script"""
    parser = argparse.ArgumentParser()
    parser.add_argument('l2pixc_annotation_file', help="PixC annotation file", type=str)
    parser.add_argument('output_dir', help="output directory", type=str)
    parser.add_argument('parameter_riverobs', help="param file for RiverObs", type=str)
    parser.add_argument('--riverobs_path', help="full path to RiverObs software (if not in os.environ)", type=str)
    parser.add_argument("--noshp", 
        help="if true, don't produce shapefiles", 
        nargs='?', type=bool, default=False, const=True)
    parser.add_argument(
        '-f', '--force', action='store_true', dest='force', default=False,
        help='force overwrite existing outputs; default is to quit')
    args = vars(parser.parse_args())

    print("===== l2pixc_to_rivertile = BEGIN =====")
    print("")

    # Load rdf files
    if os.path.isfile(args['l2pixc_annotation_file']):
        # Unite file case
        rdf_files = glob.glob(args['l2pixc_annotation_file'])
    else:
        # multi files case
        rdf_files = glob.glob(os.path.join(args['l2pixc_annotation_file'],"pixc*.rdf"))
    if len(rdf_files) == 0:
        print("> NO PixC annotation file to deal with")
    else:
        print("> %d PixC annotation file(s) to deal with" % len(rdf_files))
    print()
   
    for pixc_ann_file in rdf_files:
        
        print("> Dealing with PixC annotation file %s" % pixc_ann_file)
        
        # Open and read RDF file
        rdf = my_rdf.myRdfReader(os.path.abspath(pixc_ann_file))
        ann_cfg = rdf.parameters
        
        # Get value of orbit cycle, orbit number and latitude tile coordinates
        
        try:
            pixc_ann_split = os.path.basename(pixc_ann_file).split(".")[0].split("_")
            num_cycle = pixc_ann_split[2]
            num_orbit = pixc_ann_split[3]
            lat_tile = pixc_ann_split[4]
            print("Cycle number = %s - Orbit number = %s - Tile ref = %s" % (num_cycle, num_orbit, lat_tile))
            print()
        except:
            print("No information in pixc annotation file, assuming num_cycle = 0, num_orbit=0, lat_tile=0")
            num_cycle = '0'
            num_orbit = '0'
            lat_tile = '0'
            
                        
        #~ ann_cfg = rdf.parse(
            #~ os.path.abspath(args['l2pixc_annotation_file']), comment='!')
    
        # Clone the tail of proc directory structure
        river_dir = os.path.abspath(os.path.join(args['output_dir'], 'pixc'))
        if not os.path.isdir(river_dir):
            os.makedirs(river_dir)
    
        # Prepare args for l2pixc_to_rivertile.py
        # File path
        pixc_file = os.path.abspath(ann_cfg['l2pixc file'])
        output_riverobs = os.path.join(river_dir, "rivertile_" + num_cycle + "_" + num_orbit + "_" + lat_tile + ".nc")
        output_pixcvec = os.path.join(river_dir, "pixcvec_" + num_cycle + "_" + num_orbit + "_" + lat_tile + ".nc")
        river_ann_file = os.path.join(args['output_dir'],'river-annotation_' + num_cycle + '_' + num_orbit + '_' + lat_tile + '.rdf')
        
        # Make rivertile sas with pixelcloud
        # Path to RiverObs
        if args['riverobs_path'] is not None:
            path_to_riverobs = args['riverobs_path']
        else:
            path_to_riverobs = os.environ['RIVEROBS']
        # Build RiverObs main lib
        prog_name = os.path.join(
            path_to_riverobs,
            'src','bin','swot_pixc2rivertile.py')
        # Build command
        cmd = "{} {} {} {} {}".format(prog_name,
                                         pixc_file,
                                         output_riverobs,
                                         output_pixcvec,
                                         os.path.abspath(args['parameter_riverobs']))
        # Excute command
        print ("> executing:", cmd) 
        #subprocess.check_call(cmd, shell=True)
        execute(cmd)
        print("== Execution OK")
        print()
  
        # Make the desired shape files
        if not args["noshp"]:
            print("> Converting .nc files to shapefiles...")
            pixcvec_vars = ["height_vectorproc", 
                           "node_index", 
                           "reach_index", 
                           "azimuth_index", 
                           "range_index"]
            node_vars = ["time",
                         "wse",
                         "width",
                         "area_total",
                         "node_id", 
                         "reach_id"]
            reach_vars = ["reach_id", "width", "wse", "slope", "area_total"]

            # write node shape file
            print(". Node shapefile")
            cnes.modules.geoloc.lib.pixc_to_shp.pixc_to_shp(
                output_riverobs, 
                output_riverobs.replace(".nc","_node.shp"), 
                "latitude", "longitude", node_vars, group_name="nodes")

            # write reach shape file
            print(". Reach shapefile")
            cnes.modules.geoloc.lib.pixc_to_shp.pixc_to_shp(
                output_riverobs, 
                output_riverobs.replace(".nc","_reach.shp"), 
                "p_latitud", "p_longitud", reach_vars, group_name="reaches")

            # write pixcvec shape file
            print(". PIXCVec shapefile")
            cnes.modules.geoloc.lib.pixc_to_shp.pixc_to_shp(
                output_pixcvec, 
                output_pixcvec.replace(".nc",".shp"), 
                "latitude_vectorproc", "longitude_vectorproc", 
                pixcvec_vars, group_name=None)

        
        write_annotation_file(
            river_ann_file, 
            pixc_file,
            output_pixcvec)

    print("")
    print("===== l2pixc_to_rivertile = END =====")


#######################################


if __name__ == "__main__":
    main()
