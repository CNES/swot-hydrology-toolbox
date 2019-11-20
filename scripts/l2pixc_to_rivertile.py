#!/usr/bin/env python
'''
This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

'''
import argparse
import glob
import os

import tools.my_filenames as my_names
import tools.my_rdf as my_rdf

import cnes.modules.geoloc.lib.pixc_to_shp

import os
import ast
import argparse

import logging

import RDF as RDF
import SWOTRiver.Estimate as Estimate
from SWOTRiver.products.pixcvec import L2PIXCVector

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

def l2pixc_to_rivertile(pixc_file, out_riverobs_file, out_pixc_vector_file, rdf_file, shpbasedir=None, log_level="info", gdem_file=None):
    LOGGER = logging.getLogger('swot_pixc2rivertile')

    level = {'debug': logging.DEBUG, 'info': logging.INFO,
             'warning': logging.WARNING, 'error': logging.ERROR}[log_level]
    format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=level, format=format)

    config = RDF.RDF()
    config.rdfParse(rdf_file)
    config = dict(config)

    # typecast most config values with eval since RDF won't do it for me
    # (excluding strings)
    for key in config.keys():
        if key in ['geolocation_method', 'reach_db_path', 'height_agg_method',
                   'area_agg_method']:
            continue
        config[key] = ast.literal_eval(config[key])

    if gdem_file is not None:
        import fake_pixc_from_gdem
        import tempfile
        pixc_file = tempfile.mktemp()
        fake_pixc_from_gdem.fake_pixc_from_gdem(gdem_file, pixc_file, pixc_file)

    l2pixc_to_rivertile = Estimate.L2PixcToRiverTile(pixc_file, out_pixc_vector_file)

    l2pixc_to_rivertile.load_config(config)

    # generate empty output file on errors
    try:
        l2pixc_to_rivertile.do_river_processing()
        l2pixc_to_rivertile.match_pixc_idx()
        l2pixc_to_rivertile.do_improved_geolocation()
        l2pixc_to_rivertile.flag_lakes_pixc()

    except Exception as exception:
        LOGGER.error(
            'Unable to continue river processing: {}'.format(exception))

    l2pixc_to_rivertile.build_products()

    # rewrite index file to make it look like an SDS one
    L2PIXCVector.from_ncfile(l2pixc_to_rivertile.index_file
                             ).to_ncfile(l2pixc_to_rivertile.index_file)

    l2pixc_to_rivertile.rivertile_product.to_ncfile(out_riverobs_file)
    if shpbasedir is not None:
        if not os.path.isdir(shpbasedir):
            os.mkdir(shpbasedir)
        l2pixc_to_rivertile.rivertile_product.nodes.write_shapes(
            os.path.join(shpbasedir, 'nodes.shp'))
        l2pixc_to_rivertile.rivertile_product.reaches.write_shapes(
            os.path.join(shpbasedir, 'reaches.shp'))

    if gdem_file is not None:
        os.remove(pixc_file)


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
    print()
    
    # Process per annotation file   
    for pixc_ann_file in rdf_files:
        
        print("***** Dealing with PixC annotation file %s *****" % pixc_ann_file)
        
        # Open and read RDF file
        rdf = my_rdf.myRdfReader(os.path.abspath(pixc_ann_file))
        ann_cfg = rdf.parameters
        
        # Create output directories name
        # For RiverTile products
        river_dir = os.path.abspath(os.path.join(args['output_dir'], 'rivertile'))
        if not os.path.isdir(river_dir):
            os.makedirs(river_dir)
        # For PIXCVec products
        pixcvec_dir = os.path.abspath(os.path.join(args['output_dir'], 'pixcvec'))
        if not os.path.isdir(pixcvec_dir):
            os.makedirs(pixcvec_dir)
    
        # Prepare args for l2pixc_to_rivertile.py
        # File path
        pixc_file = os.path.abspath(ann_cfg['l2pixc file'])
        river_filenames = my_names.riverTileFilenames(IN_pixc_file=pixc_file)
        output_riverobs = os.path.join(river_dir, river_filenames.rivertile_file)
        output_pixcvec = os.path.join(pixcvec_dir, river_filenames.pixc_vec_river_file)
        river_ann_file = os.path.join(args['output_dir'], river_filenames.annot_file)


        print("== Run RiverObs ==")
        l2pixc_to_rivertile(pixc_file, output_riverobs, output_pixcvec, os.path.abspath(args['parameter_riverobs']), shpbasedir=river_dir,
                            log_level="info", gdem_file=None)
        print("== Run RiverObs OK ==")

        print()
        
        # Rename the shapefiles
        # Node files
        new_nodes_files = glob.glob(os.path.join(river_dir, "nodes*"))
        for node_file in new_nodes_files:
            ext = node_file.split(".")[-1]
            new_filename = os.path.join(river_dir, river_filenames.rivertile_nodes_file+"."+ext)
            if os.path.exists(new_filename):
                os.remove(new_filename)
            os.rename(node_file, new_filename)
        # Reach files
        new_reaches_files = glob.glob(os.path.join(river_dir, "reaches*"))
        for reach_file in new_reaches_files:
            ext = reach_file.split(".")[-1]
            new_filename = os.path.join(river_dir, river_filenames.rivertile_reaches_file+"."+ext)
            if os.path.exists(new_filename):
                os.remove(new_filename)
            os.rename(reach_file, new_filename)
  
        # Convert PIXCVecRiver into shapefile if wanted
        if not args["noshp"]:
            
            # write pixcvec shapefile
            print("> Converting PIXCVecRiver .nc file to shapefile...")
            pixcvec_vars = ["azimuth_index", 
                            "range_index",
                            "pixc_index",
                            "height_vectorproc", 
                            "lake_flag",
                            "node_id",
                            "reach_id"]
            cnes.modules.geoloc.lib.pixc_to_shp.pixc_to_shp(
                output_pixcvec, 
                output_pixcvec.replace(".nc",".shp"), 
                "latitude_vectorproc", "longitude_vectorproc", 
                pixcvec_vars, group_name=None)

        # Write annotation file
        write_annotation_file(
            river_ann_file, 
            pixc_file,
            output_pixcvec)
        
        print()

    print("")
    print("===== l2pixc_to_rivertile = END =====")


#######################################


if __name__ == "__main__":
    main()
