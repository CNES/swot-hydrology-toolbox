#!/usr/bin/env python
'''
This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

'''
import argparse
import ast
import datetime
import glob
import logging
import os
import multiprocessing as mp
from copy import deepcopy

import tools.my_filenames as my_names
import tools.my_rdf as my_rdf

import cnes.modules.geoloc.lib.pixc_to_shp

import RDF as RDF
import SWOTRiver.Estimate as Estimate
from SWOTRiver.products.pixcvec import L2PIXCVector


#######################################


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
        LOGGER.info("0 - Produce fake PIXC from GDEM")
        LOGGER.info("fake_pixc_from_gdem.fake_pixc_from_gdem(%s, %s, %s)" % (gdem_file, pixc_file, pixc_file))
        fake_pixc_from_gdem.fake_pixc_from_gdem(gdem_file, pixc_file, pixc_file)
        LOGGER.info()

    LOGGER.info("1 - [RiverObs] Init class for running RiverObs on a SWOT L2 PixelCloud data product")
    LOGGER.info("> Estimate.L2PixcToRiverTile(%s, %s)" % (pixc_file, out_pixc_vector_file))
    l2pixc_to_rivertile = Estimate.L2PixcToRiverTile(pixc_file, out_pixc_vector_file)
    LOGGER.info("- end -")

    LOGGER.info("2 - [RiverObs] Copies config object into self's storage from main")
    LOGGER.info("> Estimate.L2PixcToRiverTile.load_config(config)")
    l2pixc_to_rivertile.load_config(config)
    LOGGER.info("- end -")

    # generate empty output file on errors
    LOGGER.info("3 - Run RiverObs")
    try:
        LOGGER.info("3.1 - [RiverObs] Estimate.L2PixcToRiverTile.do_river_processing()")
        l2pixc_to_rivertile.do_river_processing()
        LOGGER.info("3.2 - [RiverObs] Estimate.L2PixcToRiverTile.match_pixc_idx()")
        l2pixc_to_rivertile.match_pixc_idx()
        LOGGER.info("3.3 - [RiverObs] Estimate.L2PixcToRiverTile.do_improved_geolocation()")
        l2pixc_to_rivertile.do_improved_geolocation()
        LOGGER.info("3.4 - [RiverObs] Estimate.L2PixcToRiverTile.flag_lakes_pixc()")
        l2pixc_to_rivertile.flag_lakes_pixc()

    except Exception as exception:
        LOGGER.error(
            'Unable to continue river processing: {}'.format(exception))
    LOGGER.info("- end -")

    LOGGER.info("4 - [RiverObs] build_products()")
    l2pixc_to_rivertile.build_products()
    LOGGER.info("- end -")

    # rewrite index file to make it look like an SDS one
    LOGGER.info("5 - [RiverObs] Generate PIXCVecRiver file")
    LOGGER.info("> L2PIXCVector.from_ncfile(%s).to_ncfile(%s)" % (l2pixc_to_rivertile.index_file, l2pixc_to_rivertile.index_file))
    L2PIXCVector.from_ncfile(l2pixc_to_rivertile.index_file
                             ).to_ncfile(l2pixc_to_rivertile.index_file)
    LOGGER.info("- end -")

    LOGGER.info("6 - [RiverObs] Generate output files")
    LOGGER.info("> Estimate.L2PixcToRiverTile.rivertile_product.to_ncfile(%s)" % out_riverobs_file)
    l2pixc_to_rivertile.rivertile_product.to_ncfile(out_riverobs_file)
    if shpbasedir is not None:
        if not os.path.isdir(shpbasedir):
            os.mkdir(shpbasedir)
        LOGGER.info("> Estimate.L2PixcToRiverTile.rivertile_product.nodes.write_shapes(%s)" % os.path.join(shpbasedir, 'nodes.shp'))
        l2pixc_to_rivertile.rivertile_product.nodes.write_shapes(
            os.path.join(shpbasedir, 'nodes.shp'))
        LOGGER.info("> Estimate.L2PixcToRiverTile.rivertile_product.reaches.write_shapes(%s)" % os.path.join(shpbasedir, 'reaches.shp'))
        l2pixc_to_rivertile.rivertile_product.reaches.write_shapes(
            os.path.join(shpbasedir, 'reaches.shp'))
    LOGGER.info("- end -")

    if gdem_file is not None:
        os.remove(pixc_file)


def run_river_tile( work):
    args, pixc_ann_file = work
    logging.info("***** Dealing with PixC annotation file %s *****" % pixc_ann_file)

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

    logging.info("== Run RiverObs ==")
    try :
        l2pixc_to_rivertile(pixc_file, output_riverobs, output_pixcvec, os.path.abspath(args['parameter_riverobs']), shpbasedir=river_dir,
                        log_level="info", gdem_file=None)
        logging.info("== Run RiverObs OK ==")
        logging.info("")
    except:
        logging.warning("== Run RiverObs NOK ==")
        logging.info("")
        return None

    logging.info("")

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
        logging.info("> Converting PIXCVecRiver .nc file to shapefile...")
        pixcvec_vars = ["azimuth_index",
                        "range_index",
                        "pixc_index",
                        "height_vectorproc",
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

    logging.info("")

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
    parser.add_argument("--verbose", help="Verbose level (debug or info=default)", nargs="?", type=str, default="info")
    parser.add_argument("--writelog", help="if true, write riverobs output to log file", nargs='?', type=bool, default=False, const=True)
    parser.add_argument(
        '-f', '--force', action='store_true', dest='force', default=False,
        help='force overwrite existing outputs; default is to quit')
    args = vars(parser.parse_args())

    # Set logger
    level = {'debug': logging.DEBUG, 'info': logging.INFO,
             'warning': logging.WARNING, 'error': logging.ERROR}[args["verbose"]]
    format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    if args["writelog"]:
        logFile = os.path.join(args['output_dir'], "RiverTile_" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + ".log")
        line_to_write = "> Log file = {}".format(logFile)
        logging.basicConfig(filename=logFile, level=level, format=format)
    else:
        line_to_write = "> No log file; print info on screen"
        logging.basicConfig(level=level, format=format)

    logging.info("===== l2pixc_to_rivertile = BEGIN =====")
    logging.info("")
    logging.info(line_to_write)
    logging.info("")

    # Load rdf files
    if os.path.isfile(args['l2pixc_annotation_file']):
        # Unite file case
        rdf_files = glob.glob(args['l2pixc_annotation_file'])
    else:
        # multi files case
        rdf_files = glob.glob(os.path.join(args['l2pixc_annotation_file'],"pixc*.rdf"))
    if len(rdf_files) == 0:
        logging.info("> NO PixC annotation file to deal with")
    else:
        logging.info("> %d PixC annotation file(s) to deal with" % len(rdf_files))
    logging.info("")
    logging.info("")

    n_cores = int(mp.cpu_count() / 2)
    pool = mp.Pool(n_cores, print('Initializing process {}'.format(os.getpid())))
    with pool:
        print('Running map')
        work = [(deepcopy(args), pixc_ann_file) for pixc_ann_file in rdf_files]
        tmp_result = pool.map(run_river_tile, work)
        print(tmp_result)
        pool.close()
        pool.join()


    logging.info("")
    logging.info("===== l2pixc_to_rivertile = END =====")


#######################################


if __name__ == "__main__":
    main()
