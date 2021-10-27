# -*- coding: utf8 -*-
'''
Convert shapefile to ply file

Copyright (c) 2018, CNES
'''

import argparse
import os
import sys
import logging
import traceback
import shapefile
import utm
import pandas as pd
import geopandas as gpd

import floodplain.io.ply as ply


def get_parser():
    """
    Parser from command line arguments

    :return myparser: Parser
    :rtype: ArgumentParser

    """
    myparser = argparse.ArgumentParser(os.path.basename(__file__))
    myparser.add_argument("-i", type=str, required=True,
                          help="Shapefile")
    myparser.add_argument("-o", type=str, required=True,
                          help="Ply file")
    return myparser


def process_file(inputfile: str, outputfile: str):
    '''
    Convert shapefile to ply file

    :param inputfile: Input file 
    :param outputfile: Output file
    '''

    try:
        data = gpd.read_file(inputfile)
        # Convert to utm
        lambdafunc = lambda row: pd.Series([*utm.from_latlon(row['latitude'],
                                                   row['longitude'])[0:2],
                                  row['height']]) 
        data[['x','y','z']] = data.apply(lambdafunc,axis=1)
        ply.gdf_to_file(outputfile,data)
    except:
        traceback.print_exc()
        logging.error(f"Writing output file")
    else:
        logging.info("Write output file: OK")


# Main program
if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    level = getattr(logging, "INFO")
    logging.basicConfig(filename=None,format='%(asctime)s [%(levelname)s] %(message)s', level=level)
    process_file(args.i, args.o)

# End
