# -*- coding: utf8 -*-
'''
Convert water mask for swot-hydrology-toolbox to shapefile

Copyright (c) 2018, CNES
'''

import argparse
import os
import sys
import logging
import traceback
import shapefile

import floodplain.io.shp as shp


def get_parser():
    """
    Parser from command line arguments

    :return myparser: Parser
    :rtype: ArgumentParser

    """
    myparser = argparse.ArgumentParser(os.path.basename(__file__))
    myparser.add_argument("-i", type=str, required=True,
                          help="Water mask file")
    myparser.add_argument("-o", type=str, required=True,
                          help="Output file")
    return myparser


def process_file(inputfile: str, outputfile: str):
    '''
    Create water mask for swot-hydrology-toolbox to shapefile

    :param inputfile: Input file 
    :param outputfile: Output file
    '''

    try:
        polygons = shp.from_file(inputfile)
        heights = []
        with shapefile.Reader(inputfile) as dbf:
            for record in dbf.records():
                heights.append(record.as_dict()['height'])
        points = []
        for polygon, height in zip(polygons,heights):
            points += [(*coord, height) for coord in polygon.exterior.coords]
        shp.points_to_file(outputfile,points)
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
