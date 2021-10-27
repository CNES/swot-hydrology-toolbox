# -*- coding: utf8 -*-
'''
Convert ply file to shapefile

Copyright (c) 2018, CNES
'''

import argparse
import os
import sys
import logging
import traceback
import shapefile
import utm
import numpy as np
from plyfile import PlyData, PlyElement

import floodplain.io.shp as shp


def get_parser():
    """
    Parser from command line arguments

    :return myparser: Parser
    :rtype: ArgumentParser

    """
    myparser = argparse.ArgumentParser(os.path.basename(__file__))
    myparser.add_argument("-i", type=str, required=True,
                          help="Ply file")
    myparser.add_argument("-o", type=str, required=True,
                          help="Shapefile")
    return myparser


def process_file(inputfile: str, outputfile: str):
    '''
    Convert ply file to shapefile

    :param inputfile: Input file 
    :param outputfile: Output file
    '''

    try:
        plydata = PlyData.read(inputfile) 
        nb = plydata.elements[0].count
        comment = plydata.comments[0].split("UTM")[-1].strip()
        zone_number = int(comment[0:-1])
        zone_letter = comment[-1]
        points_utm = np.ndarray((nb,3))
        points_utm[:,0] = plydata['vertex']['x']
        points_utm[:,1] = plydata['vertex']['y']
        points_utm[:,2] = plydata['vertex']['z']
        # Convert to latlon
        points = [ [*utm.to_latlon(point_utm[0],point_utm[1],zone_number,zone_letter),point_utm[2]] for point_utm in points_utm]
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
