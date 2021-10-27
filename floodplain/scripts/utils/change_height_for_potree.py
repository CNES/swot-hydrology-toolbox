# -*- coding: utf8 -*-
'''
Change height in plyfile for visualisation

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

from floodplain.io.ply import MyPlyElement


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
                          help="Output file")
    myparser.add_argument("-f", type=int, required=True,
                          help="Factor")
    return myparser


def process_file(inputfile: str, outputfile: str, factor: int):
    '''
    Mutiply height by a factor

    :param inputfile: Input file 
    :param outputfile: Output file
    :param factor: Multiplicative factor
    '''

    try:
        plydata = PlyData.read(inputfile) 
        nb = plydata.elements[0].count
        comment = plydata.comments[0].split("UTM")[-1].strip()
        zone_number = int(comment[0:-1])
        zone_letter = comment[-1]
        points = np.ndarray((nb,3))
        points[:,0] = plydata['vertex']['x']
        points[:,1] = plydata['vertex']['y']
        points[:,2] = plydata['vertex']['z']
        points[:,2] = [ h * factor for h in points[:,2]]
        #ply.points_to_file(outputfile,points)
        # Utm coordinates to numpy array
        vertex = np.array([tuple(value) for value in points],
                dtype=[('x', 'f4'), ('y', 'f4'), ('z', 'f4')])
        el = MyPlyElement.describe(vertex, 'vertex',
                    comments=[f'projection: UTM {zone_number}{zone_letter}'])
        PlyData([el], text=True).write(outputfile)
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
    process_file(args.i, args.o, args.f)

# End
