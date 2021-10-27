# -*- coding: utf8 -*-
'''
Change height of water mask for swot-hydrology-toolbox

Copyright (c) 2018, CNES
'''

import argparse
import os
import sys
import logging
import traceback

import geopandas as gpd

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
    myparser.add_argument("-v", type=float, required=True,
                          help="Height value")
    return myparser


def process_file(inputfile: str, outputfile: str, height: float):
    '''
    Change height of water mask for swot-hydrology-toolbox

    :param inputfile: Input file 
    :param outputfile: Output file
    '''

    try:
        data = gpd.read_file(inputfile)
        data.at[0,'height'] = height
        data.to_file(outputfile)
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
    process_file(args.i, args.o, args.v)

# End
