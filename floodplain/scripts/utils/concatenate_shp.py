# -*- coding: utf8 -*-
'''
Concatenate a list of shapefiles

Copyright (c) 2018, CNES
'''

import argparse
import os
import sys
import logging
import traceback
import subprocess
import glob
import pathlib
import geopandas as gpd
import pandas as pd

import floodplain.io.shp as shp


def get_parser():
    """
    Parser from command line arguments

    :return myparser: Parser
    :rtype: ArgumentParser

    """
    myparser = argparse.ArgumentParser(os.path.basename(__file__))
    myparser.add_argument("--dir", type=str, required=True,
                          help="Parent directory")
    myparser.add_argument("--out", type=str, required=True,
                          help="Output directory")
    return myparser


# Main program
if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    level = getattr(logging, "INFO")
    logging.basicConfig(filename=None,format='%(asctime)s [%(levelname)s] %(message)s', level=level)

    shpfiles = list(pathlib.Path(args.dir).glob('*.shp'))
    frames = [gpd.read_file(shpfile) for shpfile in shpfiles]
    result = pd.concat(frames)
    shp.gdf_to_file(args.out,result)                                                                                                                     
# End
