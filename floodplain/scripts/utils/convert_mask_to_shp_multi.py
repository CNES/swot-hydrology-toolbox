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
import subprocess


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

    cmd = f"find {args.dir} -maxdepth 1 -type d"
    directories = subprocess.check_output(cmd, shell=True)
    directories = directories.decode('ascii', 'ignore').split('\n')[1:-2]
    os.makedirs(args.out)

    for directory in directories:
        name = os.path.basename(directory)
        inputfilename = os.path.join(args.dir, name, '0_water_mask', name + ".shp")
        #os.makedirs(os.path.join(args.out,name, '0_water_mask')
        #outputfilename = os.path.join(args.out, name, '0_water_mask', name + ".shp") 
        outputfilename = os.path.join(args.out, name + ".shp") 
        cmd = f"python {os.path.join(os.path.dirname(__file__),'convert_mask_to_shp.py')} -i {inputfilename} -o {outputfilename}"
        subprocess.check_call(cmd, shell=True)

# End
