# -*- coding: utf8 -*-
'''
Print height

Copyright (c) 2018, CNES
'''

import argparse
import os
import sys
import logging
import traceback
import shapefile

def get_parser():
    """
    Parser from command line arguments

    :return myparser: Parser
    :rtype: ArgumentParser

    """
    myparser = argparse.ArgumentParser(os.path.basename(__file__))
    myparser.add_argument("-i", type=str, required=True,
                          help="Water mask file")
    return myparser


# Main program
if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    level = getattr(logging, "INFO")
    logging.basicConfig(filename=None,format='%(asctime)s [%(levelname)s] %(message)s', level=level)
    heights = []
    with shapefile.Reader(args.i) as dbf:
        for record in dbf.records():
            heights.append(record.as_dict()['height'])
    print(f"File {args.i}: {heights[0]}")

# End
