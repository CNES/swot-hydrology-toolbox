#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
.. module sisimp.py
    :synopsis: SWOT large scale simulator

.. module author: D.Blumstein + Capgemini

    Copyright (c) 2016 CNES/LEGOS/CTOH. All rights reserved.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
import datetime
import os

import lib.my_api as my_api
import lib.my_timer as my_timer
import lib.my_tools as my_tools
import sisimp_processing as sisimp_ps


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("parameter_file", type=str)
    parser.add_argument("-v", "--verbose", help="Verbose level", nargs="?", type=str, default="INFO")
    parser.add_argument("-l", "--logfile", help="Write prints to a logfile", nargs='?', type=bool, default=False, const=True)
    args = parser.parse_args()

    print("===== sisimp = BEGIN =====")
    timer = my_timer.Timer()
    timer.start()

    # Parameter file
    my_tools.testFile(args.parameter_file)

    # Verbose level
    verbose_level = my_api.setVerbose(args.verbose)
    print("> Verbose level = {}".format(verbose_level))

    # Log file
    if args.logfile:
        logFile = os.path.join(os.path.dirname(args.parameter_file), "sisimp_" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + ".log")
        my_api.initLogger(logFile, verbose_level)
        print("> Log file = {}".format(logFile))
    else:
        print("> No log file ; print info on screen")
        print()

    # 1 - Initialisation
    mySimu = sisimp_ps.Processing()
    my_api.printInfo(timer.info(0))
    
    # 2 - Run preprocessing
    mySimu.run_preprocessing(args.parameter_file)
    my_api.printInfo(timer.info(0))
    
    # 3 - Run processing
    mySimu.run_processing()
    my_api.printInfo(timer.info(0))

    # 4 - Run post-processing
    mySimu.run_postprocessing()
    my_api.printInfo(timer.info(0))
        
    my_api.printInfo("")
    my_api.printInfo(timer.stop())
    
    # Close logger
    if args.logfile:
        my_api.closeLogger()
    
    
    print("")
    print(timer.stop())
    print("===== sisimp = END =====")