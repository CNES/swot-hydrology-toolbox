# -*- coding: utf8 -*-
"""
.. module select_orbit_cnes.py
    :synopsis:Basic main to run in simulations
    Created on 21 sept. 2012
    
.. moduleauthor: Capgemini

 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
 
    $Id: module_select_orbit_cnes.py 1093 2015-04-23 09:50:42Z nestival $
"""

from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import argparse

from module_processing import Processing
import ressources.utils.my_timer as my_timer


if __name__ == '__main__':
    
    # 0 - Parse inline parameters
    parser = argparse.ArgumentParser(description="Compute orbit files specific to the studied area")
    parser.add_argument("param_file", help="full path to the parameter file (*.rdf)")
    parser.add_argument("output_dir", help="full path to the output directory")
    args = parser.parse_args()
    
    print("===== select_orbit_cnes = BEGIN =====")
    print("")
    timer = my_timer.Timer()
    timer.start()
    
    try:
        
        # 1 - Init processing
        process = Processing(args.param_file, args.output_dir)
        
        # 2 - Run preprocessing
        ret = process.run_preprocessing()
            
        if (ret == 0):
            
            # 3 - Run processing
            ret = process.run_processing()
            if (ret == 0):
                # 4 - Run postprocessing
                ret = process.run_postprocessing()
                    
    except (BaseException) as bex:
        print("Uncaught exception in module processing")
        raise
        
    print(timer.stop())
    print("===== select_orbit_cnes = END =====")
    
