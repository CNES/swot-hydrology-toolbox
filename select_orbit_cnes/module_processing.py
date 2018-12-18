# -*- coding: utf8 -*-
'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


"""
.. module processing.py
    :synopsis: Basic processing class to run in simulations
    Created on 21 sept. 2012
    
.. moduleauthor: Capgemini

    $Id: module_processing.py 1093 2015-04-23 09:50:42Z nestival $
"""

from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import os

from findOrbit import findOrbit
from ressources.utils.input_reader import Input_Reader
import ressources.utils.passplan as lib_passplan


class Processing(object):

    def __init__(self, parameter_file, output_directory):
        """
        Constructor, initializes the module and read the command file
        
        Args:
            context(str) : The path to the module execution context
        """
        print("[select_orbit_cnes] == INIT ==")
        print("[select_orbit_cnes] Parameter file = %s" % parameter_file)
        print("[select_orbit_cnes] Output directory = %s" % output_directory)
        
        self.parameter_file = parameter_file
        self.output_directory = output_directory
        self.north= None
        self.south= None
        self.east= None
        self.west= None
        self.near_range = None
        self.swath_length = None
        self.orbit_directory = None
        self.start_mission_time = None
        self.makePassPlan = None
        self.simulation_start = None
        self.simulation_stop = None
        
        print()
    
    def run_preprocessing(self):
        """
        Preprocessing, commonly used to perform some checking, read the configuration
        and initialize structures that will be used during the processing.
        
        Returns:
            int. return code
        """
        print("[select_orbit_cnes] == PRE-PROCESSING... ==")
        
        # Once API initialized in the __init__ function, it is possible to call API functions
        # to get elements from the command.
        
        # To get a module input, just use its ID as follows:
        # self.api.get_input(id) => get a single file
        # self.api.get_input_list(id) => get a list of files
        #
        # The ID is the same as the one declared in the module_interface_descriptor.par file
        return_code = 0
        
        # getting parameters from the chain parameter file
        param_reader = Input_Reader(self.parameter_file)
        # reading the parameter file
        param_reader.read_param_file()

        # Find orbit directory path
        try:
            self.orbit_directory = param_reader.get_orbit_repository()
        except:
            self.api.error("orbits folder not found")
            return_code = 101
        self.north = float(param_reader.get_north_latitude())
        self.south = float(param_reader.get_south_latitude())
        self.east = float(param_reader.get_east_longitude())
        self.west = float(param_reader.get_west_longitude())

        self.near_range = float(param_reader.get_near_range())
        self.swath_length = float(param_reader.get_swath())
        self.gdem_prefix = param_reader.get_gdem_prefix()
        self.start_mission_time = param_reader.get_start_mission()
        self.makePassPlan = param_reader.get_make_pass_plan()
        self.simulation_start_time = param_reader.get_simulation_start()
        self.simulation_stop_time = param_reader.get_simulation_stop()
        self.cycle_duration = float(param_reader.get_cycle_duration())
        
        print()
        return return_code
    
    def run_processing(self):
        """
        Main process, computations are done here
        
        Returns:
            int. return code
        """
        print("[select_orbit_cnes] == PROCESSING... ==")
        
        # The log file is handled automatically
        # In order to log a message use the API corresponding functions:
        # self.api.info(message)
        # or
        # self.api.error(message)
        # There are 5 levels of messages : trace, debug, info, warning, error
       
        # If you need to create a temporary file, use the dedicated function:
        # file = self.api.create_tmp_file()
        
        # If an external binary has to be called from this source code, use
        # the following procedure :
        #     1) prepare everything the binary needs (maybe a particular input
        # in order to determine what to do ?)
        #     2) start the binary from command line:
        # import subprocess
        # command = ["<path/to/binary>", "<argument_1>", "<argument_2>", "<argument_3>", ...]
        # return_code = subprocess.call(command)
        #     3) the "call" function is blocked until the binaray execution ends
        #     4) when call ends, control the binary outputs and return code
        
        # Modify coordinates if not coherent between each others
        if self.north < self.south:
            self.north, self.south = self.south, self.north
        if self.west < self.east:
            self.west, self.east = self.east, self.west
            
        # Init orbit class
        gdem_orbit = findOrbit(self.north, self.south, self.east, self.west, self.near_range, self.swath_length)
        
        # Compute orbit files specific to studied area
        prefix = os.path.join(self.output_directory, self.gdem_prefix)
        orbit_over_dem = gdem_orbit.orbit_over_dem(self.orbit_directory, prefix, start_mission_time = self.start_mission_time)       
        
        # Compute pass plan if asked
        if self.makePassPlan == "yes":
            passplan = lib_passplan.Passplan(self.output_directory, self.start_mission_time, self.cycle_duration, self.simulation_start_time, self.simulation_stop_time)
            passplan.run_preprocessing()
            passplan.run_processing()
            
        print()
        return 0
    
    def run_postprocessing(self):
        """
        Postprocessing, at this point, the output product is written, and memory structures 
        are freed, file are closed. 
        
        Returns:
            int. return code
        """
        print("[select_orbit_cnes] == POST-PROCESSING... ==")
        
        # If an output file has to be written, use the API function to know
        # where you can write it:
        # output_location = self.api.get_output_folder()
        
        # Keep this as the last line of this function
        #self.api.end_module(True)
        
        print()
        return 0
        
