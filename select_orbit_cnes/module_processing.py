# -*- coding: utf8 -*-
"""
.. module_processing.py
    :synopsis: Basic processing class to run in simulations
    Created on 21 sept. 2012
    
.. moduleauthor: Capgemini

 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
 
"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import os

from find_orbit import findOrbit
from ressources.utils.input_reader import Input_Reader
import ressources.utils.passplan as lib_passplan
from ressources.utils import my_api

class Processing(object):

    def __init__(self, in_parameter_file, in_output_directory):
        """
        Constructor: initializes the self attributes
        
        :param in_parameter_file: parameter file full path
        :type in_parameter_file: str
        :param in_output_directory: output directory full path
        :type in_output_directory: str
        """
        my_api.printInfo("[select_orbit_cnes] == INIT ==")
        my_api.printInfo("[select_orbit_cnes] Parameter file = %s" % in_parameter_file)
        my_api.printInfo("[select_orbit_cnes] Output directory = %s" % in_output_directory)
        my_api.printInfo("")
        
        self.parameter_file = in_parameter_file
        self.output_directory = in_output_directory
        
        self.mission_name = None
        self.mission_start_time = None
        self.orbit_directory = None
        self.south = None
        self.north = None
        self.west = None
        self.east = None
        self.azimuth_spacing = None
        self.swath_width = None
        self.near_range = None
        self.makePassPlan = None
        self.simulation_start = None
        self.simulation_stop = None        
    
    #----------------------------------
    
    def run_preprocessing(self):
        """
        Preprocessing: read the parameter file and initialize variables
        
        :return: return code (0=OK - 101=error)
        :rtype: int
        """
        my_api.printInfo("[select_orbit_cnes] == PRE-PROCESSING... ==")
        
        return_code = 0
        
        # 1 - Read the RDF parameter file
        param_reader = Input_Reader(self.parameter_file)  # Init reader
        param_reader.read_param_file()  # Read the file

        # 2 - Find orbit directory path
        try:
            self.orbit_directory = param_reader.get_orbit_repository()
            if not os.path.exists(self.orbit_directory):
                raise FileNotFoundError("Orbit repository doesn't exist: %s" % self.orbit_directory)
        except:
            raise FileNotFoundError("Orbit repository not populated in the configuration file")
            return_code = 101
        
        # 3 - Set class attributes values
        
        self.mission_name = param_reader.get_mission_name()
        self.mission_start_time = param_reader.get_mission_start_time()        
        
        self.north = float(param_reader.get_north_latitude())
        self.south = float(param_reader.get_south_latitude())
        self.east = float(param_reader.get_east_longitude())
        self.west = float(param_reader.get_west_longitude())
        
        self.azimuth_spacing = float(param_reader.get_azimuth_spacing())
        self.near_range = float(param_reader.get_near_range())
        self.swath_width = float(param_reader.get_swath())
        self.gdem_prefix = param_reader.get_gdem_prefix()
        self.makePassPlan = param_reader.get_make_pass_plan()
        self.simulation_start_time = param_reader.get_simulation_start()
        self.simulation_stop_time = param_reader.get_simulation_stop()
        
        my_api.printInfo("")
        return return_code
    
    def run_processing(self):
        """
        Main process, computations are done here
        
        Returns:
            int. return code
        """
        my_api.printInfo("[select_orbit_cnes] == PROCESSING... ==")
        
        # 1 - Modify coordinates if not coherent between each others
        if self.north < self.south:
            self.north, self.south = self.south, self.north
        if self.west < self.east:
            self.west, self.east = self.east, self.west
            
        # 2 - Init orbit class
        gdem_orbit = findOrbit(self.south, self.north, self.west, self.east, self.swath_width, self.near_range)
        
        # 3 - Compute orbit files specific to studied area
        prefix = os.path.join(self.output_directory, self.gdem_prefix)
        cycle_duration = gdem_orbit.orbit_over_dem(self.orbit_directory, prefix, self.azimuth_spacing, self.swath_width, in_mission_start_time=self.mission_start_time)       
        
        # Compute pass plan if asked
        if self.makePassPlan == "yes":
            passplan = lib_passplan.Passplan(self.output_directory, self.mission_start_time, cycle_duration, self.simulation_start_time, self.simulation_stop_time)
            passplan.run_preprocessing()
            passplan.run_processing()
            
        my_api.printInfo("")
        return 0
    
    def run_postprocessing(self):
        """
        Postprocessing, at this point, the output product is written, and memory structures 
        are freed, file are closed. 
        
        Returns:
            int. return code
        """
        my_api.printInfo("[select_orbit_cnes] == POST-PROCESSING... ==")
        
        # If an output file has to be written, use the API function to know
        # where you can write it:
        # output_location = self.api.get_output_folder()
        
        # Keep this as the last line of this function
        #self.api.end_module(True)
        
        my_api.printInfo("")
        return 0
        
