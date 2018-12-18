# -*- coding: utf8 -*-
"""
.. module processing.py
    :synopsis: Process PGE_L2_HR_LakeTile, i.e. generate L2_HR_LakeTile product from one tile of L2_HR_PIXC product and associated L2_HR_PIXCVec product
                ** Main program for SAM environment **
    02/27/2017 - Creation

.. module author: Claire POTTIER - CNES DSO/SI/TR

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

from sam.framework.utils.api.module_common import ModuleCommon
from swot_hr.common.rdf.rdf_reader import RdfReader
from swot_hr.common.rdf.rdf_enums import RDF_DEFAULT

import lib.my_api as my_api
import pge_lake_tile


class Processing(object):

    def __init__(self, context):
        """
        Constructor, initializes the module and read the command file
        
        Args:
            context(str) : The path to the module execution context
        """
        # DO NOT delete this line, API will not be usable if self.api is
        # not initialized properly.
        self.api = ModuleCommon(context)
        my_api.initApi(self.api, "DEBUG") # Choose DEBUG or INFO verbose level
        
        # Init PGE_LAKE_TILE object
        self.objPgeLakeTile = pge_lake_tile.Processing()
    
    def run_preprocessing(self):
        """
        Preprocessing, commonly used to perform some checking, read the configuration
        and initialize structures that will be used during the processing.
        
        Returns:
            int. return code
        """
        
        # Once API initialized in the __init__ function, it is possible to call API functions
        # to get elements from the command.
        
        # To get a module input, just use its ID as follows:
        # self.api.get_input(id) => get a single file
        # self.api.get_input_list(id) => get a list of files
        #
        # The ID is the same as the one declared in the module_interface_descriptor.par file
        
        return_code = 0
        
        my_api.printInfo("")
        my_api.printInfo("[lakeTileProcessing] PRE-PROCESSING...")
        my_api.printInfo("")
        
        # 1 - Read the parameter file
        my_api.printInfo("[lakeTileProcessing] > 0 - Reading parameter file...")
        try:
            
            # Read parameter file
            parameters = RdfReader(str(self.api.get_input("parameter_file")))
            
            # Get working directories
            self.objPgeLakeTile.pixc_main_dir = str(parameters.get_parameter(RDF_DEFAULT, "PIXC_main directory"))
            self.objPgeLakeTile.pixc_sensor_dir = str(parameters.get_parameter(RDF_DEFAULT, "PIXC_sensor directory"))
            self.objPgeLakeTile.pixc_vec_river_dir = str(parameters.get_parameter(RDF_DEFAULT, "PIXCVecRiver directory"))
            self.objPgeLakeTile.output_dir = self.api.get_output_folder()
            
            # Get tile(s) information
            self.objPgeLakeTile.cycle = str(self.api.get_param("cycle")) # Cycle number
            self.objPgeLakeTile.orbit = str(self.api.get_param("orbit")) # Orbit number
            self.objPgeLakeTile.tile = str(self.api.get_param("tile")) # Tile id

            # Get config parameters
            # General
            self.objPgeLakeTile.lakeDb_file = str(parameters.get_parameter(RDF_DEFAULT, "Lake a priori database"))
            self.objPgeLakeTile.minSize = float(parameters.get_parameter(RDF_DEFAULT, "Min size for lake"))
            self.objPgeLakeTile.impGeoloc = int(parameters.get_parameter(RDF_DEFAULT, "Improve geolocation"))
            self.objPgeLakeTile.hullMethod = int(parameters.get_parameter(RDF_DEFAULT, "Hull computation"))
            self.objPgeLakeTile.flagProdShp = int(parameters.get_parameter(RDF_DEFAULT, "Produce shp"))
            self.objPgeLakeTile.minSize = int(parameters.get_parameter(RDF_DEFAULT, "Min size for lake"))
            # Flags to process
            self.objPgeLakeTile.classif_flags = str(parameters.get_parameter(RDF_DEFAULT, "Classif flags to keep"))
            self.objPgeLakeTile.ice_flags = [int(flag) for flag in str(parameters.get_parameter(RDF_DEFAULT, "Ice flags to keep")).split(";")]
            self.objPgeLakeTile.layover_flags = [int(flag) for flag in str(parameters.get_parameter(RDF_DEFAULT, "Layover flags to keep")).split(";")]
            self.objPgeLakeTile.dark_water_flags = [int(flag) for flag in str(parameters.get_parameter(RDF_DEFAULT, "Dark water flags to keep")).split(";")]
            # Big lake congig
            self.objPgeLakeTile.geoloc_biglake_min_size = float(parameters.get_parameter(RDF_DEFAULT, "Geoloc biglake min size"))
            self.objPgeLakeTile.geoloc_biglake_grid_spacing = float(parameters.get_parameter(RDF_DEFAULT, "Geoloc biglake grid spacing"))
            self.objPgeLakeTile.geoloc_biglake_grid_resolution = float(parameters.get_parameter(RDF_DEFAULT, "Geoloc biglake grid resolution"))

        except IOError: 
            my_api.exitWithError("[lakeTileProcessing]   Parameter file not found")
        my_api.printInfo("")
        
        # 2 - Run PGE_LakeTile pre-processing
        self.objPgeLakeTile.run_preprocessing()
        
        return return_code
    
    def run_processing(self):
        """
        Main process, computations are done here
        
        Returns:
            int. return code
        """
        
        # The log file is handled automatically
        # In order to log a message use the API corresponding functions:
        # my_api.printInfo(message)
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
        
        self.objPgeLakeTile.run_processing()
        
        return 0
    
    def run_postprocessing(self):
        """
        Postprocessing, at this point, the output product is written, and memory structures 
        are freed, file are closed. 
        
        Returns:
            int. return code
        """
        
        # If an output file has to be written, use the API function to know
        # where you can write it:
        # output_location = self.api.get_output_folder()
        
        self.objPgeLakeTile.run_postprocessing()
        
        # Keep this as the last line of this function
        self.api.end_module(True)
        return 0
