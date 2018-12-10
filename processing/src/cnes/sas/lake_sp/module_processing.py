# -*- coding: utf8 -*-
"""
.. module:: module_processing.py
    :synopsis: Process PGE_L2_HR_LakeSP, i.e. generate L2_HR_LakeSP shp and update L2_HR_PIXC_VEC NetCDF products from files of produced by PGE_L2_HR_LakeTile
                ** Main program for SAM **
    Created on 02/02/2018

.. moduleauthor: Cécile Cazals - CS

Copyright (c) 2017 CNES. All rights reserved.
"""


from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division


from sam.framework.utils.api.module_common import ModuleCommon
from swot_hr.common.rdf.rdf_reader import RdfReader
from swot_hr.common.rdf.rdf_enums import RDF_DEFAULT

import lib.my_api as my_api
import pge_lake_sp


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
        my_api.initApi(self.api, "DEBUG")  # Choose DEBUG or INFO verbose level

        # Init PGE_LakeSP object
        self.objPgeLakeSP = pge_lake_sp.Processing()

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
        my_api.printInfo("[lakeSPProcessing] PRE-PROCESSING...")
        my_api.printInfo("")

        # 1 - Read the parameter file
        my_api.printInfo("[lakeSPProcessing] > 0 - Reading parameter file...")
        try:

            # Read parameter file
            parameters = RdfReader(str(self.api.get_input("parameter_file")))

            # Get working directories
            self.objPgeLakeSP.lake_tile_dir = str(parameters.get_parameter(RDF_DEFAULT, "LakeTile directory"))
            self.objPgeLakeSP.lake_sp_dir = self.api.get_output_folder()

            # Get tile(s) information
            self.objPgeLakeSP.cycle_number = str(self.api.get_param("cycle"))  # Cycle number
            self.objPgeLakeSP.pass_number = str(self.api.get_param("orbit"))  # Orbit number

            # Get config parameters

            # Get geographical area to process
            self.objPgeLakeSP.south_lat = int(parameters.get_parameter(RDF_DEFAULT, "South latitude"))
            self.objPgeLakeSP.north_lat = int(parameters.get_parameter(RDF_DEFAULT, "North latitude"))

            # General
            self.objPgeLakeSP.lakeDb_file = str(parameters.get_parameter(RDF_DEFAULT, "Lake a priori database"))
            self.objPgeLakeSP.minSize = float(parameters.get_parameter(RDF_DEFAULT, "Min size for lake"))
            self.objPgeLakeSP.impGeoloc = int(parameters.get_parameter(RDF_DEFAULT, "Improve geolocation"))
            self.objPgeLakeSP.flagProdShp = str(parameters.get_parameter(RDF_DEFAULT, "Produce shp"))
            # Flags to process
            self.objPgeLakeSP.ice_flags = [int(flag) for flag in str(parameters.get_parameter(RDF_DEFAULT, "Ice flags to keep")).split(";")]
            self.objPgeLakeSP.layover_flags = [int(flag) for flag in str(parameters.get_parameter(RDF_DEFAULT, "Layover flags to keep")).split(";")]
            self.objPgeLakeSP.dark_water_flags = [int(flag) for flag in str(parameters.get_parameter(RDF_DEFAULT, "Dark water flags to keep")).split(";")]
            # Big lake congig
            self.objPgeLakeSP.geoloc_biglake_min_size = float(parameters.get_parameter(RDF_DEFAULT, "Geoloc biglake min size"))
            self.objPgeLakeSP.geoloc_biglake_grid_spacing = float(parameters.get_parameter(RDF_DEFAULT, "Geoloc biglake grid spacing"))
            self.objPgeLakeSP.geoloc_biglake_grid_resolution = float(parameters.get_parameter(RDF_DEFAULT, "Geoloc biglake grid resolution"))

        except IOError:
            my_api.exitWithError("[lakeSPProcessing]   Parameter file not found")
        my_api.printInfo("")

        # 2 - Run PGE_LakeSP pre-processing
        self.objPgeLakeSP.run_preprocessing()

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

        self.objPgeLakeSP.run_processing()

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

        self.objPgeLakeSP.run_postprocessing()

        # Keep this as the last line of this function
        self.api.end_module(True)
        return 0
