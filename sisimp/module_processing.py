# -*- coding: utf8 -*-
"""
.. module processing.py
    :synopsis: SAM basic processing class to run in large scale simulations
    Created on 21 sept. 2012
    
.. module author: D.Blumstein + Capgemini

    $Id: module_processing.py 1465 2016-07-01 10:05:12Z nestival $
    Copyright (c) 2016 CNES/LEGOS/CTOH. All rights reserved.
"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import numpy as np
import os
import glob
from osgeo import ogr, osr

from sam.framework.utils.api.module_common import ModuleCommon
from swot_hr.common.rdf.rdf_reader import RdfReader
from swot_hr.common.rdf.rdf_enums import RDF_DEFAULT

import lib.my_api as my_api
import proc_sisimp as sisimp

from lib.my_variables import NOISE_MULTIPLIER_FACTOR


class Processing(object):
    """
    classdocs
    """

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
        
        # Init PGE_LAKE_TILE object
        self.objSisimp = sisimp.Processing()
    
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
        my_api.printInfo("[LargeScaleSimulator] PRE-PROCESSING...")
        my_api.printInfo("")

        # Init
        parameters = None
        noise_file_path = None
        run_directory_for_orbits = None
        path_to_orbit_file = None
        
        # 1 - Read the parameter file
        my_api.printInfo("[LargeScaleSimulator] > 0 - Reading parameter file...")
        try:

            # Read parameter file
            parameters = RdfReader(str(self.api.get_input("parameter_file")))

            self.objSisimp.out_dir = self.api.get_output_folder()

            # Get tile(s) information
            self.objSisimp.cycle = int(self.api.get_param("cycle"))  # Cycle number
            self.objSisimp.orbit_number = int(self.api.get_param("orbit"))  # Orbit number

            # Instrument parameters
            # self.objSisimp.swath_width = float(parameters.get_parameter(RDF_DEFAULT, "Swath width"))
            # self.objSisimp.nr_cross_track = float(parameters.get_parameter(RDF_DEFAULT, "NR cross track"))
            # self.objSisimp.baseline = float(parameters.get_parameter(RDF_DEFAULT, "Baseline"))
            # self.objSisimp.sensor_wavelength = float(parameters.get_parameter(RDF_DEFAULT, "Sensor wavelength"))
            # self.objSisimp.range_sampling = float(parameters.get_parameter(RDF_DEFAULT, "Range sampling"))
            # self.objSisimp.nb_pix_range = int(parameters.get_parameter(RDF_DEFAULT, "Number of pixels in range"))

            # Noise parameters
            noise_file_path = parameters.get_parameter(RDF_DEFAULT, "Noise file path")
            # self.objSisimp.HEIGHT_BIAS_STD = float(parameters.get_parameter(RDF_DEFAULT, "Height bias std"))
            self.objSisimp.geolocalisation_improvement = (parameters.get_parameter(RDF_DEFAULT, "Geolocalisation improvement")).lower()

            # Orbit parameters
            run_directory_for_orbits = parameters.get_parameter(RDF_DEFAULT, "Run directory for orbits")
            # self.objSisimp.ORBIT_JITTER = float(parameters.get_parameter(RDF_DEFAULT, "Orbit jitter"))
            # self.objSisimp.CYCLE_DURATION = float(parameters.get_parameter(RDF_DEFAULT, "Cycle duration"))

            # Water bodies parameters
            self.objSisimp.shapefile_path = parameters.get_parameter(RDF_DEFAULT, "Shapefile path")
            # self.objSisimp.HEIGHT_MODEL_A = float(parameters.get_parameter(RDF_DEFAULT, "Height model A"))
            # self.objSisimp.HEIGHT_MODEL_t0 = float(parameters.get_parameter(RDF_DEFAULT, "Height model t0"))
            # self.objSisimp.HEIGHT_MODEL_PERIOD = float(parameters.get_parameter(RDF_DEFAULT, "Height model period"))

            # Water flag
            self.objSisimp.water_flag = int(parameters.get_parameter(RDF_DEFAULT, "Water flag"))

        except IOError: 
            my_api.exitWithError("[LargeScaleSimulator]   Parameter file not found")
        my_api.printInfo("")
        
        # Create shapefile  
        try:
            self.objSisimp.create_shapefile = parameters.get_parameter(RDF_DEFAULT, "Create shapefile").lower()            
            self.objSisimp.create_shapefile = self.objSisimp.create_shapefile in ['oui', 'yes', 'yep']
        except Exception: 
            self.objSisimp.create_shapefile = False
            my_api.printInfo("No Create shapefile parameter set, no shapefile will be created")
        
        # Create dummy L2_HR_PIXC_VEC_RIVER product, associated to pixel cloud  
        try:
            self.objSisimp.create_pixc_vec_river = parameters.get_parameter(RDF_DEFAULT, "Create dummy pixc vec river file").lower()            
            self.objSisimp.create_pixc_vec_river = self.objSisimp.create_pixc_vec_river in ['oui', 'yes', 'yep']
        except Exception: 
            self.objSisimp.create_pixc_vec_river = False
            my_api.printInfo("No Create dummy pixc vec river file parameter set, no L2_HR_PIXC_VEC_RIVER file will be created")
        
        # Load the noise tab    
        try:
            self.objSisimp.noise_height = np.loadtxt(noise_file_path, skiprows=1)
            # noise_multiplier = float(parameters.get_parameter(RDF_DEFAULT, "Noise multiplier factor"))
            self.objSisimp.noise_height[:, 1] = NOISE_MULTIPLIER_FACTOR * self.objSisimp.noise_height[:, 1]
        except IOError:
            my_api.exitWithError("Noise file %s not found" % noise_file_path)
       
        # Build the orbit_path to get the orbit_file   
        try:
            path_to_orbit_file = os.path.join(run_directory_for_orbits, "makeGdemOrbit_*/outputs/")
            path_to_orbit_file += "*pass_" + str(self.objSisimp.orbit_number).zfill(4) + ".nc"
            orbit_file = glob.glob(path_to_orbit_file)[0]
            self.objSisimp.read_orbit(orbit_file)
        except IndexError:
            my_api.printError("Orbit file not found = %s" % path_to_orbit_file)
            my_api.exitWithError("Please check that orbit number %d has been generated" % self.objSisimp.orbit_number)

        # Check input shapefile
        # Check if all shapefiles are there   
        file_missing = False
        if not os.path.isfile(self.objSisimp.shapefile_path + ".dbf"):
            my_api.printError("The file " + self.objSisimp.shapefile_path + ".dbf is missing.")
            file_missing = True
        if not os.path.isfile(self.objSisimp.shapefile_path + ".shp"):
            my_api.printError("The file " + self.objSisimp.shapefile_path + ".shp is missing.")
            file_missing = True
        if not os.path.isfile(self.objSisimp.shapefile_path + ".shx"):
            my_api.printError("The file " + self.objSisimp.shapefile_path + ".shx is missing.")
            file_missing = True
        if file_missing:
            raise IOError("One or several shapefile files are missing, check logs to know which one")
        # Loading shapefile
        driver = ogr.GetDriverByName(str("ESRI Shapefile"))         
        da_shape_file = driver.Open(self.objSisimp.shapefile_path + ".shp", 0)  # 0 means read-only. 1 means writeable.
        wb_layer = da_shape_file.GetLayer()
        # Check if the informations in the shapefile are right   
        shp_srs = wb_layer.GetSpatialRef()
        lonlat_srs = osr.SpatialReference()
        lonlat_srs.ImportFromEPSG(4326)
        if not lonlat_srs.IsSame(shp_srs):
            raise IOError("This is not a shapefile in lon/lat WGS84 projection")
        # self.objSisimp.compute_pixc_vec_river to True only if self.objSisimp.create_pixc_vec_river is True and RIV_FLAG field is here
        self.objSisimp.compute_pixc_vec_river = False
        if self.objSisimp.create_pixc_vec_river and (wb_layer.FindFieldIndex(str("RIV_FLAG"), True) != -1):
            self.objSisimp.compute_pixc_vec_river = True

        # # Random seed
        # try:
        #     self.objSisimp.random_seed = int(parameters.get_parameter(RDF_DEFAULT, "Random seed"))
        # except Exception:
        #     my_api.printInfo("No random seed parameter set, used default parameter instead")

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
        
        self.objSisimp.run_processing()
        
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
        
        self.objSisimp.run_postprocessing()
        
        # Keep this as the last line of this function
        self.api.end_module(True)
        return 0
