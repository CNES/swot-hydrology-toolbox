#!/usr/bin/env python
# -*- coding : utf8 -*-
"""
module sisimp_processing.py
module author Capgemini
This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import os
from osgeo import osr, ogr
import re
import zipfile

import lib.my_api as my_api
import lib.my_filenames as my_names
import lib.my_passplan as my_plan
import lib.my_rdf_file as my_rdf
import lib.my_tools as my_tools

import sisimp_function as sisimp_fct
from write_polygons import orbitAttributes

from lib.my_variables import SWATH_WIDTH, NR_CROSS_TRACK, SENSOR_WAVELENGTH, NB_PIX_RANGE, ORBIT_JITTER, HEIGHT_MODEL_A, \
        HEIGHT_MODEL_t0, HEIGHT_MODEL_PERIOD, HEIGHT_MODEL_MIN_AREA, HEIGHT_BIAS_STD, NOISE_MULTIPLIER_FACTOR, RANGE_SAMPLING, \
        WATER_FLAG, MULTIPLE_ORBIT, COEFF_X2, COEFF_Y2, COEFF_X, COEFF_Y, COEFF_XY, COEFF_CST, GEOLOCATION_IMPROVEMENT, \
        FACT_ECHELLE, HEIGHT_MODEL, HEIGHT_MODEL_STDV, GEN_APPROX_RAD_EARTH, RAD2DEG, DEG2RAD, FACT_ECHELLE_DW, DW_PERCENT, DARKWATER_FLAG, \
        SCALE_FACTOR_NON_DETECTED_DW, DW_DETECTED_PERCENT, DW_DETECTED_NOISE_FACTOR, DW_CORRELATION_LENGTH


def read_parameter(IN_rdf_reader, IN_instrument_name, IN_instrument_default_value, read_type):
    try:
        OUT_instrument_param = read_type(IN_rdf_reader.getValue(IN_instrument_name))
        my_api.printInfo("%s : %s" % (IN_instrument_name ,str(OUT_instrument_param)))
    except:
        OUT_instrument_param = IN_instrument_default_value
        my_api.printInfo("Default value for %s : %s" % (IN_instrument_name ,str(OUT_instrument_param)))
    return OUT_instrument_param

class Processing(object):

    def __init__(self):
        """Constructor, initializes the module and read the command file"""
        my_api.printInfo("[sisimp_processing] == INIT ==")

        self.my_attributes = orbitAttributes()
        
    #--------------------------------------------

    def run_preprocessing(self, IN_paramFile):
        """
        Preprocessing, commonly used to perform some checking, read the configuration
        and initialize structures that will be used during the processing.

        :param IN_paramFile: parameter filename ; file in RDF format (<key> = <value>)
        :type IN_paramFile: string
        """
        my_api.printInfo("")
        my_api.printInfo("")
        my_api.printInfo("[sisimp_processing] PRE-PROCESSING...")
        my_api.printInfo("")

        # Init param
        parameters = None
        run_directory_for_orbits = None
        path_to_orbit_file = None

        # Read param file
        try:
            
            # Get the parameters
            parameters = my_rdf.myRdfReader(IN_paramFile)

            # Directory of orbit files (computed by select_orbit_cnes)
            run_directory_for_orbits = os.path.expandvars(parameters.getValue("Run directory for orbits"))
            my_tools.testDir(run_directory_for_orbits)

            # Shapefile of water bodies
            self.my_attributes.shapefile_path = os.path.expandvars(parameters.getValue("Shapefile path"))
            
            # Output directory
            self.my_attributes.out_dir = os.path.expandvars(str(parameters.getValue("Output directory")))
            # Create output dir if doesn't exist
            if os.path.exists(self.my_attributes.out_dir):
                if not os.path.isdir(self.my_attributes.out_dir):
                    my_api.exitWithError("ERROR = %s is not a directory" % self.my_attributes.out_dir)
            else:
                os.makedirs(self.my_attributes.out_dir)
                
            # Cross-over residual roll error
            try:
                self.my_attributes.roll_repo = str(parameters.getValue("roll_repository_name"))
                my_api.printInfo("Roll repository : %s " % self.my_attributes.roll_repo)
            except:
                my_api.printInfo("roll_repo_name not set, roll error won't be applied")

            # Orbit parameters
            self.my_attributes.multi_orbit_option = read_parameter(parameters, "Multiple orbit", MULTIPLE_ORBIT, str).lower()
            if self.my_attributes.multi_orbit_option not in ['no', 'yes', 'passplan']:
                self.my_attributes.multi_orbit_option = 'yes'

            if self.my_attributes.multi_orbit_option == 'no':
                try:
                    self.my_attributes.orbit_number = int(parameters.getValue("Orbit"))
                    my_api.printInfo("Orbit number : %d" % self.my_attributes.orbit_number)
                except:
                     my_api.exitWithError("Multiple orbit = no => Orbit number should be set")

                try:
                    self.my_attributes.cycle_number = int(parameters.getValue("Cycle number"))
                    my_api.printInfo("Orbit number : %d" % self.my_attributes.orbit_number)
                except:
                    my_api.printInfo("Multiple orbit = no => Cycle number should be set. Set to default value : 1")
                    self.my_attributes.cycle_number = 1
                     
                                 
            if self.my_attributes.multi_orbit_option == 'passplan':
                try:
                    self.my_attributes.passplan_path = str(parameters.getValue("Passplan path"))
                    my_api.printInfo("Passplan path : %s" % self.my_attributes.passplan_path)
                except:
                    try:
                        self.my_attributes.passplan_path = os.path.join(run_directory_for_orbits, "passplan.txt")
                    except:
                        my_api.exitWithError("Multiple orbit = passplan => no passplan in orbit repo")
                        
            # Simulation parameters
            self.my_attributes.swath_width = read_parameter(parameters, "Swath width", SWATH_WIDTH, float) 
            self.my_attributes.nr_cross_track = read_parameter(parameters, "NR cross track", NR_CROSS_TRACK, float)
            self.my_attributes.sensor_wavelength = read_parameter(parameters, "Sensor wavelength", SENSOR_WAVELENGTH, float)
            self.my_attributes.range_sampling = read_parameter(parameters, "Range sampling", RANGE_SAMPLING, float)
            self.my_attributes.nb_pix_range = read_parameter(parameters, "Number of pixels in range", NB_PIX_RANGE, int)
            self.my_attributes.orbit_jitter = read_parameter(parameters, "Orbit jitter", ORBIT_JITTER, float)

            # Height model parameter
            self.my_attributes.height_model = read_parameter(parameters, "Height model", HEIGHT_MODEL, str)
    
            # True height file
            try:
                self.my_attributes.trueheight_file = os.path.expandvars(parameters.getValue("True height file"))
                my_api.printInfo("True height file : %s" % self.my_attributes.trueheight_file)
            except:
                self.my_attributes.trueheight_file = None
                my_api.printInfo("True height file not set, True height model won't be applied")
            # Dark water
            self.my_attributes.dark_water = read_parameter(parameters, "Dark water",'No' , str)
            self.my_attributes.dw_pourcent = read_parameter(parameters, "Dark water percentage",DW_PERCENT,float)
            self.my_attributes.darkwater_flag = read_parameter(parameters, "Dark water flag", DARKWATER_FLAG, int)
            self.my_attributes.dw_seed=read_parameter(parameters, "Dark water seed", None, float)
            self.my_attributes.scale_factor_non_detected_dw = read_parameter(parameters, "Scale factor non detected dw", SCALE_FACTOR_NON_DETECTED_DW,float)
            self.my_attributes.dw_detected_percent = read_parameter(parameters, "Dark water detected percentage", DW_DETECTED_PERCENT, float)
            self.my_attributes.dw_detected_noise_factor = read_parameter(parameters, "Dark water detected noise factor", DW_DETECTED_NOISE_FACTOR, float)
            self.my_attributes.dw_correlation_length = read_parameter(parameters, "Dark water correlation length", DW_CORRELATION_LENGTH, int)

            # Water flag
            self.my_attributes.water_flag = read_parameter(parameters, "Water flag", WATER_FLAG, float)

            # Noise parameters
            self.my_attributes.height_bias_std = read_parameter(parameters, "Height bias std", HEIGHT_BIAS_STD, float)
            self.my_attributes.noise_multiplier_factor = read_parameter(parameters, "Noise multiplier factor", NOISE_MULTIPLIER_FACTOR, float)
            self.my_attributes.geolocalisation_improvement = read_parameter(parameters, "Geolocalisation improvement", GEOLOCATION_IMPROVEMENT, str)
            
            # Tropo model
            self.my_attributes.tropo_model = read_parameter(parameters, "Tropo model", None, str)
            self.my_attributes.tropo_error_correlation = read_parameter(parameters, "Tropo error correlation", None, int)
                
            if self.my_attributes.tropo_model == 'gaussian':
                self.my_attributes.tropo_error_stdv = read_parameter(parameters, "Tropo error stdv", None, float)
                self.my_attributes.tropo_error_mean = read_parameter(parameters, "Tropo error mean", None, float)
            if self.my_attributes.tropo_model == 'map':
                self.my_attributes.tropo_error_map_file = os.path.expandvars(parameters.getValue("Tropo error map file")  )           
            # Height model parameters
            
            # More complex model
            self.my_attributes.height_model = read_parameter(parameters, "Height model", None, str)
            
            if self.my_attributes.height_model is None:
                my_api.printInfo("Height only given by a simple model A/t0/T")
                # Simple model
                self.my_attributes.height_model_a = read_parameter(parameters, "Constant height model A", HEIGHT_MODEL_A, float)
                self.my_attributes.height_model_t0 = read_parameter(parameters, "Constant height model t0", HEIGHT_MODEL_t0, float)
                self.my_attributes.height_model_period = read_parameter(parameters, "Constant height model period", HEIGHT_MODEL_PERIOD, float)
            
            else:
                
                if self.my_attributes.height_model in ["polynomial", "gaussian"]:
                    # Simple model
                    self.my_attributes.height_model_a = read_parameter(parameters, "Constant height model A", HEIGHT_MODEL_A, float)
                    self.my_attributes.height_model_t0 = read_parameter(parameters, "Constant height model t0", HEIGHT_MODEL_t0, float)
                    self.my_attributes.height_model_period = read_parameter(parameters, "Constant height model period", HEIGHT_MODEL_PERIOD, float)
                    # Specific attributes
                    self.my_attributes.height_model_min_area = read_parameter(parameters, "Height 2d model min area", HEIGHT_MODEL_MIN_AREA, float)
                    # Only for gaussian model
                    if self.my_attributes.height_model == "gaussian":
                        self.my_attributes.height_model_stdv = read_parameter(parameters, "Height 2d model stdv", HEIGHT_MODEL_STDV, float)
                    
                elif self.my_attributes.height_model == "reference_height":
                    self.my_attributes.height_name = read_parameter(parameters, "Height shp attribute name", "HEIGHT", str)
                    
                elif self.my_attributes.height_model == "reference_file": # True height file
                    self.my_attributes.trueheight_file = os.path.expandvars(parameters.getValue("True height file"))
                    if self.my_attributes.trueheight_file is None:
                        my_api.exitWithError("True height file not filled")
                    my_tools.testFile(self.my_attributes.trueheight_file, IN_extent=".nc")
                        
                else:
                    my_api.exitWithError("Height model value UNKNOWN => should be one of polynomial / gaussian / reference_height / reference_file")
                
            
        except IOError:
            my_api.exitWithError("[sisimp_processing/run_preprocessing] Parameter file not found = %s" % IN_paramFile)

        # Check input shapefile
        # Check if all needed files of shapefile are there   
        file_missing = False
        if not os.path.isfile(self.my_attributes.shapefile_path + ".dbf"):
            my_api.printError("The file " + self.my_attributes.shapefile_path + ".dbf is missing.")
            file_missing = True
        if not os.path.isfile(self.my_attributes.shapefile_path + ".shp"):
            my_api.printError("The file " + self.my_attributes.shapefile_path + ".shp is missing.")
            file_missing = True
        if not os.path.isfile(self.my_attributes.shapefile_path + ".shx"):
            my_api.printError("The file " + self.my_attributes.shapefile_path + ".shx is missing.")
            file_missing = True
        if file_missing:
            raise IOError("One or several shapefile files are missing, check logs to know which one")
        # Loading shapefile
        driver = ogr.GetDriverByName(str("ESRI Shapefile"))
        da_shape_file = driver.Open(self.my_attributes.shapefile_path + ".shp", 0)  # 0 means read-only. 1 means writeable.
        wb_layer = da_shape_file.GetLayer()
        # Check if the informations in the shapefile are right   
        shp_srs = wb_layer.GetSpatialRef()
        lonlat_srs = osr.SpatialReference()
        lonlat_srs.ImportFromEPSG(4326)
        if not lonlat_srs.IsSame(shp_srs):
            raise IOError("This is not a shapefile in lon/lat WGS84 projection")
        # Name of height variable if used
        if self.my_attributes.height_model == "reference_height":
            field_names = [field.name for field in wb_layer.schema]
            if self.my_attributes.height_name not in field_names:
                my_api.exitWithError("%s attribute not found in %s; reference_height model option can't be applied" % (self.my_attributes.height_name, self.my_attributes.shapefile_path+".shp"))


        # Load the noise tab
        noise_file_path = ""
        try:
            self.my_attributes.noise_height = np.loadtxt(os.path.expandvars(parameters.getValue("Noise file path")), skiprows=1)
            self.my_attributes.noise_height[:, 1] = self.my_attributes.noise_multiplier_factor * self.my_attributes.noise_height[:, 1]
            self.my_attributes.dw_detected_noise_height = np.loadtxt(os.path.expandvars(parameters.getValue("Noise file path")), skiprows=1)
            self.my_attributes.dw_detected_noise_height[:,1]=self.my_attributes.dw_detected_noise_factor*self.my_attributes.noise_height[:,1]
            
        except IOError:
            my_api.exitWithError("Noise file %s not found in the folder %s " % (noise_file_path, os.path.dirname(__file__)))


        # Load the tile database file
        noise_file_path = ""
        try:
            archive = zipfile.ZipFile(os.path.expandvars(parameters.getValue("Tile database path")), "r")
            imgfile = archive.open("tiles_full.txt")
            self.my_attributes.tile_database = np.loadtxt(imgfile, skiprows=1)
        except IOError:
            my_api.exitWithError("Tile database not found ")
            
            
        # Orbit processing
        # For loop on orbit_file generated
        if self.my_attributes.multi_orbit_option == 'yes' or self.my_attributes.multi_orbit_option == 'passplan':
            try:
                # Retrieve all orbit files
                for path, subdirs, files in (os.walk(run_directory_for_orbits)):
                    for file in files:
                        file_name = os.path.join(path, file)
                        if len(re.findall("cycle_[0-9]+_pass_[0-9]+", file_name)) > 0:
                            orbit_number = int(file_name[-6:-3])
                            orbit_cycle = int(file_name[-16:-13])
                            self.my_attributes.orbit_list.append([orbit_cycle, orbit_number, os.path.join(run_directory_for_orbits, file_name)])
                            
                # When passplan option is selected: self.orbit_files corresponds to the passplan
                if self.my_attributes.multi_orbit_option == 'passplan':
                    try:
                        plan = my_plan.orbitPassplan(self.my_attributes.passplan_path)
                        TMP_new_orbit_list = [None] * len(plan.cycle_orbit_pairs)
                        for orbit_file in self.my_attributes.orbit_list:
                            for ind in plan.plan_orbit[orbit_file[1]]:
                                TMP_new_orbit_list[ind] = [plan.cycle_orbit_pairs[ind][0], orbit_file[1], orbit_file[2]]
                        self.my_attributes.orbit_list = TMP_new_orbit_list
                    except FileNotFoundError:
                        my_api.exitWithError("Passplan file not found: %s" % self.my_attributes.passplan_path)
                
            except IndexError:
                my_api.printError("Orbit file not found")
                my_api.exitWithError("Check orbit files present in the orbit folder")

        else:
            # Build the orbit_path to get the orbit_file   
            try:
                for path, subdirs, files in (os.walk(run_directory_for_orbits)):
                    for file in files:
                        file_name = os.path.join(path, file)
                        if ("pass_%04d" % self.my_attributes.orbit_number) in file_name:
                            self.my_attributes.orbit_list.append([self.my_attributes.cycle_number, self.my_attributes.orbit_number, file_name])

            except IndexError:
                my_api.printError("Orbit file not found = %s" % path_to_orbit_file)
                my_api.exitWithError("Please check that orbit number %d has been generated" % self.my_attributes.orbit_number)
        my_api.printInfo("")
        my_api.printInfo("List of orbit files to process =")
        
        for elem in self.my_attributes.orbit_list:
            my_api.printInfo("Cycle=%03d - Pass=%03d - Orbit file=%s" % (elem[0], elem[1], os.path.basename(elem[2])))

        # Create shapefile
        try:
            self.my_attributes.create_shapefile = parameters.getValue("Create shapefile").lower()
            self.my_attributes.create_shapefile = self.my_attributes.create_shapefile in ['oui', 'yes', 'yep']
        except Exception:
            self.my_attributes.create_shapefile = False
            my_api.printInfo("No Create shapefile parameter set, no shapefile will be created")

        # Create dummy L2_HR_PIXCVecRiver product, associated to pixel cloud  
        try:
            self.my_attributes.create_pixc_vec_river = parameters.getValue("Create dummy pixc vec river file").lower()
            self.my_attributes.create_pixc_vec_river = self.my_attributes.create_pixc_vec_river in ['oui', 'yes', 'yep']
        except Exception:
            self.my_attributes.create_pixc_vec_river = False
            my_api.printInfo("No Create dummy pixc vec river file parameter set, no L2_HR_PIXCVecRiver file will be created")

        # self.my_attributes.compute_pixc_vec_river to True only if self.my_attributes.create_pixc_vec_river is True and RIV_FLAG field is here
        self.my_attributes.compute_pixc_vec_river = False
        if self.my_attributes.create_pixc_vec_river and (wb_layer.FindFieldIndex(str("RIV_FLAG"), True) != -1):
            self.my_attributes.compute_pixc_vec_river = True


    def run_processing(self):
        """Main process, computations are done here"""
        my_api.printInfo("")
        my_api.printInfo("")
        my_api.printInfo("[sisimp_processing] PROCESSING...")
        my_api.printInfo("")

        for elem in self.my_attributes.orbit_list:  # Process per element in orbit list = triplet (cycle_number, orbit_number, orbit_file)
            
            my_api.printInfo(">>> CYCLE %03d and ORBIT %03d <<<" % (elem[0], elem[1]))
            my_api.printInfo("")
            
            # 1 - Read orbit file
            self.my_attributes = sisimp_fct.read_orbit(elem[2], elem[0], self.my_attributes)
            my_api.printInfo("")
            
            # 2 - Init SISIMP filenames object
            self.my_attributes.sisimp_filenames = my_names.sisimpFilenames(self.my_attributes.out_dir, self.my_attributes.mission_start_time, self.my_attributes.cycle_duration, elem[0], elem[1])
            
            # 3 - Process right swath
            self.my_attributes = sisimp_fct.make_pixel_cloud("Right", elem[0], elem[1],self.my_attributes)
            my_api.printInfo("")
            
            # 4 - Process left swath
            self.my_attributes = sisimp_fct.make_pixel_cloud("Left", elem[0], elem[1], self.my_attributes)
            my_api.printInfo("")
            
            # 5 - Write swath polygons shapefile
            sisimp_fct.write_swath_polygons(self.my_attributes)
            my_api.printInfo("")
            my_api.printInfo("")

    def run_postprocessing(self):
        """
        Run post-processing
        """
        my_api.printInfo("")
        my_api.printInfo("")
        my_api.printInfo("[sisimp_processing] POST-PROCESSING...")
        my_api.printInfo("Nothing to do...")
        my_api.printInfo("")
