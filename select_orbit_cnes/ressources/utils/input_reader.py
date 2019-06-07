"""
.. module input_reader.py
    :synopsis: handle RDF parameter file for select_orbit_cnes processing
    Created on 21 sept. 2015

.. moduleauthor: Capgemini
    
This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import os

from ressources.const_input_reader import LON_WEST, LON_EAST, LAT_NORTH, LAT_SOUTH, AZIMUTH_SPACING, SWATH, NEAR_RANGE, GDEM_PREFIX, ORBIT_REPOSITORY, MISSION_NAME, MISSION_START_TIME, MAKE_PASS_PLAN, SIMULATION_START, SIMULATION_STOP
from ressources.rdf.rdf_reader import RdfReader
from ressources.rdf.rdf_enums import RDF_DEFAULT


class Input_Reader(object):
    """
    Handle RDF parameter file for select_orbit_cnes processing
    Getters are provided for convenience access.
    """

    def __init__(self, parameter_file):
        """
        Constructor.
        Initialize internal attributes

        Args:
            parameter filename(str): the full path to the parameter file

        Raises:
            FileNotFoundError
        """
        
        # Init self attributes
        self.mission_name = None
        self.mission_start_time = None
        self.longitude_west = None
        self.longitude_east = None
        self.latitude_north = None
        self.latitude_south = None
        self.azimuth_spacing = None
        self.swath = None
        self.near_range = None
        self.dem_file = None
        self.gdem_orbit_prefix = None
        self.simulation_start = None
        self.simulation_stop = None

        # Save parameter file path
        if os.path.exists(parameter_file):
            self.parameter_file = parameter_file
        else:
            raise FileNotFoundError("Parameter file doesn't exist: %s" % parameter_file)
    
    #----------------------------------

    def read_param_file(self):
        """
        Read the parameter file
        """
        
        # 1 - Pointer to parameter file (RDF)
        parameter_file = RdfReader(self.parameter_file)

        # 2 - Getting parameters
        
        # 2.1 - Mission specific parameters
        self.mission_name = parameter_file.get_parameter(RDF_DEFAULT, MISSION_NAME)
        self.mission_start_time = parameter_file.get_parameter_or_default(RDF_DEFAULT, MISSION_START_TIME, None)
        
        # 2.2 - Orbit parameters
        self.orbit_repository = os.path.expandvars(parameter_file.get_parameter_or_default(RDF_DEFAULT, ORBIT_REPOSITORY, None))

        self.latitude_south = parameter_file.get_parameter(RDF_DEFAULT, LAT_SOUTH)
        self.latitude_north = parameter_file.get_parameter(RDF_DEFAULT, LAT_NORTH)
        self.longitude_west = parameter_file.get_parameter(RDF_DEFAULT, LON_WEST)
        self.longitude_east = parameter_file.get_parameter(RDF_DEFAULT, LON_EAST)

        self.azimuth_spacing = parameter_file.get_parameter_or_default(RDF_DEFAULT, AZIMUTH_SPACING, None)
        self.swath = parameter_file.get_parameter_or_default(RDF_DEFAULT, SWATH, None)
        self.near_range = parameter_file.get_parameter_or_default(RDF_DEFAULT, NEAR_RANGE, None)
        
        # 2.3 - Pass plan parameters
        self.make_pass_plan = parameter_file.get_parameter_or_default(RDF_DEFAULT, MAKE_PASS_PLAN, None)
        self.simulation_start = parameter_file.get_parameter_or_default(RDF_DEFAULT, SIMULATION_START, None)
        self.simulation_stop = parameter_file.get_parameter_or_default(RDF_DEFAULT, SIMULATION_STOP, None)
        
        # 2.4 - Output parameters
        self.gdem_orbit_prefix = parameter_file.get_parameter_or_default(RDF_DEFAULT, GDEM_PREFIX, None)
    
    #----------------------------------
    # Getters
    #----------------------------------

    def get_mission_name(self):
        return self.mission_name

    def get_mission_start_time(self):
        return self.mission_start_time
        
    def get_orbit_repository(self):
        return self.orbit_repository

    def get_south_latitude(self):
        return self.latitude_south

    def get_north_latitude(self):
        return self.latitude_north

    def get_west_longitude(self):
        return self.longitude_west

    def get_east_longitude(self):
        return self.longitude_east

    def get_azimuth_spacing(self):
        return self.azimuth_spacing

    def get_swath(self):
        return self.swath

    def get_near_range(self):
        return self.near_range

    def get_make_pass_plan(self):
        return self.make_pass_plan

    def get_simulation_start(self):
        return self.simulation_start

    def get_simulation_stop(self):
        return self.simulation_stop

    def get_gdem_prefix(self):
        return self.gdem_orbit_prefix
