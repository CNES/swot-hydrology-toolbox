"""
.. module input_reader.py
    :synopsis: SAR Geo. equation module
    Created on 21 sept. 2015

.. moduleauthor: Capgemini

    $Id: input_reader.py 1465 2016-07-01 10:05:12Z nestival $
This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""

from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division


import os

from ressources.const_input_reader import LON_WEST, LON_EAST, LAT_NORTH, LAT_SOUTH, AZIMUTH_SPACING, SWATH, NEAR_RANGE, GDEM_PREFIX, DEM_FILE, ORBIT_REPOSITORY, START_MISSION_TIME, MAKE_PASS_PLAN, SIMULATION_START, SIMULATION_STOP, CYCLE_DURATION

from ressources.rdf.rdf_reader import RdfReader
from ressources.rdf.rdf_enums import RDF_DEFAULT

class Input_Reader(object):
    """
    This class permit to read all config files given by SAM to prepare inputs for this module

    This class reads the 1 file given during the process :
        - parameter file (.par)

    Getters are provided for convenience accesses.
    """

    def __init__(self, parameter_file):
        """
        Constructor.
        Initialize internal attributes

        Args:
            parameter filename(str): the full path to the parameter file

            api(object): instance to the api top level for enabling logging messages in this module

        Raises:
            ModuleException
        """
        self.radar_file = None
        self.longitude_west = None
        self.longitude_east = None
        self.latitude_north = None
        self.latitude_south = None
        self.azimuth_spacing = None
        self.swath = None
        self.near_range = None
        self.dem_file = None
        self.gdem_orbit_prefix = None
        self.start_mission_time = None
        self.simulation_start = None
        self.simulation_stop = None
        self.cycle_duration = None

        # enable message logging in this class

        # save files path
        if os.path.exists(parameter_file):
            self.parameter_file = parameter_file
        else:
            raise FileNotFoundError("Parameter file doesn't exist: %s" % parameter_file)


    def read_param_file(self):
        """
        This function reads the parameter file

        Raises:
            ModuleException
        """
        
        # Reading parameter file (RDF)
        parameter_file = RdfReader(self.parameter_file)

        # Getting parameters

        self.longitude_west = parameter_file.get_parameter(RDF_DEFAULT, LON_WEST)
        self.longitude_east = parameter_file.get_parameter(RDF_DEFAULT, LON_EAST)
        self.latitude_north = parameter_file.get_parameter(RDF_DEFAULT, LAT_NORTH)
        self.latitude_south = parameter_file.get_parameter(RDF_DEFAULT, LAT_SOUTH)

        self.azimuth_spacing = parameter_file.get_parameter_or_default(RDF_DEFAULT, AZIMUTH_SPACING, None)
        self.swath = parameter_file.get_parameter_or_default(RDF_DEFAULT, SWATH, None)
        self.near_range = parameter_file.get_parameter_or_default(RDF_DEFAULT, NEAR_RANGE, None)
        self.dem_file = parameter_file.get_parameter_or_default(RDF_DEFAULT, DEM_FILE, None)
        self.gdem_orbit_prefix = parameter_file.get_parameter_or_default(RDF_DEFAULT, GDEM_PREFIX, None)

        self.orbit_repository = parameter_file.get_parameter_or_default(RDF_DEFAULT, ORBIT_REPOSITORY, None)
        self.start_mission_time = parameter_file.get_parameter_or_default(RDF_DEFAULT, START_MISSION_TIME, None)
        self.make_pass_plan = parameter_file.get_parameter_or_default(RDF_DEFAULT, MAKE_PASS_PLAN, None)
        self.simulation_start = parameter_file.get_parameter_or_default(RDF_DEFAULT, SIMULATION_START, None)
        self.simulation_stop = parameter_file.get_parameter_or_default(RDF_DEFAULT, SIMULATION_STOP, None)
        self.cycle_duration = parameter_file.get_parameter_or_default(RDF_DEFAULT, CYCLE_DURATION, None)

    def get_north_latitude(self):
        return self.latitude_north

    def get_south_latitude(self):
        return self.latitude_south

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

    def get_gdem_prefix(self):
        return self.gdem_orbit_prefix
        
    def get_orbit_repository(self):
        return self.orbit_repository

    def get_start_mission(self):
        return self.start_mission_time

    def get_make_pass_plan(self):
        return self.make_pass_plan

    def get_simulation_start(self):
        return self.simulation_start

    def get_simulation_stop(self):
        return self.simulation_stop

    def get_cycle_duration(self):
        return self.cycle_duration
