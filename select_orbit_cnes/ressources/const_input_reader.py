"""
.. module const_input_reader.py
    :synopsis: Constants for input reader part
    Created on 21 sept. 2015

.. moduleauthor: Capgemini

    $Id: const_input_reader.py 1465 2016-07-01 10:05:12Z nestival $
    
This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division


# Error messages

# Log messages
TRACE_START_READING_PARAM = "Starting reading global parameter file."
TRACE_END_READING_PARAM = "The global parameter file was successfully read."

# RDF parameters name
MISSION_NAME = "Mission name"
MISSION_START_TIME = "Mission start time"
ORBIT_REPOSITORY = "Orbit repository"
LAT_SOUTH = "DEM south latitude"
LAT_NORTH = "DEM north latitude"
LON_WEST = "DEM west longitude"
LON_EAST = "DEM east longitude"
AZIMUTH_SPACING = "Azimuth spacing"
SWATH = "Swath width"
NEAR_RANGE = "NR cross track"
MAKE_PASS_PLAN = "passplan"
SIMULATION_START = "simulation_start_time"
SIMULATION_STOP = "simulation_stop_time"
GDEM_PREFIX = "GDEM Orbit prefix"
