#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# $Id:
#
# ======================================================
#
# Project : SWOT
# Program : common python modules
# Produit par Capgemini.
#
# ======================================================
# HISTORIQUE
#
# FIN-HISTORIQUE
# ======================================================
'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


"""
.. module:: service_error.py
    :synopsis: basic error class for Swot

"""

####################################################################
# Base class SwotError
####################################################################
class SwotError(Exception):
    """
        Base class for exceptions in SWOT
    """
    pass
####################################################################
# Main class definition
####################################################################
####################################################################
# Error class definition configFileError inherits the SwotError class
####################################################################
class ConfigFileError(SwotError):
    """
        Base subclass for exception in the configuration file
        :param msg: explanation of the error
    """
    def __init__(self, msg):
        SwotError.__init__(self, msg)
        self.msg = msg

    def __str__(self):
        return repr(self.msg)

####################################################################
# Error class definition processingError inherits the SwotError class
####################################################################
class ProcessingError(SwotError):
    """
        Base subclass for exception in the main processing
        :param msg: explanation of the error
    """
    def __init__(self, msg, logger):
        SwotError.__init__(self, msg)
        self.msg = msg
        logger.error(self.msg)

    def __str__(self):
        return repr(self.msg)

####################################################################
# Sub class definition
####################################################################
####################################################################
# List of error class definition for processing
# inherits the ProcessingError class
####################################################################

####################################################################
# SAS lake_tile error
####################################################################
class SASLakeTileError(ProcessingError):
    """
        Exception raised for errors during execution of
        the SAS Lake tile
        :param msg: explanation of the error
        :param logger: instance to the logger
    """
    def __init__(self, msg, logger):
        self.msg = "Error in SAS Lake tile " + msg
        ProcessingError.__init__(self, self.msg, logger)
    def __str__(self):
        return self.msg

####################################################################
# SAS lake_sp error
####################################################################
class SASLakeSpError(ProcessingError):
    """
        Exception raised for errors during execution of
        the SAS Lake sp
        :param msg: explanation of the error
        :param logger: instance to the logger
    """
    def __init__(self, msg, logger):
        self.msg = "Error in SAS Lake sp " + msg
        ProcessingError.__init__(self, self.msg, logger)
    def __str__(self):
        return self.msg

####################################################################
# SAS lake_avg error
####################################################################
class SASLakeAvgError(ProcessingError):
    """
        Exception raised for errors during execution of
        the SAS Lake avg
        :param msg: explanation of the error
        :param logger: instance to the logger
    """
    def __init__(self, msg, logger):
        self.msg = "Error in SAS Lake avg " + msg
        ProcessingError.__init__(self, self.msg, logger)
    def __str__(self):
        return self.msg

####################################################################
# List of error class definition for the configuration file,
# inherits the configFileError class
####################################################################
class ParameterError(ConfigFileError):
    """
        Exception raised for errors in a parameter in the configuration file
        (like missing mandatory variable)
        :param section: name of the section
        :param msg: explanation of the error
    """
    def __init__(self, section, msg):
        self.section = section
        self.msg = "Error: In section " + repr(self.section) + ", " + repr(msg)
        ConfigFileError.__init__(self, self.msg)
    def __str__(self):
        return self.msg

class DirFileError(ConfigFileError):
    """
        Exception raised for errors in mandatory directory
        :param directory: name of the directory
    """
    def __init__(self, directory):
        self.directory = directory
        self.msg = "Error: " + repr(self.directory) + " doesn't exist"
        ConfigFileError.__init__(self, self.msg)
    def __str__(self):
        return self.msg

class ConfigError(ConfigFileError):
    """
        Exception raised for configuration errors in the configuration file
        (like incompatible parameters, bad value)
        :param msg: explanation of the error
    """
    def __init__(self, msg):
        self.msg = "Error: " + repr(msg)
        ConfigFileError.__init__(self, self.msg)
    def __str__(self):
        return self.msg

class FileError(ConfigFileError):
    """ Exception raised for errors inside an input file
        (like a bad format or missing variable)
        :param msg: explanation of the error
    """
    def __init__(self, msg):
        self.msg = "Error: " + repr(msg)
        ConfigFileError.__init__(self, self.msg)
    def __str__(self):
        return self.msg

####################################################################
####################################################################


