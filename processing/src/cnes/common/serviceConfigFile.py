#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# $Id:
#
# ======================================================
#
# Project : SWOT
# Program : common python module
# Produit par Capgemini.
#
# ======================================================
# HISTORIQUE
#
# FIN-HISTORIQUE
# ======================================================
"""
   module:: serviceConfigFile.py
    :synopsis: Manage the configuration file for SWOT
"""

import sys
import configparser
from configparser import RawConfigParser
import cnes.common.serviceError as serviceError

# this is a pointer to the module object instance itself.
THIS = sys.modules[__name__]

# declaration of path_conf and cfg variables
THIS.path_conf = None
THIS.cfg = None

# method to clear the config
def clear_config():
    """
        This function clear the configuration previously loadded
    """
    if THIS.path_conf is not None:
        # also in local function scope. no scope specifier like global is needed
        THIS.path_conf = None
        THIS.cfg = None

# method to get instance of ServiceConfigFile
def get_instance():
    """
        This function return the instance of service_config_file
    """
    return THIS.cfg

class ServiceConfigFile(RawConfigParser):
    """
        The class ServiceConfigFile defines all methods to access to the
        configuration file and to check the variables.
    """

    # Constructor
    def __init__(self, path_conf):
        """
            Init class ServiceConfigFile
            :param path_conf: string path of the config file
        """
        if THIS.path_conf is None:
            # first set of path_conf and cfg
            if path_conf is None:
                raise Exception("First call to ServiceConfigFile: path_conf_name is not define")
            THIS.path_conf = path_conf
            # we call the constructor of mother class
            RawConfigParser.__init__(self)
            # we load the configuration file
            self.read(path_conf)
            # we save instance of class
            THIS.cfg = self
        #initialize_config(self, path_conf)
        self.path_conf = THIS.path_conf

    # Default message
    def __repr__(self):
        """
            Default message print
        """
        return "Configuration file : " + self.path_conf

    def test_var_config_file(self, section, variable, var_type, valeurs="",
                             val_defaut=None, valid_min=None, valid_max=None):
        """
            This function check if variable is in obj
            and if it has var_type type.
            Optionnaly it can check if variable has values in valeurs
            Exit the code if any error are detected
            :param section: section name of the obj where to find
            :param variable: string name of the variable
            :param var_type: type type of the variable for verification
            :param valeurs: string list of the possible value of variable
            :param val_defaut: value by default if variable is not in the configuration file
            :param valid_min: minimum value possible for this variable
            :param valid_max: maximum value possible for this variable
        """
        try:
            # get variable in function of type
            tmp_var = self.get_var_from_type(section, variable, var_type)

        except configparser.NoSectionError:
        # error section not exist
            raise serviceError.ConfigFileError("Section '" + str(section)
                                               + "' is not in the configuration file")
        except configparser.NoOptionError:
            if val_defaut is not None:
                # Variable not exist, but a default value has been set
                # create variable with this value
                self.set(section, variable, str(val_defaut))
            else:
                # Error variable not exist
                message = "mandatory variable '" + str(variable) +\
                "' is missing in the configuration file"
                raise serviceError.ParameterError(section, message)
        except ValueError:
            # wrong type for variable
            message = "variable '" + str(variable) +\
            "' has a wrong type expected: " + str(var_type)
            raise serviceError.ParameterError(section, message)

        # variable exist get value
        # if a list of value is given, tests it
        if valeurs != "":
            val_ok = 0
            for index in range(len(valeurs)):
                if tmp_var == valeurs[index]:
                    val_ok = 1
            if val_ok == 0:
                message = "bad value for '" + variable +\
                "' variable. Value accepted: " + str(valeurs) +\
                " Value read: " + str(tmp_var)
                raise serviceError.ParameterError(section, message)
        # if valid_min is given, check if variable is not lower than this value
        if (valid_min is not None) and (tmp_var < valid_min):
            message = "Value too small for '" + variable +\
            "' variable. Limits: " + str(valid_min) +\
            " to " + str(valid_max) + " Value read: " + str(tmp_var)
            raise serviceError.ParameterError(section, message)
        # if valid_max is given, check if variable is not greater than this value
        if (valid_max is not None) and (tmp_var > valid_max):
            message = "Value too big for '" + variable +\
            "' variable. Limits: " + str(valid_min) +\
            " to " + str(valid_max) + " Value read: " + str(tmp_var)
            raise serviceError.ParameterError(section, message)

    def get_var_from_type(self, section, variable, var_type):
        """
            This function return variable reading in function of type requested
            :param section: section name of the obj where to find
            :param variable: string name of the variable
            :param var_type: type type of the variable for verification
        """
        if var_type == int:
            tmp_var = self.getint(section, variable)
        elif var_type == float:
            tmp_var = self.getfloat(section, variable)
        elif var_type == bool:
            tmp_var = self.getboolean(section, variable)
        else:
            tmp_var = self.get(section, variable)

        return tmp_var

