"""
.. module swot_hr.common.rdf.rdf_reader.py
    :synopsis:Class responsible for reading RDF files from JPL simulator
    Created on 24 avr. 2013

.. moduleauthor: Capgemini

    $Id: rdf_reader.py 1465 2016-07-01 10:05:12Z nestival $
This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""

from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import os

from ressources.rdf.rdf_exception import RdfException
from ressources.rdf.rdf_enums import RDF_WRONG_PATH, \
    RDF_DEFAULT, RDF_EOL, RDF_EQUAL, RDF_UNIT, RDF_COMMENT, RDF_WRONG_SECTION,\
    RDF_WRONG_PARAM, RDF_NULL_PARAM


import re

class RdfReader(object):
    """
    Class responsible for reading RDF files from JPL simulator.
    RDF files from the simulator defines the following grammar
    
        SECTION NAME
        param name        (unit) = value;    ! comment
    
    Comment and unit are optional.
    This reader provides functions to assess the following goal: 
    get a parameter value from a section.
    
    If no section name is given, it is saved in RDF_DEFAULT enumerations 
    in rdf_enums.py. 
    """


    def __init__(self, filename):
        """
        Constructor, reads automatically the file and 
        fill a dictionnary with the sections and parameters
        for further direct accesses.
        
        Args:
            filename(str): the full path to the file to read
            
        Raises:
            RdfException
        """
        if os.path.exists(filename):
            # Open the file
            self._filename = filename
            
            # Initialize the dictionary for sections and params
            self._content = dict()
            self._readfile()
        else:
            exception = RdfException(RDF_WRONG_PATH)
            exception.add_data("Path", filename)
            raise exception
        
    def _readfile(self):
        """
        This function reads the file automatically and fills.
        The parsing of the file is designed as follows:
            
            * If a line is not empty
            * If the line contains an equal symbol then it is a param, we get its value
            * Else it is a section name
            
        Raises:
            RdfException
        """
        with open(self._filename) as rdf:
            lines = rdf.readlines()
            
        found_section = RDF_DEFAULT
        for line in lines:
            if len(line.strip()) > 0:
                equal = line.find(RDF_EQUAL)
                # There is info
                if equal == -1:
                    # It has to be a section
                    if line.isupper():
                        # Ok, this is our section
                        found_section = line.strip()
                    else:
                        # It is very bad but the JPL's files contain this kind of lines...
                        pass
                else:
                    # It should be a param line
                    # Break down the line : param name (unit) = value    ! comment
                    # param name
                    parenthesis = line.find(RDF_UNIT, 0, equal) 
                    if parenthesis != -1:
                        # there is a unit, param name is the concatenation of all other words with spaces
                        param = line[0:parenthesis].strip()
                    else:
                        param = line[0:equal].strip()
                        
                    # value, ignore comment
                    exclamation = line.find(RDF_COMMENT, equal)
                    if exclamation != -1:
                        value = line[equal + 1:exclamation].strip()
                    else:
                        value = line[equal + 1:].strip()
                    
                    # Remove eventual end of line
                    if len(value) > 0: 
                        if value[-1] == RDF_EOL:
                            value = value[:-1]
                    else:
                        value = None
                        
                    # now save the param and its value
                    if not (found_section in self._content):
                        self._content[found_section] = dict()
                    
                    self._content[found_section][param] = value
        
    def get_sections(self):
        """
        Get the list of sections for a RDF file.
        If the RDF_DEFAULT section is used, it is also returned.
        
        Returns:
            list(str). The list of sections
        """        
        return self._content.keys()
    
    
    def get_parameters(self, section):
        """
        Returns the list of parameters in a section.
        
        Args:
            section(str): The section name
        
        Returns:
            list(str) The list of parameters
            
        Raises:
            RdfException
        """
        if section in self._content:
            return self._content[section].keys()
        else:
            exception = RdfException(RDF_WRONG_SECTION)
            exception.add_data("Section", section)
            raise exception
    
    def get_parameter(self, section, parameter):
        """
        Get the parameter value of a named parameter and in
        a specific section.
        
        Args:
            section(str): The section name
            
            parameter(str): The parameter name
            
        Returns:
            str. The parameter value, always as string
        
        Raises:
            RdfException
        """
        if section in self._content:
            if parameter in self._content[section]:
                return self._content[section][parameter]
            else:
                exception = RdfException(RDF_WRONG_PARAM)
                exception.add_data("Section", section)
                exception.add_data("Parameter", parameter)
                raise exception
        else:
            exception = RdfException(RDF_WRONG_SECTION)
            exception.add_data("Section", section)
            raise exception
    
    def get_not_null_parameter(self, section, parameter):
        """
        Get the parameter value of a named parameter and in
        a specific section if the value of the parameter is
        not set raise a RdfException.
        
        Args:
            section(str): The section name
            
            parameter(str): The parameter name
            
        Returns:
            str. The parameter value, always as string
        
        Raises:
            RdfException
        """
        if section in self._content:
            if parameter in self._content[section]:
                value = self._content[section][parameter]
                if value == None:
                    exception = RdfException(RDF_NULL_PARAM)
                    exception.add_data("Section", section)
                    exception.add_data("Parameter", parameter)
                    raise exception
                else:
                    return value
            else:
                exception = RdfException(RDF_WRONG_PARAM)
                exception.add_data("Section", section)
                exception.add_data("Parameter", parameter)
                raise exception
        else:
            exception = RdfException(RDF_WRONG_SECTION)
            exception.add_data("Section", section)
            raise exception
    
    def get_parameter_or_default(self, section, parameter, defaultValue):
        """
        Get the parameter value of a named parameter and in
        a specific section if the value is not found return the defaultValue
        passed in arguments.
        
        Args:
            section(str): The section name
            
            parameter(str): The parameter name
            
            defaultValue(str): The default value to return if parameter not set
            
        Returns:
            str. The parameter value, always as string

        """
        if section in self._content:
            if parameter in self._content[section]:
                value = self._content[section][parameter]
                if value == None:
                    return defaultValue
                else:
                    return value
            else:
                return defaultValue
        else:
            return defaultValue
        
        
    def get_parameter_at_index(self, section, index, parameter_regex):
        """
        Get the parameter value of a named parameter and in
        a specific section.
        
        Args:
            section(str): The section name
            
            index(str): 
            
            parameter_regex(str): Regular expression defining the pattern of the parameter name
            
        Returns:
            str. The parameter value, always as string
        
        Raises:
            RdfException
        """
        if section in self._content:
            
            reg = re.compile(parameter_regex)
            
            count = 0
            for param in sorted(self._content[section].keys()):
                found = reg.match(param) 
                if found != None:
                    if count == index:
                        return found.group(1), self._content[section][param]
                    else:
                        count += 1
            else:
                exception = RdfException(RDF_WRONG_PARAM)
                exception.add_data("Section", section)
                exception.add_data("Parameter index", index)
                exception.add_data("Parameter regex", parameter_regex)
                raise exception
        else:
            exception = RdfException(RDF_WRONG_SECTION)
            exception.add_data("Section", section)
            raise exception
        
        
