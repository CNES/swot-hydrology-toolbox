# -*- coding: utf8 -*-
"""
.. module:: my_rdf_file.py
    :synopsis: Deals with RDF parameter files (reader)
    Created on 30/01/2018

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

Copyright (c) 2018 CNES. All rights reserved.
"""
from __future__ import absolute_import, division, print_function, unicode_literals 

import cnes.common.lib.my_api as my_api 


class myRdfReader(object):

    def __init__(self, IN_filename):
        """
        Constructor
        
        :param IN_filename: full path of the RDF file
        :type IN_filename: string
        
        Variables:
            filename / string: full path of the RDF file
            parameters / dict: parser of the RDF file
        """
        
        # 1 - Set filename attribute
        self.filename = IN_filename
        
        # 2 - Init parameters list as a dictionary
        self.parameters = {}
        
        # 3 - Read file and fill dictionary
        self.setParams()
    
    # ----------------------------------------
    
    def setParams(self):
        """
        Read RDF file and fill dictionary with key=value
        """
        
        # 1 - Read RDF file
        rdf_reader = open(self.filename, 'r')  # Open file in reading mode
        rdf_lines = rdf_reader.read().splitlines()  # Read file and store line / line
        
        # 2 - Assign each line as a key of dictionary
        for line in rdf_lines:
            
            # Consider line only if contain "=" (not empty, not a comment, ...)
            if ("=" in line) and (not line.startswith("!")):
                
                # 2.1 - Remove comment if exists
                TMP_line = line.split("!")[0]
                
                # 2.2 - Assign key and value
                TMP_split = TMP_line.split("=")
                self.parameters[TMP_split[0].strip()] = TMP_split[1].strip()  # Remove blank at begin and end of the string
                
        # 3 - Close file
        rdf_reader.close()
    
    # ----------------------------------------
        
    def getValue(self, IN_key):
        """
        Get the value associated to the key named IN_key
        
        :param IN_key: name of the key
        :type IN_key: string
        
        :return: value associated to the key, None if doesn't exist
        :rtype: string
        """
        if IN_key in self.parameters:
            return self.parameters.get(IN_key)
        else:
            for key in self.parameters.keys():
                if IN_key in key:
                    return self.parameters.get(key)
            # Key doesn't exist
            my_api.printInfo("[my_rdf_file/getValue] Unkown key = %s" % IN_key)
            return None

    def getFloatOrNone(self, key):
        """
        Get the value associated to the key named IN_key and convert it as a float
        
        :param key: name of the key
        :type key: string
        
        :return: value associated to the key, None if doesn't exist
        :rtype: float
        """
        if self.getValue(key) == "None":
            return None
        else:
            return float(self.getValue(key))


#######################################


if __name__ == "__main__":
    
    paramFile = "D:\\Utilisateurs\\pottierc\\Documents\\workspace_qgis\\sisimp_tests_unitaires\\FT2_Cas2\\3_lake_tile_pc\\parameter_lake_tile.rdf"
    parameters = myRdfReader(paramFile)
        
    # Get working directories
    print(str(parameters.getValue("PIXC_MAIN directory")))
    print(str(parameters.getValue("PIXC_SENSOR directory")))
    print(str(parameters.getValue("PIXC_VEC_RIVER directory")))
    print(str(parameters.getValue("Output directory")))
            
    # Get tile(s) information
    # Cycle number
    TMP_param = parameters.getValue("Cycle")
    if TMP_param:
        print(int(TMP_param))
    # Orbit number
    TMP_param = parameters.getValue("Orbit")
    if TMP_param:
        print(int(TMP_param))
    # Tile id
    TMP_param = parameters.getValue("Tile")
    if TMP_param:
        print(str(TMP_param))
            
    # Get config parameters
    print(str(parameters.getValue("Lake a priori database")))
    print(str(parameters.getValue("Classif flags to keep")))
    print(float(parameters.getValue("Min size for lake")))
    print(int(parameters.getValue("Improve geolocation")))
    
    # Get inexistant value
    print(parameters.getValue("I don't exist"))
