# -*- coding: utf8 -*-
'''
.. module my_api.py
    :synopsis: Use appropriate functions as the run environment is stand-alone or SAM  
    Created on 02/01/2018

.. module author: Claire POTTIER - CNES DSO/SI/TR

Copyright (c) 2018 CNES. All rights reserved.
'''
from __future__ import absolute_import, division, print_function, unicode_literals   
  
import logging
import sys

# Global variables
GEN_PRINT_LEVEL = "INFO" # Level of print info; values = DEBUG - INFO
GEN_ENV = 2 # Run environment: 0=SAM - 1=stand-alone (with logging) - 2(default)=stand-alone (with print)
GEN_API = None # Variable of SAM api if used (ie GEN_ENV=0); None=stand-alone (default)


#----------------------------------------
# Initialization of environment
#----------------------------------------

def initApi(IN_api, IN_verbose_level):
    '''
    Initialize api with SAM api and verbose level
    
    :param IN_api: SAM api if use of SAM environment
    :type IN_api: ModuleCommon
    :param IN_verbose_level: verbose level; level = DEBUG or INFO (default)
    :type IN_verbose_level: string
    '''
    
    # Set SAM flag
    global GEN_ENV    
    GEN_ENV = 0
    
    # Set SAM api
    global GEN_API
    GEN_API = IN_api
    
    # Set verbose level
    global GEN_PRINT_LEVEL
    if ( IN_verbose_level == "DEBUG" ):
        GEN_PRINT_LEVEL = "DEBUG"
    else:
        GEN_PRINT_LEVEL = "INFO"
    
def initLogger(IN_logFile, IN_verbose_level):
    '''
    Initialize logger and verbose level
    
    :param IN_log_file: full path of the log file
    :type IN_log_file: string
    :param IN_verbose_level: verbose level; level = DEBUG or INFO (default)
    :type IN_verbose_level: string
    '''
    
    # Unset SAM flag
    global GEN_ENV    
    GEN_ENV = 1
    
    # Set SAM api
    global GEN_API
    GEN_API = None
    
    # Set verbose level
    global GEN_PRINT_LEVEL
    if ( IN_verbose_level == "DEBUG" ):
        GEN_PRINT_LEVEL = "DEBUG"
        vLevel = logging.DEBUG
    else:
        GEN_PRINT_LEVEL = "INFO"
        vLevel = logging.INFO
    
    # Configure logger
    logging.basicConfig(filename=IN_logFile, format='%(asctime)s [%(levelname)s] %(message)s', level=vLevel)
    
def setVerbose(IN_verbose_level):
    '''
    Setter verbose level = GEN_PRINT_LEVEL
    
    :param IN_verbose_level: verbose level; level = DEBUG or INFO (default)
    :type IN_verbose_level: string
    
    :return GEN_PRINT_LEVEL: the true verbose level
    :rtype: string
    '''
    
    # Set verbose level
    global GEN_PRINT_LEVEL
    if ( IN_verbose_level == "DEBUG" ):
        GEN_PRINT_LEVEL = "DEBUG"
    else:
        GEN_PRINT_LEVEL = "INFO"
        
    return GEN_PRINT_LEVEL
            
        
#----------------------------------------
# Printing function
#----------------------------------------

def printInfo(IN_msg):
    '''
    Print information message
    
    :param IN_msg: message to print
    :type IN_msg: string
    '''
    
    if ( GEN_ENV == 0 ):
        GEN_API.info(IN_msg)
    elif ( GEN_ENV == 1 ):
        logging.info(IN_msg)
    else:
        print("[I] %s" % ( IN_msg ))

def printDebug(IN_msg):
    '''
    Print debug message
    
    :param IN_msg: message to print
    :type IN_msg: string
    '''
    
    if ( GEN_PRINT_LEVEL == "DEBUG" ):
        if ( GEN_ENV == 0 ):
            GEN_API.info(IN_msg)
        elif ( GEN_ENV == 1 ):
            logging.debug(IN_msg)
        else:
            print("[D] %s" % ( IN_msg )) 

def printError(IN_msg):
    '''
    Print error message
    
    :param IN_msg: message to print
    :type IN_msg: string
    '''
    
    if ( GEN_ENV == 0 ):
        GEN_API.error(IN_msg)
    elif ( GEN_ENV == 1 ):
        logging.error(IN_msg)
    else:
        print("[E] %s" % ( IN_msg ))

        
#----------------------------------------
# Exit function
#----------------------------------------

def exitWithError(IN_msg):
    '''
    Exit program with an error message
    
    :param IN_msg: message to print
    :type IN_msg: string
    '''
    
    printError(IN_msg)
    
    if ( GEN_ENV == 0 ):
        GEN_API.end_module(False, 101)
    else:
        sys.exit(IN_msg)