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
.. module:: serviceLogger.py
    :synopsis: logging class for Swot

"""

import sys
import logging
import logging.handlers
import cnes.common.serviceConfigFile as serviceConfigFile

# pointer to the module object instance itself.
THIS = sys.modules[__name__]

THIS.klass = None

# Used in unitary test
def get_instance():
    """
        This function return the instance of the class
        ServiceLogger
    """
    return THIS.klass

class ServiceLogger(logging.getLoggerClass()):
    """
        The class ServiceLogger defines all logging parameter.
        It's an interface to python logging class.
        It's define as a singleton
    """
    instance = None
    def __new__(cls, name):
        if cls.instance is None:
            cls.instance = object.__new__(cls)
        return cls.instance

    def __init__(self, name):
        """
            Init class ServiceLogger
            :param name: name of the logger
        """
        THIS.klass = self
        cfg = serviceConfigFile.get_instance()
        # logging format
        # LEVEL : DEBUG, INFO, WARNING, ERROR
        # log format :
        # YYYY-MM-DDThh:mm:ss.mmm     LEVEL:ClassName:FunctionName: message
        self.logFormatter = logging.Formatter(fmt='%(asctime)s.%(msecs)03d     %(levelname)s:%(name)s::%(funcName)s: %(message)s', datefmt='%Y-%m-%dT%H:%M:%S')

        # set the name of the class in log messages
        self.rootLogger = logging.getLogger()
        # set the logging level from the configuration file
        self.rootLogger.setLevel(cfg.get('LOGGING', 'logFileLevel'))
        if not hasattr(self, 'first'):
            # First call to ServiceLogger
            self.first = True
            # create a log file
            self.fileHandler = logging.FileHandler(cfg.get('LOGGING', 'logFile'), mode='w')
            self.fileHandler.setFormatter(self.logFormatter)
            self.fileHandler.setLevel(cfg.get('LOGGING', 'logFileLevel'))
            # create a memory Handler to bufferize SAS log message
            self.memoryHandler = logging.handlers.MemoryHandler(1000, target=self.fileHandler)
            # add logger
            self.rootLogger.addHandler(self.memoryHandler)

            if cfg.get('LOGGING', 'logConsole') == 'True':
                # logging in console
                self.consoleHandler = logging.StreamHandler()
                self.consoleHandler.setFormatter(self.logFormatter)
                self.consoleHandler.setLevel(cfg.get('LOGGING', 'logConsoleLevel'))
                self.rootLogger.addHandler(self.consoleHandler)

    def flush_log(self):
        """
            This method flush the log file
        """
        self.memoryHandler.flush()

    def close(self):
        """
            This method close the logging file
        """
        self.fileHandler.close()
        self.rootLogger.removeHandler(self.fileHandler)
        self.consoleHandler.close()
        self.rootLogger.removeHandler(self.consoleHandler)
####################################################################
####################################################################


