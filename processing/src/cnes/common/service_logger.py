#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2019/05/17:version initiale.
# VERSION:2.0.0:DM:#91:2020/07/03:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: service_logger.py
    :synopsis: logging class for Swot
.. moduleauthor:: capgemini
..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
"""

import datetime
import os, sys
import logging
import logging.handlers
import cnes.common.service_config_file as service_config_file

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
    def __new__(cls):
        """
            __new__ method for class ServiceLogger
        """
        if cls.instance is None:
            cls.instance = object.__new__(cls)
        return cls.instance

    def __init__(self):
        """
            Init class ServiceLogger
        """
        THIS.klass = self
        cfg = service_config_file.get_instance()
        # logging level
        # LEVEL : DEBUG, INFO, WARNING, ERROR, SIGMSG
        # log format :
        # YYYY-MM-DDThh:mm:ss.mmm     LEVEL:ClassName:FunctionName: message
        self.log_formatter = logging.Formatter(fmt='%(asctime)s.%(msecs)03d     %(levelname)s | %(name)s::%(funcName)s | %(message)s',
                                               datefmt='%Y-%m-%dT%H:%M:%S')

        # set the name of the class in log messages
        self.root_logger = logging.getLogger()

        # Define SIGMSG level
        SIGMSG_LEVEL_NUM = 45
        logging.SIGMSG = SIGMSG_LEVEL_NUM
        level_name = "SIGMSG"
        logging.addLevelName(SIGMSG_LEVEL_NUM, level_name)
        def log_sigmsg(self, message, *args, **kwargs):
            """
               Local function sigmsg log
               :param message: log message
               :type message: string
               :param *args: number of argument
               :type *args: *int
               :param **kwargs: list of argument
               :type **kwargs: **type
            """
            self._log(logging.SIGMSG, message, args, **kwargs)
        logging.Logger.sigmsg = log_sigmsg

        # set the logging level from the configuration file
        self.root_logger.setLevel(cfg.get('LOGGING', 'logFileLevel'))
        if not hasattr(self, 'first'):
            # First call to ServiceLogger
            self.first = True
            # create a log file
            self.file_handler = logging.FileHandler(cfg.get('LOGGING', 'logFile'), mode='w')
            self.file_handler.setFormatter(self.log_formatter)
            self.file_handler.setLevel(cfg.get('LOGGING', 'logFileLevel'))
            # create the error log file
            try:
                # Try to get errorFile
                log_error_file = cfg.get('LOGGING', 'errorFile')
            except:
                # default name error.log
                log_error_filename = "error.log"
                log_error_file = os.path.join(cfg.get('PATHS', 'Output directory'), log_error_filename)
            self.file_handler_error = logging.FileHandler(log_error_file, mode='w')
            self.file_handler_error.setFormatter(self.log_formatter)
            self.file_handler_error.setLevel("ERROR")

            # create a memory Handler to bufferize SAS log message
            self.memory_handler = logging.handlers.MemoryHandler(1000, target=self.file_handler)
            # add logger
            self.root_logger.addHandler(self.memory_handler)
            self.root_logger.addHandler(self.file_handler_error)

            if cfg.get('LOGGING', 'logConsole') == 'True':
                # logging in console
                self.console_handler = logging.StreamHandler()
                self.console_handler.setFormatter(self.log_formatter)
                self.console_handler.setLevel(cfg.get('LOGGING', 'logConsoleLevel'))
                self.root_logger.addHandler(self.console_handler)
            else:
                self.console_handler = None
                    

    def flush_log(self):
        """
            This method flush the log file
        """
        self.memory_handler.flush()

    def close(self):
        """
            This method close the logging file
        """
        self.memory_handler.close()
        self.root_logger.removeHandler(self.memory_handler)
        self.file_handler.close()
        self.file_handler_error.close()
        self.root_logger.removeHandler(self.file_handler)
        self.root_logger.removeHandler(self.file_handler_error)
        if self.console_handler is not None:
            self.console_handler.close()
            self.root_logger.removeHandler(self.console_handler)
        del self.first
        self.instance = None
        THIS.klass = None

####################################################################
####################################################################


