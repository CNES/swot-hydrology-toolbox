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
# VERSION:4.0.0:DM:#91:2022/05/05:Poursuite industrialisation
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
import multiprocessing

# pointer to the module object instance itself.
THIS = sys.modules[__name__]

THIS.klass = None

logcounter = { 'SIGMSG':0, 'ERROR':0, 'WARNING':0, 'INFO':0, 'DEBUG':0}

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
    def __new__(cls, nb_proc = 1):
        """
            __new__ method for class ServiceLogger
        """
        if cls.instance is None:
            cls.instance = object.__new__(cls)
        return cls.instance

    def __init__(self, nb_proc = 1):
        """
            Init class ServiceLogger
        """
        THIS.klass = self
        cfg = service_config_file.get_instance()
        self.nb_proc = nb_proc
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
            global logcounter
            logcounter['SIGMSG'] += 1
        logging.Logger.sigmsg = log_sigmsg

        def log_error(self, message: str, *args, **kwargs):
            """
               Local function error log

               :param message: log message
               :type message: string
               :param *args: number of argument
               :type *args: *int
               :param **kwargs: list of argument
               :type **kwargs: **type
            """
            self._log(logging.ERROR, message, args, **kwargs)
            global logcounter
            logcounter['ERROR'] += 1
        logging.Logger.error = log_error

        def log_warning(self, message: str, *args, **kwargs):
            """
               Local function warning log

               :param message: log message
               :type message: string
               :param *args: number of argument
               :type *args: *int
               :param **kwargs: list of argument
               :type **kwargs: **type
            """
            self._log(logging.WARNING, message, args, **kwargs)
            global logcounter
            logcounter['WARNING'] += 1
        logging.Logger.warning = log_warning

        def log_info(self, message: str, *args, **kwargs):
            """
               Local function info log

               :param message: log message
               :type message: string
               :param *args: number of argument
               :type *args: *int
               :param **kwargs: list of argument
               :type **kwargs: **type
            """
            self._log(logging.INFO, message, args, **kwargs)
            global logcounter
            logcounter['INFO'] += 1
        logging.Logger.info = log_info

        def log_debug(self, message: str, *args, **kwargs):
            """
               Local function debug log

               :param message: log message
               :type message: string
               :param *args: number of argument
               :type *args: *int
               :param **kwargs: list of argument
               :type **kwargs: **type
            """
            self._log(logging.DEBUG, message, args, **kwargs)
            global logcounter
            logcounter['DEBUG'] += 1
        logging.Logger.debug = log_debug


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
        get_stat(self.nb_proc)
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


def get_stat(nb_proc):
    """
        This method get and log stats

        :param nb_proc: Number of proc used
        :type message: int
    """

    logger = logging.getLogger("service_logger")
    stat_file = "/proc/self/status"
    stat_fp = open(stat_file, 'r')
    for ligne in stat_fp:
        if ("VmPeak:" in ligne):
            logger.sigmsg(ligne.replace("VmPeak:", "peak_vm").rstrip('\n'))
        if ("VmHWM:" in ligne):
            logger.sigmsg(ligne.replace("VmHWM:", "max_rss").rstrip('\n'))
    stat_fp.close()

    logger.sigmsg("nb_core " + str(nb_proc))
    logger.sigmsg("sigmsg " + str(logcounter['SIGMSG']))
    logger.sigmsg("error " + str(logcounter['ERROR']))
    logger.sigmsg("warning " + str(logcounter['WARNING']))
    logger.sigmsg("info " + str(logcounter['INFO']))
    logger.sigmsg("debug " + str(logcounter['DEBUG']))


def _listener_process(queue, logger):
    """
        This function is the listener process

        :param queue: Number of proc used
        :type queue: int
        :param logger: logger instance
        :type logger: logger class
    """

    while True:
        try:
            record = queue.get()
            if record is None:
                break
            logger.handle(record)
        except (KeyboardInterrupt, SystemExit):
            pass
        except Exception:  # as exc:
            #  import sys, traceback
            #  print('Whoops! LOG Problem:', file=sys.stderr)
            #  traceback.print_exc(file=sys.stderr)
            break

class MPLogger(object):
    """
        The class MPLogger manage parallel logging.
    """
    def __init__(self, logger):
        """
            __init__ method of MPLogger class

            :param logger: logger instance
            :type logger: logger class
        """

        self._queue = multiprocessing.Queue(-1)
        if not logger.handlers:
            h = logging.handlers.QueueHandler(self._queue)
            logger.addHandler(h)
            logger.setLevel(logger.level)

        self._loger = logger
        self._log_listener = multiprocessing.Process(
            target=_listener_process,
            args=(self._queue, self._loger),
        )

    def start_listenner(self):
        """
            start_listenner method
        """

        if self._log_listener is None:
            raise RuntimeError("Listenner has already been started and stopped !")
        self._log_listener.start()

    def stop_listenner(self):
        """
            stop_listenner method
        """

        self._queue.put_nowait(None)
        self._log_listener.join()
        self._log_listener = None

    def worker_configurer(self, name, main=False):
        """
            Create a new log worker
            If main is True, the name of the log worker is NAME.PROCESS
                NAME:  MPLogger's name (self.name)
                PROCESS: Name of the process.
            Else the the name of the log worker is self.name

            :param name: name of the worker
            :type name: string
            :param main: logger instance
            :type main: boolean
        """
        if not main:
            pname = multiprocessing.current_process().name
            pname = "%s.%s" % (name, pname)
        else:
            pname = name

        root = logging.getLogger(pname)

        if not root.handlers:
            h = logging.handlers.QueueHandler(self._queue)
            root.addHandler(h)
            root.setLevel(self._loger.level)
        return root




####################################################################
####################################################################


