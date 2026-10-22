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
# VERSION:5.0.0:DM:#91:2022/12/08:Poursuite industrialisation
# VERSION:6.2.1:FA:#822:2025/04/08:Nouveau patch Lake_AVG_V4.2.1
# VERSION:6.2.4_SSALTO-330:FA:#831:2026/02/26: SSALTO-330 [PGE][SAS LAKES] Mise à jour des TU (SSALPGE-14313)
# VERSION:6.2.6:DM:#900:2026/07/06: [PGE] [LAKES] Ajout d'une interface niveau moyen des lacs
# VERSION:6.2.6:DM:#900:2026/07/06:SSALTO-565 [PGE] [LAKES] Ajout d'une interface coeff NN pour calcul ice_dyn_f
# VERSION:6.2.7:DM:#906:2026/07/28:SSALTO-872 [PGE][LAKES] Gestion des erreurs/warning GDAL
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
import multiprocessing
import logging
import logging.handlers
import os
from osgeo import gdal, ogr, osr
import resource
import sys

import cnes.common.service_config_file as service_config_file


# pointer to the module object instance itself.
THIS = sys.modules[__name__]

THIS.klass = None

logcounter = { 'SIGMSG':0, 'ERROR':0, 'WARNING':0, 'INFO':0, 'DEBUG':0}

STACK_LEVEL = 2

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
        
        # INSERT ADT: to catch GDAL messages
        gdal.UseExceptions()
        ogr.UseExceptions()
        osr.UseExceptions()
        # INSERT ADT - End

        THIS.klass = self
        cfg = service_config_file.get_instance()
        self.nb_proc = nb_proc
        self.sys_time = 0
        self.user_time = 0
        self.total_mem = 0
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
            kwargs["stacklevel"] = STACK_LEVEL
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
            kwargs["stacklevel"] = STACK_LEVEL
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
            kwargs["stacklevel"] = STACK_LEVEL
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
            kwargs["stacklevel"] = STACK_LEVEL
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
            kwargs["stacklevel"] = STACK_LEVEL
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
                dir_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "TEST/output")
                log_error_filename = "error.log"
                log_error_file = os.path.join(dir_path, cfg.get('PATHS', 'Output directory'), log_error_filename)
            self.file_handler_error = logging.FileHandler(log_error_file, mode='w')
            self.file_handler_error.setFormatter(self.log_formatter)
            self.file_handler_error.setLevel("ERROR")

            # create a memory Handler to bufferize SAS log message
            # Pb memory handler self.memory_handler = logging.handlers.MemoryHandler(1000, target=self.file_handler)
            # add logger
            # Pb memory handler self.root_logger.addHandler(self.memory_handler)
            self.root_logger.addHandler(self.file_handler)
            self.root_logger.addHandler(self.file_handler_error)

            if cfg.get('LOGGING', 'logConsole') == 'True':
                # logging in console
                self.console_handler = logging.StreamHandler()
                self.console_handler.setFormatter(self.log_formatter)
                self.console_handler.setLevel(cfg.get('LOGGING', 'logConsoleLevel'))
                self.root_logger.addHandler(self.console_handler)
            else:
                self.console_handler = None

    def set_stat(self, sys_time, user_time, total_mem):
        """
           Local function debug log

           :param sys_time: sys_time of subprocess
           :type sys_time: double
           :param user_time: user_time of subprocess
           :type user_time: double
           :param total_mem: total_mem of subprocess
           :type total_mem: double
        """
        self.sys_time = sys_time
        self.user_time = user_time
        self.total_mem = total_mem
 

    def flush_log(self):
        """
            This method flush the log file
        """
        # Pb memory handler self.memory_handler.flush()
        pass

    def close(self):
        """
            This method close the logging file
        """
        try:
            get_stat(self.nb_proc, self.sys_time, self.user_time, self.total_mem)
        except:
            pass
        # Pb memory handler self.memory_handler.close()
        # Pb memory handler self.root_logger.removeHandler(self.memory_handler)
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


def get_stat(nb_proc, in_sys_time, in_user_time, in_total_mem):
    """
        This method get and log stats

        :param nb_proc: Number of proc used
        :type message: int
    """

    logger = logging.getLogger("service_logger")
    # lecture de /proc/self/status pour avoir l'occupation mémoire
    stat_file = "/proc/self/status"
    stat_fp = open(stat_file, 'r')
    for ligne in stat_fp:
        if ("VmPeak:" in ligne):
            logger.sigmsg(ligne.replace("VmPeak:", "peak_vm:").rstrip('\n'))
#        if ("VmHWM:" in ligne):
#            logger.sigmsg(ligne.replace("VmHWM:", "max_rss").rstrip('\n'))
    stat_fp.close()

    # lecture de /proc/self/stat pour avoir les temps user et cpu
    stat_file = "/proc/self/stat"
    stat_fp = open(stat_file, 'r')
    content = stat_fp.readlines()
    ligne = content[0].split()
    stat_fp.close()
    # Calcul user et sys time
    user_time = float(ligne[13])/100.0
    sys_time = float(ligne[14])/100.0
    # Calcul date start time
    process_start_time = float(ligne[21])/100.0
    # lecture de /proc/uptime pour récupérer le temps système
    stat_file = "/proc/uptime"
    stat_fp = open(stat_file, 'r')
    content = stat_fp.readlines()
    sys_uptime_sec = float(content[0].split()[0])
    stat_fp.close()
    # log
    resource_self = resource.getrusage(resource.RUSAGE_SELF)
    resource_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    sys_time = resource_self.ru_stime + resource_children.ru_stime + in_sys_time
    user_time = resource_self.ru_utime + resource_children.ru_utime + in_user_time
    total_mem = resource_self.ru_maxrss + resource_children.ru_maxrss + in_total_mem
    logger.sigmsg("max_rss: " + str(total_mem) + " kB")
    logger.sigmsg("sys_cpu: " + str(sys_time) + " s")
    logger.sigmsg("user_cpu: " + str(user_time) + " s")
    logger.sigmsg("total_run_time: " + str(round(float(sys_uptime_sec) - process_start_time, 2)) + " s")
    logger.sigmsg("nb_core: " + str(nb_proc))
    logger.sigmsg("sigmsg: " + str(logcounter['SIGMSG']))
    logger.sigmsg("error: " + str(logcounter['ERROR']))
    logger.sigmsg("warning: " + str(logcounter['WARNING']))
    logger.sigmsg("info: " + str(logcounter['INFO']))
    logger.sigmsg("debug: " + str(logcounter['DEBUG']))





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
        self._logger = logger
        self._log_listener = multiprocessing.Process(
            target=self.listener_process,
            args=(self._queue, self._logger),
        )

    def listener_process(self, queue, logger):
        """
            This function is the listener process

            :param queue: Number of proc used
            :type queue: int
            :param logger: logger instance
            :type logger: logger class
        """
        root = logging.getLogger()
        handler = logging.StreamHandler()
        root.addHandler(handler)
        root.setLevel(logger.level)

        while True:
            try:
                record_data = queue.get()
                if record_data is None:
                    break

                record = logging.makeLogRecord(record_data) # ajout
                logger = logging.getLogger(record.name) # ajout

                logger.handle(record)
            except (KeyboardInterrupt, SystemExit):
                #pass
                raise # ajout
            except Exception:  # as exc:
                import sys, traceback
                print('Error in logger process', file=sys.stderr)
                traceback.print_exc(file=sys.stderr)
                break

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

    def worker_configurer(self, name, main=False):
        """
            Create a new log worker
            If main is True, the name of the log worker is NAME.PROCESS
            NAME:  MPLogger's name (self.name)
            PROCESS: Name of the process.
            Else the name of the log worker is self.name

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
        root.setLevel(self._logger.level)
        queue_handler = logging.handlers.QueueHandler(multiprocessing.Queue(-1))
        root.addHandler(queue_handler)

        return root




####################################################################
####################################################################


