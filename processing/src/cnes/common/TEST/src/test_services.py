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

import sys
import os
import os.path
import shutil
import logging
import unittest

import cnes.common.serviceConfigFile as serviceConfigFile
import cnes.common.serviceError as serviceError
import cnes.common.serviceLogger as serviceLogger

class TestServices(unittest.TestCase):

  def setUp(self):
    # define configuration file for test
    #fileConfig='../data/testConfig.cfg'
    fileConfig=sys.path[0]+'/../data/testConfig.cfg'
    self.assertTrue(os.path.isfile(fileConfig))

    # Read the configuration file
    cfg = serviceConfigFile.ServiceConfigFile(fileConfig)

  def testserviceConfigFile(self):

    """
        unitary test for serviceConfigFile
    """

    cfg = serviceConfigFile.get_instance()
    # test values read in the file
    # Chain test
    self.assertEqual(cfg.get('main', 'chaineCaractere'), 'Une chaine de caractere')
    # Word test
    self.assertEqual(cfg.get('main', 'unMot'), 'mot')
    # integer test
    self.assertEqual(cfg.getint('main', 'entier'), 1)
    # float test
    self.assertEqual(cfg.getfloat('main', 'flottant'), 1.2)
    #boolean test
    #booleen1=True
    self.assertTrue(cfg.getboolean('main', 'booleen1'))
    #booleen2=Yes
    self.assertTrue(cfg.getboolean('main', 'booleen2'))
    #booleen3=Y
    with self.assertRaises(ValueError):
        cfg.getboolean('main', 'booleen3')
    #booleen4=1
    self.assertTrue(cfg.getboolean('main', 'booleen4'))


    # Test check methods
    # Test wrong section
    with self.assertRaises(serviceError.ConfigFileError):
        cfg.test_var_config_file('badSection', 'entier', int)
    # Test wrong variable
    with self.assertRaises(serviceError.ParameterError):
        cfg.test_var_config_file('main', 'badVariable', int)
    # Test good and wrong type
    cfg.test_var_config_file('main', 'entier', int)
    with self.assertRaises(serviceError.ParameterError):
        cfg.test_var_config_file('main', 'flottant', int)
    # Test list of value (good and bad)
    cfg.test_var_config_file('main', 'entier', int, valeurs=[1, 2, 3, 4])
    with self.assertRaises(serviceError.ParameterError):
        cfg.test_var_config_file('main', 'entier', int, valeurs=[2, 3, 4])
    # Test test missing variable with default value
    cfg.test_var_config_file('main', 'missingVariable', int, val_defaut=42)
    self.assertEqual(cfg.getint('main', 'missingVariable'), 42)
    # Test valid_min and valid_max (good and bad)
    # test valeur dans l'interval
    cfg.test_var_config_file('main', 'flottant', float, valid_min=0.5, valid_max=1.5)
    # test valeur trop petite
    with self.assertRaises(serviceError.ParameterError):
        cfg.test_var_config_file('main', 'flottant', float, valid_min=1.3, valid_max=1.5)
    # test valeur trop grande
    with self.assertRaises(serviceError.ParameterError):
        cfg.test_var_config_file('main', 'flottant', float, valid_min=0.5, valid_max=0.9)

  def testserviceLogging(self):

    """
        unitary test for serviceLogging
    """

    # get instance of configuration file for test
    cfg = serviceConfigFile.get_instance()
    ## Starting of logging service
    # First call, set name, level and type of log
    serviceLogger.ServiceLogger(self.__class__.__name__)

    ## Local instanciation of logging
    # Not a first call, only name can be overlord
    logger = logging.getLogger(self.__class__.__name__)
    # Write some log messages
    logger.info("START of testserviceLogging")
    # DEBUG log
    logger.debug("This is DEBUG log")
    # INFO log
    logger.info("This is INFO log")
    # WARNING log
    logger.warning("This is WARNING log")
    # ERROR log
    logger.error("This is ERROR log")
    logger.info("END of testserviceLogging")

    # Test if logfile is initiate
    self.assertTrue(os.path.isfile(cfg.get('LOGGING', 'logFile')))

    # Flush of log file
    instance = serviceLogger.get_instance()
    instance.memoryHandler.flush()
    # Openning log file
    fichierLog = open(cfg.get('LOGGING', 'logFile'), "r")
    # Verify containt of log file
    lines = fichierLog.readlines()
    # Test message format in log file
    self.assertTrue(("INFO:TestServices::testserviceLogging: START of testserviceLogging" in lines[0]))
    self.assertTrue(("DEBUG:TestServices::testserviceLogging: This is DEBUG log" in lines[1]))
    self.assertTrue(("INFO:TestServices::testserviceLogging: This is INFO log" in lines[2]))
    self.assertTrue(("WARNING:TestServices::testserviceLogging: This is WARNING log" in lines[3]))
    self.assertTrue(("ERROR:TestServices::testserviceLogging: This is ERROR log" in lines[4]))
    self.assertTrue(("INFO:TestServices::testserviceLogging: END of testserviceLogging" in lines[5]))
    

  def testserviceError(self):

    """
        unitary test for serviceError
    """
    # useless ??
    logger = logging.getLogger(self.__class__.__name__)
    # Test if we can raise SASLakeTileError
    with self.assertRaises(serviceError.SASLakeTileError):
       msg = "Test crash with SASLakeTileError"
       raise serviceError.SASLakeTileError(msg, logger)

    # Test if we can raise SASLakeSpError
    with self.assertRaises(serviceError.SASLakeSpError):
       msg = "Test crash with SASLakeSpError"
       raise serviceError.SASLakeSpError(msg, logger)

#  def tearDown(self):
 #   self.assertTrue(False)


