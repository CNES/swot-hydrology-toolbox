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
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: sas_lake_tile
   :synopsis: Process SAS_L2_HR_LakeTile, i.e. generate L2_HR_LakeTile product from one tile of L2_HR_PIXC product and its associated L2_HR_PIXCVecRiver product
    Created on 2017/02/27

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR
..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import logging

import cnes.common.service_error as service_error

import cnes.common.lib.my_timer as my_timer


class SASLakeTile(object):
    """
    Class handling LakeTile SAS
    Dedicated to manage main processing steps
    """
    
    def __init__(self, in_obj_pixc, in_obj_pixc_vec, in_obj_lake_db, in_obj_lake):
        """
        Constructor: initialize variables

        :param in_obj_pixc: PIXC object
        :type in_obj_pixc: lake_tile.proc_pixc
        :param in_obj_pixc_vec: PIXCVecRiver object
        :type in_obj_pixc_vec: lib_lake.proc_pixc_vec
        :param in_obj_lake_db: Prior Lake Database (PLD) object
        :type in_obj_lake_db: lib_lake.lake_db
        :param in_obj_lake: lake product (LakeTile) object
        :type in_obj_lake: lib_lake.proc_lake
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("")

        # Objects
        self.obj_lake_db = in_obj_lake_db  # Lake DB object
        self.obj_pixc = in_obj_pixc  # PIXC object
        self.obj_pixc_vec = in_obj_pixc_vec  # PIXCVecRiver object
        self.obj_lake = in_obj_lake  # LakeTile object

    def run_preprocessing(self):

        """
        Process LakeTile preprocessing =
        - reshape PIXCVecRiver arrays to the same format as PIXC arrays
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.sigmsg("")
        logger.sigmsg("**************************")
        logger.sigmsg("***** PRE-PROCESSING *****")
        logger.sigmsg("**************************")
        logger.sigmsg("")

        try:
            # 1 - Reshape PIXCVecRiver arrays
            logger.info("> 1 - Reshape PIXCVecRiver arrays...")
            self.obj_pixc_vec.reshape(self.obj_pixc)
            
        except:
            message = "[lakeTileProcessing]   Something wrong happened in run_preprocessing"
            raise service_error.SASLakeTileError(message, logger)

    def run_processing(self):
        """
        Process SAS_L2_HR_LakeTile
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.sigmsg("")
        logger.sigmsg("**********************")
        logger.sigmsg("***** PROCESSING *****")
        logger.sigmsg("**********************")
        logger.sigmsg("")
        timer_proc = my_timer.Timer()
        timer_proc.start()

        try:
            # Processing only if PixC pixels are selected
            if self.obj_pixc.nb_selected != 0:

                # 2 - F2-F3-F3b = Identify all separate entities in the water mask
                logger.info("1 - Identifying all separate entities in the water mask...")
                self.obj_pixc.compute_separate_entities()
                logger.info("" + timer_proc.info(0))
                logger.info("")


                # 3 - F4 = Retrieve pixels corresponding to lakes and unknown entirely inside the tile
                logger.info("2 - Getting pixels corresponding to lakes and unknown entirely inside the tile...")
                self.obj_pixc.compute_obj_inside_tile()
                logger.info("" + timer_proc.info(0))
                logger.info("")

            
                # 4 - F5 = Retrieve pixels indices and associated label of objects at the top/bottom edge of the tile
                logger.info("3 - Getting pixels corresponding to objects at the top/bottom edge of the tile...")
                self.obj_pixc.compute_edge_indices_and_label()
                logger.info("" + timer_proc.info(0))
                logger.info("")

                # 5 - F6 = Fill lake product
                logger.info("4 - Filling LakeTile product...")
                self.obj_lake.compute_lake_products(self.obj_pixc.labels_inside)
                logger.info("" + timer_proc.info(0))
                logger.info("")

            else:
                logger.info("NO selected PixC => empty lake tile product generated")
                logger.info("")

        except:
            message = "Something wrong happened in run_processing"
            raise service_error.SASLakeTileError(message, logger)

    def run_postprocessing(self):
        """
        Process LakeTile postprocessing = 
        - Nothing to do in lake_tile
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.sigmsg("")
        logger.sigmsg("***************************")
        logger.sigmsg("***** POST-PROCESSING *****")
        logger.sigmsg("***************************")
        logger.sigmsg("")
        