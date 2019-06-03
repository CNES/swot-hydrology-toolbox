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
.. module:: sas_lake_sp
   :synopsis: Process SAS_L2_HR_LakeSP, generate L2_HR_LakeSP and L2_HR_PIXCVec product from all tiles of L2_HR_LakeTile
    Created on 2017/02/27

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR
                  Cécile Cazals - C-S
..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import cnes.common.lib.my_timer as my_timer
import cnes.common.service_error as service_error
import numpy as np


class SASLakeSP(object):
    """
        SAS lake sp class
    """

    def __init__(self,  in_obj_pixc_sp, in_obj_pixc_vec_sp, in_obj_lake_db, in_obj_lake_l, in_obj_lake_r):
        """
        Constructor: initialize variables
        TODO : redefine inputs
        :param in_obj_pixc_sp: PIXC_SP object
        :type in_obj_pixc_sp: lake_sp.PixC_SP
        :param in_obj_pixc_vec_sp: PIXCVec_SP object
        :type in_obj_pixc_vec_sp: lake_sp.PixC_Vec_SP
        :param in_obj_lake_db: lake_db object
        :type in_obj_lake_db: lib_lake.lake_db
        :param in_obj_lake: proc_lake object
        :type in_obj_lake: lib_lake.proc_lake
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("")

        # Objects
        self.obj_lake_db = in_obj_lake_db  # Lake DB object
        self.obj_pixc_sp = in_obj_pixc_sp  # PIXC object
        self.obj_pixc_vec_sp = in_obj_pixc_vec_sp  # PIXCVecRiver object
        self.obj_lake_l = in_obj_lake_l  # LakeTile object
        self.obj_lake_r = in_obj_lake_r  # LakeTile object

    def run_preprocessing(self):

        """
        preprocessing lake_sp
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.sigmsg("")
        logger.sigmsg("**************************")
        logger.sigmsg("***** PRE-PROCESSING *****")
        logger.sigmsg("**************************")
        logger.sigmsg("")

        logger.info("NOTHING TO DO")


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

        try :


            if self.obj_pixc_sp.pixc_edge_l.nb_pixels > 0:
                logger.info("***** Processing Left Swath ******")

                # 3.1.1 - Gather pixels by entities for all the tiles of this swath
                logger.info("1 - Relabeling edge entities ...")
                self.obj_pixc_sp.pixc_edge_l.swath_global_relabeling()

                # 3.1.2 - Compute lake products
                logger.info("2 - Filling LakeTile product...")
                self.obj_lake_l.computeLakeProducts(np.unique(self.obj_pixc_sp.pixc_edge_l.labels))

            else:
                logger.info("No pixel to process")

            if self.obj_pixc_sp.pixc_edge_r.nb_pixels > 0:
                logger.info("***** Processing Right Swath ******")

                # 3.1.1 - Gather pixels by entities for all the tiles of this swath
                logger.info("1 - Relabeling edge entities ...")
                self.obj_pixc_sp.pixc_edge_r.swath_global_relabeling()

                # 3.1.2 - Compute lake products
                logger.info("2 - Filling LakeTile product...")
                self.obj_lake_r.computeLakeProducts(np.unique(self.obj_pixc_sp.pixc_edge_r.labels))

            else:

                logger.info("No pixel to process")


        except:
            message = "Something wrong happened in run_processing"
            raise service_error.SASLakeSpError(message, logger)

    def run_postprocessing(self):
        """
            Process PGE_L2_HR_LakeSP OUT
            Nothing to do in lake_SP
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.sigmsg("")
        logger.sigmsg("***************************")
        logger.sigmsg("***** POST-PROCESSING *****")
        logger.sigmsg("***************************")
        logger.sigmsg("")

#######################################
