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
# VERSION:3.0.0:DM:#91:2021/03/12:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: sas_lake_sp
   :synopsis: Process SAS_L2_HR_LakeSP, i.e. generate L2_HR_LakeSP and L2_HR_PIXCVec product from all tiles of L2_HR_LakeTile
    Created on 2017/02/27

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR
                  Cécile CAZALS - C-S
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
    Class handling LakeSP SAS
    Dedicated to manage main processing steps
    """

    def __init__(self,  in_obj_pixc_sp, in_obj_pixcvec_sp, in_obj_lake_db, in_obj_lake):
        """
        Constructor: initialize variables
        
        :param in_obj_pixc_sp: object grouping LakeTile_edge data
        :type in_obj_pixc_sp: proc_pixc_sp.PixCEdge
        :param in_obj_pixcvec_sp: object grouping LakeTile_pixcvec data
        :type in_obj_pixcvec_sp: proc_pixc_vec_sp.PixCVecSP
        :param in_obj_lake_db: Prior Lake Database (PLD) object
        :type in_obj_lake_db: lib_lake.lake_db
        :param in_obj_lake: LakeSP product object
        :type in_obj_lake: lib_lake.proc_lake
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("")

        # Objects
        self.obj_lake_db = in_obj_lake_db  # Lake DB object
        self.obj_pixc_sp = in_obj_pixc_sp  # LakeTile_edge object
        self.obj_pixcvec_sp = in_obj_pixcvec_sp  # LakeTile_pixcvec object
        self.obj_lake = in_obj_lake  # LakeSP object

    def run_preprocessing(self):
        """
        Process LakeSP pre-processing = nothing to do
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.sigmsg("")
        logger.sigmsg("==========================")
        logger.sigmsg("===== PRE-PROCESSING =====")
        logger.sigmsg("==========================")
        logger.sigmsg("")
        logger.info("NOTHING TO DO")

    def run_processing(self):
        """
        Process SAS_L2_HR_LakeSP
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.sigmsg("")
        logger.sigmsg("======================")
        logger.sigmsg("===== PROCESSING =====")
        logger.sigmsg("======================")
        logger.sigmsg("")
        timer_proc = my_timer.Timer()
        timer_proc.start()

        try :
            
            # 1 - Processing left swath if there are pixels
            logger.info("***** Processing Left Swath ******")
            
            if self.obj_pixc_sp.pixc_edge_l.nb_pixels > 0:

                # 1.1 - Gather edge pixels in separate entities for all the tiles of this swath
                logger.info("1 - Gathering edge pixels in separate entities for all the tiles of this swath...")
                self.obj_pixc_sp.pixc_edge_l.swath_global_relabeling()
                logger.info("" + timer_proc.info(0))
                logger.info("")

                # 1.2 - Compute LakeSP product for this swath
                logger.info("2 - Computing LakeSP features for this swath...")
                self.obj_lake.swath_l.compute_lake_features(np.unique(self.obj_pixc_sp.pixc_edge_l.labels))
                logger.info("" + timer_proc.info(0))
                logger.info("")

            else:
                logger.info("No pixel to process for this swath")

            # 2 - Processing right swath if there are pixels
            logger.info("***** Processing Right Swath ******")
            
            if self.obj_pixc_sp.pixc_edge_r.nb_pixels > 0:

                # 2.1 - Gather edge pixels in separate entities for all the tiles of this swath
                logger.info("1 - Gathering edge pixels in separate entities for all the tiles of this swath...")
                self.obj_pixc_sp.pixc_edge_r.swath_global_relabeling()
                logger.info("" + timer_proc.info(0))
                logger.info("")

                # 2.2 - Compute LakeSP product for this swath
                logger.info("2 - Computing LakeSP features for this swath...")
                self.obj_lake.swath_r.compute_lake_features(np.unique(self.obj_pixc_sp.pixc_edge_r.labels))
                logger.info("" + timer_proc.info(0))
                logger.info("")

            else:
                logger.info("No pixel to process for this swath")

        except:
            message = "Something wrong happened in run_processing"
            raise service_error.SASLakeSpError(message, logger)

    def run_postprocessing(self):
        """
        Process LakeSP post-processing = nothing to do
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.sigmsg("")
        logger.sigmsg("===========================")
        logger.sigmsg("===== POST-PROCESSING =====")
        logger.sigmsg("===========================")
        logger.sigmsg("")
        logger.info("NOTHING TO DO")
