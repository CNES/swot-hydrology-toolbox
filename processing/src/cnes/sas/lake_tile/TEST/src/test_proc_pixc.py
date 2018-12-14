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

import cnes.common.lib_lake.proc_lake as proc_lake
import cnes.common.lib_lake.locnes_filenames as locnes_filenames
import cnes.common.lib_lake.proc_pixc_vec as proc_pixc_vec
import cnes.sas.lake_tile.proc_pixc as proc_pixc
import cnes.common.lib_lake.lake_db as lake_db
import cnes.common.lib_lake.locnes_variables as locnes_variables

import cnes.common.serviceConfigFile as serviceConfigFile
import cnes.common.serviceLogger as serviceLogger
import cnes.common.serviceError as serviceError

class TestPixelCloud(unittest.TestCase):
  """
      class for unitary test of PixelCloud module
  """
  def test_initPixelCloud(self):
    """
        unitary test for __init__ from PixelCloud class
    """


    # Prepare data
    pixc_file = "../data/SWOT_L2_HR_PIXC_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01.nc"
    pixc_vec_river_file = "../data/SWOT_L2_HR_PIXCVecRiver_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01.nc"
    # prepare output directory
    out_dir = "../output/"
    # get objetlakeTileFilenames
    objetlakeTileFilenames = locnes_filenames.lakeTileFilenames(pixc_file, pixc_vec_river_file, out_dir)

    locnes_variables.CONTINENT_FILE = os.environ['REF_PATH'] + "/../BD/major_basins/FAO/major_hydrobasins.shp"
    # load objPixcVec
    objPixcVec = proc_pixc_vec.PixelCloudVec("TILE", pixc_vec_river_file)
    # load objPixc
    objPixc = proc_pixc.PixelCloud(pixc_file, objPixcVec.reject_idx)
    print(objPixc.labels)




