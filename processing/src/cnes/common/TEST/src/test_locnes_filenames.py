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



import unittest

import cnes.common.lib_lake.locnes_filenames as locnes_filenames
import cnes.common.serviceError as serviceError

class TestLocnesFilenames(unittest.TestCase):
    """
        class for unitary test of locnes_filenames
    """
    def test_getInfoFromFilename(self):
        """
            unitary test for getInfoFromFilename
        """

        # Test PIXC
        filename = "SWOT_L2_HR_PIXC_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01.nc"
        filetype = "PIXC"
        retour = locnes_filenames.getInfoFromFilename(filename, filetype)
        # Test values returned
        ref_value = {'pass': '366', 'stop_date': '20140114T011056', 'tile_ref': '43N-R', 'crid': 'Dx0000', 'start_date': '20140114T011055', 'cycle': '000', 'counter': '01'}
        self.assertEqual(retour, ref_value)

        # Test PIXCVecRiver
        filename = "SWOT_L2_HR_PIXCVecRiver_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01.nc"
        filetype = "PIXCVecRiver"
        retour = locnes_filenames.getInfoFromFilename(filename, filetype)

        # Test LakeTile
        filename = "SWOT_L2_HR_LakeTile_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01_pixcvec.nc"
        filetype = "LakeTile"
        retour = locnes_filenames.getInfoFromFilename(filename, filetype)

        # Test error
        # Test bad key
        filename = "SWOT_L2_HR_LakeTile_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01_pixcvec.nc"
        filetype = "BadKey"
        with self.assertRaises(serviceError.SASLakeTileError):
            retour = locnes_filenames.getInfoFromFilename(filename, filetype)


    def test_lakeTileFilenames(self):
        """
            unitary test for lakeTileFilenames
        """

        # prepare input data
        pixc_file = "../data/SWOT_L2_HR_PIXC_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01.nc"
        pixc_vec_river_file = "../data/SWOT_L2_HR_PIXCVecRiver_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01.nc"
        # prepare output directory
        out_dir = "../output/"
        # Instanciate class (call constructor)
        objet = locnes_filenames.lakeTileFilenames(pixc_file, pixc_vec_river_file, out_dir)
        # Test values read in the filename
        self.assertEqual(objet.cycle_num, 0)
        self.assertEqual(objet.pass_num, 366)
        self.assertEqual(objet.tile_ref, "43N-R")
        self.assertEqual(objet.start_date, "20140114T011055")
        self.assertEqual(objet.stop_date, "20140114T011056")

        # test computeLakeTileFilename_shp method
        objet.computeLakeTileFilename_shp()
        self.assertEqual(objet.lake_tile_shp_file, "../output/SWOT_L2_HR_LakeTile_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01.shp")

        # test computeLakeTileFilename_edge
        objet.computeLakeTileFilename_edge()
        self.assertEqual(objet.lake_tile_edge_file, "../output/SWOT_L2_HR_LakeTile_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01_edge.nc")

        # test computeLakeTileFilename_pixcvec
        objet.computeLakeTileFilename_pixcvec()
        self.assertEqual(objet.lake_tile_pixcvec_file, "../output/SWOT_L2_HR_LakeTile_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01_pixcvec.nc")

        # test computeLakeIdPrefix
        objet.computeLakeIdPrefix()
        self.assertEqual(objet.lake_id_prefix, "5_000_366_43N-R_")

        # test testPixcVecRiverFilename
        retour = objet.testPixcVecRiverFilename()
        self.assertTrue(retour)

        # ERROR TEST
        # prepare input data with bad filename
        pixc_file = "../data/SWOT_L2_HR_PIXC_366_43N-R_20140114T011055_20140114T011056_Dx0000_01.nc"
        pixc_vec_river_file = "../data/SWOT_L2_HR_PIXCVecRiver_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01.nc"
        # prepare output directory
        out_dir = "../output/"
        # Instanciate class (call constructor)
        with self.assertRaises(serviceError.ProcessingError):
            objet = locnes_filenames.lakeTileFilenames(pixc_file, pixc_vec_river_file, out_dir)


    def test_lakeSPFilenames(self):
        """
            unitary test for lakeSPFilenames
        """

        # prepare input data
        lake_tile_list = ["../data/SWOT_L2_HR_LakeTile_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01_pixcvec.nc"]
        continent_value = "EU"
        # prepare output directory
        out_dir = "../output/"
        # Instanciate class (call constructor)
        objet = locnes_filenames.lakeSPFilenames(lake_tile_list, continent_value, out_dir)
        # Test values read in the filename
        self.assertEqual(objet.cycle_num, 0)
        self.assertEqual(objet.pass_num, 366)
        self.assertEqual(objet.start_date, "20140114T011055")
        self.assertEqual(objet.stop_date, "20140114T011056")

        # test
        objet.computeStartStopDates()
        self.assertEqual(objet.stop_date, "20140114T011056")
        self.assertEqual(objet.start_date, "20140114T011055")

        # test
        objet.computeLakeSPFilename()
        self.assertEqual(objet.lake_sp_file, "../output/SWOT_L2_HR_LakeSP_000_366_EU_20140114T011055_20140114T011056_Dx0000_01.shp")

        # test
        objet.computeLakeIdPrefix()
        self.assertEqual(objet.lake_id_prefix, "5_000_366_EU_")


        # ERROR TEST
        # prepare input data with bad filename
        lake_tile_list = ["../data/SWOT_L2_HR_LakeTile_366_43N-R_20140114T011055_20140114T011056_Dx0000_01_pixcvec.nc"]
        continent_value = "EU"
        # prepare output directory
        out_dir = "../output/"
        # Instanciate class (call constructor)
        with self.assertRaises(serviceError.ProcessingError):
            objet = locnes_filenames.lakeSPFilenames(lake_tile_list, continent_value, out_dir)



    def test_computePixcvecFilename(self):
        """
            unitary test for computePixcvecFilename
        """

        # prepare input data
        pixcvec_file = "../data/SWOT_L2_HR_LakeTile_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01_pixcvec.nc"
        # prepare output directory
        out_dir = "../output/"
        # Call
        out_file = locnes_filenames.computePixcvecFilename(pixcvec_file, out_dir)
        self.assertEqual(out_file, "../output/SWOT_L2_HR_PIXCVec_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01.nc")

        # ERROR TEST
        # prepare input data with bad filename
        pixcvec_file = "../data/SWOT_L2_HR_LakeTile_366_20140114T011055_20140114T011056_Dx0000_01_pixcvec.nc"
        # prepare output directory
        out_dir = "../output/"
        # Call and error test
        with self.assertRaises(serviceError.ProcessingError):
            out_file = locnes_filenames.computePixcvecFilename(pixcvec_file, out_dir)


