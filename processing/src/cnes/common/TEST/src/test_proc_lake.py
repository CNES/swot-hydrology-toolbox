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



import os
import os.path
import unittest

import cnes.common.lib_lake.proc_lake as proc_lake
import cnes.common.lib_lake.locnes_filenames as locnes_filenames
import cnes.common.lib_lake.proc_pixc_vec as proc_pixc_vec
import cnes.sas.lake_tile.proc_pixc as proc_pixc
import cnes.common.lib_lake.lake_db as lake_db
import cnes.common.lib_lake.locnes_variables as locnes_variables

class TestProcLake(unittest.TestCase):
    """
        class for unitary test of proc_lake module
    """
    def test_LakeProduct(self):
        """
            unitary test for LakeProduct class
        """

        # Prepare data
        pixc_file = "../data/SWOT_L2_HR_PIXC_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01.nc"
        pixc_vec_river_file = "../data/SWOT_L2_HR_PIXCVecRiver_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01.nc"
        # prepare output directory
        out_dir = "../output/"
        # get objetlakeTileFilenames
        objetlakeTileFilenames = locnes_filenames.lakeTileFilenames(pixc_file, pixc_vec_river_file, out_dir)

        #locnes_variables.overwriteParameter("CONTINENT_FILE", os.environ['REF_PATH'] + "../BD/major_basins/FAO/major_hydrobasins.shp")
        locnes_variables.CONTINENT_FILE = os.environ['REF_PATH'] + "/../BD/major_basins/FAO/major_hydrobasins.shp"
        # load objPixcVec
        objPixcVec = proc_pixc_vec.PixelCloudVec("TILE", pixc_vec_river_file)
        # load objPixc
        objPixc = proc_pixc.PixelCloud(pixc_file, objPixcVec.reject_idx)
        print(objPixc.labels)
        # load objLakeDb
        lakeDb_file = os.environ['REF_PATH'] + "/../BD/BD_lakes/20180223_France/priordb_lakes_france.shp"
        objLakeDb = lake_db.LakeDb_shp(lakeDb_file, objPixc.tile_poly)

        # Data for lakeProduct
        productType = "TILE"
        layer_name = "SWOT_L2_HR_LakeTile_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01"
        id_prefix = "5_000_366_43N-R_"
        # Instanciate class
        objet = proc_lake.LakeProduct(productType, objPixc, objPixcVec, objLakeDb, layer_name, id_prefix)

        objPixc.computeSeparateEntities()
        objPixc.computeObjInsideTile()

    #    # Test writeMetadataFile
    #    filename = "../output/" + layer_name + ".xml"
    #    objet.writeMetadataFile(filename)

        # Test computeLakeProducts
        liste_label = [1]
    #    objet.computeLakeProducts(liste_label)

        # Test computeProduct
    #    computeProduct(self, IN_lake_id, IN_indices, IN_classif_dict, IN_size, IN_meanHeight, IN_imp_lon, IN_imp_lat)

        # Test sortPixelsWrtClassifFlags
    #    sortPixelsWrtClassifFlags(In_ind)

        # Test computeCtTime
    #    computeCtTime(self, IN_point):


        # TODO test output
        # TODO test error
        # TODO test other methods
        # TODO test other configuration ? for instance productType = "SP"

    def test_selectWaterDarkLayoverPixels(self):
        """
            unitary test for selectWaterDarkLayoverPixels function
        """

        # Prepare data
        classif_dict = {'layover': [10, 11, 12, 13, 14], 'water': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 'dark': [15, 16, 17, 18, 19]}
        flag_water = True
        flag_dark = False
        flag_layover = False
        # Call function
        out_array = proc_lake.selectWaterDarkLayoverPixels(classif_dict, flag_water, flag_dark, flag_layover)
        self.assertEqual(out_array, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

        # Prepare data
        flag_water = False
        flag_dark = False
        flag_layover = True
        # Call function
        out_array = proc_lake.selectWaterDarkLayoverPixels(classif_dict, flag_water, flag_dark, flag_layover)
        self.assertEqual(out_array, [10, 11, 12, 13, 14])

        # Prepare data
        flag_water = False
        flag_dark = True
        flag_layover = False
        # Call function
        out_array = proc_lake.selectWaterDarkLayoverPixels(classif_dict, flag_water, flag_dark, flag_layover)
        self.assertEqual(out_array, [15, 16, 17, 18, 19])

