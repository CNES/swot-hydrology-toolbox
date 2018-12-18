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

import cnes.common.lib.my_netcdf_file as my_netcdf_file
import numpy

class TestMyNetcdf(unittest.TestCase):
    """
        class for unitary test of my_netcdf_file module
    """
    def test_myNcReader(self):
        """
            unitary test for myNcReader
        """
        # prepare input data
        fichier = "../data/SWOT_L2_HR_PIXC_000_366_43N-R_20140114T011055_20140114T011056_Dx0000_01.nc"

        # load netCDF file
        netcdf = my_netcdf_file.myNcReader(fichier)

        # dimension test
        groupes = netcdf.content.groups
        dimension = netcdf.getListDim(IN_group=groupes['pixel_cloud'])
        print(dimension)
        netcdf.printListDim(IN_group=groupes['pixel_cloud'])
        nom_dim = 'record'
        dim_value = netcdf.getDimValue(nom_dim, IN_group=groupes['pixel_cloud'])
        self.assertEqual(dim_value, 1198)

        # attribute test
        liste_att = netcdf.getListAtt()
        netcdf.printListAtt(IN_group=groupes['pixel_cloud'])
        nom_att = "pass_number"
        att_value = netcdf.getAttValue(nom_att, IN_group=groupes['pixel_cloud'])
        self.assertEqual(att_value, 366)

        # variable test
        netcdf.printListVar(IN_group=groupes['pixel_cloud'])
        nom_var = "latitude"
        var = netcdf.getVarValue(nom_var, IN_group=groupes['pixel_cloud'])
        print(var)
        unite = netcdf.getVarUnit(nom_var, IN_group=groupes['pixel_cloud'])
        self.assertEqual(unite, "degrees_north")

        netcdf.close()


    def test_myNcWriter(self):
        """
            unitary test for myNcReader
        """
        # prepare input data
        fichier = "out.nc"

        # load netCDF file
        netcdf = my_netcdf_file.myNcWriter(fichier)
        # add group
        groupe = netcdf.add_group("toto")
        # add dimension
        netcdf.add_dimension('record', 10, IN_group=groupe)
        # add global attribute
        netcdf.add_global_attribute('batteur', 'Porcaro', IN_group=groupe)
        # add variable
        netcdf.add_variable('africa', int, 'record', IN_fill_value=-2147483647,  IN_group=groupe)
        netcdf.add_variable_attribute('africa', 'style', 'rock', IN_group=groupe)
        numpy.random.seed(10)
        data = numpy.random.randint(100, size=10)
        netcdf.fill_variable('africa', data, IN_group=groupe)
        # close
        netcdf.close()

        # reopen and test
        netcdf = my_netcdf_file.myNcReader(fichier)
        print(dir(netcdf.content))
        groupes = netcdf.content.groups
        # test variable
        var = netcdf.getVarValue('africa', IN_group=groupes['toto'])
        numpy.testing.assert_array_equal(var, data)
        # test variable attribute
        att_value = groupes['toto'].variables['africa'].style
        self.assertEqual(att_value, 'rock')
        # test global attribute
        att_value = netcdf.getAttValue('batteur', IN_group=groupes['toto'])
        self.assertEqual(att_value, 'Porcaro')
        # close
        netcdf.close()


        



