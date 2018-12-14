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

import unittest
import numpy

import cnes.common.serviceError as serviceError
import cnes.common.lib_lake.locnes_variables as locnes_variables
import cnes.common.lib.my_hull as my_hull

def prepare_data(input_data):
    """
        This function convert POLYGON into list of value for assert
    """
    output_data = list(map(float, input_data.replace("POLYGON ((", "").replace("))", "").replace(", ", " ").replace(",", " ").split(' ')))
    return output_data


class TestMyHull(unittest.TestCase):
    """
        class for unitary test of my_hull module
    """
    def test_compute_lake_boundaries(self):
        """
            unitary test for compute_lake_boundaries
        """

        # Prepare input data
        nb_pix_range = 3117
        v_long = numpy.array(list(map(float, open("../data/imp_lon", "r").read().replace(" \n", "").split(' '))))
        v_lat = numpy.array(list(map(float, open("../data/imp_lat", "r").read().replace(" \n", "").split(' '))))
        v_range = numpy.array(list(map(int, open("../data/range_idx", "r").read().replace(" \n", "").split(' '))))
        v_azimuth = numpy.array(list(map(int, open("../data/azimuth_idx", "r").read().replace(" \n", "").split(' '))))

        retour = my_hull.compute_lake_boundaries(v_long, v_lat, v_range, v_azimuth, nb_pix_range)

        ref = open("../data/polygon_computeLakeBoundaries", "r").read().replace("\n", "")
        self.assertAlmostEqual(prepare_data(str(retour)), prepare_data(ref), places=8)

        # Test bad HULL_METHOD
        locnes_variables.HULL_METHOD = 12
        with self.assertRaises(serviceError.ProcessingError):
            my_hull.compute_lake_boundaries(v_long, v_lat, v_range, v_azimuth, nb_pix_range)

    def test_alpha_shape(self):
        """
            unitary test for alpha_shape
        """

        # Prepare input data
        nb_pix_range = 3117
        v_long = numpy.array(list(map(float, open("../data/imp_lon", "r").read().replace(" \n", "").split(' '))))
        v_lat = numpy.array(list(map(float, open("../data/imp_lat", "r").read().replace(" \n", "").split(' '))))
        v_range = numpy.array(list(map(int, open("../data/range_idx", "r").read().replace(" \n", "").split(' '))))
        coords = numpy.zeros((v_long.size, 2))
        coords[:, 0] = v_long
        coords[:, 1] = v_lat
        alpha = (2000 + 2000 * v_range / nb_pix_range).astype('int')

        # call function
        retour = my_hull.alpha_shape(coords, alpha)
        # get ref data
        ref = open("../data/polygon_alpha_shape", "r").read().replace("\n", "")
        # assert
        self.assertAlmostEqual(prepare_data(str(retour)), prepare_data(ref), places=8)

    def test_get_concav_hull_bis(self):
        """
            unitary test for get_concav_hull_bis
        """

        # Prepare input data
        v_long = numpy.array(list(map(float, open("../data/imp_lon", "r").read().replace(" \n", "").split(' '))))
        v_lat = numpy.array(list(map(float, open("../data/imp_lat", "r").read().replace(" \n", "").split(' '))))
        coords = numpy.zeros((v_long.size, 2))
        coords[:, 0] = v_long
        coords[:, 1] = v_lat

        # call function
        retour = my_hull.get_concav_hull_bis(coords)
        # get ref data
        ref = open("../data/polygon_concavHullBis", "r").read().replace("\n", "")
        # assert
        self.assertAlmostEqual(prepare_data(str(retour)), prepare_data(ref), places=8)

    def test_get_concave_hull_from_radar_vectorisation(self):
        """
            unitary test for get_concave_hull_from_radar_vectorisation
        """

        # Prepare input data
        v_long = numpy.array(list(map(float, open("../data/imp_lon", "r").read().replace(" \n", "").split(' '))))
        v_lat = numpy.array(list(map(float, open("../data/imp_lat", "r").read().replace(" \n", "").split(' '))))
        v_range = numpy.array(list(map(int, open("../data/range_idx", "r").read().replace(" \n", "").split(' '))))
        v_azimuth = numpy.array(list(map(int, open("../data/azimuth_idx", "r").read().replace(" \n", "").split(' '))))

        # call function
        retour = my_hull.get_concave_hull_from_radar_vectorisation(v_range, v_azimuth, v_long, v_lat)
        # get ref data
        ref = open("../data/polygon_concavHullFromRadarVectorisation", "r").read().replace("\n", "")
        # assert
        self.assertAlmostEqual(prepare_data(str(retour)), prepare_data(ref), places=8)

    def test_get_circum_ratio(self):
        """
            unitary test for get_circum_ratio
        """

        # Prepare input data
        point_a = (1.0, 40.0)
        point_b = (3.0, 40.0)
        point_c = (1.0, 20.0)

        # call function
        retour = my_hull.get_circum_ratio(point_a, point_b, point_c)
        # assert
        self.assertAlmostEqual(retour, 10.04987562, places=8)

    def test_get_max_segment(self):
        """
            unitary test for get_max_segment
        """

        # Prepare input data
        point_a = (1.0, 40.0)
        point_b = (3.0, 40.0)
        point_c = (1.0, 20.0)

        # call function
        retour = my_hull.get_max_segment(point_a, point_b, point_c)
        # assert
        self.assertAlmostEqual(retour, 20.09975124, places=8)

