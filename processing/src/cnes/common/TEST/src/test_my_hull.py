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
import numpy

from osgeo import ogr, osr
import cnes.common.serviceError as serviceError
import cnes.common.lib_lake.locnes_variables as locnes_variables
import cnes.common.lib.my_hull as my_hull

def prepare_data(inputData):
  """
      This function convert POLYGON into list of value for assert
  """
  outputData=list(map(float, inputData.replace("POLYGON ((","").replace("))","").replace(", "," ").replace(","," ").split(' ')))
  return outputData


class TestMyHull(unittest.TestCase):
  """
      class for unitary test of my_hull module
  """
  def test_computeLakeBoundaries(self):
    """
        unitary test for computeLakeBoundaries
    """

    # Prepare input data
    nb_pix_range = 3117
    vLong = numpy.array(list(map(float, open("../data/imp_lon", "r").read().replace(" \n","").split(' '))))
    vLat = numpy.array(list(map(float, open("../data/imp_lat", "r").read().replace(" \n","").split(' '))))
    vRange = numpy.array(list(map(int, open("../data/range_idx", "r").read().replace(" \n","").split(' '))))
    vAzimuth = numpy.array(list(map(int, open("../data/azimuth_idx", "r").read().replace(" \n","").split(' '))))
    
    retour = my_hull.computeLakeBoundaries(vLong, vLat, vRange, vAzimuth, nb_pix_range)
   
    ref = open("../data/polygon_computeLakeBoundaries", "r").read().replace("\n","")
    self.assertAlmostEqual(prepare_data(str(retour)), prepare_data(ref), places=8)

    # Test bad HULL_METHOD
    locnes_variables.HULL_METHOD = 12
    with self.assertRaises(serviceError.ProcessingError):
       objet = my_hull.computeLakeBoundaries(vLong, vLat, vRange, vAzimuth, nb_pix_range)

  def test_alpha_shape(self):
    """
        unitary test for alpha_shape
    """

    # Prepare input data
    nb_pix_range = 3117
    vLong = numpy.array(list(map(float, open("../data/imp_lon", "r").read().replace(" \n","").split(' '))))
    vLat = numpy.array(list(map(float, open("../data/imp_lat", "r").read().replace(" \n","").split(' '))))
    vRange = numpy.array(list(map(int, open("../data/range_idx", "r").read().replace(" \n","").split(' '))))
    coords = numpy.zeros((vLong.size, 2))
    coords[:, 0] = vLong
    coords[:, 1] = vLat
    alpha = (2000 + 2000 * vRange / nb_pix_range).astype('int')

    # call function
    retour = my_hull.alpha_shape(coords, alpha)
    # get ref data
    ref = open("../data/polygon_alpha_shape", "r").read().replace("\n","")
    # assert
    self.assertAlmostEqual(prepare_data(str(retour)), prepare_data(ref), places=8)

  def test_getConcavHull_bis(self):
    """
        unitary test for getConcavHull_bis
    """

    # Prepare input data
    vLong = numpy.array(list(map(float, open("../data/imp_lon", "r").read().replace(" \n","").split(' '))))
    vLat = numpy.array(list(map(float, open("../data/imp_lat", "r").read().replace(" \n","").split(' '))))
    coords = numpy.zeros((vLong.size, 2))
    coords[:, 0] = vLong
    coords[:, 1] = vLat

    # call function
    retour = my_hull.getConcavHull_bis(coords)
    # get ref data
    ref = open("../data/polygon_concavHullBis", "r").read().replace("\n","")
    # assert
    self.assertAlmostEqual(prepare_data(str(retour)), prepare_data(ref), places=8)

  def test_getConcaveHullFromRadarVectorisation(self):
    """
        unitary test for getConcaveHullFromRadarVectorisation
    """

    # Prepare input data
    vLong = numpy.array(list(map(float, open("../data/imp_lon", "r").read().replace(" \n","").split(' '))))
    vLat = numpy.array(list(map(float, open("../data/imp_lat", "r").read().replace(" \n","").split(' '))))
    vRange = numpy.array(list(map(int, open("../data/range_idx", "r").read().replace(" \n","").split(' '))))
    vAzimuth = numpy.array(list(map(int, open("../data/azimuth_idx", "r").read().replace(" \n","").split(' '))))

    # call function
    retour = my_hull.getConcaveHullFromRadarVectorisation(vRange, vAzimuth, vLong, vLat)
    # get ref data
    ref = open("../data/polygon_concavHullFromRadarVectorisation", "r").read().replace("\n","")
    # assert
    self.assertAlmostEqual(prepare_data(str(retour)), prepare_data(ref), places=8)

  def test_getCircumRatio(self):
    """
        unitary test for getCircumRatio
    """

    # Prepare input data
    pa = (1.0, 40.0)
    pb = (3.0, 40.0)
    pc = (1.0, 20.0)

    # call function
    retour = my_hull.getCircumRatio(pa, pb, pc)
    # assert
    self.assertAlmostEqual(retour, 10.04987562, places=8)

  def test_getMaxSegment(self):
    """
        unitary test for getMaxSegment
    """

    # Prepare input data
    pa = (1.0, 40.0)
    pb = (3.0, 40.0)
    pc = (1.0, 20.0)

    # call function
    retour = my_hull.getMaxSegment(pa, pb, pc)
    # assert
    self.assertAlmostEqual(retour, 20.09975124, places=8)

