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

from osgeo import ogr, osr
import cnes.common.lib_lake.locnes_variables as locnes_variables
import cnes.common.lib.my_basins as my_basins

class TestMyBassins(unittest.TestCase):
  """
      class for unitary test of my_bassins module
  """
  def test_linkPolyToContinent(self):
    """
        unitary test for linkPolyToContinent
    """

    # Prepare name of continent file
    locnes_variables.CONTINENT_FILE = os.environ['REF_PATH'] + "/../BD/major_basins/FAO/major_hydrobasins.shp"
    # prepare coordinates of tile
    ring = ogr.Geometry(ogr.wkbLinearRing)
    inner_first_lon = -0.441402512973483
    inner_first_lat = 43.5240335977192
    outer_first_lon = -0.45041238828001945
    outer_first_lat = 43.53357124344608
    outer_last_lon = -0.4385241863532871
    outer_last_lat = 43.52509800108395
    inner_last_lon = -0.4537197700431493
    inner_last_lat = 43.531366731706214
    # add point for polygon
    ring.AddPoint(inner_first_lon, inner_first_lat)
    ring.AddPoint(outer_first_lon, outer_first_lat)
    ring.AddPoint(outer_last_lon, outer_last_lat)
    ring.AddPoint(inner_last_lon, inner_last_lat)
    ring.AddPoint(inner_first_lon, inner_first_lat)
    tile_poly = ogr.Geometry(ogr.wkbPolygon)
    tile_poly.AddGeometry(ring)
    # Launch test
    continent = my_basins.linkPolyToContinent(tile_poly)
    # assert 
    self.assertEqual(continent,"EU")

    # prepare coordinates of tile
    ring = ogr.Geometry(ogr.wkbLinearRing)
    inner_first_lon = -100.441402512973483
    inner_first_lat = 43.5240335977192
    outer_first_lon = -100.45041238828001945
    outer_first_lat = 43.53357124344608
    outer_last_lon = -100.4385241863532871
    outer_last_lat = 43.52509800108395
    inner_last_lon = -100.4537197700431493
    inner_last_lat = 43.531366731706214
    # add point for polygon
    ring.AddPoint(inner_first_lon, inner_first_lat)
    ring.AddPoint(outer_first_lon, outer_first_lat)
    ring.AddPoint(outer_last_lon, outer_last_lat)
    ring.AddPoint(inner_last_lon, inner_last_lat)
    ring.AddPoint(inner_first_lon, inner_first_lat)
    tile_poly = ogr.Geometry(ogr.wkbPolygon)
    tile_poly.AddGeometry(ring)
    # Launch test
    continent = my_basins.linkPolyToContinent(tile_poly)
    # assert 
    self.assertEqual(continent,"NA")

