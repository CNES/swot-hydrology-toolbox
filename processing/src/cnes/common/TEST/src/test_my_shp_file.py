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

import os
import os.path
import unittest

from osgeo import ogr, osr
import cnes.common.lib.my_shp_file as my_shp_file
import cnes.common.lib.my_tools as my_tools

def prepare_data(inputData):
    """
        This function convert POLYGON into list of value for assert
    """
    outputData = list(map(float, inputData.replace("POLYGON ((", "").replace("))", "").replace(", ", " ").replace(",", " ").split(' ')))
    return outputData


class TestMyVariables(unittest.TestCase):
    """
        class for unitary test of my_shp_file module
    """
    def test_merge_mem_layer_with_shp(self):
        """
            unitary test for merge_mem_layer_with_shp
        """
         # TODO check ordre des points des polygones d'entrée, peut-être une erreur ???
        # prepare input data
        # create spatial reference
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        # create first shapefile
        driver = ogr.GetDriverByName(str('ESRI Shapefile'))
        out_shapefile = "layer_1.shp"
        if os.path.exists(out_shapefile):
            driver.DeleteDataSource(out_shapefile)
        data_source1 = driver.CreateDataSource(out_shapefile)
        layer = data_source1.CreateLayer("polygone", srs, ogr.wkbPolygon)
        id_field = ogr.FieldDefn("id", ogr.OFTInteger)
        layer.CreateField(id_field)
        feature = ogr.Feature(layer.GetLayerDefn())
        polygone = ogr.Geometry(ogr.wkbPolygon)
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(0.0, 40.0)
        ring.AddPoint(1.0, 40.0)
        ring.AddPoint(1.0, 41.0)
        ring.AddPoint(0.0, 41.0)
        ring.AddPoint(0.0, 40.0)
        polygone.AddGeometry(ring)
        feature.SetGeometry(polygone)
        feature.SetField("id", 1)
        layer.CreateFeature(feature)

        data_source1.CopyLayer(layer, "polygone")

        # create second shapefile
        out_shapefile2 = "layer_2.shp"
        if os.path.exists(out_shapefile2):
            driver.DeleteDataSource(out_shapefile2)
        data_source2 = driver.CreateDataSource(out_shapefile2)
        layer = data_source2.CreateLayer("polygone", srs, ogr.wkbPolygon)
        id_field = ogr.FieldDefn("id", ogr.OFTInteger)
        layer.CreateField(id_field)
        feature = ogr.Feature(layer.GetLayerDefn())
        polygone = ogr.Geometry(ogr.wkbPolygon)
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(1.0, 40.0)
        ring.AddPoint(2.0, 40.0)
        ring.AddPoint(2.0, 41.0)
        ring.AddPoint(1.0, 41.0)
        ring.AddPoint(1.0, 40.0)
        polygone.AddGeometry(ring)
        feature.SetGeometry(polygone)
        feature.SetField("id", 1)
        layer.CreateFeature(feature)

        data_source2.CopyLayer(layer, "polygone")

        # create memory layer
        mem_driver = ogr.GetDriverByName(str('MEMORY'))
        data_source = mem_driver.CreateDataSource('memData')
        layer = data_source.CreateLayer("polygone", srs, ogr.wkbPolygon)
        id_field = ogr.FieldDefn("id", ogr.OFTInteger)
        layer.CreateField(id_field)

        liste = ["layer_1.shp", "layer_2.shp"]
        # call function
        datasource, layer = my_shp_file.merge_mem_layer_with_shp(liste, layer)
        # ref data
        ref_geom = []
        ref_geom.append("POLYGON ((0 40,0 41,1 41,1 40,0 40))")
        ref_geom.append("POLYGON ((1 40,1 41,2 41,2 40,1 40))")

        # assert
        indice = 0
        for feature in layer:
            geom = feature.GetGeometryRef()
            self.assertEqual(str(geom), ref_geom[indice])
            print(geom)
            print(ref_geom[indice])
            indice = indice + 1

    def test_merge_2_layers(self):
        """
            unitary test for merge_2_layers
        """
        # prepare input data
        # create spatial reference
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        # create first layer
        mem_driver = ogr.GetDriverByName(str('MEMORY'))
        data_source1 = mem_driver.CreateDataSource('memData')
        layer1 = data_source1.CreateLayer("polygone", srs, ogr.wkbPolygon)
        id_field = ogr.FieldDefn("id", ogr.OFTInteger)
        layer1.CreateField(id_field)
        feature = ogr.Feature(layer1.GetLayerDefn())
        polygone = ogr.Geometry(ogr.wkbPolygon)
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(0.0, 40.0)
        ring.AddPoint(1.0, 40.0)
        ring.AddPoint(1.0, 41.0)
        ring.AddPoint(0.0, 41.0)
        ring.AddPoint(0.0, 40.0)
        polygone.AddGeometry(ring)
        feature.SetGeometry(polygone)
        feature.SetField("id", 1)
        layer1.CreateFeature(feature)

        # create second layer
        data_source2 = mem_driver.CreateDataSource('memData')
        layer2 = data_source2.CreateLayer("polygone", srs, ogr.wkbPolygon)
        id_field = ogr.FieldDefn("id", ogr.OFTInteger)
        layer2.CreateField(id_field)
        feature = ogr.Feature(layer2.GetLayerDefn())
        polygone = ogr.Geometry(ogr.wkbPolygon)
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(1.0, 40.0)
        ring.AddPoint(2.0, 40.0)
        ring.AddPoint(2.0, 41.0)
        ring.AddPoint(1.0, 41.0)
        ring.AddPoint(1.0, 40.0)
        polygone.AddGeometry(ring)
        feature.SetGeometry(polygone)
        feature.SetField("id", 1)
        layer2.CreateFeature(feature)

        # call function
        data_source, layer = my_shp_file.merge_2_layers(layer1, layer2)

        # reference data
        ref_geom = []
        # FIXME pourquoi l'inversion des points ici ??
        ref_geom.append("POLYGON ((1 40 0,0 40 0,0 41 0,1 41 0,1 40 0))")
        ref_geom.append("POLYGON ((1 40 0,2 40 0,2 41 0,1 41 0,1 40 0))")

        layer = data_source.GetLayer()
        # assert
        indice = 0
        for feature in layer:
            geom = feature.GetGeometryRef()
            self.assertEqual(str(geom), ref_geom[indice])
            print(geom)
            print(ref_geom[indice])
            indice = indice + 1

    def test_write_mem_layer_as_shp(self):
        """
            unitary test for write_mem_layer_as_shp
        """
        # prepare input data
        # create spatial reference
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        # create layer
        mem_driver = ogr.GetDriverByName(str('MEMORY'))
        data_source = mem_driver.CreateDataSource('memData')
        layer = data_source.CreateLayer("polygone", srs, ogr.wkbPolygon)
        id_field = ogr.FieldDefn("id", ogr.OFTInteger)
        layer.CreateField(id_field)
        feature = ogr.Feature(layer.GetLayerDefn())
        polygone = ogr.Geometry(ogr.wkbPolygon)
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(0.0, 40.0)
        ring.AddPoint(1.0, 40.0)
        ring.AddPoint(1.0, 41.0)
        ring.AddPoint(0.0, 41.0)
        ring.AddPoint(0.0, 40.0)
        polygone.AddGeometry(ring)
        feature.SetGeometry(polygone)
        feature.SetField("id", 1)
        layer.CreateFeature(feature)
        file_name = "layer_test_writeMemLayer_asShp.shp"

        # call function
        my_shp_file.write_mem_layer_as_shp(layer, file_name)

        # assert
        my_tools.testFile(file_name)


