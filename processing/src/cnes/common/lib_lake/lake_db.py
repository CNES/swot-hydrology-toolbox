# -*- coding: utf8 -*-
"""
.. module:: lake_db.py
    :synopsis: Deal with operationnal lake a priori database; consider shapefile and sqlite formats
    Created on 08/27/2018

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""

import numpy as np
from osgeo import ogr
from shapely.geometry import Point
from shapely.wkt import loads
import logging

import cnes.common.lib.my_api as my_api
import cnes.common.lib_lake.locnes_variables as my_var


class LakeDb_shp(object):
    
    def __init__(self, IN_filename, IN_poly=None):
        """
        Constructor
        
        :param IN_filename: full path of the lake a priori database
        :type IN_filename: string
        :param IN_poly: polygon to spatially select lakes from DB
        :type IN_poly: ogr.Polygon

        Variables of the object:
        filename / string: full path of the lake a priori database
        dataSource / osgeo.ogr.DataSource: reader of lake database
        layer / osgeo.ogr.Layer: layer of the lake database
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("[LakeDb_shp] == INIT ==")
        logger.info("[LakeDb_shp] Lake DB = %s" % IN_filename)
        
        # Init with values
        self.filename = IN_filename  # Full path of the lake a priori database
        
        # Open database
        self.open_db(IN_poly)

    # ----------------------------------------
        
    def open_db(self, IN_poly=None):
        """
        Open database, optionnally spatially select polygons and copy layer to memory
        
        :param IN_poly: polygon to spatially select lakes from DB
        :type IN_poly: ogr.Polygon
        """
        logger = logging.getLogger(self.__class__.__name__)
        # 1 - Open shapefile in read-only access
        shpDriver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Shapefile driver
        shpDataSource = shpDriver.Open(self.filename, 0)
        
        # 2 - Get the layer
        layer = shpDataSource.GetLayer()
        logger.info("[LakeDb_shp] %d lakes stored in database" % layer.GetFeatureCount())
        
        # 3 - Select some lakes among BD using IN_poly
        if IN_poly is not None:
            layer.SetSpatialFilter(IN_poly)
            logger.info("[LakeDb_shp] %d lakes after focus over studied area" % layer.GetFeatureCount())
        
        # 4 - Create an output datasource in memory
        memDriver = ogr.GetDriverByName('MEMORY')  # Memory driver
        self.dataSource = memDriver.CreateDataSource('memData')

        # 5 - Open the memory datasource with write access
        tmp = memDriver.Open('memData', 1)

        # 6 - Copy the layer to memory
        db_mem = self.dataSource.CopyLayer(layer,'lake_db')
        
        # 7 - Get memory layer
        self.layer = self.dataSource.GetLayer('lake_db')
        
        # 8 - Close shapefile
        shpDataSource.Destroy()
        
    def close_db(self):
        """
        Close database
        """
        self.dataSource.Destroy()

    # ----------------------------------------
    
    def linkToDb(self, IN_poly, IN_lon, IN_lat):
        """
        Links polygon IN_poly to a priori database, i.e. returns, when available, 
        the list of the ID(s) of the a priori lake(s) intersecting the input polygon,
        and the reference ID array (corresponding one-to-one with the L2_HR_PIXC) 
        by giving the ID of the closest priori lake
        
        If IN_poly corresponds to no a priori lake: return None, None
        
        If IN_poly corresponds to 1 a priori lake: return prior_id, [prior_id * size_IN_lon]
        
        If IN_poly corresponds to 2 or more a priori lakes:
            -OUT_prior_id contains all prior ID sorted by intersection areas
            -OUT_pixc_vec_tag returns the prior ID of the closest priori lake for each pixels cloud point
        
        :param IN_poly: polygon delineating a water body
        :type IN_poly: OGRPolygon
        :param IN_lon: improved longitude of PixC related to IN_poly
        :type IN_lon: 1D array of float
        :param IN_lat: improved latitude of PixC related to IN_poly
        :type IN_lat: 1D array of float

        :return: OUT_prior_id = list of lake identifiers from the a priori database that intersect IN_poly
        :rtype: OUT_prior_id = list of string
        :return: OUT_pixc_vec_tag = lake identifiers from the a priori database for each point of the PixC corresponding one-to-one with (IN_lon, IN_lat)
        :rtype: OUT_pixc_vec_tag = list of string
        """
        my_api.printDebug("[lake_db] == linkToDb ==")

        # 1 - Spatial filter of the lake DB over the area covered by the studied polygon
        self.layer.SetSpatialFilter(IN_poly)
        
        # 2 - Processing according to the number of a priori lakes intersecting polygon
        nb_lakes = self.layer.GetFeatureCount()
        
        if nb_lakes == 0:  # Polygon matches no a priori lake
            self.layer.SetSpatialFilter(None)  # Delete spatial filter
            return None, None
        
        elif nb_lakes == 1:  # Easy match: polygon matches only one a priori lake
            for curLake_bd in self.layer:
                # Get the a priori identifier
                curCode = curLake_bd.GetField(my_var.LAKE_DB_ID)  
                # Test but should not occur...
                if curCode is None:
                    return None, None
                # Return output variables
                self.layer.SetSpatialFilter(None)  # Delete spatial filter
                # Compute PIXCVec_tag
                OUT_pixc_vec_tag = np.empty(IN_lon.shape, dtype=object)
                OUT_pixc_vec_tag[:] = str(curCode)
                return curCode, OUT_pixc_vec_tag
            
        else:  # Many matches: polygon matches 2 or more a priori lakes
            
            # 2.1 - Init
            OUT_prior_id = []  # List of prior id
            OUT_pixc_vec_tag = np.empty(IN_lon.shape, dtype=object)  # Array of prior id
            prior_geom_list = []  # List of prior geometries
            area_intersection = []  # List of area intersection between IN poly and each polygon of prior database (prior_geom_list)
            
            # 2.2 - List area of intersections with a priori geometries and id of these geometries
            for curLake_bd in self.layer:
                # Get the a priori identifier
                curId = curLake_bd.GetField(my_var.LAKE_DB_ID)  
                # Test but should not occur...
                if curId is None:
                    return None, None
                # Get geometry
                curGeom = curLake_bd.GetGeometryRef()
                # Compute exact area of intersection
                intersection = IN_poly.Intersection(curGeom)
                if intersection is not None:
                    prior_geom_list.append(loads(curGeom.Clone().ExportToWkt()))  # Add prior lake in shapely format to list prior_geom_list
                    area_intersection.append(intersection.GetArea())  # Add 
                    OUT_prior_id.append(str(curId))  # Add prior ID to list OUT_prior_id
                    
            self.layer.SetSpatialFilter(None)  # Delete spatial filter
            
            # If no true intersection
            if len(OUT_prior_id) == 0:
                return None, None
            
            # 2.3 - Put output in good format
            
            # Compute closest prior lake of each point of pixel cloud
            OUT_pixc_vec_tag = [OUT_prior_id[computeIdxOfClosestPolygonFromPoint(Point(IN_lon[idx], IN_lat[idx]), prior_geom_list)] for idx in range(len(IN_lon))]
            # Sort OUT_prior_id by decreasing area intersection
            sorted_idx = sorted(range(len(area_intersection)), key=lambda k: area_intersection[k], reverse=True)
            OUT_prior_id = [OUT_prior_id[idx] for idx in sorted_idx]

            return ';'.join(OUT_prior_id), OUT_pixc_vec_tag
        

#######################################


class LakeDb_sqlite(object):
    
    def __init__(self, IN_filename, IN_poly=None):
        """
        Constructor
        
        :param IN_filename: full path of the lake a priori database
        :type IN_filename: string
        :param IN_poly: polygon to spatially select lakes from DB
        :type IN_poly: ogr.Polygon

        Variables of the object:
        filename / string: full path of the lake a priori database
        driver / string: driver to access to file
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("[LakeDb_sqlite] == INIT ==")
        logger.info("[LakeDb_sqlite] Lake DB = %s" % IN_filename)
        
        # Init with values
        self.filename = IN_filename  # Full path of the lake a priori database
        
        # Open database
        self.open_db(IN_poly)

    # ----------------------------------------
        
    def open_db(self, IN_poly=None):
        """
        Open database, optionnally spatially select polygons and copy layer to memory
        
        :param IN_poly: polygon to spatially select lakes from DB
        :type IN_poly: ogr.Polygon
        """
        pass
    
    def close_db(self):
        """
        Close database
        """
        pass
        

#######################################


def computeIdxOfClosestPolygonFromPoint(IN_point, IN_geom_list):
    """
    Get index in list IN_geom_list of the polygon which is the closest to IN_point
    
    :param IN_point: Point of PIXCVec
    :type IN_point: Shapely point
    :param IN_geom_list: List of polygons of prior lake database
    :type IN_geom_list: list of shapely polygones
    
    :return: indice of IN_geom_list
    :rtype: int
    """
    distance_vector = []

    for geom_idx, geom in enumerate(IN_geom_list):

        # If the point is contained inside of the polygon
        if geom.contains(IN_point):
            return geom_idx
        else:
            # Compute distance between polygon and point
            distance_vector.append(geom.exterior.distance(IN_point))

    # Sort by increasing distance and return the first element
    return np.argsort(distance_vector)[0]


#######################################


if __name__ == '__main__':
    
    # Create ring
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(1.2, 46.7)
    ring.AddPoint(1.3, 46.7)
    ring.AddPoint(1.3, 46.8)
    ring.AddPoint(1.2, 46.8)
    ring.AddPoint(1.2, 46.7)

    # Create polygon
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    
    lake_db = LakeDb_shp("C:\\Users\\pottierc\\Documents\\01_SWOT\\ADT\\WG10_HydroProducts\\Support\\2016-2018\\CSF_20180430\\Lot2_BDs_apriori\\Tache_2.1_BD_lacs\\2.1.4_Implementation_BD_lacs\\inputs\\priordb_lakes_france.shp", poly)
    lake_db.layer.SetSpatialFilter(None)
    lake_db.close_db()
