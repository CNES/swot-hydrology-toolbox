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

import json
import logging
import numpy as np
from scipy.spatial import KDTree
from osgeo import ogr
import sqlite3

import cnes.common.lib_lake.locnes_variables as my_var


class LakeDb_shp(object):
    
    def __init__(self, in_filename, in_poly=None):
        """
        Constructor
        
        :param in_filename: full path of the lake a priori database
        :type in_filename: string
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon

        Variables of the object:
        filename / string: full path of the lake a priori database
        dataSource / osgeo.ogr.DataSource: reader of lake database
        layer / osgeo.ogr.Layer: layer of the lake database
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Lake DB = %s", in_filename)
        
        # Init with values
        self.filename = in_filename  # Full path of the lake a priori database
        
        # Open database
        self.open_db(in_poly)

    # ----------------------------------------
        
    def open_db(self, in_poly=None):
        """
        Open database, optionnally spatially select polygons and copy layer to memory
        
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Open shapefile in read-only access
        shpDriver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Shapefile driver
        shpDataSource = shpDriver.Open(self.filename, 0)
        
        # 2 - Get the layer
        layer = shpDataSource.GetLayer()
        logger.info("%d lakes stored in database", layer.GetFeatureCount())
        
        # 3 - Select some lakes among BD using in_poly
        if in_poly is not None:
            layer.SetSpatialFilter(in_poly)
            logger.info("%d lakes after focus over studied area", layer.GetFeatureCount())
        
        # 4 - Create an output datasource in memory
        memDriver = ogr.GetDriverByName('MEMORY')  # Memory driver
        self.dataSource = memDriver.CreateDataSource('memData')

        # 5 - Open the memory datasource with write access
        tmp = memDriver.Open('memData', 1)

        # 6 - Copy the layer to memory
        db_mem = self.dataSource.CopyLayer(layer, 'lake_db')
        
        # 7 - Get memory layer
        self.layer = self.dataSource.GetLayer('lake_db')

        # 8 - Close shapefile
        shpDataSource.Destroy()
        
    def close_db(self):
        """
        Close database
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")        
        self.dataSource.Destroy()

    # ----------------------------------------
    
    def getRefValues(self, in_id):
        """
        Getter of reference height and area given the lake identifier
        
        :return out_ref_height: reference height
        :rtype out_ref_height: float
        :return out_ref_area: reference area
        :rtype out_ref_area: float
        """
        
        # Select feature given its id
        self.layer.SetAttributeFilter("%s = '%s'" % (my_var.LAKE_DB_ID, in_id))
        prior_lake = self.layer.GetNextFeature()
        
        # Get reference height
        try:
            out_ref_height = prior_lake.GetField(str("ref_height"))
        except:
            out_ref_height = None
        # Get reference area
        try:
            out_ref_area = prior_lake.GetField(str("ref_area"))
        except:
            out_ref_area = None
        
        # Release filter
        self.layer.SetAttributeFilter(None)
        
        # Output
        return out_ref_height, out_ref_area

    # ----------------------------------------
    def linkToDb(self, in_poly, in_lon, in_lat):
        """
        Links polygon in_poly to a priori database, i.e. returns, when available, 
        the list of the ID(s) of the a priori lake(s) intersecting the input polygon,
        and the reference ID array (corresponding one-to-one with the L2_HR_PIXC) 
        by giving the ID of the closest priori lake
        
        If in_poly corresponds to no a priori lake: return None, None
        
        If in_poly corresponds to 1 a priori lake: return prior_id, [prior_id * size_in_lon]
        
        If in_poly corresponds to 2 or more a priori lakes:
            -out_prior_id contains all prior ID sorted by intersection areas
            -out_pixc_vec_tag returns the prior ID of the closest priori lake for each pixels cloud point
        
        :param in_poly: polygon delineating a water body
        :type in_poly: OGRPolygon
        :param in_lon: improved longitude of PixC related to in_poly
        :type in_lon: 1D array of float
        :param in_lat: improved latitude of PixC related to in_poly
        :type in_lat: 1D array of float

        :return: out_prior_id = list of lake identifiers from the a priori database that intersect in_poly
        :rtype: out_prior_id = list of string
        :return: out_pixc_vec_tag = lake identifiers from the a priori database for each point of the PixC corresponding one-to-one with (in_lon, in_lat)
        :rtype: out_pixc_vec_tag = list of string
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")

        # 1 - Spatial filter of the lake DB over the area covered by the studied polygon
        self.layer.SetSpatialFilter(in_poly)
        
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

                out_pixc_vec_tag = np.empty(in_lon.shape, dtype=object)
                out_pixc_vec_tag[:] = str(curCode)
                return str(curCode), out_pixc_vec_tag
            
        else:  # Many matches: polygon matches 2 or more a priori lakes

            # 2.1 - Init

            set_prior_id = set()  # Set of prior id
            out_pixc_vec_tag = np.empty(in_lon.shape, dtype=object)  # Array of prior id
            prior_geom_coords = []  # List of point coordinates from prior database
            area_intersection = []  # List of area intersection between IN poly and each polygon of prior database

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
                intersection = in_poly.Intersection(curGeom)
                if intersection is not None:

                    area_intersection.append(intersection.GetArea())  # Add the intersection area
                    set_prior_id.add(str(curId))  # Add prior ID to set_prior_id

                    # Get the coordinates with json library
                    json_geom = curGeom.Clone().ExportToJson()
                    multi_polygon = json.loads(json_geom)['coordinates']

                    multi_polygon_array = np.array([])
                    # Loop on multi_polygon in case of multi polygon geometry
                    for poly in multi_polygon:
                        multi_polygon_array = np.append(multi_polygon_array, poly[:-1])

                    # Reshape the 1-D array into 2-D array
                    multi_polygon_array = multi_polygon_array.reshape(int(multi_polygon_array.shape[0]/2), 2)
                    prior_geom_coords.append(multi_polygon_array) # Add the coordinates of the curGeom to the list

            self.layer.SetSpatialFilter(None)  # Delete spatial filter

            # If no true intersection
            if len(set_prior_id) == 0:
                return None, None

            # 2.3 - Put output in good format

            # Compute closest prior lake of each point of pixel cloud
            out_prior_id = list(set_prior_id)
            out_pixc_vec_tag = computeClosestPolygonWithKDTree(in_lon, in_lat, prior_geom_coords, out_prior_id)

            # Sort out_prior_id by decreasing area intersection
            sorted_idx = sorted(range(len(area_intersection)), key=lambda k: area_intersection[k], reverse=True)
            out_prior_id = [out_prior_id[idx] for idx in sorted_idx]

            return out_prior_id, out_pixc_vec_tag
        

#######################################


class LakeDb_sqlite(LakeDb_shp):
    
    def __init__(self, in_filename, in_poly=None):
        """
        Constructor
        
        :param in_filename: full path of the lake a priori database
        :type in_filename: string
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon

        Variables of the object:
        filename / string: full path of the lake a priori database
        driver / string: driver to access to file
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Lake DB = %s", in_filename)
        
        # Init with values
        self.filename = in_filename  # Full path of the lake a priori database
        self.lake_db = None # store the lake database in SQLite format
        self.dataSource = None # store the memory layer (shp file) from the lake database

        # Open database
        self.open_db(in_poly)

    # ----------------------------------------

    def open_db(self, IN_poly=None):
        """
        Open database, optionnally spatially select polygons and copy layer to memory
        
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # #######################################################################################
        # # Open the SQLite database with ogr
        # inDriver = ogr.GetDriverByName('SQLite')
        #
        # inDB = inDriver.Open(self.filename, 0)
        #
        # # Create the output ogr layer (in memory)
        # memDriver = ogr.GetDriverByName('MEMORY')  # Memory driver
        # self.dataSource = memDriver.CreateDataSource('memData')
        #
        # # Open the memory datasource with write access
        # tmp = memDriver.Open('memData', 1)
        #
        # # Copy the SQLite layer to memory
        # pipes_mem = self.dataSource.CopyLayer(inDB.GetLayer('lake'), 'lake', ['OVERWRITE=YES'])
        #
        # self.layer = self.dataSource.GetLayer('lake')
        # #######################################################################################

        # Transform 3D geometry into 2D geometry (necessary for spatialite query)
        IN_poly.FlattenTo2D()

        # Open the SQLite database and define the connector
        self.db_conn = sqlite3.connect(self.filename, timeout=10)

        # Load spatialite extension
        self.db_conn.enable_load_extension(True)
        self.db_conn.execute('SELECT load_extension("mod_spatialite")')

        # Define the cursor
        self.db_cur = self.db_conn.cursor()

        """
        if my_api.GEN_PRINT_LEVEL == 'DEBUG':
            (lakes_nb,) = self.db_cur.execute('SELECT count(*) from lake').fetchone()
            my_api.printDebug("[LakeDb_sqlite] {} lakes stored in database".format(lakes_nb))
        """

        # Create an output datasource in memory
        memDriver = ogr.GetDriverByName('MEMORY')  # Memory driver

        # Open the memory datasource with write access
        self.dataSource = memDriver.CreateDataSource('memData')

        # Set spatial projection
        srs = ogr.osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        # Creating memory layer
        lyr = self.dataSource.CreateLayer(str('lake_db'), srs=srs, geom_type=ogr.wkbPolygon)
        lyr.CreateField(ogr.FieldDefn(my_var.LAKE_DB_ID, ogr.OFTString))

        # Define the layer
        lyr_defn = lyr.GetLayerDefn()

        if IN_poly is not None:
            
            self.db_cur.execute("SELECT {}, AsText(geometry) \
                                 FROM lake \
                                 WHERE MBRIntersects(GeomFromText('{}'), lake.geometry);".format(my_var.LAKE_DB_ID, IN_poly.ExportToWkt()))

            for row in self.db_cur:
                
                # Create empty feature/entity
                out_feat = ogr.Feature(lyr_defn)

                # Fill feature with ID and geometry from SQLite request
                out_feat.SetField(my_var.LAKE_DB_ID, str(row[0]))
                multi_poly = ogr.CreateGeometryFromWkt(row[1])
                out_feat.SetGeometry(multi_poly)

                lyr.CreateFeature(out_feat)

                out_feat.Destroy()

            # Get memory layer
            self.layer = self.dataSource.GetLayer('lake_db')

            logger.info("%d lakes after focus over studied area" % self.layer.GetFeatureCount())

        # Close spatialite database
        self.db_conn.close()


#######################################


class CollectCoordinates:

    def __init__(self, multi_polygon):
        """
        Constructor.

        :param multi_polygon: list of lists containing polygons
        :type multi_polygon: list
        :param list_coords: numpy 2D-array of lon/lat coordinates of the multi polygon

        """
        
        self.list_coords = np.array([])
        # self.multi_polygon = multi_polygon

        # Transform the list of lists into 1D array
        self.multipolygon2coordinates(multi_polygon)

        # Reshape the 1D array into 2D array (lon,lat)
        self.list_coords = self.list_coords.reshape(int(self.list_coords.shape[0] / 2), 2)

    def multipolygon2coordinates(self, list_polygon):
        """
        Recursive function to transform a multipolygon into a list of coordinates.

        :param list_polygon: list containing multi polygons, polygons or coordinates.
        :type list_polygon: list
        """

        # Index of the first comma
        idx = str(list_polygon).index(',')
        # Number of list imbrication
        nb_level = str(list_polygon)[:idx].count('[') - str(list_polygon)[:idx].count(']')

        if nb_level > 3:
            # Call recursively the function
            for elmt in list_polygon:
                self.multipolygon2coordinates(elmt)
        else:
            # Extract [lon, lat] coordinates and store them into self.list_coords variable
            for xy in list_polygon:
                self.list_coords = np.append(self.list_coords, xy[:-1])


#######################################
                

def computeClosestPolygonWithKDTree(IN_lon, IN_lat, prior_geom_coords, prior_id):
    """
    Associate to each pixc_vec coordinate (IN_lon, IN_lat) the closest prior lake and its id

    :param IN_lon: improved longitude of PixC related to IN_poly
    :type IN_lon: 1D array of float
    :param IN_lat: improved latitude of PixC related to IN_poly
    :type IN_lat: 1D array of float
    :param prior_geom_coords: List of coordinates of polygons of prior lake database selected
    :type prior_geom_coords: 2D array of float
    :param prior_id: List of prior ID from lake DB
    :type prior_id: 1D array of str

    :return: list of the closest prior_id associated to the (IN_lon, IN_lat) points
    :rtype: list of str
    """
    
    lon_coords, lat_coords, prior_id_list = np.array([]), np.array([]), np.array([], dtype=object)
    for coords, id in zip(prior_geom_coords, prior_id):
        lon_coords = np.append(lon_coords, coords[:,0])
        lat_coords = np.append(lat_coords, coords[:,1])
        prior_id_list = np.append(prior_id_list, len(coords[:,0])*[id]) # Fill the associated prior_id list of the lon/lat coordinates

    # Build the K-d tree
    tree = KDTree(list(zip(lon_coords, lat_coords)))

    # Built the list point for the query
    points_list = np.vstack((IN_lon, IN_lat)).T

    # Apply K-d tree and get the result: distance and index
    _, kd_tree_idx = tree.query(points_list)

    return prior_id_list[kd_tree_idx]


#######################################


if __name__ == '__main__':
    
    # Create ring
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(1.9, 45.1)
    ring.AddPoint(2.1, 45.1)
    ring.AddPoint(2.1, 45.4)
    ring.AddPoint(1.9, 45.4)
    ring.AddPoint(1.9, 45.1)

    # Create polygon
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    
    lake_db = LakeDb_shp("C:\\Users\\pottierc\\Documents\\no_save\\data\\BD_lacs\\20181002_EU\\apriori_db_lakes_EU.shp", poly)
    print(lake_db.layer.GetFeatureCount())
    print()
    for lake in lake_db.layer:
        print(lake.GetField(str("lake_id")))
    print()
    print(lake_db.getRefValues("44008001651"))
    lake_db.close_db()
