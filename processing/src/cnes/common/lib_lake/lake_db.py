# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2019/05/17:version initiale.
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: lake_db.py
    :synopsis: Deal with operationnal lake a priori database; consider shapefile and sqlite formats
     Created on 2018/08/27

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""

import json
import logging
import numpy as np
from osgeo import ogr
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist
import sqlite3

import cnes.common.service_config_file as service_config_file
import cnes.common.lib.my_tools as my_tools

class LakeDbShp(object):
    """
        class LakeDbShp
    """
    def __init__(self, in_lake_db_filename, in_influence_lake_db_filename=None, in_poly=None):
        """
        Constructor
        
        :param in_lake_db_filename: full path of the lake a priori database
        :type in_lake_db_filename: string
        :param in_influence_lake_db_filename: full path of the influence lake a priori table
        :type in_influence_lake_db_filename: string
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon


        Variables of the object:
        
         - filename / string: full path of the lake a priori database
         - data_source / osgeo.ogr.DataSource: reader of lake database
         - layer / osgeo.ogr.Layer: layer of the lake database
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Lake DB = %s", in_lake_db_filename)
        # Get config file
        cfg = service_config_file.get_instance()
        # get lake_db_id parameter
        self.lake_db_id = cfg.get("DATABASES", "LAKE_DB_ID")

        # Set influence map flag to True if influence map path is given
        self.influence_map_flag = in_influence_lake_db_filename is not None

        # Init with values
        self.filename = in_lake_db_filename  # Full path of the lake a priori database
        
        # Open database
        self.lake_ds, self.lake_layer = self.open_db(self.filename, in_poly) # Open Lake database
        if self.influence_map_flag : # open influence map database
            logger.info("Influence Lake DB = %s", in_influence_lake_db_filename)
            self.influence_lake_ds, self.influence_lake_layer = self.open_db(in_influence_lake_db_filename, in_poly)

    # ----------------------------------------
        
    def open_db(self, in_file_path, in_poly=None):
        """
        Open database, optionnally spatially select polygons and copy layer to memory

        :param in_file_path in_poly: full path to DB
        :type in_file_path in_poly: string
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Start")

        # 1 - Open shapefile in read-only access
        shp_driver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Shapefile driver
        shp_data_source = shp_driver.Open(in_file_path, 0)
        
        # 2 - Get the layer
        layer = shp_data_source.GetLayer()
        logger.info("%d lakes stored in database", layer.GetFeatureCount())
        
        # 3 - Select some lakes among BD using in_poly
        if in_poly is not None:
            layer.SetSpatialFilter(in_poly)
            logger.info("%d lakes after focus over studied area", layer.GetFeatureCount())
        
        # 4 - Create an output datasource in memory
        mem_driver = ogr.GetDriverByName('MEMORY')  # Memory driver
        data_source = mem_driver.CreateDataSource('memData')

        # 5 - Open the memory datasource with write access
        mem_driver.Open('memData', 1)

        # 6 - Copy the layer to memory
        data_source.CopyLayer(layer, 'lake_db')
        
        # 7 - Get memory layer
        layer = data_source.GetLayer('lake_db')

        # 8 - Close shapefile
        shp_data_source.Destroy()
        logger.info("Stop")
        return data_source, layer
        
    def close_db(self):
        """
        Close database
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- close database -")
        self.lake_ds.Destroy()
        if self.influence_map_flag :
            self.influence_lake_ds.Destroy()

    # ----------------------------------------
    
    def get_prior_values(self, in_id):
        """
        Getter of name, GRanD identifier, reference height and area given the lake identifier
        
        :return: out_name: prior name
        :rtype: string
        :return: out_grand: GRanD identifier
        :rtype: string
        :return: out_ref_height: reference height
        :rtype: float
        :return: out_ref_area: reference area
        :rtype: float
        """
        
        # Select feature given its id
        self.lake_layer.SetAttributeFilter("%s = '%s'" % (self.lake_db_id, in_id))
        prior_lake_feat = self.lake_layer.GetNextFeature()
            
        # Get name
        try:
            out_name = prior_lake_feat.GetField(str("name"))
        except:
            out_name = None
            
        # Get GRanD identifier : Dam ID from Global Reservoir and Dam (GranD) database
        try:
            out_grand = prior_lake_feat.GetField(str("grand_id"))
        except:
            out_grand = None

        # Get reference height
        try:
            out_ref_height = prior_lake_feat.GetField(str("ref_height"))
        except:
            out_ref_height = None
            
        # Get reference area
        try:
            out_ref_area = prior_lake_feat.GetField(str("ref_area"))
        except:
            out_ref_area = None
        
        # Release filter
        self.lake_layer.SetAttributeFilter(None)
        
        # Output
        return out_name, out_grand, out_ref_height, out_ref_area

    # ----------------------------------------
    
    def link_to_db(self, in_poly, in_lon, in_lat):
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
        :rtype: list of string
        :return: out_pixc_vec_tag = lake identifiers from the a priori database for each point of the PixC corresponding one-to-one with (in_lon, in_lat)
        :rtype: list of string
        """
        logger = logging.getLogger(self.__class__.__name__)


        # 1 - Spatial filter of the lake DB over the area covered by the studied polygon
        self.lake_layer.SetSpatialFilter(in_poly)

        # 2 - Processing according to the number of a priori lakes intersecting polygon
        nb_lakes = self.lake_layer.GetFeatureCount()
        logger.debug("Current lake matched with %d lakes from a priori lake database" %(nb_lakes))

        retour_1 = None
        retour_2 = None

        if nb_lakes == 0:  # Polygon matches no a priori lake
            self.lake_layer.SetSpatialFilter(None)  # Delete spatial filter
            retour_1 = None
            retour_2 = None

        elif nb_lakes == 1:  # Easy match: polygon matches only one a priori lake

            cur_lake_bd = self.lake_layer.GetNextFeature()
            cur_id = cur_lake_bd.GetField(self.lake_db_id)

            logger.debug("A priori lakedb_id is : %s" %(cur_id))
            # Test but should not occur...
            if cur_id is None:
                retour_1 = None
                retour_2 = None
            else:
                # Return output variables
                self.lake_layer.SetSpatialFilter(None)  # Delete spatial filter
                # Compute PIXCVec_tag
                out_pixc_vec_tag = np.empty(in_lon.shape, dtype=object)
                out_pixc_vec_tag[:] = str(cur_id)
                retour_1 = str(cur_id)
                retour_2 = out_pixc_vec_tag
            
        else:  # Many matches: polygon matches 2 or more a priori lakes

            # 2.1 - Init
            prior_id = []  # Set of prior id
            area_intersection = []  # List of area intersection between IN poly and each polygon of prior database
            prior_geoms = []

            # 2.2 - List area of intersections with a priori geometries and id of these geometries
            for cur_lake_bd in self.lake_layer:
                # Get the a priori identifier
                cur_id = cur_lake_bd.GetField(self.lake_db_id)

                # Test but should not occur...
                if cur_id is None:
                    retour_1 = None
                    retour_2 = None
                else:
                    # Get geometry
                    cur_geom = cur_lake_bd.GetGeometryRef().Clone()
                    # Compute exact area of intersection
                    intersection = in_poly.Intersection(cur_geom)
                    if intersection is not None:
    
                        area_intersection.append(intersection.GetArea())  # Add the intersection area
                        prior_id.append(str(cur_id))  # Add prior ID to set_prior_id
                        prior_geoms.append(cur_geom)

            self.lake_layer.SetSpatialFilter(None)  # Delete spatial filter

            # 2.3 - Put output in good format
            if len(prior_id) == 0:
                # If no true intersection
                retour_1 = None
                retour_2 = None
            else:
                # Computation time : compute_pixc_vec_tag_with_influence_area_map * 2,5 = compute_closest_polygon_with_KDtree
                if self.influence_map_flag:
                    out_pixc_vec_tag = self.compute_pixc_vec_tag_with_influence_area_map(in_lon, in_lat)
                else:
                    out_pixc_vec_tag = compute_closest_polygon_with_KDtree(in_lon, in_lat, prior_geoms, prior_id)

                # ATTENTION : Different results !!
                # print(self.compute_pixc_vec_tag_with_influence_area_map(in_lon, in_lat) == compute_closest_polygon_with_KDtree(in_lon, in_lat, prior_geoms, prior_id))
                # compute_pixc_vec_tag_with_influence_area_map is more precise.
                # compute_closest_polygon_with_KDtree less precise because computes the distance between pixels and polygon coordinates and not polygon edges.

                # Sort out_prior_id by decreasing area intersection
                sorted_idx = sorted(range(len(area_intersection)), key=lambda k: area_intersection[k], reverse=True)
                out_prior_id = [prior_id[idx] for idx in sorted_idx]

                # Print number of pixels and lakedb_id
                unique, counts = np.unique(out_pixc_vec_tag, return_counts = True)
                for i, unique_val in enumerate(unique):
                    logger.debug("%d pixels of current lake belong to a priori lake id %s " % (counts[i], unique_val))

                retour_1 = out_prior_id
                retour_2 = out_pixc_vec_tag

        return retour_1, retour_2
        
    def getPointPriorID(self, p_lon, p_lat):
        """
        Compute lakedb_id of pixel of coordinate (p_lon, p_lat)

        :param p_lon: pixels longitude
        :type p_lon: float
        :param p_lat: pixels latitude
        :type p_lat: float
        :return: lakedb_id
        :type in_poly: string
        """
        logger = logging.getLogger(self.__class__.__name__)

        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(p_lon, p_lat)
        self.influence_lake_layer.SetSpatialFilter(point)

        if self.influence_lake_layer.GetFeatureCount() == 1 :
            feat = self.influence_lake_layer.GetNextFeature()
            lakedb_id = feat.GetField(self.lake_db_id)
        else :
            logger.debug("Wrong lakedb_id from influence area map database for pixel of coordinates %f %f" %(p_lon, p_lat))
            lakedb_id = ""

        return lakedb_id


    def compute_pixc_vec_tag_with_influence_area_map(self, in_lon, in_lat):
        """
        Compute lakedb_id for every pixel of concerned lake in the case of more than one match with a priori database.

        :param in_lon: array of longitudes
        :type in_lon: numpy array of floats
        :param in_lat: array of latitude
        :type in_lat: numpy array of floats
        :return: list of lakedb_id
        :type p_lon: list of string
        """

        prior_id_list = [self.getPointPriorID(p_lon, p_lat) for p_lon, p_lat in zip(in_lon, in_lat)]

        return prior_id_list
#######################################


class LakeDbSqlite(LakeDbShp):
    """
        class LakeDbSqlite
    """
    def __init__(self, in_lake_db_filename, in_poly=None):
        """
        Constructor
        
        :param in_lake_db_filename: full path of the lake a priori database
        :type in_lake_db_filename: string
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon


        Variables of the object:
        
         - filename / string: full path of the lake a priori database
         - driver / string: driver to access to file
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Lake DB = %s", in_lake_db_filename)
        # Get config file
        cfg = service_config_file.get_instance()
        # get lake_db_id parameter
        self.lake_db_id = cfg.get("DATABASES", "LAKE_DB_ID")
        
        # Init with values
        self.filename = in_lake_db_filename  # Full path of the lake a priori database
        self.lake_db = None # store the lake database in SQLite format
        self.data_source = None # store the memory layer (shp file) from the lake database

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
        
        # #######################################################################################
        # # Open the SQLite database with ogr
        # inDriver = ogr.GetDriverByName('SQLite')
        #
        # inDB = inDriver.Open(self.filename, 0)
        #
        # # Create the output ogr layer (in memory)
        # memDriver = ogr.GetDriverByName('MEMORY')  # Memory driver
        # self.data_source = memDriver.CreateDataSource('memData')
        #
        # # Open the memory datasource with write access
        # tmp = memDriver.Open('memData', 1)
        #
        # # Copy the SQLite layer to memory
        # pipes_mem = self.data_source.CopyLayer(inDB.GetLayer('lake'), 'lake', ['OVERWRITE=YES'])
        #
        # self.layer = self.data_source.GetLayer('lake')
        # #######################################################################################

        logger = logging.getLogger(self.__class__.__name__)
        cfg = service_config_file.get_instance()
        # Transform 3D geometry into 2D geometry (necessary for spatialite query)
        in_poly.FlattenTo2D()

        # Open the SQLite database and define the connector
        self.db_conn = sqlite3.connect(self.filename, timeout=10)

        # Load spatialite extension
        self.db_conn.enable_load_extension(True)
        self.db_conn.execute('SELECT load_extension("mod_spatialite")')

        # Define the cursor
        self.db_cur = self.db_conn.cursor()

        if cfg.get('LOGGING', 'logFileLevel') == 'DEBUG':
            (lakes_nb,) = self.db_cur.execute('SELECT count(*) from lake').fetchone()
            logger.debug("[LakeDbSqlite] {} lakes stored in database".format(lakes_nb))

        # Create an output datasource in memory
        mem_driver = ogr.GetDriverByName('MEMORY')  # Memory driver

        # Open the memory datasource with write access
        self.data_source = mem_driver.CreateDataSource('memData')

        # Set spatial projection
        srs = ogr.osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        # Creating memory layer
        lyr = self.data_source.CreateLayer(str('lake_db'), srs=srs, geom_type=ogr.wkbPolygon)
        lyr.CreateField(ogr.FieldDefn(self.lake_db_id, ogr.OFTString))

        # Define the layer
        lyr_defn = lyr.GetLayerDefn()

        if in_poly is not None:
            
            self.db_cur.execute("SELECT {}, AsText(geometry) \
                                 FROM lake \
                                 WHERE MBRIntersects(GeomFromText('{}'), lake.geometry);".format(self.lake_db_id, in_poly.ExportToWkt()))

            for row in self.db_cur:
                
                # Create empty feature/entity
                out_feat = ogr.Feature(lyr_defn)

                # Fill feature with ID and geometry from SQLite request
                out_feat.SetField(self.lake_db_id, str(row[0]))
                multi_poly = ogr.CreateGeometryFromWkt(row[1])
                out_feat.SetGeometry(multi_poly)

                lyr.CreateFeature(out_feat)

                out_feat.Destroy()

            # Get memory layer
            self.layer = self.data_source.GetLayer('lake_db')

            logger.info("%d lakes after focus over studied area" % self.layer.GetFeatureCount())

        # Close spatialite database
        self.db_conn.close()




#######################################


class CollectCoordinates:
    """
        class CollectCoordinates
    """
    def __init__(self, multi_polygon):
        """
        Constructor.

        :param multi_polygon: list of lists containing polygons
        :type multi_polygon: list
        :param list_coords: numpy 2D-array of lon/lat coordinates of the multi polygon

        """
        
        self.list_coords = np.array([])

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
                

def compute_closest_polygon_with_KDtree(in_lon, in_lat, prior_geoms, prior_id):
    """
    Associate to each pixc_vec coordinate (in_lon, in_lat) the closest prior lake and its id

    :param in_lon: improved longitude of PixC related to in_poly
    :type in_lon: 1D array of float
    :param in_lat: improved latitude of PixC related to in_poly
    :type in_lat: 1D array of float
    :param prior_geom_coords: List of coordinates of polygons of prior lake database selected
    :type prior_geom_coords: 2D array of float
    :param prior_id: List of prior ID from lake DB
    :type prior_id: 1D array of str

    :return: list of the closest prior_id associated to the (in_lon, in_lat) points
    :rtype: list of str
    """
    prior_geom_coords = []
    for geom in prior_geoms:
        # Get the coordinates with json library
        json_geom = geom.ExportToJson()
        multi_polygon = json.loads(json_geom)['coordinates']

        multi_polygon_array = np.array([])
        # Loop on multi_polygon in case of multi polygon geometry
        for poly in multi_polygon:
            multi_polygon_array = np.append(multi_polygon_array, poly[:-1])

        # Reshape the 1-D array into 2-D array
        multi_polygon_array = multi_polygon_array.reshape(int(multi_polygon_array.shape[0] / 2), 2)
        prior_geom_coords.append(multi_polygon_array)  # Add the coordinates of the cur_geom to the list

    lon_coords, lat_coords, prior_id_list = np.array([]), np.array([]), np.array([], dtype=object)
    for coords, id in zip(prior_geom_coords, prior_id):
        lon_coords = np.append(lon_coords, coords[:,0])
        lat_coords = np.append(lat_coords, coords[:,1])
        prior_id_list = np.append(prior_id_list, len(coords[:,0])*[id]) # Fill the associated prior_id list of the lon/lat coordinates

    # Project coordinates to UTM before compute distances
    X_coords_utm, Y_coords_utm, utm_code = my_tools.get_utm_coords_from_lonlat(lon_coords, lat_coords)
    X_point, Y_point, utm_code = my_tools.get_utm_coords_from_lonlat(in_lon, in_lat)

    # Cdist computation
    # dist_mat = cdist(np.stack((X_coords_utm, Y_coords_utm), axis=-1), np.stack((X_point, Y_point), axis=-1))
    # cdist_idx = np.argmin(dist_mat, axis = 0)

    # Build the K-d tree
    tree_utm = KDTree(list(zip(X_coords_utm, Y_coords_utm)))

    # Built the list point for the query
    points_list_utm = np.vstack((X_point, Y_point)).T

    # Apply K-d tree and get the result: distance and index
    _, kd_tree_idx_utm = tree_utm.query(points_list_utm)

    # return prior_id_list[cdist_idx]
    return prior_id_list[kd_tree_idx_utm]


#######################################



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
    
    lake_db = LakeDbShp("C:\\Users\\pottierc\\Documents\\no_save\\data\\BD_lacs\\20181002_EU\\apriori_db_lakes_EU.shp", poly)
    print(lake_db.lake_layer.GetFeatureCount())
    print()
    for lake in lake_db.lake_layer:
        print(lake.GetField(str("lake_id")))
    print()
    print(lake_db.get_ref_values("44008001651"))
    lake_db.close_db()
