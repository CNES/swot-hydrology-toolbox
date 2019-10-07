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
    :synopsis: Deal with operationnal prior lake database; consider shapefile and sqlite formats
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
import sqlite3

import cnes.common.service_config_file as service_config_file
import cnes.common.lib.my_variables as my_var
import cnes.common.lib.my_tools as my_tools

class LakeDb(object):
    """
        This class deals with lakedb.
    """

    def __init__(self):
        """
        Constructor: set the general values

        Variables of the object:
        - lake_db_id / String: Fieldname of lake id in lakedb
        - basin_db_id / String: Fieldname of basin id in basindb
        - lake_layer / osgeo.ogr.Layer: lake_layer of a priori lake database
        - lake_ds / osgeo.ogr.DataSource: datasource of a priori lake database
        - influence_lake_layer / osgeo.ogr.Layer: lake_influence_layer of a priori lake database
        - influence_lake_ds / osgeo.ogr.DataSource: datasource of influence of  a priori lake database
        - basin_layer / osgeo.ogr.Layer: basin_layer of a priori lake database
        - basin_ds / osgeo.ogr.DataSource: datasource of basin of a priori lake database
        """

        self.lake_db_id = ""
        self.basin_db_id = ""

        self.lake_layer = None
        self.lake_ds = None

        self.influence_lake_layer = None
        self.influence_lake_ds = None

        self.basin_layer = None
        self.basin_ds = None

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

        if self.lake_layer:
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

        else :
            out_name = None
            out_grand = None
            out_ref_height = None
            out_ref_area = None

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

        :return: out_prior_id_list = list of lake identifiers from the a priori database that intersect in_poly
        :rtype: list of string
        :return: out_pixc_vec_tag = lake identifiers from the a priori database for each point of the PixC corresponding one-to-one with (in_lon, in_lat)
        :rtype: list of string
        """
        logger = logging.getLogger(self.__class__.__name__)
        out_prior_id_list = []
        out_pixc_vec_tag = np.empty(in_lon.shape, dtype=object)
        out_pixc_vec_tag[:] = ""

        if self.lake_layer:
            # 1 - Spatial filter of the lake DB over the area covered by the studied polygon
            self.lake_layer.SetSpatialFilter(in_poly)

            # 2 - Processing according to the number of a priori lakes intersecting polygon
            nb_lakes = self.lake_layer.GetFeatureCount()
            logger.debug("Current lake matched with %d lakes from a priori lake database" % (nb_lakes))

            if nb_lakes == 1:  # Easy match: polygon matches only one a priori lake

                cur_lake_bd = self.lake_layer.GetNextFeature()
                cur_id = cur_lake_bd.GetField(self.lake_db_id)

                logger.debug("A priori lakedb_id is : %s" % (cur_id))
                # Test but should not occur...
                if cur_id :
                    # Compute PIXCVec_tag
                    out_pixc_vec_tag[:] = str(cur_id)
                    out_prior_id_list.append(str(cur_id))
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
                    if cur_id :
                        # Get geometry
                        cur_geom = cur_lake_bd.GetGeometryRef().Clone()
                        # Compute exact area of intersection
                        intersection = in_poly.Intersection(cur_geom)
                        if intersection is not None:
                            area_intersection.append(intersection.GetArea())  # Add the intersection area
                            prior_id.append(str(cur_id))  # Add prior ID to set_prior_id
                            prior_geoms.append(cur_geom)

                # 2.3 - Put output in good format
                if len(prior_id) > 0:
                    # Computation time : compute_pixc_vec_tag_with_influence_area_map * 2,5 = compute_closest_polygon_with_kdtree
                    if self.influence_lake_layer:
                        out_pixc_vec_tag = self.compute_pixc_vec_tag_with_influence_area_map(in_lon, in_lat, prior_geoms,
                                                                                             prior_id)
                    else:
                        out_pixc_vec_tag = compute_closest_polygon_with_kdtree(in_lon, in_lat, prior_geoms, prior_id)

                    # ATTENTION : Different results !!
                    # print(self.compute_pixc_vec_tag_with_influence_area_map(in_lon, in_lat) == compute_closest_polygon_with_kdtree(in_lon, in_lat, prior_geoms, prior_id))
                    # compute_pixc_vec_tag_with_influence_area_map is more precise.
                    # compute_closest_polygon_with_kdtree less precise because computes the distance between pixels and polygon coordinates and not polygon edges.

                    # Sort out_prior_id by decreasing area intersection
                    sorted_idx = sorted(range(len(area_intersection)), key=lambda k: area_intersection[k], reverse=True)
                    out_prior_id_list = [prior_id[idx] for idx in sorted_idx]

                    # Print number of pixels and lake_db_id
                    unique, counts = np.unique(out_pixc_vec_tag, return_counts=True)
                    for i, unique_val in enumerate(unique):
                        logger.debug("%d pixels of current lake belong to a priori lake id %s " % (counts[i], unique_val))

            self.lake_layer.SetSpatialFilter(None)  # Delete spatial filter

        return out_prior_id_list, out_pixc_vec_tag

    def get_point_prior_id(self, p_lon, p_lat, prior_geom_coords, prior_id_list):
        """
        Compute lakedb_id of pixel of coordinate (p_lon, p_lat)

        :param p_lon: pixels longitude
        :type p_lon: float
        :param p_lat: pixels latitude
        :type p_lat: float
        :param prior_geom_coords: List of coordinates of polygons of prior lake database selected
        :type prior_geom_coords: 2D array of float
        :param prior_id_list: List of prior ID from lake DB
        :type prior_id_list: 1D array of str
        :return: lakedb_id
        :type in_poly: string
        """
        logger = logging.getLogger(self.__class__.__name__)

        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(p_lon, p_lat)

        self.influence_lake_layer.SetSpatialFilter(point)

        if self.influence_lake_layer.GetFeatureCount() == 1:
            feat = self.influence_lake_layer.GetNextFeature()
            lakedb_id = str(int(feat.GetField(self.lake_db_id)))
        else:
            logger.debug(
                "Wrong lakedb_id from influence area map database for pixel of coordinates %f %f" % (p_lon, p_lat))
            lakedb_id = ""

        if lakedb_id not in prior_id_list:
            lakedb_id = compute_closest_polygon_with_kdtree([p_lon], [p_lat], prior_geom_coords, prior_id_list)[0]

        return lakedb_id

    def compute_pixc_vec_tag_with_influence_area_map(self, in_lon, in_lat, prior_geom_coords, prior_id_list):
        """
        Compute lakedb_id for every pixel of concerned lake in the case of more than one match with a priori database.

        :param in_lon: array of longitudes
        :type in_lon: numpy array of floats
        :param in_lat: array of latitude
        :type in_lat: numpy array of floats
        :param prior_geom_coords: List of coordinates of polygons of prior lake database selected
        :type prior_geom_coords: 2D array of float
        :param prior_id_array: List of prior ID from lake DB
        :type prior_id_array: 1D array of str
        :return: list of lakedb_id
        :type p_lon: list of string
        """
        prior_id_list = [self.get_point_prior_id(p_lon, p_lat, prior_geom_coords, prior_id_list) for p_lon, p_lat in
                         zip(in_lon, in_lat)]
        return prior_id_list

    def link_poly_to_basin(self, in_poly):
        """
        Link a polygon to a list of basins(s) by considering intersection of both

        :param in_poly: polygon to link to a basin
        :type in_poly: ogr.Polygon

        :return: list of basis(s) associated to polygon
        :rtype: list of string
        """
        continent_list = []
        if self.basin_layer :
            # 1 - Compute intersection
            self.basin_layer.SetSpatialFilter(in_poly)

            # 2 - Get continent name

            if self.basin_layer.GetFeatureCount() == 0:
                continent_list.append("000")
            elif self.basin_layer.GetFeatureCount() == 1:
                feat = self.basin_layer.GetNextFeature()
                basin_id = str(feat.GetField(self.basin_db_id))
                continent_list.append(basin_id)
            else :
                continent_list = []
                area_intersection = []
                for feat in self.basin_layer:
                    print(self.basin_db_id)
                    basin_id = str(feat.GetField(self.basin_db_id))
                    inter = feat.GetGeometryRef().Intersection(in_poly)
                    area_intersection.append(inter.GetArea())
                    continent_list.append(basin_id)
                # Sort out_prior_id by decreasing area intersection
                sorted_idx = sorted(range(len(area_intersection)), key=lambda k: area_intersection[k], reverse=True)
                continent_list = [continent_list[idx] for idx in sorted_idx]
        else :
            continent_list.append("010") # If no input continent is given, code is 010 => no continent used

        return continent_list

    def link_poly_to_continent(self, in_poly):
        """
        Link a polygon to a list of continent(s) by considering intersection of both

        :param in_poly: polygon to link to a continent
        :type in_poly: ogr.Polygon

        :return: list of continent(s) associated to polygon
        :rtype: list of string
        """

        basins_list = self.link_poly_to_basin(in_poly)

        # Get continent name
        out_continent = set()
        for basin_id in basins_list:
            if compute_continent_from_basin_id(basin_id):
                out_continent.add(compute_continent_from_basin_id(basin_id))

        return ';'.join(out_continent)

    def close_db(self):
        pass

#######################################

class LakeDbShp(LakeDb):
    """
        class LakeDbShp
    """

    def __init__(self, in_lake_db_filename, in_poly=None):
        """
        Constructor

        :param in_lake_db_filename: full path of the prior lake database
        :type in_lake_db_filename: string
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon

        Variables of the object:
            - lakedb_id / String: Fieldname of lake id in lakedb
            - lake_layer / osgeo.ogr.Layer: lake_layer of a priori lake database
            - lake_ds / osgeo.ogr.DataSource: datasource of a priori lake database
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Lake DB = %s", in_lake_db_filename)

        # 1 - Init LakeDb
        super().__init__()

        # 2 - Get config file
        cfg = service_config_file.get_instance()
        # Get lakedb_id parameter
        self.lake_db_id = cfg.get("DATABASES", "LAKE_DB_ID")

        # 3 - Open database
        self.lake_ds, self.lake_layer = self.open_shp(in_lake_db_filename, in_poly)  # Open Lake database

    # ----------------------------------------

    def open_shp(self, in_file_path, in_poly=None):
        """
        Open database, optionnally spatially select polygons and copy lake_layer to memory

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

        # 2 - Get the lake_layer
        layer = shp_data_source.GetLayer()
        logger.info("%d features stored in database %s" % (layer.GetFeatureCount(), layer.GetName()))

        # 3 - Select some lakes among BD using in_poly
        if in_poly is not None:
            layer.SetSpatialFilter(in_poly)
            logger.info("%d features after focus over studied area", layer.GetFeatureCount())

        # 4 - Create an output datasource in memory
        mem_driver = ogr.GetDriverByName('MEMORY')  # Memory driver
        data_source = mem_driver.CreateDataSource('memData')

        # 5 - Open the memory datasource with write access
        mem_driver.Open('memData', 1)

        # 6 - Copy the lake_layer to memory
        data_source.CopyLayer(layer, 'lake_db')

        # 7 - Get memory lake_layer
        layer = data_source.GetLayer()

        # 8 - Close shapefile
        shp_data_source.Destroy()
        logger.info("Stop")
        return data_source, layer

    # ----------------------------------------

    def close_db(self):
        """
        Close database
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- close database -")
        self.lake_ds.Destroy()


#######################################

class LakeDbSqlite(LakeDb):
    """
        class LakeDbSqlite
    """

    def __init__(self, in_lake_db_filename, in_poly=None):
        """
        Constructor

        :param in_lake_db_filename: full path of the prior lake database
        :type in_lake_db_filename: string
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon

        Variables of the object:
            - lake_db_id / String: Fieldname of lake id in lakedb
            - lake_layer / osgeo.ogr.Layer: lake_layer of a priori lake database
            - lake_ds / osgeo.ogr.DataSource: datasource of a priori lake database
            - influence_lake_layer / osgeo.ogr.Layer: lake_influence_layer of a priori lake database
            - influence_lake_ds / osgeo.ogr.DataSource: datasource of influence of  a priori lake database
            - influence_map_flag / Bool: flag that determine if LakeDb uses influence map
            - basin_layer / osgeo.ogr.Layer: basin_layer of a priori lake database
            - basin_ds / osgeo.ogr.DataSource: datasource of basin of a priori lake database
            - basin_flag / Bool: flag that determine if LakeDb uses basin
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Lake DB = %s", in_lake_db_filename)

        # 1. Get tables and fields names
        lakedb_table_name = my_var.lakedb_table_name
        lake_infl_table_name = my_var.lake_infl_table_name
        basin_table_name = my_var.basin_table_name
        self.lake_db_id = my_var.lake_db_id
        self.basin_db_id = my_var.basin_db_id

        # 2 - Open database
        self.lake_ds, self.lake_layer = self.open_db(in_lake_db_filename, lakedb_table_name, self.lake_db_id, in_poly)
        self.influence_lake_ds, self.influence_lake_layer = self.open_db(in_lake_db_filename, lake_infl_table_name, self.lake_db_id,
                                                                         in_poly)
        self.basin_ds, self.basin_layer = self.open_db(in_lake_db_filename, basin_table_name, self.basin_db_id, in_poly)

        # 3 - Set db flags
        self.influence_map_flag = True
        self.basin_flag = True

    # ----------------------------------------

    def open_db(self, in_lake_db_filename, table_name, field_name, in_poly=None):
        """
        Open database, optionnally spatially select polygons and copy layer to memory

        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon
        :param table_name: table_name to load from file
        :type table_name: String
        :param field_name: field_name to load from table
        :type field_name: String
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon
        """
        logger = logging.getLogger(self.__class__.__name__)
        cfg = service_config_file.get_instance()
        
        # Transform 3D geometry into 2D geometry (necessary for spatialite query)
        in_poly.FlattenTo2D()

        # Open the SQLite database and define the connector
        self.db_conn = sqlite3.connect(in_lake_db_filename, timeout=10)

        # Load spatialite extension
        self.db_conn.enable_load_extension(True)
        self.db_conn.execute('SELECT load_extension("mod_spatialite")')

        # Define the cursor
        self.db_cur = self.db_conn.cursor()

        if cfg.get('LOGGING', 'logFileLevel') == 'DEBUG':
            (lakes_nb,) = self.db_cur.execute('SELECT count(*) from lake').fetchone()
            logger.debug(" %d features stored in table %s of database" % (lakes_nb, table_name))

        # Create an output datasource in memory
        mem_driver = ogr.GetDriverByName('MEMORY')  # Memory driver

        # Open the memory datasource with write access
        ds = mem_driver.CreateDataSource('memData')

        # Set spatial projection
        srs = ogr.osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        # Creating memory layer
        lyr = ds.CreateLayer(str('layer'), srs=srs, geom_type=ogr.wkbPolygon)
        lyr.CreateField(ogr.FieldDefn(field_name, ogr.OFTString))

        # Define the layer
        lyr_defn = lyr.GetLayerDefn()

        if in_poly is not None:
            cmd = "SELECT %s, AsText(geometry) FROM %s WHERE MBRIntersects(GeomFromText('%s'), %s.geometry);" % (
                field_name, table_name, in_poly.ExportToWkt(), table_name)
            self.db_cur.execute(cmd)
        else:
            cmd = "SELECT %s, AsText(geometry) FROM %s ;" % (field_name, table_name)
            self.db_cur.execute(cmd)

        for row in self.db_cur:
            # Create empty feature/entity
            out_feat = ogr.Feature(lyr_defn)

            # Fill feature with ID and geometry from SQLite request
            out_feat.SetField(field_name, str(row[0]))
            multi_poly = ogr.CreateGeometryFromWkt(row[1])
            out_feat.SetGeometry(multi_poly)

            lyr.CreateFeature(out_feat)

            out_feat.Destroy()

        lyr.ResetReading()
        # Get memory lake_layer
        logger.info("%d features after focus over studied area" % lyr.GetFeatureCount())

        # Close spatialite database
        self.db_conn.close()

        return ds, lyr

    # ----------------------------------------

    def close_db(self):
        """
        Close database
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- close database -")
        self.lake_ds.Destroy()
        if self.influence_map_flag:
            self.influence_lake_ds.Destroy()
        if self.basin_flag:
            self.basin_ds.Destroy()

#######################################


def compute_closest_polygon_with_kdtree(in_lon, in_lat, prior_geoms, prior_id):
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
        lon_coords = np.append(lon_coords, coords[:, 0])
        lat_coords = np.append(lat_coords, coords[:, 1])
        prior_id_list = np.append(prior_id_list, len(coords[:, 0]) * [
            id])  # Fill the associated prior_id list of the lon/lat coordinates

    # Project coordinates to UTM before compute distances
    x_coords_utm, y_coords_utm, utm_code = my_tools.get_utm_coords_from_lonlat(lon_coords, lat_coords)
    x_point, y_point, utm_code = my_tools.get_utm_coords_from_lonlat(in_lon, in_lat)

    # Cdist computation
    # dist_mat = cdist(np.stack((x_coords_utm, y_coords_utm), axis=-1), np.stack((x_point, y_point), axis=-1))
    # cdist_idx = np.argmin(dist_mat, axis = 0)

    # Build the K-d tree
    tree_utm = KDTree(list(zip(x_coords_utm, y_coords_utm)))

    # Built the list point for the query
    points_list_utm = np.vstack((x_point, y_point)).T

    # Apply K-d tree and get the result: distance and index
    _, kd_tree_idx_utm = tree_utm.query(points_list_utm)

    # code commenté return prior_id_list[cdist_idx]
    return prior_id_list[kd_tree_idx_utm]


def compute_continent_from_basin_id(in_basin_id):
    """
    Compute continent from basin ID (HydroBASINS nomenclature)
    * 0 for Ocean
    * 1 for Africa
    * 2 for Europe
    * 3 for Siberia
    * 4 for Asia
    * 5 for Australia
    * 6 for South America
    * 7 for North America
    * 8 for Arctic (North America)
    * 9 for Greenland
    * 10 for no continent

    :param in_basin_id: basin identifier
    :type in_basin_id: string

    :return: 2 letters related to the continent
    :rtype: string
    """

    retour = ""

    if in_basin_id == "010":
        retour = ""
    elif in_basin_id == "000":
        retour = "OCEAN"
    elif in_basin_id.startswith("1"):
        retour = "AF"
    elif in_basin_id.startswith("2"):
        retour = "EU"
    elif in_basin_id.startswith("3"):
        retour = "SI"
    elif in_basin_id.startswith("4"):
        retour = "AS"
    elif in_basin_id.startswith("5"):
        retour = "AU"
    elif in_basin_id.startswith("6"):
        retour = "SA"
    elif in_basin_id.startswith("7"):
        retour = "NA"
    elif in_basin_id.startswith("8"):
        retour = "AR"
    elif in_basin_id.startswith("9"):
        retour = "GR"
    else:
        retour = "xx"

    return retour


def compute_basin_id_from_continent(in_continent):
    """
    Compute continent from basin ID (HydroBASINS nomenclature)
    * 0 for Ocean
    * 1 for Africa
    * 2 for Europe
    * 3 for Siberia
    * 4 for Asia
    * 5 for Australia
    * 6 for South America
    * 7 for North America
    * 8 for Arctic (North America)
    * 9 for Greenland
    * 10 for no continent

    :param in_continent: continent code
    :type in_continent: string

    :return: pfaf_continent id
    :rtype: string
    """

    retour = ""
    if in_continent == "OCEAN":
        retour = "0"
    elif in_continent == "AF":
        retour = "1"
    elif in_continent == "EU":
        retour = "2"
    elif in_continent == "SI":
        retour = "3"
    elif in_continent == "AS":
        retour = "4"
    elif in_continent == "AU":
        retour = "5"
    elif in_continent == "SA":
        retour = "6"
    elif in_continent == "NA":
        retour = "7"
    elif in_continent == "AR":
        retour = "8"
    elif in_continent == "GR":
        retour = "9"
    elif in_continent == "":
        retour = "10"
    else:
        retour = ""

    return retour
