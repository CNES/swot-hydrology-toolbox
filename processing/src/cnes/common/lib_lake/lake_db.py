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
    :synopsis: Deal with Prior Lake Database (PLD); consider shapefile and SQLite formats
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

import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_variables as my_var


class LakeDb(object):
    """
    This class is the parent class managing Prior Lake Database
    """

    def __init__(self):
        """
        Constructor: set the general values

        Variables of the object:
        - lakedb_id_name / String: Fieldname of identifier in lake table
        - basindb_id_name / String: Fieldname of identifier in basin table
        - lake_layer / osgeo.ogr.Layer: layer of lake table in PLD
        - lake_ds / osgeo.ogr.DataSource: associated DataSource
        - influence_lake_flag / boolean: flag indicating if the influence lake geometries are used or not
        - influence_lake_layer / osgeo.ogr.Layer: layer of lake influence area table in PLD
        - influence_lake_ds / osgeo.ogr.DataSource: associated DataSource
        - basin_flag / boolean: flag indicating if the basin geometries are used or not
        - basin_layer / osgeo.ogr.Layer: layer of basin table in PLD
        - basin_ds / osgeo.ogr.DataSource: associated DataSource
        """

        # Fieldnames
        self.lakedb_id_name = ""  # Fieldname for lake identifier
        self.lakedb_id_type = ""  # Fieldname for lake identifier
        self.basindb_id_name = ""  # Fieldname for basin identifier

        # Lake table layer and DataSource
        self.lake_layer = None
        self.lake_ds = None

        # Influence area table layer and DataSource
        self.influence_lake_flag = False  # Use of influence lake geometries (default = False)
        self.influence_lake_layer = None
        self.influence_lake_ds = None

        # Basin table layer and dataSource
        self.basin_flag = False  # Use of basin geometries (default = False)
        self.basin_layer = None
        self.basin_ds = None

        self.pld_names = None
        self.pld_grand = None
        self.pld_height = None
        self.pld_area = None


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
            if self.lakedb_id_type == "String":
                self.lake_layer.SetAttributeFilter("%s like '%s'" % (self.lakedb_id_name, str(in_id)))
            else :
                self.lake_layer.SetAttributeFilter("%s = %s" % (self.lakedb_id_name, str(in_id)))

            prior_lake_feat = self.lake_layer.GetNextFeature()
            # Get name
            try:
                out_name = prior_lake_feat.GetField(self.pld_names)
            except:
                out_name = None

            # Get GRanD identifier : Dam ID from Global Reservoir and Dam (GranD) database
            try:
                out_grand = prior_lake_feat.GetField(self.pld_grand)
            except:
                out_grand = None

            # Get reference height
            try:
                out_ref_height = prior_lake_feat.GetField(self.pld_height)
            except:
                out_ref_height = None

            # Get reference area
            try:
                out_ref_area = prior_lake_feat.GetField(self.pld_area)
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
                cur_id = str(int(cur_lake_bd.GetField(self.lakedb_id_name)))

                logger.debug("A priori lakedb_id is : %s" % (cur_id))
                # Test but should not occur...
                if cur_id :
                    # Compute PIXCVec_tag
                    out_pixc_vec_tag[:] = cur_id
                    out_prior_id_list.append(cur_id)
            else:  # Many matches: polygon matches 2 or more a priori lakes

                # 2.1 - Init
                prior_id = []  # Set of prior id
                area_intersection = []  # List of area intersection between IN poly and each polygon of prior database
                prior_geoms = []

                # 2.2 - List area of intersections with a priori geometries and id of these geometries
                for cur_lake_bd in self.lake_layer:
                    # Get the a priori identifier
                    cur_id = str(int(cur_lake_bd.GetField(self.lakedb_id_name)))

                    # Test but should not occur...
                    if cur_id :
                        # Get geometry
                        cur_geom = cur_lake_bd.GetGeometryRef().Clone()
                        # Compute exact area of intersection
                        intersection = in_poly.Intersection(cur_geom)
                        if intersection is not None:
                            area_intersection.append(intersection.GetArea())  # Add the intersection area
                            prior_id.append(cur_id)  # Add prior ID to set_prior_id
                            prior_geoms.append(cur_geom)

                # 2.3 - Put output in good format
                if len(prior_id) > 0:
                    # Computation time : compute_pixc_vec_tag_with_influence_area_map * 2,5 = compute_closest_polygon_with_kdtree
                    if self.influence_lake_layer:
                        logger.debug("Compute pixel cloud lakedb_id with influence area map")
                        out_pixc_vec_tag = self.compute_pixc_vec_tag_with_influence_area_map(in_lon, in_lat, prior_geoms,
                                                                                             prior_id)
                    else:
                        logger.debug("Compute pixel cloud lakedb_id with kdtree")
                        out_pixc_vec_tag = compute_closest_polygon_with_kdtree(in_lon, in_lat, prior_geoms, prior_id)

                    # ATTENTION : Different results !!
                    # print(self.compute_pixc_vec_tag_with_influence_area_map(in_lon, in_lat) == compute_closest_polygon_with_kdtree(in_lon, in_lat, prior_geoms, prior_id))
                    # compute_pixc_vec_tag_with_influence_area_map is more precise.
                    # compute_closest_polygon_with_kdtree less precise because computes the distance between pixels and polygon coordinates and not polygon edges.

                    # Sort out_prior_id by decreasing area intersection
                    sorted_idx = sorted(range(len(area_intersection)), key=lambda k: area_intersection[k], reverse=True)
                    out_prior_id_list = [prior_id[idx] for idx in sorted_idx]

                    # Print number of pixels and lake_id
                    unique, counts = np.unique(out_pixc_vec_tag, return_counts=True)
                    for i, unique_val in enumerate(unique):
                        logger.debug("%d pixels of current lake belong to a priori lake id %s " % (counts[i], unique_val))

            self.lake_layer.SetSpatialFilter(None)  # Delete spatial filter

        return out_prior_id_list, out_pixc_vec_tag

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
        lakedb_id_pixcvec = np.zeros(in_lat.size, dtype=object)
        lakedb_id_pixcvec[:] = ""

        # 1. Filter Influence Area following attributs
        logger = logging.getLogger(self.__class__.__name__)

        request = "%s = '%s'" %(self.lakedb_id_name, prior_id_list[0])
        for lakedb_id in prior_id_list[1:]:
            request += " or %s = '%s'" %(self.lakedb_id_name, lakedb_id)
        self.influence_lake_layer.SetAttributeFilter(request)
        logger.debug("Filter %d influence Area with request %s" % (self.influence_lake_layer.GetFeatureCount(), request))

        # 2. Create memory layer to store pixcvec points
        lyr, ds = load_pixcvec_to_memory_layer(in_lon, in_lat)

        # 3. Filter point of pixcvec with influence area
        for infl_area_feat in self.influence_lake_layer:
            lyr.SetSpatialFilter(infl_area_feat.GetGeometryRef())
            for feat_point in lyr:
                lakedb_id_pixcvec[feat_point.GetFID()] = infl_area_feat.GetField(self.lakedb_id_name)
            lyr.SetSpatialFilter(None)
        self.influence_lake_layer.SetAttributeFilter(None)
        ds.Destroy()

        # 4. Compute lakedb_id for pixels located out of the influence area
        unassigned_pixels = np.where(lakedb_id_pixcvec == "")
        nb_pt_ass_kd = unassigned_pixels[0].size
        nb_pt_ass_infl = lakedb_id_pixcvec.size - nb_pt_ass_kd

        if nb_pt_ass_infl > 0 :
            logger.debug("%d pixels assigned using influence map" %(nb_pt_ass_infl))
        if nb_pt_ass_kd > 0 :
            lakedb_id_pixcvec[unassigned_pixels] = compute_closest_polygon_with_kdtree(in_lon[unassigned_pixels], in_lat[unassigned_pixels], prior_geom_coords, prior_id_list)
            logger.debug("%d pixels assigned using kd tree" %(nb_pt_ass_kd))

        return lakedb_id_pixcvec

    def link_poly_to_basin(self, in_poly):
        """
        Link a polygon to a list of basins(s) by considering intersection of both

        :param in_poly: polygon to link to a basin
        :type in_poly: ogr.Polygon

        :return: out_basin_list = list of basin(s) associated to the input polygon
        :rtype: list of string
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Init output
        out_basin_list = []
        
        if self.basin_layer:
            
            # 1 - Compute intersection
            self.basin_layer.SetSpatialFilter(in_poly)

            # 2 - Get continent name
            if self.basin_layer.GetFeatureCount() == 0:
                out_basin_list.append("000")  # Ocean
                
            elif self.basin_layer.GetFeatureCount() == 1:
                feat = self.basin_layer.GetNextFeature()
                basin_id = str(feat.GetField(self.basindb_id_name))
                out_basin_list.append(basin_id)
                
            else :
                area_intersection = []
                for feat in self.basin_layer:
                    basin_id = str(feat.GetField(self.basindb_id_name))
                    inter = feat.GetGeometryRef().Intersection(in_poly)
                    area_intersection.append(inter.GetArea())
                    out_basin_list.append(basin_id)
                # Sort out_basin_list by area intersection decreasing 
                sorted_idx = sorted(range(len(area_intersection)), key=lambda k: area_intersection[k], reverse=True)
                out_basin_list = [out_basin_list[idx] for idx in sorted_idx]
                
        else :
            out_basin_list.append("010")  # If no input continent is given, code is 010 => no continent used

        logger.info(out_basin_list)
        return out_basin_list

    def link_poly_to_continent(self, in_poly):
        """
        Link a polygon to a list of continent(s) by considering intersection of both

        :param in_poly: polygon to link to a continent
        :type in_poly: ogr.Polygon

        :return: out_continent_list = list of continent(s) associated to the input polygon
        :rtype: list of string
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Init output
        out_continent_list = set()
        
        # Compute list of basins related to the input polygon
        basin_list = self.link_poly_to_basin(in_poly)

        # Compute list of continents
        for basin_id in basin_list:
            continent = compute_continent_from_basin_id(basin_id)
            if continent:
                out_continent_list.add(continent)

        logger.info(';'.join(out_continent_list))
        return ';'.join(out_continent_list)

    def close_db(self):
        """
        Close database
        """
        pass
#######################################


class LakeDbShp(LakeDb):
    """
    This class heritates from the main class to manage PLD in shapefile format
    """

    def __init__(self, in_lakedb_filename, in_poly=None):
        """
        Constructor

        :param in_lakedb_filename: full path of PLD
        :type in_lakedb_filename: string
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Lake DB in SHAPEFILE format = %s", in_lakedb_filename)

        # 1 - Init LakeDb object
        super().__init__()

        # 2 - Get config file
        cfg = service_config_file.get_instance()
        # Get lakedb_id parameter
        self.lakedb_id_name = cfg.get("DATABASES", "LAKE_DB_ID")

        # 3 - Open database
        self.open_db(in_lakedb_filename, in_poly)

        layer_defn = self.lake_layer.GetLayerDefn()
        layer_field = {}
        for i in range(layer_defn.GetFieldCount()) :
            field_name = layer_defn.GetFieldDefn(i).GetName()
            field_type = layer_defn.GetFieldDefn(i).GetFieldTypeName(layer_defn.GetFieldDefn(i).GetType())
            layer_field[field_name] = field_type

        # Init field names and types
        self.lakedb_id_type = layer_field[self.lakedb_id_name]
        if my_var.PLD_FIELD_NAMES in layer_field.keys():
            self.pld_names = my_var.PLD_FIELD_NAMES
        if my_var.PLD_FIELD_GRAND_ID in layer_field.keys():
            self.pld_grand = my_var.PLD_FIELD_GRAND_ID
        if my_var.PLD_FIELD_REF_HEIGHT in layer_field.keys():
            self.pld_height = my_var.PLD_FIELD_REF_HEIGHT
        if my_var.PLD_FIELD_REF_AREA in layer_field.keys():
            self.pld_area = my_var.PLD_FIELD_REF_AREA

# ----------------------------------------

    def open_db(self, in_lakedb_filename, in_poly=None):
        """
        Open PLD, optionnally spatially select polygons, and copy lake_layer to memory

        :param in_lakedb_filename: full path of PLD
        :type in_lakedb_filename: string
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")

        # 1 - Open shapefile in read-only access
        shp_driver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Shapefile driver
        shp_data_source = shp_driver.Open(in_lakedb_filename, 0)

        # 2 - Get the lake_layer
        layer = shp_data_source.GetLayer()
        logger.info("%d lakes stored in PLD %s" % (layer.GetFeatureCount(), layer.GetName()))

        # 3 - Select subset among PLD lakes using in_poly
        if in_poly is not None:
            layer.SetSpatialFilter(in_poly)
            logger.info("%d lakes after focus over studied area", layer.GetFeatureCount())

        # 4 - Create an output DataSource in memory
        mem_driver = ogr.GetDriverByName('MEMORY')  # Memory driver
        self.lake_ds = mem_driver.CreateDataSource('memData')

        # 5 - Open the memory DataSource with write access
        mem_driver.Open('memData', 1)

        # 6 - Copy the lake_layer to memory
        self.lake_ds.CopyLayer(layer, 'lake_db')

        # 7 - Get memory lake_layer
        self.lake_layer = self.lake_ds.GetLayer()

        # 8 - Close shapefile
        shp_data_source.Destroy()

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
    This class heritates from the main class to manage PLD in SQLite format
    """

    def __init__(self, in_lakedb_filename, in_poly=None):
        """
        Constructor

        :param in_lakedb_filename: full path of PLD
        :type in_lakedb_filename: string
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Lake DB in SQLITE format = %s", in_lakedb_filename)

        # 1. Get tables and fields names
        lake_table_name = my_var.PLD_TABLE_LAKE
        lake_infl_table_name = my_var.PLD_TABLE_LAKE_INFL
        basin_table_name = my_var.PLD_TABLE_BASIN
        self.lakedb_id_name = my_var.PLD_FIELD_LAKE_ID
        self.basindb_id_name = my_var.PLD_FIELD_BASIN_ID

        self.pld_names = my_var.PLD_FIELD_NAMES
        self.pld_grand = my_var.PLD_FIELD_GRAND_ID
        self.pld_height = my_var.PLD_FIELD_REF_HEIGHT
        self.pld_area = my_var.PLD_FIELD_REF_AREA

        # 2 - Open database
        # 2.1 - Table lake
        self.lake_ds, self.lake_layer = self.open_db(in_lakedb_filename,
                                                     lake_table_name,
                                                     [self.lakedb_id_name, self.pld_names, self.pld_grand, self.pld_height, self.pld_area],
                                                     in_poly)
        # 2.2 - Table influence area
        self.influence_lake_ds, self.influence_lake_layer = self.open_db(in_lakedb_filename, 
                                                                         lake_infl_table_name, 
                                                                         [self.lakedb_id_name],
                                                                         in_poly)
        # 2.3 - Table basin
        self.basin_ds, self.basin_layer = self.open_db(in_lakedb_filename, 
                                                       basin_table_name, 
                                                       [self.basindb_id_name], 
                                                       in_poly)

        # 3 - Set flags indicating influence area and basin tables are used
        self.influence_lake_flag = True
        self.basin_flag = True

        # Init field type
        layer_defn = self.lake_layer.GetLayerDefn()
        layer_field = {}
        for i in range(layer_defn.GetFieldCount()):
            field_name = layer_defn.GetFieldDefn(i).GetName()
            field_type = layer_defn.GetFieldDefn(i).GetFieldTypeName(layer_defn.GetFieldDefn(i).GetType())
            layer_field[field_name] = field_type

        self.lakedb_id_type = layer_field[self.lakedb_id_name]


    # ----------------------------------------

    def open_db(self, in_lakedb_filename, in_table_name, in_field_name_list, in_poly=None):
        """
        Open database, optionnally spatially select polygons and copy layer to memory

        :param in_lakedb_filename: full path of PLD
        :type in_lakedb_filename: str
        :param in_table_name: name of table to load from DB
        :type in_table_name: str
        :param in_field_name_list: list of fieldnames to load from table, first element is the identifier
        :type in_field_name_list: list of str
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon
        
        :return: out_data_source = DataSource of the specified table
        :rtype: osgeo.ogr.DataSource
        :return: out_layer = layer associated to the specified table
        :rtype: osgeo.ogr.Layer
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        cfg = service_config_file.get_instance()
        
        # 0 - Init output in memory
        mem_driver = ogr.GetDriverByName('MEMORY')  # Memory driver
        # 0.1 - Open the memory DataSource with write access
        out_data_source = mem_driver.CreateDataSource('memData')
        # 0.2 - Set spatial projection
        srs = ogr.osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        # 0.3 - Create memory layer
        out_layer = out_data_source.CreateLayer(str('layer'), srs=srs, geom_type=ogr.wkbPolygon)
        # 0.4 - Create needed fields
        out_layer.CreateField(ogr.FieldDefn(in_field_name_list[0], ogr.OFTString))
        if in_table_name == my_var.PLD_TABLE_LAKE:
            out_layer.CreateField(ogr.FieldDefn(my_var.PLD_FIELD_NAMES, ogr.OFTString))
            out_layer.CreateField(ogr.FieldDefn(my_var.PLD_FIELD_GRAND_ID, ogr.OFTInteger))
            out_layer.CreateField(ogr.FieldDefn(my_var.PLD_FIELD_REF_HEIGHT, ogr.OFTReal))
            out_layer.CreateField(ogr.FieldDefn(my_var.PLD_FIELD_REF_AREA, ogr.OFTReal))
        # Retrieve layer definition
        lyr_defn = out_layer.GetLayerDefn()

        # 1 - Open the SQLite database
        # 1.1 - Define the connector
        db_connector = sqlite3.connect(in_lakedb_filename, timeout=10)
        # 1.2 - Load spatialite extension
        db_connector.enable_load_extension(True)
        db_connector.execute('SELECT load_extension("mod_spatialite")')
        # 1.3 - Define the cursor
        db_cursor = db_connector.cursor()
        # Print info
        if cfg.get('LOGGING', 'logFileLevel') == 'DEBUG':
            (lakes_nb,) = db_cursor.execute('SELECT count(*) from %s' %(in_table_name)).fetchone()
            logger.debug("%d features stored in table <%s>" % (lakes_nb, in_table_name))

        # 2 - Select subset among PLD lakes using in_poly
        if in_poly is not None:
            in_poly.FlattenTo2D()  # Transform 3D geometry into 2D geometry (necessary for spatialite query)
            cmd = "SELECT %s, AsText(geometry) FROM %s WHERE MBRIntersects(GeomFromText('%s'), %s.geometry);" % (
                ",".join(in_field_name_list), in_table_name, in_poly.ExportToWkt(), in_table_name)
            db_cursor.execute(cmd)
        else:
            cmd = "SELECT %s, AsText(geometry) FROM %s ;" % (",".join(in_field_name_list), in_table_name)
            db_cursor.execute(cmd)

        # 3 - Copy selected lakes to output memory layer
        for lake in db_cursor:
            
            # 3.1 - Create empty feature
            tmp_feat = ogr.Feature(lyr_defn)

            # 3.2 - Fill feature with ID and geometry from SQLite request
            tmp_feat.SetField(in_field_name_list[0], str(lake[0]))
            if in_table_name == my_var.PLD_TABLE_LAKE:
                tmp_feat.SetField(my_var.PLD_FIELD_NAMES, str(lake[1]))
                tmp_feat.SetField(my_var.PLD_FIELD_GRAND_ID, str(lake[2]))
                tmp_feat.SetField(my_var.PLD_FIELD_REF_HEIGHT, str(lake[3]))
                tmp_feat.SetField(my_var.PLD_FIELD_REF_AREA, str(lake[4]))
            poly = ogr.CreateGeometryFromWkt(lake[-1])
            tmp_feat.SetGeometry(poly)

            # 3.3 - Add feature to output layer
            out_layer.CreateFeature(tmp_feat)

            # 3.4 - Close temporary feature
            tmp_feat.Destroy()

        # 4 - Reset reading pointer
        out_layer.ResetReading()

        # 5 - Close spatialite database
        db_connector.close()

        logger.info("%d features after focus over studied area" % out_layer.GetFeatureCount())
        return out_data_source, out_layer

    # ----------------------------------------

    def close_db(self):
        """
        Close database
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Close lake table
        logger.info("- close lake table -")
        self.lake_ds.Destroy()
        
        # Close influence area table
        logger.info("- close influence area table -")
        self.influence_lake_ds.Destroy()
        
        # Close basin table
        logger.info("- close lake table -")
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

def load_pixcvec_to_memory_layer(in_lon, in_lat):
    # Create Memory layer
    mem_driver = ogr.GetDriverByName('MEMORY')  # Memory driver
    # Open the memory datasource with write access
    ds = mem_driver.CreateDataSource('memData')
    # Set spatial projection
    srs = ogr.osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    # Creating memory layer
    lyr = ds.CreateLayer(str('point'), srs=srs, geom_type=ogr.wkbPoint)

    for (lon, lat) in zip(in_lon, in_lat):
        p = ogr.Geometry(ogr.wkbPoint)
        p.AddPoint(lon, lat)
        point = ogr.Feature(lyr.GetLayerDefn())
        point.SetGeometry(p)
        lyr.CreateFeature(point)

    return lyr, ds
