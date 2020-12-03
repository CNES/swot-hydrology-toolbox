# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2019/05/17:version initiale.
# VERSION:2.0.0:DM:#91:2020/07/03:Poursuite industrialisation
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
import os

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
        - db_format / str: database format (=None if no DB used; ="shp" for shapefile format; ="sqlite" for SQLite format)
        - lakedb_id_name / str: fieldname of identifier (=primary key) of lake table
        - lakedb_id_type / str: type of identifier (=primary key) in lake table
        - basindb_id_name / str: fieldname of identifier (=primary key) in basin table
        - pld_names / str: fieldname of lake names in lake DB
        - pld_grand / str: fieldname of GRanD identifier in lake DB
        - pld_max_wse / str: fieldname of maximum water surface elevation for lakes in DB
        - pld_max_area / str: fieldname of maximum area for lakes in DB
        - pld_ref_date / str: fieldname of reference date for storage change in DB
        - pld_ref_ds / str: fieldname of reference storage change in DB
        - pld_storage / str: fieldname of storage in DB
        - lake_layer / osgeo.ogr.Layer: layer of lake table in PLD
        - lake_ds / osgeo.ogr.DataSource: associated DataSource
        - list_lakeid / set: list of lake_id of the PLD lakes located over the retrieved subset of PLD
        - influence_lake_flag / boolean: flag indicating if the influence lake geometries are used or not
        - influence_lake_layer / osgeo.ogr.Layer: layer of lake influence area table in PLD
        - influence_lake_ds / osgeo.ogr.DataSource: associated DataSource
        - basin_flag / boolean: flag indicating if the basin geometries are used or not
        - basin_layer / osgeo.ogr.Layer: layer of basin table in PLD
        - basin_ds / osgeo.ogr.DataSource: associated DataSource
        """
        
        # Type of DB
        self.db_format = None

        # Fieldnames
        self.lakedb_id_name = my_var.PLD_FIELD_LAKE_ID  # Fieldname for lake identifier
        self.lakedb_id_type = ""  # Type of lake identifier
        self.basindb_id_name = my_var.PLD_FIELD_BASIN_ID  # Fieldname for basin identifier
        self.pld_names = my_var.PLD_FIELD_LAKE_NAMES
        self.pld_grand = my_var.PLD_FIELD_LAKE_GRAND_ID
        self.pld_max_wse = my_var.PLD_FIELD_LAKE_MAX_WSE
        self.pld_max_area = my_var.PLD_FIELD_LAKE_MAX_AREA
        self.pld_ref_date = my_var.PLD_FIELD_LAKE_REF_DATE
        self.pld_ref_ds = my_var.PLD_FIELD_LAKE_REF_DS
        self.pld_storage = my_var.PLD_FIELD_LAKE_STORAGE
        self.set_list_pld_infos()

        # Lake table layer and DataSource
        self.lake_layer = None
        self.lake_ds = None
        self.list_lakeid = set()

        # Influence area table layer and DataSource
        self.influence_lake_flag = False  # Use of influence lake geometries (default = False)
        self.influence_lake_layer = None
        self.influence_lake_ds = None

        # Basin table layer and dataSource
        self.basin_flag = False  # Use of basin geometries (default = False)
        self.basin_layer = None
        self.basin_ds = None

    def close_db(self):
        """
        Close database
        """
        pass

    # ----------------------------------------
    
    def set_list_pld_infos(self):
        """
        Group lake DB available attribute names in a list
        """
        
        self.pld_infos = []
        
        if self.pld_names is not None:
            self.pld_infos.append(self.pld_names)
        
        if self.pld_grand is not None:
            self.pld_infos.append(self.pld_grand)
        
        if self.pld_max_wse is not None:
            self.pld_infos.append(self.pld_max_wse)
        
        if self.pld_max_area is not None:
            self.pld_infos.append(self.pld_max_area)
        
        if self.pld_ref_date is not None:
            self.pld_infos.append(self.pld_ref_date)
        
        if self.pld_ref_ds is not None:
            self.pld_infos.append(self.pld_ref_ds)
        
        if self.pld_storage is not None:
            self.pld_infos.append(self.pld_storage)

    # ----------------------------------------
    
    def init_fields_name_and_type(self, test_fields=True):
        """
        Init type of lake identifier and enentually test existence of needed fieldnames
        
        :param test_fields: flag to test existence of needed fieldnames
        :type: boolean
        """
        
        # 1 - Retrieve list of available fieldnames and their type
        dict_field_type = my_tools.get_layer_fields_name_and_type(self.lake_layer)
        
        # 2 - Set type of lake identifier
        self.lakedb_id_type = dict_field_type[self.lakedb_id_name]
        
        # 3 - Test existence of fieldnames and set to None if they don't exist
        if test_fields:
            self.pld_names = my_tools.test_key(dict_field_type, self.pld_names)
            self.pld_grand = my_tools.test_key(dict_field_type, self.pld_grand)
            self.pld_max_wse = my_tools.test_key(dict_field_type, self.pld_max_wse)
            self.pld_max_area = my_tools.test_key(dict_field_type, self.pld_max_area)
            self.pld_ref_date = my_tools.test_key(dict_field_type, self.pld_ref_date)
            self.pld_ref_ds = my_tools.test_key(dict_field_type, self.pld_ref_ds)
            self.pld_storage = my_tools.test_key(dict_field_type, self.pld_storage)
            self.set_list_pld_infos()

    # ----------------------------------------
    
    def set_list_lakeid(self):
        """
        Set the list of lake_id of the PLD lakes located over the retrieved subset of PLD
        """
        self.lake_layer.ResetReading()
        for cur_lake in self.lake_layer:
            self.list_lakeid.add(cur_lake.GetField(self.lakedb_id_name))
        self.lake_layer.ResetReading()

    # ----------------------------------------
    
    def get_prior_values(self, in_id):
        """
        Getter of prior geometry and wanted infos given the lake identifier
        
        :param in_id: identifier of the lake
        :type in_id: string
        
        :return: out_geom = geometry of the prior lake
        :rtype: out_geom = OGRPolygon
        :return: out_name = name of the prior lake
        :rtype: out_name = string
        :return: out_grand = GRanD identifier 
        :rtype: out_grand = int
        :return: out_max_wse = PLD lake maximum water surface elevation
        :rtype: out_max_wse = float
        :return: out_max_area = PLD lake maximum area
        :rtype: out_max_area = float
        :return: out_ref_date = reference date for storage change computation
        :rtype: out_ref_date = string
        :return: out_ref_ds = reference storage change for storage change computation
        :rtype: out_ref_ds = float
        :return: out_storage = PLD lake maximum storage value
        :rtype: out_storage = float
        """
        
        # 0 - Init output variables
        out_geom = None
        out_name = None
        out_grand = None
        out_max_wse = None
        out_max_area = None
        out_ref_date = None
        out_ref_ds = None
        out_storage = None

        if self.lake_layer: # In case a PLD is used
            
            # Init dictionary of available PLD infos
            dict_pld_info = {}
            for item in self.pld_infos:
                dict_pld_info[item] = None
            
            # 1 - Select feature given its identifier
            if self.lakedb_id_type == "String":
                self.lake_layer.SetAttributeFilter("%s like '%s'" % (self.lakedb_id_name, str(in_id)))
            else :
                self.lake_layer.SetAttributeFilter("%s = %s" % (self.lakedb_id_name, str(in_id)))
            pld_lake_feat = self.lake_layer.GetNextFeature()
            
            # 2.1 - Retrieve geometry
            out_geom = pld_lake_feat.GetGeometryRef().Clone()
            # 2.2 - Retrieve information when exists
            for item in self.pld_infos:
                dict_pld_info[item] = pld_lake_feat.GetField(item)

            # 3 - Release filter
            self.lake_layer.SetAttributeFilter(None)
            
            # 4 - Format output
            # 4.1 - List of names
            out_name = my_tools.get_value(dict_pld_info, self.pld_names)
            if (out_name is not None) and (out_name in ["", my_var.FV_STRING_SHP]):
                out_name = None
            # 4.2 - GRanD identifier
            out_grand = my_tools.get_value(dict_pld_info, self.pld_grand)
            if (out_grand is not None) and (out_grand < 0):
                out_grand = None
            # 4.3 - Max water surface elevation
            out_max_wse = my_tools.get_value(dict_pld_info, self.pld_max_wse)
            if (out_max_wse is not None) and (out_max_wse < 0):
                out_max_wse = None
            # 4.4- Max area
            out_max_area = my_tools.get_value(dict_pld_info, self.pld_max_area)
            if (out_max_area is not None) and (out_max_area < 0):
                out_max_area = None
            # 4.5 - Reference date
            out_ref_date = my_tools.get_value(dict_pld_info, self.pld_ref_date)
            if (out_ref_date is not None) and (out_ref_date in ["", my_var.FV_STRING_SHP]):
                out_ref_date = None
            # 4.6 - Reference data storage
            out_ref_ds = my_tools.get_value(dict_pld_info, self.pld_ref_ds)
            if (out_ref_ds is not None) and (out_ref_ds == my_var.FV_REAL):
                out_ref_ds = None
            # 4.7 - Absolute water storage
            out_storage = my_tools.get_value(dict_pld_info, self.pld_storage)
            if (out_storage is not None) and (out_storage < 0):
                out_storage = None

        return out_geom, out_name, out_grand, out_max_wse, out_max_area, out_ref_date, out_ref_ds, out_storage

    # ----------------------------------------
    
    def get_influence_area_poly(self, in_lakeid):
        """
        Get influence area polygon related to in_lakeid
        or None if the influence_area layer doesn't exist
        
        :return: out_poly = influence area polygon of input PLD lake
        :rtype: out_poly = OGRPolygon
        """
        
        out_poly = None
        
        if self.influence_lake_layer is not None:
            request = "%s = '%s'" %(self.lakedb_id_name, in_lakeid)
            self.influence_lake_layer.SetAttributeFilter(request)
            cur_feature = self.influence_lake_layer.GetNextFeature()
            out_poly = cur_feature.GetGeometryRef().Clone()
            
        return out_poly

    # ----------------------------------------
    
    def link_to_db(self, in_poly, in_lon, in_lat):
        """
        Links polygon in_poly to PLD, i.e. returns, when available, 
        the list of the ID(s) of the PLD lake(s) intersecting the input polygon,
        the list of the fraction of observed lake covered by each PLD lake,
        and the reference ID array (corresponding one-to-one with the L2_HR_PIXC) 
        by giving the ID of the closest PLD lake
        
        If in_poly corresponds to no PLD lake: return None, None, None
        
        If in_poly corresponds to 1 PLD lake: return prior_id, fraction, [prior_id * size_in_lon]
        
        If in_poly corresponds to 2 or more PLD lakes:
            - out_list_prior_id contains all PLD ID (=lake_id) sorted by overlap fraction
            - out_list_pld_overlap contrains the list of fractions of observed feature covered by each PLD
            - out_pixcvec_lakeid returns the prior ID of the closest PLD lake for each PIXC
        
        :param in_poly: polygon delineating a water body
        :type in_poly: OGRPolygon
        :param in_lon: improved longitude of PixC related to in_poly
        :type in_lon: 1D array of float
        :param in_lat: improved latitude of PixC related to in_poly
        :type in_lat: 1D array of float

        :return: out_list_prior_id = list of lake identifiers from the a priori database that intersect in_poly, ordered by overlaping
                                     area decreasing
        :rtype: list of string
        :return: out_list_pld_overlap = list of fractions of observed feature covered by each PLD lake identified in out_prior_id_list, 
                                        ordered by fraction decreasing
        :rtype: list of string
        :return: out_pixcvec_lakeid = lake identifiers from the a priori database for each point of the PixC corresponding one-to-one 
                                      with (in_lon, in_lat)
        :rtype: list of string
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Init output variables
        out_list_prior_id = []
        out_list_pld_overlap = []
        out_pixcvec_lakeid = np.empty(in_lon.shape, dtype=object)
        out_pixcvec_lakeid[:] = ""

        if self.lake_layer:  # Processing only if a PLD is used
            
            # 1 - Spatial filter of the PLD over the area covered by the studied polygon
            self.lake_layer.SetSpatialFilter(in_poly)

            # 2 - Processing according to the number of PLD lakes intersecting polygon
            nb_lakes = self.lake_layer.GetFeatureCount()
            logger.debug("Current observed lake matched with %d lakes from Prior Lake Database" % (nb_lakes))

            if nb_lakes == 1:  # Easy match: polygon matches only one PLD lake

                # 2.1 - Retrieve PLD lake info
                cur_lake_bd = self.lake_layer.GetNextFeature()  # PLD lake feature
                cur_id = str(int(cur_lake_bd.GetField(self.lakedb_id_name)))  # PLD lake identifier
                logger.debug("Associated PLD identifier = %s" % (cur_id))
                
                if cur_id:  # Test but should not occur...
                    
                    # 2.2 - Save PLD identifier
                    out_list_prior_id.append(cur_id)
                
                    # 2.3 - Compute associated fraction of observed lake covered by PLD lake
                    area_obs = my_tools.get_area(in_poly)
                    geom_inter = in_poly.Intersection(cur_lake_bd.GetGeometryRef())
                    if geom_inter is not None:
                        area_inter = my_tools.get_area(geom_inter)
                        frac_inter = str(round(area_inter/area_obs*100.))
                        out_list_pld_overlap.append(frac_inter)
                
                    # 2.4 - Set lake_id attributes in PIXCVec product
                    out_pixcvec_lakeid[:] = cur_id
                    
                else:
                    logger.error("Something wrong happened in PLD: no identifier for this PLD lake!", exc_info=True)
                    raise
                    
            else:  # Many matches: polygon matches 2 or more PLD lakes

                # Init variables
                tmp_list_prior_id = []  # List of PLD identifiers
                tmp_list_pld_overlap = []  # List of fractions of observed lake covered by PLD lakes
                prior_geoms = []

                for cur_lake_bd in self.lake_layer:
                    
                    # 2.1 - Retrieve PLD lake identifier
                    cur_id = str(int(cur_lake_bd.GetField(self.lakedb_id_name)))
                    logger.debug("Associated PLD identifier = %s" % (cur_id))

                    if cur_id:  # Test but should not occur...
                    
                        # 2.2 - Save PLD identifier
                        tmp_list_prior_id.append(cur_id)
                        
                        # 2.3 - Compute associated fraction of observed lake covered by PLD lake
                        area_obs = my_tools.get_area(in_poly)
                        cur_geom = cur_lake_bd.GetGeometryRef().Clone()  # PLD lake geometry
                        geom_inter = in_poly.Intersection(cur_geom)
                        if geom_inter is not None:
                            area_inter = my_tools.get_area(geom_inter)
                            frac_inter = round(area_inter/area_obs*100.)
                            tmp_list_pld_overlap.append(frac_inter)
                            prior_geoms.append(cur_geom)
                            
                    else:
                        logger.error("Something wrong happened in PLD: no identifier for this PLD lake!", exc_info=True)
                        raise

                # Compute PIXCVec_tag and format output lists
                if len(tmp_list_prior_id) > 0:
                    
                    # 2.4 - Compute PIXCVec_tag
                    # Computation time: compute_pixcvec_lakeid_with_influence_area_map * 2,5 = compute_closest_polygon_with_kdtree
                    if self.influence_lake_layer:
                        logger.debug("Compute pixel cloud lake_id with influence area map")
                        out_pixcvec_lakeid = self.compute_pixcvec_lakeid_with_influence_area_map(in_lon, in_lat, 
                                                                                                 prior_geoms, tmp_list_prior_id)
                    else:
                        logger.debug("Compute pixel cloud lake_id with kdtree")
                        out_pixcvec_lakeid = compute_closest_polygon_with_kdtree(in_lon, in_lat, 
                                                                                 prior_geoms, tmp_list_prior_id)

                    # ATTENTION : Different results !!
                    # print(self.compute_pixcvec_lakeid_with_influence_area_map(in_lon, in_lat) == compute_closest_polygon_with_kdtree(in_lon, \
                    #                in_lat, prior_geoms, prior_id))
                    # compute_pixcvec_lakeid_with_influence_area_map is more precise.
                    # compute_closest_polygon_with_kdtree less precise because computes the distance between pixels and polygon
                    # coordinates and not polygon edges.

                    # 2.5 - Sort output prior ID and overlap fractions lists by decreasing area intersection
                    sorted_idx = np.argsort(tmp_list_pld_overlap)[::-1]
                    out_list_prior_id = [tmp_list_prior_id[idx] for idx in sorted_idx]
                    out_list_pld_overlap = [str(tmp_list_pld_overlap[idx]) for idx in sorted_idx]

                    # Print number of pixels and lake_id
                    unique, counts = np.unique(out_pixcvec_lakeid, return_counts=True)
                    for ind, unique_val in enumerate(unique):
                        logger.debug("%d pixels of current lake belong to lake_id %s " % (counts[ind], unique_val))

            self.lake_layer.SetSpatialFilter(None)  # Delete spatial filter

        return out_list_prior_id, out_list_pld_overlap, out_pixcvec_lakeid

    def compute_pixcvec_lakeid_with_influence_area_map(self, in_lon, in_lat, prior_geom_coords, list_prior_id):
        """
        Compute lake_id for each PIXC of associated observed lake in the case of more than one match with PLD.

        :param in_lon: longitudes of pixels
        :type in_lon: 1D-array of float
        :param in_lat: latitudes of pixels
        :type in_lat: 1D-array of float
        :param prior_geom_coords: coordinates of polygons of overlapping PLD lakes
        :type prior_geom_coords: list of OGRPolygon
        :param list_prior_id: prior ID of overlapping PLD lakes
        :type list_prior_id: list of str
        
        :return: out_lakedb_id_pixcvec = list of lake_id associated to each PIXC constituting observed lake
        :rtype: list of str
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 0. Init output variable
        out_lakedb_id_pixcvec = np.zeros(in_lat.size, dtype=object)
        out_lakedb_id_pixcvec[:] = ""

        # 1. Filter Influence Area following attributes
        request = "%s = '%s'" %(self.lakedb_id_name, list_prior_id[0])
        for lakedb_id in list_prior_id[1:]:
            request += " or %s = '%s'" %(self.lakedb_id_name, lakedb_id)
        self.influence_lake_layer.SetAttributeFilter(request)
        logger.debug("Filter %d influence Area with request %s" % (self.influence_lake_layer.GetFeatureCount(), request))

        # 2. Create memory layer to store pixcvec points
        lyr, ds = my_tools.load_pixels_to_mem_layer(in_lon, in_lat)

        # 3. Filter point of pixcvec with influence area
        for infl_area_feat in self.influence_lake_layer:
            lyr.SetSpatialFilter(infl_area_feat.GetGeometryRef())
            for feat_point in lyr:
                out_lakedb_id_pixcvec[feat_point.GetFID()] = infl_area_feat.GetField(self.lakedb_id_name)
            lyr.SetSpatialFilter(None)
        self.influence_lake_layer.SetAttributeFilter(None)
        ds.Destroy()

        # 4. Compute lakedb_id for pixels located out of the influence area
        unassigned_pixels = np.where(out_lakedb_id_pixcvec == "")
        nb_pt_ass_kd = unassigned_pixels[0].size
        nb_pt_ass_infl = out_lakedb_id_pixcvec.size - nb_pt_ass_kd

        if nb_pt_ass_infl > 0 :
            logger.debug("%d pixels assigned using influence map" %(nb_pt_ass_infl))
        if nb_pt_ass_kd > 0 :
            out_lakedb_id_pixcvec[unassigned_pixels] = compute_closest_polygon_with_kdtree(in_lon[unassigned_pixels], in_lat[unassigned_pixels], 
                                                                                           prior_geom_coords, list_prior_id)
            logger.debug("%d pixels assigned using kd tree" %(nb_pt_ass_kd))

        return out_lakedb_id_pixcvec

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
        
        # 2 - Set DB type
        self.db_format = "shp"

        # 3 - Open database
        self.open_db(in_lakedb_filename, in_poly)

        # 4 - Init fields name and type
        self.init_fields_name_and_type()
        
        # 5 - Set list of selected lake_id
        self.set_list_lakeid()

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

        # 1 - Init LakeDb object
        super().__init__()
        
        # 2 - Set DB type
        self.db_format = "sqlite"

        # 3 - Open database
        # 3.1 - Table lake
        self.lake_ds, self.lake_layer = self.open_db(in_lakedb_filename, my_var.PLD_TABLE_LAKE,
                                                     [self.lakedb_id_name, self.pld_names, self.pld_grand, self.pld_max_wse, \
                                                      self.pld_max_area, self.pld_ref_date, self.pld_ref_ds, self.pld_storage],
                                                     ['text', 'text', 'int9', 'float', 'float', 'text', 'float', 'float'],
                                                     in_poly=in_poly)
        # 3.2 - Table influence area
        self.influence_lake_flag = True  # Set flag indicating influence area table is used
        self.influence_lake_ds, self.influence_lake_layer = self.open_db(in_lakedb_filename, 
                                                                         my_var.PLD_TABLE_LAKE_INFL, 
                                                                         [self.lakedb_id_name],
                                                                         ['text'],
                                                                         in_poly=in_poly)
        # 3.3 - Table basin
        self.basin_flag = True  # Set flag indicating basin table is used
        self.basin_ds, self.basin_layer = self.open_db(in_lakedb_filename, 
                                                       my_var.PLD_TABLE_BASIN, 
                                                       [self.basindb_id_name],
                                                       ['text'],
                                                       in_poly=in_poly)

        # 4 - Init fields name and type
        self.init_fields_name_and_type(test_fields=False)
        
        # 5 - Set list of selected lake_id
        self.set_list_lakeid()

    # ----------------------------------------

    def open_db(self, in_lakedb_filename, in_table_name, in_field_name_list, in_field_type_list=None, in_poly=None):
        """
        Open database, optionnally spatially select polygons and copy layer to memory

        :param in_lakedb_filename: full path of PLD
        :type in_lakedb_filename: str
        :param in_table_name: name of table to load from DB
        :type in_table_name: str
        :param in_field_name_list: list of fieldnames to load from table, first element is the identifier
        :type in_field_name_list: list of str
        :param in_field_type_list: list of type of each in_field_name_list, except the identifier (=None if no attribute in addition to identifier)
        :type in_field_type_list: list of str
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon
        
        :return: out_data_source = DataSource of the specified table
        :rtype: osgeo.ogr.DataSource
        :return: out_layer = layer associated to the specified table
        :rtype: osgeo.ogr.Layer
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Loading fields %s of table %s from file %s" %(" ".join(in_field_name_list), in_table_name, in_lakedb_filename))
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
        for ind, field_type in enumerate(in_field_type_list):
            out_layer.CreateField(ogr.FieldDefn(in_field_name_list[ind], my_var.FORMAT_OGR[field_type]))
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

        # 3 - Copy selected features to output memory layer
        for db_feature in db_cursor:
            
            # 3.1 - Create empty feature
            tmp_feat = ogr.Feature(lyr_defn)

            # 3.2 - Fill feature with attributes and geometry from SQLite request
            poly = ogr.CreateGeometryFromWkt(db_feature[-1])
            for ind, fieldname in enumerate(in_field_name_list):
                tmp_feat.SetField(fieldname, str(db_feature[ind]))
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
        logger.info("- close basin table -")
        self.basin_ds.Destroy()
        

#######################################

class LakeDbDirectory(LakeDb):
    """
        class LakeDbDirectory
    """

    def __init__(self, in_lake_db_directory, in_poly=None):
        """
        Constructor

        :param in_lake_db_directory: full path of the prior lake database directory
        :type in_lake_db_directory: string
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
        logger.info("PLD directory = %s", in_lake_db_directory)

        # 1 - Init LakeDb object
        super().__init__()

        # 2 - Set DB type
        self.db_format = "sqlite"

        # 3 - Get list of PLD files
        pld_path_list = get_list_of_pld_path(in_lake_db_directory)

        # 4 - Open table basin
        self.basin_flag = True  # Set flag indicating basin table is used
        self.basin_ds, self.basin_layer = self.open_db(pld_path_list[0],
                                                       my_var.PLD_TABLE_BASIN,
                                                       [self.basindb_id_name],
                                                       ['text'],
                                                       in_poly=in_poly)

        # 5 - Select pld path corrponding to continents containted in self.basin_layer
        pld_path_list = select_pld_file_from_polygon(pld_path_list, self.basin_layer)

        # 2 - Open database
        self.lake_ds, self.lake_layer = self.open_db_multipld(pld_path_list, my_var.PLD_TABLE_LAKE,
                                                     [self.lakedb_id_name, self.pld_names, self.pld_grand, self.pld_max_wse, \
                                                      self.pld_max_area, self.pld_ref_date, self.pld_ref_ds, self.pld_storage],
                                                     ['text', 'text', 'int9', 'float', 'float', 'text', 'float', 'float'],
                                                     in_poly=in_poly)

        # 3.2 - Table influence area
        self.influence_lake_flag = True  # Set flag indicating influence area table is used
        self.influence_lake_ds, self.influence_lake_layer = self.open_db_multipld(pld_path_list,
                                                                         my_var.PLD_TABLE_LAKE_INFL,
                                                                         [self.lakedb_id_name],
                                                                         ['text'],
                                                                         in_poly=in_poly)


        # 4 - Init fields name and type
        self.init_fields_name_and_type(test_fields=False)

        # 5 - Set list of selected lake_id
        self.set_list_lakeid()

        # ----------------------------------------

    def open_db_multipld(self, pld_path_list, in_table_name, in_field_name_list, in_field_type_list=None, in_poly=None):

        """
        Open several databases, optionnally spatially select polygons and copy layer to memory

        :param pld_path_list: list of full path of PLD
        :type pld_path_list: list of str
        :param in_table_name: name of table to load from DB
        :type in_table_name: str
        :param in_field_name_list: list of fieldnames to load from table, first element is the identifier
        :type in_field_name_list: list of str
        :param in_field_type_list: list of type of each in_field_name_list, except the identifier (=None if no attribute in addition to identifier)
        :type in_field_type_list: list of str
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon

        :return: out_data_source = DataSource of the specified table
        :rtype: osgeo.ogr.DataSource
        :return: out_layer = layer associated to the specified table
        :rtype: osgeo.ogr.Layer
        """
        logger = logging.getLogger(self.__class__.__name__)

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
        for ind, field_type in enumerate(in_field_type_list):
            out_layer.CreateField(ogr.FieldDefn(in_field_name_list[ind], my_var.FORMAT_OGR[field_type]))
        # Retrieve layer definition
        lyr_defn = out_layer.GetLayerDefn()

        for pld_path in pld_path_list:
            logger.debug("Loading fields %s of table %s from folder %s" % (
            " ".join(in_field_name_list), in_table_name, pld_path))

            # 1 - Open the SQLite database
            # 1.1 - Define the connector
            db_connector = sqlite3.connect(pld_path, timeout=10)
            # 1.2 - Load spatialite extension
            db_connector.enable_load_extension(True)
            db_connector.execute('SELECT load_extension("mod_spatialite")')
            # 1.3 - Define the cursor
            db_cursor = db_connector.cursor()
            # Print info
            if cfg.get('LOGGING', 'logFileLevel') == 'DEBUG':
                (lakes_nb,) = db_cursor.execute('SELECT count(*) from %s' % (in_table_name)).fetchone()
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

            logger.debug("Loading features to memory layer")
            # 3 - Copy selected features to output memory layer
            for db_feature in db_cursor:

                # 3.1 - Create empty feature
                tmp_feat = ogr.Feature(lyr_defn)

                # 3.2 - Fill feature with attributes and geometry from SQLite request
                poly = ogr.CreateGeometryFromWkt(db_feature[-1])
                for ind, fieldname in enumerate(in_field_name_list):
                    tmp_feat.SetField(fieldname, str(db_feature[ind]))
                tmp_feat.SetGeometry(poly)

                # 3.3 - Add feature to output layer
                out_layer.CreateFeature(tmp_feat)

                # 3.4 - Close temporary feature
                tmp_feat.Destroy()

            # 5 - Close spatialite database
            db_connector.close()

        # 4 - Reset reading pointer
        out_layer.ResetReading()


        logger.info("%d features after focus over studied area" % out_layer.GetFeatureCount())
        return out_data_source, out_layer

    def open_db(self, in_lakedb_filename, in_table_name, in_field_name_list, in_field_type_list=None, in_poly=None):
        """
        Open database, optionnally spatially select polygons and copy layer to memory

        :param in_lakedb_filename: full path of PLD
        :type in_lakedb_filename: str
        :param in_table_name: name of table to load from DB
        :type in_table_name: str
        :param in_field_name_list: list of fieldnames to load from table, first element is the identifier
        :type in_field_name_list: list of str
        :param in_field_type_list: list of type of each in_field_name_list, except the identifier (=None if no attribute in addition to identifier)
        :type in_field_type_list: list of str
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon

        :return: out_data_source = DataSource of the specified table
        :rtype: osgeo.ogr.DataSource
        :return: out_layer = layer associated to the specified table
        :rtype: osgeo.ogr.Layer
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Loading fields %s of table %s from file %s" % (
        " ".join(in_field_name_list), in_table_name, in_lakedb_filename))
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
        for ind, field_type in enumerate(in_field_type_list):
            out_layer.CreateField(ogr.FieldDefn(in_field_name_list[ind], my_var.FORMAT_OGR[field_type]))
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
            (lakes_nb,) = db_cursor.execute('SELECT count(*) from %s' % (in_table_name)).fetchone()
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

        # 3 - Copy selected features to output memory layer
        for db_feature in db_cursor:

            # 3.1 - Create empty feature
            tmp_feat = ogr.Feature(lyr_defn)

            # 3.2 - Fill feature with attributes and geometry from SQLite request
            poly = ogr.CreateGeometryFromWkt(db_feature[-1])
            for ind, fieldname in enumerate(in_field_name_list):
                tmp_feat.SetField(fieldname, str(db_feature[ind]))
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

def get_list_of_pld_path(pld_directory_path):
    """
    Get list of pld file path from folder

    :pld_directory_path : path if directory containing pld files
    :pld_directory_path : list of string

    :return: list of pld path contained in pld_directory_path
    :rtype: list of str
    """
    pld_path_list = []

    for file in os.listdir(pld_directory_path):
        if file.endswith(".sqlite") :
            pld_path_list.append(os.path.join(pld_directory_path, file))

    return pld_path_list



def select_pld_file_from_polygon( pld_path_list, basin_lyr):
    """
    Select files for pld_path_list that are related to the continents contained in basin_lyr

    :pld_path_list: list of pld path related to the continents contained in basin_lyr
    :pld_path_list: list of str

    :return: list of pld path covering continents contained in basin_lyr
    :rtype: list of str
    """
    continent_set = set()
    for feat in basin_lyr :
       basin_id = feat.GetField("basin_id")
       continent_set.add(compute_continent_from_basin_id(str(basin_id)))

    basin_lyr.ResetReading()

    pld_file_list_selected = []
    for cont_id in continent_set:
        pld_file_list_selected += [pld_file for pld_file in pld_path_list if cont_id in pld_file]

    return pld_file_list_selected


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
