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
# VERSION:3.0.0:DM:#91:2021/03/12:Poursuite industrialisation
# VERSION:3.1.0:DM:#91:2021/05/21:Poursuite industrialisation
# VERSION:3.2.0:DM:#91:2021/10/27:Poursuite industrialisation
# VERSION:4.0.0:DM:#91:2022/05/05:Poursuite industrialisation
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
import cnes.common.lib_lake.storage_change as storage_change


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
        - tiledb_id_name / str: fieldname for tile reference
        - tiledb_cont / str: fieldname of list of IDs of continents intersecting the tile (case of LakeTile processing)
        - tiledb_az_0_line / str: fieldname of azimuth 0 line in pass tile table
        - tiledb_az_max_line / str: fieldname of azimuth max line in pass tile table
        - pld_names / str: fieldname of lake names in lake DB
        - pld_res_id / str: fieldname of reservoir identifier in lake DB
        - pld_lon / float: fieldname of longitude of the centroid of the lake
        - pld_lat / float: fieldname of latitude of the centroid of the lake
        - pld_max_wse / str: fieldname of maximum water surface elevation for lakes in DB
        - pld_max_wse_u / str: fieldname of uncertainty over maximum water surface elevation for lakes in DB
        - pld_max_area / str: fieldname of maximum area for lakes in DB
        - pld_max_area_u / str: fieldname of uncertainty over maximum area for lakes in DB
        - pld_ref_date / str: fieldname of reference date for storage change in DB
        - pld_ref_ds / str: fieldname of reference storage change in DB
        - pld_storage / str: fieldname of storage in DB
        - pld_ice_clim / str: fieldname of ice_clim flag in DB
        - pld_ice_dyn / str: fieldname of ice_dyn flag in DB
        - pld_reach / str: fieldname of list of Prior River Database reach_id intersecting the lake
        - pld_overlap_avg_l / str: fieldname of average overlap of lake and tile of current pass and swath L
        - pld_nb_tiles_l / str: fieldname of number of tiles covered by a lake over swath L
        - pld_overlap_avg_r / str: fieldname of average overlap of lake and tile of current pass and swath R
        - pld_nb_tiles_r / str: fieldname of number of tiles covered by a lake over swath R
        - big_lake_geom = list of string : list of lakes_id considered as too big to be process in LakeSP
        - lake_layer / osgeo.ogr.Layer: layer of lake table in PLD
        - lake_ds / osgeo.ogr.DataSource: associated DataSource
        - list_lakeid / set: list of lake_id of the PLD lakes located over the retrieved subset of PLD
        - flag_ocean / boolean: (case of LakeTile processing) =True if tile is entirely over ocean
        - influence_lake_flag / boolean: flag indicating if the influence lake geometries are used or not
        - influence_lake_layer / osgeo.ogr.Layer: layer of lake influence area table in PLD
        - influence_lake_ds / osgeo.ogr.DataSource: associated DataSource
        - basin_flag / boolean: flag indicating if the basin geometries are used or not
        - basin_layer / osgeo.ogr.Layer: layer of basin table in PLD
        - basin_ds / osgeo.ogr.DataSource: associated DataSource
        - az_0_geom / osgeo.ogr.Geometry: geometries of lakes located at the bottom of tile
        - az_max_geom / osgeo.ogr.Geometry: geometries of lakes located at the top of tile
        - az_0_and_max_geom / osgeo.ogr.Geometry: geometries of lakes located at the top and bottom of tile
        """
        try:
            # Get instance of service config file
            cfg = service_config_file.get_instance()
            # Save min_overlap parameter
            self.min_overlap = cfg.getfloat('CONFIG_PARAMS', 'MIN_OVERLAP')
        except:
            self.min_overlap = 2.0
        
        # Type of DB
        self.db_format = None

        # Fieldnames
        self.lakedb_id_name = my_var.PLD_FIELD_LAKE_ID  # Fieldname for lake identifier
        self.lakedb_id_type = ""  # Type of lake identifier
        self.basindb_id_name = my_var.PLD_FIELD_BASIN_ID  # Fieldname for basin identifier
        self.tiledb_id_name = my_var.PLD_FIELD_TILE_REF  # Fieldname for tile reference
        self.tiledb_cont = my_var.PLD_FIELD_TILE_CONT
        self.tiledb_az_0_line = my_var.PLD_FIELD_TILE_AZ_0_LINE
        self.tiledb_az_max_line = my_var.PLD_FIELD_TILE_AZ_MAX_LINE
        self.pld_names = my_var.PLD_FIELD_LAKE_NAMES
        self.pld_res_id = my_var.PLD_FIELD_LAKE_RES_ID
        self.pld_lon = my_var.PLD_FIELD_LAKE_LON
        self.pld_lat = my_var.PLD_FIELD_LAKE_LAT
        self.pld_max_wse = my_var.PLD_FIELD_LAKE_MAX_WSE
        self.pld_max_wse_u = my_var.PLD_FIELD_LAKE_MAX_WSE_U
        self.pld_max_area = my_var.PLD_FIELD_LAKE_MAX_AREA
        self.pld_max_area_u = my_var.PLD_FIELD_LAKE_MAX_AREA_U
        self.pld_ref_date = my_var.PLD_FIELD_LAKE_REF_DATE
        self.pld_ref_ds = my_var.PLD_FIELD_LAKE_REF_DS
        self.pld_storage = my_var.PLD_FIELD_LAKE_STORAGE
        self.pld_ice_clim = my_var.PLD_FIELD_ICE_CLIM
        self.pld_ice_dyn = my_var.PLD_FIELD_ICE_DYN
        self.pld_reach = my_var.PLD_FIELD_LAKE_REACH
        self.pld_overlap_avg_l = my_var.PLD_FIELD_LAKE_OVERLAP_AVG_L
        self.pld_nb_tiles_l = my_var.PLD_FIELD_LAKE_NB_TILES_L
        self.pld_overlap_avg_r = my_var.PLD_FIELD_LAKE_OVERLAP_AVG_R
        self.pld_nb_tiles_r = my_var.PLD_FIELD_LAKE_NB_TILES_R
        self.set_list_pld_infos()
        self.big_lake_geom = []

        # Lake table layer and DataSource
        self.lake_layer = None
        self.lake_ds = None
        self.list_lakeid = set()
        self.flag_ocean = False

        # Influence area table layer and DataSource
        self.influence_lake_flag = False  # Use of influence lake geometries (default = False)
        self.influence_lake_layer = None
        self.influence_lake_ds = None

        # Basin table layer and dataSource
        self.basin_flag = False  # Use of basin geometries (default = False)
        self.basin_layer = None
        self.basin_ds = None

        # Geometries of lakes located at the top or bottom of tile
        self.az_0_geom = ogr.Geometry(ogr.wkbMultiPolygon)
        self.az_max_geom = ogr.Geometry(ogr.wkbMultiPolygon)
        self.az_0_and_max_geom = ogr.Geometry(ogr.wkbMultiPolygon)

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
        
        if self.pld_res_id is not None:
            self.pld_infos.append(self.pld_res_id)
        
        if self.pld_lon is not None:
            self.pld_infos.append(self.pld_lon)
        
        if self.pld_lat is not None:
            self.pld_infos.append(self.pld_lat)
        
        if self.pld_max_wse is not None:
            self.pld_infos.append(self.pld_max_wse)
        
        if self.pld_max_wse_u is not None:
            self.pld_infos.append(self.pld_max_wse_u)
        
        if self.pld_max_area is not None:
            self.pld_infos.append(self.pld_max_area)
        
        if self.pld_max_area_u is not None:
            self.pld_infos.append(self.pld_max_area_u)
        
        if self.pld_ref_date is not None:
            self.pld_infos.append(self.pld_ref_date)
        
        if self.pld_ref_ds is not None:
            self.pld_infos.append(self.pld_ref_ds)
        
        if self.pld_storage is not None:
            self.pld_infos.append(self.pld_storage)
        
        if self.pld_ice_clim is not None:
            self.pld_infos.append(self.pld_ice_clim)
        
        if self.pld_ice_dyn is not None:
            self.pld_infos.append(self.pld_ice_dyn)
        
        if self.pld_reach is not None:
            self.pld_infos.append(self.pld_reach)

        if self.pld_overlap_avg_l is not None:
            self.pld_infos.append(self.pld_overlap_avg_l)

        if self.pld_nb_tiles_l is not None:
            self.pld_infos.append(self.pld_nb_tiles_l)

        if self.pld_overlap_avg_r is not None:
            self.pld_infos.append(self.pld_overlap_avg_r)

        if self.pld_nb_tiles_r is not None:
            self.pld_infos.append(self.pld_nb_tiles_r)
    
    def init_fields_name_and_type(self, test_fields=True):
        """
        Init type of lake identifier and eventually test existence of needed fieldnames
        
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
            self.pld_res_id = my_tools.test_key(dict_field_type, self.pld_res_id)
            self.pld_lon = my_tools.test_key(dict_field_type, self.pld_lon)
            self.pld_lat = my_tools.test_key(dict_field_type, self.pld_lat)
            self.pld_max_wse = my_tools.test_key(dict_field_type, self.pld_max_wse)
            self.pld_max_wse_u = my_tools.test_key(dict_field_type, self.pld_max_wse_u)
            self.pld_max_area = my_tools.test_key(dict_field_type, self.pld_max_area)
            self.pld_max_area_u = my_tools.test_key(dict_field_type, self.pld_max_area_u)
            self.pld_ref_date = my_tools.test_key(dict_field_type, self.pld_ref_date)
            self.pld_ref_ds = my_tools.test_key(dict_field_type, self.pld_ref_ds)
            self.pld_storage = my_tools.test_key(dict_field_type, self.pld_storage)
            self.pld_ice_clim = my_tools.test_key(dict_field_type, self.pld_ice_clim)
            self.pld_ice_dyn = my_tools.test_key(dict_field_type, self.pld_ice_dyn)
            self.pld_reach = my_tools.test_key(dict_field_type, self.pld_reach)
            self.pld_overlap_avg_l = my_tools.test_key(dict_field_type, self.pld_overlap_avg_l)
            self.pld_nb_tiles_l = my_tools.test_key(dict_field_type, self.pld_nb_tiles_l)
            self.pld_overlap_avg_r = my_tools.test_key(dict_field_type, self.pld_overlap_avg_r)
            self.pld_nb_tiles_r = my_tools.test_key(dict_field_type, self.pld_nb_tiles_r)

        # 4 - Compute again self.pld_info, LakeTile/SP ot Avg pld_infos are differentes
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

    def build_border_geometry(self, in_az_0_line, in_az_max_line):
        """
        Compute Multipolygon from all polygon intersecting in_az_0_line or in_az_max_line and a buffer of in_az_[0-max]_line. Border polygones are stored
        in self.az_0_geom, self.az_max_geom and self.az_0_and_max_geom

        :param in_az_0_line: Line with azimuth = 0
        :type: OGR Geometry Linestring
        :param in_az_max_line: Line with azimuth = azimuth max
        :type: OGR Geometry Linestring
        """
        logger = logging.getLogger(self.__class__.__name__)

        self.az_0_geom = in_az_0_line.Buffer(0.000001)
        self.az_max_geom = in_az_max_line.Buffer(0.000001)
        self.az_0_and_max_geom = self.az_0_geom.Union(self.az_max_geom)

        if in_az_0_line:
            self.lake_layer.ResetReading()
            self.lake_layer.SetSpatialFilter(in_az_0_line)

            for feat in self.lake_layer:
                geom = feat.GetGeometryRef().Buffer(0)
                lake_id = feat.GetField(self.lakedb_id_name)
                self.az_0_geom = self.az_0_geom.Union(geom)
                logger.debug("Lake_id %s will be processed in LakeSP because pld_lake intersects tile border" % (lake_id))
                self.list_lakeid.remove(lake_id)

            self.lake_layer.SetSpatialFilter(None)
            self.lake_layer.ResetReading()

        if in_az_max_line:
            self.lake_layer.ResetReading()
            self.lake_layer.SetSpatialFilter(in_az_max_line)

            for feat in self.lake_layer:
                geom = feat.GetGeometryRef().Buffer(0)
                lake_id = feat.GetField(self.lakedb_id_name)

                if geom.Intersects(self.az_0_geom):
                    self.az_0_and_max_geom = self.az_0_and_max_geom.Union(geom)
                else :
                    self.az_max_geom = self.az_max_geom.Union(geom)
                logger.debug("Lake_id %s will be processed in LakeSP because pld_lake intersects tile border" % (lake_id))
                if lake_id in self.list_lakeid:  # in the case of a lake crossing in_az_0_line and in_az_max_line, lake is already removed
                    self.list_lakeid.remove(lake_id)

            self.lake_layer.SetSpatialFilter(None)
            self.lake_layer.ResetReading()

    # ----------------------------------------

    def get_influence_area_poly(self, in_lakeid):
        """
        Get influence area polygon related to in_lakeid
        or None if the influence_area layer doesn't exist
        
        :return: out_poly = influence area polygon of input PLD lake
        :rtype: out_poly = OGRPolygon
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        out_poly = None
        
        if self.influence_lake_layer is not None:
            request = "%s = '%s'" %(self.lakedb_id_name, in_lakeid)
            self.influence_lake_layer.SetAttributeFilter(request)
            
            if self.influence_lake_layer.GetFeatureCount() == 1:
                cur_feature = self.influence_lake_layer.GetNextFeature()
                out_poly = cur_feature.GetGeometryRef().Clone()
                
            else:
                ## TODO : this case should never happend : a lake in PLD should ALWAYS have an influence area
                logger.warning("lake_id %s corresponds to %d feature in lake_influence table" % (in_lakeid, self.influence_lake_layer.GetFeatureCount()))
                logger.warning("lake_id %s corresponds to %d feature in lake_influence table" \
                               % (in_lakeid, self.influence_lake_layer.GetFeatureCount()))
                out_poly = None
                
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

            if nb_lakes != 0:
                # Init variables
                tmp_list_prior_id = []  # List of PLD identifiers
                tmp_list_pld_overlap = []  # List of fractions of observed lake covered by PLD lakes
                prior_geoms = []

                # 3 - Check number of string matches
                for cur_lake_bd in self.lake_layer:

                    # 3.1 - Retrieve PLD lake identifier and geometry
                    cur_id = str(cur_lake_bd.GetField(self.lakedb_id_name))
                    logger.debug("Associated PLD identifier = %s" % cur_id)
                    cur_geom = cur_lake_bd.GetGeometryRef().Clone()

                    if cur_id:  # Test but should not occur...
                        # 3.2 - Compute associated fraction of observed lake covered by PLD lake
                        frac_inter = self.compute_frac_inter(in_poly, cur_geom)

                        # 3.3 - Save info to output lists
                        # only if overlap percentage is important enough
                        if frac_inter > self.min_overlap:
                            tmp_list_prior_id.append(cur_id)  # PLD identifier
                            tmp_list_pld_overlap.append(str(frac_inter))  # Overlap
                            prior_geoms.append(cur_geom)
                        else:
                            logger.debug(
                                "Overlap = {}% (<{}%) => PLD lake not linked with observed lake".format(frac_inter,
                                                                                                        self.min_overlap))
                    else:
                        logger.error("Something wrong happened in PLD: no identifier for this PLD lake!", exc_info=True)
                        raise

                # 4 - Compute out_list_prior_id, out_list_pld_overlap and out_pixcvec_lakeid
                if len(tmp_list_prior_id) == 1: # Easy match: polygon strongly matches only one PLD lake
                    out_list_prior_id = tmp_list_prior_id  # PLD identifier
                    out_list_pld_overlap = tmp_list_pld_overlap  # Overlap
                    out_pixcvec_lakeid[:] = tmp_list_prior_id[0]  # Set lake_id attributes in PIXCVec product
                elif len(tmp_list_prior_id) > 1:  # In the case where in_poly intersects several lakes from PLD and more than one lake have an overlap bigger than 2%
                    # 4.2 - Compute PIXCVec_tag
                    # Computation time: compute_pixcvec_lakeid_with_influence_area_map * 2,5 = compute_closest_polygon_with_kdtree
                    if self.influence_lake_layer:
                        logger.debug("Compute PIXCVec lake_id with influence area map")
                        out_pixcvec_lakeid = self.compute_pixcvec_lakeid_with_influence_area_map(in_lon, in_lat,
                                                                                                 prior_geoms, tmp_list_prior_id)
                    else:
                        logger.debug("Compute PIXCVec lake_id with kdtree")
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
                        logger.debug("%d pixels of current observed lake belong to lake_id %s " % (counts[ind], unique_val))

            self.lake_layer.SetSpatialFilter(None)  # Delete spatial filter

        return out_list_prior_id, out_list_pld_overlap, out_pixcvec_lakeid

    def compute_frac_inter(self, lake_geom, pld_geom):
        """
        Compute intersection between lake_geom and pld_geom.

        :param lake_geom: polygon delineating a water body
        :type lake_geom: OGRPolygon
        :param pld_geom: polygon delineating a water body in Prior Lake Database
        :type pld_geom: OGRPolygon

        :return: out_frac_inter = relative overlapping area = area(intersection)/area(pld_geom)
        :rtype: out_frac_inter = float
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Compute the geometry of the intersection between the observed lake and the PLD lake
        geom_inter = lake_geom.Intersection(pld_geom)

        # 2 - Compute relative intersection area wrt observed area
        if geom_inter is not None:
            area_obs = my_tools.get_area(lake_geom)
            area_inter = my_tools.get_area(geom_inter)
            out_frac_inter = round(area_inter / area_obs * 100.)
        else:
            logger.warning("Something wrong happened in PLD: this lake doesn't intersect PLD!")
            out_frac_inter = 0.0

        return out_frac_inter

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
        request = "%s IN ('%s'" %(self.lakedb_id_name, list_prior_id[0])
        for lakedb_id in list_prior_id[1:]:
            request += ", '%s'" %(lakedb_id)
        request += ")"
        self.influence_lake_layer.SetAttributeFilter(request)

        if self.influence_lake_layer.GetFeatureCount() == 0:
            ## TODO : this case should never happend : a lake in PLD should ALWAYS have an influence area
            logger.warning("No influence area found for lake_id %s" % str(list_prior_id))
        else:
            logger.debug("Filter %d influence Area with request %s" % (self.influence_lake_layer.GetFeatureCount(), request))

        # 2. Create memory layer to store pixcvec points
        ds, lyr = my_tools.load_pixels_to_mem_layer(in_lon, in_lat)

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

        if nb_pt_ass_infl > 0:
            logger.debug("%d pixels assigned using influence map" % nb_pt_ass_infl)
        if nb_pt_ass_kd > 0:
            out_lakedb_id_pixcvec[unassigned_pixels] = compute_closest_polygon_with_kdtree(in_lon[unassigned_pixels], in_lat[unassigned_pixels], 
                                                                                           prior_geom_coords, list_prior_id)
            logger.debug("%d pixels assigned using kd tree" % nb_pt_ass_kd)

        return out_lakedb_id_pixcvec

    def link_poly_to_basin(self, in_poly):
        """
        Link a polygon to a list of basins(s) by considering intersection of both

        :param in_poly: polygon to link to a basin
        :type in_poly: ogr.Polygon

        :return: out_basin_list = list of IDs of the basin(s) associated to the input polygon
        :rtype: list of string
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Init output
        out_basin_list = []
        
        if self.basin_layer:
            
            # 1 - Compute intersection
            self.basin_layer.SetSpatialFilter(in_poly)

            # 2 - Get basin code 
            if self.basin_layer.GetFeatureCount() == 0:
                out_basin_list.append("000")  # Ocean
                
            elif self.basin_layer.GetFeatureCount() == 1:
                feat = self.basin_layer.GetNextFeature()
                basin_id = str(feat.GetField(self.basindb_id_name))
                out_basin_list.append(basin_id)
                
            else:
                area_intersection = []
                for basin_feat in self.basin_layer:
                    basin_id = str(basin_feat.GetField(self.basindb_id_name))
                    basin_geom = basin_feat.GetGeometryRef()
                    if basin_geom.IsValid():
                        inter = in_poly.Intersection(basin_geom)
                    else:
                        logger.warning("PLD basin table contains an invalid geometry with basin_id %s" % basin_id)
                        inter = in_poly.Intersection(basin_geom.Buffer(0))
                    area_intersection.append(inter.GetArea())
                    out_basin_list.append(basin_id)
                # Sort out_basin_list by area intersection decreasing 
                sorted_idx = sorted(range(len(area_intersection)), key=lambda k: area_intersection[k], reverse=True)
                out_basin_list = [out_basin_list[idx] for idx in sorted_idx]
                
        else:
            out_basin_list.append("010")  # If no input continent is given, code is 010 => no continent used

        logger.debug(out_basin_list)
        return out_basin_list

    def link_poly_to_continent_and_basin(self, in_poly):
        """
        Link a polygon to a list of continent(s) and basin(s) by considering intersection of both

        :param in_poly: polygon to link to continent(s) and basin(s)
        :type in_poly: ogr.Polygon

        :return: out_continent_id_list = two-letter identifier of the continent(s) associated to the input polygon, ordered by increasing continent_code
        :rtype: string
        :return: out_continent_code_list = 1-digit code of the continent(s) associated to the input polygon, ordered by increasing continent_code
        :rtype: string
        :return: out_basin_code_list = code of the basin(s) associated to the input polygon, ordered in increasing identifier
        :rtype: string
        """
        
        # 1 - Compute code of the basins related to the input polygon
        out_basin_code_list = "None"
        basin_code_list = self.link_poly_to_basin(in_poly) # Basin codes are ordered by descending intersection area between tile and basins
        for code in ["000", "010"]:
            if code in basin_code_list:
                basin_code_list.remove(code)
        if len(basin_code_list) > 0:
            out_basin_code_list = ';'.join(basin_code_list)
            
        if out_basin_code_list == "None":
            out_continent_id_list = "None"
            out_continent_code_list = "None"
            
        else:

            # 2 - Compute lists of continents
            # 2.1 - Compute continent code
            continent_code_list = []
            for basin_code in basin_code_list:
                if basin_code[0] not in continent_code_list:
                    continent_code_list.append(basin_code[0])
            out_continent_code_list = ';'.join(continent_code_list)
            # 2.2 - Compute continent identifier
            continent_id_list = []
            for continent_code in continent_code_list:
                continent_id = compute_continent_id_from_basin_code(continent_code)
                if continent_id:
                    continent_id_list.append(continent_id)
            out_continent_id_list = ';'.join(continent_id_list)
                
        return out_continent_id_list, out_continent_code_list, out_basin_code_list

    def remove_from_list_lake_id(self, in_poly):
        """
        Remove from list_lake_id id of lakes intersecting in_poly. In_poly contains bounding box of water bodies at the top or bottom of tile.

        :param in_poly: bounding box of water bodies at the top or bottom of tile
        :type in_poly: ogr.Polygon
        """

        logger = logging.getLogger(self.__class__.__name__)

        self.lake_layer.SetSpatialFilter(in_poly)
        for feat in self.lake_layer:
            lake_id = feat.GetField(self.lakedb_id_name)
            if lake_id in self.list_lakeid:
                logger.debug("Lake_id %s will be processed in LakeSP because obs_lake intersects tile border" %(lake_id))
                self.list_lakeid.remove(lake_id)
        self.lake_layer.SetSpatialFilter(None)

    
#######################################


class LakeDbShp(LakeDb):
    """
    This class heritates from the main class to manage PLD in shapefile format
    """

    def __init__(self, in_lakedb_filename, 
                 in_poly=None, in_az_0_line=None, in_az_max_line=None):
        """
        Constructor

        :param in_lakedb_filename: full path of lake DB shapefile
        :type in_lakedb_filename: string
        
        :param in_poly: polygon to spatially select lakes
        :type in_poly: OGRPolygon
        :param in_az_0_line: line delineating the bottom of the tile
        :type in_az_0_line: OGRLineString
        :param in_az_max_line: line delineating the top of the tile
        :type in_az_max_line: OGRLineString
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Lake DB in SHAPEFILE format = %s", in_lakedb_filename)

        # 1 - Init LakeDb object
        super().__init__()
        
        # 2 - Set DB type
        self.db_format = "shp"

        # 3 - Open database
        self.open_db(in_lakedb_filename, in_poly=in_poly)

        # 4 - Init fields name and type
        self.init_fields_name_and_type()
        
        # 5 - Set list of selected lake_id
        self.set_list_lakeid()

        if in_az_0_line and in_az_max_line:
            # 6 - Build border geometry
            self.build_border_geometry(in_az_0_line, in_az_max_line)

    # ----------------------------------------

    def open_db(self, in_lakedb_filename, in_poly=None):
        """
        Open PLD, optionnally spatially select polygons or with attribute lake_id, and copy lake_layer to memory

        :param in_lakedb_filename: full path of PLD
        :type in_lakedb_filename: string
        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")

        # 1 - Open shapefile in read-only access
        shp_driver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Shapefile driver
        shp_data_source = shp_driver.Open(in_lakedb_filename, 0)

        # 2 - Get the lake_layer
        layer = shp_data_source.GetLayer()
        logger.debug("%d lakes stored in PLD %s" % (layer.GetFeatureCount(), layer.GetName()))

        # 3 - Select subset among PLD lakes using in_poly
        if in_poly is not None:
            layer.SetSpatialFilter(in_poly)
            logger.debug("%d lakes after focus over studied area", layer.GetFeatureCount())

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
        logger.debug("- close database -")
        self.lake_ds.Destroy()


#######################################
        

class LakeDbSqlite(LakeDb):
    """
    This class heritates from the main class to manage PLD in SQLite format
    """

    def __init__(self, in_lakedb_filename, 
                 in_tile_ref=None, in_continent_id=None, in_swath=None, in_basin_id=None,
                 in_poly=None, in_az_0_line=None, in_az_max_line=None):
        """
        Constructor

        :param in_lakedb_filename: full path of PLD file
        :type in_lakedb_filename: string
        :param in_tile_ref: reference of tile on which to focus (specific to LakeTile processing, in case of use of operational PLD); format = <3-digit_tile_number><swath_side>
        :type in_tile_ref: string
        :param in_continent_id: 1-digit continent identifier on which to focus (specific to LakeSP processing, in case of use of operational PLD)
        :type in_continent_id: int
        :param in_swath: "R" or "L" for Right ou Left swath
        :type in_swath: String
        :param in_basin_id: 2-digit basin code on which to focus (specific to LakeAvg processing, in case of use of operational PLD)
        :type in_basin_id: int
        :param in_poly: polygon to spatially select lake features (in case of user lake DB) 
        :type in_poly: OGRPolygon
        :param in_az_0_line: line delineating the bottom of the tile (specific to LakeTile processing, in case of use of user lake DB)
        :type in_az_0_line: OGRLineString
        :param in_az_max_line: line delineating the top of the tile (specific to LakeTile processing, in case of use of user lake DB)
        :type in_az_max_line: OGRLineString
        
        Possible combinations:
            - in_tile_ref=not None, others=None => LakeTile processing with operational PLD
            - in_continent_id=not None, others=None => LakeSP processing with operational PLD
            - in_basin_id=not None, others=None => LakeAvg processing with operational PLD
            - in_tile_poly=not None, others=None
                - with in_az_0_line=None, in_az_max_line=None => LakeSP processing with user lake DB
                - with in_az_0_line=not None, in_az_max_line=not None => LakeTile processing with user lake DB
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Lake DB in SQLITE format = %s", in_lakedb_filename)

        # 1 - Init LakeDb object
        super().__init__()
        
        # 2 - Set DB type
        self.db_format = "sqlite"

        # 3 - Store database tables
        start_value = None
        if in_continent_id is not None:
            start_value = in_continent_id
        elif in_basin_id is not None:
            start_value = in_basin_id
        
        # 3.1 - Table pass_tile
        if in_tile_ref is not None: # LakeTile processing ; select tile geometry
            # Read table pass_tile
            pass_tile_ds, pass_tile_layer = self.open_db(in_lakedb_filename, 
                                                         my_var.PLD_TABLE_PASS_TILE,
                                                         [self.tiledb_id_name, self.tiledb_cont, self.tiledb_az_0_line, self.tiledb_az_max_line],
                                                         ['text', 'text', 'text', 'text'],
                                                         field_id=in_tile_ref)
            nb_tile_features = pass_tile_layer.GetFeatureCount()
            if nb_tile_features != 1:
                logger.error("%d features correspond to tile %s" % (nb_tile_features, in_tile_ref))
                raise
            # Get feature
            feat = pass_tile_layer.GetNextFeature()
            # Retrieve continent
            continent_value = feat.GetField(self.tiledb_cont)
            if continent_value == "0":
                self.flag_ocean = True
            # Retrieve the polygon of the tile
            poly_4selection = feat.GetGeometryRef()
            # poly_4selection can be multipolygon in the case if a tile with longitude=+=180.
            if poly_4selection.GetGeometryName() not in ["POLYGON", "MULTIPOLYGON"]:
                logger.warning("poly_4selection variable type %s invalid" %poly_4selection.GetGeometryName())
                poly_4selection = None
            # Retreive az_0 and az_max_lines
            az_0_line = ogr.CreateGeometryFromWkt(feat.GetField(self.tiledb_az_0_line))
            az_max_line = ogr.CreateGeometryFromWkt(feat.GetField(self.tiledb_az_0_line))

        elif in_swath is not None: # LakeSP processing ; select current swath pass
            pass_tile_ds, pass_tile_layer = self.open_db(in_lakedb_filename,
                                                         my_var.PLD_TABLE_PASS_TILE,
                                                         [self.tiledb_id_name],
                                                         ['text'],
                                                         field_id='%'+str(in_swath))
            nb_tile_features = pass_tile_layer.GetFeatureCount()
            if nb_tile_features < 1:
                logger.error("%d features correspond to swath %s" % (nb_tile_features, in_swath))
                raise
            poly_4selection = ogr.Geometry(ogr.wkbMultiPolygon)
            for pass_tile_feat in pass_tile_layer:
                for geom in pass_tile_feat.GetGeometryRef():
                    poly_4selection.AddGeometry(geom)
            poly_4selection = poly_4selection.UnionCascaded()
            self.pld_overlap_avg_l = None
            self.pld_nb_tiles_l = None
            self.pld_overlap_avg_r = None
            self.pld_nb_tiles_r = None
        else: # LakeAvg processing
            poly_4selection = in_poly
            # Fields pld_overlap_avg_l, pld_nb_tiles_l, pld_overlap_avg_r, pld_nb_tiles_r do not extsits in PLD by Basin. Set to None
            self.pld_overlap_avg_l = None
            self.pld_nb_tiles_l = None
            self.pld_overlap_avg_r = None
            self.pld_nb_tiles_r = None

        # 3.2 - Table lake
        field_to_load_list = [self.lakedb_id_name, self.pld_names, self.pld_res_id, \
                              self.pld_lon, self.pld_lat, \
                              self.pld_max_wse, self.pld_max_wse_u, \
                              self.pld_max_area, self.pld_max_area_u, \
                              self.pld_ref_date, self.pld_ref_ds, self.pld_storage, \
                              self.pld_ice_clim, self.pld_ice_dyn, \
                              self.pld_reach]

        field_type_to_load_list = ['text', 'text', 'int9', \
                                                      'float', 'float', \
                                                      'float', 'float', \
                                                      'float', 'float', \
                                                      'text', 'float', 'float', \
                                                      'int4', 'int4', \
                                                      'text']
        if in_tile_ref : # Only for LakeTile processing
            field_to_load_list += [self.pld_overlap_avg_l, self.pld_nb_tiles_l,
                              self.pld_overlap_avg_r, self.pld_nb_tiles_r]
            field_type_to_load_list += ['float', 'int', 'float', 'int']

        self.lake_ds, self.lake_layer = self.open_db(in_lakedb_filename, my_var.PLD_TABLE_LAKE, field_to_load_list, field_type_to_load_list,
                                                     poly=poly_4selection, start_value=start_value)
        
        # 3.2 - Table influence area
        self.influence_lake_flag = True  # Set flag indicating influence area table is used
        try:
            self.influence_lake_ds, self.influence_lake_layer = self.open_db(in_lakedb_filename, 
                                                                             my_var.PLD_TABLE_LAKE_INFL, 
                                                                             [self.lakedb_id_name],
                                                                             ['text'],
                                                                             poly=poly_4selection, start_value=start_value)
        except:
            logger.warning("No table named %s in the database" % my_var.PLD_TABLE_LAKE_INFL)
            self.influence_lake_flag = False
        
        # 3.3 - Table basin
        self.basin_flag = True  # Set flag indicating basin table is used
        self.basin_ds, self.basin_layer = self.open_db(in_lakedb_filename, 
                                                       my_var.PLD_TABLE_BASIN, 
                                                       [self.basindb_id_name, \
                                                        my_var.PLD_FIELD_BASIN_LON_MIN, my_var.PLD_FIELD_BASIN_LON_MAX, \
                                                        my_var.PLD_FIELD_BASIN_LAT_MIN, my_var.PLD_FIELD_BASIN_LAT_MAX],
                                                       ['text', 'float', 'float', 'float', 'float'],
                                                       poly=poly_4selection, start_value=in_continent_id)
        
        # 4 - Init fields name and type
        self.init_fields_name_and_type(test_fields=False)
        
        # 5 - Set list of selected lake_id
        self.set_list_lakeid()

        # 6 - Build border geometry (only for LakeTile processing)
        if in_tile_ref is not None:
            self.build_border_geometry(az_0_line, az_max_line)
        
        # 7 - Save CB basin coordinates (only for LakeAvg processing)
        self.basin_coords = None
        if in_basin_id is not None:
            cbb_basin = self.basin_layer.GetNextFeature()
            self.basin_coords = dict()
            basin_coords_keys = ["lon_min", "lon_max", "lat_min", "lat_max"]
            
            for key in basin_coords_keys:
                self.basin_coords[key] = cbb_basin.GetField(key)
            
            self.basin_layer.ResetReading()

        # 8 - Identify lakes with pld_overlap_avg_[l/r]*pld_nb_tiles_[l/r] > MAX_NB_TILES_FULL_AZ
        # Get instance of service config file
        cfg = service_config_file.get_instance()

        # only for LakeTile processing, -1 for MAX_NB_TILES_FULL_AZ parameter disable selection of big lakes
        if in_tile_ref is not None:
            MAX_NB_TILES_FULL_AZ = cfg.getint("CONFIG_PARAMS", "MAX_NB_TILES_FULL_AZ")
            if MAX_NB_TILES_FULL_AZ != -1:
                if "L" in in_tile_ref:
                    self.lake_layer.SetAttributeFilter("overlap_avg_l*nb_tiles_l > %d " % MAX_NB_TILES_FULL_AZ)
                else:
                    self.lake_layer.SetAttributeFilter("overlap_avg_r*nb_tiles_r > %d " % MAX_NB_TILES_FULL_AZ)

                for lake in self.lake_layer:
                    lake_id = lake.GetField("lake_id")
                    logger.info("Lake %s covers %d tiles -> processed in LakeTile" % (lake_id, lake.GetField("overlap_avg_%s"%(in_tile_ref[-1]))))
                    self.big_lake_geom.append((lake_id, lake.GetGeometryRef().Clone()))

                self.lake_layer.SetAttributeFilter(None)
                self.lake_layer.ResetReading()

    # ----------------------------------------

    def open_db(self, in_lakedb_filename, in_table_name, in_field_name_list, in_field_type_list, 
                field_id=None, start_value=None, poly=None):
        """
        Open database, optionnally select features depending on their identifier, location or related continent
        and copy layer to memory

        :param in_lakedb_filename: full path of PLD
        :type in_lakedb_filename: string
        :param in_table_name: name of table to load from DB
        :type in_table_name: string
        :param in_field_name_list: list of fieldnames to load from table, first element is the identifier
        :type in_field_name_list: list of string
        :param in_field_type_list: list of type of each in_field_name_list
        :type in_field_type_list: list of string
        :param field_id: to select only the feature having this identifier
        :type field_id: string
        :param start_value: to select features having their identifier beginning with it
        :type start_value: string
        :param poly: polygon to spatially select features
        :type poly: ogr.Polygon
        
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
        out_layer = out_data_source.CreateLayer(str('layer'), srs=srs, geom_type=ogr.wkbMultiPolygon)
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

        # 2 - Select subset among PLD lakes using field_id, poly or having identifier starting with continent_id
        cmd = "SELECT %s, AsText(ST_Multi(geometry)) FROM %s " %(",".join(in_field_name_list), in_table_name)
        cmd_where = []
        if field_id is not None:
            cmd_where.append("%s LIKE '%s' " %(in_field_name_list[0], field_id))
        if start_value is not None:
            cmd_where.append(" {} LIKE '{}%' ".format(in_field_name_list[0], start_value))
        if poly is not None:
            poly.FlattenTo2D()  # Transform 3D geometry into 2D geometry (necessary for spatialite query)
            cmd_where.append(" ST_Intersects(GeomFromText('%s'), %s.geometry) " % (poly.ExportToWkt(), in_table_name))
        if cmd_where:
            cmd += " WHERE " + " AND ".join(cmd_where) + ";"
        else:
            cmd += ";"
        #logger.debug(cmd)
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

        logger.debug("%d features of table %s after focus over studied area" % (out_layer.GetFeatureCount(), in_table_name))
        return out_data_source, out_layer

    # ----------------------------------------

    def close_db(self):
        """
        Close database
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Close lake table
        logger.debug("- close lake table -")
        self.lake_ds.Destroy()
        
        # Close influence area table
        if self.influence_lake_flag == True:
            logger.debug("- close influence area table -")
            self.influence_lake_ds.Destroy()
        
        # Close basin table
        if self.basin_flag == True:
            logger.debug("- close basin table -")
            self.basin_ds.Destroy()
        

#######################################
        
        
def compute_tile_borders(in_tile_ref, in_tile_geom):
    """
    Compute bottom and top line (in the along-track direction) of input tile
    
    :param in_tile_ref: tile reference; format = <3-digit_tile_number><swath_side>
    :type in_tile_ref: string
    :param in_tile_geom: polygon of tile
    :type in_tile_geom: OGRPolygon
    
    :return: out_line_bottom = bottom line of input tile
    :rtype: out_line_bottom = OGRLineString
    :return: out_line_top = top line of input tile
    :rtype: out_line_top = OGRLineString
    """
    logger = logging.getLogger(__name__)

    # 1 - Retrieve ring inside OGRPolygon
    if in_tile_geom.GetGeometryName() == "POLYGON":
        ring = in_tile_geom.GetGeometryRef(0)
    elif in_tile_geom.GetGeometryName() == "MULTIPOLYGON":
        ring = in_tile_geom.GetGeometryRef(0).GetGeometryRef(0)
    else:
        logger.error("Unknown geometry %s " % in_tile_geom.GetGeometryName())

    # 2 - Define indices of wanted points in the polygon, depending on the tile swath side
    if in_tile_ref.endswith("L"):
        ind_bottom_1 = 0
        ind_bottom_2 = 3
        ind_top_1 = 1
        ind_top_2 = 2
    else:
        ind_bottom_1 = 0
        ind_bottom_2 = 1
        ind_top_1 = 2
        ind_top_2 = 3
        
    # 3 - Compute bottom line
    out_line_bottom = ogr.Geometry(ogr.wkbLineString)
    point1 = ring.GetPoint(ind_bottom_1)
    out_line_bottom.AddPoint(point1[0], point1[1])
    point2 = ring.GetPoint(ind_bottom_2)
    out_line_bottom.AddPoint(point2[0], point2[1])
    
    # 4 - Compute top line
    out_line_top = ogr.Geometry(ogr.wkbLineString)
    point1 = ring.GetPoint(ind_top_1)
    out_line_top.AddPoint(point1[0], point1[1])
    point2 = ring.GetPoint(ind_top_2)
    out_line_top.AddPoint(point2[0], point2[1])
    
    return out_line_bottom, out_line_top
        

#######################################


def compute_closest_polygon_with_kdtree(in_lon, in_lat, prior_geoms, prior_id):
    """
    Associate to each PIXCVec coordinate (in_lon, in_lat) the closest prior lake and its id

    :param in_lon: improved longitude of PixC related to in_poly
    :type in_lon: 1D array of float
    :param in_lat: improved latitude of PixC related to in_poly
    :type in_lat: 1D array of float
    :param prior_geoms: List of coordinates of polygons of prior lake database selected
    :type prior_geoms: 2D array of float
    :param prior_id: List of prior ID from lake DB
    :type prior_id: 1D array of str

    :return: list of the closest prior_id associated to the (in_lon, in_lat) points
    :rtype: list of str
    """
    lon_coords, lat_coords, prior_id_list = np.array([]), np.array([]), np.array([], dtype=object)
    for geom, id in zip(prior_geoms, prior_id):

        # Get the coordinates with json library
        json_geom = json.loads(geom.ExportToJson())

        if json_geom["type"] == 'Polygon':
            coords = list(json_geom['coordinates'])
            coords_array = np.array(coords).ravel()
        elif json_geom["type"] == 'MultiPolygon':
            coords_array = np.array([poly for poly in json_geom['coordinates']]).ravel()
            while coords_array.dtype == 'object':  # array of list : depending if multipolygon contains polygons with holes
                new_coords = np.array([])
                for x in coords_array:
                    new_coords = np.append(new_coords, np.array(x).ravel())
                coords_array = np.array(new_coords).ravel()
        else:
            coords_array = np.array([])

        lon_coords = np.append(lon_coords, coords_array[0::2])
        lat_coords = np.append(lat_coords, coords_array[1::2])
        prior_id_list = np.append(prior_id_list, [id] * coords_array[1::2].size)

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


def compute_continent_id_from_basin_code(in_basin_code):
    """
    Compute continent identifier from basin code (HydroBASINS nomenclature)
    * 1 for Africa
    * 2 for Europe
    * 3 for Siberia
    * 4 for Asia
    * 5 for Australia
    * 6 for South America
    * 7 for North America
    * 8 for Arctic (North America)
    * 9 for Greenland
    Added for technical purposes:
    * 0 for Ocean
    * 10 for no continent

    :param in_basin_code: basin code
    :type in_basin_code: string

    :return: 2-letter identifier of the continent
    :rtype: string
    """

    retour = "None"

    if in_basin_code.startswith("1"):
        retour = "AF"
    elif in_basin_code.startswith("2"):
        retour = "EU"
    elif in_basin_code.startswith("3"):
        retour = "SI"
    elif in_basin_code.startswith("4"):
        retour = "AS"
    elif in_basin_code.startswith("5"):
        retour = "AU"
    elif in_basin_code.startswith("6"):
        retour = "SA"
    elif in_basin_code.startswith("7"):
        retour = "NA"
    elif in_basin_code.startswith("8"):
        retour = "AR"
    elif in_basin_code.startswith("9"):
        retour = "GR"

    return retour


def compute_continent_code(in_continent_id):
    """
    Compute continent code from continent short name (HydroBASINS nomenclature)
    * 1 for Africa
    * 2 for Europe
    * 3 for Siberia
    * 4 for Asia
    * 5 for Australia
    * 6 for South America
    * 7 for North America
    * 8 for Arctic (North America)
    * 9 for Greenland
    Added for technical purposes:
    * 0 for Ocean
    * 10 for no continent

    :param in_continent_id: 2-letter continent identifier
    :type in_continent_id: string

    :return: 1-digit continent code
    :rtype: string
    """

    retour = "None"
    
    if in_continent_id == "AF":
        retour = "1"
    elif in_continent_id == "EU":
        retour = "2"
    elif in_continent_id == "SI":
        retour = "3"
    elif in_continent_id == "AS":
        retour = "4"
    elif in_continent_id == "AU":
        retour = "5"
    elif in_continent_id == "SA":
        retour = "6"
    elif in_continent_id == "NA":
        retour = "7"
    elif in_continent_id == "AR":
        retour = "8"
    elif in_continent_id == "GR":
        retour = "9"

    return retour
    
    
#######################################
    

class PriorLake(object):
    """
    This class manages one lake from the Lake A priori Database (PLD)
    """

    def __init__(self, in_obj_lakedb, in_lakeid):
        """
        Constructor: set the general values

        Variables of the object:
        - lake_id / str: identifier of the prior lake
        - geom / OGRPolygon: polygon of the prior lake
        - name / str: name of the prior lake
        - res_id / int: reservoir identifier if the lake is a reservoir identified in a database (GRanD for the moment)
        - lon / float: longitude of the centroid of the lake
        - lat / float: latitude of the centroid of the lake
        - max_wse / float: maximum water surface elevation used for storage change computation (in m)
        - max_wse_u / float: uncertainty over maximum water surface elevation used for storage change computation (in m)
        - max_area / float: maximum area used for storage change computation (in km2)
        - max_area_u / float: uncertainty over maximum area used for storage change computation (in km2)
        - ref_date / str: date of reference storage change value, used as 0
        - ref_ds / str: reference storage change value, used as 0 (in km3)
        - storage / float: highest measured storage change value
        - ice_clim / int: ice clim flag
        - ice_dyn / int: ice dyn flag
        - reach_id / str: list of Prior Reach Database reach_id intersecting the lake
        - ok_to_compute_stocc / boolean: =True if storage change is computable, =False otherwise
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # 0 - Init variables
        self.lake_id = in_lakeid
        self.geom = None
        self.name = None
        self.res_id = None
        self.lon = None
        self.lat = None
        self.max_wse = None
        self.max_wse_u = None
        self.max_area = None
        self.max_area_u = None
        self.ref_date = None
        self.ref_ds = None
        self.storage = None
        self.ice_clim = None
        self.ice_dyn = None
        self.reach_id = None

        if in_obj_lakedb.lake_layer: # In case a PLD is used
            
            # Init dictionary of available PLD infos
            dict_pld_info = {}
            for item in in_obj_lakedb.pld_infos:
                dict_pld_info[item] = None
            
            # 1 - Select feature given its identifier
            if in_obj_lakedb.lakedb_id_type == "String":
                in_obj_lakedb.lake_layer.SetAttributeFilter("%s like '%s'" % (in_obj_lakedb.lakedb_id_name, str(in_lakeid)))
            else:
                in_obj_lakedb.lake_layer.SetAttributeFilter("%s = %s" % (in_obj_lakedb.lakedb_id_name, str(in_lakeid)))
            pld_lake_feat = in_obj_lakedb.lake_layer.GetNextFeature()
            
            # 2.1 - Retrieve geometry
            self.geom = pld_lake_feat.GetGeometryRef().Clone()
            # 2.2 - Retrieve information when exists
            for item in in_obj_lakedb.pld_infos:
                dict_pld_info[item] = pld_lake_feat.GetField(item)

            # 3 - Release filter
            in_obj_lakedb.lake_layer.SetAttributeFilter(None)
            
            # 4 - Format output
            self.check_value(dict_pld_info, in_obj_lakedb)
 
        # 5 - Set flag to compute storage change if reference data are available
        self.ok_to_compute_stocc = False
        if (self.max_wse is not None) and (self.max_area is not None):
            self.ok_to_compute_stocc = True

    # ----------------------------------------

    def check_value(self, dict_pld_info, in_obj_lakedb):
        """
        Check PLD value

        :param dict_pld_info: dictionary of available PLD infos
        :type dict_pld_info: dictionary
        :param in_obj_lakedb: lakedb object
        :type in_obj_lakedb: lakedb object
        """
        # 4.1 - List of names
        self.name = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_names)
        if (self.name is not None) and (self.name in ["", my_var.FV_STRING_SHP]):
            self.name = None
        # 4.2 - GRanD identifier
        self.res_id = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_res_id)
        if (self.res_id is not None) and (self.res_id < 0):
            self.res_id = None
        # 4.3 - Longitude of the centroid
        self.lon = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_lon)
        if (self.lon is not None) and (self.lon < -180.0):
            self.lon = None
        # 4.4 - Latitude of the centroid
        self.lat = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_lat)
        if (self.lat is not None) and (self.lat < -90.0):
            self.lat = None
        # 4.5 - Max water surface elevation
        self.max_wse = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_max_wse)
        if (self.max_wse is not None) and (self.max_wse < 0):
            self.max_wse = None
        # 4.6 - Uncertainty over max water surface elevation
        self.max_wse_u = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_max_wse_u)
        if (self.max_wse_u is not None) and (self.max_wse_u < 0):
            self.max_wse_u = None
        # 4.7 - Max area
        self.max_area = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_max_area)
        if (self.max_area is not None) and (self.max_area < 0):
            self.max_area = None
        # 4.8 - Uncertainty over max area
        self.max_area_u = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_max_area_u)
        if (self.max_area_u is not None) and (self.max_area_u < 0):
            self.max_area_u = None
        # 4.9 - Reference date
        self.ref_date = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_ref_date)
        if (self.ref_date is not None) and (self.ref_date in ["", my_var.FV_STRING_SHP]):
            self.ref_date = None
        # 4.10 - Reference data storage
        self.ref_ds = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_ref_ds)
        if (self.ref_ds is not None) and (self.ref_ds < -9e7):
            self.ref_ds = None
        # 4.11 - Absolute water storage
        self.storage = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_storage)
        if (self.storage is not None) and (self.storage < 0):
            self.storage = None
        # 4.12 - List of PRD intersecting reach_id
        self.reach_id = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_reach)
        if (self.reach_id is not None) and (self.reach_id in ["", my_var.FV_STRING_SHP]):
            self.reach_id = None

        # 4.13 - Average overlap of lake and tile of current pass and swath L
        self.overlap_avg_l = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_overlap_avg_l)
        if (self.overlap_avg_l is not None) and (self.overlap_avg_l < 0):
            self.overlap_avg_l = None

        # 4.14 - Number of tiles covered by a lake over swath L
        self.nb_tiles_l = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_nb_tiles_l)
        if (self.nb_tiles_l is not None) and (self.nb_tiles_l < 0):
            self.nb_tiles_l = None

        # 4.15 - Average overlap of lake and tile of current pass and swath R
        self.overlap_avg_r = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_overlap_avg_r)
        if (self.overlap_avg_r is not None) and (self.overlap_avg_r < 0):
            self.overlap_avg_r = None

        # 4.16 - Number of tiles covered by a lake over swath R
        self.nb_tiles_r = my_tools.get_value(dict_pld_info, in_obj_lakedb.pld_nb_tiles_r)
        if (self.nb_tiles_r is not None) and (self.nb_tiles_r < 0):
            self.nb_tiles_r = None
      
    def format_attributes(self):
        """
        Format PLD attributes to the format expected in the LakeProduct
        i.e. convert None to correct _FillValue
        
        :return: out_attributes = attributes correctly formatted
        :rtype: out_attributes = dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # 0 - Init output dictionary
        out_attributes = dict()
        
        # 1 - Set lake_id
        out_attributes["lake_id"] = self.lake_id
        
        # 2 - Compute values to write
        
        # 2.1 - Set PLD name
        if self.name is not None:
            out_attributes["lake_name"] = self.name
        else:
            out_attributes["lake_name"] = my_var.FV_STRING_SHP
                
        # 2.2 - Set reservoir identifier
        if self.res_id is not None:
            out_attributes["p_res_id"] = self.res_id
        else:
            out_attributes["p_res_id"] = my_var.FV_INT9_SHP
                
        # 2.3 - Set longitude of centroid
        if self.lon is not None:
            out_attributes["p_lon"] = self.lon
        else:
            out_attributes["p_lon"] = my_var.FV_REAL
                
        # 2.4 - Set latitude of centroid
        if self.lat is not None:
            out_attributes["p_lat"] = self.lat
        else:
            out_attributes["p_lat"] = my_var.FV_REAL

        # 2.5 - Set ref water surface elevation
        if self.max_wse is not None:
            out_attributes["p_ref_wse"] = self.max_wse
        else:
            out_attributes["p_ref_wse"] = my_var.FV_REAL

        # 2.6 - Set ref area
        if self.max_area is not None:
            out_attributes["p_ref_area"] = self.max_area
        else:
            out_attributes["p_ref_area"] = my_var.FV_REAL

        # 2.7 - Set date t0
        if self.ref_date is not None:
            out_attributes["p_date_t0"] = self.ref_date
        else:
            out_attributes["p_date_t0"] = my_var.FV_STRING_SHP

        # 2.8 - Set storage change at t0
        if self.ref_ds is not None:
            out_attributes["p_ds_t0"] = self.ref_ds
        else:
            out_attributes["p_ds_t0"] = my_var.FV_REAL
                
        # 2.9 - Update maximum water storage value
        if self.storage is not None:
            out_attributes["p_storage"] = self.storage
        else:
            out_attributes["p_storage"] = my_var.FV_REAL
                
        # 2.10 - Update ice clim flag value
        if self.ice_clim is not None:
            out_attributes["ice_clim_f"] = self.ice_clim
        else:
            out_attributes["ice_clim_f"] = my_var.FV_INT_SHP
                
        # 2.11 - Update ice dyn flag value
        if self.ice_dyn is not None:
            out_attributes["ice_dyn_f"] = self.ice_dyn
        else:
            out_attributes["ice_dyn_f"] = my_var.FV_INT_SHP
        
        # 2.12 - Set reach identifier
        if self.reach_id is not None:
            out_attributes["reach_id"] = self.reach_id
        else:
            out_attributes["reach_id"] = my_var.FV_STRING_SHP
            
        return out_attributes

    # ----------------------------------------
    
    def run_stocc(self, in_list_obs):
        """
        Call the different functions of storage change computation
        
        :param in_list_obs: list of observed features intersecting the prior lake
                            - in_list_obs[obs_id]["area"] = area of observed lake obs_id
                            - in_list_obs[obs_id]["area_u"] = uncertainty over area of observed lake obs_id
                            - in_list_obs[obs_id]["wse"] = water surface elevation of observed lake obs_id
                            - in_list_obs[obs_id]["wse_u"] = uncertainty over water surface elevation of observed lake obs_id
                            - in_list_obs[obs_id]["alpha"] = proprotionnal coefficient related to area of observed lake wrt
                                                                all observed lakes linked to the prior lake
        :type in_list_obs: dict
        
        :return: out_delta_s_l = linear storage change value
        :rtype: out_delta_s_l = float
        :return: out_ds_l_u = linear storage change error
        :rtype: out_ds_l_u = float
        :return: out_delta_s_q = quadratic storage change value
        :rtype: out_delta_s_q = float
        :return: out_ds_q_u = qradratic storage change error
        :rtype: out_ds_q_u = float
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # 1 - Linear storage change
        # 1.1 - Compute raw linear storage change
        stoc_val, stoc_u = storage_change.stocc_linear_basic(in_list_obs,
                                                             self.max_area,
                                                             self.max_area_u,
                                                             self.max_wse, 
                                                             self.max_wse_u)
        # 1.2 - Deduce linear storage change value
        out_delta_s_l = my_var.FV_REAL
        if stoc_val is not None:
            if self.ref_ds is not None:
                out_delta_s_l = stoc_val - self.ref_ds
            else:
                out_delta_s_l = stoc_val
        # 1.3 - Deduce linear storage change error
        out_ds_l_u = my_var.FV_REAL
        if stoc_u is not None:
            out_ds_l_u = stoc_u
        
        # 2 - Quadratic storage change
        # 2.1 - Compute raw linear storage change
        stoc_val, stoc_u = storage_change.stocc_quadratic_basic(in_list_obs, 
                                                                 self.max_area,
                                                                 self.max_area_u,
                                                                 self.max_wse, 
                                                                 self.max_wse_u)
        # 2.2 - Deduce quadratic storage change value
        out_delta_s_q = my_var.FV_REAL
        if stoc_val is not None:
            if self.ref_ds is not None:
                out_delta_s_q = stoc_val - self.ref_ds
            else:
                out_delta_s_q = stoc_val
        # 2.3 - Deduce quadratic storage change error
        out_ds_q_u = my_var.FV_REAL
        if stoc_u is not None:
            out_ds_q_u = stoc_u
        
        return out_delta_s_l, out_ds_l_u, out_delta_s_q, out_ds_q_u
    
