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
# VERSION:5.0.0:DM:#91:2022/12/08:Poursuite industrialisation
# VERSION:5.0.1:DM:#91:2023/02/24:Poursuite industrialisation
# VERSION:5.2.0:DM:#91:2023/05/23:Poursuite industrialisation
# VERSION:5.2.1:DM:#91:2023/07/24:Poursuite industrialisation
# VERSION:5.3.0:DM:#91:2023/10/04:Poursuite industrialisation
# VERSION:6.1.0:DM:#91:2024/09/20:Poursuite industrialisation
# VERSION:6.2.0:DM:#91:2024/11/29:Poursuite industrialisation
# VERSION:6.2.1:FA:#823:2025/04/08:Nouveau patch Lake_AVG_V4.2.1
# VERSION:6.2.4:FA:#854:2026/02/18:#SSALTO-446 [PGE][LAKE_SP] doublon de SUBSTR dans la requête INSERT INTO
# VERSION:6.2.4:FA:#831:2026/02/26:SSALTO-330 [PGE][SAS LAKES] Mise à jour des TU (SSALPGE-14313)
# VERSION:6.2.6:DM:#900:2026/07/06: [PGE] [LAKES] Ajout d'une interface niveau moyen des lacs
# VERSION:6.2.6:DM:#900:2026/07/06:SSALTO-565 [PGE] [LAKES] Ajout d'une interface coeff NN pour calcul ice_dyn_f
# VERSION:6.2.7:DM:#906:2026/07/28:SSALTO-872 [PGE][LAKES] Gestion des erreurs/warning GDAL
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
        - filename / str: database file full path
        - db_format / str: database format (=None if no DB used; ="shp" for shapefile format; ="sqlite" for SQLite format)
        - lakedb_id_name / str: fieldname of identifier (=primary key) of lake table
        - lakedb_id_type / str: type of identifier (=primary key) in lake table
        - basindb_id_name / str: fieldname of identifier (=primary key) in basin table
        - pld_field_lake_verified / dict: all elements of my_var.PLD_FIELD_LAKE truly in the input PLD
        - big_lake_geom / list: list of pairs (lake_id, polygon) of PLD lakes considered too big to be process in LakeSP
        - lake_layer / osgeo.ogr.Layer: layer of lake table in PLD
        - lake_ds / osgeo.ogr.DataSource: associated DataSource
        - flag_ocean / boolean: (case of LakeTile processing) =True if tile is entirely over ocean
        - lake_assignment_layer_name / str : name of layer of lake assignment table in PLD (lake_influence or TopoCat.lake_catchment)
        - lake_assignment_flag / boolean: flag indicating if the assignment lake geometries are used or not
        - lake_assignment_layer / osgeo.ogr.Layer: layer of lake assignment table in PLD (lake_influence or TopoCat.lake_catchment)
        - lake_assignment_ds / osgeo.ogr.DataSource: associated DataSource
        - basin_flag / boolean: flag indicating if the basin geometries are used or not
        - basin_layer / osgeo.ogr.Layer: layer of basin table in PLD
        - basin_ds / osgeo.ogr.DataSource: associated DataSource
        - db_connector / sqlite3.Connection: connector to a memory SQLite database (when operationnal SQLite format; case of LakeSP and LakeAvg processing)
        - db_cursor / sqlite3.Cursor: cursor to the memory SQLite database (when operationnal SQLite format; case of LakeSP and LakeAvg processing)
        - flag_prior_init / boolean: =True if the LakeSP_Prior or LakeAvg layer has been init with prior lakes, =False otherwise
        - az_0_geom / osgeo.ogr.Geometry: geometries of lakes located at the bottom of tile
        - az_max_geom / osgeo.ogr.Geometry: geometries of lakes located at the top of tile
        - az_0_and_max_geom / osgeo.ogr.Geometry: geometries of lakes located at the top and bottom of tile
        """
        try:
            # Get instance of service config file
            cfg = service_config_file.get_instance()
            # Save min_overlap parameter
            self.min_overlap = cfg.getfloat('CONFIG_PARAMS', 'MIN_OVERLAP')
            # Layer name to use for lake assignment
            self.lake_assignment_layer_name = cfg.get("CONFIG_PARAMS", "ASSIGNMENT_LAYER")
        except:
            self.min_overlap = 2.0
            self.lake_assignment_layer_name = "lake_catchment"
        
        # DB filename and type
        self.filename = None
        self.db_format = None

        # Fieldnames
        self.lakedb_id_name = my_var.PLD_FIELD_LAKE["id"]  # Fieldname for lake identifier
        self.lakedb_id_type = ""  # Type of lake identifier
        self.basindb_id_name = my_var.PLD_FIELD_BASIN["id"]  # Fieldname for basin identifier
        self.pld_field_lake_verified = dict()  # All elements of my_var.PLD_FIELD_LAKE truly in the input PLD
        
        # List of (lake_id, related geometries) of PLD lakes in the studied tile
        self.big_lake_geom = []

        # Lake table layer and DataSource
        self.lake_layer = None
        self.lake_ds = None
        self.flag_ocean = False

        # Lake assignment table layer and DataSource
        self.lake_assignment_flag = False  # Use of lake assignment geometries (default = False)
        self.lake_assignment_layer = None
        self.lake_assignment_ds = None

        # Basin table layer and dataSource
        self.basin_flag = False  # Use of basin geometries (default = False)
        self.basin_layer = None
        self.basin_ds = None
        
        # Connector and cursor to a memory SQLite database, if used
        self.db_connector = None
        self.db_cursor = None
        self.flag_prior_init = False  # Indicates if the LakeSP_Prior or LakeAvg layer has been init with prior lakes 

        # Union of geometries of PLD lakes located at the top and/or bottom of tile
        self.az_0_geom = None
        self.az_max_geom = None
        self.az_0_and_max_geom = None

    def close_db(self):
        """
        Close database
        """
        pass

    # ----------------------------------------        
        
    def init_prior_layer(self, in_obj_lake, in_date=None):
        """
        Init layer of input lake object with prior lakes over the granule area
        
        :param in_obj_lake: lake object
        :type in_obj_lake: LakeSPPriorProduct or LakeAvgProduct
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("WARNING layer CANNOT be init from lake database because there is none")

    # ----------------------------------------
    
    def init_fields_name_and_type(self):
        """
        Init type of lake identifier and eventually test existence of needed fieldnames
        """
        
        # 1 - Retrieve list of fieldnames and their type available in input PLD
        dict_field_type = my_tools.get_layer_fields_name_and_type(self.lake_layer)
        
        # 2 - Save type of identifier
        self.lakedb_id_type = dict_field_type[self.lakedb_id_name]
        
        # 3 - Compare to a priori fieldnames
        for key, pld_fieldname in my_var.PLD_FIELD_LAKE.items():
            if pld_fieldname in dict_field_type.keys():
                self.pld_field_lake_verified[key] = pld_fieldname

    # ----------------------------------------

    def build_border_geometry(self, in_az_0_line, in_az_max_line, tile_geom=None):
        """
        Compute Multipolygon from all polygons intersecting in_az_0_line or in_az_max_line and a buffer of in_az_[0-max]_line. 
        Border polygons are stored in self.az_0_geom, self.az_max_geom and self.az_0_and_max_geom.
        If tile_geom is given, border geometries are cut by tile_geom.

        :param in_az_0_line: Line with azimuth = 0
        :type: OGR Geometry Linestring
        :param in_az_max_line: Line with azimuth = azimuth max
        :type: OGR Geometry Linestring
        :param tile_geom: current tile geometry
        :type: OGR Geometry MultiPolygon
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        list_border_lake = set()
        
        az_0_line_buffer = in_az_0_line.Buffer(0.000001)
        az_max_line_buffer = in_az_max_line.Buffer(0.000001)

        if in_az_0_line:
            self.lake_layer.ResetReading()
            self.lake_layer.SetSpatialFilter(in_az_0_line)

            for feat in self.lake_layer:
                lake_id = feat.GetField(self.lakedb_id_name)
                geom = feat.GetGeometryRef()

                # If tile_geom is given, PLD.lake.geometry is cut by tile to deal with big geometries.
                if tile_geom:
                    try:
                        geom = geom.Intersection(tile_geom)
                    except Exception as e:
                        logger.warning(f"PLD lake contains an invalid geometry at lake_id : {lake_id}, causing the following error : {e}")
                        geom = geom.Buffer(0).Intersection(tile_geom)

                if geom.Intersects(az_max_line_buffer):
                    if self.az_0_and_max_geom is None:
                        self.az_0_and_max_geom = geom.Clone()
                    else:
                        self.az_0_and_max_geom = self.az_0_and_max_geom.Union(geom)
                    logger.debug("PLD lake_id {} = az_0_and_max_geom".format(lake_id))
                else :
                    if self.az_0_geom is None:
                        self.az_0_geom = geom.Clone()
                    else:
                        self.az_0_geom = self.az_0_geom.Union(geom)
                    logger.debug("PLD lake_id {} = az_0_geom".format(lake_id))

                list_border_lake.add(lake_id)

            self.lake_layer.SetSpatialFilter(None)
            self.lake_layer.ResetReading()

        if in_az_max_line:
            self.lake_layer.ResetReading()
            self.lake_layer.SetSpatialFilter(in_az_max_line)

            for feat in self.lake_layer:
                lake_id = feat.GetField(self.lakedb_id_name)
                geom = feat.GetGeometryRef()

                # If tile_geom is given, PLD.lake.geometry is cut by tile to deal with big geometries.
                if tile_geom:
                    try:
                        geom = geom.Intersection(tile_geom)
                    except Exception as e:
                        logger.warning(f"PLD lake contains an invalid geometry at lake_id : {lake_id}, causing the following error : {e}")
                        geom = geom.Buffer(0).Intersection(tile_geom)

                if geom.Intersects(az_0_line_buffer):
                    if lake_id not in list_border_lake:
                        if self.az_0_and_max_geom is None:
                            self.az_0_and_max_geom = geom.Clone()
                        else:
                            self.az_0_and_max_geom = self.az_0_and_max_geom.Union(geom)
                        logger.debug("PLD lake_id {} = az_0_and_max_geom".format(lake_id))
                else :
                    if self.az_max_geom is None:
                        self.az_max_geom = geom.Clone()
                    else:
                        self.az_max_geom = self.az_max_geom.Union(geom)
                    logger.debug("PLD lake_id {} = az_max_geom".format(lake_id))
                list_border_lake.add(lake_id)

            self.lake_layer.SetSpatialFilter(None)
            self.lake_layer.ResetReading()
            
        # Print
        logger.debug("%s PLD lakes will be processed in LakeSP because they intersect tile border" % len(list_border_lake))
        logger.debug(";".join(list_border_lake))

    # ----------------------------------------

    def lake_catchment_intersection(self, in_poly, in_lon, in_lat, obs_area, pld_inter_dict, sum_area_intersections):
        """
        Attribute pixels to pld lake depending on the lake catchment feature.

        :param in_poly: polygon delineating a water body
        :type in_poly: OGRPolygon
        :param in_lon: improved longitude of PixC related to in_poly
        :type in_lon: 1D array of float
        :param in_lat: improved latitude of PixC related to in_poly
        :type in_lat: 1D array of float
        :param obs_area: lake observed area
        :param obs_area: float
        :param pld_inter_dict: dict containing lake_id, intersection area and pld geometry for each PLD lake intersecting observed lake
        :type pld_inter_dict: dict
        :param sum_area_intersections: sum of intersections areas
        :type sum_area_intersections: float

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

        # 0 - Init variables
        tmp_list_prior_id = []  # List of PLD identifiers
        tmp_list_pld_overlap = []  # List of fractions of observed lake covered by PLD lakes
        prior_geoms = []

        # 1 - Compute relative intersection area wrt observed area minus intersections area
        for lake_id, (inter_area, pld_geom) in pld_inter_dict.items():

            # 1.1 - Compute associated fraction of observed lake covered by PLD lake
            frac_inter = round( inter_area / obs_area * 100.)

            # 1.2 - Save info to temporary lists
            tmp_list_prior_id.append(lake_id)  # PLD identifier
            tmp_list_pld_overlap.append(str(frac_inter))  # Overlap
            prior_geoms.append(pld_geom)
            
        # 2 - Compute out_list_prior_id, out_list_pld_overlap and out_pixcvec_lakeid
        out_list_prior_id = None
        out_list_pld_overlap = None
        out_pixcvec_lakeid = None

        if len(tmp_list_prior_id) != 0:  

            # 2.1 - Compute PIXCVec_tag

            # 2.1.0 - Init out_lakedb_id_pixcvec to empty string
            out_pixcvec_lakeid = np.zeros(in_lat.size, dtype=object)
            out_pixcvec_lakeid[:] = ""

            # 2.1.1 - Pixels within PLD.lake polygons are directly affected to its corresponding lakes
            out_pixcvec_lakeid = self.affect_pixels_within_lake_polygons(in_lon, in_lat, tmp_list_prior_id, out_pixcvec_lakeid)

            # 2.1.2 - Pixels outside of PLD.lake polygons are affected to its corresponding lakes using functions
            out_pixcvec_lakeid = self.affect_pixels_outside_lake_polygons(in_lon, in_lat, tmp_list_prior_id, prior_geoms, out_pixcvec_lakeid)

            # 2.2 - Sort output prior ID and overlap fractions lists by decreasing area intersection
            sorted_idx = np.argsort(tmp_list_pld_overlap)[::-1]
            out_list_prior_id = [tmp_list_prior_id[idx] for idx in sorted_idx]
            out_list_pld_overlap = [str(tmp_list_pld_overlap[idx]) for idx in sorted_idx]

            # Print number of pixels and lake_id
            unique, counts = np.unique(out_pixcvec_lakeid, return_counts=True)
            for ind, unique_val in enumerate(unique):
                logger.debug("%d pixels of current observed lake belong to lake_id %s" % (counts[ind], unique_val))

        return out_list_prior_id, out_list_pld_overlap, out_pixcvec_lakeid

    def lake_influence_intersection(self, in_poly, in_lon, in_lat, obs_area, pld_inter_dict, sum_area_intersections):
        """
        Attribute pixels to pld lake depending on the lake influence feature and overlap percentage.

        :param in_poly: polygon delineating a water body
        :type in_poly: OGRPolygon
        :param in_lon: improved longitude of PixC related to in_poly
        :type in_lon: 1D array of float
        :param in_lat: improved latitude of PixC related to in_poly
        :type in_lat: 1D array of float
        :param obs_area: lake observed area
        :param obs_area: float
        :param pld_inter_dict: dict containing lake_id, intersection area and pld geometry for each PLD lake intersecting observed lake
        :type pld_inter_dict: dict
        :param sum_area_intersections: sum of intersections areas
        :type sum_area_intersections: float

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

        # 0 - Init variables
        tmp_list_prior_id = []  # List of PLD identifiers
        tmp_list_pld_overlap = []  # List of fractions of observed lake covered by PLD lakes
        prior_geoms = []

        # 1 - Compute relative intersection area wrt observed area minus intersections area
        for lake_id, (inter_area, pld_geom) in pld_inter_dict.items():

            # 1.1 - Compute observed area that to not intersects other pld lakes
            unlinked_area = obs_area - (sum_area_intersections - inter_area)

            # 1.2 - Compute ratio of current intersections area and unassociated area.
            inter_ratio = round(inter_area / unlinked_area * 100.)

            # 1.3 - Save info to output lists
            # only if overlap percentage is important enough
            if inter_ratio > self.min_overlap:
                # Compute associated fraction of observed lake covered by PLD lake
                frac_inter = round( inter_area / obs_area * 100.)
                tmp_list_prior_id.append(lake_id)  # PLD identifier
                tmp_list_pld_overlap.append(str(frac_inter))  # Overlap
                prior_geoms.append(pld_geom)
            else:
                logger.debug("Overlap with lake {} = {}% (<{}%) => PLD lake not linked with observed lake".format(lake_id, inter_ratio,self.min_overlap))

        # 2 - Compute out_list_prior_id, out_list_pld_overlap and out_pixcvec_lakeid
        out_list_prior_id = None
        out_list_pld_overlap = None
        out_pixcvec_lakeid = None
        
        # 2.1 - Easy match: polygon strongly matches only one PLD lake
        if len(tmp_list_prior_id) == 1: 
            out_list_prior_id = tmp_list_prior_id  # PLD identifier
            out_list_pld_overlap = tmp_list_pld_overlap  # Overlap
            out_pixcvec_lakeid[:] = tmp_list_prior_id[0]  # Set lake_id attributes in PIXCVec product
        
        # 2.2 - In the case where in_poly intersects several lakes from PLD and more than one lake have an overlap bigger than 2%
        if len(tmp_list_prior_id) > 1:  

            # 2.2.0 - Init out_lakedb_id_pixcvec to empty string
            out_pixcvec_lakeid = np.zeros(in_lat.size, dtype=object)
            out_pixcvec_lakeid[:] = ""

            # 2.2.1 - Pixels within PLD.lake polygons are directly affected to its corresponding lakes
            out_pixcvec_lakeid = self.affect_pixels_within_lake_polygons(in_lon, in_lat, tmp_list_prior_id, out_pixcvec_lakeid)

            # 2.2.2 - Pixels outside of PLD.lake polygons are affected to its corresponding lakes using functions
            out_pixcvec_lakeid = self.affect_pixels_outside_lake_polygons(in_lon, in_lat, tmp_list_prior_id, prior_geoms, out_pixcvec_lakeid)

            # 2.2.3 - Sort output prior ID and overlap fractions lists by decreasing area intersection
            sorted_idx = np.argsort(tmp_list_pld_overlap)[::-1]
            out_list_prior_id = [tmp_list_prior_id[idx] for idx in sorted_idx]
            out_list_pld_overlap = [str(tmp_list_pld_overlap[idx]) for idx in sorted_idx]

            # Print number of pixels and lake_id
            unique, counts = np.unique(out_pixcvec_lakeid, return_counts=True)
            for ind, unique_val in enumerate(unique):
                logger.debug("%d pixels of current observed lake belong to lake_id %s" % (counts[ind], unique_val))

        return out_list_prior_id, out_list_pld_overlap, out_pixcvec_lakeid

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

                # Init dict to store lake_id, intersection area and pld geometry for each PLD lake intersecting observed lake
                # key : lake_id, value : tuple (inter_area, pld_geometry)
                pld_inter_dict = {}

                # 3 - Check PLD lake intersecting in_poly and compute intersections
                for cur_lake_bd in self.lake_layer:

                    # 3.1 - Retrieve PLD lake identifier and geometry
                    cur_id = str(cur_lake_bd.GetField(self.lakedb_id_name))
                    #logger.debug("Associated PLD identifier = %s" % cur_id)
                    cur_geom = cur_lake_bd.GetGeometryRef().Clone()

                    if cur_id:  # Tested but opposite should not occur...
                        # 3.2 - Compute intersection of observed lake and PLD lake
                        # Compute the geometry of the intersection between the observed lake and the PLD lake
                        try:
                            geom_inter = in_poly.Intersection(cur_geom)
                        except Exception as e:
                            logger.warning(f"PLD lake contains an invalid geometry at lake_id : {cur_id}, causing the following error : {e}")
                            geom_inter = in_poly.Intersection(cur_geom.Buffer(0))

                        if geom_inter is not None:
                            pld_inter_dict[cur_id] = (my_tools.get_area(geom_inter), cur_geom)
                        else:
                            logger.warning("Something wrong happened with PLD lake %s: in fact, it doesn't intersect the current observed lake!" %\
                                           cur_id)

                # Compute sum of intersections areas
                sum_area_intersections = sum([inter_area[0] for inter_area in pld_inter_dict.values()])
                # Compute observed lake area
                obs_area = my_tools.get_area(in_poly)

                if self.lake_assignment_layer_name == "lake_catchment": 
                    out_list_prior_id, out_list_pld_overlap, out_pixcvec_lakeid = self.lake_catchment_intersection( 
                    in_poly, in_lon, in_lat, obs_area, pld_inter_dict, sum_area_intersections)
                else:
                    out_list_prior_id, out_list_pld_overlap, out_pixcvec_lakeid = self.lake_influence_intersection( 
                    in_poly, in_lon, in_lat, obs_area, pld_inter_dict, sum_area_intersections)

            self.lake_layer.SetSpatialFilter(None)  # Delete spatial filter

        return out_list_prior_id, out_list_pld_overlap, out_pixcvec_lakeid

    def affect_pixels_within_lake_polygons(self, in_lon, in_lat, in_list_prior_id, out_pixcvec_lakeid):
        """
        Affect pixels from PIXCVec located within lakes from in_list_prior_id

        :param in_lon: longitudes of pixels
        :type in_lon: 1D-array of float
        :param in_lat: latitudes of pixels
        :type in_lat: 1D-array of float
        :param in_list_prior_id: prior ID of overlapping PLD lakes
        :type in_list_prior_id: numpy array of str
        :param: out_pixcvec_lakeid = list of lake_id associated to each PIXC constituting observed lake. Unaffected pixels are empty str
        :type out_pixcvec_lakeid: numpy array of str

        :return: out_pixcvec_lakeid = list of lake_id associated to each PIXC constituting observed lake. Unaffected pixels are empty str
        :rtype: numpy array of str
        """
        logger = logging.getLogger(self.__class__.__name__)
        # 1. Create memory layer to store pixcvec points
        ds, lyr = my_tools.load_pixels_to_mem_layer(in_lon, in_lat)

        # 2. Filter PLD.lake with lake_id in in_list_prior_id
        self.lake_layer.ResetReading()
        request = "%s IN ('%s'" % (self.lakedb_id_name, in_list_prior_id[0])
        for lakedb_id in in_list_prior_id[1:]:
            request += ", '%s'" % (lakedb_id)
        request += ")"
        self.lake_layer.SetAttributeFilter(request)

        # 3. Filter point of pixcvec with PLD.lake geometry
        for lake_feat in self.lake_layer:
            lyr.SetSpatialFilter(lake_feat.GetGeometryRef())
            for feat_point in lyr:
                out_pixcvec_lakeid[feat_point.GetFID()] = lake_feat.GetField(self.lakedb_id_name)
            lyr.SetSpatialFilter(None)
        ds.Destroy()

        logger.debug("%d / %d pixels within PLD.lake polygons" % (np.sum(out_pixcvec_lakeid != ""), out_pixcvec_lakeid.size))
        self.lake_layer.SetAttributeFilter(None)
        return out_pixcvec_lakeid

    def affect_pixels_outside_lake_polygons(self, in_lon, in_lat, in_list_prior_id, in_prior_geoms, out_pixcvec_lakeid):
        """
        Affect pixels from PIXCVec located outside lakes from in_list_prior_id.

        :param in_lon: longitudes of pixels
        :type in_lon: 1D-array of float
        :param in_lat: latitudes of pixels
        :type in_lat: 1D-array of float
        :param in_list_prior_id: prior ID of overlapping PLD lakes
        :type in_list_prior_id: list of str
        :param: out_pixcvec_lakeid = list of lake_id associated to each PIXC constituting observed lake. Unaffected pixels are empty str
        :type out_pixcvec_lakeid: numpy array of str

        :return: out_pixcvec_lakeid = list of lake_id associated to each PIXC constituting observed lake. Unaffected pixels are empty str
        :rtype: numpy array of str
        """
        logger = logging.getLogger(self.__class__.__name__)
        # Retrieve unaffected pixels
        pixc_out_of_lake_geom = np.where(out_pixcvec_lakeid == "")

        sub_lon = in_lon[pixc_out_of_lake_geom]
        sub_lat = in_lat[pixc_out_of_lake_geom]
        # Computation time: compute_pixcvec_lakeid_with_assignment_layer * 2,5 = compute_closest_polygon_with_kdtree
        if self.lake_assignment_layer:
            logger.debug(
                "%d / %d pixels affected using lake assignment geometries" % (np.sum(out_pixcvec_lakeid == ""), out_pixcvec_lakeid.size))
            out_pixcvec_lakeid[pixc_out_of_lake_geom] = self.compute_pixcvec_lakeid_with_assignment_layer(sub_lon, sub_lat,
                                                                                                            in_prior_geoms, in_list_prior_id)
        elif np.any(out_pixcvec_lakeid == ""):
            logger.debug(
                "%d / %d pixels affected using kd-tree" % (np.sum(out_pixcvec_lakeid == ""), out_pixcvec_lakeid.size))
            out_pixcvec_lakeid[pixc_out_of_lake_geom] = compute_closest_polygon_with_kdtree(sub_lon, sub_lat,
                                                                                            in_prior_geoms, in_list_prior_id)
        # ATTENTION : Different results !!
        # print(self.compute_pixcvec_lakeid_with_assignment_layer(in_lon, in_lat) == compute_closest_polygon_with_kdtree(in_lon, \
        #                in_lat, prior_geoms, prior_id))
        # compute_pixcvec_lakeid_with_assignment_layer is more precise.
        # compute_closest_polygon_with_kdtree less precise because computes the distance between pixels and polygon
        # coordinates and not polygon edges.

        return out_pixcvec_lakeid

    def compute_pixcvec_lakeid_with_assignment_layer(self, in_lon, in_lat, prior_geom_coords, list_prior_id):
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

        # 1. Filter assignment features following attributes
        request = "%s IN ('%s'" %(self.lakedb_id_name, list_prior_id[0])
        for lakedb_id in list_prior_id[1:]:
            request += ", '%s'" %(lakedb_id)
        request += ")"
        self.lake_assignment_layer.SetAttributeFilter(request)

        if self.lake_assignment_layer.GetFeatureCount() == 0:
            ## TODO :
            #   - this case should never happend : a lake in PLD should ALWAYS have an influence area
            #   - in Topocat.lake_catchment, a PLD lake may not have a lake catchment
            logger.debug("WARNING No assignment feature found for lake_id %s" % str(list_prior_id))
        else:
            logger.debug("Filter %d assignment features with request %s" % (self.lake_assignment_layer.GetFeatureCount(), request))

        # 2. Create memory layer to store pixcvec points
        ds, lyr = my_tools.load_pixels_to_mem_layer(in_lon, in_lat)

        # 3. Filter point of pixcvec with assignment features
        for infl_area_feat in self.lake_assignment_layer:
            lyr.SetSpatialFilter(infl_area_feat.GetGeometryRef())
            for feat_point in lyr:
                out_lakedb_id_pixcvec[feat_point.GetFID()] = infl_area_feat.GetField(self.lakedb_id_name)
            lyr.SetSpatialFilter(None)
        self.lake_assignment_layer.SetAttributeFilter(None)
        ds.Destroy()

        # 4. Compute lakedb_id for pixels located out of the assignment features
        unassigned_pixels = np.where(out_lakedb_id_pixcvec == "")
        nb_pt_ass_kd = unassigned_pixels[0].size
        nb_pt_ass_infl = out_lakedb_id_pixcvec.size - nb_pt_ass_kd

        if nb_pt_ass_infl > 0:
            logger.debug("%d pixels assigned using assignment layer" % nb_pt_ass_infl)
        if (self.lake_assignment_layer_name != "lake_catchment") & (nb_pt_ass_kd > 0):
            out_lakedb_id_pixcvec[unassigned_pixels] = compute_closest_polygon_with_kdtree(in_lon[unassigned_pixels], in_lat[unassigned_pixels], 
                                                                                           prior_geom_coords, list_prior_id)
            logger.debug("%d pixels assigned using kd tree" % nb_pt_ass_kd)

        return out_lakedb_id_pixcvec

    # ----------------------------------------

    def link_point_to_basin(self, in_point):
        """
        Link a point to a basin by considering intersection of both

        :param in_point: point to link to a basin
        :type in_point: ogr.Point

        :return: out_basin = IDs of the basin associated to the input point
        :rtype: string
        """
        
        # Init output
        # If no input continent is given, code is 010 => no continent used
        out_basin = my_var.BASIN_DEFAULT
        
        if self.basin_layer:
            
            # 1 - Compute intersection
            self.basin_layer.SetSpatialFilter(in_point)

            # 2 - Get basin code 
            if self.basin_layer.GetFeatureCount() == 0:
                out_basin = my_var.BASIN_4_OCEAN
                
            else:
                feat = self.basin_layer.GetNextFeature()
                basin_id = str(feat.GetField(self.basindb_id_name))
                out_basin = basin_id
                
            self.basin_layer.SetSpatialFilter(None)

        return out_basin

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
                out_basin_list.append(my_var.BASIN_4_OCEAN)  # Ocean
                
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
            out_basin_list.append(my_var.BASIN_DEFAULT)  # If no input continent is given, code is 010 => no continent used

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
        for code in [my_var.BASIN_4_OCEAN, my_var.BASIN_DEFAULT]:
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

    # ----------------------------------------
    
    def is_coastal(self, in_longitude, in_latitude):
        """
        Determine if a set of points falls in ocean or land, by looking at its 
        4 extreme points. If one falls in no basin, then the feature is considered as coastal

        :param in_longitude: longitude of the pixels of the water feature
        :type in_longitude: 1D-array of float
        :param in_latitude: latitude of the pixels of the water feature
        :type in_latitude: 1D-array of float

        :return: out_flag_coastal =True if the point falls in ocean, =False otherwise
        :rtype: boolean
        """
        
        out_flag_coastal = False
        
        if self.basin_layer:
        
            # 1 - Get extreme points
            # 1.1 - Longitude min            
            ind_long_min = np.argmin(in_longitude)
            point1 = ogr.Geometry(ogr.wkbPoint)
            point1.AddPoint(in_longitude[ind_long_min], in_latitude[ind_long_min])
            # 1.2 - Longitude max
            ind_long_max = np.argmax(in_longitude)
            point2 = ogr.Geometry(ogr.wkbPoint)
            point2.AddPoint(in_longitude[ind_long_max], in_latitude[ind_long_max])
            # 1.3 - Latitude min
            ind_lat_min = np.argmin(in_latitude)
            point3 = ogr.Geometry(ogr.wkbPoint)
            point3.AddPoint(in_longitude[ind_lat_min], in_latitude[ind_lat_min])
            # 1.4 - Latitude max
            ind_lat_max = np.argmax(in_latitude)
            point4 = ogr.Geometry(ogr.wkbPoint)
            point4.AddPoint(in_longitude[ind_lat_max], in_latitude[ind_lat_max])

            for point in [point1, point2, point3, point4]:
                basin = self.link_point_to_basin(point)
                if basin == my_var.BASIN_4_OCEAN:
                    out_flag_coastal = True
                    break
            
        return out_flag_coastal

    # ----------------------------------------
    
    def get_lakeid_with_mls(self):
        """
        Get the list of PLD lakes that have a Mean Lake Surface 
        
        :return: out_list_lakeid = list of PLD lakes that have a Mean Lake Surface 
        :rtype: out_list_lakeid = list()
        """
    
        out_list_lakeid = list()

        if "is_mls" in self.pld_field_lake_verified.keys():
            # 1 - Filter the lake layer over the features having is_mls=1
            self.lake_layer.SetAttributeFilter("is_mls = 1")
            
            # 2 - Retrieve their lake_id
            for cur_feat in self.lake_layer:
                out_list_lakeid.append(cur_feat.GetFieldAsString("lake_id"))
                
            # 3 - Reset the filter
            self.lake_layer.SetAttributeFilter(None)
        
        return out_list_lakeid

    
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
        
        # 2 - Set DB filename and type
        self.filename = in_lakedb_filename
        self.db_format = "shp"

        # 3 - Open database
        self.open_db(in_poly=in_poly)

        # 4 - Init fields name and type
        self.init_fields_name_and_type()

        if in_az_0_line and in_az_max_line:
            # 5 - Build border geometry
            self.build_border_geometry(in_az_0_line, in_az_max_line)

    # ----------------------------------------

    def open_db(self, in_poly=None):
        """
        Open PLD, optionnally spatially select polygons or with attribute lake_id, and copy lake_layer to memory

        :param in_poly: polygon to spatially select lakes from DB
        :type in_poly: ogr.Polygon
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")

        # 1 - Open shapefile in read-only access
        shp_driver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Shapefile driver
        shp_data_source = shp_driver.Open(self.filename, 0)

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

    # ----------------------------------------        
        
    def init_prior_layer(self, in_obj_lake, in_date=None):
        """
        Init layer of input lake object with prior lakes over the granule area
        
        :param in_obj_lake: lake object
        :type in_obj_lake: LakeSPPriorProduct or LakeAvgProduct
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("WARNING Layer CANNOT be init from lake database because it's in shapefile format")


#######################################
        

class LakeDbSqlite(LakeDb):
    """
    This class heritates from the main class to manage PLD in SQLite format
    """

    def __init__(self, in_lakedb_filename, 
                 in_tile_ref=None, in_continent_id=None, in_swath=None, in_basin_id=None,
                 in_poly=None, in_az_0_line=None, in_az_max_line=None, in_date=None):
        """
        Constructor

        :param in_lakedb_filename: full path of lake database file
        :type in_lakedb_filename: string
        :param in_tile_ref: reference of tile on which to focus (specific to LakeTile processing, in case of use of operational PLD); format = <3-digit_tile_number><swath_side>
        :type in_tile_ref: string
        :param in_continent_id: 1-digit continent identifier on which to focus (specific to LakeSP processing, in case of use of operational PLD)
        :type in_continent_id: int
        :param in_swath: "R" or "L" for Right ou Left swath
        :type in_swath: string
        :param in_basin_id: 2-digit basin code on which to focus (specific to LakeAvg processing, in case of use of operational PLD)
        :type in_basin_id: int
        :param in_poly: polygon to spatially select lake features (in case of user lake DB) 
        :type in_poly: OGRPolygon
        :param in_az_0_line: line delineating the bottom of the tile (specific to LakeTile processing, in case of use of user lake DB)
        :type in_az_0_line: OGRLineString
        :param in_az_max_line: line delineating the top of the tile (specific to LakeTile processing, in case of use of user lake DB)
        :type in_az_max_line: OGRLineString
        :param in_date: date of lake product acquisition
        :type in_date: datetime.datetime
        
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
        # Get instance of service config file
        cfg = service_config_file.get_instance()

        # 1 - Init LakeDb object
        super().__init__()
        
        # 2 - Set DB filename and type
        self.filename = in_lakedb_filename
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
            pass_tile_ds, pass_tile_layer = self.open_db(my_var.PLD_TABLE_PASS_TILE,
                                                         list(my_var.PLD_FIELD_TILE.values()),
                                                         ['text'] * len(my_var.PLD_FIELD_TILE),
                                                         field_id=in_tile_ref)
            nb_tile_features = pass_tile_layer.GetFeatureCount()
            if nb_tile_features != 1:
                logger.error("%d features correspond to tile %s" % (nb_tile_features, in_tile_ref))
                raise
            # Get feature
            feat = pass_tile_layer.GetNextFeature()
            # Retrieve continent
            continent_value = feat.GetField(my_var.PLD_FIELD_TILE["cont"])
            if continent_value == "0":
                self.flag_ocean = True
            # Retrieve the polygon of the tile
            poly_4selection = feat.GetGeometryRef()
            # poly_4selection can be multipolygon in the case if a tile with longitude=+=180.
            if poly_4selection.GetGeometryName() not in ["POLYGON", "MULTIPOLYGON"]:
                logger.warning("poly_4selection variable type %s invalid" % poly_4selection.GetGeometryName())
                poly_4selection = None
            # Retreive az_0 and az_max_lines
            az_0_line = ogr.CreateGeometryFromWkt(feat.GetField(my_var.PLD_FIELD_TILE["az_0_line"]))
            az_max_line = ogr.CreateGeometryFromWkt(feat.GetField(my_var.PLD_FIELD_TILE["az_max_line"]))
            # Close temporary layer
            pass_tile_ds.Destroy()

        elif in_swath is not None: # LakeSP processing ; select current swath pass
            pass_tile_ds, pass_tile_layer = self.open_db(my_var.PLD_TABLE_PASS_TILE,
                                                         [my_var.PLD_FIELD_TILE["id"]],
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
            # Close temporary layer
            pass_tile_ds.Destroy()
            
        else: # LakeAvg processing
            poly_4selection = in_poly
            # Fields pld_overlap_avg_l, pld_nb_tiles_l, pld_overlap_avg_r, pld_nb_tiles_r do not exists in PLD by Basin. Set to None
            self.pld_overlap_avg_l = None
            self.pld_nb_tiles_l = None
            self.pld_overlap_avg_r = None
            self.pld_nb_tiles_r = None

        # 3.2.1 - Build list of fields to load
        # Retrieve from my_variables
        self.list_lake_fields_to_load = my_var.PLD_FIELD_LAKE.copy()
        self.list_lake_field_type_to_load = my_var.PLD_FIELD_LAKE_TYPE.copy()
        # Remove some when LakeSP or LakeAvg processing
        if not in_tile_ref:  # For all processes except LakeTile, ie LakeSP+LakeAvg
            for key in ["overlap_avg_l", "nb_tiles_l", "overlap_avg_r", "nb_tiles_r"]:
                if key in self.list_lake_fields_to_load.keys():
                    self.list_lake_fields_to_load.pop(key)
                    self.list_lake_field_type_to_load.pop(key)
            if in_basin_id is not None:  # For LakeAvg process
                for key in ["ice_start", "ice_end", "is_mls"]:
                    if key in self.list_lake_fields_to_load.keys():
                        self.list_lake_fields_to_load.pop(key)
                        self.list_lake_field_type_to_load.pop(key)
        # 3.2.2 - Load table lake
        self.lake_ds, self.lake_layer = self.open_db(my_var.PLD_TABLE_LAKE,
                                                     list(self.list_lake_fields_to_load.values()), 
                                                     list(self.list_lake_field_type_to_load.values()),
                                                     poly=poly_4selection, 
                                                     start_value=start_value, date=in_date)

        # Save criteria
        self.poly_4selection = poly_4selection
        self.start_value = start_value
        
        # 3.3 - Table lake assignment
        try:
            self.lake_assignment_flag = True  # Set flag indicating if lake assignment table is used
            self.lake_assignment_ds, self.lake_assignment_layer = self.open_db(self.lake_assignment_layer_name,
                                                                               [self.lakedb_id_name],
                                                                               ['text'],
                                                                               poly=poly_4selection, start_value=start_value)
        except:
            logger.warning("No table named %s in the database" % my_var.PLD_TABLE_LAKE_INFL)
            self.lake_assignment_flag = False
        
        # 3.4 - Table basin
        self.basin_flag = True  # Set flag indicating basin table is used
        if in_basin_id is None: # LakeTile and LakeSP processing: <lon|lat>_<min|max> are not needed
            self.basin_ds, self.basin_layer = self.open_db(my_var.PLD_TABLE_BASIN, 
                                                           [self.basindb_id_name],
                                                           ['text'],
                                                           poly=poly_4selection, start_value=start_value)
        else: # LakeAvg processing: <lon|lat>_<min|max> are needed
            self.basin_ds, self.basin_layer = self.open_db(my_var.PLD_TABLE_BASIN, 
                                                           [self.basindb_id_name, \
                                                            my_var.PLD_FIELD_BASIN["lon_min"], my_var.PLD_FIELD_BASIN["lon_max"], \
                                                            my_var.PLD_FIELD_BASIN["lat_min"], my_var.PLD_FIELD_BASIN["lat_max"]],
                                                           ['text', 'float', 'float', 'float', 'float'],
                                                           poly=poly_4selection, start_value=start_value)
        
        # 4 - Init fields name and type
        self.init_fields_name_and_type()

        # 5 - Build border geometry (only for LakeTile processing)
        if in_tile_ref is not None:
            self.build_border_geometry(az_0_line, az_max_line, tile_geom=poly_4selection)
        
        # 6 - Save CB basin coordinates (only for LakeAvg processing)
        self.basin_coords = None
        if in_basin_id is not None:
            cbb_basin = self.basin_layer.GetNextFeature()
            self.basin_coords = dict()
            basin_coords_keys = ["lon_min", "lon_max", "lat_min", "lat_max"]
            
            for key in basin_coords_keys:
                self.basin_coords[key] = cbb_basin.GetField(key)
            
            self.basin_layer.ResetReading()

        # 7 - ONLY FOR LakeTile PROCESSING
        # Identify lakes with pld_overlap_avg_[l/r]*pld_nb_tiles_[l/r] > MAX_NB_TILES_FULL_AZ
        # -1 for MAX_NB_TILES_FULL_AZ parameter disable selection of big lakes
        if in_tile_ref is not None:
            MAX_NB_TILES_FULL_AZ = cfg.getfloat("CONFIG_PARAMS", "MAX_NB_TILES_FULL_AZ")
            if MAX_NB_TILES_FULL_AZ != -1:
                if "L" in in_tile_ref:
                    self.lake_layer.SetAttributeFilter("overlap_avg_l*nb_tiles_l > %d " % MAX_NB_TILES_FULL_AZ)
                else:
                    self.lake_layer.SetAttributeFilter("overlap_avg_r*nb_tiles_r > %d " % MAX_NB_TILES_FULL_AZ)

                for lake in self.lake_layer:
                    lake_id = lake.GetField("lake_id")
                    logger.info("Lake %s covers %d tiles -> processed in LakeTile" % (lake_id, lake.GetField("nb_tiles_%s" % (in_tile_ref[-1]))))
                    self.big_lake_geom.append((lake_id, lake.GetGeometryRef().Clone()))

                self.lake_layer.SetAttributeFilter(None)
                self.lake_layer.ResetReading()

    # ----------------------------------------

    def open_db(self, in_table_name, 
                in_field_name_list, in_field_type_list, 
                field_id=None, start_value=None, poly=None, date=None):
        """
        Open database, optionnally select features depending on their identifier, location or related continent
        and copy result to shapefile memory layer

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
        :param date: date of measurement for LakeTile and LakeSP processing; list of 2 elements (start/stop time) for LakeAvg processing
        :type date: datetime.datetime or list
        
        :return: out_data_source = DataSource of the specified table
        :rtype: osgeo.ogr.DataSource
        :return: out_layer = layer associated to the specified table
        :rtype: osgeo.ogr.Layer
        """
        logger = logging.getLogger(self.__class__.__name__)
        #logger.debug("Loading fields %s of table %s from file %s" %(" ".join(in_field_name_list), in_table_name, self.filename))
        cfg = service_config_file.get_instance()

        # 1 - Open lake database in SQLite format
        db_connector, db_cursor = open_sqlite_db(self.filename, in_table_name,
                                                 verbose_level=cfg.get('LOGGING', 'logFileLevel'))

        # 2 - Check that requested fields exist in table
        field_name_list, field_type_list = check_fields_existence(in_table_name, in_field_name_list, in_field_type_list, db_cursor)

        # 3 - Select subset using field_id, poly or having identifier starting with continent_id
        sql_request = build_select_request(in_table_name, field_name_list, field_id=field_id, start_value=start_value,
                                           poly=poly, date=date)
        #logger.debug(sql_request)
        db_cursor.execute(sql_request)

        # 4 - Copy the result to a shapefile memory layer
        out_data_source, out_layer = copy_sqlite_cursor_to_mem_shp(field_name_list, field_type_list, db_cursor)

        # 5 - Close SQLite database
        logger.debug("%d features of table %s after focus over studied area" % (out_layer.GetFeatureCount(), in_table_name))
        db_connector.close()

        return out_data_source, out_layer

    # ----------------------------------------

    def open_db_to_mem_sqlite(self, in_table_name, 
                              in_field_name_list, in_field_type_list,
                              in_field_metadata, 
                              field_id=None, start_value=None, poly=None, date=None):
        """
        Open database, optionnally select features depending on their identifier, location or related continent
        and copy result to SQLite memory database

        :param in_table_name: name of table to load from DB
        :type in_table_name: string
        :param in_field_name_list: list of fieldnames to load from table, first element is the identifier
        :type in_field_name_list: list of string
        :param in_field_type_list: list of type of each in_field_name_list
        :type in_field_type_list: list of string
        :param in_field_metadata: dictionary having key=fieldname and 
            value=OrderedDict with (key, value)=(name, value) of each field
        :type in_field_metadata: OrderedDict
        :param field_id: to select only the feature having this identifier
        :type field_id: string
        :param start_value: to select features having their identifier beginning with it
        :type start_value: string
        :param poly: polygon to spatially select features
        :type poly: ogr.Polygon
        :param date: date of lake product acquisition
        :type date: datetime.datetime
        """
        # log logger = logging.getLogger(self.__class__.__name__)
        # log logger.debug("Loading fields %s of table %s from file %s" %(" ".join(in_field_name_list), in_table_name, self.filename))
        cfg = service_config_file.get_instance()
        
        # 1 - Create memory SQLite database
        self.db_connector, self.db_cursor = open_sqlite_db(":memory:", in_table_name, 
                                                           verbose_level=cfg.get('LOGGING', 'logFileLevel'))

        # 2 - Copy the result to a memory SQLite database
        self.copy_sqlite_cursor_to_mem_sqlite(in_table_name, 
                                              in_field_name_list, in_field_type_list,
                                              in_field_metadata, 
                                              field_id=field_id, start_value=start_value, poly=poly, date=date)

    def copy_sqlite_cursor_to_mem_sqlite(self, in_table_name, 
                                         in_field_name_list, in_field_type_list,
                                         in_field_metadata, 
                                         field_id=None, start_value=None, poly=None, date=None):
        """
        Copy the content of a SQLite request given by in_db_cursor in a shapefile memory layer
        
        :param in_table_name: name of table to load from DB
        :type in_table_name: string
        :param in_field_name_list: list of fieldnames to load from table, first element is the identifier
        :type in_field_name_list: list of string
        :param in_field_type_list: list of type of each in_field_name_list
        :type in_field_type_list: list of string
        :param in_field_metadata: dictionary having key=fieldname and 
            value=OrderedDict with (key, value)=(name, value) of each field
        :type in_field_metadata: OrderedDict
        :param field_id: to select only the feature having this identifier
        :type field_id: string
        :param start_value: to select features having their identifier beginning with it
        :type start_value: string
        :param poly: polygon to spatially select features
        :type poly: ogr.Polygon
        :param date: date of lake product acquisition
        :type date: datetime.datetime
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Attach source DB = the lake database
        sql_request = "ATTACH DATABASE '%s' AS source;" % self.filename
        #logger.debug(sql_request)
        self.db_cursor.execute(sql_request)
        
        # 2 - Create table in the memory SQLite database
        
        # 2.1 - Build defitinions of fields
        table_def = "lake_id INTEGER PRIMARY KEY"
        for key, value in in_field_metadata.items():
            
            if key != "lake_id":
            
                if value['type'] in ["int", "int4", "int9", "float", "text"]:
                    field_type = my_var.FORMAT_SQL_STR[value['type']]
                    field_default = my_var.FV_SHP[value['type']]
                    table_def += ", %s %s DEFAULT %s" % (key, field_type, field_default)
                    
                else:
                    logger.error("Format %s unknown for SQLite => field %s CANNOT be created" % (value["type"], key))
                    raise
            
        # 2.2 - Build and execute request
        sql_request = "CREATE TABLE %s(%s);" % (in_table_name, table_def)
        #logger.debug(sql_request)
        self.db_cursor.execute(sql_request)
        
        # 3 - Insert objects from source to target DB, ie from lake database to memory SQLite database
        
        # 3.0 - Build the real list of prior attributes to retrieve from PLD fields
        list_of_prior_att = list(my_var.MAPPING_PRIOR_ATT_WRT_PLD.keys())
        for fieldname in my_var.MAPPING_PRIOR_ATT_WRT_PLD.keys():
            if fieldname not in in_field_metadata.keys():
                list_of_prior_att.remove(fieldname)
        
        # 3.1 - Build list of PLD fieldnames to retrieve
        list_pld_fields_to_retrieve = list()
        for field in list_of_prior_att:
            list_pld_fields_to_retrieve.append(my_var.MAPPING_PRIOR_ATT_WRT_PLD[field])
        
        # 3.2 - Select subset using field_id, poly or having identifier starting with continent_id
        select_request = build_select_request("source.%s" % in_table_name, list_pld_fields_to_retrieve,
                                              flag_no_geom=True,
                                              field_id=field_id, start_value=start_value, poly=poly, date=date)

        # 3.3 - Build and execute request
        sql_request = "INSERT INTO %s(%s) %s" % (in_table_name, ','.join(list_of_prior_att), select_request)
        #logger.debug(sql_request)
        self.db_cursor.execute(sql_request)
        
        # 3.4 - Commit
        self.db_connector.commit()
        
        # 4 - Detach source DB
        sql_request = "DETACH source;"
        #logger.debug(sql_request)
        self.db_cursor.execute(sql_request)
        
        # Print info
        (lakes_nb,) = self.db_cursor.execute('SELECT count(*) from %s' % (in_table_name)).fetchone()
        logger.debug("%d features retrieved from table <%s>" % (lakes_nb, in_table_name))

    # ----------------------------------------

    def close_db(self):
        """
        Close database
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Close lake table
        logger.debug("- close lake table -")
        self.lake_ds.Destroy()
        
        # Close assgnment table
        if self.lake_assignment_flag == True:
            logger.debug("- close assignment table -")
            self.lake_assignment_ds.Destroy()
        
        # Close basin table
        if self.basin_flag == True:
            logger.debug("- close basin table -")
            self.basin_ds.Destroy()

    # ----------------------------------------
    
    def init_prior_layer(self, in_obj_lake, in_date=None):
        """
        Init layer of input lake object with prior lakes over the granule area
        
        :param in_obj_lake: lake object
        :type in_obj_lake: LakeSPPriorProduct or LakeAvgProduct
        :param in_date: date of lake product acquisition
        :type in_date: datetime.datetime
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Init layer from lake database")
        
        # 1 - Open lake table of lake database and store it in a memory SQLite database
        self.open_db_to_mem_sqlite(my_var.PLD_TABLE_LAKE, 
                                   self.list_lake_fields_to_load, self.list_lake_field_type_to_load,
                                   in_obj_lake.attribute_metadata,
                                   poly=self.poly_4selection, start_value=self.start_value, date=in_date)
        
        # 2 - Copy selected features to input memory layer
        res = self.db_cursor.execute('SELECT * from %s;' % my_var.PLD_TABLE_LAKE).fetchall()
        for db_feature in res:
            
            # 2.1 - Create empty feature
            tmp_feat = ogr.Feature(in_obj_lake.layer_defn)
    
            # 2.2 - Fill feature with attributes and geometry from SQLite request
            for ind, fieldname in enumerate(in_obj_lake.attribute_metadata.keys()):
                fieldvalue = db_feature[ind]
                
                # On remplace None par _FillValue
                if fieldvalue is None:
                    fieldvalue = my_var.FV_SHP[in_obj_lake.attribute_metadata[fieldname]["type"]]
                    
                else:
                    
                    # TODO: to update when default values in PLD have been updated
                    if (fieldname in ["lake_name", "reach_id", "ice_clim_f"]) and (str(fieldvalue) == ""):
                        fieldvalue = my_var.FV_SHP[in_obj_lake.attribute_metadata[fieldname]["type"]]
                        
                    elif (fieldname in ["p_res_id", "ice_dyn_f"]) and (int(fieldvalue) <= 0):
                        fieldvalue = my_var.FV_SHP[in_obj_lake.attribute_metadata[fieldname]["type"]]
                    
                    elif (fieldname in ["p_ref_wse", "p_date_t0", "p_ds_t0"]) and (str(fieldvalue).startswith("-999")):
                        fieldvalue = my_var.FV_SHP[in_obj_lake.attribute_metadata[fieldname]["type"]]
                    #==TODO-END

                if (in_obj_lake.attribute_metadata[fieldname]["type"] in ['str', 'text', 'str32', 'string']) and (len(str(fieldvalue))>254):
                    fieldvalue = fieldvalue[:254]
                    
                tmp_feat.SetField(fieldname, fieldvalue)
                
            # 2.3 - Add feature to output layer
            in_obj_lake.layer.CreateFeature(tmp_feat)
    
            # 2.4 - Close temporary feature
            tmp_feat.Destroy()
    
        # 3 - Reset reading pointer
        in_obj_lake.layer.ResetReading()
        
        # 4 - Close connection to lake database
        self.db_connector.close()
        
        # 5 - Set flag to true
        self.flag_prior_init = True
        
    
def open_sqlite_db(in_filename, in_table_name, verbose_level="DEBUG"):
    """
    Open a SQLite database
    
    :param in_filename: full path of SQLite database
    :type in_filename: string
    :param in_table_name: name of table to load from DB
    :type in_table_name: string
    :param verbose_level: verbose level (default=DEBUG)
    :type verbose_level: string
    
    :return: db_connector = connection object of the SQLite database with full path in input
    :rtype: db_connector = sqlite3.Connection
    :return: db_cursor = cursor object of the SQLite database
    :rtype: db_cursor = sqlite3.Cursor
    """
    logger = logging.getLogger(__name__)
    
    # 1 - Define the connector
    db_connector = sqlite3.connect(in_filename, timeout=10)
    
    # 2 - Load spatialite extension
    db_connector.enable_load_extension(True)
    db_connector.execute('SELECT load_extension("mod_spatialite")')
    
    # 3 - Define the cursor
    db_cursor = db_connector.cursor()
    
    # Print info
    if (verbose_level == 'DEBUG') and (in_filename != ":memory:"):
        (lakes_nb,) = db_cursor.execute('SELECT count(*) from %s' % (in_table_name)).fetchone()
        logger.debug("%d features stored in table <%s>" % (lakes_nb, in_table_name))
        
    return db_connector, db_cursor
        
    
def copy_sqlite_cursor_to_mem_shp(in_field_name_list, in_field_type_list, in_db_cursor):
    """
    Copy the content of a SQLite request given by in_db_cursor in a shapefile memory layer
    
    :param in_field_name_list: list of fieldnames to load from table, first element is the identifier
    :type in_field_name_list: list of string
    :param in_field_type_list: list of type of each in_field_name_list
    :type in_field_type_list: list of string
    :param in_db_cursor: lake database cursor 
    :type in_db_cursor: sqlite3.Cursor
    
    :return: out_data_source = DataSource of the specified table
    :rtype: osgeo.ogr.DataSource
    :return: out_layer = layer associated to the specified table
    :rtype: osgeo.ogr.Layer
    """
    
    # 1 - Init output memory layer
    mem_driver = ogr.GetDriverByName('MEMORY')  # Memory driver
    # 1.1 - Open the memory DataSource with write access
    out_data_source = mem_driver.CreateDataSource('memData')
    # 1.2 - Set spatial projection
    srs = ogr.osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    # 1.3 - Create memory layer
    out_layer = out_data_source.CreateLayer(str('layer'), srs=srs, geom_type=ogr.wkbMultiPolygon)
    # 1.4 - Create needed fields

    for ind, field_type in enumerate(in_field_type_list):
        out_layer.CreateField(ogr.FieldDefn(in_field_name_list[ind], my_var.FORMAT_OGR[field_type]))
    # Retrieve layer definition
    lyr_defn = out_layer.GetLayerDefn()
    
    # 2 - Copy selected features to output memory layer
    for db_feature in in_db_cursor:
        
        # 2.1 - Create empty feature
        tmp_feat = ogr.Feature(lyr_defn)

        # 2.2 - Fill feature with attributes and geometry from SQLite request
        for ind, fieldname in enumerate(in_field_name_list):
            if fieldname == "ice_clim_flag":
                try:
                    tmp_feat.SetField(fieldname, int(max(list(db_feature[ind]))))
                except:
                    tmp_feat.SetField(fieldname, db_feature[ind])
            else:
                tmp_feat.SetField(fieldname, db_feature[ind])
        poly = ogr.CreateGeometryFromWkt(db_feature[-1])
        tmp_feat.SetGeometry(poly)

        # 2.3 - Add feature to output layer
        out_layer.CreateFeature(tmp_feat)

        # 2.4 - Close temporary feature
        tmp_feat.Destroy()

    # 3 - Reset reading pointer
    out_layer.ResetReading()
    
    return out_data_source, out_layer


def check_fields_existence(in_table_name, in_field_name_list, in_field_type_list, db_cursor):
    """
    Check that requested fields listed in in_field_name_list exist in table in_table_name read by db_cursor. 
    Return updated in_field_name_list without unavailable fields.

    :param in_table_name: name of table to load from DB
    :type in_table_name: string
    :param in_field_name_list: list of fieldnames to load from table, first element is the identifier
    :type in_field_name_list: list of string
    :param in_field_type_list: list of type of each in_field_name_list
    :type in_field_type_list: list of string
    :param in_db_cursor: lake database cursor
    :type in_db_cursor: sqlite3.Cursor

    :return: out_field_name_list = list of fieldnames to load from table, first element is the identifier
    :rtype: out_field_name_list = list of string
    :return: out_field_type_list = list of type of each out_field_name_list
    :rtype: out_field_type_list = list of string
    """
    logger = logging.getLogger(__name__)

    # 1 - Retrieve name of fields in current table
    list_fields_in_table = []
    db_cursor.execute("PRAGMA table_info(%s);" % in_table_name)
    for db_feature in db_cursor:
        id, field_name, field_type, _, _ , _= db_feature
        list_fields_in_table.append(field_name)

    # 2 - Check that requested fields exist in current table
    out_field_name_list = list()
    out_field_type_list = list()
    for field_name, field_type in zip(in_field_name_list, in_field_type_list):
        if field_name in list_fields_in_table:
            out_field_name_list.append(field_name)
            out_field_type_list.append(field_type)
        else:
            logger.warning("Field %s is not available in table %s ... will be set to None" % (field_name, in_table_name))

    return out_field_name_list, out_field_type_list


def build_select_request(in_table_name, in_field_name_list, flag_no_geom=False, field_id=None, start_value=None,
                         poly=None, date=None):
    """
    Build the SQL request from input parameters
    
    :param in_table_name: name of table to load from DB
    :type in_table_name: string
    :param in_field_name_list: list of fieldnames to load from table, first element is the identifier
    :type in_field_name_list: list of string
    :param flag_no_geom: =True to also retrieve the geometry, =False otherwise (default)
    :type flag_no_geom: boolean
    :param field_id: to select only the feature having this identifier
    :type field_id: string
    :param start_value: to select features having their identifier beginning with it
    :type start_value: string
    :param poly: polygon to spatially select features
    :type poly: ogr.Polygon
    :param date: date of lake product acquisition
    :type date: datetime.datetime

    :return: out_request = the SQL request
    :rtype: out_request = string
    """
    
    # 1 - Init SELECT
    if flag_no_geom:
        out_request = "SELECT %s FROM %s " %(",".join(in_field_name_list), in_table_name)
    else:
        out_request = "SELECT %s, AsText(ST_Multi(geometry)) FROM %s " %(",".join(in_field_name_list), in_table_name)
    
    # 2 - Build the WHERE options
    cmd_where = []
    
    if field_id is not None:
        cmd_where.append("%s LIKE '%s' " % (in_field_name_list[0], field_id))
        
    if start_value is not None:
        cmd_where.append(" {} LIKE '{}%' ".format(in_field_name_list[0], start_value))
        
    if poly is not None:
        poly.FlattenTo2D()  # Transform 3D geometry into 2D geometry (necessary for spatialite query)
        cmd_where.append(" ST_Intersects(GeomFromText('%s'), %s.geometry) " % (poly.ExportToWkt(), in_table_name))
        
    # 3 - Build the final request
    if cmd_where:
        out_request += " WHERE " + " AND ".join(cmd_where) + ";"
    else:
        out_request += ";"
        
    # 4 - Specific case of ice_clim_flag
    if date is not None:
        if isinstance(date, list):  # In case of LakeAvg
            out_request = out_request.replace("ice_clim_flag", "SUBSTR(ice_clim_flag || ice_clim_flag2, %d, %d)" % \
                                              (date[0].timetuple().tm_yday, date[1].timetuple().tm_yday-date[0].timetuple().tm_yday+1))
        else:  # In case of LakeTile/LakeSP
            out_request = out_request.replace("ice_clim_flag", "SUBSTR(ice_clim_flag || ice_clim_flag2, %d, 1)" % \
                                              (date.timetuple().tm_yday))
        
    return out_request
        

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
    logger = logging.getLogger(__name__)

    lon_coords, lat_coords, prior_id_list = np.array([]), np.array([]), np.array([], dtype=object)
    for geom, id in zip(prior_geoms, prior_id):
        if geom.GetGeometryName() == "POLYGON":
            lon_coords, lat_coords, prior_id_list = append_coords_from_polygon(geom, id, lon_coords, lat_coords, prior_id_list)
        elif geom.GetGeometryName() == "MULTIPOLYGON":
            for poly in geom:
                lon_coords, lat_coords, prior_id_list = append_coords_from_polygon(poly, id, lon_coords, lat_coords, prior_id_list)
        else :
            logger.error("Unknown geometry %s " % geom.GetGeometryName())

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

def append_coords_from_polygon(polygon, id, lon_coords, lat_coords, prior_id_list):
    """
        Append coordinates and prior id from polygon to lon_coords, lat_coords, prior_id_list.

        :param polygon: polygon
        :type polygon: ogr.Geometry of type ogr.wkbPolygon
        :param id: lake_id of polygon
        :type id: str
        :param lon_coords: array of longitudes
        :type lon_coords: 1D array of floats
        :param lat_coords: array of latitudes
        :type lat_coords: 1D array of floats
        :param prior_id: List of prior ID from lake DB
        :type prior_id: 1D array of str

        :return: tuple of array of longitudes, array of latitudes and array of lake_id
        :rtype: tuple of (1D array of floats, 1D array of floats, 1D array of str)
    """
    for ring in polygon:
        coords_array = get_coords_from_ring(ring)
        lon_coords = np.append(lon_coords, coords_array[0::2])
        lat_coords = np.append(lat_coords, coords_array[1::2])
        prior_id_list = np.append(prior_id_list, [id] * coords_array[1::2].size)
    return lon_coords, lat_coords, prior_id_list

def get_coords_from_ring(ring):
    """
        Get longitude and latitude coordinates from ring

        :param ring: ring
        :type ring: ogr.Geometry of type ogr.wkbLineString or ogr.wkbMultiLineString

        :return: array of longitude and latitudes
        :rtype: 1D array of floats
    """
    # Get the coordinates with json library
    json_geom = json.loads(ring.ExportToJson())

    if json_geom["type"] == 'LineString':
        coords = list(json_geom['coordinates'])
        coords_array = np.array(coords).ravel()
    elif json_geom["type"] == 'MultiLineString':
        coords_array = np.array([poly for poly in json_geom['coordinates']]).ravel()
        while coords_array.dtype == 'object':  # array of list : depending if multipolygon contains polygons with holes
            new_coords = np.array([])
            for x in coords_array:
                new_coords = np.append(new_coords, np.array(x).ravel())
            coords_array = np.array(new_coords).ravel()
    else:
        coords_array = np.array([])
    return coords_array

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
        - pld_info / dict: info from PLD for the prior lake; key=attribute name
        - ok_to_compute_stocc / boolean: =True if storage change is computable, =False otherwise
        """
        
        # 0 - Init variables
        self.lake_id = in_lakeid
        self.geom = None
        self.pld_info = dict()  # Dictionary of available PLD infos

        if in_obj_lakedb.lake_layer: # In case a PLD is used
            
            # 1 - Select feature given its identifier
            if in_obj_lakedb.lakedb_id_type == "String":
                in_obj_lakedb.lake_layer.SetAttributeFilter("%s like '%s'" % (in_obj_lakedb.lakedb_id_name, str(in_lakeid)))
            else:
                in_obj_lakedb.lake_layer.SetAttributeFilter("%s = %s" % (in_obj_lakedb.lakedb_id_name, str(in_lakeid)))
            pld_lake_feat = in_obj_lakedb.lake_layer.GetNextFeature()

            # TODO: Remove test as soon as possible : issue should be solved with PLD_V103
            if pld_lake_feat:

                # 2.1 - Retrieve geometry
                self.geom = pld_lake_feat.GetGeometryRef().Clone()

                # 2.2 - Retrieve prior information when exists
                for pld_fieldname in in_obj_lakedb.pld_field_lake_verified.values():
                    fieldvalue = pld_lake_feat.GetField(pld_fieldname)

                    if fieldvalue is not None:

                        # TODO: to update when default values in PLD have been set to NULL
                        if pld_fieldname in [my_var.PLD_FIELD_LAKE["names"], my_var.PLD_FIELD_LAKE["reach"]]:
                            if fieldvalue != "":
                                self.pld_info[pld_fieldname] = fieldvalue

                        elif pld_fieldname == my_var.PLD_FIELD_LAKE["res_id"]:
                            if int(fieldvalue) > 0:
                                self.pld_info[pld_fieldname] = fieldvalue

                        elif pld_fieldname.startswith("ice_"):
                            if int(fieldvalue) >= 0:
                                self.pld_info[pld_fieldname] = fieldvalue

                        elif pld_fieldname in [my_var.PLD_FIELD_LAKE["max_wse"], my_var.PLD_FIELD_LAKE["ref_date"]]:
                            if not str(fieldvalue).startswith("-999"):
                                self.pld_info[pld_fieldname] = fieldvalue

                        elif type(fieldvalue) == float:
                            if float(fieldvalue) > -9e6:
                                self.pld_info[pld_fieldname] = fieldvalue

                        else:
                            self.pld_info[pld_fieldname] = fieldvalue
                        #==TODO-end

            else:
                logging.debug("WARNING lake_id %s not found in PLD => Attributes set to None" % str(in_lakeid))

            # 3 - Release filter
            in_obj_lakedb.lake_layer.SetAttributeFilter(None)
 
        # 4 - Set flag to compute storage change if reference data are available
        self.ok_to_compute_stocc = False
        if (my_var.PLD_FIELD_LAKE["max_wse"] in self.pld_info.keys()) and (my_var.PLD_FIELD_LAKE["max_area"] in self.pld_info.keys()):
            self.ok_to_compute_stocc = True

    # ----------------------------------------
      
    def format_attributes(self):
        """
        Map PLD attributes to fields expected in the LakeProduct
        - Key in self.pld_info is the PLD fieldname
        - Key in out_attributes is the LakeProduct fieldname
        
        :return: out_attributes = attributes correctly formatted
        :rtype: out_attributes = dict
        """
        
        # Init output dictionary
        out_attributes = dict()
        
        for product_fieldname, pld_fieldname in my_var.MAPPING_PRIOR_ATT_WRT_PLD.items():
            if product_fieldname == "lake_id":
                out_attributes[product_fieldname] = self.lake_id
            elif pld_fieldname in self.pld_info.keys():
                out_attributes[product_fieldname] = self.pld_info[pld_fieldname]
            else:  # Empty fields as _FillValue
                out_attributes[product_fieldname] = my_var.FV_SHP[my_var.PLD_FIELD_LAKE_TYPE[pld_fieldname]]
            
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
                            - in_list_obs[obs_id]["alpha"] = proprotionnal coefficient related to area of observed lake wrt all observed lakes linked to the prior lake
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
        
        # 0 - Init variables
        stocc_vars = dict()
        for key in [my_var.PLD_FIELD_LAKE["max_area"],my_var.PLD_FIELD_LAKE["max_area_u"],
                    my_var.PLD_FIELD_LAKE["max_wse"],my_var.PLD_FIELD_LAKE["max_wse_u"],my_var.PLD_FIELD_LAKE["ref_ds"]]:
            stocc_vars[key] = None
            if key in self.pld_info.keys():
                stocc_vars[key] = self.pld_info[key]                
        
        # 1 - Linear storage change
        # 1.1 - Compute raw linear storage change
        stoc_val, stoc_u = storage_change.stocc_linear_basic(in_list_obs,
                                                             stocc_vars[my_var.PLD_FIELD_LAKE["max_area"]],
                                                             stocc_vars[my_var.PLD_FIELD_LAKE["max_area_u"]],
                                                             stocc_vars[my_var.PLD_FIELD_LAKE["max_wse"]],
                                                             stocc_vars[my_var.PLD_FIELD_LAKE["max_wse_u"]])
        # 1.2 - Deduce linear storage change value
        out_delta_s_l = my_var.FV_REAL
        if stoc_val is not None:
            if stocc_vars[my_var.PLD_FIELD_LAKE["ref_ds"]] is not None:
                out_delta_s_l = stoc_val - stocc_vars[my_var.PLD_FIELD_LAKE["ref_ds"]]
            else:
                out_delta_s_l = stoc_val
        # 1.3 - Deduce linear storage change error
        out_ds_l_u = my_var.FV_REAL
        if stoc_u is not None:
            out_ds_l_u = stoc_u
        
        # 2 - Quadratic storage change
        # 2.1 - Compute raw linear storage change
        stoc_val, stoc_u = storage_change.stocc_quadratic_basic(in_list_obs, 
                                                                stocc_vars[my_var.PLD_FIELD_LAKE["max_area"]],
                                                                stocc_vars[my_var.PLD_FIELD_LAKE["max_area_u"]],
                                                                stocc_vars[my_var.PLD_FIELD_LAKE["max_wse"]],
                                                                stocc_vars[my_var.PLD_FIELD_LAKE["max_wse_u"]])
        # 2.2 - Deduce quadratic storage change value
        out_delta_s_q = my_var.FV_REAL
        if stoc_val is not None:
            if stocc_vars[my_var.PLD_FIELD_LAKE["ref_ds"]] is not None:
                out_delta_s_q = stoc_val - stocc_vars[my_var.PLD_FIELD_LAKE["ref_ds"]]
            else:
                out_delta_s_q = stoc_val
        # 2.3 - Deduce quadratic storage change error
        out_ds_q_u = my_var.FV_REAL
        if stoc_u is not None:
            out_ds_q_u = stoc_u
        
        return out_delta_s_l, out_ds_l_u, out_delta_s_q, out_ds_q_u
