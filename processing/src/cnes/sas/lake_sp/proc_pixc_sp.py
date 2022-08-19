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
# VERSION:4.0.0:DM:#91:2022/05/05:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: proc_pixc_sp.py
   :synopsis: Deals with subset of pixel cloud, for objects located at top or bottom edge of a tile; i.e. gather pixels retrieved from all 
              L2_HR_LakeTile_Edge files 
    Created on 09/27/2017

.. moduleauthor:: Cécile Cazals - CS

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import absolute_import, division, print_function, unicode_literals 

import datetime
from dateutil import parser
import logging
import numpy as np
from osgeo import ogr

import cnes.common.lib.my_netcdf_file as my_nc
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_segmentation as my_segmentation
import cnes.common.lib.my_variables as my_var
import cnes.common.lib_lake.lake_db as lake_db

import cnes.sas.lake_tile.proc_pixc as proc_pixc

import cnes.common.service_config_file as service_config_file


class PixCEdgeSwath(proc_pixc.PixelCloud):
    """
    class PixCEdgeSwath
    This class is designed to process all L2_HR_LakeTile_Edge files of one single swath. 
    LakeTile_Edge files contain pixel cloud information (from L2_HR_PIXC) for pixels involved in lakes located 
    at the edges of tiles in the along-track direction.
    """

    def __init__(self, in_cycle_num, in_pass_num, in_continent_id, in_swath_side):
        """
        Constructor: init aggregation of all LakeTile_Edge files of one half-swath over a continent-pass
        
        The LakeTile_Edge file contains only needed PixC variables information, but also additional information like:

            - edge_loc field, the edge location of the lake: top of the tile, bottom of the tile or both top and bottom (0=bottom 1=top 2=both)
            - edge_label contains the object labels retrieved from PGE_LakeTile labeling process.
            - edge_index contains the L2_HR_PIXC initial pixels indices. This information is needed to update the improved geoloc and tag fields
              of LakeTile_PIXCVec files.

        This class processes all LakeTile_Edge files of a single swath to gather all entities at the along-track edge of tile into lake entities.

        :param in_cycle_num: cycle number
        :type in_cycle_num: int
        :param in_pass_num: pass number
        :type in_pass_num: int
        :param in_continent_id: 2-letter continent identifier
        :type in_continent_id: string
        :param in_swath_side: 1-letter swath side; R=right or L=left
        :type in_swath_side: string

        Variables of the object:
            
        - Global PixCEdgeSwath infos:
            - cycle_num / int: cycle number (= global attribute named cycle_number in LakeTile_Edge)
            - pass_num / int: pass number (= global attribute named pass_number in LakeTile_Edge)
            - swath_side / string: R=Right L=Left swath side
            - date / int: date of acquisition (ex: 20000101)
            - nb_pix_range / int: number of pixels in range dimension (= global attribute named nb_pix_range in LakeTile_Edge)
            - nb_pixels / int: total number of pixels
            - pixc_metadata / dict: processing metadata
            - tile_poly / ogr.Polygon: polygon of the swath delineated by the PixC tiles

        - Variables specific to processing:
            - tile_num / list of int: List of tile number to process ex: [76, 77, 78]
            - tile_index / list of int: Tile_num reference of each pixel ex: [0, 0, 0, 1, 2, 2, 2, 2]
            - labels / 1D-array of int: arrays of new labels
            - is_boundary_pix / list of bool: if pixel belong to the first / last azimuth of single pass
            - near_range / list of float: store the near range of each tiles
            - slant_range_spacing / list of int: store the slant range samplig of each tiles (supposed to be 0.75)
            - classif_full_water / 1D-array of int: classification as if all pixels were interior water pixels
            - classif_without_dw / 1D-array of int: classification of all water pixels (ie with dark water pixels removed)
            - classif_dict/ dict: selection (by True or False) of pixels wrt their classification, for different processes; key=water/dark/4hull/4wse/4area
            - interferogram_flattened / 1D-array of complex: flattened interferogram
            - inundated_area / 1D-array of int: area of pixels covered by water
            - height_std_pix / 1D-array of float: height std
            - corrected_height / 1D-array of float: height corrected from geoid and other tide corrections
            
        - From L2_HR_LakeTile_Edge file:
            
            - Edge objects info:
                - edge_loc
                - edge_label
                - edge_index

            - From PIXC/pixel_group variables:
                - nb_pix_range / int: number of pixels in range dimension (= global attribute named interferogram_size_range in L2_HR_PIXC)
                - nb_pix_azimuth / int: number of pixels in azimuth dimension (= global attribute named interferogram_size_azimuth in L2_HR_PIXC)
                - classif / 1D-array of byte: classification value of pixels (= variable named classification in L2_HR_PIXC)
                - range_index / 1D-array of int: range indices of water pixels (= variable named range_index in L2_HR_PIXC)
                - azimuth_index / 1D-array of int: azimuth indices of water pixels (= variable named azimuth_index in L2_HR_PIXC)
                - interferogram / 2D-array of float: rare interferogram
                - power_plus_y / 1D-array of float: power for plus_y channel
                - power_minus_y / 1D-array of float: power for minus_y channel
                - water_frac / 1D array of float: water fraction
                - water_frac_uncert / 1D array of float: water fraction uncertainty
                - false_detection_rate / 1D array of float: alse detection rate
                - missed_detection_rate / 1D array of float: missed detection rate
                - bright_land_flag / 1D array of byte: bright land flag
                - layover_impact /1D array of float: layover impact
                - eff_num_rare_looks / 1D array of byte: number of rare looks
                - latitude / 1D-array of float: latitude of water pixels
                - longitude / 1D-array of float: longitude of water pixels
                - height / 1D-array of float: height of water pixels
                - cross_track / 1D-array of float: cross-track distance from nadir to center of water pixels
                - pixel_area / 1D-array of int: area of water pixels
                - dlatitude_dphase / 1D-array of float: sensitivity of latitude estimate to interferogram phase
                - dlongitude_dphase / 1D-array of float: sensitivity of longitude estimate to interferogram phase
                - dheight_dphase / 1D array of float: sensitivity of height estimate to interferogram phase
                - darea_dheight / 1D array of float: sensitivity of pixel area to reference height
                - eff_num_medium_looks / 1D array of int: number of medium looks
                - model_dry_tropo_cor / 1D array of float: dry troposphere vertical correction
                - model_wet_tropo_cor / 1D array of float: wet troposphere vertical correction
                - iono_cor_gim_ka / 1D array of float: ionosphere vertical correction
                - height_cor_xover / 1D array of float: crossover calibration height correction
                - geoid / 1D array of float: geoid
                - solid_earth_tide / 1D array of float: solid earth tide
                - load_tide_fes / 1D array of float: load tide height (FES2014)
                - load_tide_got / 1D array of float: load tide height (GOT4.10)
                - pole_tide / 1D array of float: pole tide height
                - classification_qual / 1D array of uint: flag that indicates the quality of the classification quantities 
                - geolocation_qual / 1D array of uint: flag that indicates the quality of the geolocation quantities 
                - wavelength / float: wavelength corresponding to the effective radar carrier frequency 
                - looks_to_efflooks / float: ratio between the number of actual samples and the effective number of independent samples during 
                  spatial averaging over a large 2-D area

            - From PIXC/tvp group variables:
                - nadir_time[_tai] / 1D-array of float: observation UTC [TAI] time of each nadir pixel 
                  (= variable named time[tai] in L2_HR_PIXC file)
                - nadir_longitude / 1D-array of float: longitude of each nadir pixel (= variable named longitude in L2_HR_PIXC file)
                - nadir_latitude / 1D-array of float: latitude of each nadir pixel (= variable named latitude in L2_HR_PIXC file)
                - nadir_[x|y|z] / 1D-array of float: [x|y|z] cartesian coordinates of each nadir pixel (= variables named [x|y|z] in L2_HR_PIXC file)
                - nadir_[vx|vy|vz] / 1D-array of float: velocity vector of each nadir pixel in cartesian coordinates 
                  (= variables named velocity_unit_[x|y|z] in L2_HR_PIXC file)
                - nadir_plus_y_antenna_[x|y|z] / 1D-array of float: position vector of the +y KaRIn antenna phase center in ECEF coordinates
                  (= variables named plus_y_antenna_[x|y|z] in L2_HR_PIXC file)
                - nadir_minus_y_antenna_[x|y|z] / 1D-array of float: position vector of the -y KaRIn antenna phase center in ECEF coordinates 
                  (= variables named minus_y_antenna_[x|y|z] in L2_HR_PIXC file)
        """
        # Get instance of service config file
        self.cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # 0 - Retrieve flag variables from configuration file
        self.classif_list = [int(i) for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_LIST").split(";")]
        self.classif_water = [i == "True" for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_WATER").split(";")]
        self.classif_4hull = [i == "True" for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_4HULL").split(";")]
        self.classif_4wse = [i == "True" for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_4WSE").split(";")]
        self.classif_4area = [i == "True" for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_4AREA").split(";")]

        # 1. Init Global PixCEdgeSwath infos
        self.cycle_num = in_cycle_num # Cycle number
        self.pass_num = in_pass_num # Orbit number
        self.swath_side = in_swath_side # swath R or L
        self.continent_id = in_continent_id
        self.continent_code = lake_db.compute_continent_code(in_continent_id)
        
        # 2. Init LakeTile_Edge variables
        # 2.1. Edge info
        self.edge_loc =  np.array(()).astype('int')  # Localization (top/bottom/both)
        self.edge_label = np.array(()).astype('int')  # Label in tile
        self.edge_index = np.array(()).astype('int')  # Index in original L2_HR_PIXC product
        # 2.2. From PIXC/pixel_cloud
        self.classif = np.array(()).astype('byte')
        self.range_index = np.array(()).astype('int')
        self.azimuth_index = np.array(()).astype('int')
        self.interferogram = np.array(())
        self.power_plus_y = np.array(())
        self.power_minus_y = np.array(())
        self.water_frac = np.array(())
        self.water_frac_uncert = np.array(())
        self.false_detection_rate = np.array(())
        self.missed_detection_rate = np.array(())
        self.bright_land_flag = np.array(()).astype('byte')
        self.layover_impact = np.array(())
        self.eff_num_rare_looks = np.array(())
        self.latitude = np.array(())
        self.longitude = np.array(())
        self.height = np.array(())
        self.cross_track = np.array(())
        self.pixel_area = np.array(())
        self.dlatitude_dphase = np.array(())
        self.dlongitude_dphase = np.array(())
        self.dheight_dphase = np.array(())
        self.darea_dheight = np.array(())
        self.eff_num_medium_looks = np.array(())
        self.model_dry_tropo_cor = np.array(())
        self.model_wet_tropo_cor = np.array(())
        self.iono_cor_gim_ka = np.array(())
        self.height_cor_xover = np.array(())
        self.geoid = np.array(())
        self.solid_earth_tide = np.array(())
        self.load_tide_fes = np.array(())
        self.load_tide_got = np.array(())
        self.pole_tide = np.array(())
        self.classification_qual = np.array(()).astype('uint32')
        self.geolocation_qual = np.array(()).astype('uint32')
        self.wavelength = -9999.0
        self.looks_to_efflooks = -9999.0
        # 2.3. From PIXC/tvp group
        self.nadir_time = np.array(())  # Time
        self.nadir_time_tai = np.array(())
        self.nadir_longitude = np.array(())  # Longitude
        self.nadir_latitude = np.array(())  # Latitude
        self.nadir_x = np.array(())  # X cartesian coordinate
        self.nadir_y = np.array(())  # Y cartesian coordinate
        self.nadir_z = np.array(())  # Z cartesian coordinate
        self.nadir_vx = np.array(())  # Velocity in X coordinate
        self.nadir_vy = np.array(())  # Velocity in Y coordinate
        self.nadir_vz = np.array(())  # Velocity in Z coordinate
        self.nadir_plus_y_antenna_x = np.array(())
        self.nadir_plus_y_antenna_y = np.array(())
        self.nadir_plus_y_antenna_z = np.array(())
        self.nadir_minus_y_antenna_x = np.array(())
        self.nadir_minus_y_antenna_y = np.array(())
        self.nadir_minus_y_antenna_z = np.array(())
        
        # 3. Init dictionary of SP_Edge metadata
        self.pixc_metadata = {}
        self.pixc_metadata["cycle_number"] = -9999
        self.pixc_metadata["pass_number"] = -9999
        self.pixc_metadata["continent_id"] = self.continent_id
        self.pixc_metadata["continent_code"] = self.continent_code
        self.pixc_metadata["basin_code"] = ""
        self.pixc_metadata["time_granule_start"] = ""
        self.pixc_metadata["time_granule_end"] = ""
        self.pixc_metadata["time_coverage_start"] = ""
        self.pixc_metadata["time_coverage_end"] = ""
        self.pixc_metadata["geospatial_lon_min"] = -9999
        self.pixc_metadata["geospatial_lon_max"] = -9999
        self.pixc_metadata["geospatial_lat_min"] = -9999
        self.pixc_metadata["geospatial_lat_max"] = -9999
        self.pixc_metadata["left_first_longitude"] = -9999
        self.pixc_metadata["left_first_latitude"] = -9999
        self.pixc_metadata["left_last_longitude"] = -9999
        self.pixc_metadata["left_last_latitude"] = -9999
        self.pixc_metadata["right_first_longitude"] = -9999
        self.pixc_metadata["right_first_latitude"] = -9999
        self.pixc_metadata["right_last_longitude"] = -9999
        self.pixc_metadata["right_last_latitude"] = -9999
        self.pixc_metadata["ellipsoid_semi_major_axis"] = ""
        self.pixc_metadata["ellipsoid_flattening"] = ""

        # 4. Init variables specific to processing
        self.date = None  # Date of acquisition
        self.nb_pix_range = 0  # Number of pixels in range
        self.nb_pixels = 0  # Number of pixels to process
        self.tile_poly = ogr.Geometry(ogr.wkbMultiPolygon)  # Tiles polygon union
        self.tile_num = []  # List of tile number to process ex: [76, 77, 78]
        self.tile_index = []  # Tile reference of each pixel
        self.labels = np.array(()).astype('int')  # Init labels to 0
        self.is_boundary_pix = []  # If pixel belongs to the first / last azimuth of single pass
        self.near_range = []
        self.slant_range_spacing = []
        self.list_classif_keys = ["water", "interior_water", "as_full_water", "dark", "without_dw", "4hull", "4wse", "4area"]
        self.classif_dict = dict()
        for key in self.list_classif_keys:
            self.classif_dict[key] = np.array((), dtype='bool')
        self.interferogram_flattened = np.array(())
        self.inundated_area = np.array(())
        self.height_std_pix = np.array(())
        self.corrected_height = np.array(())

    # ----------------------------------------

    def set_from_laketile_edge_file(self, in_laketile_edge_path_list):
        """
        This function loads NetCDF information from all tiles and merge all PIXC_edge info into 1D vectors

        :param in_laketile_edge_path_list: list of sorted full path of the NetCDF file to load
        :type in_laketile_edge_path_list: list of string
        """
        logger = logging.getLogger(self.__class__.__name__)

        in_laketile_edge_path_list.sort()
        # store number of pixels in azimuth
        full_path_azimuth_nb_pix = [0]
        # store azimith boundaries : begining and end of path + edges in case of missing tiles
        full_path_azimuth_boundary = [0]

        for lake_tile_edge_path in in_laketile_edge_path_list:  # For each LakeTile_Edge file

            logger.debug("Loading LakeTile_Edge file = %s" % lake_tile_edge_path)

            # Load data
            nb_pix_loaded, current_tile_number, current_nb_pix_azimuth, current_near_range, current_slant_range_spacing = \
                 self.load_laketile_edge_data(lake_tile_edge_path)

            if not current_tile_number:
                continue

            logger.debug("%d pixels loaded" % nb_pix_loaded)

            # Update nb_pixel
            self.nb_pixels += nb_pix_loaded

            # check if tile is neighboring the previous tile
            if len(self.tile_num) >= 1:
                previous_tile_number = self.tile_num[-1]

                cpt = 0
                while current_tile_number != previous_tile_number + 1:
                    logger.warning("Tile %d%s is missing, an empty tile is added" % (previous_tile_number + 1, self.swath_side))
                    # if current tile is not adjacent to previous tile, add an empty tile to tile_ref
                    self.tile_num.append(previous_tile_number + 1)
                    self.near_range.append(0)
                    self.slant_range_spacing.append(0)
                    previous_tile_number = self.tile_num[-1]

                    # in the case of missing tile, add an azimuth boundary for adjacent tile azimuth min and max
                    full_path_azimuth_boundary.append(np.sum(full_path_azimuth_nb_pix)-1)
                    full_path_azimuth_nb_pix.append(current_nb_pix_azimuth)
                    full_path_azimuth_boundary.append(np.sum(full_path_azimuth_nb_pix))

                    #  stop process if more than 50 tiles are missing
                    cpt += 1
                    if cpt > 50:
                        break

            # add current tile information
            self.tile_index += [current_tile_number] * nb_pix_loaded
            self.tile_num.append(current_tile_number)
            self.near_range.append(current_near_range)
            self.slant_range_spacing.append(current_slant_range_spacing)

            full_path_azimuth_nb_pix.append(current_nb_pix_azimuth)

        # Convert list to numpy array
        self.tile_index = np.array(self.tile_index)
        self.labels = np.zeros((self.nb_pixels))

        # add last azimuth index to boundaries
        full_path_azimuth_boundary.append(np.sum(full_path_azimuth_nb_pix)-1)

        # compute azimuth index of the beginning of each tile
        full_path_tile_azimuth_index = np.cumsum(full_path_azimuth_nb_pix)

        # Compute if pixels belong to the beging or end od path or missing tile egde
        self.is_boundary_pix = np.zeros((self.nb_pixels))

        # compute full path azimuth from the begining of path to the end of path
        full_path_azimuth_index = np.zeros((self.nb_pixels)).astype('int')
        for i, tile in enumerate(self.tile_num):
            full_path_azimuth_index[np.where(self.tile_index == tile)] = self.azimuth_index[np.where(self.tile_index == tile)].astype('int') +\
                                                                         full_path_tile_azimuth_index[i]

        # fill is boundary array if azimuth of pixel is a path boundary
        for az_boundary in full_path_azimuth_boundary :
            self.is_boundary_pix[full_path_azimuth_index == az_boundary] = 1

        logger.debug("%d PixC loaded for current swath" % self.nb_pixels)

    # ----------------------------------------

    def load_laketile_edge_data(self, in_laketile_edge_filename):
        """
        This function loads NetCDF information.
        
        :param in_laketile_edge_filename: full path of the NetCDF file to load
        :type in_laketile_edge_filename: string
        
        :return: out_nb_pix = number of pixels loaded
        :rtype: out_nb_pix = integer
        :return: out_tile_ref = reference of loaded tile
        :rtype: out_tile_ref = string
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Open input NetCDF file in reading mode
        pixc_edge_reader = my_nc.MyNcReader(in_laketile_edge_filename)

        # 2 - Get and check tile references (cycle, pass, swath)
        out_tile_number = int(pixc_edge_reader.get_att_value("tile_number"))

        current_cycle_num = int(pixc_edge_reader.get_att_value("cycle_number"))
        if current_cycle_num != self.cycle_num:
            logger.error("Cycle of tile %d do not match with SP product %d" %(current_cycle_num, self.cycle_num))

        current_pass_number = int(pixc_edge_reader.get_att_value("pass_number"))
        if current_pass_number != self.pass_num:
            logger.error("Pass of tile %d do not match with SP product %d" %(current_pass_number, self.pass_num))

        current_swath_side = str(pixc_edge_reader.get_att_value("swath_side"))
        if current_swath_side != self.swath_side:
            logger.error("Swath of tile %s do not match with PixCEdgeSwath %s" %(current_swath_side, self.swath_side))
        current_date = parser.parse(pixc_edge_reader.get_att_value("time_granule_start"))
        if not self.date:
            self.date = current_date
        # Stop the process if acquisition dates delta > 24h
        if (current_date - self.date) > datetime.timedelta(days=1):
            logger.error("Input LakeTile_Edge file do not correspond to the same aquisition date")

        current_continent_id = str(pixc_edge_reader.get_att_value("continent_id"))
        if not self.continent_id in current_continent_id:
            # If cur_continent_id do not belong to the EDGE SP product, do not add pixc info
            logger.error("Input LakeTile_Edge %s file do not correspond to the same continent %s" % (in_laketile_edge_filename, self.continent_id))
            retour = None, None, None, None, None

        else:
    
            # 4 - Update tile_poly with new tile
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(float(pixc_edge_reader.get_att_value("inner_first_longitude")),
                          float(pixc_edge_reader.get_att_value("inner_first_latitude")))
            ring.AddPoint(float(pixc_edge_reader.get_att_value("outer_first_longitude")),
                          float(pixc_edge_reader.get_att_value("outer_first_latitude")))
            ring.AddPoint(float(pixc_edge_reader.get_att_value("outer_last_longitude")),
                          float(pixc_edge_reader.get_att_value("outer_last_latitude")))
            ring.AddPoint(float(pixc_edge_reader.get_att_value("inner_last_longitude")),
                          float(pixc_edge_reader.get_att_value("inner_last_latitude")))
            ring.AddPoint(float(pixc_edge_reader.get_att_value("inner_first_longitude")),
                          float(pixc_edge_reader.get_att_value("inner_first_latitude")))
            current_tile_poly = ogr.Geometry(ogr.wkbPolygon)
            current_tile_poly.AddGeometry(ring)

            if not current_tile_poly.IsValid():
                logger.warning("Polygon tile of file %s is not valid" % in_laketile_edge_filename)
            else:
                self.tile_poly = self.tile_poly.Union(current_tile_poly)
    
            # 3 - Initialization of object variables if not already done
            if not self.tile_num:
                
                # 3.1 - Values for metadata
                self.pixc_metadata["cycle_number"] = self.cycle_num
                self.pixc_metadata["pass_number"] = self.pass_num
                self.pixc_metadata["time_granule_start"] = str(pixc_edge_reader.get_att_value("time_granule_start"))
                self.pixc_metadata["time_granule_end"] = str(pixc_edge_reader.get_att_value("time_granule_end"))
                self.pixc_metadata["time_coverage_start"] = str(pixc_edge_reader.get_att_value("time_coverage_start"))
                self.pixc_metadata["time_coverage_end"] = str(pixc_edge_reader.get_att_value("time_coverage_end"))
                self.pixc_metadata["geospatial_lon_min"] = np.double(pixc_edge_reader.get_att_value("geospatial_lon_min"))
                self.pixc_metadata["geospatial_lon_max"] = np.double(pixc_edge_reader.get_att_value("geospatial_lon_max"))
                self.pixc_metadata["geospatial_lat_min"] = np.double(pixc_edge_reader.get_att_value("geospatial_lat_min"))
                self.pixc_metadata["geospatial_lat_max"] = np.double(pixc_edge_reader.get_att_value("geospatial_lat_max"))
                if self.swath_side == "L":
                    self.pixc_metadata["left_first_longitude"] = np.double(pixc_edge_reader.get_att_value("outer_first_longitude"))
                    self.pixc_metadata["left_first_latitude"] = np.double(pixc_edge_reader.get_att_value("outer_first_latitude"))
                    self.pixc_metadata["left_last_longitude"] = np.double(pixc_edge_reader.get_att_value("outer_last_longitude"))
                    self.pixc_metadata["left_last_latitude"] = np.double(pixc_edge_reader.get_att_value("outer_last_latitude"))
                    self.pixc_metadata["right_first_longitude"] = np.double(pixc_edge_reader.get_att_value("inner_first_longitude"))
                    self.pixc_metadata["right_first_latitude"] = np.double(pixc_edge_reader.get_att_value("inner_first_latitude"))
                    self.pixc_metadata["right_last_longitude"] = np.double(pixc_edge_reader.get_att_value("inner_last_longitude"))
                    self.pixc_metadata["right_last_latitude"] = np.double(pixc_edge_reader.get_att_value("inner_last_latitude"))
                else:
                    self.pixc_metadata["left_first_longitude"] = np.double(pixc_edge_reader.get_att_value("inner_first_longitude"))
                    self.pixc_metadata["left_first_latitude"] = np.double(pixc_edge_reader.get_att_value("inner_first_latitude"))
                    self.pixc_metadata["left_last_longitude"] = np.double(pixc_edge_reader.get_att_value("inner_last_longitude"))
                    self.pixc_metadata["left_last_latitude"] = np.double(pixc_edge_reader.get_att_value("inner_last_latitude"))
                    self.pixc_metadata["right_first_longitude"] = np.double(pixc_edge_reader.get_att_value("outer_first_longitude"))
                    self.pixc_metadata["right_first_latitude"] = np.double(pixc_edge_reader.get_att_value("outer_first_latitude"))
                    self.pixc_metadata["right_last_longitude"] = np.double(pixc_edge_reader.get_att_value("outer_last_longitude"))
                    self.pixc_metadata["right_last_latitude"] = np.double(pixc_edge_reader.get_att_value("outer_last_latitude"))
                self.pixc_metadata["continent_id"] = self.continent_id
                self.pixc_metadata["continent_code"] = self.continent_code
                self.pixc_metadata["ellipsoid_semi_major_axis"] = str(pixc_edge_reader.get_att_value("ellipsoid_semi_major_axis"))
                self.pixc_metadata["ellipsoid_flattening"] = str(pixc_edge_reader.get_att_value("ellipsoid_flattening"))
                
                # 3.2 - Other variables
                self.nb_pix_range = int(pixc_edge_reader.get_att_value("interferogram_size_range"))
                self.wavelength = np.float(pixc_edge_reader.get_att_value("wavelength"))  # Wavelength
                # Ratio between the number of actual samples and the effective number of independent samples
                self.looks_to_efflooks = np.float(pixc_edge_reader.get_att_value("looks_to_efflooks"))  
    
            # 5 - Get number of pixels
            out_nb_pix = pixc_edge_reader.get_dim_value('points')
            out_nb_pix_azimuth = pixc_edge_reader.get_att_value("interferogram_size_azimuth")
            out_near_range = pixc_edge_reader.get_att_value("near_range")
            out_slant_range_spacing = pixc_edge_reader.get_att_value("nominal_slant_range_spacing")
            
            # 6 - Update metadata from PixC info
            self.pixc_metadata["time_coverage_start"] = min(self.pixc_metadata["time_coverage_start"], 
                                                            str(pixc_edge_reader.get_att_value("time_coverage_start")))
            self.pixc_metadata["time_coverage_end"] = max(self.pixc_metadata["time_coverage_end"], 
                                                          str(pixc_edge_reader.get_att_value("time_coverage_end")))
            self.pixc_metadata["geospatial_lon_min"] = min(self.pixc_metadata["geospatial_lon_min"], 
                                                           np.double(pixc_edge_reader.get_att_value("geospatial_lon_min")))
            self.pixc_metadata["geospatial_lon_max"] = max(self.pixc_metadata["geospatial_lon_max"], 
                                                           np.double(pixc_edge_reader.get_att_value("geospatial_lon_max")))
            self.pixc_metadata["geospatial_lat_min"] = min(self.pixc_metadata["geospatial_lat_min"], 
                                                           np.double(pixc_edge_reader.get_att_value("geospatial_lat_min")))
            self.pixc_metadata["geospatial_lat_max"] = max(self.pixc_metadata["geospatial_lat_max"], 
                                                           np.double(pixc_edge_reader.get_att_value("geospatial_lat_max")))
            if self.swath_side == "L":
                if self.pass_num%2 == 1:
                    self.pixc_metadata["left_first_longitude"] = min(self.pixc_metadata["left_first_longitude"], 
                                                                     np.double(pixc_edge_reader.get_att_value("outer_first_longitude")))
                    self.pixc_metadata["left_first_latitude"] = min(self.pixc_metadata["left_first_latitude"], 
                                                                    np.double(pixc_edge_reader.get_att_value("outer_first_latitude")))
                    self.pixc_metadata["left_last_longitude"] = max(self.pixc_metadata["left_last_longitude"], 
                                                                    np.double(pixc_edge_reader.get_att_value("outer_last_longitude")))
                    self.pixc_metadata["left_last_latitude"] = max(self.pixc_metadata["left_last_latitude"], 
                                                                   np.double(pixc_edge_reader.get_att_value("outer_last_latitude")))
                    self.pixc_metadata["right_first_longitude"] = min(self.pixc_metadata["right_first_longitude"], 
                                                                      np.double(pixc_edge_reader.get_att_value("inner_first_longitude")))
                    self.pixc_metadata["right_first_latitude"] = min(self.pixc_metadata["right_first_latitude"], 
                                                                     np.double(pixc_edge_reader.get_att_value("inner_first_latitude")))
                    self.pixc_metadata["right_last_longitude"] = max(self.pixc_metadata["right_last_longitude"], 
                                                                     np.double(pixc_edge_reader.get_att_value("inner_last_longitude")))
                    self.pixc_metadata["right_last_latitude"] = max(self.pixc_metadata["right_last_latitude"], 
                                                                    np.double(pixc_edge_reader.get_att_value("inner_last_latitude")))
                else:
                    self.pixc_metadata["left_first_longitude"] = min(self.pixc_metadata["left_first_longitude"], 
                                                                     np.double(pixc_edge_reader.get_att_value("outer_first_longitude")))
                    self.pixc_metadata["left_first_latitude"] = max(self.pixc_metadata["left_first_latitude"], 
                                                                    np.double(pixc_edge_reader.get_att_value("outer_first_latitude")))
                    self.pixc_metadata["left_last_longitude"] = max(self.pixc_metadata["left_last_longitude"], 
                                                                    np.double(pixc_edge_reader.get_att_value("outer_last_longitude")))
                    self.pixc_metadata["left_last_latitude"] = min(self.pixc_metadata["left_last_latitude"], 
                                                                   np.double(pixc_edge_reader.get_att_value("outer_last_latitude")))
                    self.pixc_metadata["right_first_longitude"] = min(self.pixc_metadata["right_first_longitude"], 
                                                                      np.double(pixc_edge_reader.get_att_value("inner_first_longitude")))
                    self.pixc_metadata["right_first_latitude"] = max(self.pixc_metadata["right_first_latitude"], 
                                                                     np.double(pixc_edge_reader.get_att_value("inner_first_latitude")))
                    self.pixc_metadata["right_last_longitude"] = max(self.pixc_metadata["right_last_longitude"], 
                                                                     np.double(pixc_edge_reader.get_att_value("inner_last_longitude")))
                    self.pixc_metadata["right_last_latitude"] = min(self.pixc_metadata["right_last_latitude"], 
                                                                    np.double(pixc_edge_reader.get_att_value("inner_last_latitude")))
            else:
                if self.pass_num%2 == 1:
                    self.pixc_metadata["left_first_longitude"] = min(self.pixc_metadata["left_first_longitude"], 
                                                                     np.double(pixc_edge_reader.get_att_value("inner_first_longitude")))
                    self.pixc_metadata["left_first_latitude"] = min(self.pixc_metadata["left_first_latitude"], 
                                                                    np.double(pixc_edge_reader.get_att_value("inner_first_latitude")))
                    self.pixc_metadata["left_last_longitude"] = max(self.pixc_metadata["left_last_longitude"], 
                                                                    np.double(pixc_edge_reader.get_att_value("inner_last_longitude")))
                    self.pixc_metadata["left_last_latitude"] = max(self.pixc_metadata["left_last_latitude"], 
                                                                   np.double(pixc_edge_reader.get_att_value("inner_last_latitude")))
                    self.pixc_metadata["right_first_longitude"] = min(self.pixc_metadata["right_first_longitude"], 
                                                                      np.double(pixc_edge_reader.get_att_value("outer_first_longitude")))
                    self.pixc_metadata["right_first_latitude"] = min(self.pixc_metadata["right_first_latitude"], 
                                                                     np.double(pixc_edge_reader.get_att_value("outer_first_latitude")))
                    self.pixc_metadata["right_last_longitude"] = max(self.pixc_metadata["right_last_longitude"], 
                                                                     np.double(pixc_edge_reader.get_att_value("outer_last_longitude")))
                    self.pixc_metadata["right_last_latitude"] = max(self.pixc_metadata["right_last_latitude"], 
                                                                    np.double(pixc_edge_reader.get_att_value("outer_last_latitude")))
                else:
                    self.pixc_metadata["left_first_longitude"] = min(self.pixc_metadata["left_first_longitude"], 
                                                                     np.double(pixc_edge_reader.get_att_value("inner_first_longitude")))
                    self.pixc_metadata["left_first_latitude"] = max(self.pixc_metadata["left_first_latitude"], 
                                                                    np.double(pixc_edge_reader.get_att_value("inner_first_latitude")))
                    self.pixc_metadata["left_last_longitude"] = max(self.pixc_metadata["left_last_longitude"], 
                                                                    np.double(pixc_edge_reader.get_att_value("inner_last_longitude")))
                    self.pixc_metadata["left_last_latitude"] = min(self.pixc_metadata["left_last_latitude"], 
                                                                   np.double(pixc_edge_reader.get_att_value("inner_last_latitude")))
                    self.pixc_metadata["right_first_longitude"] = min(self.pixc_metadata["right_first_longitude"], 
                                                                      np.double(pixc_edge_reader.get_att_value("outer_first_longitude")))
                    self.pixc_metadata["right_first_latitude"] = max(self.pixc_metadata["right_first_latitude"], 
                                                                     np.double(pixc_edge_reader.get_att_value("outer_first_latitude")))
                    self.pixc_metadata["right_last_longitude"] = max(self.pixc_metadata["right_last_longitude"], 
                                                                     np.double(pixc_edge_reader.get_att_value("outer_last_longitude")))
                    self.pixc_metadata["right_last_latitude"] = min(self.pixc_metadata["right_last_latitude"], 
                                                                    np.double(pixc_edge_reader.get_att_value("outer_last_latitude")))
            
            # 7 - Update vectors if there are pixels in the current LakeTile_Edge file
            if out_nb_pix > 0:
    
                # 7.1 - Add edge objects info
                self.edge_label = np.concatenate((self.edge_label, pixc_edge_reader.get_var_value("edge_label")))
                self.edge_index = np.concatenate((self.edge_index, pixc_edge_reader.get_var_value("edge_index")))
                self.edge_loc = np.concatenate((self.edge_loc, pixc_edge_reader.get_var_value("edge_loc")))
    
                # 7.2 - Add variables from PIXC/pixel_cloud
                tmp_classif = pixc_edge_reader.get_var_value("classification")
                self.classif = np.concatenate((self.classif, tmp_classif))
                
                classif_dict = dict()
                for key in self.list_classif_keys:
                    classif_dict[key] = np.full((out_nb_pix), False)
                
                # Simulate classification of edge + full water pixels
                # All PIXC are set to INTERIOR_WATER
                # LAND_EDGE + WATER_EDGE PIXC are set to WATER_EDGE
                classif_dict["as_full_water"] = np.zeros(out_nb_pix) + my_var.CLASSIF_INTERIOR_WATER
                classif_dict["as_full_water"][tmp_classif == my_var.CLASSIF_LAND_EDGE] = my_var.CLASSIF_WATER_EDGE
                classif_dict["as_full_water"][tmp_classif == my_var.CLASSIF_WATER_EDGE] = my_var.CLASSIF_WATER_EDGE
                # Keep only classification of water pixels (ie remove dark water flags)
                classif_dict["without_dw"] = np.copy(tmp_classif)
                classif_dict["without_dw"][tmp_classif == my_var.CLASSIF_DARK] = 0
            
                # PIXC indices wrt category of process
                # Boolean arrays of classification flags (=True if PIXC is of classif "key")
                tmp_classif_dict = dict()
                for classif_flag in self.classif_list:
                    tmp_classif_dict[classif_flag] = np.full((out_nb_pix), False)
                    tmp_classif_dict[classif_flag][tmp_classif == classif_flag] = True
                    if classif_flag == my_var.CLASSIF_INTERIOR_WATER:
                        classif_dict["interior_water"] = tmp_classif_dict[classif_flag]
                # PIXC indices of water and dark water pixels
                for classif_flag, ok_water in zip(self.classif_list, self.classif_water):
                    if ok_water:
                        classif_dict["water"] += tmp_classif_dict[classif_flag]
                    else:
                        classif_dict["dark"] += tmp_classif_dict[classif_flag]
                # PIXC indices of pixels to be used for hull computation
                for classif_flag, ok_4hull in zip(self.classif_list, self.classif_4hull):
                    if ok_4hull:
                        classif_dict["4hull"] += tmp_classif_dict[classif_flag]
                # PIXC indices of pixels to be used for wse computation
                for classif_flag, ok_4wse in zip(self.classif_list, self.classif_4wse):
                    if ok_4wse:
                        classif_dict["4wse"] += tmp_classif_dict[classif_flag]
                # PIXC indices of pixels to be used for area_total attribute computation
                for classif_flag, ok_4area in zip(self.classif_list, self.classif_4area):
                    if ok_4area:
                        classif_dict["4area"] += tmp_classif_dict[classif_flag]
                        
                # Concatenate classif_dict lists
                for key in self.list_classif_keys:
                    self.classif_dict[key] = np.concatenate((self.classif_dict[key], classif_dict[key]))
                
                self.range_index = np.concatenate((self.range_index, pixc_edge_reader.get_var_value("range_index")))
                self.azimuth_index = np.concatenate((self.azimuth_index, pixc_edge_reader.get_var_value("azimuth_index")))
                
                interferogram_value = pixc_edge_reader.get_var_value("interferogram")
                tmp_interferogram = interferogram_value[:,0] + 1j*interferogram_value[:,1]
                self.interferogram = np.concatenate((self.interferogram, tmp_interferogram))
                tmp_interferogram_flattened = 0 * tmp_interferogram
                self.interferogram_flattened = np.concatenate((self.interferogram_flattened, tmp_interferogram_flattened))
                self.power_plus_y = np.concatenate((self.power_plus_y, pixc_edge_reader.get_var_value("power_plus_y")))
                self.power_minus_y = np.concatenate((self.power_minus_y, pixc_edge_reader.get_var_value("power_minus_y")))
                
                tmp_water_frac = pixc_edge_reader.get_var_value("water_frac")
                self.water_frac = np.concatenate((self.water_frac, tmp_water_frac))
                self.water_frac_uncert = np.concatenate((self.water_frac_uncert, pixc_edge_reader.get_var_value("water_frac_uncert")))
                self.false_detection_rate = np.concatenate((self.false_detection_rate, pixc_edge_reader.get_var_value("false_detection_rate")))
                self.missed_detection_rate = np.concatenate((self.missed_detection_rate, pixc_edge_reader.get_var_value("missed_detection_rate")))
                self.bright_land_flag = np.concatenate((self.bright_land_flag, pixc_edge_reader.get_var_value("bright_land_flag")))
                self.layover_impact = np.concatenate((self.layover_impact, pixc_edge_reader.get_var_value("layover_impact")))
                self.eff_num_rare_looks = np.concatenate((self.eff_num_rare_looks, pixc_edge_reader.get_var_value("eff_num_rare_looks")))
                
                self.latitude = np.concatenate((self.latitude, pixc_edge_reader.get_var_value("latitude")))
                self.longitude = np.concatenate((self.longitude, pixc_edge_reader.get_var_value("longitude")))
                tmp_height = pixc_edge_reader.get_var_value("height")
                self.height = np.concatenate((self.height, tmp_height))
                
                self.cross_track = np.concatenate((self.cross_track, pixc_edge_reader.get_var_value("cross_track")))
                tmp_pixel_area = pixc_edge_reader.get_var_value("pixel_area")
                self.pixel_area = np.concatenate((self.pixel_area, tmp_pixel_area))
                tmp_inundated_area = np.copy(tmp_pixel_area)
                ind_ok = np.where(tmp_water_frac < my_var.FV_FLOAT)
                if len(ind_ok) > 0:
                    tmp_inundated_area[ind_ok] = tmp_pixel_area[ind_ok] * tmp_water_frac[ind_ok]
                self.inundated_area = np.concatenate((self.inundated_area, tmp_inundated_area))
                
                tmp_phase_noise_std = pixc_edge_reader.get_var_value("phase_noise_std")
                self.dlatitude_dphase = np.concatenate((self.dlatitude_dphase, pixc_edge_reader.get_var_value("dlatitude_dphase")))
                self.dlongitude_dphase = np.concatenate((self.dlongitude_dphase, pixc_edge_reader.get_var_value("dlongitude_dphase")))
                tmp_dheight_dphase = pixc_edge_reader.get_var_value("dheight_dphase")
                self.dheight_dphase = np.concatenate((self.dheight_dphase, tmp_dheight_dphase))
                self.darea_dheight = np.concatenate((self.darea_dheight, pixc_edge_reader.get_var_value("darea_dheight")))
                self.eff_num_medium_looks = np.concatenate((self.eff_num_medium_looks, pixc_edge_reader.get_var_value("eff_num_medium_looks")))
                
                self.model_dry_tropo_cor = np.concatenate((self.model_dry_tropo_cor, pixc_edge_reader.get_var_value("model_dry_tropo_cor")))
                self.model_wet_tropo_cor = np.concatenate((self.model_wet_tropo_cor, pixc_edge_reader.get_var_value("model_wet_tropo_cor")))
                self.iono_cor_gim_ka = np.concatenate((self.iono_cor_gim_ka, pixc_edge_reader.get_var_value("iono_cor_gim_ka")))
                self.height_cor_xover = np.concatenate((self.height_cor_xover, pixc_edge_reader.get_var_value("height_cor_xover")))
                
                tmp_geoid = pixc_edge_reader.get_var_value("geoid")
                self.geoid = np.concatenate((self.geoid, tmp_geoid))
                tmp_solid_earth_tide = pixc_edge_reader.get_var_value("solid_earth_tide")
                self.solid_earth_tide = np.concatenate((self.solid_earth_tide, tmp_solid_earth_tide))
                tmp_load_tide_fes = pixc_edge_reader.get_var_value("load_tide_fes")
                self.load_tide_fes = np.concatenate((self.load_tide_fes, tmp_load_tide_fes))
                tmp_load_tide_got = pixc_edge_reader.get_var_value("load_tide_got")
                self.load_tide_got = np.concatenate((self.load_tide_got, tmp_load_tide_got))
                tmp_pole_tide = pixc_edge_reader.get_var_value("pole_tide")
                self.pole_tide = np.concatenate((self.pole_tide, tmp_pole_tide))
                
                tmp_classification_qual = pixc_edge_reader.get_var_value("classification_qual")
                self.classification_qual = np.concatenate((self.classification_qual, tmp_classification_qual))
                tmp_geolocation_qual = pixc_edge_reader.get_var_value("geolocation_qual")
                self.geolocation_qual = np.concatenate((self.geolocation_qual, tmp_geolocation_qual))
    
                # 7.3 - Info of the nadir point associated to the PixC
                self.nadir_time = np.concatenate((self.nadir_time, pixc_edge_reader.get_var_value("nadir_time")))
                self.nadir_time_tai = np.concatenate((self.nadir_time, pixc_edge_reader.get_var_value("nadir_time_tai")))
                self.nadir_longitude = np.concatenate((self.nadir_longitude, pixc_edge_reader.get_var_value("nadir_longitude")))
                self.nadir_latitude = np.concatenate((self.nadir_latitude, pixc_edge_reader.get_var_value("nadir_latitude")))
                self.nadir_x = np.concatenate((self.nadir_x, pixc_edge_reader.get_var_value("nadir_x")))
                self.nadir_y = np.concatenate((self.nadir_y, pixc_edge_reader.get_var_value("nadir_y")))
                self.nadir_z = np.concatenate((self.nadir_z, pixc_edge_reader.get_var_value("nadir_z")))
                self.nadir_vx = np.concatenate((self.nadir_vx, pixc_edge_reader.get_var_value("nadir_vx")))
                self.nadir_vy = np.concatenate((self.nadir_vy, pixc_edge_reader.get_var_value("nadir_vy")))
                self.nadir_vz = np.concatenate((self.nadir_vz, pixc_edge_reader.get_var_value("nadir_vz")))
                self.nadir_plus_y_antenna_x = np.concatenate((self.nadir_plus_y_antenna_x, pixc_edge_reader.get_var_value("nadir_plus_y_antenna_x")))
                self.nadir_plus_y_antenna_y = np.concatenate((self.nadir_plus_y_antenna_y, pixc_edge_reader.get_var_value("nadir_plus_y_antenna_y")))
                self.nadir_plus_y_antenna_z = np.concatenate((self.nadir_plus_y_antenna_z, pixc_edge_reader.get_var_value("nadir_plus_y_antenna_z")))
                self.nadir_minus_y_antenna_x = np.concatenate((self.nadir_minus_y_antenna_x, \
                                                               pixc_edge_reader.get_var_value("nadir_minus_y_antenna_x")))
                self.nadir_minus_y_antenna_y = np.concatenate((self.nadir_minus_y_antenna_y, \
                                                               pixc_edge_reader.get_var_value("nadir_minus_y_antenna_y")))
                self.nadir_minus_y_antenna_z = np.concatenate((self.nadir_minus_y_antenna_z, \
                                                               pixc_edge_reader.get_var_value("nadir_minus_y_antenna_z")))
                
                # 7.4 - Set bad PIXC height std to high number to deweight 
                # instead of giving infs/nans
                tmp_height_std_pix = np.abs(tmp_phase_noise_std * tmp_dheight_dphase)
                bad_num = 1.0e5
                tmp_height_std_pix[tmp_height_std_pix<=0] = bad_num
                tmp_height_std_pix[np.isinf(tmp_height_std_pix)] = bad_num
                tmp_height_std_pix[np.isnan(tmp_height_std_pix)] = bad_num
                self.height_std_pix = np.concatenate((self.height_std_pix, tmp_height_std_pix))
            
                # 7.5 - Compute height wrt the geoid and apply tide corrections
                # Compute indices of PIXC for which corrections are all valid
                valid_geoid = np.where(np.isfinite(tmp_geoid))[0]
                valid_solid_earth_tide = np.where(np.isfinite(tmp_solid_earth_tide))[0]
                valid_pole_tide = np.where(np.isfinite(tmp_pole_tide))[0]
                valid_load_tide_fes = np.where(np.isfinite(tmp_load_tide_fes))[0]
                inter1 = np.intersect1d(valid_geoid, valid_solid_earth_tide)
                inter2 = np.intersect1d(valid_pole_tide, inter1)
                ind_valid_corr = np.intersect1d(valid_load_tide_fes, inter2)
                # Compute corrected height for these PIXC
                tmp_corrected_height = np.zeros(out_nb_pix) + np.nan
                tmp_corrected_height[ind_valid_corr] = tmp_height[ind_valid_corr] \
                                                       - tmp_geoid[ind_valid_corr] \
                                                       - tmp_solid_earth_tide[ind_valid_corr] \
                                                       - tmp_pole_tide[ind_valid_corr] \
                                                       - tmp_load_tide_fes[ind_valid_corr]
                self.corrected_height = np.concatenate((self.corrected_height, tmp_corrected_height))
    
            # 8 - Close file
            pixc_edge_reader.close()

            retour = out_nb_pix, out_tile_number, out_nb_pix_azimuth, out_near_range, out_slant_range_spacing

        return retour
        
    # ----------------------------------------

    def swath_global_relabeling(self):
        """
        This function gives new labels to entities gathered at tile edges.
        """
        cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Deal with all edges

        for i_edge, tile_num in enumerate(self.tile_num[:-1]):  # Loop over tile edges

            logger.debug(" ***** Processing edge of tiles %d and %d *****" % (tile_num, tile_num + 1))

            # 1.1 - Get indices of pixels processed at the current edge
            tile_idx1 = np.where(np.logical_and(self.tile_index == tile_num, np.logical_or(self.edge_loc == 1, self.edge_loc == 2)))[0]
            tile_idx2 = np.where(np.logical_and(self.tile_index == tile_num + 1, np.logical_or(self.edge_loc == 0, self.edge_loc == 2)))[0]


            # 1.2 - If one tile does not have pixel to process, continue to next iteration
            if (tile_idx1.size == 0) or (tile_idx2.size == 0):
                logger.debug("")
                continue
            logger.debug(" %d pixels of tile %d will be match with %d pixels of tile %d" % (tile_idx1.size, tile_num, tile_idx2.size, tile_num + 1))

            # 1.3 - New labels for pixels at the edge of current edge
            # compute near range variation
            delta_range = self.compute_range_variation_between_tiles(tile_num, tile_num+1)
            new_labels_subset = self.gather_edge_entities(tile_idx1, tile_idx2,  delta_range)
            logger.debug("Range variation between tile is %d pixels" % delta_range)

            # 1.4 - Link old labels to new ones
            self.label_matching(tile_idx1, tile_idx2, new_labels_subset)
            
            logger.debug("")

        # 2 - Deal with edge at beginning or end of pass
        if (self.labels == 0).any():
            
            # 2.1 - Get tiles with unprocessed labels (still initialized to zero)
            tiles_to_process = np.unique(self.tile_index[np.where(self.labels == 0)])

            for tile in tiles_to_process:

                # 2.2 - Get associated tile index
                tile_idx = np.where(self.tile_index == tile)[0]

                # 2.3 - Get old labels from PGE_LakeTile
                old_labels = np.unique(self.edge_label[tile_idx][np.where(self.labels[tile_idx] == 0)])

                # 2.4 - Compute a global new label
                new_labels = np.arange(old_labels.size) + np.max(self.labels) + 1

                # 2.5 - Update global labels
                for idx in range(old_labels.size):
                    self.labels[tile_idx[np.where(self.edge_label[tile_idx] == old_labels[idx])]] = new_labels[idx]

        self.labels = self.labels.astype('int')
        # 3 - Relabel Lake Using Segmentation Heigth
        # For each label : check if only one lake is in each label and relabels if necessary

        # 4.0. Check if lake segmentation following height needs to be run
        seg_method = cfg.getint('CONFIG_PARAMS', 'SEGMENTATION_METHOD')
        min_object_size = cfg.getfloat('CONFIG_PARAMS', 'MIN_SIZE') * 1e6

        # If SEGMENTATION_METHOD is 0, function unactivated
        if seg_method == 0:
            logger.debug("Lake segmentation following height unactivated")
        # 4.1. If SEGMENTATION_METHOD not 0, relabel lake following segmentation height
        else:

            labels_tmp = np.zeros(self.labels.shape, dtype=self.labels.dtype)
            labels, count = np.unique(self.labels, return_counts=True)

            for i, label in enumerate(labels):
                idx = np.where(self.labels == label)
                subset_pixel_area = self.pixel_area[idx]

                if np.sum(subset_pixel_area) > min_object_size:
                    min_rg = min(self.get_range_of_lake(idx))
                    min_az = min(self.get_azimuth_of_lake(idx))
                    subset_range = self.get_range_of_lake(idx) - min_rg
                    subset_azimuth = self.get_azimuth_of_lake(idx) - min_az
                    subset_height = self.height[idx]

                    relabel_obj = my_segmentation.relabel_lake_using_segmentation_heigth(subset_range, subset_azimuth,
                                                                                         subset_height, subset_pixel_area,
                                                                                         min_object_size, seg_method)

                    labels_tmp[self.labels == label] = np.max(labels_tmp) + relabel_obj
                else:
                    labels_tmp[self.labels == label] = np.max(labels_tmp) + 1

            self.labels = labels_tmp

    # ----------------------------------------

    def gather_edge_entities(self, in_tile_idx1, in_tile_idx2, delta_range):
        """
        This function gives new labels for pixels at current tile edge.
            1. Pixels within a buffer around tile edge are selected.
            2. A water mask is computed
            3. The mask is labeled
            
        :param in_tile_idx1: indices of pixels edge of tile 1
        :type in_tile_idx1: 1D-array of int
        :param in_tile_idx2: indices of pixels edge of tile 2
        :type in_tile_idx2: 1D-array of int
        :param delta_range: variation of near_range between tiles 1 and 2
        :type delta_range: int
        
        :return: out_new_labels_subset = new labels given to pixels of edge entities
        :rtype: 1D-array of int
        """
        logger = logging.getLogger(self.__class__.__name__)

        # 1 - Pixels at top / bottom of tile are loaded
        rg, az = self.select_edge_pixels(in_tile_idx1, in_tile_idx2, delta_range)

        # 2 - Equivalent matrix size in azimuth and range
        # 2.1 - Equivalent matrix size in azimuth and range
        nb_pix_range = max(rg) + 1
        nb_pix_azimuth = max(az) + 1
        # 2.2 - Compute water mask over the subset0
        water_mask = my_tools.compute_bin_mat(nb_pix_range, nb_pix_azimuth, rg, az)

        # 3 - Label entities over the subset of PixC at the edge tile
        sep_entities, nb_obj = my_tools.label_region(water_mask)

        # 4 - Convert into 1D-array
        out_new_labels_subset = my_tools.convert_2d_mat_in_1d_vec(rg, az, sep_entities)

        logger.debug("%d separate entities located at the edge tile" % np.unique(out_new_labels_subset).size)

        return out_new_labels_subset

    # ----------------------------------------

    def select_edge_pixels(self, in_tile_idx1, in_tile_idx2, delta_range):
        """
        This function selects pixels at top and bottom of tile 1 and 2
        
        :param in_tile_idx1: indices of pixels edge of tile 1
        :type in_tile_idx1: 1D-array of int
        :param in_tile_idx2: indices of pixels edge of tile 2
        :type in_tile_idx2: 1D-array of int
        :param delta_range: variation of near_range between tiles 1 and 2
        :type delta_range: int
        
        :return: out_rg = range indices of pixels at the edge
        :rtype: 1D-array of int
        :return: out_az = azimuth indices of pixels at the edge
        :rtype: 1D-array of int
        """

        # 1 - Distinguish top / bottom edge
        rg1, az1 = self.get_edge_pixels(in_tile_idx1, "top")
        rg2, az2 = self.get_edge_pixels(in_tile_idx2, "bottom")

        # # 2 - Concatenate pixels range in a numpy array
        if delta_range > 0:
            rg = np.concatenate((rg1 + delta_range, rg2))
        else:
            rg = np.concatenate((rg1, rg2 - delta_range))

        # 3 - Concatenate pixels azimuth in a numpy array
        az = np.concatenate((az1, az2 + max(az1)+1))

        # 4 - Reduce range and azimuth values in order to reduce the size of generated water mask
        out_az = az - min(az)
        out_rg = rg - min(rg)

        # Return reduced range and azimuth indices
        return out_rg, out_az

    # ----------------------------------------

    def get_edge_pixels(self, in_tile_idx, in_edge_loc_str):
        """
        The function returns range and azimuth of pixels located within a buffer around the tile edge specified in edge_loc_str.
            
        :param in_tile_idx: indices of pixels edge of tile
        :type in_tile_idx: 1D array of int
        :param in_edge_loc_str: edge location = "top" or "bottom"
        :type in_edge_loc_str: string
        
        :return: range and azimuth indices of pixels
        :rtype: 1D array of int
        """
        logger = logging.getLogger(self.__class__.__name__)
        idx_edge_buf = None
        
        # Associated location
        if in_edge_loc_str == "bottom":
            # bottom of tile, get indices of all pixels with azimuth zero
            idx_edge_buf = np.where(self.azimuth_index[in_tile_idx] == 0)[0]
            
        elif in_edge_loc_str == "top":
            # top of tile, get indices of all pixels with maximal azimuth
            az_max = max(self.azimuth_index[in_tile_idx])
            idx_edge_buf = np.where(self.azimuth_index[in_tile_idx] == az_max)[0]

        else:
            logger.error("in_edge_loc_str input variable has to be 'top' or 'bottom'")

        return self.range_index[in_tile_idx][idx_edge_buf], self.azimuth_index[in_tile_idx][idx_edge_buf]

    # ----------------------------------------

    def label_matching(self, in_tile_idx1, in_tile_idx2, in_new_labels_subset):
        """
        This function matches labels computed in LakeTile_Edge file with labels computed in the gatherEdgeEntities() function.
        Pixels belonging to the same entity but cut by the tiling process are gathered with a single new label.
            
        :param in_tile_idx1: indices of pixels in tile 1
        :type in_tile_idx1: 1D-array of int
        :param in_tile_idx2: indices of pixels in tile 2
        :type in_tile_idx2: 1D-array of int
        :param in_new_labels_subset: new labels computed in gatherEdgeEntities()
        :type in_new_labels_subset: 1D-array of int
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Matching old labels from LakeTile_Edge with new labels")

        # 1 - Get old labels computed in LakeTile_Edge
        old_labels_subset1, old_labels_subset2 = self.select_edge_labels(in_tile_idx1, in_tile_idx2)

        # 2 - Get new labels
        # NB: the first part of in_new_labels_subset contains labels of tile 1, the end of the array contains labels of tile 2
        in_new_labels_subset1 = in_new_labels_subset[:old_labels_subset1.size]
        in_new_labels_subset2 = in_new_labels_subset[old_labels_subset1.size:]

        # the structure of in_new_labels_subset{1-2} and old_labels_subset{1-2} is the exactly the same,
        # it contains the labels of all pixels within the subset defined the azimuth buffer

        in_new_labels_subset_unique = np.unique(in_new_labels_subset)

        # correspondance contains a list of tuples. Each tuple will correspond to a new entity and will contains the old labels of tiles 1 and 2
        correspondance = []

        # Matching the labels of the subset of tile 1 and 2
        for new_l in in_new_labels_subset_unique:
            # get old labels (lake tile labels) of tiles 1 et 2
            old_l1 = np.unique(old_labels_subset1[np.where(in_new_labels_subset1 == new_l)])
            old_l2 = np.unique(old_labels_subset2[np.where(in_new_labels_subset2 == new_l)])

            # Most current case: at one label of tile1 corresponds one label of tile 2
            if old_l1.size == 1 and old_l2.size == 1:
                # appending to correspondance a tuple with old label 1 and old label 2
                correspondance.append((old_l1[0], old_l2[0]))

            # At one label of tile 1 do not corresponds a label of tile 2. 
            # The case happens when the lake is entierly located at the boundary of on tile bot does not cross the border.
            elif old_l2.size == 0:
                for idx1 in np.arange(old_l1.size):
                    correspondance.append((old_l1[idx1], None))

            # At one label of tile 2 do not corresponds a label of tile 1.
            elif old_l1.size == 0:
                for idx2 in np.arange(old_l2.size):
                    correspondance.append((None, old_l2[idx2]))

            # Case that rarely occurs: the lake meanders along the border between two tiles, in this case severals labels 
            # from tile1 matches with several labels of tile2.
            else:
                for idx1 in np.arange(old_l1.size):
                    for idx2 in np.arange(old_l2.size):
                        correspondance.append((old_l1[idx1], old_l2[idx2]))

        # To give an explicit example, labels of tile1 of replaced by letters belonging to ['a', ..] and labels of tile2 
        # are replaces by lettres belonging to [ 'r', ...]
        # correspondance: [(a, r), (a, s), (b, s), (b, t), (c, u), (d, u)]
        # labels a, b and r, s, t belong to the same entity, labels c, d, u belong to a separate entity.

        # The fonction lib.matchLabels returns reorganised correspondance labels:
        # label_matched: [([a, b], [r, s, t]), ([c, d], [u])]
        label_matched = match_labels(correspondance)

        logger.debug(" %d labels of first tile matched with %d labels of second tile into %d entities" % (np.unique(old_labels_subset1).size,
                                                                                                          np.unique(old_labels_subset2).size,
                                                                                                          len(label_matched)))

        # need new labels for global lake_sp processings
        unique_label = np.arange(len(label_matched)).astype('int') + max(self.labels) + 1
        # local labels from half edge tiles 1 and 2 are moved into new labels global only for Lake_sp processings
        # Ex: label_matched: [([a, b],[r,s,t]),([c,d],[u])]
        # The for loop iterates over entities at current tile edge
        for idx, (label_tile_1, label_tile_2) in enumerate(label_matched):
            # l contains the tuple corresponding to the idx^th entity
            # Ex: l: ([a, b], [r, s, t])

            # new_labels[idx] contains a label specefic for lake_SP processings
            new_label = unique_label[idx]
            if label_tile_1:
                for old_l1 in label_tile_1:
                    # At one label of tile 2 do not corresponds a label of tile 1. The case happens when the lake is entierly located at
                    # the boundary of on tile bot does not cross the border.
                    # In this case, old_l1 is setted to None, then, it is not processed.
                    if old_l1:
                        # Get label of entity already computed in the case of a lake covering more than two tiles
                        labels_concerned = np.unique(self.labels[in_tile_idx1][np.where(np.logical_and(self.edge_label[in_tile_idx1] == old_l1,
                                                                                                       self.edge_loc[in_tile_idx1] == 2))])
                        # Deleting label = 0 as those label are not already computed
                        labels_concerned = np.delete(labels_concerned, np.where(labels_concerned == 0))
                        # For each already processed and more than two tile lake, update the global labels
                        for label_to_relabel in labels_concerned:
                            self.labels[np.where(self.labels == label_to_relabel)] = new_label

                        # Set global label
                        self.labels[in_tile_idx1[np.where(self.edge_label[in_tile_idx1] == old_l1)]] = new_label

            if label_tile_2:
                for old_l2 in label_tile_2:
                    # At one label of tile 1 do not corresponds a label of tile 2. The case happens when the lake is entierly 
                    # located at the boundary of on tile bot does not cross the border.
                    # In this case, old_l2 is setted to None, then, it is not processed.
                    if old_l2:
                        # Get label of entity already computed in the case of a lake covering more than two tiles
                        labels_concerned = (np.unique(self.labels[in_tile_idx2][np.where(np.logical_and(self.edge_label[in_tile_idx2] == old_l2,
                                                                                                        self.edge_loc[in_tile_idx2] == 2))]))
                        # Deleting label = 0 as those label are not already computed
                        labels_concerned = np.delete(labels_concerned, np.where(labels_concerned == 0))

                        # For each already processed and more than two tile lake, update the global labels
                        for label_to_relabel in labels_concerned:
                            self.labels[np.where(self.labels == label_to_relabel)] = new_label

                        # Set global label
                        self.labels[in_tile_idx2[np.where(self.edge_label[in_tile_idx2] == old_l2)]] = new_label

            logger.debug("Labels %s of tile N matched with labels %s of tiles N+1 into new labels %s" % (str(label_tile_1), str(label_tile_2), str(new_label)))

        nb_edge_entities = np.unique(np.concatenate((self.labels[in_tile_idx1], self.labels[in_tile_idx2]))).size

        logger.debug(" %d separate entities are located at tile edge" % nb_edge_entities)

    # ----------------------------------------

    def select_edge_labels(self, in_tile_idx1, in_tile_idx2):
        """
        This function selects old labels from LakeTile_Edge at top and bottom of tile 1 and 2

        :param in_tile_idx1: indices of edge pixels in tile 1
        :type in_tile_idx1: 1D-array of int
        :param in_tile_idx2: indices of edge pixels in tile 2
        :type in_tile_idx2: 1D-array of int
        
        :return: labels of edge pixels at edge 1 and 2
        :rtype: 1D-array of int
        """

        label1 = self.get_edge_labels(in_tile_idx1, "top")
        label2 = self.get_edge_labels(in_tile_idx2, "bottom")

        return label1, label2

    # ----------------------------------------

    def get_edge_labels(self, in_tile_idx, in_edge_loc_str):
        """
        This function returns the LakeTile_Edge labels of pixels within the buffer zone
        
        :param in_tile_idx: indices of pixels at the edge of tile
        :type in_tile_idx: 1D array of int
        :param in_edge_loc_str: edge location = "top" or "bottom"
        :type in_edge_loc_str: string
        
        :return: LakeTile_Edge labels of pixels within the buffer zone
        :rtype: 1D-array of int
        """
        logger = logging.getLogger(self.__class__.__name__)
        idx_edge_buf = None
        
        # Associated location
        if in_edge_loc_str == "bottom":
            # Bottom of tile: get indices of all pixels with azimuth zero
            idx_edge_buf = np.where(self.azimuth_index[in_tile_idx] == 0)[0]
            
        elif in_edge_loc_str == "top":
            # Top of tile, get indices of all pixels with maximal azimuth
            az_max = max(self.azimuth_index[in_tile_idx])
            idx_edge_buf = np.where(self.azimuth_index[in_tile_idx] == az_max)[0]

        else:
            logger.error("in_edge_loc_str input variable has to be 'top' or 'bottom'")

        return self.edge_label[in_tile_idx[idx_edge_buf]]
        
    # ----------------------------------------

    def get_azimuth_of_lake(self, in_indices):
        """
            This function returns a re-computed azimuth index in order to have a continous azimuth along a lake at the edge of tiles

        :param in_indices: indices of pixels of a lake
        :type in_indices: 1D-array of int

        :return: recomputed azimuth_idx of the lake
        :rtype: 1D-array of int
        """

        tiles = np.unique(self.tile_index[in_indices])
        lake_tile_idx = self.tile_index[in_indices]
        lake_azimuth_idx = self.azimuth_index[in_indices]

        for tile in tiles:
            lake_azimuth_idx[np.where(lake_tile_idx > tile)] = lake_azimuth_idx[np.where(lake_tile_idx > tile)] +\
                                                               max(self.azimuth_index[in_indices][np.where(lake_tile_idx == tile)]) + 1

        return lake_azimuth_idx

    # ----------------------------------------

    def get_range_of_lake(self, in_indices):
        """
            This function returns a re-computed range index in order to have a continous range across a lake at the edge of tiles

        :param in_indices: indices of pixels of a lake
        :type in_indices: 1D-array of int

        :return: recomputed range_idx of the lake
        :rtype: 1D-array of int
        """

        out_range = []
        range_tmp = self.range_index[in_indices]
        tile_index_tmp = self.tile_index[in_indices]
        concerned_tiles, counts = np.unique(tile_index_tmp, return_counts=True)
        d_rg = 0
        for i, tile_num in enumerate(concerned_tiles):

            if i == 0:
                out_range += list(range_tmp[np.where(tile_index_tmp == tile_num)])
            else:
                d_rg += self.compute_range_variation_between_tiles(concerned_tiles[i-1], concerned_tiles[i])
                out_range += list(range_tmp[np.where(tile_index_tmp == tile_num)] - d_rg)
        if min(out_range) != 0:
            out_range = out_range + min(out_range)

        return np.array(out_range)


    def get_majority_pixels_tile_ref(self, in_lake_tile_label):
        """
        This fuction returns the tile reference of the tile containing the larger number of pixels with the given label.
            
        :param in_lake_tile_label: labels of lake to process
        :type in_lake_tile_label: int
        
        :return: tile reference
        :rtype: string
        """
        
        # 1 - Get unique values and counts of tiles for pixels with edge_label in_lake_tile_label
        unique, counts = np.unique(self.tile_index[np.where(self.edge_label == in_lake_tile_label)], return_counts = True)

        # 2 - Get tile ref corresponding to the max number of pixels
        out_tile_max_pix = str(unique[np.where(counts == max(counts))][0]).rjust(3, str('0')) + self.swath_side

        return out_tile_max_pix
        
    # ----------------------------------------

    def get_lake_tile_label(self, in_new_label):
        """
        This function is designed to retrieve old labels of PGE_LakeTile. 
        The given new label corresponds to a global label, corresponding to several old labels.
        The old label involving the largest number of pixels is return.
            
        :param in_new_label: global new label
        :type in_new_label: int
        
        :return: LakeTile_Edge label involving the largest number of pixels
        :rtype: string
        """

        # 1 - Get tiles concerned by the current new label
        tiles_concerned = np.unique(self.tile_index[np.where(self.labels == in_new_label)])

        nb_max_pix = 0
        out_final_label = 0

        for tile in tiles_concerned:

            # 2 - Get indices of in_new_label pixels in tile
            label_idx = np.where(self.labels == in_new_label)

            # 3 - Get lake_tile label value and number of pixels
            # unique, count = np.unique(self.edge_label[tile_idx], return_counts=True)
            unique, idx = np.unique(self.edge_label[label_idx], return_inverse=True)
            count = np.bincount(idx)
            np_max_pix_tmp = np.max(count)
            label_max_pix_tmp = unique[np.where(count == np_max_pix_tmp)][0]

            if np_max_pix_tmp > nb_max_pix:
                nb_max_pix = np_max_pix_tmp
                out_final_label = label_max_pix_tmp

        # 4 - Returns the lake tile label involving the largest number of pixels
        return str(out_final_label)

    # ----------------------------------------

    def compute_range_variation_between_tiles(self, tile_num1, tile_num2):
        """
        This function is designed to compute the variation of the first pixels in range between tiles 1 and 2.

        :param tile_num1: number of tile 1
        :type tile_num1: int
        :param tile_num2: number of tile 2
        :type tile_num2: int

        :return: variation of near range in pixels.
        :rtype: int
        """
        idx1 = self.tile_num.index(tile_num1)
        idx2 = self.tile_num.index(tile_num2)
        near_range1 = self.near_range[idx1]
        near_range2 = self.near_range[idx2]
        slant_range_spacing = self.slant_range_spacing[idx1]

        delta_near_range = np.rint((near_range1-near_range2)/slant_range_spacing)
        return int(delta_near_range)

    # ----------------------------------------

    def get_near_range(self, in_indices):
        """
        This function is designed to return a near range array, with the near range of each tiles corresponding to each pixels of in_indices.

        :param in_indices: indices of pixels of a lake
        :type in_indices: 1D-array of int

        :return: near rang of tiles covering indices of pixels of lales
        :rtype: 1D-array of float
        """
        near_range = []
        concerned_tiles, counts = np.unique(self.tile_index[in_indices], return_counts=True)
        for i, tile_num in enumerate(concerned_tiles):
            idx = self.tile_num.index(tile_num)
            near_range_tmp = self.near_range[idx]
            near_range += [near_range_tmp] * counts[i]

        return np.array(near_range)


#######################################


def match_labels(in_liste):
    """
    This function reorganise labels in order to group labels by entities
        
    :param IN_liste: ex: [(a, r), (a, s), (b, s), (b, t), (c, u), (d, u)]. Labels a, b and r, s, t belong to the same entity, labels c, d, u belong to a separate entity.
    
    :return: labels gathered by entities
             ex: [({a, b}, {r, s, t])), ({c, d}, {u})] <=> [([a, b], [r, s, t]), ([c, d], [u])]
    """
    labels_list = group_by_second(group_by_first(in_liste))
    labels_match_out = []
    for (set_label_tile1, set_label_tile2) in labels_list:

        if set_label_tile1 == {None}:
            for label in set_label_tile2:
                labels_match_out.append((set_label_tile1, {label}))
        elif set_label_tile2 == {None}:
            for label in set_label_tile1:
                labels_match_out.append(({label}, None))
        else:
            set_label_tile1.discard(None)
            set_label_tile2.discard(None)
            labels_match_out.append((set_label_tile1, set_label_tile2))

    return labels_match_out


def group_by_first(in_list):
    """
    This function take a list of tuples. The list is grouped by the first element of tuple.
        
    :param IN_list: la list of tuple. Ex: [(a, r), (b, r),  (c, s), (c, t)]
    
    :return out_list: ex: [(a,{r}), (b,{r}), (c,{s,t})]
    """

    l1 = []
    l2 = []

    for (ind_i, ind_j) in in_list:

        if ind_i is None:
            l1.append(ind_i)
            l2.append({ind_j})
        else :
            if ind_i in l1:
                idx = l1.index(ind_i)
                l2[idx].add(ind_j)
            else :
                l1.append(ind_i)
                l2.append({ind_j})

    out_list = []
    for (elem1, elem2) in zip(l1, l2):
        if len(elem2) > 1 and None in elem2:
            elem2.remove(None)
        out_list.append((elem1, elem2))

    out_list = [(l1[i], l2[i]) for i in range(len(l1))]

    return out_list


def group_by_second(in_list):
    """
    This function take a dictionary. The dictionary is grouped by second elements and returned as a list of tuple of set.
        
    :param in_list: result of group by first function: [(a,{r}), (b,{r}), (c,{s,t})]
    
    :return: a list of tuple of set. ex: [({a,b},{r}), ({c},{s,t})]
    """

    # Init
    out_results = []
    
    for (in_elem1, in_elem2) in in_list:
        has_been_found = False
        if in_elem2 == {None}:
            has_been_found = False
        else :
            for ind, (out_elem1, out_elem2) in enumerate(out_results):

                if checklist_inter(in_elem2, (out_elem1, out_elem2)):
                    if in_elem1 is None:
                        in_elem2 = in_elem2.difference(checklist_inter(in_elem2, (out_elem1, out_elem2)))
                    else:
                        has_been_found = True
                        out_results[ind] = merge( (out_elem1, out_elem2), in_elem1, in_elem2)

        if not has_been_found:
            out_results.append(({in_elem1}, in_elem2))

    return out_results


def checklist_inter(in_value, in_result):
    """
    Check if an element of value is contained in result
        
    :param in_value: a set of elements ex: {r}
    :param in_result: a tuple of set ex: ({a}, {r})
    
    :return: {r}
    """
    return in_result[1].intersection(in_value)


def merge(in_result, in_key, in_value):
    """
    Merge value with element of tuple key in result.
        
    :param in_result: a tuple of set ex: ({a}, {r})
    :param in_key: key of first element of tuple. ex: b
    :param in_value: set to merge ex: {r}
    
    :return: merged tuple of set. ex:{a, b]), {r])
    """
    return (in_result[0].union({in_key}), in_result[1].union(in_value))
