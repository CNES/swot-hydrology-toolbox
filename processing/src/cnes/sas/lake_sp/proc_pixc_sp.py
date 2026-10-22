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
# VERSION:5.0.0:DM:#91:2022/12/08:Poursuite industrialisation
# VERSION:5.0.1:DM:#91:2023/02/24:Poursuite industrialisation
# VERSION:5.1.0:DM:#91:2023/03/31:Poursuite industrialisation
# VERSION:5.2.0:DM:#91:2023/05/23:Poursuite industrialisation
# VERSION:5.3.0:DM:#91:2023/10/04:Poursuite industrialisation
# VERSION:6.1.0:DM:#91:2024/09/20:Poursuite industrialisation
# VERSION:6.2.0:DM:#91:2024/11/29:Poursuite industrialisation
# VERSION:6.3.1:FA:#855:2026/02/26: L2_HR_LAKE_SP fails on many CT (ex C024P232CT2) with message "ValueError: not enough values to unpack (expected 8, got 6)"
# VERSION:6.3.1:FA:#831:2026/02/26:SSALTO-330 [PGE][SAS LAKES] Mise à jour des TU (SSALPGE-14313)
# VERSION:6.3.2:DM:#900:2026/07/06: [PGE] [LAKES] Ajout d'une interface niveau moyen des lacs
# VERSION:6.3.2:DM:#900:2026/07/06:SSALTO-565 [PGE] [LAKES] Ajout d'une interface coeff NN pour calcul ice_dyn_f
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

import cnes.common.service_config_file as service_config_file

import cnes.common.lib.my_timer as my_timer
import cnes.common.lib.my_netcdf_file as my_nc
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_segmentation as my_segmentation
import cnes.common.lib.my_variables as my_var
import cnes.common.lib_lake.lake_db as lake_db
import cnes.common.lib_lake.proc_pixc as proc_pixc


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
            - edge_index contains the L2_HR_PIXC initial pixels indices. \
            This information is needed to update the improved geoloc and tag fields of LakeTile_PIXCVec files.
            - pixel_at_window_edge flags the pixels corresponding to the edge of the data window in the range direction (=1; =0 otherwise).
            - ice_dyn_f contains the dynamical ice cover flag.
            - cross_track_refined contains the refined cross-track location of the pixels.

        This class processes all LakeTile_Edge files of a single swath to gather all water regions at the along-track edge of tile into lake water regions.

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
            - labels / 1D-array of int: array of new labels
            - labels_biglakes / list of int: array of new labels corresponding to biglakes edge patches
            - labels_to_process / 1D-array of int: array of labels to be processed by the LakeSP processing
            - is_boundary_pix / list of bool: if pixel belong to the first / last azimuth of single pass
            - big_lakeid / set of str: store big lakes covered by all the tiles
            - near_range / list of float: store the near range of each tile
            - polarization / list of str: store the radar signal polarization of each tile
            - slant_range_spacing / list of int: store the slant range sampling of each tile (supposed to be 0.75)
            - pixc_first_line_index / list of int: store the azimuth index corresponding to the first slant range interferogram line of each tile
            - pixc_last_line_index / list of int: store the azimuth index corresponding to the last slant range interferogram line of each tile
            - classif_dict / dict: selection (by True or False) of pixels wrt their classification, for different processes \
                + already init arrays useful for aggregate.py lib; \
                key=water/open_water/interior_water/water_edge/dark/low_coherence/mandatory/specular_ringing/not_specular_ringing/4hull/4wse/4area/Idw_in/Idw_in_no_dw/Idw/Idw_no_dw/Ide/I/I_no_dw
            - interferogram_flattened / 1D-array of complex: flattened interferogram
            - inundated_area / 1D-array of int: area of pixels covered by water
            - height_std_pix / 1D-array of float: height std
            - corrected_height / 1D-array of float: height corrected from geoid and other tide corrections
            
        - From L2_HR_LakeTile_Edge file:
            
            - Edge objects info:
                - edge_loc
                - edge_label
                - edge_index
                - pixel_at_window_edge
                - ice_dyn_f
                - cross_track_refined

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
                - sig0 / 1D array of float: sigma0
                - sig0_u / 1D array of float: sigma0 uncertainty (= variable named sig0_uncert in L2_HR_PIXC)
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
                - prev_xover / 1D array of float: time difference from the previous crossover (= variable named time_from_prev_xover in L2_HR_PIXC file)
                - next_xover / 1D array of float: time difference to the next crossover (= variable named time_to_next_xover in L2_HR_PIXC file)
                - wavelength / float: wavelength corresponding to the effective radar carrier frequency 
                - looks_to_efflooks / float: ratio between the number of actual samples and the effective number of independent \
                                                samples during spatial averaging over a large 2-D area

            - From PIXC/tvp group variables:
                - nadir_time[_tai] / 1D-array of float: observation UTC [TAI] time of each nadir pixel (= variable named time[tai] in L2_HR_PIXC file)
                - nadir_longitude / 1D-array of float: longitude of each nadir pixel (= variable named longitude in L2_HR_PIXC file)
                - nadir_latitude / 1D-array of float: latitude of each nadir pixel (= variable named latitude in L2_HR_PIXC file)
                - nadir_[x|y|z] / 1D-array of float: [x|y|z] cartesian coordinates of each nadir pixel (= variables named [x|y|z] in L2_HR_PIXC file)
                - nadir_[vx|vy|vz] / 1D-array of float: velocity vector of each nadir pixel in cartesian coordinates \
                (= variables named velocity_unit_[x|y|z] in L2_HR_PIXC file)
                - nadir_plus_y_antenna_[x|y|z] / 1D-array of float: position vector of the +y KaRIn antenna phase center in ECEF coordinates \
                (= variables named plus_y_antenna_[x|y|z] in L2_HR_PIXC file)
                - nadir_minus_y_antenna_[x|y|z] / 1D-array of float: position vector of the -y KaRIn antenna phase center in ECEF coordinates \
                (= variables named minus_y_antenna_[x|y|z] in L2_HR_PIXC file)
        """
        # Get instance of service config file
        self.cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        self.product_type = "SP"
        
        # 0 - Retrieve CLASSIF_LIST from configuration file
        self.classif_list = [int(i) for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_LIST").split(";")]

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
        self.laketile_labels = np.array(()).astype('object')  # Label in tile with tile reference (ex : '261R_23')
        self.edge_index = np.array(()).astype('int')  # Index in original L2_HR_PIXC product
        self.pixel_at_window_edge = np.array(()).astype('int')  # Flag pixels at the window edge in the range direction
        self.ice_dyn_f = np.array(()).astype('int')  # Dynamical ice cover flag
        self.cross_track_refined = np.array(())  # Refined cross-track location of the pixels
        
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
        self.sig0 = np.array(())
        self.sig0_u = np.array(())
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
        self.prev_xover = np.array(())
        self.next_xover = np.array(())
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
        self.pixc_metadata["big_lakeid"] = ""
        self.pixc_metadata["polarization"] = ""
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
        self.labels_biglakes = []
        self.labels_to_process = 0
        self.is_boundary_pix = []  # If pixel belongs to the first / last azimuth of single pass
        self.big_lakeid = set()
        self.near_range = []
        self.polarization = []
        self.slant_range_spacing = []
        self.pixc_first_line_index = []
        self.pixc_last_line_index = []
        self.list_classif_keys = ["water", "open_water", "interior_water", "water_edge", "dark", "low_coherence", 
                                  "mandatory", "specular_ringing", "not_specular_ringing",
                                  "4hull", "4wse", "4area", 
                                  "Idw_in", "Idw_in_no_dw", "Idw", "Idw_no_dw", "Ide", "I", "I_no_dw"]
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

        # 0 - Init variables
        in_laketile_edge_path_list.sort()
        # Store number of pixels in azimuth
        nb_pix_azimuth_in_full_pass = [0]
        # Store azimuth boundaries: beginning and end of pass + edges in case of missing tiles
        azimuth_boundary_in_full_pass = [0]

        for lake_tile_edge_path in in_laketile_edge_path_list:  # For each LakeTile_Edge file

            logger.debug(". Loading LakeTile_Edge file = %s" % lake_tile_edge_path)

            # 2 - Load data
            nb_pix_loaded, current_tile_number, \
                current_pixc_first_line_index, current_pixc_last_line_index, \
                current_pixc_first_line_is_start_hr, current_pixc_last_line_is_stop_hr, \
                current_big_lakeid, \
                current_near_range, current_polarization, current_slant_range_spacing = \
                self.load_laketile_edge_data(lake_tile_edge_path)

            # Stop process file if the tile doesn't correspond to the current continent
            if not current_tile_number:
                continue

            logger.debug("%d pixels loaded" % nb_pix_loaded)

            # 3 - Update number of pixels over the entire pass
            self.nb_pixels += nb_pix_loaded

            # 4 - Check if tile is consecutive to the previous tile
            if len(self.tile_num) >= 1:
                previous_tile_number = self.tile_num[-1]

                cpt = 0  # Counter of missing tiles
                while current_tile_number != previous_tile_number + 1:
                    # If current tile is not consecutive to previous tile, add an empty tile to tile_ref
                    logger.warning("Tile %d%s is missing => a fake tile is added" % (previous_tile_number + 1, self.swath_side))
                    
                    # 4.1 - Add fake tile info
                    self.tile_num.append(previous_tile_number + 1)
                    self.near_range.append(0)
                    self.polarization.append("")
                    self.slant_range_spacing.append(0)
                    self.pixc_first_line_index.append(current_pixc_first_line_index)
                    self.pixc_last_line_index.append(current_pixc_last_line_index)
                    previous_tile_number = self.tile_num[-1]

                    # 4.2 - Add an azimuth boundary for adjacent tile azimuth min and max
                    azimuth_boundary_in_full_pass.append(np.sum(nb_pix_azimuth_in_full_pass)-1)
                    nb_pix_azimuth_in_full_pass.append(current_pixc_last_line_index-current_pixc_first_line_index+1)
                    azimuth_boundary_in_full_pass.append(np.sum(nb_pix_azimuth_in_full_pass))

                    # Stop process if more than 50 tiles are missing
                    cpt += 1
                    if cpt > 50:
                        break

            # 5 - Add current tile information
            self.tile_index += [current_tile_number] * nb_pix_loaded
            self.tile_num.append(current_tile_number)
            self.big_lakeid = self.big_lakeid.union(set(current_big_lakeid.split(";")))
            self.near_range.append(current_near_range)
            self.polarization.append(current_polarization)
            self.slant_range_spacing.append(current_slant_range_spacing)
            self.pixc_first_line_index.append(current_pixc_first_line_index)
            self.pixc_last_line_index.append(current_pixc_last_line_index)
            if current_pixc_first_line_is_start_hr == 1:
                logger.debug("WARNING - Current tile is a start of an HR segment => features with pixels at top of previous tile are considered as partial")
                azimuth_boundary_in_full_pass.append(np.sum(nb_pix_azimuth_in_full_pass)-1)
            nb_pix_azimuth_in_full_pass.append(current_pixc_last_line_index-current_pixc_first_line_index+1)
            if current_pixc_last_line_is_stop_hr == 1:
                logger.debug("WARNING - Current tile is a stop of an HR segment => features with pixels at bottom of next tile are considered as partial")
                azimuth_boundary_in_full_pass.append(np.sum(nb_pix_azimuth_in_full_pass))

        # 6.1 - Add last azimuth index to list of azimuth indices of boundaries
        azimuth_boundary_in_full_pass.append(np.sum(nb_pix_azimuth_in_full_pass)-1)
        # 6.2 - Retrieve unique indices and sort
        azimuth_boundary_in_full_pass = np.sort(np.unique(np.array(azimuth_boundary_in_full_pass)))
        
        # 7 - Find if pixels belong to the beginning or end of pass or missing tile egde

        # 7.1 - Compute azimuth index of the beginning of each tile
        full_pass_tile_azimuth_index = np.cumsum(nb_pix_azimuth_in_full_pass)

        # 7.2 - Compute equivalent azimuth indices of each edge pixel within the full pass
        self.tile_index = np.array(self.tile_index)  # Convert to numpy array
        self.azimuth_index_in_full_pass = np.zeros((self.nb_pixels)).astype('int')
        for i, tile in enumerate(self.tile_num):
            idx = np.where(self.tile_index == tile)[0]
            self.azimuth_index_in_full_pass[idx] = self.azimuth_index[idx].astype('int') \
                                                    - self.pixc_first_line_index[i] \
                                                    + full_pass_tile_azimuth_index[i]

        # 7.3 - Fill boundary array if azimuth of pixels is a pass boundary
        self.is_boundary_pix = np.zeros((self.nb_pixels))
        for az_boundary in azimuth_boundary_in_full_pass:
            self.is_boundary_pix[self.azimuth_index_in_full_pass == az_boundary] = 1

        # 8 - Init labels of water regions
        self.labels = np.zeros((self.nb_pixels), dtype=int)

        logger.debug("%d PixC loaded for current swath" % self.nb_pixels)
        if len(self.big_lakeid) == 0:
            logger.debug("NO big lake considered for this swath")
        else:
            logger.debug("{} big lakes considered for this swath = {}".format(len(self.big_lakeid), self.big_lakeid))

    # ----------------------------------------

    def load_laketile_edge_data(self, in_laketile_edge_filename):
        """
        This function loads NetCDF information.
        
        :param in_laketile_edge_filename: full path of the NetCDF file to load
        :type in_laketile_edge_filename: string
        
        :return: out_nb_pix = number of pixels loaded
        :rtype: out_nb_pix = integer
        :return: out_tile_number = number of loaded tile
        :rtype: out_tile_number = int
        :return: out_pixc_first_line_index = azimuth index corresponding to the first slant range interferogram line that is inside the pixel cloud tile
        :rtype: out_pixc_first_line_index = int
        :return: out_pixc_last_line_index = azimuth index corresponding to the last slant range interferogram line that is inside the pixel cloud tile
        :rtype: out_pixc_last_line_index = int
        :return: out_pixc_first_line_is_start_hr = =True if the pixel cloud tile contains the start of a new HR segment. =False otherwise
        :rtype: out_pixc_first_line_is_start_hr = bool
        :return: out_pixc_last_line_is_stop_hr = =True if the pixel cloud tile contains the end of an HR segment. =False otherwise
        :rtype: out_pixc_last_line_is_stop_hr = bool
        :return: out_big_lakeid = list of PLD lake_id of big lakes covered by the tile
        :rtype: out_big_lakeid = str
        :return: out_near_range = slant range for the first image pixel
        :rtype: out_near_range = int
        :return: out_polarization = radar polarization for the tile
        :rtype: out_polarization = str
        :return: out_slant_range_spacing = range spacing corresponding to the 200 MHz sampling frequency
        :rtype: out_slant_range_spacing = int
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Open input NetCDF file in reading mode
        pixc_edge_reader = my_nc.MyNcReader(in_laketile_edge_filename)

        # 2 - Get and check tile references (cycle, pass, swath)
        out_tile_number = int(pixc_edge_reader.get_att_value("tile_number"))

        current_cycle_num = int(pixc_edge_reader.get_att_value("cycle_number"))
        if current_cycle_num != self.cycle_num:
            logger.error("Cycle of tile %d doesn't match with SP product %d" %(current_cycle_num, self.cycle_num))

        current_pass_number = int(pixc_edge_reader.get_att_value("pass_number"))
        if current_pass_number != self.pass_num:
            logger.error("Pass of tile %d doesn't match with SP product %d" %(current_pass_number, self.pass_num))

        current_swath_side = str(pixc_edge_reader.get_att_value("swath_side"))
        if current_swath_side != self.swath_side:
            logger.error("Swath of tile %s doesn't match with PixCEdgeSwath %s" %(current_swath_side, self.swath_side))
        current_date = parser.parse(pixc_edge_reader.get_att_value("time_granule_start"))
        if not self.date:
            self.date = current_date
        # Stop the process if acquisition dates delta > 24h
        if (current_date - self.date) > datetime.timedelta(days=1):
            logger.error("Input LakeTile_Edge file doesn't correspond to the same aquisition date")

        current_continent_id = str(pixc_edge_reader.get_att_value("continent_id"))
        if not self.continent_id in current_continent_id:
            # If cur_continent_id doesn't belong to the EDGE SP product, do not add pixc info
            logger.error("Input LakeTile_Edge %s file doesn't correspond to the same continent %s" % (in_laketile_edge_filename, self.continent_id))
            retour = None, None, None, None, None, None, None, None, None, None

        else:
    
            # 3 - Update tile_poly with new tile
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
    
            # 4 - Initialization of object variables if not already done
            if not self.tile_num:
                
                # 4.1 - Values for metadata
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
                
                # 4.2 - Other variables
                self.nb_pix_range = int(pixc_edge_reader.get_att_value("interferogram_size_range"))
                self.wavelength = np.float64(pixc_edge_reader.get_att_value("wavelength"))  # Wavelength
                # Ratio between the number of actual samples and the effective number of independent samples
                self.looks_to_efflooks = np.float64(pixc_edge_reader.get_att_value("looks_to_efflooks"))
    
            # 5 - Retrieve global attributes
            out_nb_pix = pixc_edge_reader.get_dim_value('points')
            out_pixc_first_line_index  = pixc_edge_reader.get_att_value("pixc_first_line_index")
            out_pixc_last_line_index  = pixc_edge_reader.get_att_value("pixc_last_line_index")
            out_pixc_first_line_is_start_hr = bool(pixc_edge_reader.get_att_value("pixc_first_line_is_start_hr"))
            out_pixc_last_line_is_stop_hr = bool(pixc_edge_reader.get_att_value("pixc_last_line_is_stop_hr"))
            #code commente out_nb_pix_azimuth = pixc_edge_reader.get_att_value("interferogram_size_azimuth")
            out_big_lakeid = pixc_edge_reader.get_att_value("big_lakeid")
            out_near_range = pixc_edge_reader.get_att_value("near_range")
            out_polarization = pixc_edge_reader.get_att_value("polarization")
            out_slant_range_spacing = pixc_edge_reader.get_att_value("nominal_slant_range_spacing")
            
            # 6 - Update metadata from LakeTile_Edge info
            self.update_metadata_from_edge_file(pixc_edge_reader)
            
            # 7 - Update vectors if there are pixels in the current LakeTile_Edge file
            if out_nb_pix > 0:
    
                # 7.1 - Add edge objects info
                tmp_edge_label = pixc_edge_reader.get_var_value("edge_label")
                self.edge_label = np.concatenate((self.edge_label, tmp_edge_label))
                tmp_laketile_labels = np.empty(len(tmp_edge_label), dtype='object')
                tmp_laketile_labels[:] = ["%d%s_%d" % (out_tile_number, self.swath_side, l) for l in tmp_edge_label]
                self.laketile_labels = np.concatenate((self.laketile_labels, tmp_laketile_labels))
                self.edge_index = np.concatenate((self.edge_index, pixc_edge_reader.get_var_value("edge_index")))
                self.edge_loc = np.concatenate((self.edge_loc, pixc_edge_reader.get_var_value("edge_loc")))
                self.pixel_at_window_edge = np.concatenate((self.pixel_at_window_edge, pixc_edge_reader.get_var_value("pixel_at_window_edge")))
                self.ice_dyn_f = np.concatenate((self.ice_dyn_f, pixc_edge_reader.get_var_value("ice_dyn_f")))
                self.cross_track_refined = np.concatenate((self.cross_track, pixc_edge_reader.get_var_value("cross_track_refined")))
    
                # 7.2 - Add variables from PIXC/pixel_cloud
                tmp_classif = pixc_edge_reader.get_var_value("classification")
                self.classif = np.concatenate((self.classif, tmp_classif))
                
                tmp_classification_qual = pixc_edge_reader.get_var_value("classification_qual")
                self.classification_qual = np.concatenate((self.classification_qual, tmp_classification_qual))
                tmp_geolocation_qual = pixc_edge_reader.get_var_value("geolocation_qual")
                self.geolocation_qual = np.concatenate((self.geolocation_qual, tmp_geolocation_qual))
                
                classif_dict = self.init_classif_dict(self.classif_list, out_nb_pix, tmp_classif, classif_qual=tmp_classification_qual)
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
                
                self.sig0 = np.concatenate((self.sig0, pixc_edge_reader.get_var_value("sig0")))
                self.sig0_u = np.concatenate((self.sig0_u, pixc_edge_reader.get_var_value("sig0_uncert")))
                
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
                
                self.prev_xover = np.concatenate((self.prev_xover, pixc_edge_reader.get_var_value("time_from_prev_xover")))
                self.next_xover = np.concatenate((self.next_xover, pixc_edge_reader.get_var_value("time_to_next_xover")))
    
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
                tmp_corrected_height = np.full(out_nb_pix, np.nan, dtype=np.float32)
                tmp_corrected_height[ind_valid_corr] = tmp_height[ind_valid_corr] \
                                                       - tmp_geoid[ind_valid_corr] \
                                                       - tmp_solid_earth_tide[ind_valid_corr] \
                                                       - tmp_pole_tide[ind_valid_corr] \
                                                       - tmp_load_tide_fes[ind_valid_corr]
                self.corrected_height = np.concatenate((self.corrected_height, tmp_corrected_height))
    
            # 8 - Close file
            pixc_edge_reader.close()
            
            retour = out_nb_pix, out_tile_number, \
                        out_pixc_first_line_index, out_pixc_last_line_index, \
                        out_pixc_first_line_is_start_hr, out_pixc_last_line_is_stop_hr, \
                        out_big_lakeid, \
                        out_near_range, out_polarization, out_slant_range_spacing

        return retour
    
    def update_metadata_from_edge_file(self, pixc_edge_reader):
        """
        Update metadata from LakeTile_Edge file
        
        :param pixc_edge_reader: reader of LakeTile_Edge NetCDF file
        :type pixc_edge_reader: my_netcdf_file.MyNcReader
        """
        
        # 1 - Global variables
        self.pixc_metadata["time_granule_start"] = min(self.pixc_metadata["time_granule_start"],
                                                        str(pixc_edge_reader.get_att_value("time_granule_start")))
        self.pixc_metadata["time_granule_end"] = max(self.pixc_metadata["time_granule_end"],
                                                      str(pixc_edge_reader.get_att_value("time_granule_end")))
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
        
        # 2 - Bounding box
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
        
    # ----------------------------------------

    def swath_global_relabeling(self):
        """
        Give a new label to water regions gathered at tile edges.
        """
        cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        timer = my_timer.Timer()
        timer.start()

        # 1 - Gather regions split by PIXC tiling for the whole swath
        labels_at_edges = []

        for tile_num1, tile_num2 in zip(self.tile_num, self.tile_num[1:]):  # Loop over tile edges
            tile_ref1 = "%d%s" % (tile_num1, self.swath_side)
            tile_ref2 = "%d%s" % (tile_num2, self.swath_side)
            logger.debug(" ***** Processing top of tile %s and bottom of tile %s *****" % (tile_ref1, tile_ref2))
            
            # 1.1 - Get indices of pixels processed at the current edge
            # Pixels related to regions at the top of tile_1 and crossing tile_1
            idx_tile1 = np.where(np.logical_and(self.tile_index == tile_num1,
                                                np.logical_or(np.logical_or(self.edge_loc == 1, self.edge_loc == 2), self.edge_loc == 11)))[0]
            # Pixels related to regions at the bottom of tile_2 and crossing tile_2
            idx_tile2 = np.where(np.logical_and(self.tile_index == tile_num2,
                                                np.logical_or(np.logical_or(self.edge_loc == 0, self.edge_loc == 2), self.edge_loc == 10)))[0]

            # Process only if there is more than one pixel at both tiles edges
            tmp_labels_at_edges = None
            if idx_tile1.size != 0 and idx_tile2.size != 0:
                logger.debug("%d pixels of tile %s will be linked with %d pixels of tile %s" % \
                                 (idx_tile1.size, tile_ref1, idx_tile2.size, tile_ref2))

                # 1.2 - Compute near range variation between both tiles
                delta_range = self.compute_range_variation_between_tiles(tile_num1, tile_num2)
                logger.debug("Range variation between both tiles is %d pixels" % delta_range)
                
                # 1.3 - Match LakeTile labels on both sides of the current edge
                # New labels are computed and stored in labels_at_edges
                idx1 = self.tile_num.index(tile_num1)
                idx2 = self.tile_num.index(tile_num2)
                tmp_labels_at_edges = self.gather_regions_at_edge(idx_tile1, idx_tile2,
                                                                  delta_range,
                                                                  self.pixc_last_line_index[idx1], 
                                                                  self.pixc_first_line_index[idx2])
            if tmp_labels_at_edges is not None:
                labels_at_edges += tmp_labels_at_edges
        
        # 1.4. Match all labels together, at the swath scale
        self.gather_regions_of_swath(labels_at_edges)
        logger.debug("%d features found in current swath" % len(np.unique(self.labels)))
        logger.debug("")

        # 2 - Compute lake segmentation on each edge region
        
        # 2.1 - Check if lake segmentation following height needs to be run
        seg_method = cfg.getint('CONFIG_PARAMS', 'SEGMENTATION_METHOD')
        min_object_size = cfg.getfloat('CONFIG_PARAMS', 'MIN_SIZE') * 1e6
        min_object_nb_pixc = cfg.getint('CONFIG_PARAMS', 'MIN_NB_PIXC')

        # If SEGMENTATION_METHOD is 0, function unactivated
        if seg_method == 0:
            logger.debug("Lake segmentation following height unactivated")
            
        # 2.2 - If SEGMENTATION_METHOD not 0, relabel lake following segmentation height
        else:
            logger.debug("Lake segmentation following height activated")

            labels_tmp = np.zeros(self.labels.shape, dtype=self.labels.dtype)
            labels, count = np.unique(self.labels, return_counts=True)
            
            labels_biglakes_tmp = []

            for i, label in enumerate(labels):
                
                # Specific case of big lake
                if label in self.labels_biglakes:
                    new_label = np.max(labels_tmp) + 1
                    labels_tmp[self.labels == label] = new_label
                    labels_biglakes_tmp.append(new_label)
                    continue
                    
                idx = np.where(self.labels == label)
                subset_pixel_area = self.pixel_area[idx]

                if np.sum(subset_pixel_area) > min_object_size:
                    
                    # Subset of range indices
                    subset_range = self.get_range_of_lake(idx)
                    min_rg = min(subset_range)
                    subset_range -= min_rg  # Reduce the indices
                    
                    # Subset of azimuth indices
                    subset_azimuth = self.azimuth_index_in_full_pass[idx]
                    min_az = min(subset_azimuth)
                    subset_azimuth -= min_az  # Reduce the indices
                    
                    # Subset of height
                    subset_height = self.height[idx]

                    # Compute segmentation height
                    relabel_obj = my_segmentation.relabel_lake_using_segmentation_height(subset_range, subset_azimuth,
                                                                                         subset_height, subset_pixel_area,
                                                                                         min_object_size, min_object_nb_pixc,
                                                                                         seg_method)

                    labels_tmp[self.labels == label] = np.max(labels_tmp) + relabel_obj
                else:
                    labels_tmp[self.labels == label] = np.max(labels_tmp) + 1

            self.labels = labels_tmp
            self.labels_biglakes = labels_biglakes_tmp
            
        # 3 - Clean labels
        self.clean_labels()
        
        logger.debug(timer.stop(message="swath_global_relabeling computation time : "))

    def compute_range_variation_between_tiles(self, tile_num1, tile_num2):
        """
        This function is designed to compute the variation of the first pixels in range between tiles 1 and 2.

        :param tile_num1: number of tile_1
        :type tile_num1: int
        :param tile_num2: number of tile_2
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

    def gather_regions_at_edge(self, idx_tile1, idx_tile2, delta_range, pixc_last_line_index1, pixc_first_line_index2):
        """
        Match labels at top of tile_1 with those at bottom of tile_2 

        :param idx_tile1: indices of pixels related to regions located at the "top" edge of (or crossing) tile_1
        :type idx_tile1: 1D-array of int
        :param idx_tile2: indices of pixels related to regions located at the "bottom" edge of (or crossing) tile_2
        :type idx_tile2: 1D-array of int
        :param delta_range: variation of near_range between tile_1 and tile_2
        :type delta_range: int
        :param pixc_last_line_index1: azimuth index corresponding to the last slant range interferogram line that is inside tile_1 (= azimuth_index "max")
        :type pixc_last_line_index1: int
        :param pixc_first_line_index2: azimuth index corresponding to the first slant range interferogram line that is inside tile_2 (= azimuth_index "min")
        :type pixc_first_line_index2: int

        :return: out_label_matching = list of list of old LakeTile labels corresponding to the same new label (for each new label); =None if there is None
        :rtype: out_label_matching = list of list of str
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 0.1 - Compute indices of pixels (within idx_tile1) at the very top edge of tile_1
        idx_edge_in_tile1 = np.where(self.azimuth_index[idx_tile1] == pixc_last_line_index1)[0]
        # 0.2 - Compute indices of pixels (within idx_tile2) at the very bottom edge of tile_2
        idx_edge_in_tile2 = np.where(self.azimuth_index[idx_tile2] == pixc_first_line_index2)[0]
        
        out_label_matching = None  # Init output list
        if (len(idx_edge_in_tile1) + len(idx_edge_in_tile2)) > 0:

            # 1 - Compute labels of regions crossing the tile edge
            # NB: consider only "1 range line pixels" at the edge for each tile
            new_labels_at_edge = self.compute_labels_at_edge(idx_tile1[idx_edge_in_tile1],
                                                             idx_tile2[idx_edge_in_tile2], 
                                                             delta_range)
    
            # 2 - Retrieve LakeTile labels at the very edge of each tile
            label1 = self.laketile_labels[idx_tile1][idx_edge_in_tile1]
            label2 = self.laketile_labels[idx_tile2][idx_edge_in_tile2]
            laketile_labels_at_edge = np.concatenate((label1, label2))
            
            # 3 - Link old LakeTile labels to new labels at edge
            out_label_matching = self.link_old_to_new_labels(laketile_labels_at_edge, new_labels_at_edge)
            
        else:
            logger.debug("No water feature is truly located at the tile edge")

        return out_label_matching

    def compute_labels_at_edge(self, idx_edge1, idx_edge2, delta_range):
        """
        This function gives new labels for pixels in a 1-pixel buffer around tile edge
            1. Compute the corresponding water mask
            2. Compute labels associated to these "water features"
            
        :param idx_edge1: indices of pixels at the very "top" edge of tile_1
        :type idx_edge1: 1D-array of int
        :param idx_edge2: indices of pixels at the very "bottom" edge of tile_2
        :type idx_edge2: 1D-array of int
        :param delta_range: variation of near_range between tiles 1 and 2
        :type delta_range: int
        
        :return: out_new_labels_at_edge = new labels given to pixels of edge water regions
        :rtype: out_new_labels_at_edge = 1D-array of int
        """
        logger = logging.getLogger(self.__class__.__name__)

        # 1.1 - Retrieve range and azimuth indices of pixels at the very top edge of tile_1
        rg1 = self.range_index[idx_edge1]
        az1 = self.azimuth_index_in_full_pass[idx_edge1]
        # 1.2 - Retrieve range and azimuth indices of pixels at the very bottom edge of tile_2
        rg2 = self.range_index[idx_edge2]
        az2 = self.azimuth_index_in_full_pass[idx_edge2]

        # 2 - Compute water mask over tile edges
        # 2.1 - Concatenate pixels of two tiles into a single azimuth and range array
        if delta_range > 0:
            rg = np.concatenate((rg1 + delta_range, rg2))
        else:
            rg = np.concatenate((rg1, rg2 - delta_range))
        az = np.concatenate((az1, az2))

        # 2.2 - Reduce range and azimuth values in order to reduce the size of generated water mask
        az -= min(az)
        rg -= min(rg)

        # 2.3 - Compute water mask over the subset0
        water_mask = my_tools.compute_bin_mat(max(rg) + 1, max(az) + 1, rg, az)

        # 3 - Label water regions over the subset of PixC at the edge tile
        sep_regions, nb_obj = my_tools.label_region(water_mask)

        # 4 - Convert into 1D-array
        out_new_labels_at_edge = my_tools.convert_2d_mat_in_1d_vec(rg, az, sep_regions)

        logger.debug("%d separate water regions located at the edge tile" % np.unique(out_new_labels_at_edge).size)

        return out_new_labels_at_edge

    def link_old_to_new_labels(self, in_old_labels, in_new_labels):
        """
        Link old LakeTile labels with new labels at edge

        :param in_old_labels: LakeTile labels for each pixel at the very edge of two consecutive tiles
        :type in_old_labels: 1D-array of int
        :param in_new_labels: new edge labels for each pixel at the very edge of two consecutive tiles
        :type in_new_labels: 1D-array of int

        :return: out_label_matching = list of list of old LakeTile labels corresponding to the same new label (for each new label)
        :rtype: out_label_matching = list of list of str
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        out_label_matching = []
        
        for new_label in np.unique(in_new_labels):
            
            # Retrieve LakeTile labels corresponding to new_label
            linked_laketile_labels = np.unique(in_old_labels[np.where(in_new_labels == new_label)])
            
            # Add to output list
            out_label_matching.append(list(linked_laketile_labels))
            
        logger.debug("%d matching have been found at this edge" % len(out_label_matching))
        return out_label_matching

    def gather_regions_of_swath(self, label_matching):
        """
        Reorganize all LakeTile labels at the scale of the whole swath using label matching computed by gather_regions_at_edge.

        :param label_matching: list of list of old LakeTile labels corresponding to the same new label (for each new label)
        :type label_matching: list of list of str
        """

        # 1 - Gather LakeTile labels across current swath
        laketile_labels_sorted_by_object = []
        # laketile_labels_sorted_by_object is a list of sets gathering laketile_labels into waterbody. ex :
        #   edge_correspondence :['261R_23'], ['261R_26', '262R_1'], ['261R_13', '262R_5'], ['262R_5', '264R_10']]
        #   laketile_labels_sorted_by_object :[{'261R_23'}, {'261R_26', '262R_1'}, {'261R_13', '262R_5', '264R_10'}]
        for obj_labels in label_matching:
            added = False
            for i, obj in enumerate(laketile_labels_sorted_by_object):
                if set(obj_labels) & obj:
                    if not added:
                        laketile_labels_sorted_by_object[i] = obj | set(obj_labels)
                        added = True

            if not added:
                laketile_labels_sorted_by_object.append(set(obj_labels))
            else :
                for i, obj_labels_1 in enumerate(laketile_labels_sorted_by_object):
                    for j, obj_labels_2 in enumerate(laketile_labels_sorted_by_object[i+1:]):
                        if obj_labels_1 & obj_labels_2:
                            laketile_labels_sorted_by_object[i] = set(obj_labels_1 | obj_labels_2)
                            laketile_labels_sorted_by_object.remove(obj_labels_2)

        # 2 - Add LakeTile labels to laketile_labels_correspondence if no pixel is found at tiles borders
        for laketile_label in list(np.unique(self.laketile_labels)):
            if laketile_label not in [l for obj_labels in laketile_labels_sorted_by_object for l in obj_labels]:
                laketile_labels_sorted_by_object.append(set([laketile_label]))

        # 3 - Relabel swath labels self.labels using laketile_labels_sorted_by_object
        for i, obj_labels in enumerate(laketile_labels_sorted_by_object):
            lakesp_label = i+1
            for l in obj_labels:
                idx = np.where(self.laketile_labels == l)[0]
                self.labels[idx] = lakesp_label 
                if self.edge_loc[idx[0]] >= 10:
                    self.labels_biglakes.append(lakesp_label)
        
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
        
    # ----------------------------------------

    def get_majority_pixels_tile_ref(self, in_pixc_index, in_lake_tile_label):
        """
        This function returns the tile reference of the tile containing the largest number of pixels with the given label.
        
        :param in_pixc_index: indices of pixels defining the current feature
        :type in_pixc_index: 1D-array of int
        :param in_lake_tile_label: labels of lake to process
        :type in_lake_tile_label: int
        
        :return: tile reference
        :rtype: string
        """
        
        # 1 - Get unique values and counts of tiles for pixels with edge_label in_lake_tile_label
        unique, counts = np.unique(self.tile_index[in_pixc_index][np.where(self.labels[in_pixc_index] == in_lake_tile_label)], return_counts = True)

        # 2 - Get tile ref corresponding to the max number of pixels
        out_tile_max_pix = str(unique[np.where(counts == max(counts))][0]).rjust(3, str('0')) + self.swath_side

        return out_tile_max_pix
        
    # ----------------------------------------

    def get_lake_tile_label(self, in_new_label):
        """
        This function retrieves all old edge_label for pixels having labels=in_new_label.
        The given new label corresponds to a global continent-pass label, corresponding to several old edge_labels.
        The old label involving the largest number of pixels is returned.
            
        :param in_new_label: global new label
        :type in_new_label: int
        
        :return: LakeTile_Edge label involving the largest number of pixels
        :rtype: string
        """
        # 1 - Retrieve LakeTile labels corresponding to current in_new_label label
        label_idx = np.where(self.labels == in_new_label)

        # 2 - Get unique values and counts
        unique, counts = np.unique(self.laketile_labels[label_idx], return_counts=True)

        # 3 - Find major value
        out_final_label = unique[np.argsort(-counts)][0]

        return str(out_final_label[5:])

    # ----------------------------------------

    def get_near_range(self, in_indices):
        """
        This function is designed to return a near range array, with the near range of each tile corresponding to each pixels of in_indices.

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

