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
# VERSION:5.2.1:DM:#91:2023/07/24:Poursuite industrialisation
# VERSION:5.3.0:DM:#91:2023/10/04:Poursuite industrialisation
# VERSION:6.1.0:DM:#91:2024/09/20:Poursuite industrialisation
# VERSION:6.2.0:DM:#91:2024/11/29:Poursuite industrialisation
# VERSION:6.2.2:DM:#822:2025/06/30:Nouveau patch Lake_TILE_V6.2.1
# VERSION:6.2.3:FA:#839:2025/11/10:[LAKE_TILE] L2_HR_LAKETILE fails on C032P429T305R with message ValueError: list.remove(x): x not in list
# VERSION:6.2.6:DM:#900:2026/07/06: [PGE] [LAKES] Ajout d'une interface niveau moyen des lacs
# VERSION:6.2.6:DM:#900:2026/07/06:SSALTO-565 [PGE] [LAKES] Ajout d'une interface coeff NN pour calcul ice_dyn_f
# VERSION:6.2.7:DM:#906:2026/07/28:SSALTO-872 [PGE][LAKES] Gestion des erreurs/warning GDAL
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: proc_pixc.py
   :synopsis: Deals with SWOT pixel cloud product
    Created on 2017/02/28

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import numpy as np
import os
from osgeo import ogr, osr
from scipy import interpolate
from scipy.ndimage.morphology import binary_erosion
from dateutil import parser

import cnes.common.service_config_file as service_config_file
import cnes.common.service_error as service_error

import cnes.common.lib.my_hull as my_hull
import cnes.common.lib.my_timer as my_timer
import cnes.common.lib.my_netcdf_file as my_nc
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_segmentation as my_segmentation
import cnes.common.lib.my_variables as my_var
import cnes.common.lib_lake.locnes_products_netcdf as nc_file
import cnes.modules.geoloc.lib.geoloc as my_geoloc

import jpl.modules.aggregate_cnes as jpl_aggregate

import udes.pixel_ice as pixel_ice


class PixelCloud(object):
    """
    class PixelCloud
    Deal with SWOT pixel cloud product
    """
    
    def __init__(self):
        """
        Constructor: init pixel cloud object
        
        Variables of the object:

        - From pixel_cloud group in L2_HR_PIXC file
        
            Having the same number of PIXC as in the L2_HR_PIXC file
            - origin_classif / 1D-array of byte: classification value of pixels (= variable named classification in L2_HR_PIXC)
            - origin_range_index / 1D-array of int: range indices of water pixels (= variable named range_index in L2_HR_PIXC)
            - origin_azimuth_index / 1D-array of int: azimuth indices of water pixels (= variable named azimuth_index in L2_HR_PIXC)
            - origin_latitude / 1D-array of float: latitude of water pixels (= variable named longitude in L2_HR_PIXC)
            - origin_longitude / 1D-array of float: longitude of water pixels (= variable named latitude in L2_HR_PIXC)
            - origin_cross_track / 1D-array of float: cross-track distance of water pixels (= variable named cross_track in L2_HR_PIXC)
            - origin_cross_track_refined / 1D-array of float: refined cross-track distance from nadir to center of water pixels
            
            Subset of L2_HR_PIXC (after rejection of some PIXC, as river pixels except those associated to reservoirs)
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
            - false_detection_rate / 1D array of float: false detection rate
            - missed_detection_rate / 1D array of float: missed detection rate
            - bright_land_flag / 1D array of byte: bright land flag
            - layover_impact /1D array of float: layover impact
            - eff_num_rare_looks / 1D array of byte: number of rare looks
            - latitude / 1D-array of float: latitude of water pixels
            - longitude / 1D-array of float: longitude of water pixels
            - height / 1D-array of float: height of water pixels
            - cross_track / 1D-array of float: cross-track distance from nadir to center of water pixels
            - pixel_area / 1D-array of int: area of water pixels
            - inc / 1D-array of float: incidence angle
            - phase_noise_std / 1D-array of float: phase noise standard deviation
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
            - interferogram_qual / 1D array of uint: flag that indicates the quality of the interferogram quantities 
            - classification_qual / 1D array of uint: flag that indicates the quality of the classification quantities 
            - geolocation_qual / 1D array of uint: flag that indicates the quality of the geolocation quantities 
            - sig0_qual / 1D array of uint: quality flag for sig0
            - prev_xover / 1D array of float: time difference from the previous crossover (= variable named time_from_prev_xover in L2_HR_PIXC file)
            - next_xover / 1D array of float: time difference to the next crossover (= variable named time_to_next_xover in L2_HR_PIXC file)
            - wavelength / float: wavelength corresponding to the effective radar carrier frequency 
            - looks_to_efflooks / float: ratio between the number of actual samples and the effective number of independent samples during spatial averaging over a large 2-D area
        - From tvp group in L2_HR_PIXC file
            - nadir_time[_tai] / 1D-array of float: observation UTC [TAI] time of each nadir pixel (= variable named time[tai] in L2_HR_PIXC file)
            - nadir_time_all / 1D-array of float: observation UTC time at nadir corresponding to all PIXC (= variable named time in L2_HR_PIXC file)
            - nadir_longitude / 1D-array of float: longitude of each nadir pixel (= variable named longitude in L2_HR_PIXC file)
            - nadir_latitude / 1D-array of float: latitude of each nadir pixel (= variable named latitude in L2_HR_PIXC file)   
            - nadir_[x|y|z] / 1D-array of float: [x|y|z] cartesian coordinates of each nadir pixel (= variables named [x|y|z] in L2_HR_PIXC file)
            - nadir_[vx|vy|vz] / 1D-array of float: velocity vector of each nadir pixel in cartesian coordinates (= variables 
              named velocity_unit_[x|y|z] in L2_HR_PIXC file)
            - nadir_plus_y_antenna_[x|y|z] / 1D-array of float: position vector of the +y KaRIn antenna phase center in ECEF coordinates (= variables named plus_y_antenna_[x|y|z] in L2_HR_PIXC file)
            - nadir_minus_y_antenna_[x|y|z] / 1D-array of float: position vector of the -y KaRIn antenna phase center in ECEF coordinates (= variables named minus_y_antenna_[x|y|z] in L2_HR_PIXC file)
        - From processing
            - product_type / string: type of product among "SP"=LakeSP and "TILE"=LakeTile
            - tile_poly / ogr.Polygon: polygon of the PixC tile
            - az_0_line / ogr.LineString: line with azimuth = 0
            - az_max_line / ogr.LineString: line with azimuth = az_max
            - date / datetime.datetime : date of acquisition
            - continent / string: continent covered by the tile (if global var CONTINENT_FILE exists)
            - cross_track_refined / 1D-array of float: refined cross-track distance from nadir to center of water pixels
            - selected_index / 1D-array of int: indices from original 1D-arrays of not rejected pixels with specified classification indices
            - nb_selected / int: number of selected pixels (=selected_index.size)
            - nb_water_pix / int: number of water pixels
            - list_classif_keys / list: list of keys of classif_dict
            - classif_dict / dict: selection (by True or False) of pixels wrt their classification, for different processes + already init arrays useful for aggregate.py lib; key=water/open_water/interior_water/water_edge/dark/low_coherence/mandatory/specular_ringing/not_specular_ringing/4hull/4wse/4area/Idw_in/Idw_in_no_dw/Idw/Idw_no_dw/Ide/I/I_no_dw
            - interferogram_flattened / 1D-array of complex: flattened interferogram
            - inundated_area / 1D-array of int: area of pixels covered by water
            - pixc_first_line_index / int: azimuth index corresponding to the first line in tile
            - pixc_last_line_index / int: azimuth index corresponding to the last line in tile
            - seg_hr_start / bool: =True if tile is a start for an HR segment
            - seg_hr_stop / bool: =False if tile is a stop for an HR segment
            - origin_pixel_at_window_edge / 1D-array of int: =1 for pixels having their range_index == pixel_cloud.data_window_first|last_valid; =0 otherwise (same number of PIXC as in the L2_HR_PIXC file)
            - pixel_at_window_edge / 1D-array of boolean: =1 for pixels having their range_index == pixel_cloud.data_window_firts|last_valid; =0 otherwise 
            - is_boundary_pix / 1D-array of boolean: =1 for pixels belonging to a start or a stop of an HR segment
            - height_std_pix / 1D-array of float: height std
            - corrected_height / 1D-array of float: height corrected from geoid and other tide corrections
            - ice_dyn_f / 1D-array of int: ice dynamical flag computed from SWOT observations for each pixel
            - labels / 1D-array of int: labelled regions associated to each PixC water pixel; pixels of this vector correspond one-to-one to the pixels of data from L2_HR_PIXC and L2_HR_PIXCVecRiver
            - labels_unique / 1D-array of int: list of unique labels in labels
            - labels_counts / 1D-array of int: list of number of pixels associated to labels_unique
            - labels_coastal_zones / set of int: list of water regions associated to coastal zones
            - labels_big_lakes / set of int: list of water regions associated to big lakes
            - big_lakeid / set of string: list of PLD lake_id considered as big for the tile
            - labels_coastal_zones / / 1D-array of int: list of water regions associated to coastal zones
            - nb_obj / int : number of separate regions in the PixC tile
            - labels_inside / 1D-array of int: label of objects entirely inside the tile
            - nb_obj_inside / int : number of these objects
            - labels_[at_top|at_bottom|at_both]_edge / 1D-array of int: label of objects at the top/bottom/both edges of the tile
            - nb_obj_[at_top|at_bottom|at_both]_edge / int : number of these objects 
            - edge_index / 1D-array of int: indices of pixels contained in objects at top/bottom edges
            - edge_label / 1D-array of int: object label for each pixel contained in objects at top/bottom edges
            - edge_loc / 1D-array of int: object edge location (0=bottom 1=top 2=both) for each pixel contained in objects at top/bottom edges
            - nb_edge_pix / int: number of pixels contained in objects at top/bottom edges
        """
        # Get instance of service config file
        self.cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        self.product_type = "TILE"

        # Init PIXC variables
        self.origin_classif = None  
        self.origin_range_index = None
        self.origin_pixel_at_window_edge = None
        self.origin_azimuth_index = None
        self.origin_longitude = None
        self.origin_latitude = None
        self.origin_cross_track = None
        self.origin_cross_track_refined = None
        self.classif = None  
        self.range_index = None
        self.pixel_at_window_edge = None
        self.azimuth_index = None
        self.interferogram = None
        self.power_plus_y = None
        self.power_minus_y = None
        self.water_frac = None
        self.water_frac_uncert = None
        self.false_detection_rate = None
        self.missed_detection_rate = None
        self.bright_land_flag = None
        self.layover_impact = None
        self.eff_num_rare_looks = None
        self.latitude = None
        self.longitude = None
        self.height = None
        self.cross_track = None
        self.cross_track_refined = None
        self.pixel_area = None
        self.inc = None
        self.phase_noise_std = None
        self.dlatitude_dphase = None
        self.dlongitude_dphase = None
        self.dheight_dphase = None
        self.darea_dheight = None
        self.eff_num_medium_looks = None
        self.sig0 = None
        self.sig0_u = None
        self.model_dry_tropo_cor = None
        self.model_wet_tropo_cor = None
        self.iono_cor_gim_ka = None
        self.height_cor_xover = None
        self.geoid = None
        self.solid_earth_tide = None
        self.load_tide_fes = None
        self.load_tide_got = None
        self.pole_tide = None
        self.interferogram_qual = None
        self.classification_qual = None
        self.geolocation_qual = None
        self.sig0_qual = None
        self.prev_xover = None
        self.next_xover = None
        self.wavelength = -9999.0
        self.looks_to_efflooks = -9999.0
        self.nadir_time = None
        self.nadir_time_all = None
        self.nadir_time_tai = None
        self.nadir_longitude = None
        self.nadir_latitude = None
        self.nadir_x = None
        self.nadir_y = None
        self.nadir_z = None
        self.nadir_vx = None
        self.nadir_vy = None
        self.nadir_vz = None
        self.nadir_plus_y_antenna_x = None
        self.nadir_plus_y_antenna_y = None
        self.nadir_plus_y_antenna_z = None
        self.nadir_minus_y_antenna_x = None
        self.nadir_minus_y_antenna_y = None
        self.nadir_minus_y_antenna_z = None
        
        # Init dictionary of PIXC metadata
        self.pixc_metadata = {}
        self.pixc_metadata["source"] = ""
        self.pixc_metadata["cycle_number"] = -9999
        self.pixc_metadata["pass_number"] = -9999
        self.pixc_metadata["tile_number"] = -9999
        self.pixc_metadata["swath_side"] = ""
        self.pixc_metadata["tile_name"] = ""
        self.pixc_metadata["time_granule_start"] = ""
        self.pixc_metadata["time_granule_end"] = ""
        self.pixc_metadata["time_coverage_start"] = ""
        self.pixc_metadata["time_coverage_end"] = ""
        self.pixc_metadata["geospatial_lon_min"] = -9999
        self.pixc_metadata["geospatial_lon_max"] = -9999
        self.pixc_metadata["geospatial_lat_min"] = -9999
        self.pixc_metadata["geospatial_lat_max"] = -9999
        self.pixc_metadata["inner_first_latitude"] = -9999.0
        self.pixc_metadata["inner_first_longitude"] = -9999.0
        self.pixc_metadata["inner_last_latitude"] = -9999.0
        self.pixc_metadata["inner_last_longitude"] = -9999.0
        self.pixc_metadata["outer_first_latitude"] = -9999.0
        self.pixc_metadata["outer_first_longitude"] = -9999.0
        self.pixc_metadata["outer_last_latitude"] = -9999.0
        self.pixc_metadata["outer_last_longitude"] = -9999.0
        self.pixc_metadata["continent_id"] = "None"
        self.pixc_metadata["continent_code"] = "None"
        self.pixc_metadata["wavelength"] = -9999.0
        self.pixc_metadata["near_range"] = -9999.0
        self.pixc_metadata["nominal_slant_range_spacing"] = -9999.0
        self.pixc_metadata["polarization"] = ""
        self.pixc_metadata["interferogram_size_range"] = -9999
        self.pixc_metadata["interferogram_size_azimuth"] = -9999
        self.pixc_metadata["looks_to_efflooks"] = -9999.0
        self.pixc_metadata["ellipsoid_semi_major_axis"] = ""
        self.pixc_metadata["ellipsoid_flattening"] = ""

        # Variables specific to processing
        self.nb_pix_range = 0  # Number of pixels in range dimension
        self.nb_pix_azimuth = 0  # Number of pixels in azimuth dimension
        self.tile_poly = None  # Polygon representing the PIXC tile
        self.az_0_line = None  # Line with azimuth = 0
        self.az_max_line  = None  # Line with azimuth = az_max
        self.date = None  # Date of acquisition
        self.continent_id = None  # 2-letter identifier of the continent covered by the tile
        self.continent_code = None  # 1-digit code of the continent covered by the tile
        self.selected_index = None  # Indices of selected pixels
        self.nb_selected = 0  # Number of selected pixels
        self.nb_water_pix = 0  # Number of water pixels
        self.list_classif_keys = ["water", "open_water", "water_edge", "dark", "low_coherence", 
                                  "mandatory", "specular_ringing", "not_specular_ringing",
                                  "4hull", "4wse", "4area", 
                                  "Idw_in", "Idw_in_no_dw", "Idw", "Idw_no_dw", "Ide", "I", "I_no_dw"]
        self.classif_dict = None  # Selection of pixels per process, wrt their classification
        self.interferogram_flattened = None  # Flattened interferogram
        self.inundated_area = None  # Area of pixel where water
        self.pixc_first_line_index = None  # Azimuth index corresponding to the first line in tile
        self.pixc_last_line_index = None   # Azimuth index corresponding to the last line in tile
        self.seg_hr_start = False  # Tile being a start for an HR segment (ie pixels aren't considered at bottom_edge)
        self.seg_hr_stop = False  # Tile being a stop for an HR segment (ie pixels aren't considered at top_edge)
        self.is_boundary_pix = []  # =1 for pixels belonging to a start or a stop of an HR segment
        self.height_std_pix = None  # Height std
        self.corrected_height = None  # Height corrected from geoid and other tide corrections
        self.ice_dyn_f = None  # Ice dynamical flag computed from SWOT observations for each pixel
        self.labels = None  # Vector of entity labels associated to each pixel
        self.labels_unique = None  # Vector of unique labels associated in self.labels
        self.labels_counts = None  # Vector of number of pixels associated to labels_unique
        self.labels_coastal_zones = set()  # Vector of water regions associated to coastal zones
        self.labels_big_lakes = set()  # Vector of water regions associated to big lakes
        self.big_lakeid = set()  # Vector of PLD lake_id considered as big for the tile
        self.nb_obj = None  # Number of separate water regions
        self.labels_inside = None  # Labels of water regions entirely inside the tile
        self.nb_obj_inside = None  # Number of water regions inside the tile
        self.labels_at_top_edge = set()  # Labels of water regions at the top edge of the tile
        self.nb_obj_at_top_edge = None  # Number of water regions at the top edge of the tile
        self.labels_at_bottom_edge = set()  # Labels of water regions at the bottom edge of the tile
        self.nb_obj_at_bottom_edge = None  # Number of water regions at the bottom edge of the tile
        self.labels_at_both_edges = set()  # Labels of water regions at the top and bottom edges of the tile
        self.nb_obj_at_both_edges = None  # Number of water regions at the top and bottom edges of the tile
        self.edge_index = np.array([], dtype=int)  # Indices of pixels contained in objects at top/bottom edges
        self.edge_label = np.array([], dtype=int)  # Object label for each pixel contained in objects at top/bottom edges
        self.edge_loc = np.array([], dtype=int)  # Object edge location (0=bottom 1=top 2=both) for each pixel contained in objects at top/bottom edges
        self.nb_edge_pix = 0  # Number of pixels contained in objects at top/bottom edges

    def get_params_from_config_file(self):
        """
        Get parameters from configuration file
        
        :return: classif_list = CLASSIF_LIST parameter from configuration file
        :rtype: classif_list = list
        :return: use_fractional_inundation = USE_FRACTIONAL_INUNDATION parameter from configuration file
        :rtype: use_fractional_inundation = list
        :return: pixc_qual_list = PIXC_QUAL_LIST parameter from configuration file
        :rtype: pixc_qual_list = list
        :return: pixc_qual_mask = PIXC_QUAL_MASK parameter from configuration file
        :rtype: pixc_qual_mask = list
        :return: exclude_bright_land = EXCLUDE_BRIGHT_LAND parameter from configuration file
        :rtype: exclude_bright_land = list
        """

        # List of classif flags to keep for processing
        classif_list = [int(i) for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_LIST").split(";")]
        
        # True/False list of classif flags to use water fraction for inundated computation
        use_fractional_inundation = [i == "True" for i in self.cfg.get("CONFIG_PARAMS", "USE_FRACTIONAL_INUNDATION").split(";")]
        
        # List of quality flags to consider for PixC selection
        if self.cfg.get("CONFIG_PARAMS", "PIXC_QUAL_LIST") == "":
            pixc_qual_list = None
        else:
            pixc_qual_list = self.cfg.get("CONFIG_PARAMS", "PIXC_QUAL_LIST").split(";")
            
        # Mask of quality flags to use for PixC selection
        if self.cfg.get("CONFIG_PARAMS", "PIXC_QUAL_MASK") == "":
            pixc_qual_mask = None
        else:
            pixc_qual_mask = self.cfg.get("CONFIG_PARAMS", "PIXC_QUAL_MASK").split(";")
        if (pixc_qual_list is not None) and (len(pixc_qual_list) != len(pixc_qual_mask)):
            message = "PIXC_QUAL_LIST and PIXC_QUAL_MASK must have the same number of elements"
            raise service_error.ConfigError(message)
            
        # If EXCLUDE_BRIGHT_LAND in lake_tile_param.cfg: retrieved as a string => getboolean needed to convert as a boolean
        # If not: init as a boolean => must use get and not getboolean
        try:
            exclude_bright_land = self.cfg.getboolean("CONFIG_PARAMS", "EXCLUDE_BRIGHT_LAND")
        except:
            exclude_bright_land = self.cfg.get("CONFIG_PARAMS", "EXCLUDE_BRIGHT_LAND")

        return classif_list, use_fractional_inundation, pixc_qual_list, pixc_qual_mask, exclude_bright_land
        
    def extract_from_pixc(self, in_pixc_file, in_index_reject, flag_simu=False, ice_model=None):
        """
        Retrieve needed data from pixel cloud
        
        :param in_pixc_file: full path of L2_HR_PIXC file
        :type in_pixc_file: string
        :param in_index_reject: list of indices to reject before all processing
        :type in_index_reject: 1D-array of int
        :param flag_simu: =True if input PIXC and PIXCVecRiver files are simulated data; =False otherwise (default)
        :type flag_simu: boolean
        :param ice_model: pre-trained classification model used to compute dynamical ice flags
        :type ice_model: 
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("L2_HR_PIXC file = %s" % in_pixc_file)
        
        # Get flag variables from configuration file
        classif_list, use_fractional_inundation, pixc_qual_list, pixc_qual_mask, exclude_bright_land = self.get_params_from_config_file()
        logger.debug("Keep pixels with classification flags = %s" % classif_list)
        
        # 1 - Open pixel cloud file in reading mode
        pixc_reader = my_nc.MyNcReader(in_pixc_file)
        pixc_group = pixc_reader.content.groups['pixel_cloud']
        sensor_group = pixc_reader.content.groups['tvp']
        
        # 2 - Retrieve needed global attributes
        pixc_keys = pixc_reader.get_list_att()
        for key in self.pixc_metadata.keys():
            if key in pixc_keys:
                self.pixc_metadata[key] = pixc_reader.get_att_value(key)
        # pixel_cloud attributes
        pixc_keys = pixc_reader.get_list_att(in_group=pixc_group)
        for key in self.pixc_metadata.keys():
            if key in pixc_keys:
                self.pixc_metadata[key] = pixc_reader.get_att_value(key, in_group=pixc_group)
                    
        # Save specific variables
        self.nb_pix_range = int(self.pixc_metadata["interferogram_size_range"])  # Number of pixels in range dimension
        self.nb_pix_azimuth = int(self.pixc_metadata["interferogram_size_azimuth"])  # Number of pixels in azimuth dimension
        self.near_range = np.double(self.pixc_metadata["near_range"])  # Slant range for the 1st image pixel
        self.wavelength = np.float64(self.pixc_metadata["wavelength"])  # Wavelength
        self.polarization = str(self.pixc_metadata["polarization"])  # Radar signal polarization of the tile
        # Ratio between the number of actual samples and the effective number of independent samples
        self.looks_to_efflooks = np.float64(self.pixc_metadata["looks_to_efflooks"])
        self.date = parser.parse(self.pixc_metadata["time_granule_start"])  # Date of acquisition

        # Azimuth indices, corresponding to begin and end of the tile (without SLC overlapping)
        # Try/except because pixc_line_qual not in PIXC generated by Large Scale Simulator
        
        try:
            # Case when there is pixel_cloud.pixc_line_qual variable
            pixc_line_to_tvp = pixc_reader.get_var_value("pixc_line_qual", in_group=pixc_group)
            pixc_line_to_tvp_mod2 = pixc_line_to_tvp % 2  # =1 where not_in_tile; =0 otherwise
            if pixc_line_to_tvp_mod2[0] == 0:
                self.pixc_first_line_index = 0
                if not flag_simu:
                    self.seg_hr_start = True
            if pixc_line_to_tvp_mod2[-1] == 0:
                self.pixc_last_line_index = len(pixc_line_to_tvp) - 1
                if not flag_simu:
                    self.seg_hr_stop = True
            pixc_line_to_tvp_roll1 = np.roll(pixc_line_to_tvp_mod2, 1)  # Roll of 1 index to right
            pixc_line_to_tvp_var = (pixc_line_to_tvp_mod2 + pixc_line_to_tvp_roll1) % 2  # =1 for tile edges
            index_start_stop_tile = np.where(pixc_line_to_tvp_var)[0]
            if self.pixc_first_line_index is None:
                self.pixc_first_line_index = index_start_stop_tile[0]
            if self.pixc_last_line_index is None:
                self.pixc_last_line_index = index_start_stop_tile[1] - 1
            logger.debug("Use of pixc_line_qual variable to define the tile edge")
            
        except:
            # Case when there is no pixel_cloud.pixc_line_qual variable
            logger.debug("There is no pixel_cloud.pixc_line_qual variable in the PIXC \
                         => use interferogram_size_azimuth attribute to define the tile edge")
            self.pixc_first_line_index = 0
            self.pixc_last_line_index = self.nb_pix_azimuth - 1
        logger.debug("> pixc_first_line_index = %d - pixc_last_line_index = %d" % (self.pixc_first_line_index, self.pixc_last_line_index))
        if self.seg_hr_start:
            logger.debug("WARNING - Start of HR segment => no bottom_edge pixel")
        if self.seg_hr_stop:
            logger.debug("WARNING - Stop of HR segment => no top_edge pixel")

        # 3 - Create polygon of tile from global attributes
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(my_tools.convert_to_m180_180(float(pixc_reader.get_att_value("inner_first_longitude"))),
                                                   float(pixc_reader.get_att_value("inner_first_latitude")))
        ring.AddPoint(my_tools.convert_to_m180_180(float(pixc_reader.get_att_value("outer_first_longitude"))),
                                                   float(pixc_reader.get_att_value("outer_first_latitude")))
        ring.AddPoint(my_tools.convert_to_m180_180(float(pixc_reader.get_att_value("outer_last_longitude"))), 
                                                   float(pixc_reader.get_att_value("outer_last_latitude")))
        ring.AddPoint(my_tools.convert_to_m180_180(float(pixc_reader.get_att_value("inner_last_longitude"))),
                                                   float(pixc_reader.get_att_value("inner_last_latitude")))
        ring.AddPoint(my_tools.convert_to_m180_180(float(pixc_reader.get_att_value("inner_first_longitude"))),
                                                   float(pixc_reader.get_att_value("inner_first_latitude")))
        self.tile_poly = ogr.Geometry(ogr.wkbPolygon)
        self.tile_poly.AddGeometry(ring)

        if not self.tile_poly.IsValid():
            logging.warning("Tile polygon is invalid : geometry is fixed using ConvexHull")
            lon, lat, z = zip(*ring.GetPoints())
            self.tile_poly = my_hull.get_convex_hull(np.array(lon), np.array(lat))

        # Build line with azimuth = 0
        self.az_0_line = ogr.Geometry(ogr.wkbLineString)
        self.az_0_line.AddPoint(my_tools.convert_to_m180_180(float(pixc_reader.get_att_value("inner_first_longitude"))),
                                                   float(pixc_reader.get_att_value("inner_first_latitude")))
        self.az_0_line.AddPoint(my_tools.convert_to_m180_180(float(pixc_reader.get_att_value("outer_first_longitude"))),
                                                   float(pixc_reader.get_att_value("outer_first_latitude")))

        # Build line with azimuth max
        self.az_max_line = ogr.Geometry(ogr.wkbLineString)
        self.az_max_line.AddPoint(my_tools.convert_to_m180_180(float(pixc_reader.get_att_value("inner_last_longitude"))),
                                                   float(pixc_reader.get_att_value("inner_last_latitude")))
        self.az_max_line.AddPoint(my_tools.convert_to_m180_180(float(pixc_reader.get_att_value("outer_last_longitude"))),
                                                   float(pixc_reader.get_att_value("outer_last_latitude")))

        # 4 - Retrieve variables to help the PixC selection
        # 4.1 - Classification flag
        self.origin_classif = pixc_reader.get_var_value("classification", in_group=pixc_group)
        nb_pixels_origin = int(pixc_reader.get_dim_value("points", in_group=pixc_group))
        # 4.2 - Range indices of water pixels
        self.origin_range_index = pixc_reader.get_var_value("range_index", in_group=pixc_group)
        # 4.3 - Azimuth indices of water pixels
        self.origin_azimuth_index = pixc_reader.get_var_value("azimuth_index", in_group=pixc_group)
        # 4.4 - Longitude
        origin_longitude = my_tools.convert_to_m180_180(pixc_reader.get_var_value("longitude", in_group=pixc_group))
        # 4.5 - Latitude
        origin_latitude = pixc_reader.get_var_value("latitude", in_group=pixc_group)
        # 4.6 - Height
        origin_height = pixc_reader.get_var_value("height", in_group=pixc_group)
        # 4.7 - Approximate cross-track  location
        self.origin_cross_track = pixc_reader.get_var_value("cross_track", in_group=pixc_group)
        # 4.8 - Illumination time
        origin_illumination_time = pixc_reader.get_var_value("illumination_time", in_group=pixc_group)
        # 4.9 - PixC quality flags to consider (geolocation_qual and classification_qual obligatory)
        # Build list of quality flags to retrieve
        if pixc_qual_list is None:
            tmp_list = list()
        else:
            tmp_list = pixc_qual_list.copy()
        for tmp_flag in ["classification_qual", "geolocation_qual"]:
            if tmp_flag not in tmp_list:
                tmp_list.append(tmp_flag)
        # Retrieve quality flags
        origin_qual = dict()
        for flag_name in tmp_list:
            origin_qual[flag_name] = pixc_reader.get_var_value(flag_name, in_group=pixc_group)
        # 4.10 - bright_land_flag
        origin_bright_land_flag = pixc_reader.get_var_value("bright_land_flag", in_group=pixc_group)
        # 4.11 - Data window first and last range indices
        try:
            origin_data_window_first_valid = pixc_reader.get_var_value("data_window_first_valid", in_group=pixc_group)
        except:  # For Large Scale Simulator PIXC
            origin_data_window_first_valid = np.zeros(self.nb_pix_azimuth, dtype=int)
        try:
            origin_data_window_last_valid = pixc_reader.get_var_value("data_window_last_valid", in_group=pixc_group)
        except:  # For Large Scale Simulator PIXC
            origin_data_window_last_valid = np.zeros(self.nb_pix_azimuth, dtype=int) + self.nb_pix_range-1
        self.origin_pixel_at_window_edge = np.zeros(nb_pixels_origin, dtype=int)
        
        # 5 - Compute true cross track  
        
        # 5.1 - Interpolate nadir_time wrt illumination time
        origin_nadir_time = pixc_reader.get_var_value("time", in_group=sensor_group)  # Read nadir_time values

        # 5.2 - Link TVP and pixel_cloud indices
        # Try/except because pixc_line_to_tvp not in PIXC generated by Large Scale Simulator
        try:
            # Case when there is pixel_cloud.pixc_line_to_tvp variable
            pixc_line_to_tvp = pixc_reader.get_var_value("pixc_line_to_tvp", in_group=pixc_group)
            nadir_index_all = pixc_line_to_tvp[self.origin_azimuth_index].astype(int)
            logger.debug("Use of pixel_cloud.pixc_line_to_tvp variable to link TVP and pixel_cloud indices")

        except:
            # Case when there is no pixel_cloud.pixc_line_to_tvp variable
            logger.warning("There is no pixc_line_to_tvp variable in the PIXC \
                           => use illumination_time variable to link TVP and pixel_cloud indices")
            if len(origin_nadir_time) > 1:
                f = interpolate.interp1d(origin_nadir_time, range(len(origin_nadir_time)))  # Interpolator
                nadir_index_all = (np.rint(f(origin_illumination_time))).astype(int)  # Link between all PixC and nadir pixels
            else:
                nadir_index_all = np.zeros(len(origin_illumination_time), dtype=int)

        # Nadir cartesian coordinates
        nadir_x = pixc_reader.get_var_value("x", in_group=sensor_group)[nadir_index_all]
        nadir_y = pixc_reader.get_var_value("y", in_group=sensor_group)[nadir_index_all]
        nadir_z = pixc_reader.get_var_value("z", in_group=sensor_group)[nadir_index_all]
        # Nadir velocity in cartesian coordinates
        nadir_vx = pixc_reader.get_var_value("vx", in_group=sensor_group)[nadir_index_all]
        nadir_vy = pixc_reader.get_var_value("vy", in_group=sensor_group)[nadir_index_all]
        nadir_vz = pixc_reader.get_var_value("vz", in_group=sensor_group)[nadir_index_all]
                                     
        # Convert geodetic coordinates (lat, lon, height) to cartesian coordinates (x, y, z)
        x, y, z = my_geoloc.convert_llh2ecef(origin_latitude, origin_longitude,
                                             origin_height, my_var.GEN_RAD_EARTH_EQ,
                                             my_var.GEN_RAD_EARTH_POLE)

        # 4.7.2 - Cross track estimated from SWOT PIXC positions
        self.origin_cross_track_refined = my_geoloc.compute_true_cross_track_wrt_swot_convention(np.transpose(np.array([x, y, z])),
                                                                                                 np.transpose(np.array([nadir_x, nadir_y, nadir_z])),
                                                                                                 np.transpose(np.array([nadir_vx, nadir_vy, nadir_vz])))
    
        # 5 - Select PixC given criteria (Non-river, Finite values, quality flags)
        logger.debug("There are %d pixels in the input PIXC" % nb_pixels_origin)
        self.select_pixc_from_criteria(in_index_reject,
                                       origin_longitude, origin_latitude, origin_height, origin_illumination_time,
                                       origin_data_window_first_valid, origin_data_window_last_valid,
                                       origin_bright_land_flag, origin_qual, 
                                       classif_list, pixc_qual_list, pixc_qual_mask, exclude_bright_land)

        # 6 - Retrieve PixC data only for selected pixels
        self.nadir_time_all = origin_nadir_time[nadir_index_all]
        
        if self.nb_selected != 0:
            
            # 6.1 - In PixC group
            
            # Classification flags
            self.classif = self.origin_classif[self.selected_index]
            
            # Init dictionnary of classification wrt to selected classification flags set as parameters
            if "classification_qual" in origin_qual.keys():
                self.classif_dict = self.init_classif_dict(classif_list, self.nb_selected, self.classif, 
                                                           classif_qual=origin_qual["classification_qual"][self.selected_index])
            else:
                self.classif_dict = self.init_classif_dict(classif_list, self.nb_selected, self.classif)
            
            # Range indices of water pixels
            self.range_index = self.origin_range_index[self.selected_index]
            # Number of water pixels
            self.nb_water_pix = self.range_index.size
            # Range indices of water pixels
            self.azimuth_index = self.origin_azimuth_index[self.selected_index]
            # Localisation of pixels wrt window edge
            self.pixel_at_window_edge = self.origin_pixel_at_window_edge[self.selected_index]
            # Find pixels belonging to a start or a stop of the HR mask
            self.is_boundary_pix = np.zeros((self.nb_selected))
            if self.seg_hr_start:
                self.is_boundary_pix[self.azimuth_index == self.pixc_first_line_index] = 1
            if self.seg_hr_stop:
                self.is_boundary_pix[self.azimuth_index == self.pixc_last_line_index] = 1
            
            # Interferogram
            interferogram = pixc_reader.get_var_value("interferogram", in_group=pixc_group)[self.selected_index]
            self.interferogram = interferogram[:,0] + 1j*interferogram[:,1]
            self.interferogram_flattened = 0 * self.interferogram 
            # Sensitivity of height estimate to interferogram phase
            self.power_plus_y = pixc_reader.get_var_value("power_plus_y", in_group=pixc_group)[self.selected_index]            
            # Sensitivity of height estimate to interferogram phase
            self.power_minus_y = pixc_reader.get_var_value("power_minus_y", in_group=pixc_group)[self.selected_index]   
            
            # Water fraction
            self.water_frac = pixc_reader.get_var_value("water_frac", in_group=pixc_group)[self.selected_index]
            # Water fraction uncertainty
            self.water_frac_uncert = pixc_reader.get_var_value("water_frac_uncert", in_group=pixc_group)[self.selected_index]
            # False detection rate
            self.false_detection_rate = pixc_reader.get_var_value("false_detection_rate", in_group=pixc_group)[self.selected_index]
            # Missed detection rate
            self.missed_detection_rate = pixc_reader.get_var_value("missed_detection_rate", in_group=pixc_group)[self.selected_index]
            # Bright land flag
            self.bright_land_flag = pixc_reader.get_var_value("bright_land_flag", in_group=pixc_group)[self.selected_index]
            # Layover impact
            self.layover_impact = pixc_reader.get_var_value("layover_impact", in_group=pixc_group)[self.selected_index]
            # Number of rare looks
            self.eff_num_rare_looks = pixc_reader.get_var_value("eff_num_rare_looks", in_group=pixc_group)[self.selected_index]
            
            # Latitude
            self.latitude = origin_latitude[self.selected_index]
            # Longitude
            self.longitude = origin_longitude[self.selected_index]
            # Height
            self.height = origin_height[self.selected_index]
            
            # Cross-track distance
            self.cross_track = self.origin_cross_track[self.selected_index]
            self.cross_track_refined = self.origin_cross_track_refined[self.selected_index]
            # Pixel area
            self.pixel_area = pixc_reader.get_var_value("pixel_area", in_group=pixc_group)[self.selected_index]
            # Inundated area
            self.inundated_area = np.copy(self.pixel_area)
            for classif_flag, frac_ok in zip(classif_list, use_fractional_inundation):
                if frac_ok:
                    ind_ok = np.where((self.classif == classif_flag))
                    if len(ind_ok) > 0:
                        logger.debug("=> Use water fraction to compute pixel area for classif flag %d" % classif_flag)
                        self.inundated_area[ind_ok] = self.pixel_area[ind_ok] * self.water_frac[ind_ok]
            # Incidence angle
            self.inc = pixc_reader.get_var_value("inc", in_group=pixc_group)[self.selected_index]
            # Phase noise standard deviation
            self.phase_noise_std = pixc_reader.get_var_value("phase_noise_std", in_group=pixc_group)[self.selected_index]
            
            # Sensitivity of latitude estimate to interferogram phase
            self.dlatitude_dphase = pixc_reader.get_var_value("dlatitude_dphase", in_group=pixc_group)[self.selected_index]
            # Sensitivity of longitude estimate to interferogram phase
            self.dlongitude_dphase = pixc_reader.get_var_value("dlongitude_dphase", in_group=pixc_group)[self.selected_index]
            # Sensitivity of height estimate to interferogram phase
            self.dheight_dphase = pixc_reader.get_var_value("dheight_dphase", in_group=pixc_group)[self.selected_index]
            # Sensitivity of pixel area to reference height
            self.darea_dheight = pixc_reader.get_var_value("darea_dheight", in_group=pixc_group)[self.selected_index]
            
            # Time of illumination of each pixel
            illumination_time = origin_illumination_time[self.selected_index]
            # Number of medium looks
            self.eff_num_medium_looks = pixc_reader.get_var_value("eff_num_medium_looks", in_group=pixc_group)[self.selected_index]
            # Sigma0
            self.sig0 = pixc_reader.get_var_value("sig0", in_group=pixc_group)[self.selected_index]
            # Sigma0 uncertainty
            self.sig0_u = pixc_reader.get_var_value("sig0_uncert", in_group=pixc_group)[self.selected_index]
                        
            # Dry troposphere vertical correction
            self.model_dry_tropo_cor = pixc_reader.get_var_value("model_dry_tropo_cor", in_group=pixc_group)[self.selected_index]
            # Wet troposphere vertical correction
            self.model_wet_tropo_cor = pixc_reader.get_var_value("model_wet_tropo_cor", in_group=pixc_group)[self.selected_index]
            # Ionosphere vertical correction
            self.iono_cor_gim_ka = pixc_reader.get_var_value("iono_cor_gim_ka", in_group=pixc_group)[self.selected_index]
            # Crossover calibration height correction
            self.height_cor_xover = pixc_reader.get_var_value("height_cor_xover", in_group=pixc_group)[self.selected_index]
            # Geoid
            self.geoid = pixc_reader.get_var_value("geoid", in_group=pixc_group)[self.selected_index]
            # Solid earth tide
            self.solid_earth_tide = pixc_reader.get_var_value("solid_earth_tide", in_group=pixc_group)[self.selected_index]
            # Load tide height (FES2014)
            self.load_tide_fes = pixc_reader.get_var_value("load_tide_fes", in_group=pixc_group)[self.selected_index]
            # Load tide height (GOT4.10)
            self.load_tide_got = pixc_reader.get_var_value("load_tide_got", in_group=pixc_group)[self.selected_index]
            # Pole tide height
            self.pole_tide = pixc_reader.get_var_value("pole_tide", in_group=pixc_group)[self.selected_index]
            
            # Quality flag for the interferogram quantities
            self.interferogram_qual = pixc_reader.get_var_value("interferogram_qual", in_group=pixc_group)[self.selected_index]
            # Quality flag for the classification quantities 
            self.classification_qual = origin_qual["classification_qual"][self.selected_index]
            # Quality flag for the geolocation quantities
            self.geolocation_qual = origin_qual["geolocation_qual"][self.selected_index]
            # Quality flag for sigma0
            self.sig0_qual = pixc_reader.get_var_value("sig0_qual", in_group=pixc_group)[self.selected_index]
            
            # Time difference from(to) the previous(next) crossover
            self.prev_xover = np.full(self.nb_selected, np.nan, dtype=np.float32)
            self.next_xover = np.full(self.nb_selected, np.nan, dtype=np.float32)
            # Try/except because time_from|to_prev|next_xover neither in PIXC generated by Large Scale Simulator nor in PIXC version < D2
            try:  # PIXC from version D2
                # Read variables in the PIXC
                # NB: data is given per azimuth_index
                time_from_prev_xover = pixc_reader.get_var_value("time_from_prev_xover", in_group=pixc_group)[self.selected_index]
                time_to_next_xover = pixc_reader.get_var_value("time_to_next_xover", in_group=pixc_group)[self.selected_index]
                # Set value for each selected pixels
                self.prev_xover = time_from_prev_xover[self.azimuth_index]
                self.next_xover = time_to_next_xover[self.azimuth_index]
            except:  # PIXC until version D0 + Large Scale Simulator
                logger.warning("There is neither time_from_prev_xover nor time_to_next_xover variable in the PIXC => set to NaN")

            # 6.2 - In TVP group
            
            # Link TVP and pixel_cloud indices
            # Try/except because pixc_line_to_tvp not in PIXC generated by Large Scale Simulator
            try:
                # Case when there is pixel_cloud.pixc_line_to_tvp variable
                nadir_index = pixc_line_to_tvp[self.azimuth_index].astype(int)
                
            except:
                # Case when there is no pixel_cloud.pixc_line_to_tvp variable
                if len(origin_nadir_time) > 1:
                    nadir_index = (np.rint(f(illumination_time))).astype(int)  # Link between PixC and nadir pixels
                else:
                    nadir_index = np.zeros(len(illumination_time), dtype=int)
                    
            # Nadir time
            self.nadir_time = origin_nadir_time[nadir_index]
            # Nadir time TAI
            self.nadir_time_tai = pixc_reader.get_var_value("time_tai", in_group=sensor_group)[nadir_index]
            # Nadir longitude
            self.nadir_longitude = my_tools.convert_to_m180_180(pixc_reader.get_var_value("longitude", in_group=sensor_group)[nadir_index])
            # Nadir latitude
            self.nadir_latitude = pixc_reader.get_var_value("latitude", in_group=sensor_group)[nadir_index]
            # Nadir cartesian coordinates
            self.nadir_x = pixc_reader.get_var_value("x", in_group=sensor_group)[nadir_index]
            self.nadir_y = pixc_reader.get_var_value("y", in_group=sensor_group)[nadir_index]
            self.nadir_z = pixc_reader.get_var_value("z", in_group=sensor_group)[nadir_index]
            # Nadir velocity in cartesian coordinates
            self.nadir_vx = pixc_reader.get_var_value("vx", in_group=sensor_group)[nadir_index]
            self.nadir_vy = pixc_reader.get_var_value("vy", in_group=sensor_group)[nadir_index]
            self.nadir_vz = pixc_reader.get_var_value("vz", in_group=sensor_group)[nadir_index]
            # Coordinates of plus_y antenna phase center in the ECEF frame
            self.nadir_plus_y_antenna_x = pixc_reader.get_var_value("plus_y_antenna_x", in_group=sensor_group)[nadir_index]
            self.nadir_plus_y_antenna_y = pixc_reader.get_var_value("plus_y_antenna_y", in_group=sensor_group)[nadir_index]
            self.nadir_plus_y_antenna_z = pixc_reader.get_var_value("plus_y_antenna_z", in_group=sensor_group)[nadir_index]
            # Coordinates of minus_y antenna phase center in the ECEF frame
            self.nadir_minus_y_antenna_x = pixc_reader.get_var_value("minus_y_antenna_x", in_group=sensor_group)[nadir_index]
            self.nadir_minus_y_antenna_y = pixc_reader.get_var_value("minus_y_antenna_y", in_group=sensor_group)[nadir_index]
            self.nadir_minus_y_antenna_z = pixc_reader.get_var_value("minus_y_antenna_z", in_group=sensor_group)[nadir_index]
            
            # 6.3 - Set bad PIXC height std to high number to deweight 
            # instead of giving infs/nans
            self.height_std_pix = np.abs(self.phase_noise_std * self.dheight_dphase)
            bad_num = 1.0e5
            self.height_std_pix[self.height_std_pix<=0] = bad_num
            self.height_std_pix[np.isinf(self.height_std_pix)] = bad_num
            self.height_std_pix[np.isnan(self.height_std_pix)] = bad_num
            
            # 6.4 - Compute height wrt the geoid and apply tide corrections
            # Compute indices of PIXC for which all corrections are valid
            valid_geoid = np.where(np.isfinite(self.geoid))[0]
            valid_solid_earth_tide = np.where(np.isfinite(self.solid_earth_tide))[0]
            valid_pole_tide = np.where(np.isfinite(self.pole_tide))[0]
            valid_load_tide_fes = np.where(np.isfinite(self.load_tide_fes))[0]
            inter1 = np.intersect1d(valid_geoid, valid_solid_earth_tide)
            inter2 = np.intersect1d(valid_pole_tide, inter1)
            ind_valid_corr = np.intersect1d(valid_load_tide_fes, inter2)
            # Compute corrected height for these PIXC
            self.corrected_height = np.full(self.nb_water_pix, np.nan, dtype=np.float32)
            self.corrected_height[ind_valid_corr] = self.height[ind_valid_corr] \
                                                    - self.geoid[ind_valid_corr] \
                                                    - self.solid_earth_tide[ind_valid_corr] \
                                                    - self.pole_tide[ind_valid_corr] \
                                                    - self.load_tide_fes[ind_valid_corr]
                    
            # 6.5 - Compute ice_dyn_f at pixel level
            self.ice_dyn_f = pixel_ice.compute_by_pixel(ice_model, self.inc, self.sig0, self.phase_noise_std, self.interferogram_qual, self.sig0_qual)
        
        # 7 - Close file
        pixc_reader.close()
        
    def select_pixc_from_criteria(self, in_index_reject,
                                  origin_longitude, origin_latitude, origin_height, origin_illumination_time,
                                  origin_data_window_first_valid, origin_data_window_last_valid,
                                  origin_bright_land_flag, origin_qual, 
                                  classif_list, pixc_qual_list, pixc_qual_mask, exclude_bright_land):
        """
        Select pixels to process given criteria:
            - Exclude river-only pixels
            - Exclude pixels with longitude, latitude, height, cross_track and illumination_time = numpy.nan
            - Exclude pixels with classification, range_index and azimuth_index = _FillValue
            - Exclude pixels that belong to the other swath side
            - Exclude pixels outside the data window
            - Exclude pixels corresponding to the unwanted quality masks
            - Exclude pixels having bright_land_flag = 1
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # PixC excluded from LakeTile processing will have their classification flag set to 100
        tmp_classif = self.origin_classif.copy()  # Init temporary classif vector
        
        if in_index_reject is not None:
            
            # 1.1 - pixc_index values from PIXCVECRiver MUST be less than PIXC.size-1 => exit with error if not the case
            index_error = np.where(in_index_reject > self.origin_classif.size-1)[0]
            nb_index_error = len(index_error)
            if nb_index_error > 0:
                logger.error("%d pixc_index values from PIXCVECRiver are larger than PIXC length (%d)" % (nb_index_error, self.origin_classif.size-1))
                logger.error(in_index_reject[index_error])
                logger.error("for PIXCVecRiver indices {}".format(index_error))
                raise service_error.ConfigError("Some pixc_index in input L2_HR_PIXCVecRiver are not compatible with input L2_HR_PIXC")
                
            # 1.2 - Exclude river-only pixels
            logger.debug("%d pixels are river (but not connected lakes) pixels => will be rejected" % in_index_reject.size)
            tmp_classif[in_index_reject] = 100
            
        # 2.1 - Exclude pixels with longitude, latitude, height and illumination_time = numpy.nan
        vars_to_look_for_nans = {}
        vars_to_look_for_nans["longitude"] = origin_longitude
        vars_to_look_for_nans["latitude"] = origin_latitude
        vars_to_look_for_nans["height"] = origin_height
        vars_to_look_for_nans["cross_track"] = self.origin_cross_track
        vars_to_look_for_nans["cross_track_refined"] = self.origin_cross_track_refined
        vars_to_look_for_nans["time"] = origin_illumination_time
        for var_name, var_value in vars_to_look_for_nans.items():
            nan_idx = np.argwhere(np.isnan(var_value))
            nb_nan = len(nan_idx)
            if nb_nan > 0:
                logger.debug("%d pixels have NaN in %s variable => will be rejected" % (nb_nan, var_name))
                tmp_classif[nan_idx] = 100 
        # 2.2 - Exclude pixels with classification, range_index and azimuth_index = _FillValue
        vars_to_look_for_fillvalue = {}
        vars_to_look_for_fillvalue["classification"] = self.origin_classif
        vars_to_look_for_fillvalue["range_index"] = self.origin_range_index
        vars_to_look_for_fillvalue["azimuth_index"] = self.origin_azimuth_index
        for var_name, var_value in vars_to_look_for_fillvalue.items():
            if var_name == "classification":
                max_value = 10.
            else:
                max_value = 999999.
            fv_idx = np.argwhere(var_value > max_value)
            nb_fv = len(fv_idx)
            if nb_fv > 0:
                logger.debug("%d pixels have _FillValue in %s variable => will be rejected" % (nb_fv, var_name))
                tmp_classif[fv_idx] = 100 
                
        # 3.1 - Exclude pixels that belong to the other swath side
        if self.pixc_metadata["swath_side"] == "L":  # Left swath
            nok_idx1 = np.where(self.origin_cross_track > 0)[0]  # Wrt cross_track variable
            nb_nok1 = len(nok_idx1)
            nok_idx2 = np.where(self.origin_cross_track_refined > 0)[0]  # Wrt cross_track_refined variable
            nb_nok2 = len(nok_idx2)
            nok_idx = np.unique(np.concatenate((nok_idx1, nok_idx2)))
            nb_nok = len(nok_idx)
            if nb_nok > 0:
                if nb_nok1 > 0:
                    logger.debug("%d pixels have cross_track > 0 while swath_side=L" % nb_nok1)
                else:
                    logger.debug("NO pixel has cross_track > 0 while swath_side=L")
                if nb_nok2 > 0:
                    logger.debug("%d pixels have cross_track_refined > 0 while swath_side=L" % nb_nok2)
                else:
                    logger.debug("NO pixel has cross_track_refined > 0 while swath_side=L")
                logger.debug("=> %d pixels will be rejected" % nb_nok)
                tmp_classif[nok_idx] = 100 
        elif self.pixc_metadata["swath_side"] == "R":  # Right swath
            nok_idx1 = np.where(self.origin_cross_track < 0)[0]  # Wrt cross_track variable
            nb_nok1 = len(nok_idx1)
            nok_idx2 = np.where(self.origin_cross_track_refined < 0)[0]  # Wrt cross_track_refined variable
            nb_nok2 = len(nok_idx2)
            nok_idx = np.unique(np.concatenate((nok_idx1, nok_idx2)))
            nb_nok = len(nok_idx)
            if nb_nok > 0:
                if nb_nok1 > 0:
                    logger.debug("%d pixels have cross_track < 0 while swath_side=R" % nb_nok1)
                else:
                    logger.debug("NO pixel has cross_track < 0 while swath_side=R")
                if nb_nok2 > 0:
                    logger.debug("%d pixels have cross_track_refined < 0 while swath_side=R" % nb_nok2)
                else:
                    logger.debug("NO pixel has cross_track_refined < 0 while swath_side=R")
                logger.debug("=> %d pixels will be rejected" % nb_nok)
                tmp_classif[nok_idx] = 100 
        # 3.2 - Exclude pixels outside the data window
        for cur_az_index in np.arange(len(origin_data_window_first_valid)):
            pixc_idx = np.where(self.origin_azimuth_index == cur_az_index)[0]
            if len(pixc_idx) > 0:
                
                # First range index
                first_range = origin_data_window_first_valid[cur_az_index]
                # Outside the data window
                nok_idx = np.where(self.origin_range_index[pixc_idx] < first_range)[0]
                nb_nok = len(nok_idx)
                if nb_nok > 0:
                    logger.debug("%d pixels are outside the data window (range_index < %d for azimuth_index = %d) => will be rejected" % (nb_nok, first_range, cur_az_index))
                    tmp_classif[pixc_idx[nok_idx]] = 100 
                # At the edge of the data window
                edge_idx = np.where(self.origin_range_index[pixc_idx] == first_range)[0]
                if edge_idx.size > 0:
                    #logger.debug("Pixel at the NR edge of the data window (range_index_min = %d for azimuth_index = %d)" % (first_range, cur_az_index))
                    self.origin_pixel_at_window_edge[pixc_idx[edge_idx]] = 1
                
                # Last range index
                last_range = origin_data_window_last_valid[cur_az_index]
                # Outside the data window
                nok_idx = np.where(self.origin_range_index[pixc_idx] > last_range)[0]
                nb_nok = len(nok_idx)
                if nb_nok > 0:
                    logger.debug("%d pixels are outside the data window (range_index > %d for azimuth_index = %d) => will be rejected" % (nb_nok, last_range, cur_az_index))
                    tmp_classif[pixc_idx[nok_idx]] = 100 
                # At the edge of the data window
                edge_idx = np.where(self.origin_range_index[pixc_idx] == last_range)[0]
                if edge_idx.size > 0:
                    #logger.debug("Pixel at the FR edge of the data window (range_index_max = %d for azimuth_index = %d)" % (last_range, cur_az_index))
                    self.origin_pixel_at_window_edge[pixc_idx[edge_idx]] = 1
        logger.debug("%d pixels are at the edge of the data window" % np.sum(self.origin_pixel_at_window_edge))
                    
        # 4 - Exclude pixels having classification_qual.detected_water_but_no_prior_water=1
        # and classification_qual.specular_ringing_degraded=1
        tmp_comp = np.bitwise_and(origin_qual["classification_qual"], 524292)
        nok_idx = np.where(tmp_comp == 524292)[0]
        nb_nok = len(nok_idx)
        if nb_nok == 0:
            logger.debug("No pixel has classification_qual.detected_water_but_no_prior_water=1 AND classification_qual.specular_ringing_degraded=1 => ALL pixels KEPT")
        else:
            logger.debug("%d pixels have classification_qual.detected_water_but_no_prior_water=1 AND classification_qual.specular_ringing_degraded=1 => will be rejected" % nb_nok)
            tmp_classif[nok_idx] = 100 

        # 5 - Exclude pixels having classification_qual.detected_water_but_no_prior_water=1
        # and classification_qual.water_false_detection_rate_suspect=1 and pixel 5, 6 or 7
        tmp_comp = np.bitwise_and(origin_qual["classification_qual"], 20)
        nok_idx = np.where((tmp_comp == 20) & ((tmp_classif == 5) | (tmp_classif == 6) | (tmp_classif == 7)))[0]
        nb_nok = len(nok_idx)
        if nb_nok == 0:
            logger.debug("No pixel has classification_qual.detected_water_but_no_prior_water=1 AND classification_qual.water_false_detection_rate_suspect=1 AND is of classif 5, 6 or 7 => ALL pixels KEPT")
        else:
            logger.debug("%d pixels have classification_qual.detected_water_but_no_prior_water=1 AND classification_qual.water_false_detection_rate_suspect=1 AND is of classif 5, 6 or 7 => will be rejected" % nb_nok)
            tmp_classif[nok_idx] = 100 
                
        # 6 - Exclude pixels corresponding to the unwanted quality masks
        if pixc_qual_list is not None:
            for flag_name, flag_mask in zip(pixc_qual_list, pixc_qual_mask):
                if int(flag_mask) == 0:
                    nok_idx = np.argwhere(origin_qual[flag_name] > 0)
                else:
                    mask_nok = np.uint32(int(flag_mask, 2))
                    tmp_comp = np.bitwise_and(origin_qual[flag_name], mask_nok)
                    nok_idx = np.argwhere(tmp_comp > 0)
                nb_nok = len(nok_idx)
                if nb_nok == 0:
                    logger.debug("ALL pixels are kept after filtering over %s quality flag" % flag_name)
                else:
                    if int(flag_mask) == 0:
                        logger.debug("%d pixels aren't totally good [value!=0] for %s quality flag => will be rejected" % (nb_nok, flag_name))
                    else:
                        logger.debug("%d pixels correspond to the unwanted quality mask [%s=%d] for %s quality flag => will be rejected" \
                                     % (nb_nok, flag_mask, mask_nok, flag_name))
                    tmp_classif[nok_idx] = 100 
                    
        # 7 - Exclude pixels having bright_land_flag = 1
        if exclude_bright_land:
            nok_idx = np.argwhere(origin_bright_land_flag == 1)
            nb_nok = len(nok_idx)
            if nb_nok == 0:
                logger.debug("No pixel has bright_land_flag=1 => ALL pixels KEPT")
            else:
                logger.debug("%d pixels have bright_land_flag=1 => will be rejected" % nb_nok)
                tmp_classif[nok_idx] = 100 
        
        # 8 - Get indices of selected pixels
        self.selected_index = None  # Init wanted indices vector
        for classif_flag in classif_list:
            v_ind = np.where(tmp_classif == classif_flag)[0]
            logger.debug("%d pixels with classification flag = %d" % (v_ind.size, classif_flag))
            if v_ind.size > 0:
                if self.selected_index is None:
                    self.selected_index = v_ind
                else:
                    self.selected_index = np.concatenate((self.selected_index, v_ind))
        if self.selected_index is None:
            self.nb_selected = 0
        else:
            self.nb_selected = self.selected_index.size
        logger.debug("=> %d pixels to keep" % self.nb_selected)
        
    def init_classif_dict(self, in_classif_list, in_nb_pixc, in_classif, classif_qual=None):
        """
        Initialize self.classif_dict wrt classification flags set as input parameters
        
        :param in_classif_list: list of classif flags to keep for processing
        :type in_classif_list: list of int
        :param in_nb_pixc: number of selected pixels
        :type in_nb_pixc: int
        :param in_classif: classes of selected pixels
        :type in_classif: 1D-array of int
        :param classif_qual: classification quality flags of selected pixels
        :type classif_qual: 1D-array of int
        
        :return: out_classif_dict = True/False flags wrt key info
        :rtype: out_classif_dict = dict of 1D-arrays of boolean
        """
        
        # PIXC indices wrt category of process
        out_classif_dict = dict()
        for key in self.list_classif_keys:
            out_classif_dict[key] = np.full((in_nb_pixc), False)
        
        # Boolean arrays of classification flags (=True if PIXC is of classif "key")
        tmp_classif_dict = dict()
        mask_interior_water = np.full((in_nb_pixc), False)
        classif_mandatory = set([int(i) for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_MANDATORY").split(";")])  # Set of classif flags a lake feature must contain at least
        for classif_flag in in_classif_list:
            tmp_classif_dict[classif_flag] = np.full((in_nb_pixc), False)
            tmp_classif_dict[classif_flag][in_classif == classif_flag] = True
            if classif_flag in set(my_var.CLASSIF_INTERIOR_WATER).intersection(set(in_classif_list)):
                mask_interior_water += tmp_classif_dict[classif_flag]
            if classif_flag in set(my_var.CLASSIF_WATER_EDGE).intersection(set(in_classif_list)):
                out_classif_dict["water_edge"] += tmp_classif_dict[classif_flag] 
            if classif_flag in set(my_var.CLASSIF_OPEN_WATER).intersection(set(in_classif_list)):
                out_classif_dict["open_water"] += tmp_classif_dict[classif_flag]
            if classif_flag in set(my_var.CLASSIF_LOW).intersection(set(in_classif_list)):
                out_classif_dict["low_coherence"] += tmp_classif_dict[classif_flag]
            if classif_flag in classif_mandatory.intersection(set(in_classif_list)):
                out_classif_dict["mandatory"] += tmp_classif_dict[classif_flag]
                
        # Boolean array for specular_ringing flag activated (=True if PIXC has classification_qual bit 19 is set to 1)
        mask_specular_ringing = np.full((in_nb_pixc), False)
        if classif_qual is not None:
            mask_specular_ringing = np.uint32(2**19)
            tmp_comp = np.bitwise_and(classif_qual, mask_specular_ringing)
            mask_specular_ringing = tmp_comp.astype(bool)
        out_classif_dict["specular_ringing"] = mask_specular_ringing
        out_classif_dict["not_specular_ringing"] = np.invert(mask_specular_ringing)
                
        # PIXC indices of water and dark water pixels
        classif_water = [i == "True" for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_WATER").split(";")]
        for classif_flag, ok_water in zip(in_classif_list, classif_water):
            if ok_water:
                out_classif_dict["water"] += tmp_classif_dict[classif_flag]
            else:
                out_classif_dict["dark"] += tmp_classif_dict[classif_flag]
                
        # PIXC indices of pixels to be used for hull computation
        classif_4hull = [i == "True" for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_4HULL").split(";")]
        for classif_flag, ok_4hull in zip(in_classif_list, classif_4hull):
            if ok_4hull:
                out_classif_dict["4hull"] += tmp_classif_dict[classif_flag]
                
        # PIXC indices of pixels to be used for wse computation
        classif_4wse = [i == "True" for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_4WSE").split(";")]
        for classif_flag, ok_4wse in zip(in_classif_list, classif_4wse):
            if ok_4wse:
                out_classif_dict["4wse"] += tmp_classif_dict[classif_flag]
        # Remove specular_ringing pixels
        out_classif_dict["4wse"] = out_classif_dict["4wse"] & out_classif_dict["not_specular_ringing"]
                
        # PIXC indices of pixels to be used for area_total attribute computation
        classif_4area = [i == "True" for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_4AREA").split(";")]
        for classif_flag, ok_4area in zip(in_classif_list, classif_4area):
            if ok_4area:
                out_classif_dict["4area"] += tmp_classif_dict[classif_flag]
                
        # Binary detected-water mask for pixels detected as interior water (to be further used by JPL aggregate.py library)
        # Without dark water
        out_classif_dict["Idw_in_no_dw"] = np.copy(mask_interior_water) 
        # With dark water
        out_classif_dict["Idw_in"] = np.copy(out_classif_dict["Idw_in_no_dw"])  
        for classif_flag in set(my_var.CLASSIF_DARK).intersection(set(in_classif_list)):
            out_classif_dict["Idw_in"] += tmp_classif_dict[classif_flag]
        
        # Binary detected-water mask representing pixels detected as water (to be further used by JPL aggregate.py library)
        # Without dark water
        out_classif_dict["Idw_no_dw"] = np.full((in_nb_pixc), False)  
        for classif_flag in set(my_var.CLASSIF_INTERIOR_WATER + my_var.CLASSIF_WATER_EDGE).intersection(set(in_classif_list)):
            out_classif_dict["Idw_no_dw"] += tmp_classif_dict[classif_flag] 
        # With dark water
        out_classif_dict["Idw"] = np.copy(out_classif_dict["Idw_no_dw"])  
        for classif_flag in set(my_var.CLASSIF_DARK).intersection(set(in_classif_list)):
            out_classif_dict["Idw"] += tmp_classif_dict[classif_flag]
        
        # Binary detected-water mask for pixels detected as edge water (to be further used by JPL aggregate.py library)
        out_classif_dict["Ide"] = np.full((in_nb_pixc), False)  
        for classif_flag in set(my_var.CLASSIF_WATER_EDGE + my_var.CLASSIF_LAND_EDGE).intersection(set(in_classif_list)):
            out_classif_dict["Ide"] += tmp_classif_dict[classif_flag] 
        
        # Binary mask for pixels detected as water or land-near-water (to be further used by JPL aggregate.py library)
        #out_classif_dict["I"] = out_classif_dict["Idw_in"] + out_classif_dict["Idw"] + out_classif_dict["Ide"] => no need Idw?
        out_classif_dict["I"] = out_classif_dict["Idw_in"] + out_classif_dict["Ide"]
        # Idw_no_dw out_classif_dict["I_no_dw"] = out_classif_dict["Idw_in_no_dw"] + out_classif_dict["Idw_no_dw"] + out_classif_dict["Ide"]
        out_classif_dict["I_no_dw"] = out_classif_dict["Idw_in_no_dw"] + out_classif_dict["Ide"]
        
        # Convert from boolean to int
        for key in ["Idw_in", "Idw_in_no_dw", "Idw", "Idw_no_dw", "Ide", "I", "I_no_dw"]:
            out_classif_dict[key] = out_classif_dict[key]*1
            
        return out_classif_dict

    # ----------------------------------------
    
    def remove_pixels_outside_window(self, in_out_pixc_index):
        """
        Considering the pixels related to the water feature given by in_pixc_index,
        remove pixels outside given min and max cross-track distance and pixel_cloud.data_window_first|last_valid variables
        and compute partial flag depending if the water feature crosses right or left swath edges.
        
        :param in_out_pixc_index: indices of pixels defining the studied feature
        :type in_out_pixc_index: 1D-array of int
        
        :return: in_out_pixc_index = indices of remaining pixels defining the studied feature
        :rtype: in_out_pixc_index = 1D-array of int
        :return: out_flag_partial =1 if the studied feature crosses right or left swath edges; =0 otherwise
        :rtype: out_flag_partial = int
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 0 - Init output variable
        out_flag_partial = 0
        
        # 1 - Near range crosstrack distance
        min_xtrack = self.cfg.getfloat('CONFIG_PARAMS', 'MIN_XTRACK')
        if min_xtrack > 0.:
            min_xtrack_idx1 = np.where(np.abs(self.cross_track[in_out_pixc_index]) < min_xtrack)[0]
            nb_min_xtrack1 = len(min_xtrack_idx1)
            min_xtrack_idx2 = np.where(np.abs(self.cross_track_refined[in_out_pixc_index]) < min_xtrack)[0]
            nb_min_xtrack2 = len(min_xtrack_idx2)
            min_xtrack_idx = np.unique(np.concatenate((min_xtrack_idx1, min_xtrack_idx2)))
            nb_min_xtrack = len(min_xtrack_idx)
            if nb_min_xtrack > 0:
                if nb_min_xtrack1 > 0:
                    logger.debug("%d pixels have cross_track < %d m" % (nb_min_xtrack1, round(min_xtrack)))
                else:
                    logger.debug("NO pixel has cross_track < %d m" % round(min_xtrack))
                if nb_min_xtrack2 > 0:
                    logger.debug("%d pixels have cross_track_refined < %d m" % (nb_min_xtrack2, round(min_xtrack)))
                else:
                    logger.debug("NO pixel has cross_track_refined < %d m" % round(min_xtrack))
                logger.debug("=> %d pixels will be rejected" % nb_min_xtrack)
                out_flag_partial = 1
                in_out_pixc_index = np.delete(in_out_pixc_index, min_xtrack_idx)
                
        # 2 - Far range crosstrack distance
        max_xtrack = self.cfg.getfloat('CONFIG_PARAMS', 'MAX_XTRACK')
        if max_xtrack > 0.:
            max_xtrack_idx1 = np.where(np.abs(self.cross_track[in_out_pixc_index]) > max_xtrack)[0]
            nb_max_xtrack1 = len(max_xtrack_idx1)
            max_xtrack_idx2 = np.where(np.abs(self.cross_track_refined[in_out_pixc_index]) > max_xtrack)[0]
            nb_max_xtrack2 = len(max_xtrack_idx2)
            max_xtrack_idx = np.unique(np.concatenate((max_xtrack_idx1, max_xtrack_idx2)))
            nb_max_xtrack = len(max_xtrack_idx)
            if nb_max_xtrack > 0:
                if nb_max_xtrack1 > 0:
                    logger.debug("%d pixels have cross_track > %d m" % (nb_max_xtrack1, round(max_xtrack)))
                else:
                    logger.debug("NO pixel has cross_track > %d m" % round(max_xtrack))
                if nb_max_xtrack2 > 0:
                    logger.debug("%d pixels have cross_track_refined > %d m" % (nb_max_xtrack2, round(max_xtrack)))
                else:
                    logger.debug("NO pixel has cross_track_refined > %d m" % round(max_xtrack))
                logger.debug("=> %d pixels will be rejected" % nb_max_xtrack)
                out_flag_partial = 1
                in_out_pixc_index = np.delete(in_out_pixc_index, max_xtrack_idx)
                
        # 3 - Consider the range indices of the valid data window
        if np.sum(self.pixel_at_window_edge[in_out_pixc_index]) > 0:
            out_flag_partial = 1
                
        return in_out_pixc_index, out_flag_partial

    # ----------------------------------------

    def compute_water_mask(self):
        """
        Create the water mask (i.e. a 2D binary matrix) in radar geometry,
        from the pixel cloud (1D-array layers of azimuth_index, range_index, classification and continuous classification)
        
        :return: water mask in radar geometry, i.e. a 2D matrix with "1" for each pixel in input
        :rtype: 2D binary matrix of int 0/1
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")

        return my_tools.compute_bin_mat(self.nb_pix_range, self.nb_pix_azimuth, self.range_index, self.azimuth_index)

    def compute_separate_regions(self):
        """
        Identify all separate water regions in the water mask
        """
        cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        timer = my_timer.Timer()
        timer.start()

        # 1 - Create the water mask
        water_mask = self.compute_water_mask()
        if cfg.getboolean("OPTIONS", "Save intermediate data"):
            out_filename = "INTER_water_mask_%03d_%s.npy" % (self.pixc_metadata["cycle_number"], self.pixc_metadata["tile_name"])
            out_file = os.path.join(cfg.get("PATHS", "Output directory"), out_filename)
            np.save(out_file, water_mask)

        # 2 - Identify all separate water regions in a 2D binary mask
        sep_regions, self.nb_obj = my_tools.label_region(water_mask)
        if cfg.getboolean("OPTIONS", "Save intermediate data"):
            out_filename = "INTER_label_region_%03d_%s.npy" % (self.pixc_metadata["cycle_number"], self.pixc_metadata["tile_name"])
            out_file = os.path.join(cfg.get("PATHS", "Output directory"), out_filename)
            np.save(out_file, sep_regions)

        # 3 - Convert 2D labelled mask in 1D-array layer of the same size of the L2_HR_PIXC layers
        self.labels = my_tools.convert_2d_mat_in_1d_vec(self.range_index, self.azimuth_index, sep_regions)
        self.labels = self.labels.astype(int)  # Conversion from float to integer

        # 4 - Relabel Lake Using Segmentation Heigth
        # For each label : check if only one lake is in each label and relabels if necessary

        # 4.0. Check if lake segmentation following height needs to be run
        seg_method = cfg.getint('CONFIG_PARAMS', 'SEGMENTATION_METHOD')
        min_object_size = cfg.getfloat('CONFIG_PARAMS', 'MIN_SIZE') * 1e6
        min_object_nb_pixc = cfg.getint('CONFIG_PARAMS', 'MIN_NB_PIXC')

        # If SEGMENTATION_METHOD is 0, function unactivated
        logger.debug("Lake segmentation following height : %s" %(my_segmentation.segmentation_method_name[seg_method]))
        # 4.1. If SEGMENTATION_METHOD not 0, relabel lake following segmentation height
        if seg_method != 0:
            labels_tmp = np.zeros(self.labels.shape, dtype=self.labels.dtype)
            labels, count = np.unique(self.labels, return_counts=True)

            for i, label in enumerate(labels):
                idx = np.where(self.labels == label)
                subset_pixel_area = self.pixel_area[idx]

                if np.sum(subset_pixel_area) > min_object_size*2:
                    min_rg = min(self.range_index[idx])
                    min_az = min(self.azimuth_index[idx])
                    subset_range = self.range_index[idx] - min_rg
                    subset_azimuth = self.azimuth_index[idx] - min_az
                    subset_corrected_height = self.corrected_height[idx]

                    relabel_obj = my_segmentation.relabel_lake_using_segmentation_height(subset_range, subset_azimuth, subset_corrected_height,
                                                                                         subset_pixel_area, min_object_size, min_object_nb_pixc,
                                                                                         seg_method)

                    labels_tmp[self.labels == label] = np.max(labels_tmp) + relabel_obj
                else:
                    labels_tmp[self.labels == label] = np.max(labels_tmp) + 1

            self.labels = labels_tmp
            
        self.labels_unique, self.labels_counts = np.unique(self.labels, return_counts=True)
        self.nb_obj = self.labels_unique.size
        
        logger.debug(timer.stop(message="compute_separate_regions computation time : "))

    def localize_regions_wrt_tile(self, obj_lake_db):
        """
        Separate labels of water features entirely inside the tile, from labels of features at top and/or bottom of the tile.
        Objects intersecting PLD lakes located at the edge of the tile are not considered as 'inside the tile'.

        :param obj_lake_db: Prior Lake Database (PLD) object
        :type obj_lake_db: lib_lake.lake_db
        """
        logger = logging.getLogger(self.__class__.__name__)
        flag_exclude_coastal_zones = self.cfg.getboolean("CONFIG_PARAMS", "EXCLUDE_COASTAL_ZONES")
        min_nb_pixc_4label = self.cfg.getint("CONFIG_PARAMS", "MIN_NB_PIXC")
                
        # 1 - Identify labels related to PLD lakes intersecting az_0_geom, az_max_geom or az_0_and_max_geom
        # NB: do not consider too small features to save time computation
        labels_over_min_nb_pixc = self.labels_unique[self.labels_counts >= min_nb_pixc_4label]
        labels_linked_to_pld = set()
        for label in labels_over_min_nb_pixc.tolist():
            idx = np.where(self.labels == label)[0]

            # 1.1 - Polygon of the bounding box around the pixels related to the current label
            bb_geom = my_tools.get_bounding_box(np.min(self.latitude[idx]), np.max(self.latitude[idx]),
                                                np.min(self.longitude[idx]), np.max(self.longitude[idx]))
            
            # 1.2 - Test its intersection with the top and the bottom edge of the PIXC tile
            # NB: test simple bounding box before exact convex hull to save time computation
            #     (no need to compute convex hull for each water region)
            convex_geom = None

            flag_bottom_ok = False
            flag_top_ok = False
            flag_top_and_bottom_ok = False
            
            # /!\ The "if" below MUST be kept separated
            # If ... elif ... elif ... is not sufficient /!\
            
            # 1.2.1 - Test if current geom is at the bottom of the tile
            if (not self.seg_hr_start) and (obj_lake_db.az_0_geom is not None) and bb_geom.Intersects(obj_lake_db.az_0_geom):
                convex_geom = my_hull.get_convex_hull(self.longitude[idx], self.latitude[idx])
                if convex_geom.Intersects(obj_lake_db.az_0_geom):
                    flag_bottom_ok = True

            # 1.2.2 - Test if current geom is at the top of the tile
            if (not self.seg_hr_stop) and (obj_lake_db.az_max_geom is not None) and bb_geom.Intersects(obj_lake_db.az_max_geom):
                if convex_geom is None:
                    convex_geom = my_hull.get_convex_hull(self.longitude[idx], self.latitude[idx])
                if convex_geom.Intersects(obj_lake_db.az_max_geom):
                    flag_top_ok = True
                    
            # 1.2.3 - Test if current geom is at the top and the bottom of the tile
            if (obj_lake_db.az_0_and_max_geom is not None) and bb_geom.Intersects(obj_lake_db.az_0_and_max_geom):
                if convex_geom is None:
                    convex_geom = my_hull.get_convex_hull(self.longitude[idx], self.latitude[idx])
                if convex_geom.Intersects(obj_lake_db.az_0_and_max_geom):
                    if self.seg_hr_start:
                        if not self.seg_hr_stop:
                            flag_top_ok = True
                    else:
                        if self.seg_hr_stop:
                            flag_bottom_ok = True
                        else:
                            flag_top_and_bottom_ok = True

            # 1.3 - Add to lists depending on the result
            if flag_top_and_bottom_ok:
                self.labels_at_both_edges.add(label)
                labels_linked_to_pld.add(label)
            elif flag_bottom_ok:
                self.labels_at_bottom_edge.add(label)
                labels_linked_to_pld.add(label)
            elif flag_top_ok:
                self.labels_at_top_edge.add(label)
                labels_linked_to_pld.add(label)
                
        logger.debug("{} labels related to PLD lakes crossing the edges of the tile".format(len(labels_linked_to_pld)))
        logger.debug("{} labels related to PLD lakes crossing bottom edge of the tile = {}".format(len(self.labels_at_bottom_edge,), self.labels_at_bottom_edge))
        logger.debug("{} labels related to PLD lakes crossing top edge of the tile = {}".format(len(self.labels_at_top_edge), self.labels_at_top_edge))
        logger.debug("{} labels related to PLD lakes crossing both edges of the tile = {}".format(len(self.labels_at_both_edges), self.labels_at_both_edges))
        
        # 2 - Identify other labels stuck to first or last azimuth lines

        # 2.1 - Identify objects at azimuth = self.pixc_first_line_index
        if self.seg_hr_start:
            labels_at_az_0 = np.array([], dtype=np.int64)
        else:
            idx_at_az_0 = np.where(self.azimuth_index == self.pixc_first_line_index)[0]
            labels_at_az_0 = np.unique(self.labels[idx_at_az_0])

        # 2.2 - Identify objects at azimuth = self.pixc_last_line_index
        if self.seg_hr_stop:
            labels_at_az_max = np.array([], dtype=np.int64)
        else:
            idx_at_az_max = np.where(self.azimuth_index == self.pixc_last_line_index)[0]
            labels_at_az_max = np.unique(self.labels[idx_at_az_max])

        # 2.3.1 - Identify objects at azimuth = 0 and azimuth = max
        tmp_labels_at_both_edges = np.intersect1d(labels_at_az_0, labels_at_az_max)
        logger.debug("Labels identified at az_0+az_max = {}".format(tmp_labels_at_both_edges))
        # 2.3.2 - Remove objects at both edges from labels_at_az_0
        tmp_labels_at_bottom_edge = np.setdiff1d(labels_at_az_0, tmp_labels_at_both_edges)
        logger.debug("Labels identified at az_0 = {}".format(tmp_labels_at_bottom_edge))
        # 2.3.3 - Remove objects at both edges from labels_at_az_max
        tmp_labels_at_top_edge = np.setdiff1d(labels_at_az_max, tmp_labels_at_both_edges)
        logger.debug("Labels identified at az_max = {}".format(tmp_labels_at_top_edge))

        # 3 - Add to arrays previously computed
        # 3.1 - Merging labels intersecting az_0_geom with labels at azimuth = self.pixc_first_line_index
        self.labels_at_bottom_edge.update(set(tmp_labels_at_bottom_edge))
        # 3.2 - Merging labels intersecting az_max_geom with labels at  azimuth = self.pixc_last_line_index
        self.labels_at_top_edge.update(set(tmp_labels_at_top_edge))
        # 3.3 - Merging labels intersecting az_0_and_max_geom with labels at both at azimuth = self.pixc_first_line_index and azimuth = self.pixc_last_line_index
        self.labels_at_both_edges.update(set(tmp_labels_at_both_edges))
        
        # 4 - List edge features that should be processed in LakeTile anyway
        labels_processed = self.labels_at_bottom_edge | self.labels_at_top_edge | self.labels_at_both_edges
        labels_to_remove = set()
        labels_to_remove_from_processing = set()
        for label in labels_processed:
            
            idx = np.where(self.labels == label)[0]
            lon_tmp = self.longitude[idx]
            lat_tmp = self.latitude[idx]
            
            # TODO: remove when no remaining over-detection in PIXC
            # 4.1 - Remove coastal features which have to be processed in lakeTile
            # NB: as they generally cross azimuth edges, they are listed in self.labels_at_edge lists
            if obj_lake_db.is_coastal(lon_tmp, lat_tmp):
                if flag_exclude_coastal_zones:
                    logger.debug("Label %s corresponds to a coastal zone => totally removed from Lake processing" % label)
                    labels_to_remove_from_processing.add(label)
                else:
                    logger.debug("Label %s corresponds to a coastal zone => processed in LakeTile" % label)
                labels_to_remove.add(label)
                self.labels_coastal_zones.add(label)

            # 4.2 - Identify labels of objects linked to PLD lakes such as pld_overlap_avg_[l/r]*pld_nb_tiles_[l/r] > MAX_NB_TILES_FULL_AZ
            elif len(obj_lake_db.big_lake_geom) > 0:
                # Create memory layer to store PixC points
                ds, lyr = my_tools.load_pixels_to_mem_layer(lon_tmp, lat_tmp)
                # Compare to geometry of PLD lakes at azimuth edges
                for lake_id, geom in obj_lake_db.big_lake_geom:
                    lyr.SetSpatialFilter(geom)
                    if lyr.GetFeatureCount() > 0 :
                        logger.debug("Label %s = PLD lake %s which is too big => processed in LakeTile even if crosses azimuth border" \
                                     % (label, lake_id))
                        self.labels_big_lakes.add(label)
                        self.big_lakeid.add(lake_id)
                    lyr.SetSpatialFilter(None)
                ds.Destroy()
                
            # 4.3 - Identify labels of objects covering more than 1 million of pixels
            # and remove them from self.labels_at_edge lists if not linked to PLD lake
            elif idx.size > 1e6 :
                # TODO: remove comment when use of MAX_NB_TILES_FULL_AZ improved
                """
                if label in labels_linked_to_pld:
                    logger.debug("Label %s corresponds to a big water body (%d pixels) BUT linked to a PLD lake not considered as big "
                             "=> its pixels go to LakeTile_Edge file and it will be processed in LakeSP" % (label, idx.size))
                else:
                    logger.debug("Label %s corresponds to a big water body (%d pixels) NEITHER linked to a PLD lake NOR a coastal area " 
                                "=> processed in LakeTile" % (label, idx.size))
                """
                logger.debug("Label %s corresponds to a big water body (%d pixels) NEITHER linked to a large PLD lake NOR a coastal area " 
                                "=> processed in LakeTile" % (label, idx.size))
                labels_to_remove.add(label)
                
        # 4.5 - Remove listed labels from lists of edge features
        labels_to_remove_np = np.array(list(labels_to_remove), dtype=int)
        self.labels_at_both_edges = np.setdiff1d(np.array(list(self.labels_at_both_edges), dtype=int), labels_to_remove_np)
        self.labels_at_top_edge = np.setdiff1d(np.array(list(self.labels_at_top_edge), dtype=int), labels_to_remove_np)
        self.labels_at_bottom_edge = np.setdiff1d(np.array(list(self.labels_at_bottom_edge), dtype=int), labels_to_remove_np)

        # 5 - Count features and print info
        self.nb_obj_at_both_edges = self.labels_at_both_edges.size
        logger.debug("> %d labels at bottom (az=%d) AND top (az=%d) of the tile" % (self.nb_obj_at_both_edges, self.pixc_first_line_index, self.pixc_last_line_index))
        if self.nb_obj_at_both_edges > 0:
            logger.debug(";".join(list(self.labels_at_both_edges.astype('str'))))

        self.nb_obj_at_top_edge = self.labels_at_top_edge.size
        logger.debug("> %d labels at top of the tile (az=%d)" % (self.nb_obj_at_top_edge, self.pixc_last_line_index))
        if self.nb_obj_at_top_edge > 0:
            logger.debug(";".join(list(self.labels_at_top_edge.astype('str'))))

        self.nb_obj_at_bottom_edge = self.labels_at_bottom_edge.size
        logger.debug("> %d labels at bottom of the tile (az=%d)" % (self.nb_obj_at_bottom_edge, self.pixc_first_line_index))
        if self.nb_obj_at_bottom_edge > 0:
            logger.debug(";".join(list(self.labels_at_bottom_edge.astype('str'))))

        # 6 - Get labels of objects entirely inside the tile
        self.labels_inside = np.copy(self.labels_unique)  # Initialisation
        self.labels_inside = np.setdiff1d(self.labels_inside, self.labels_at_bottom_edge)  # Delete labels at bottom
        self.labels_inside = np.setdiff1d(self.labels_inside, self.labels_at_top_edge)  # Delete labels at top
        self.labels_inside = np.setdiff1d(self.labels_inside, self.labels_at_both_edges)  # Delete labels at top and bottom
        self.labels_inside = np.setdiff1d(self.labels_inside, np.array(list(labels_to_remove_from_processing)))  # Remove coastal features if specified
        self.labels_inside = np.concatenate((self.labels_inside, np.array(list(self.labels_big_lakes), dtype=int)))  # Process big lakes in LakeTile
        self.nb_obj_inside = self.labels_inside.size
        logger.debug("> %d objects entirely inside the tile" % self.nb_obj_inside)
        #logger.debug(";".join(list(self.labels_inside.astype('str'))))
        
        # 6 - Clean labels inside
        self.clean_labels()
        self.nb_obj_inside = self.labels_inside.size
        logger.debug("> %d objects entirely inside the tile after cleaning" % self.nb_obj_inside)
        #logger.debug(";".join(list(self.labels_inside.astype('str'))))
        
    def clean_labels(self):
        """
        For further processing, keep only labels containing at least one mandatory pixel, 
        i.e. of the mandatory class, and not flagged specular_ringing
        and containing at least 5 pixels (NB: 1 pixel = 23m*70m = 1610m2 max => 6.2 pixels min for 1ha)
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Determine if a pixel belongs to the mandatory class
        # Remainder: pixel should be of mandatory class and not flagged specular_ringing
        idx_mandatory = self.classif_dict["mandatory"] & self.classif_dict["not_specular_ringing"]
        
        # 2 - Retrieve the labels of these pixels and compute the list of single labels
        remaining_labels, count_labels = np.unique(self.labels[idx_mandatory], return_counts=True)
        logger.debug("> %d labels contain at least 1 mandatory pixel which is not flagged specular_ringing" % remaining_labels.size)
        
        # 3 - Remove the labels containing less than 5 pixels
        # NB: 1 pixel = 23m*70m = 1610m2 max => 6.2 pixels min for 1ha
        idx_min_size = np.where(count_labels >= 5)[0]
        logger.debug("> Among them, %d labels contain at least 5 pixels" % idx_min_size.size)
        
        # 4 - Compute the remaining labels_inside to be processed
        if self.product_type == "TILE":
            self.labels_inside = np.intersect1d(self.labels_inside, remaining_labels[idx_min_size])
        elif self.product_type == "SP":
            self.labels_to_process = np.intersect1d(np.unique(self.labels), remaining_labels[idx_min_size])
    
    def compute_edge_pixels(self):
        """
        Compute edge pixels indices and their associated label
        The label previously attached is replaced by a decreasing counter from '9'*NB_DIGITS
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        if (self.nb_obj_at_top_edge + self.nb_obj_at_bottom_edge + self.nb_obj_at_both_edges) == 0:
            logger.debug("NO edge pixel to deal with")
            
        else:
            
            min_xtrack = self.cfg.getfloat('CONFIG_PARAMS', 'MIN_XTRACK')
            max_xtrack = self.cfg.getfloat('CONFIG_PARAMS', 'MAX_XTRACK')
            nb_digits = self.cfg.getint('CONFIG_PARAMS', 'NB_DIGITS')  # Retrieve nb_digits from config file
            cpt_label = int('9'*nb_digits)
            
            # 0 - Masks used below => computed only once
            # 0.1 - Azimuth lines
            first_line_plus2 = self.pixc_first_line_index + 2
            last_line_minus2 = self.pixc_last_line_index - 2
            mask_bottom_line = self.azimuth_index < first_line_plus2
            mask_top_line = self.azimuth_index > last_line_minus2
            # 0.2 - Cross-track distance
            abs_xtrack = np.abs(self.cross_track)
            abs_xtrack_r = np.abs(self.cross_track_refined)
            mask_xtrack = ((abs_xtrack >= min_xtrack) & (abs_xtrack <= max_xtrack) &
                           (abs_xtrack_r >= min_xtrack) & (abs_xtrack_r <= max_xtrack))
            # 0.3 - Practical functions
            edge_idx_lst, edge_lbl_lst, edge_loc_lst = [], [], []
            def _indices(label: int, line_mask: np.ndarray) -> np.ndarray:
                """Combined mask label + azimuth line + cross‑track -> indices."""
                label_mask = (self.labels == label)
                return np.flatnonzero(label_mask & line_mask & mask_xtrack)
            def _store(indices: np.ndarray, loc_code: int) -> None:
                nonlocal cpt_label
                if indices.size == 0:
                    return
                edge_idx_lst.append(indices)
                edge_lbl_lst.append(np.full(indices.size, cpt_label, dtype=int))
                edge_loc_lst.append(np.full(indices.size, loc_code, dtype=int))
                cpt_label -= 1
            
            # 1 - Fill with bottom edge objects
            for indl in np.arange(self.nb_obj_at_bottom_edge):  # Loop over the bottom edge objects
                cur_label = self.labels_at_bottom_edge[indl]
                if cur_label in self.labels_big_lakes:
                    idx = _indices(cur_label, mask_bottom_line)
                    _store(idx, loc_code=10)  # 10 = code for big lakes at bottom edge
                else:
                    idx = np.where(self.labels == cur_label)[0]
                    _store(idx, loc_code=0)  # 0 = code for classic bottom edge objects
                    
            # 2 - Fill with top edge objects
            for indl in np.arange(self.nb_obj_at_top_edge):  # Loop over the top edge objects
                cur_label = self.labels_at_top_edge[indl]
                if cur_label in self.labels_big_lakes:
                    idx = _indices(cur_label, mask_top_line)
                    _store(idx, loc_code=11)  # 11 = code for big lakes at top edge
                else:
                    idx = np.where(self.labels == cur_label)[0]
                    _store(idx, loc_code=1)  # 0 = code for classic top edge objects
                    
            # 3 - Fill with bottom and top edges objects
            for indl in np.arange(self.nb_obj_at_both_edges):  # Loop over both edges objects
                cur_label = self.labels_at_both_edges[indl]
                if cur_label in self.labels_big_lakes:
                    # Bottom edge
                    idx_bot = _indices(cur_label, mask_bottom_line)
                    _store(idx_bot, loc_code=10)
                    # Top edge
                    idx_top = _indices(cur_label, mask_top_line)
                    _store(idx_top, loc_code=11)
                else:
                    idx = np.where(self.labels == cur_label)[0]
                    _store(idx, loc_code=2)  # 0 = code for classic both edges objects
                    
            # 4 - Concatenate lists
            self.edge_index = (np.concatenate(edge_idx_lst)
                                if edge_idx_lst else np.empty(0, dtype=int))
            self.edge_label = (np.concatenate(edge_lbl_lst)
                                if edge_lbl_lst else np.empty(0, dtype=int))
            self.edge_loc   = (np.concatenate(edge_loc_lst)
                                if edge_loc_lst else np.empty(0, dtype=int))
                    
            # 5 - Number of edge pixels
            self.nb_edge_pix = self.edge_index.size

    # ----------------------------------------
    
    def compute_azimuth_partial_flag(self, in_pixc_index):
        """
        Considering the pixels related to the water feature given by in_pixc_index,
        compute partial flag in the azimuth direction, depending:
        - for Tile: if the feature corresponds to a start or a stop of a HR segment
        - for SP: if the feature crosses the azimuth edge of the 1st or last tile of the single-pass band
        
        :param in_pixc_index: indices of pixels defining the studied feature
        :type in_pixc_index: 1D-array of int
        
        :return: out_flag_partial =1 if the studied feature crosses top or bottom of the swath edges; =0 otherwise
        :rtype: out_flag_partial = int
        """
        
        # Init output variable
        out_flag_partial = 0
                
        # Get pixels that belong to the first or last azimuth line
        az_edge_pix = np.where(self.is_boundary_pix[in_pixc_index])[0]
        if len(az_edge_pix) > 0:
            out_flag_partial = 1  
                
        return out_flag_partial

    # ----------------------------------------
    
    def compute_sig0_with_uncertainties(self, in_pixc_index, method='rare'):
        """
        Caller of JPL aggregate.py/sig0_with_uncerts function 
        which computes the aggregation of PIXC sigma0 over a feature, with corresponding uncertainty
        Compute height only over pixels having self.sig0 != numpy.nan
    
        :param in_pixc_index: indices of pixels to consider for computation
        :type in_pixc_index: 1D-array of int
        :param method: type of aggregator ('rare' or 'medium')
        :type method: string
        
        :return: out_sig0 = aggregated sigma0
        :rtype: out_sig0 = float
        :return: out_sig0_std = standard deviation of the sigma0
        :rtype: out_sig0_std = float
        :return: out_sig0_unc = sigma0 uncertainty for the feature
        :rtype: out_sig0_unc = float
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Get pixels having self.sig0 != numpy_nan

        not_nan_idx = np.argwhere(np.isfinite(self.sig0[in_pixc_index]))
        nb_not_nan = len(not_nan_idx)
        
        # Process depending their number
        if len(not_nan_idx) == 0:
            logger.debug("WARNING All pixels have sig0=NaN => output is set to NaN")
            out_sig0 = np.nan
            out_sig0_std = np.nan
            out_sig0_unc = np.nan

        else:            
            nb_pixc = len(in_pixc_index)
            
            if nb_not_nan == nb_pixc:
                # log logger.debug("There are NO pixels having sig0=NaN => aggregate sig0 over all pixels")
                # Call JPL function
                out_sig0, out_sig0_std, out_sig0_unc = jpl_aggregate.sig0_with_uncerts(
                    self.sig0[in_pixc_index], 
                    self.sig0_u[in_pixc_index],
                    method=method)
                
            else:
                logger.debug("WARNING There are %d/%d pixels having sig0=NaN => aggregate sig0 only over the others" 
                               % (nb_pixc-nb_not_nan, nb_pixc))
                selected_pixc_idx = [in_pixc_index[not_nan_idx]]
                # Call JPL function
                out_sig0, out_sig0_std, out_sig0_unc = jpl_aggregate.sig0_with_uncerts(
                    self.sig0[selected_pixc_idx], 
                    self.sig0_u[selected_pixc_idx], 
                    method=method)
        
        return out_sig0, out_sig0_std, out_sig0_unc
    
    def compute_height(self, in_pixc_index, method='weight'):
        """
        Caller of JPL aggregate.py/height_only function 
        which computes the aggregation of PIXC height over a feature
        Compute height only over pixels having self.height != numpy.nan
        
        !!! Used in input of height-constrained geolocation => MUST BE APPLIED 
            ON self.height AND NOT self.corrected_height because use of heights
            above the ellipsoid
    
        :param in_pixc_index: indices of pixels to consider for computation
        :type in_pixc_index: 1D-array of int
        :param method: type of aggregator ('weight', 'uniform', or 'median')
        :type method: string
        
        :return: out_height = aggregated height
        :rtype: out_height = float
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Get pixels having self.height != numpy_nan
        not_nan_idx = np.argwhere(np.isfinite(self.height[in_pixc_index]))
        nb_not_nan = len(not_nan_idx)
        
        # Process depending their number
        if len(not_nan_idx) == 0:
            logger.debug("WARNING All pixels have height=NaN => output height is set to NaN")
            out_height = np.nan
            
        else:            
            nb_pixc = len(in_pixc_index)
            
            if nb_not_nan == nb_pixc:
                # log logger.debug("There are NO pixels having height=NaN => aggregate height over all pixels")
                # Call JPL function
                out_height, weight_norm = jpl_aggregate.height_only(self.height[in_pixc_index], 
                                                                    height_std=self.height_std_pix[in_pixc_index],
                                                                    method=method)
                
            else:
                logger.debug("WARNING There are %d/%d pixels having height=NaN => aggregate height only over the others" % (nb_pixc-nb_not_nan, nb_pixc))
                selected_pixc_idx = [in_pixc_index[not_nan_idx]]
                # Call JPL function
                out_height, weight_norm = jpl_aggregate.height_only(self.height[selected_pixc_idx], 
                                                                    height_std=self.height_std_pix[selected_pixc_idx],
                                                                    method=method)
        
        return out_height
    
    def compute_height_with_uncertainties(self, in_pixc_index, method='weight', delta_mls=None):
        """
        Caller of JPL aggregate.py/height_with_uncerts function 
        which computes the aggregation of PIXC height over a feature, with corresponding uncertainty
        Compute height only over pixels having self.corrected_height != numpy.nan
    
        :param in_pixc_index: indices of pixels to consider for computation
        :type in_pixc_index: 1D-array of int
        :param method: type of aggregator ('weight', 'uniform', or 'median')
        :type method: string
        :param delta_mls: difference between SWOT-observed mean lake surface (MLS) and the EGM2008 geoid model (self.geoid)
        :type delta_mls: 1D-array of float
        
        :return: out_height = aggregated height
        :rtype: out_height = float
        :return: out_height_std = standard deviation of the heights
        :rtype: out_height_std = float
        :return: out_height_unc = height uncertainty for the feature
        :rtype: out_height_unc = float
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Applied delta_mls if specified
        tmp_corrected_height = np.copy(self.corrected_height[in_pixc_index])
        if delta_mls is not None:
            tmp_corrected_height -= delta_mls
        
        # Get pixels having self.corrected_height != numpy_nan

        not_nan_idx = np.argwhere(np.isfinite(tmp_corrected_height))
        nb_not_nan = len(not_nan_idx)
        
        # Process depending their number
        if len(not_nan_idx) == 0:
            logger.debug("WARNING All pixels have corrected_height=NaN => output height is set to NaN")
            out_height = np.nan
            out_height_std = np.nan
            out_height_unc = np.nan

        else:            
            nb_pixc = len(in_pixc_index)
            
            if nb_not_nan == nb_pixc:
                # log logger.debug("There are NO pixels having corrected_height=NaN => aggregate height over all pixels")
                # Call JPL function
                out_height, out_height_std, out_height_unc, lat_unc, lon_unc = jpl_aggregate.height_with_uncerts(
                    tmp_corrected_height, 
                    self.eff_num_rare_looks[in_pixc_index], 
                    self.eff_num_medium_looks[in_pixc_index],
                    self.interferogram_flattened[in_pixc_index], 
                    self.power_plus_y[in_pixc_index], 
                    self.power_minus_y[in_pixc_index], 
                    self.looks_to_efflooks,
                    self.dheight_dphase[in_pixc_index], 
                    self.dlatitude_dphase[in_pixc_index], 
                    self.dlongitude_dphase[in_pixc_index], 
                    height_std=self.height_std_pix[in_pixc_index],
                    method=method)
                
            else:
                logger.debug("WARNING There are %d/%d pixels having corrected_height=NaN => aggregate height only over the others" 
                               % (nb_pixc-nb_not_nan, nb_pixc))
                selected_pixc_idx = [in_pixc_index[not_nan_idx]]
                # Call JPL function
                out_height, out_height_std, out_height_unc, lat_unc, lon_unc = jpl_aggregate.height_with_uncerts(
                    tmp_corrected_height[not_nan_idx], 
                    self.eff_num_rare_looks[selected_pixc_idx], 
                    self.eff_num_medium_looks[selected_pixc_idx],
                    self.interferogram_flattened[selected_pixc_idx], 
                    self.power_plus_y[selected_pixc_idx], 
                    self.power_minus_y[selected_pixc_idx], 
                    self.looks_to_efflooks,
                    self.dheight_dphase[selected_pixc_idx], 
                    self.dlatitude_dphase[selected_pixc_idx], 
                    self.dlongitude_dphase[selected_pixc_idx], 
                    height_std=self.height_std_pix[selected_pixc_idx],
                    method=method)
        
        return out_height, out_height_std, out_height_unc
    
    def compute_lakeshore_erosion(self, in_pixc_index):
        """
        Erodes the contours of water pixels in radar geometry so that only half the edge pixels are considered.

        :param in_pixc_index: indices of pixels to consider for computation
        :type in_pixc_index: 1D-array of int

        :return: :param in_pixc_labels_mask_erode: indices of pixels after erosion
        :type in_pixc_labels_mask_erode: 1D-array of int
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        idx_mandatory = self.classif_dict["mandatory"][in_pixc_index]
        size_idx_mandatory = np.sum(idx_mandatory)
        idx_open_water = self.classif_dict["open_water"][in_pixc_index]
        size_idx_open_water = np.sum(idx_open_water)

        if (size_idx_open_water != 0) & (size_idx_mandatory > size_idx_open_water) :
            idx_for_edge_pixels = self.classif_dict["water_edge"][in_pixc_index]

            # Convert feature in a 2D matrix
            in_x = self.range_index[in_pixc_index]
            in_x -= int(np.min(in_x))
            nb_pix_range = int(np.max(in_x) + 1)
            in_y = self.azimuth_index[in_pixc_index]
            in_y -= int(np.min(in_y))
            nb_pix_azimuth = int(np.max(in_y) + 1)
            in_pixc_mask = my_tools.compute_bin_mat(nb_pix_range, nb_pix_azimuth, in_x, in_y)
            
            # Erosion with kernel vertical to remove only one pixel in azimuth
            in_pixc_mask_erode = binary_erosion(in_pixc_mask, structure=np.array([[1], [1]])).astype(in_pixc_mask.dtype)
            in_pixc_mask_erode_1D = my_tools.convert_2d_mat_in_1d_vec(in_x, in_y, in_pixc_mask_erode)

            # Keep erosion only if the pixel is of classif 3 or 6
            in_pixc_labels_mask_erode = [int(idx) for i, idx in enumerate(in_pixc_index) if not (in_pixc_mask_erode_1D[i] == 0 and idx_for_edge_pixels[i] == 1)]
            nb_erode = len(in_pixc_index) - len(in_pixc_labels_mask_erode)
            logger.debug("%d pixels are rejected after lakeshore erosion" % nb_erode)

            return in_pixc_labels_mask_erode
        
        logger.debug("WARNING NO water edge pixel to compute lakeshore erosion => all pixels kept")
        return in_pixc_index

    def compute_area(self, in_pixc_index, method='composite', flag_all=True):
        """
        Caller of JPL aggregate.py/area_only function 
        which computes the aggregation of PIXC area over a feature
        Compute area only over pixels having self.pixel_area != numpy.nan
    
        :param in_pixc_index: indices of pixels to consider for computation
        :type in_pixc_index: 1D-array of int
        :param method: type of aggregator ('weight', 'uniform', or 'median')
        :type method: string
        
        :return: out_area = total water area (=water + dark water pixels)
        :rtype: out_area = float
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        if len(in_pixc_index) == 0:
            logger.debug("WARNING NO pixel in input => output area is set to 0")
            out_area = 0.0
            
        else:
        
            # 1 - Select only pixels having self.pixel_area != numpy.nan
            not_nan_idx = np.argwhere(np.isfinite(self.pixel_area[in_pixc_index]))
            nb_not_nan = len(not_nan_idx)
            
            if len(not_nan_idx) == 0:
                logger.debug("WARNING All pixels have pixel_area=NaN => output area is set to NaN")
                out_area = np.nan
                
            else:            
                nb_pixc = len(in_pixc_index)
                
                if nb_not_nan == nb_pixc:
                    # log logger.debug("There are NO pixels having pixel_area=NaN => aggregate area over all pixels")
                    selected_pixc_idx = in_pixc_index
                else:
                    logger.debug("WARNING There are %d/%d pixels having pixel_area=NaN => aggregate area only over the others"
                                   % (nb_pixc-nb_not_nan, nb_pixc))
                    selected_pixc_idx = in_pixc_index[not_nan_idx]
    
                # Compute classif dict only for the selected pixels
                selected_classif = dict()
                for key in self.classif_dict.keys():
                    selected_classif[key] = self.classif_dict[key][selected_pixc_idx]
    
                # 2 - Aggregate PIXC area and uncertainty considering all PIXC (ie water + dark water)
                # Call external function common with RiverObs
                out_area, nb_pixc = jpl_aggregate.area_only(
                        self.pixel_area[selected_pixc_idx], self.water_frac[selected_pixc_idx], 
                        selected_classif, method=method)
                    
            # 4 - Handle specific case when out_area and/or out_area_detct is negative
            # due to amount of negative water_frac leading to negative inundated_area
            if (out_area < 0):
                logger.debug("WARNING negative area => output area is set to NaN")
                out_area = np.nan
        
        return out_area

    def compute_area_with_uncertainties(self, in_pixc_index, method='composite', flag_all=True):
        """
        Caller of JPL aggregate.py/area_with_uncert function 
        which computes the aggregation of PIXC area over a feature, with corresponding uncertainty
        Compute area only over pixels having self.pixel_area != numpy.nan
    
        :param in_pixc_index: indices of pixels to consider for computation
        :type in_pixc_index: 1D-array of int
        :param method: type of aggregator ('weight', 'uniform', or 'median')
        :type method: string
        :param flag_all: if True (default), compute all areas and uncertainties; else, compute only total water area and uncertainty
        :type flag_all: boolean
        
        :return: out_area = total water area (=water + dark water pixels)
        :rtype: out_area = float
        :return: out_area_unc = uncertainty in total water area
        :rtype: out_area_unc = float
        :return: out_area_detct = area of detected water pixels
        :rtype: out_area_detct = float
        :return: out_area_detct_unc = uncertainty in area of detected water
        :rtype: out_area_detct_unc = float
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        if len(in_pixc_index) == 0:
            logger.debug("WARNING NO pixel in input => output area is set to 0")
            out_area = 0.0
            out_area_unc = np.nan
            out_area_detct = np.nan
            out_area_detct_unc = np.nan
            
        else:
        
            # 1 - Select only pixels having self.pixel_area != numpy.nan
            not_nan_idx = np.argwhere(np.isfinite(self.pixel_area[in_pixc_index]))
            nb_not_nan = len(not_nan_idx)
            
            if len(not_nan_idx) == 0:
                logger.debug("WARNING All pixels have pixel_area=NaN => output area are set to NaN")
                out_area = np.nan
                out_area_unc = np.nan
                out_area_detct = np.nan
                out_area_detct_unc = np.nan
                
            else:            
                nb_pixc = len(in_pixc_index)
                
                if nb_not_nan == nb_pixc:
                    # log logger.debug("There are NO pixels having pixel_area=NaN => aggregate area over all pixels")
                    selected_pixc_idx = in_pixc_index
                else:
                    logger.debug("WARNING There are %d/%d pixels having pixel_area=NaN => aggregate area only over the others"
                                   % (nb_pixc-nb_not_nan, nb_pixc))
                    selected_pixc_idx = in_pixc_index[not_nan_idx]
                # Remove 1 pixel in top/bottom azimuth to compute lake area
                if self.cfg.getboolean("CONFIG_PARAMS", "DEL_1PIX_IN_AZ_4HULL"):
                    selected_pixc_idx = self.compute_lakeshore_erosion(selected_pixc_idx)
    
                # Compute classif dict only for the selected pixels
                selected_classif = dict()
                for key in self.classif_dict.keys():
                    selected_classif[key] = self.classif_dict[key][selected_pixc_idx]
    
                # 2 - Aggregate PIXC area and uncertainty considering all PIXC (ie water + dark water)
                # Call external function common with RiverObs
                out_area, out_area_unc, area_pcnt_uncert = jpl_aggregate.area_with_uncert(
                        self.pixel_area[selected_pixc_idx], self.water_frac[selected_pixc_idx], self.water_frac_uncert[selected_pixc_idx],
                        self.darea_dheight[selected_pixc_idx], selected_classif, self.false_detection_rate[selected_pixc_idx],
                        self.missed_detection_rate[selected_pixc_idx],
                        method=method)
                
                if flag_all:
                
                    # 3 - Aggregate PIXC area and uncertainty considering only water PIXC (ie detected water without dark water)
                    # Call external function common with RiverObs
                    out_area_detct, out_area_detct_unc, area_detct_pcnt_uncert = jpl_aggregate.area_with_uncert(
                            self.pixel_area[selected_pixc_idx], self.water_frac[selected_pixc_idx], self.water_frac_uncert[selected_pixc_idx],
                            self.darea_dheight[selected_pixc_idx], selected_classif, self.false_detection_rate[selected_pixc_idx],
                            self.missed_detection_rate[selected_pixc_idx],
                            dark_water_klasses=False,
                            method=method)
                    
                else:
                    out_area_detct = np.nan
                    out_area_detct_unc = np.nan
                    
            # 4 - Handle specific case when out_area and/or out_area_detct is negative
            # due to amount of negative water_frac leading to negative inundated_area
            if (out_area < 0) or (out_area_detct < 0):
                logger.debug("WARNING negative area => output area are set to NaN")
                out_area = np.nan
                out_area_unc = np.nan
                out_area_detct = np.nan
                out_area_detct_unc = np.nan
        
        return out_area, out_area_unc, out_area_detct, out_area_detct_unc
    
    def compute_geophysical_ref(self, in_name, in_pixc_index, method='weight', delta_mls=None):
        """
        Caller of JPL aggregate.py/height_only function 
        which computes the aggregation of PIXC geophysical reference value over a feature
        Compute geophysical ref value only over pixels having self.<geophysical_ref> != numpy.nan
    
        :param in_name: name of the geophysical ref to process ('geoid_hght', 'delta_mls', 'solid_tide', 'pole_tide', 'load_tidef', or 'load_tideg')
        :type in_name: string
        :param in_pixc_index: indices of pixels to consider for computation
        :type in_pixc_index: 1D-array of int
        :param method: type of aggregator ('weight', 'uniform', or 'median')
        :type method: string
        :param delta_mls: difference between SWOT-observed mean lake surface and the EGM2008 geoid model; same number of elements as in_pixc_index
        :type delta_mls: 1D-array of float
        
        :return: out_value = aggregated geophysical value
        :rtype: out_value = float
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Retrieve geophysical ref array to process
        geophysical_ref = None
        if in_name == 'geoid_hght':
            geophysical_ref = self.geoid[in_pixc_index]
        elif in_name == 'delta_mls':
            if delta_mls is None:
                logger.debug("MLS is set to None => return delta_mls = NaN")
                return np.nan
            geophysical_ref = delta_mls
        elif in_name == 'solid_tide':
            geophysical_ref = self.solid_earth_tide[in_pixc_index]
        elif in_name == 'pole_tide':
            geophysical_ref = self.pole_tide[in_pixc_index]
        elif in_name == 'load_tidef':
            geophysical_ref = self.load_tide_fes[in_pixc_index]
        elif in_name == 'load_tideg':
            geophysical_ref = self.load_tide_got[in_pixc_index]
        else:
            logger.error("in_name = %d => should be among geoid_hght, delta_mls, solid_tide, pole_tide, load_tidef, or load_tideg" % in_name)
            raise
        
        # Get pixels having geophysical_ref != numpy_nan
        not_nan_idx = np.argwhere(np.isfinite(geophysical_ref))
        nb_not_nan = len(not_nan_idx)
        
        # Process depending their number
        if len(not_nan_idx) == 0:
            logger.debug("WARNING All pixels have %s=NaN => output value is set to NaN" % in_name)
            out_value = np.nan
            
        else:            
            nb_pixc = len(in_pixc_index)
            
            if nb_not_nan == nb_pixc:
                #logger.debug("There are NO pixels having %s=NaN => aggregate over all pixels" % in_name)
                # Call JPL function
                out_value, weight_norm = jpl_aggregate.height_only(geophysical_ref, 
                                                                   height_std=self.height_std_pix[in_pixc_index], 
                                                                   method=method)
                
            else:
                logger.debug("WARNING There are %d/%d pixels having %s=NaN => aggregate only over the others" % (nb_pixc-nb_not_nan, nb_pixc, in_name))
                selected_pixc_idx = in_pixc_index[not_nan_idx]
                # Call JPL function
                out_value, weight_norm = jpl_aggregate.height_only(geophysical_ref[not_nan_idx], 
                                                                   height_std=self.height_std_pix[selected_pixc_idx],
                                                                   method=method)
        
        return out_value
    
    def compute_interferogram_flattened(self, in_pixc_index, in_p_final):
        """
        Flatten the interferogram for the feature definied by the PIXC of indices in_pixc_index
        
        :param in_pixc_index: indices of pixels defining the studied feature
        :type in_pixc_index: 1D-array of int
        :param in_p_final: position of pixels of the feature in geocentric coordinates (same length as in_pixc_index)
        :type in_p_final: numpy 2D-array of float; size (3=x/y/z, nb_pixc)
        """
        
        # 1 - Format variables
        # 1.1 - Cartesian coordinates of plus_y antenna
        plus_y_antenna_xyz = (self.nadir_plus_y_antenna_x , self.nadir_plus_y_antenna_y , self.nadir_plus_y_antenna_z)
        # 1.2 - Cartesian coordinates of minus_y antenna
        minus_y_antenna_xyz = (self.nadir_minus_y_antenna_x , self.nadir_minus_y_antenna_y, self.nadir_minus_y_antenna_z)
        # 1.3 - Position of PIXC in geocentric coordinates
        target_xyz = (in_p_final[:,0], in_p_final[:,1], in_p_final[:,2])
        
        # 2 - Call to external function, common with RiverObs
        self.interferogram_flattened[in_pixc_index] = jpl_aggregate.flatten_interferogram(
                self.interferogram[in_pixc_index], 
                plus_y_antenna_xyz, minus_y_antenna_xyz, 
                target_xyz, in_pixc_index, self.wavelength)          

    # ---------------------------------------- 
    
    def write_obj_pixc_pixcvec_as_shp(self, in_filename, in_obj_pixcvec, flag_vectorproc=True):
        """
        Write PIXC and PIXCVec variables of interest in a shapefile

        :param in_filename: full path of the output file
        :type in_filename: string
        :param in_obj_pixcvec: PIXCVec object related to self PIXC object
        :type in_obj_pixcvec: proc_pixc_vec.PixelCloudVec
        :param flag_vectorproc: =True to use height-constrained geolocation as points coordinates; =False to use PIXC original coordinates
        :type flag_vectorproc: boolean
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Output PIXC_PIXCVec shapefile = %s" % in_filename)

        # 1 - Initialisation du fichier de sortie
        # 1.1 - Driver
        shp_driver = ogr.GetDriverByName(str("ESRI Shapefile"))
        # 1.2 - Create file
        if os.path.exists(in_filename):
            shp_driver.DeleteDataSource(in_filename)
        out_data_source = shp_driver.CreateDataSource(in_filename)
        # 1.3 - Create layer
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)  # WGS84
        out_layer = out_data_source.CreateLayer(str(str(os.path.basename(in_filename)).replace('.shp', '')), srs, geom_type=ogr.wkbPoint)
        # 1.4 - Create attributes
        out_layer.CreateField(ogr.FieldDefn(str('az_index'), ogr.OFTInteger))
        out_layer.CreateField(ogr.FieldDefn(str('r_index'), ogr.OFTInteger))
        out_layer.CreateField(ogr.FieldDefn(str('classif'), ogr.OFTInteger))
        out_layer.CreateField(ogr.FieldDefn(str('classif_q'), ogr.OFTInteger))
        out_layer.CreateField(ogr.FieldDefn(str('geoloc_q'), ogr.OFTInteger))
        tmp_field = ogr.FieldDefn(str('xtrack'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(4)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('xtrack2'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(4)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('dif_xtrack'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(4)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('sig0'), ogr.OFTReal)
        tmp_field.SetWidth(20)
        tmp_field.SetPrecision(6)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('height'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(4)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('height_cor'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(4)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('height2'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(4)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('dif_height'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(4)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('pixel_area'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(6)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('water_frac'), ogr.OFTReal)
        tmp_field.SetWidth(20)
        tmp_field.SetPrecision(6)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('inund_area'), ogr.OFTReal)
        tmp_field.SetWidth(20)
        tmp_field.SetPrecision(6)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('xover'), ogr.OFTReal)
        tmp_field.SetWidth(10)
        tmp_field.SetPrecision(6)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('geoid'), ogr.OFTReal)
        tmp_field.SetWidth(10)
        tmp_field.SetPrecision(6)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('layover'), ogr.OFTReal)
        tmp_field.SetWidth(10)
        tmp_field.SetPrecision(6)
        out_layer.CreateField(tmp_field)
        out_layer.CreateField(ogr.FieldDefn(str('reach_id'), ogr.OFTString))
        out_layer.CreateField(ogr.FieldDefn(str('lake_id'), ogr.OFTString))
        out_layer.CreateField(ogr.FieldDefn(str('obs_id'), ogr.OFTString))
        out_layer_defn = out_layer.GetLayerDefn()

        # 2 - On traite point par point
        for ind_pixc, ind_pixcvec in zip(range(self.nb_selected), self.selected_index):
            # 2.1 - On cree l'objet dans le format de la couche de sortie
            out_feature = ogr.Feature(out_layer_defn)
            # 2.2 - On lui assigne le point
            point = ogr.Geometry(ogr.wkbPoint)
            if flag_vectorproc:
                point.AddPoint(in_obj_pixcvec.longitude_vectorproc[ind_pixcvec], in_obj_pixcvec.latitude_vectorproc[ind_pixcvec])
            else:
                point.AddPoint(self.longitude[ind_pixc], self.latitude[ind_pixc])
            out_feature.SetGeometry(point)
            # 2.3 - On lui assigne les attributs
            out_feature.SetField(str('az_index'), int(self.azimuth_index[ind_pixc]))
            out_feature.SetField(str('r_index'), int(self.range_index[ind_pixc]))
            out_feature.SetField(str('classif'), int(self.classif[ind_pixc]))
            out_feature.SetField(str('classif_q'), int(self.classification_qual[ind_pixc]))
            out_feature.SetField(str('geoloc_q'), int(self.geolocation_qual[ind_pixc]))
            out_feature.SetField(str('xtrack'), float(self.cross_track[ind_pixc]))
            out_feature.SetField(str('xtrack2'), float(self.cross_track_refined[ind_pixc]))
            out_feature.SetField(str('dif_xtrack'), float(self.cross_track_refined[ind_pixc])-float(self.cross_track[ind_pixc]))
            out_feature.SetField(str('sig0'), float(self.sig0[ind_pixc]))
            out_feature.SetField(str('height'), float(self.height[ind_pixc]))
            out_feature.SetField(str('height_cor'), float(self.corrected_height[ind_pixc]))
            tmp_height2 = float(in_obj_pixcvec.height_vectorproc[ind_pixcvec])
            if tmp_height2 > -1e6:
                out_feature.SetField(str('height2'), tmp_height2)
                out_feature.SetField(str('dif_height'), tmp_height2-float(self.height[ind_pixc]))
            out_feature.SetField(str('pixel_area'), float(self.pixel_area[ind_pixc]))
            out_feature.SetField(str('water_frac'), float(self.water_frac[ind_pixc]))
            if float(self.inundated_area[ind_pixc]) > -1e6:
                out_feature.SetField(str('inund_area'), float(self.inundated_area[ind_pixc]))
            out_feature.SetField(str('xover'), float(self.height_cor_xover[ind_pixc]))
            out_feature.SetField(str('geoid'), float(self.geoid[ind_pixc]))
            out_feature.SetField(str('layover'), float(self.layover_impact[ind_pixc]))
            if in_obj_pixcvec.reach_id[ind_pixcvec] != b'': 
                out_feature.SetField(str('reach_id'), in_obj_pixcvec.reach_id[ind_pixcvec].decode('utf-8'))
            if in_obj_pixcvec.lake_id[ind_pixcvec] != b'': 
                out_feature.SetField(str('lake_id'), in_obj_pixcvec.lake_id[ind_pixcvec].decode('utf-8'))
            out_feature.SetField(str('obs_id'), in_obj_pixcvec.obs_id[ind_pixcvec].decode('utf-8'))
            # 2.4 - On ajoute l'objet dans la couche de sortie
            out_layer.CreateFeature(out_feature)

        # 3 - Destroy the data sources to free resources
        out_data_source.Destroy()

    # ----------------------------------------

    def write_edge_file(self, in_filename, in_proc_metadata):
        """
        Save pixels related to objects at top/bottom edge of the PixC tile in a NetCDF file.
        These 1D-arrays indicate pixel index, associated label, location (top/bottom/both edges) and needed PixC variables.
        If there is no pixel at tile edge, the file is empty but exists.
        
        :param in_filename: full path of the output file
        :type in_filename: string
        :param in_proc_metadata: processing metadata
        :type in_proc_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Output L2_HR_LakeTile_Edge NetCDF file = %s" % in_filename)
        
        # 0.1 - Save pixc_metadata dict
        tmp_pixc_metadata = self.pixc_metadata.copy()
        tmp_pixc_metadata.pop("time_coverage_start")
        tmp_pixc_metadata.pop("time_coverage_end")
        # 0.2 - Estimate time_coverage_start and time_coverage_end
        time_str_dict = dict()
        if self.nb_edge_pix != 0:
            nadir_time_edge_pixels = self.nadir_time[self.edge_index]
            time_str_dict["time_coverage_start"] = my_tools.convert_utc_to_str(min(nadir_time_edge_pixels), in_format=2)
            time_str_dict["time_coverage_end"] = my_tools.convert_utc_to_str(max(nadir_time_edge_pixels), in_format=2)
        else:
            time_str_dict["time_coverage_start"] = "None"
            time_str_dict["time_coverage_end"] = "None"
        # 0.3 - Add self.pixc_first|last_line_index 
        metadata_azimuth_index_for_tile = dict()
        metadata_azimuth_index_for_tile["pixc_first_line_index"] = self.pixc_first_line_index
        metadata_azimuth_index_for_tile["pixc_last_line_index"] = self.pixc_last_line_index
        metadata_azimuth_index_for_tile["pixc_first_line_is_start_hr"] = int(self.seg_hr_start)
        metadata_azimuth_index_for_tile["pixc_last_line_is_stop_hr"] = int(self.seg_hr_stop)
        # 0.4 - Add list of big lakes processed in the tile
        metadata_big_lakes = dict()
        metadata_big_lakes["big_lakeid"] = ";".join(self.big_lakeid)
        
        # 1 - Init product
        edge_file = nc_file.LakeTileEdgeProduct(in_pixc_metadata=tmp_pixc_metadata, 
                                                in_proc_metadata={**in_proc_metadata, **time_str_dict, **metadata_azimuth_index_for_tile, **metadata_big_lakes})

        # 2 - Form dictionary with variables to write
        vars_to_write = {}
        if self.nb_selected != 0:
            vars_to_write["edge_index"] = self.selected_index[self.edge_index]
            vars_to_write["edge_label"] = self.edge_label
            vars_to_write["edge_loc"] = self.edge_loc
            vars_to_write["pixel_at_window_edge"] = self.pixel_at_window_edge[self.edge_index]
            vars_to_write["ice_dyn_f"] = self.ice_dyn_f[self.edge_index]
            vars_to_write["classification"] = self.classif[self.edge_index]
            vars_to_write["range_index"] = self.range_index[self.edge_index]
            vars_to_write["azimuth_index"] = self.azimuth_index[self.edge_index]
            vars_to_write["interferogram"] = np.stack((np.real(self.interferogram[self.edge_index]), np.imag(self.interferogram[self.edge_index]))).T
            vars_to_write["power_plus_y"] = self.power_plus_y[self.edge_index]
            vars_to_write["power_minus_y"] = self.power_minus_y[self.edge_index]
            vars_to_write["water_frac"] = self.water_frac[self.edge_index]
            vars_to_write["water_frac_uncert"] = self.water_frac_uncert[self.edge_index]
            vars_to_write["false_detection_rate"] = self.false_detection_rate[self.edge_index]
            vars_to_write["missed_detection_rate"] = self.missed_detection_rate[self.edge_index]
            vars_to_write["bright_land_flag"] = self.bright_land_flag[self.edge_index]
            vars_to_write["layover_impact"] = self.layover_impact[self.edge_index]
            vars_to_write["eff_num_rare_looks"] = self.eff_num_rare_looks[self.edge_index]
            vars_to_write["latitude"] = self.latitude[self.edge_index]
            vars_to_write["longitude"] = self.longitude[self.edge_index]
            vars_to_write["height"] = self.height[self.edge_index]
            vars_to_write["cross_track"] = self.cross_track[self.edge_index]
            vars_to_write["cross_track_refined"] = self.cross_track_refined[self.edge_index]
            vars_to_write["pixel_area"] = self.pixel_area[self.edge_index]
            vars_to_write["phase_noise_std"] = self.phase_noise_std[self.edge_index]
            vars_to_write["dlatitude_dphase"] = self.dlatitude_dphase[self.edge_index]
            vars_to_write["dlongitude_dphase"] = self.dlongitude_dphase[self.edge_index]
            vars_to_write["dheight_dphase"] = self.dheight_dphase[self.edge_index]
            vars_to_write["darea_dheight"] = self.darea_dheight[self.edge_index]
            vars_to_write["eff_num_medium_looks"] = self.eff_num_medium_looks[self.edge_index]
            vars_to_write["sig0"] = self.sig0[self.edge_index]
            vars_to_write["sig0_uncert"] = self.sig0_u[self.edge_index]
            vars_to_write["model_dry_tropo_cor"] = self.model_dry_tropo_cor[self.edge_index]
            vars_to_write["model_wet_tropo_cor"] = self.model_wet_tropo_cor[self.edge_index]
            vars_to_write["iono_cor_gim_ka"] = self.iono_cor_gim_ka[self.edge_index]
            vars_to_write["height_cor_xover"] = self.height_cor_xover[self.edge_index]
            vars_to_write["geoid"] = self.geoid[self.edge_index]
            vars_to_write["solid_earth_tide"] = self.solid_earth_tide[self.edge_index]
            vars_to_write["load_tide_fes"] = self.load_tide_fes[self.edge_index]
            vars_to_write["load_tide_got"] = self.load_tide_got[self.edge_index]
            vars_to_write["pole_tide"] = self.pole_tide[self.edge_index]
            vars_to_write["classification_qual"] = self.classification_qual[self.edge_index]
            vars_to_write["geolocation_qual"] = self.geolocation_qual[self.edge_index]
            vars_to_write["time_from_prev_xover"] = self.prev_xover[self.edge_index]
            vars_to_write["time_to_next_xover"] = self.next_xover[self.edge_index]
            vars_to_write["nadir_time"] = self.nadir_time[self.edge_index]
            vars_to_write["nadir_time_tai"] = self.nadir_time_tai[self.edge_index]
            vars_to_write["nadir_longitude"] = self.nadir_longitude[self.edge_index]
            vars_to_write["nadir_latitude"] = self.nadir_latitude[self.edge_index]
            vars_to_write["nadir_x"] = self.nadir_x[self.edge_index]
            vars_to_write["nadir_y"] = self.nadir_y[self.edge_index]
            vars_to_write["nadir_z"] = self.nadir_z[self.edge_index]
            vars_to_write["nadir_vx"] = self.nadir_vx[self.edge_index]
            vars_to_write["nadir_vy"] = self.nadir_vy[self.edge_index]
            vars_to_write["nadir_vz"] = self.nadir_vz[self.edge_index]
            vars_to_write["nadir_plus_y_antenna_x"] = self.nadir_plus_y_antenna_x[self.edge_index]
            vars_to_write["nadir_plus_y_antenna_y"] = self.nadir_plus_y_antenna_y[self.edge_index]
            vars_to_write["nadir_plus_y_antenna_z"] = self.nadir_plus_y_antenna_z[self.edge_index]
            vars_to_write["nadir_minus_y_antenna_x"] = self.nadir_minus_y_antenna_x[self.edge_index]
            vars_to_write["nadir_minus_y_antenna_y"] = self.nadir_minus_y_antenna_y[self.edge_index]
            vars_to_write["nadir_minus_y_antenna_z"] = self.nadir_minus_y_antenna_z[self.edge_index]
            
        # 3 - Write file
        edge_file.write_product(in_filename, self.nb_edge_pix, vars_to_write)

    def write_edge_file_as_shp(self, in_filename):
        """
        Write PixC subset related to edge (top/bottom) objects as a shapefile

        :param in_filename: full path of the output file
        :type in_filename: string
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Output L2_HR_LakeTile_Edge shapefile = %s" % in_filename)

        # 1 - Initialisation du fichier de sortie
        # 1.1 - Driver
        shp_driver = ogr.GetDriverByName(str("ESRI Shapefile"))
        # 1.2 - Create file
        if os.path.exists(in_filename):
            shp_driver.DeleteDataSource(in_filename)
        out_data_source = shp_driver.CreateDataSource(in_filename)
        # 1.3 - Create layer
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)  # WGS84
        out_layer = out_data_source.CreateLayer(str(str(os.path.basename(in_filename)).replace('.shp', '')), srs, geom_type=ogr.wkbPoint)
        # 1.4 - Create attributes
        out_layer.CreateField(ogr.FieldDefn(str('edge_index'), ogr.OFTInteger))
        out_layer.CreateField(ogr.FieldDefn(str('edge_label'), ogr.OFTInteger))
        out_layer.CreateField(ogr.FieldDefn(str('edge_loc'), ogr.OFTInteger))
        out_layer.CreateField(ogr.FieldDefn(str('az_index'), ogr.OFTInteger))
        out_layer.CreateField(ogr.FieldDefn(str('r_index'), ogr.OFTInteger))
        out_layer.CreateField(ogr.FieldDefn(str('classif'), ogr.OFTInteger))
        tmp_field = ogr.FieldDefn(str('xtrack'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(4)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('xtrack2'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(4)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('sig0'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(4)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('height'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(4)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('pixel_area'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(6)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('water_frac'), ogr.OFTReal)
        tmp_field.SetWidth(11)
        tmp_field.SetPrecision(6)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('nadir_t'), ogr.OFTReal)
        tmp_field.SetWidth(13)
        tmp_field.SetPrecision(3)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('nadir_long'), ogr.OFTReal)
        tmp_field.SetWidth(11)
        tmp_field.SetPrecision(6)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('nadir_lat'), ogr.OFTReal)
        tmp_field.SetWidth(10)
        tmp_field.SetPrecision(6)
        out_layer.CreateField(tmp_field)
        out_layer_defn = out_layer.GetLayerDefn()

        # 2 - On traite point par point
        for indp in range(self.nb_edge_pix):
            tmp_index = self.edge_index[indp]
            # 2.1 - On cree l'objet dans le format de la couche de sortie
            out_feature = ogr.Feature(out_layer_defn)
            # 2.2 - On lui assigne le point
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(self.longitude[tmp_index], self.latitude[tmp_index])
            out_feature.SetGeometry(point)
            # 2.3 - On lui assigne les attributs
            out_feature.SetField(str('edge_index'), int(self.selected_index[tmp_index]))
            out_feature.SetField(str('edge_label'), int(self.edge_label[indp]))
            out_feature.SetField(str('edge_loc'), int(self.edge_loc[indp]))
            out_feature.SetField(str('az_index'), int(self.azimuth_index[tmp_index]))
            out_feature.SetField(str('r_index'), int(self.range_index[tmp_index]))
            out_feature.SetField(str('classif'), int(self.classif[tmp_index]))
            out_feature.SetField(str('water_frac'), float(self.water_frac[tmp_index]))
            out_feature.SetField(str('xtrack'), float(self.cross_track[tmp_index]))
            out_feature.SetField(str('xtrack2'), float(self.cross_track_refined[tmp_index]))
            out_feature.SetField(str('sig0'), float(self.sig0[tmp_index]))
            out_feature.SetField(str('height'), float(self.height[tmp_index]))
            out_feature.SetField(str('pixel_area'), float(self.pixel_area[tmp_index]))
            out_feature.SetField(str('nadir_t'), float(self.nadir_time[tmp_index]))
            out_feature.SetField(str('nadir_long'), float(self.nadir_longitude[tmp_index]))
            out_feature.SetField(str('nadir_lat'), float(self.nadir_latitude[tmp_index]))
            # 2.4 - On ajoute l'objet dans la couche de sortie
            out_layer.CreateFeature(out_feature)

        # 3 - Destroy the data sources to free resources
        out_data_source.Destroy()
