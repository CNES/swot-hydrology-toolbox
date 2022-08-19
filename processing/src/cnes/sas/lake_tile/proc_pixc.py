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

import cnes.common.service_config_file as service_config_file

import cnes.common.lib.my_netcdf_file as my_nc
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_segmentation as my_segmentation
import cnes.common.lib.my_variables as my_var
import cnes.common.lib_lake.locnes_products_netcdf as nc_file

import jpl.modules.aggregate as jpl_aggregate


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
            - origin_latitude / 1D-array of float: latitude of water pixels
            - origin_longitude / 1D-array of float: longitude of water pixels
            
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
            - phase_noise_std / 1D-array of float: phase noise standard deviation
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
            - classification_qual / 1D array of float: flag that indicates the quality of the classification quantities 
            - geolocation_qual / 1D array of float: flag that indicates the quality of the geolocation quantities 
            - wavelength / float: wavelength corresponding to the effective radar carrier frequency 
            - looks_to_efflooks / float: ratio between the number of actual samples and the effective number of independent samples
            during spatial averaging over a large 2-D area
        - From tvp group in L2_HR_PIXC file
            - nadir_time[_tai] / 1D-array of float: observation UTC [TAI] time of each nadir pixel (= variable named time[tai] in L2_HR_PIXC file)
            - nadir_time_all / 1D-array of float: observation UTC time at nadir corresponding to all PIXC (= variable named time in L2_HR_PIXC file)
            - nadir_longitude / 1D-array of float: longitude of each nadir pixel (= variable named longitude in L2_HR_PIXC file)
            - nadir_latitude / 1D-array of float: latitude of each nadir pixel (= variable named latitude in L2_HR_PIXC file)   
            - nadir_[x|y|z] / 1D-array of float: [x|y|z] cartesian coordinates of each nadir pixel (= variables named [x|y|z] in L2_HR_PIXC file)
            - nadir_[vx|vy|vz] / 1D-array of float: velocity vector of each nadir pixel in cartesian coordinates (= variables 
              named velocity_unit_[x|y|z] in L2_HR_PIXC file)
            - nadir_plus_y_antenna_[x|y|z] / 1D-array of float: position vector of the +y KaRIn antenna phase center in ECEF coordinates 
              (= variables named plus_y_antenna_[x|y|z] in L2_HR_PIXC file)
            - nadir_minus_y_antenna_[x|y|z] / 1D-array of float: position vector of the -y KaRIn antenna phase center in ECEF coordinates 
              (= variables named minus_y_antenna_[x|y|z] in L2_HR_PIXC file)
        - From processing
            - tile_poly / ogr.Polygon: polygon of the PixC tile
            - az_0_line / ogr.LineString: line with azimuth = 0
            - az_max_line / ogr.LineString: line with azimuth = az_max
            - continent / string: continent covered by the tile (if global var CONTINENT_FILE exists)
            - selected_index / 1D-array of int: indices from original 1D-arrays of not rejected pixels with specified classification indices
            - nb_selected / int: number of selected pixels (=selected_index.size)
            - nb_water_pix / int: number of water pixels
            - list_classif_keys / list: list of keys of classif_dict
            - classif_dict/ dict: selection (by True or False) of pixels wrt their classification, for different processes; key=water/dark/as_full_water/without_dw/4hull/4wse/4area
            - interferogram_flattened / 1D-array of complex: flattened interferogram
            - inundated_area / 1D-array of int: area of pixels covered by water
            - height_std_pix / 1D-array of float: height std
            - corrected_height / 1D-array of float: height corrected from geoid and other tide corrections
            - labels / 1D-array of int: labelled regions associated to each PixC water pixel; pixels of this vector correspond one-to-one 
              to the pixels of data from L2_HR_PIXC and L2_HR_PIXCVecRiver
            - nb_obj / int : number of separate entities in the PixC tile
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

        # Init PIXC variables
        self.origin_classif = None  
        self.origin_range_index = None
        self.origin_azimuth_index = None
        self.origin_longitude = None
        self.origin_latitude = None
        self.classif = None  
        self.range_index = None
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
        self.pixel_area = None
        self.phase_noise_std = None
        self.dlatitude_dphase = None
        self.dlongitude_dphase = None
        self.dheight_dphase = None
        self.darea_dheight = None
        self.eff_num_medium_looks = None
        self.model_dry_tropo_cor = None
        self.model_wet_tropo_cor = None
        self.iono_cor_gim_ka = None
        self.height_cor_xover = None
        self.geoid = None
        self.solid_earth_tide = None
        self.load_tide_fes = None
        self.load_tide_got = None
        self.pole_tide = None
        self.classification_qual = None
        self.geolocation_qual = None
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
        self.pixc_metadata["interferogram_size_range"] = -9999
        self.pixc_metadata["interferogram_size_azimuth"] = -9999
        self.pixc_metadata["looks_to_efflooks"] = -9999.0
        self.pixc_metadata["ellipsoid_semi_major_axis"] = ""
        self.pixc_metadata["ellipsoid_flattening"] = ""

        # Variables specific to processing
        self.nb_pix_range = 0  # Number of pixels in range dimension
        self.nb_pix_azimuth = 0  # Number of pixels in azimuth dimension
        self.tile_poly = None  # Polygon representing the PIXC tile
        self.az_0_line = None # Line with azimuth = 0
        self.az_max_line  = None # Line with azimuth = az_max
        self.continent_id = None  # 2-letter identifier of the continent covered by the tile
        self.continent_code = None  # 1-digit code of the continent covered by the tile
        self.selected_index = None  # Indices of selected pixels
        self.nb_selected = 0  # Number of selected pixels
        self.nb_water_pix = 0  # Number of water pixels
        self.list_classif_keys = ["water", "interior_water", "as_full_water", "dark", "without_dw", "4hull", "4wse", "4area"]
        self.classif_dict = None  # Selection of pixels per process, wrt their classification
        self.interferogram_flattened = None  # Flattened interferogram
        self.inundated_area = None  # Area of pixel where water
        self.height_std_pix = None  # Height std
        self.labels = None  # Vector of entity labels associated to each pixel
        self.nb_obj = None  # Number of separate entities
        self.labels_inside = None  # Labels of entities entirely inside the tile
        self.nb_obj_inside = None  # Number of entities inside the tile
        self.labels_at_top_edge = None  # Labels of entities at the top edge of the tile
        self.nb_obj_at_top_edge = None  # Number of entities at the top edge of the tile
        self.labels_at_bottom_edge = None  # Labels of entities at the bottom edge of the tile
        self.nb_obj_at_bottom_edge = None  # Number of entities at the bottom edge of the tile
        self.labels_at_both_edges = None  # Labels of entities at the top and bottom edges of the tile
        self.nb_obj_at_both_edges = None  # Number of entities at the top and bottom edges of the tile
        self.edge_index = None  # Indices of pixels contained in objects at top/bottom edges
        self.edge_label = None  # Object label for each pixel contained in objects at top/bottom edges
        self.edge_loc = None  # Object edge location (0=bottom 1=top 2=both) for each pixel contained in objects at top/bottom edges
        self.nb_edge_pix = 0  # Number of pixels contained in objects at top/bottom edges
        
    def set_from_pixc_file(self, in_pixc_file, in_index_reject):
        """
        Retrieve needed data from pixel cloud
        
        :param in_pixc_file: full path of L2_HR_PIXC file
        :type in_pixc_file: string
        :param in_index_reject: list of indices to reject before all processing
        :type in_index_reject: 1D-array of int
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("L2_HR_PIXC file = %s" % in_pixc_file)
        
        # Get flag variables from configuration file
        classif_list = [int(i) for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_LIST").split(";")]
        classif_water = [i == "True" for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_WATER").split(";")]
        classif_4hull = [i == "True" for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_4HULL").split(";")]
        classif_4wse = [i == "True" for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_4WSE").split(";")]
        classif_4area = [i == "True" for i in self.cfg.get("CONFIG_PARAMS", "CLASSIF_4AREA").split(";")]
        use_fractional_inundation = [i == "True" for i in self.cfg.get("CONFIG_PARAMS", "USE_FRACTIONAL_INUNDATION").split(";")]
        pixc_qual_list = self.cfg.get("CONFIG_PARAMS", "PIXC_QUAL_LIST").split(";")
        pixc_qual_mask = self.cfg.get("CONFIG_PARAMS", "PIXC_QUAL_MASK").split(";")
        # If EXCLUDE_BRIGHT_LAND in lake_tile_param.cfg: retrieved as a string => getboolean needed to convert as a boolean
        # If not: init as a boolean => must use get and not getboolean
        try:
            exclude_bright_land = self.cfg.getboolean("CONFIG_PARAMS", "EXCLUDE_BRIGHT_LAND")
        except:
            exclude_bright_land = self.cfg.get("CONFIG_PARAMS", "EXCLUDE_BRIGHT_LAND")
        
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
                    
        # Savec of specific variables
        self.nb_pix_range = int(self.pixc_metadata["interferogram_size_range"])  # Number of pixels in range dimension
        self.nb_pix_azimuth = int(self.pixc_metadata["interferogram_size_azimuth"])  # Number of pixels in azimuth dimension
        self.near_range = np.double(self.pixc_metadata["near_range"])  # Slant range for the 1st image pixel
        self.wavelength = np.float(self.pixc_metadata["wavelength"])  # Wavelength
        # Ratio between the number of actual samples and the effective number of independent samples
        self.looks_to_efflooks = np.float(self.pixc_metadata["looks_to_efflooks"])  

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
        # 4.7 - Illumination time
        origin_illumination_time = pixc_reader.get_var_value("illumination_time", in_group=pixc_group)
        # 4.8 - PixC quality flags to consider
        if pixc_qual_list is not None:
            origin_qual = dict()
            for flag_name in pixc_qual_list:
                origin_qual[flag_name] = pixc_reader.get_var_value(flag_name, in_group=pixc_group)
        # 4.9 - bright_land_flag
        origin_bright_land_flag = pixc_reader.get_var_value("bright_land_flag", in_group=pixc_group)
                
        # 5 - Select PixC given criteria (Non-river, Finite values, quality flags)
        logger.debug("There are %d pixels in the input PIXC" % len(origin_longitude))
        
        # PixC excluded from LakeTile processing will have their classification flag set to 100
        tmp_classif = self.origin_classif  # Init temporary classif vector
        
        # 5.1 - Exclude river-only pixels
        if in_index_reject is not None:
            logger.debug("%d pixels are river (but not connected lakes) pixels => will be rejected" % in_index_reject.size)
            tmp_classif[in_index_reject] = 100
            
        # 5.2.1 - Exclude pixels with longitude, latitude, height and illumination_time = numpy.nan
        vars_to_look_for_nans = {}
        vars_to_look_for_nans["longitude"] = origin_longitude
        vars_to_look_for_nans["latitude"] = origin_latitude
        vars_to_look_for_nans["height"] = origin_height
        vars_to_look_for_nans["time"] = origin_illumination_time
        for var_name, var_value in vars_to_look_for_nans.items():
            nan_idx = np.argwhere(np.isnan(var_value))
            nb_nan = len(nan_idx)
            if nb_nan != 0:
                logger.debug("%d pixels have NaN in %s variable => will be rejected" % (nb_nan, var_name))
                tmp_classif[nan_idx] = 100 
        # 5.2.2 - Exclude pixels with classification, range_index and azimuth_index = _FillValue
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
            if nb_fv != 0:
                logger.debug("%d pixels have _FillValue in %s variable => will be rejected" % (nb_fv, var_name))
                tmp_classif[fv_idx] = 100 
                
        # 5.3 - Exclude pixels corresponding to the unwanted quality masks
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
                        logger.debug("%d pixels correspond to the unwanted quality mask [%s=%d] for %s quality flag => will be rejected" % (nb_nok, flag_mask, mask_nok, flag_name))
                    tmp_classif[nok_idx] = 100 
                    
        # 5.4 - Exclude pixels having bright_land_flag = 1
        if exclude_bright_land:
            nok_idx = np.argwhere(origin_bright_land_flag == 1)
            nb_nok = len(nok_idx)
            if nb_nok == 0:
                logger.debug("No pixel has bright_land_flag=1 => ALL pixels KEPT")
            else:
                logger.debug("%d pixels have bright_land_flag=1 => will be rejected" % nb_nok)
                tmp_classif[nok_idx] = 100 
        
        # 5.5 - Get indices of selected pixels
        self.selected_index = None  # Init wanted indices vector
        for classif_flag in classif_list:
            v_ind = np.where(tmp_classif == classif_flag)[0]
            logger.debug("%d pixels with classification flag = %d" % (v_ind.size, classif_flag))
            if v_ind.size != 0:
                if self.selected_index is None:
                    self.selected_index = v_ind
                else:
                    self.selected_index = np.concatenate((self.selected_index, v_ind))
        if self.selected_index is None:
            self.nb_selected = 0
        else:
            self.nb_selected = self.selected_index.size
        logger.debug("=> %d pixels to keep" % self.nb_selected)

        # 6 - Retrieve PixC data only for selected pixels
        if self.nb_selected == 0:
            # Interpolate nadir_time wrt illumination time
            tmp_nadir_time = pixc_reader.get_var_value("time", in_group=sensor_group)  # Read nadir_time values

            if len(tmp_nadir_time) > 1:
                f = interpolate.interp1d(tmp_nadir_time, range(len(tmp_nadir_time)))  # Interpolator
                nadir_index_all = (np.rint(f(origin_illumination_time))).astype(
                    int)  # Link between all PixC and nadir pixels
            else:
                nadir_index_all = np.zeros(len(origin_illumination_time), dtype=int)

            # Nadir time
            self.nadir_time_all = tmp_nadir_time[nadir_index_all]

        elif self.nb_selected != 0:
            
            # 6.1 - In PixC group
            
            # Classification flags
            self.classif = self.origin_classif[self.selected_index]
            
            # PIXC indices wrt category of process
            self.classif_dict = dict()
            for key in self.list_classif_keys:
                self.classif_dict[key] = np.full((self.nb_selected), False)
            
            # Simulate classification of edge + full water pixels
            # All PIXC are set to INTERIOR_WATER
            # LAND_EDGE + WATER_EDGE PIXC are set to WATER_EDGE
            self.classif_dict["as_full_water"] = np.zeros(np.shape(self.classif)) + my_var.CLASSIF_INTERIOR_WATER
            self.classif_dict["as_full_water"][self.classif == my_var.CLASSIF_LAND_EDGE] = my_var.CLASSIF_WATER_EDGE
            self.classif_dict["as_full_water"][self.classif == my_var.CLASSIF_WATER_EDGE] = my_var.CLASSIF_WATER_EDGE
            # Keep only classification of water pixels (ie remove dark water flags)
            self.classif_dict["without_dw"] = np.copy(self.classif)
            self.classif_dict["without_dw"][self.classif == my_var.CLASSIF_DARK] = 0
            
            # Boolean arrays of classification flags (=True if PIXC is of classif "key")
            tmp_classif_dict = dict()
            for classif_flag in classif_list:
                tmp_classif_dict[classif_flag] = np.full((self.nb_selected), False)
                tmp_classif_dict[classif_flag][self.classif == classif_flag] = True
                if classif_flag == my_var.CLASSIF_INTERIOR_WATER:
                    self.classif_dict["interior_water"] = tmp_classif_dict[classif_flag]
            # PIXC indices of water and dark water pixels
            for classif_flag, ok_water in zip(classif_list, classif_water):
                if ok_water:
                    self.classif_dict["water"] += tmp_classif_dict[classif_flag]
                else:
                    self.classif_dict["dark"] += tmp_classif_dict[classif_flag]
            # PIXC indices of pixels to be used for hull computation
            for classif_flag, ok_4hull in zip(classif_list, classif_4hull):
                if ok_4hull:
                    self.classif_dict["4hull"] += tmp_classif_dict[classif_flag]
            # PIXC indices of pixels to be used for wse computation
            for classif_flag, ok_4wse in zip(classif_list, classif_4wse):
                if ok_4wse:
                    self.classif_dict["4wse"] += tmp_classif_dict[classif_flag]
            # PIXC indices of pixels to be used for area_total attribute computation
            for classif_flag, ok_4area in zip(classif_list, classif_4area):
                if ok_4area:
                    self.classif_dict["4area"] += tmp_classif_dict[classif_flag]
            
            # Range indices of water pixels
            self.range_index = self.origin_range_index[self.selected_index]
            # Number of water pixels
            self.nb_water_pix = self.range_index.size
            # Range indices of water pixels
            self.azimuth_index = self.origin_azimuth_index[self.selected_index]
            
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
            self.cross_track = pixc_reader.get_var_value("cross_track", in_group=pixc_group)[self.selected_index]
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
            
            # Quality flag for the classification quantities 
            self.classification_qual = origin_qual["classification_qual"][self.selected_index]
            # Quality flag for the geolocation quantities
            self.geolocation_qual = origin_qual["geolocation_qual"][self.selected_index]

            # 6.2 - In TVP group
            
            # Interpolate nadir_time wrt illumination time
            tmp_nadir_time = pixc_reader.get_var_value("time", in_group=sensor_group)  # Read nadir_time values
   
            if len(tmp_nadir_time) > 1:
                f = interpolate.interp1d(tmp_nadir_time, range(len(tmp_nadir_time)))  # Interpolator
                nadir_index = (np.rint(f(illumination_time))).astype(int)  # Link between PixC and nadir pixels
                nadir_index_all = (np.rint(f(origin_illumination_time))).astype(int)  # Link between all PixC and nadir pixels
            else:
                nadir_index = np.zeros(len(illumination_time), dtype=int)
                nadir_index_all = np.zeros(len(origin_illumination_time), dtype=int)
                
            # Nadir time
            self.nadir_time = tmp_nadir_time[nadir_index]
            self.nadir_time_all = tmp_nadir_time[nadir_index_all]
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
            self.corrected_height = np.zeros(self.nb_water_pix) + np.nan
            self.corrected_height[ind_valid_corr] = self.height[ind_valid_corr] \
                                                    - self.geoid[ind_valid_corr] \
                                                    - self.solid_earth_tide[ind_valid_corr] \
                                                    - self.pole_tide[ind_valid_corr] \
                                                    - self.load_tide_fes[ind_valid_corr]
                    
        # 7 - Close file
        pixc_reader.close()

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

    def compute_separate_entities(self):
        """
        Identify all separate entities in the water mask
        """
        cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")

        # 1 - Create the water mask
        water_mask = self.compute_water_mask()

        # 2 - Identify all separate entities in a 2D binary mask
        sep_entities, self.nb_obj = my_tools.label_region(water_mask)

        # 3 - Convert 2D labelled mask in 1D-array layer of the same size of the L2_HR_PIXC layers
        self.labels = my_tools.convert_2d_mat_in_1d_vec(self.range_index, self.azimuth_index, sep_entities)
        self.labels = self.labels.astype(int)  # Conversion from float to integer

        # 4 - Relabel Lake Using Segmentation Heigth
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

                if np.sum(subset_pixel_area) > min_object_size*2:
                    min_rg = min(self.range_index[idx])
                    min_az = min(self.azimuth_index[idx])
                    subset_range = self.range_index[idx] - min_rg
                    subset_azimuth = self.azimuth_index[idx] - min_az
                    subset_corrected_height = self.corrected_height[idx]

                    relabel_obj = my_segmentation.relabel_lake_using_segmentation_heigth(subset_range, subset_azimuth, subset_corrected_height,
                                                                                         subset_pixel_area, min_object_size, seg_method)

                    labels_tmp[self.labels == label] = np.max(labels_tmp) + relabel_obj
                else:
                    labels_tmp[self.labels == label] = np.max(labels_tmp) + 1

            self.labels = labels_tmp
        self.nb_obj = np.unique(self.labels).size

    def compute_obj_inside_tile(self, obj_lake_db):
        """
        Separate labels of lakes and unknown objects entirely inside the tile, from labels of objects at top or bottom of the tile.
        Objects intersecting PLD lakes located at the edge of the tile are not considered as 'inside the tile'.

        :param obj_lake_db: Prior Lake Database (PLD) object
        :type obj_lake_db: lib_lake.lake_db

        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")

        # 1 - Identify objects at azimuth = 0
        idx_at_az_0 = np.where(self.azimuth_index == 0)[0]
        labels_at_az_0 = np.unique(self.labels[idx_at_az_0])

        # 2 - Identify objects at azimuth = max
        idx_at_az_max = np.where(self.azimuth_index == self.nb_pix_azimuth - 1)[0]
        labels_at_az_max = np.unique(self.labels[idx_at_az_max])

        # 3 - Identify objects at azimuth = 0 and azimuth = max
        self.labels_at_both_edges = np.intersect1d(labels_at_az_0, labels_at_az_max)
        self.nb_obj_at_both_edges = self.labels_at_both_edges.size

        # 4 - Identify labels...
        # 4.1 - Only at azimuth = 0
        self.labels_at_bottom_edge = np.setdiff1d(labels_at_az_0, self.labels_at_both_edges)
        # 4.2 - Only at azimuth = max
        self.labels_at_top_edge = np.setdiff1d(labels_at_az_max, self.labels_at_both_edges)

        # 5 - Identify labels of objects intersecting az_0_geom, az_max_geom, az_0_and_max_geom
        for l in np.unique(self.labels):
            idx = np.where(self.labels == l)

            bb_geom = my_tools.get_bounding_box(np.min(self.latitude[idx]), np.max(self.latitude[idx]),
                                                np.min(self.longitude[idx]), np.max(self.longitude[idx]))

            if bb_geom.Intersects(obj_lake_db.az_0_and_max_geom):
                obj_lake_db.remove_from_list_lake_id(bb_geom)
                if l not in self.labels_at_both_edges:
                    self.labels_at_both_edges = np.append(self.labels_at_both_edges, l)
            elif bb_geom.Intersects(obj_lake_db.az_max_geom):
                obj_lake_db.remove_from_list_lake_id(bb_geom)
                if l not in self.labels_at_top_edge:
                    self.labels_at_top_edge = np.append(self.labels_at_top_edge, l)
            elif bb_geom.Intersects(obj_lake_db.az_0_geom):
                obj_lake_db.remove_from_list_lake_id(bb_geom)
                if l not in self.labels_at_bottom_edge:
                    self.labels_at_bottom_edge = np.append(self.labels_at_bottom_edge, l)

        # 6 - Identify labels of objects with nb_tiles_az > MAX_NB_TILES_FULL_AZ
        for l in np.unique(self.labels):
            idx = np.where(self.labels == l)
            lon_tmp = self.longitude[idx]
            lat_tmp = self.latitude[idx]

            # Create memory layer to store pixc points
            ds, lyr = my_tools.load_pixels_to_mem_layer(lon_tmp, lat_tmp)
            for lake_id, geom in obj_lake_db.big_lake_geom:
                lyr.SetSpatialFilter(geom)
                if lyr.GetFeatureCount() > 0 :
                    logger.info("Lake %s is processed in LakeTile even if crosses azimuth border" % lake_id)
                    if l in self.labels_at_both_edges:
                        self.labels_at_both_edges = np.setdiff1d(self.labels_at_both_edges, [l])
                    if l in self.labels_at_top_edge:
                        self.labels_at_top_edge = np.setdiff1d(self.labels_at_top_edge, [l])
                    if l in self.labels_at_bottom_edge:
                        self.labels_at_bottom_edge = np.setdiff1d(self.labels_at_bottom_edge, [l])
                lyr.SetSpatialFilter(None)
            ds.Destroy()

        self.nb_obj_at_both_edges = self.labels_at_both_edges.size
        logger.debug("> %d labels at bottom (az=0) AND top (az=%d) of the tile" % (self.nb_obj_at_both_edges, self.nb_pix_azimuth))
        logger.debug(";".join(list(self.labels_at_both_edges.astype('str'))))

        self.nb_obj_at_top_edge = self.labels_at_top_edge.size
        logger.debug("> %d labels at top of the tile (az=%d)" % (self.nb_obj_at_top_edge, self.nb_pix_azimuth))
        logger.debug(";".join(list(self.labels_at_top_edge.astype('str'))))

        self.nb_obj_at_bottom_edge = self.labels_at_bottom_edge.size
        logger.debug("> %d labels at bottom of the tile (az=0)" % self.nb_obj_at_bottom_edge)
        logger.debug(";".join(list(self.labels_at_bottom_edge.astype('str'))))

        # 7 - Get labels of objects entirely inside the tile
        self.labels_inside = np.arange(self.nb_obj) + 1  # Initialisation
        if self.nb_obj_at_bottom_edge != 0:
            self.labels_inside = np.setdiff1d(self.labels_inside, self.labels_at_bottom_edge)  # Delete labels at bottom
        if self.nb_obj_at_top_edge != 0:
            self.labels_inside = np.setdiff1d(self.labels_inside, self.labels_at_top_edge)  # Delete labels at top
        if self.nb_obj_at_both_edges != 0:
            self.labels_inside = np.setdiff1d(self.labels_inside,
                                              self.labels_at_both_edges)  # Delete labels at top and bottom
        self.nb_obj_inside = self.labels_inside.size
        logger.debug("> %d objects entirely inside the tile" % self.nb_obj_inside)
        
    def compute_edge_indices_and_label(self):
        """
        Compute edge pixels indices and their associated label
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        if (self.nb_obj_at_top_edge + self.nb_obj_at_bottom_edge + self.nb_obj_at_both_edges) == 0:
            logger.debug("NO edge pixel to deal with")
            
        else:
            
            flag_start = False
            
            # 1 - Fill with bottom edge objects
            if self.nb_obj_at_bottom_edge != 0:
                flag_start = True
                self.edge_index = np.where(self.labels == self.labels_at_bottom_edge[0])[0]  # Get PixC pixels related to bottom edge object
                self.edge_label = np.ones(self.edge_index.size) * self.labels_at_bottom_edge[0]  # Associated label vector
                self.edge_loc = np.zeros(self.edge_index.size)  # Associated location: 0=bottom 1=top 2=both
                for indl in np.arange(1, self.nb_obj_at_bottom_edge):  # Loop over the other bottom edge objects
                    tmp_index = np.where(self.labels == self.labels_at_bottom_edge[indl])[0]  # Get pixels related to edge object
                    self.edge_index = np.append(self.edge_index, tmp_index)
                    self.edge_label = np.append(self.edge_label, np.ones(tmp_index.size) * self.labels_at_bottom_edge[indl])
                    self.edge_loc = np.append(self.edge_loc, np.zeros(tmp_index.size))
                    
            # 2 - Fill with top edge objects
            if self.nb_obj_at_top_edge != 0:
                if flag_start is False:
                    flag_start = True
                    idx_start = 1
                    self.edge_index = np.where(self.labels == self.labels_at_top_edge[0])[0]  # Get PixC pixels related to top edge object
                    self.edge_label = np.ones(self.edge_index.size) * self.labels_at_top_edge[0]  # Associated label vector
                    self.edge_loc = np.zeros(self.edge_index.size) + 1  # Associated location: 0=bottom 1=top 2=both
                else:
                    idx_start = 0
                for indl in np.arange(idx_start, self.nb_obj_at_top_edge):  # Loop over the other top edge objects
                    tmp_index = np.where(self.labels == self.labels_at_top_edge[indl])[0]  # Get pixels related to edge object
                    self.edge_index = np.append(self.edge_index, tmp_index)
                    self.edge_label = np.append(self.edge_label, np.ones(tmp_index.size) * self.labels_at_top_edge[indl])
                    self.edge_loc = np.append(self.edge_loc, np.zeros(tmp_index.size) + 1)
                    
            # 3 - Fill with bottom and top edges objects
            if self.nb_obj_at_both_edges != 0:
                if flag_start is False:
                    idx_start = 1
                    self.edge_index = np.where(self.labels == self.labels_at_both_edges[0])[0]  # Get PixC pixels related to both edges object
                    self.edge_label = np.ones(self.edge_index.size) * self.labels_at_both_edges[0]  # Associated label vector
                    self.edge_loc = np.zeros(self.edge_index.size) + 2  # Associated location: 0=bottom 1=top 2=both
                else:
                    idx_start = 0
                for indl in np.arange(idx_start, self.nb_obj_at_both_edges):  # Loop over the other both edges objects
                    tmp_index = np.where(self.labels == self.labels_at_both_edges[indl])[0]  # Get pixels related to edge object
                    self.edge_index = np.append(self.edge_index, tmp_index)
                    self.edge_label = np.append(self.edge_label, np.ones(tmp_index.size) * self.labels_at_both_edges[indl])
                    self.edge_loc = np.append(self.edge_loc, np.zeros(tmp_index.size) + 2)
                    
            # 4 - Number of edge pixels
            self.nb_edge_pix = self.edge_index.size

    # ----------------------------------------
    
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
            logger.warning("All pixels have height=NaN => output height is set to NaN")
            out_height = np.nan
            
        else:            
            nb_pixc = len(in_pixc_index)
            
            if nb_not_nan == nb_pixc:
                logger.debug("There are NO pixels having height=NaN => aggregate height only over all pixels")
                # Call JPL function
                out_height, weight_norm = jpl_aggregate.height_only(self.height, in_pixc_index, height_std=self.height_std_pix, method=method)
                
            else:
                logger.warning("There are %d/%d pixels having height=NaN => aggregate height only over these" % (nb_pixc-nb_not_nan, nb_pixc))
                # Call JPL function
                out_height, weight_norm = jpl_aggregate.height_only(self.height, in_pixc_index[not_nan_idx], height_std=self.height_std_pix,
                                                                    method=method)
        
        return out_height
    
    def compute_height_with_uncertainties(self, in_pixc_index, method='weight'):
        """
        Caller of JPL aggregate.py/height_with_uncerts function 
        which computes the aggregation of PIXC height over a feature, with corresponding uncertainty
        Compute height only over pixels having self.corrected_height != numpy.nan
    
        :param in_pixc_index: indices of pixels to consider for computation
        :type in_pixc_index: 1D-array of int
        :param method: type of aggregator ('weight', 'uniform', or 'median')
        :type method: string
        
        :return: out_height = aggregated height
        :rtype: out_height = float
        :return: out_height_std = standard deviation of the heights
        :rtype: out_height_std = float
        :return: out_height_unc = height uncertainty for the feature
        :rtype: out_height_unc = float
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Get pixels having self.corrected_height != numpy_nan
        not_nan_idx = np.argwhere(np.isfinite(self.corrected_height[in_pixc_index]))
        nb_not_nan = len(not_nan_idx)
        
        # Process depending their number
        if len(not_nan_idx) == 0:
            logger.warning("All pixels have corrected_height=NaN => output height is set to NaN")
            out_height = np.nan
            out_height_std = np.nan
            out_height_unc = np.nan

        else:            
            nb_pixc = len(in_pixc_index)
            
            if nb_not_nan == nb_pixc:
                logger.debug("There are NO pixels having corrected_height=NaN => aggregate height over all pixels")
                # Call JPL function
                out_height, out_height_std, out_height_unc, lat_unc, lon_unc = jpl_aggregate.height_with_uncerts(
                    self.corrected_height, in_pixc_index, self.eff_num_rare_looks, self.eff_num_medium_looks,
                    self.interferogram_flattened, self.power_plus_y, self.power_minus_y, self.looks_to_efflooks,
                    self.dheight_dphase, self.dlatitude_dphase, self.dlongitude_dphase, height_std=self.height_std_pix,
                    method=method)
                
            else:
                logger.warning("There are %d/%d pixels having corrected_height=NaN => aggregate height only over these" 
                               % (nb_pixc-nb_not_nan, nb_pixc))
                # Call JPL function
                out_height, out_height_std, out_height_unc, lat_unc, lon_unc = jpl_aggregate.height_with_uncerts(
                    self.corrected_height, in_pixc_index[not_nan_idx], self.eff_num_rare_looks, self.eff_num_medium_looks,
                    self.interferogram_flattened, self.power_plus_y, self.power_minus_y, self.looks_to_efflooks,
                    self.dheight_dphase, self.dlatitude_dphase, self.dlongitude_dphase, height_std=self.height_std_pix,
                    method=method)
        
        return out_height, out_height_std, out_height_unc

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
        
        # == TODO: compute the pixel assignment error?
        # call the general function
        
        # 0 - Select only pixels having self.pixel_area != numpy.nan
        not_nan_idx = np.argwhere(np.isfinite(self.pixel_area[in_pixc_index]))
        nb_not_nan = len(not_nan_idx)
        
        if len(not_nan_idx) == 0:
            logger.warning("All pixels have corrected_height=NaN => output area are set to NaN")
            out_area = np.nan
            out_area_unc = np.nan
            out_area_detct = np.nan
            out_area_detct_unc = np.nan
            
        else:            
            nb_pixc = len(in_pixc_index)
            
            if nb_not_nan == nb_pixc:
                logger.debug("There are NO pixels having corrected_height=NaN => aggregate area over all pixels")
                selected_pixc_idx = in_pixc_index
            else:
                logger.warning("There are %d/%d pixels having corrected_height=NaN => aggregate height only over these" 
                               % (nb_pixc-nb_not_nan, nb_pixc))
                selected_pixc_idx = in_pixc_index[not_nan_idx]

            # 1 - Aggregate PIXC area and uncertainty considering all PIXC (ie water + dark water)
            # Call external function common with RiverObs
            out_area, out_area_unc, area_pcnt_uncert = jpl_aggregate.area_with_uncert(
                    self.pixel_area, self.water_frac, self.water_frac_uncert,
                    self.darea_dheight, self.classif_dict["as_full_water"], self.false_detection_rate,
                    self.missed_detection_rate, selected_pixc_idx,
                    method=method)
            
            if flag_all:
            
                # 2 - Aggregate PIXC area and uncertainty considering only water PIXC (ie detected water)
                # Call external function common with RiverObs
                out_area_detct, out_area_detct_unc, area_detct_pcnt_uncert = jpl_aggregate.area_with_uncert(
                        self.pixel_area, self.water_frac, self.water_frac_uncert,
                        self.darea_dheight, self.classif_dict["without_dw"], self.false_detection_rate,
                        self.missed_detection_rate, selected_pixc_idx,
                        method=method)
                
            else:
                out_area_detct = np.nan
                out_area_detct_unc = np.nan
        
        return out_area, out_area_unc, out_area_detct, out_area_detct_unc
    
    def compute_geophysical_ref(self, in_name, in_pixc_index, method='weight'):
        """
        Caller of JPL aggregate.py/height_only function 
        which computes the aggregation of PIXC geophysical reference value over a feature
        Compute geophysical ref value only over pixels having self.<geophysical_ref> != numpy.nan
    
        :param in_name: name of the geophysical ref to process ('geoid_hght', 'solid_tide', 'pole_tide', 'load_tidef', or 'load_tideg')
        :type in_name: string
        :param in_pixc_index: indices of pixels to consider for computation
        :type in_pixc_index: 1D-array of int
        :param method: type of aggregator ('weight', 'uniform', or 'median')
        :type method: string
        
        :return: out_value = aggregated geophysical value
        :rtype: out_value = float
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Retrieve geophysical ref array to process
        geophysical_ref = None
        if in_name == 'geoid_hght':
            geophysical_ref = self.geoid
        elif in_name == 'solid_tide':
            geophysical_ref = self.solid_earth_tide
        elif in_name == 'pole_tide':
            geophysical_ref = self.pole_tide
        elif in_name == 'load_tidef':
            geophysical_ref = self.load_tide_fes
        elif in_name == 'load_tideg':
            geophysical_ref = self.load_tide_got
        else:
            logger.error("in_name = %d => should be among geoid_hght, solid_tide, pole_tide, load_tidef, or load_tideg" % in_name)
            raise
        
        # Get pixels having geophysical_ref != numpy_nan
        not_nan_idx = np.argwhere(np.isfinite(geophysical_ref[in_pixc_index]))
        nb_not_nan = len(not_nan_idx)
        
        # Process depending their number
        if len(not_nan_idx) == 0:
            logger.warning("All pixels have %s=NaN => output value is set to NaN" % in_name)
            out_value = np.nan
            
        else:            
            nb_pixc = len(in_pixc_index)
            
            if nb_not_nan == nb_pixc:
                logger.debug("There are NO pixels having %s=NaN => aggregate over all pixels" % in_name)
                # Call JPL function
                out_value, weight_norm = jpl_aggregate.height_only(geophysical_ref, in_pixc_index, height_std=self.height_std_pix, method=method)
                
            else:
                logger.warning("There are %d/%d pixels having %s=NaN => aggregate only over these" % (nb_pixc-nb_not_nan, nb_pixc, in_name))
                # Call JPL function
                out_value, weight_norm = jpl_aggregate.height_only(geophysical_ref, in_pixc_index[not_nan_idx], height_std=self.height_std_pix,
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
        
        # 1 - Init product
        edge_file = nc_file.LakeTileEdgeProduct(in_pixc_metadata=tmp_pixc_metadata, 
                                                in_proc_metadata={**in_proc_metadata, **time_str_dict})

        # 2 - Form dictionary with variables to write
        vars_to_write = {}
        if self.nb_selected != 0:
            vars_to_write["edge_index"] = self.selected_index[self.edge_index]
            vars_to_write["edge_label"] = self.edge_label
            vars_to_write["edge_loc"] = self.edge_loc
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
            vars_to_write["pixel_area"] = self.pixel_area[self.edge_index]
            vars_to_write["phase_noise_std"] = self.phase_noise_std[self.edge_index]
            vars_to_write["dlatitude_dphase"] = self.dlatitude_dphase[self.edge_index]
            vars_to_write["dlongitude_dphase"] = self.dlongitude_dphase[self.edge_index]
            vars_to_write["dheight_dphase"] = self.dheight_dphase[self.edge_index]
            vars_to_write["darea_dheight"] = self.darea_dheight[self.edge_index]
            vars_to_write["eff_num_medium_looks"] = self.eff_num_medium_looks[self.edge_index]
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
        tmp_field = ogr.FieldDefn(str('water_frac'), ogr.OFTReal)
        tmp_field.SetWidth(10)
        tmp_field.SetPrecision(6)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('crosstrack'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(4)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('pixel_area'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(6)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('height'), ogr.OFTReal)
        tmp_field.SetWidth(12)
        tmp_field.SetPrecision(4)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('nadir_t'), ogr.OFTReal)
        tmp_field.SetWidth(13)
        tmp_field.SetPrecision(3)
        out_layer.CreateField(tmp_field)
        tmp_field = ogr.FieldDefn(str('nadir_long'), ogr.OFTReal)
        tmp_field.SetWidth(10)
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
            out_feature.SetField(str('crosstrack'), float(self.cross_track[tmp_index]))
            out_feature.SetField(str('pixel_area'), float(self.pixel_area[tmp_index]))
            out_feature.SetField(str('height'), float(self.height[tmp_index]))
            out_feature.SetField(str('nadir_t'), float(self.nadir_time[tmp_index]))
            out_feature.SetField(str('nadir_long'), float(self.nadir_longitude[tmp_index]))
            out_feature.SetField(str('nadir_lat'), float(self.nadir_latitude[tmp_index]))
            # 2.4 - On ajoute l'objet dans la couche de sortie
            out_layer.CreateFeature(out_feature)

        # 3 - Destroy the data sources to free resources
        out_data_source.Destroy()
