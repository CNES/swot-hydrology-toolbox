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
.. module:: proc_pixc_sp.py
   :synopsis: Deals with subset of pixel cloud, for objects located at top or bottom edge of a tile; i.e. gather pixels involved in edge lake product retrieved from all tiles of L2_HR_LakeTile_edge files
    Created on 27/09/2017

.. moduleauthor:: Cécile Cazals - CS

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from __future__ import absolute_import, division, print_function, unicode_literals 

import numpy as np
from osgeo import ogr
import logging

import cnes.common.lib.my_basins as my_basins
import cnes.common.lib.my_netcdf_file as my_nc
import cnes.common.lib.my_tools as my_tools
import cnes.common.service_config_file as service_config_file

class PixC_Edge(object):
    """
        class PixC_Edge
    """
    def __init__(self, in_lake_tile_edge_path_list, in_cycle, in_pass, in_continent) :
        """
            This class is designed to process all L2_HR_LakeTile edge files of two swath thought 2 PixC_Edge_swath objects

            :param in_lake_tile_edge_path_list: list of LakeTile edge files full path
            :type in_lake_tile_edge_path_list: list of string
            :param in_cycle num: cycle number
            :type in_cycle: int
            :param in_pass_num: pass number
            :type in_pass_num: int
            :param in_continent: continent code
            :type in_continent: str

            Variables of the object:
                - lake_tile_edge_path_r / list of string : list of full path to laketile edge files for right swath
                - lake_tile_edge_path_l / list of string : list of full path to laketile edge files for left swath
                - pixc_edge_r / PixC_Edge_swath: object to process right swath
                - pixc_edge_l / PixC_Edge_swath: object to process left swath
                - tile_poly / ogr.Polygon: polygon of all tiles on both swath
                - nb_pixels / int: number of pixels to process in both swath

        """

        # Get instance of service config file
        self.cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)

        logger.info("- start -")

        self.lake_tile_edge_path_r = [ file for file in in_lake_tile_edge_path_list if "R_20" in file]
        self.lake_tile_edge_path_r.sort()
        self.lake_tile_edge_path_l = [ file for file in in_lake_tile_edge_path_list if "L_20" in file]
        self.lake_tile_edge_path_l.sort()

        logger.info("== INIT Right swath for continent %s with files ==" %(in_continent))
        for file in self.lake_tile_edge_path_r:
            logger.info(file)
        logger.info("===================================")
        self.pixc_edge_r = PixC_Edge_swath( in_cycle, in_pass, in_continent, "R")

        logger.info("== INIT Left swath for continent %s with files ==" %(in_continent))
        for file in self.lake_tile_edge_path_l:
            logger.info(file)
        logger.info("===================================")
        self.pixc_edge_l = PixC_Edge_swath( in_cycle, in_pass, in_continent, "L")

        self.pixc_metadata = {}  # dictionary of SP metadata

    def set_pixc_edge_from_laketile_edge_file(self):
        """
            Load pixc data from both swath
            Fill metadata

        """
        logger = logging.getLogger(self.__class__.__name__)

        # Load left swath
        self.pixc_edge_l.set_pixc_edge_swath_from_laketile_edge_file(self.lake_tile_edge_path_l)

        # Load right swath
        self.pixc_edge_r.set_pixc_edge_swath_from_laketile_edge_file(self.lake_tile_edge_path_r)

        # Update tile_poly with both swath
        self.tile_poly = self.pixc_edge_r.tile_poly.Union(self.pixc_edge_l.tile_poly)

        # Update nb_pixels with both swath
        self.nb_pixels = self.pixc_edge_r.nb_pixels + self.pixc_edge_l.nb_pixels

        # Assert that both swath belong to the same SP product (cycle, pass and continent) and fill metadata
        if self.pixc_edge_l.cycle_num == self.pixc_edge_r.cycle_num :
            self.pixc_metadata["cycle_number"] =  self.pixc_edge_l.cycle_num
        else :
            logger.error("Cycle number for swath left and right are not the same")

        if self.pixc_edge_l.pass_num == self.pixc_edge_r.pass_num :
            self.pixc_metadata["pass_number"] =  self.pixc_edge_l.pass_num
        else :
            logger.error("Pass number for swath left and right are not the same")

        if self.pixc_edge_l.continent == self.pixc_edge_r.continent :
            self.pixc_metadata["continent"] =  self.pixc_edge_l.continent
        else :
            logger.error("Continent for swath left and right are not the same")

        self.pixc_metadata["start_time"] = min(self.pixc_edge_l.pixc_metadata["start_time"], self.pixc_edge_r.pixc_metadata["start_time"])
        self.pixc_metadata["stop_time"] = max(self.pixc_edge_l.pixc_metadata["stop_time"], self.pixc_edge_r.pixc_metadata["stop_time"])

        self.pixc_metadata["polygon"] = str(self.tile_poly)

        # Assert that ellipsoid_semi_major_axis and ellipsoid_flattening are the same and fill metadat
        if self.pixc_edge_l.pixc_metadata["ellipsoid_semi_major_axis"] == self.pixc_edge_r.pixc_metadata["ellipsoid_semi_major_axis"] :
            self.pixc_metadata["ellipsoid_semi_major_axis"] =  self.pixc_edge_l.pixc_metadata["ellipsoid_semi_major_axis"]
        else :
            logger.error("Ellipsoid semi major axis for swath left and right are not the same")

        if self.pixc_edge_l.pixc_metadata["ellipsoid_flattening"] == self.pixc_edge_r.pixc_metadata["ellipsoid_flattening"] :
            self.pixc_metadata["ellipsoid_flattening"] =  self.pixc_edge_l.pixc_metadata["ellipsoid_flattening"]
        else :
            logger.error("Ellipsoid flattening for swath left and right are not the same")


class PixC_Edge_swath(object):
    """
        class PixC_Edge_swath
    """

    def __init__(self, in_cycle, in_pass, in_continent, in_swath_side):
        """
        This class is designed to process all L2_HR_LakeTile edge files of one swath. LakeTile edge files contain pixel cloud information (from L2_HR_PIXC) for pixels involved in lakes located at the top/bottom edges of tiles.

        The LakeTile edge file contains only PixC variables needed information, but also additional information like:

            - edge_loc field, the edge location of the lake : top of the tile, bottom of the tile or both top and bottom (0=bottom 1=top 2=both)
            - edge_label contains the object labels retrieved from PGE_LakeTile labeling process.
            - edge_index contains the L2_HR_PIXC initial pixels indices. This information is needed to update the improved geoloc and tag fields of LakeTile pixcvec files.

        This class processes all LakeTile edge files of a single swath to gather all entities at the edge of tile into lake entities.

        :param in_lake_tile_edge_path_list: list of LakeTile edge files full path
        :type in_lake_tile_edge_path_list: list of string
        :param in_cycle num: cycle number
        :type in_cycle: int
        :param in_pass_num: pass number
        :type in_pass_num: int
        :param in_continent: continent code
        :type in_continent: str
        :param in_swath_side: R=Right L=Left swath side
        :type in_swath_side: string

        Variables of the object:
            - Global PixC_Edge_swath infos:
                - cycle_num / int: cycle number (= global attribute named cycle_number in LakeTile_edge)
                - pass_num / int: pass number (= global attribute named pass_number in LakeTile_edge)
                - swath_side / string : R=Right L=Left swath side
                - date / int : date of acquisition (ex : 20000101 )
                - nb_pix_range / int: number of pixels in range dimension (= global attribute named nb_pix_range in LakeTile_edge)
                - nb_pixels / int : total number of pixels
                - pixc_metadata / dict : metadata processing ...
                - tile_poly / ogr.Polygon: polygon of the PixC tile

            - Variables specific to processing:
                - tile_num / list of int : List of tile number to process ex: [76, 77, 78]
                - tile_index /list of int : Tile_num reference of each pixel ex: [0, 0, 0, 1, 2, 2, 2, 2]
                - labels / 1D-array of int : arrays of new labels
                - is_boundary_pix / list of bool : if pixel belong to the first / last azimuth of single pass

            - Edge objects info:
                - edge_loc
                - edge_label
                - edge_index

            - Pixc variables from laketile edge file:
                - nb_pix_range / int: number of pixels in range dimension (= global attribute named interferogram_size_range in L2_HR_PIXC)
                - nb_pix_azimuth / int: number of pixels in azimuth dimension (= global attribute named interferogram_size_azimuth in L2_HR_PIXC)
                - [origin_]classif / 1D-array of byte: classification value of pixels (= variable named classification in L2_HR_PIXC)
                - [origin_]range_index / 1D-array of int: range indices of water pixels (= variable named range_index in L2_HR_PIXC)
                - [origin_]azimuth_index / 1D-array of int: azimuth indices of water pixels (= variable named azimuth_index in L2_HR_PIXC)
                - water_frac / 1D array of float: water fraction
                - water_frac_uncert / 1D array of float: water fraction uncertainty
                - false_detection_rate / 1D array of float: alse detection rate
                - missed_detection_rate / 1D array of float: missed detection rate
                - prior_water_prob / 1D array of float: prior water probability
                - bright_land_flag / 1D array of byte: bright land flag
                - layover_impact /1D array of float: layover impact
                - num_rare_looks / 1D array of byte: number of rare looks
                - [origin_]latitude / 1D-array of float: latitude of water pixels
                - [origin_]longitude / 1D-array of float: longitude of water pixels
                - height / 1D-array of float: height of water pixels
                - cross_track / 1D-array of float: cross-track distance from nadir to center of water pixels
                - pixel_area / 1D-array of int: area of water pixels
                - inc / 1D array of float: incidence angle
                - dheight_dphase / 1D array of float: sensitivity of height estimate to interferogram phase
                - dheight_droll / 1D array of float: sensitivity of height estimate to spacecraft roll
                - dheight_dbaseline / 1D array of float: sensitivity of height estimate to interferometric baseline
                - dheight_drange / 1D array of float: sensitivity of height estimate to range
                - darea_dheight / 1D array of float: sensitivity of pixel area to reference height
                - num_med_looks / 1D array of int: number of medium looks
                - sig0 / 1D array of float: sigma0
                - phase_unwrapping_region / 1D array of float: phase unwrapping region index
                - instrument_range_cor / 1D array of float: instrument range correction
                - instrument_phase_cor / 1D array of float: instrument phase correction
                - instrument_baseline_cor / 1D array of float: instrument baseline correction
                - instrument_attitude_cor / 1D array of float: instrument attitude correction
                - model_dry_tropo_cor / 1D array of float: dry troposphere vertical correction
                - model_wet_tropo_cor / 1D array of float: wet troposphere vertical correction
                - iono_cor_gim_ka / 1D array of float: ionosphere vertical correction
                - xover_height_cor / 1D array of float: crossover calibration height correction
                - load_tide_sol1 / 1D array of float: load tide height (GOT4.10)
                - load_tide_sol2 / 1D array of float: load tide height (FES2014)
                - pole_tide / 1D array of float: pole tide height
                - solid_earth_tide / 1D array of float: solid earth tide
                - geoid / 1D array of float: geoid
                - surface_type_flag / 1D array of byte: surface type flag

            - Info of the nadir point associated to the PixC:
                - nadir_time[_tai] / 1D-array of float: observation UTC [TAI] time of each nadir pixel (= variable named time[tai] in L2_HR_PIXC file)
                - nadir_longitude / 1D-array of float: longitude of each nadir pixel (= variable named longitude in L2_HR_PIXC file)
                - nadir_latitude / 1D-array of float: latitude of each nadir pixel (= variable named latitude in L2_HR_PIXC file)
                - nadir_[x|y|z] / 1D-array of float: [x|y|z] cartesian coordinates of each nadir pixel (= variables named [x|y|z] in L2_HR_PIXC file)
                - nadir_[vx|vy|vz] / 1D-array of float: velocity vector of each nadir pixel in cartesian coordinates (= variables named velocity_unit_[x|y|z] in L2_HR_PIXC file)
                - nadir_sc_event_flag / 1D array of byte: spacecraft event flag
                - nadir_tvp_qual / 1D array of byte: quality flag
        
        """
        # Get instance of service config file
        self.cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)

        logger.info("- start -")

        # 1. Init Global PixC_Edge_swath infos
        self.cycle_num = in_cycle # Cycle number
        self.pass_num = in_pass # Orbit number
        self.swath_side = in_swath_side # swath R or L
        self.continent = in_continent
        self.date = None # date of acquisition
        self.nb_pix_range = 0 # Number of pixel in range
        self.nb_pixels = 0 # Number of pixels to process
        self.tile_poly = ogr.Geometry(ogr.wkbMultiPolygon)  # Tiles polygon union
        self.pixc_metadata = {} # dictionary of SP edge metadata
        self.pixc_metadata["cycle_number"] = -9999
        self.pixc_metadata["pass_number"] = -9999
        self.pixc_metadata["start_time"] = ""
        self.pixc_metadata["stop_time"] = ""
        self.pixc_metadata["continent"] = self.continent
        self.pixc_metadata["ellipsoid_semi_major_axis"] = ""
        self.pixc_metadata["ellipsoid_flattening"] = ""

        # 2. Init variables specific to processing
        self.tile_num = [] # List of tile number to process ex: [76, 77, 78]
        self.tile_index = [] # Tile reference of each pixel
        self.labels =  np.array(()).astype('int') # Init labels to 0
        self.is_boundary_pix = [] # if pixel belong to the first / last azimuth of single pass

        # 3. Init edge objects info
        self.edge_loc =  np.array(()).astype('int') # Localization (top/bottom/both)
        self.edge_label = np.array(()).astype('int')  # Label in tile
        self.edge_index = np.array(()).astype('int')  # Index in original L2_HR_PIXC product

        # 4. Init Pixc variables from laketile edge file
        self.classif = np.array(())
        self.range_index = np.array(()).astype('int')
        self.azimuth_index = np.array(()).astype('int')
        self.pixel_area = np.array(())
        self.water_frac = np.array(())
        self.water_frac_uncert = np.array(())
        self.false_detection_rate = np.array(())
        self.missed_detection_rate = np.array(())
        self.prior_water_prob = np.array(())
        self.bright_land_flag = np.array(())
        self.layover_impact = np.array(())
        self.num_rare_looks = np.array(())
        self.latitude = np.array(())
        self.longitude = np.array(())
        self.height = np.array(())
        self.cross_track = np.array(())
        self.pixel_area = np.array(())
        self.inc = np.array(())
        self.dheight_dphase = np.array(())
        self.dheight_droll = np.array(())
        self.dheight_dbaseline = np.array(())
        self.dheight_drange = np.array(())
        self.darea_dheight = np.array(())
        self.num_med_looks = np.array(())
        self.sig0 = np.array(())
        self.phase_unwrapping_region = np.array(())
        self.instrument_range_cor = np.array(())
        self.instrument_phase_cor = np.array(())
        self.instrument_baseline_cor = np.array(())
        self.instrument_attitude_cor = np.array(())
        self.model_dry_tropo_cor = np.array(())
        self.model_wet_tropo_cor = np.array(())
        self.iono_cor_gim_ka = np.array(())
        self.xover_height_cor = np.array(())
        self.load_tide_sol1 = np.array(())
        self.load_tide_sol2 = np.array(())
        self.pole_tide = np.array(())
        self.solid_earth_tide = np.array(())
        self.geoid = np.array(())
        self.surface_type_flag = np.array(())
        self.pixc_qual = np.array(())

        # 5. Init info of the nadir point associated to the PixC
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
        self.nadir_sc_event_flag = np.array(())
        self.nadir_tvp_qual = np.array(())

    # ----------------------------------------

    def set_pixc_edge_swath_from_laketile_edge_file(self, in_lake_tile_edge_path_list):
        """
        This function loads NetCDF information form all tiles and merge of all PIXC_edge info into 1D vectors

        :param in_lake_tile_edge_path_list: list of full path of the NetCDF file to load
        :type in_lake_tile_edge_path_list: list of string

        """
        logger = logging.getLogger(self.__class__.__name__)

        in_lake_tile_edge_path_list.sort()
        nb_pix_azimuth_list = [0]
        for lake_tile_edge_path in in_lake_tile_edge_path_list:  # For each LakeTile edge file

            logger.info("Loading L2_HR_LakeTile edge file = %s" % lake_tile_edge_path)

            # Load data
            nb_pix_loaded, current_tile_number, current_nb_pix_azimuth = self.load_laketile_edge_data(lake_tile_edge_path)

            if not current_tile_number:
                continue

            logger.info("%d pixels loaded" % nb_pix_loaded)

            # Update nb_pixel
            self.nb_pixels += nb_pix_loaded

            # check if tile is neighboring the previous tile
            if len(self.tile_num) >= 1:
                previous_tile_number = self.tile_num[-1]

                cpt = 0
                while current_tile_number != previous_tile_number + 1:
                    logger.info("Adding empty tile %d%s" % (previous_tile_number + 1, self.swath_side))
                    # if current tile is not adjacent to previous tile, add an empty tile to tile_ref
                    self.tile_num.append(previous_tile_number + 1)
                    nb_pix_azimuth_list.append(current_nb_pix_azimuth)
                    previous_tile_number = self.tile_num[-1]
                    cpt += 1
                    if cpt > 50:
                        break

            self.tile_index += [current_tile_number] * nb_pix_loaded
            self.tile_num.append(current_tile_number)
            nb_pix_azimuth_list.append(current_nb_pix_azimuth)


        # Convert list to numpy array
        self.tile_index = np.array(self.tile_index)
        self.labels = np.zeros((self.nb_pixels))

        # Compute if pixels belong to the beging or the end of pass
        azimuth_increment_list = np.cumsum(nb_pix_azimuth_list)[:-1]
        full_path_azimuth = np.zeros((self.nb_pixels))

        for i, tile in enumerate(self.tile_num):
            full_path_azimuth[np.where(self.tile_index == tile)] = self.azimuth_index[np.where(self.tile_index == tile)] + azimuth_increment_list[i]

        self.is_boundary_pix = np.logical_or(full_path_azimuth == 0, full_path_azimuth == np.sum(nb_pix_azimuth_list)-1)

        logger.info("%d PixC loaded for current swath" % self.nb_pixels)

    # ----------------------------------------

    def load_laketile_edge_data(self, in_lake_tile_edge_filename):
        """
        This function loads NetCDF information.
        
        :param in_lake_tile_edge_filename: full path of the NetCDF file to load
        :type in_lake_tile_edge_filename: string
        
        :return out_nb_pix: number of pixel loaded
        :rtype out_nb_pix: integer
        :return OUT_tile_ref: reference of loaded tile
        :rtype OUT_tile_ref: string
        """

        # 1 - Open input NetCDF file in reading mode
        pixc_edge_reader = my_nc.myNcReader(in_lake_tile_edge_filename)

        # 2 - Get and check tile references (cycle, pass, swath)
        out_tile_number = int(pixc_edge_reader.getAttValue("tile_number"))

        current_cycle_num = int(pixc_edge_reader.getAttValue("cycle_number"))
        if current_cycle_num != self.cycle_num:
            logging.error("Cycle of tile %d do note match with SP product %d" %(current_cycle_num, self.cycle_num))

        current_pass_number = int(pixc_edge_reader.getAttValue("pass_number"))
        if current_pass_number != self.pass_num:
            logging.error("Pass of tile %d do note match with SP product %d" %(current_pass_number, self.pass_num))

        current_swath_side = str(pixc_edge_reader.getAttValue("swath_side"))
        if current_swath_side != self.swath_side:
            logging.error("Swath of tile %s do note match with PixC_Edge_swath %s" %(current_swath_side, self.swath_side))
        current_date = int(pixc_edge_reader.getAttValue("start_time")[:8])
        if not self.date :
            self.date = current_date
        # Stop the process if acquisition dates are not same day or next day
        if np.fabs(current_date - self.date) > 1 :
            logging.error("Input laketile_edge file do not correspond to the same aquisition date")

        # A supprimer à terme current_continent = my_basins.link_poly_to_continent(current_tile_poly)
        current_continent = str(pixc_edge_reader.getAttValue("continent"))
        if not current_continent == self.continent:
            # If cur_continent do not belong to the EDGE SP product, do not add pixc info
            logging.error("Input laketile_edge %s file do not correspond to the same continent %s" %(in_lake_tile_edge_filename, self.continent))
            retour = None, None, None

        else:
    
            # 4 - Update tile_poly with new tile
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(float(pixc_edge_reader.getAttValue("inner_first_longitude")),
                          float(pixc_edge_reader.getAttValue("inner_first_latitude")))
            ring.AddPoint(float(pixc_edge_reader.getAttValue("outer_first_longitude")),
                          float(pixc_edge_reader.getAttValue("outer_first_latitude")))
            ring.AddPoint(float(pixc_edge_reader.getAttValue("outer_last_longitude")),
                          float(pixc_edge_reader.getAttValue("outer_last_latitude")))
            ring.AddPoint(float(pixc_edge_reader.getAttValue("inner_last_longitude")),
                          float(pixc_edge_reader.getAttValue("inner_last_latitude")))
            ring.AddPoint(float(pixc_edge_reader.getAttValue("inner_first_longitude")),
                          float(pixc_edge_reader.getAttValue("inner_first_latitude")))
            current_tile_poly = ogr.Geometry(ogr.wkbPolygon)
            current_tile_poly.AddGeometry(ring)

            if not current_tile_poly.IsValid() :
                logging.warning("Polygon tile of file %s is not valid" %(in_lake_tile_edge_filename))
            else :
                self.tile_poly = self.tile_poly.Union(current_tile_poly)
    
            # 3 - Initialization of object variables if not already done
            if not self.tile_num:
                self.nb_pix_range = pixc_edge_reader.getAttValue("interferogram_size_range")
                # Find associated continent
                self.near_range = np.double(pixc_edge_reader.getAttValue("near_range"))  # Slant range for the 1st image pixel
    
                self.pixc_metadata["cycle_number"] = self.cycle_num
                self.pixc_metadata["pass_number"] = self.pass_num
                self.pixc_metadata["start_time"] = str(pixc_edge_reader.getAttValue("start_time"))
                self.pixc_metadata["stop_time"] = str(pixc_edge_reader.getAttValue("stop_time"))
                self.pixc_metadata["continent"] = self.continent
                self.pixc_metadata["ellipsoid_semi_major_axis"] = str(pixc_edge_reader.getAttValue("ellipsoid_semi_major_axis"))
                self.pixc_metadata["ellipsoid_flattening"] = str(pixc_edge_reader.getAttValue("ellipsoid_flattening"))
    
    
            # 5 - Get number of pixels
            out_nb_pix = pixc_edge_reader.getDimValue('points')
            out_nb_pix_azimuth = pixc_edge_reader.getAttValue("interferogram_size_azimuth")
    
            # 6 - Update vectors if there are pixels
            if out_nb_pix > 0:
    
                # 5.1 - Edge objects info
                self.edge_label = np.concatenate((self.edge_label, pixc_edge_reader.getVarValue("edge_label")))
                self.edge_index = np.concatenate((self.edge_index, pixc_edge_reader.getVarValue("edge_index")))
                self.edge_loc = np.concatenate((self.edge_loc, pixc_edge_reader.getVarValue("edge_loc")))
    
                # 5.2 - Variables from L2_HR_PIXC product
                self.classif = np.concatenate((self.classif, pixc_edge_reader.getVarValue("classification")))
                self.range_index = np.concatenate((self.range_index, pixc_edge_reader.getVarValue("range_index")))
                self.azimuth_index = np.concatenate((self.azimuth_index, pixc_edge_reader.getVarValue("azimuth_index")))
                self.pixel_area = np.concatenate((self.pixel_area, pixc_edge_reader.getVarValue("pixel_area")))
                self.water_frac = np.concatenate((self.pixel_area, pixc_edge_reader.getVarValue("water_frac")))
                self.water_frac_uncert = np.concatenate((self.pixel_area, pixc_edge_reader.getVarValue("water_frac_uncert")))
                self.false_detection_rate = np.concatenate((self.pixel_area, pixc_edge_reader.getVarValue("false_detection_rate")))
                self.missed_detection_rate = np.concatenate((self.pixel_area, pixc_edge_reader.getVarValue("missed_detection_rate")))
                self.prior_water_prob = np.concatenate((self.pixel_area, pixc_edge_reader.getVarValue("prior_water_prob")))
                self.bright_land_flag = np.concatenate((self.pixel_area, pixc_edge_reader.getVarValue("bright_land_flag")))
                self.layover_impact = np.concatenate((self.pixel_area, pixc_edge_reader.getVarValue("layover_impact")))
                self.num_rare_looks = np.concatenate((self.pixel_area, pixc_edge_reader.getVarValue("num_rare_looks")))
                self.latitude = np.concatenate((self.latitude, pixc_edge_reader.getVarValue("latitude")))
                self.longitude = np.concatenate((self.longitude, pixc_edge_reader.getVarValue("longitude")))
                self.height = np.concatenate((self.height, pixc_edge_reader.getVarValue("height")))
                self.cross_track = np.concatenate((self.cross_track, pixc_edge_reader.getVarValue("cross_track")))
                self.pixel_area = np.concatenate((self.pixel_area, pixc_edge_reader.getVarValue("pixel_area")))
                self.inc = np.concatenate((self.inc, pixc_edge_reader.getVarValue("inc")))
                self.dheight_dphase = np.concatenate((self.dheight_dphase, pixc_edge_reader.getVarValue("dheight_dphase")))
                self.dheight_droll = np.concatenate((self.dheight_droll, pixc_edge_reader.getVarValue("dheight_droll")))
                self.dheight_dbaseline = np.concatenate((self.dheight_droll, pixc_edge_reader.getVarValue("dheight_dbaseline")))
                self.dheight_drange = np.concatenate((self.dheight_drange, pixc_edge_reader.getVarValue("dheight_drange")))
                self.darea_dheight = np.concatenate((self.darea_dheight, pixc_edge_reader.getVarValue("darea_dheight")))
                self.num_med_looks = np.concatenate((self.num_med_looks, pixc_edge_reader.getVarValue("num_med_looks")))
                self.sig0 = np.concatenate((self.sig0, pixc_edge_reader.getVarValue("sig0")))
                self.phase_unwrapping_region = np.concatenate((self.phase_unwrapping_region, pixc_edge_reader.getVarValue("phase_unwrapping_region")))
                self.instrument_range_cor = np.concatenate((self.instrument_range_cor, pixc_edge_reader.getVarValue("instrument_range_cor")))
                self.instrument_phase_cor = np.concatenate((self.instrument_phase_cor, pixc_edge_reader.getVarValue("instrument_phase_cor")))
                self.instrument_baseline_cor = np.concatenate((self.instrument_baseline_cor, pixc_edge_reader.getVarValue("instrument_baseline_cor")))
                self.instrument_attitude_cor = np.concatenate((self.instrument_attitude_cor, pixc_edge_reader.getVarValue("instrument_attitude_cor")))
                self.model_dry_tropo_cor = np.concatenate((self.model_dry_tropo_cor, pixc_edge_reader.getVarValue("model_dry_tropo_cor")))
                self.model_wet_tropo_cor = np.concatenate((self.model_wet_tropo_cor, pixc_edge_reader.getVarValue("model_wet_tropo_cor")))
                self.iono_cor_gim_ka = np.concatenate((self.iono_cor_gim_ka, pixc_edge_reader.getVarValue("iono_cor_gim_ka")))
                self.xover_height_cor = np.concatenate((self.xover_height_cor, pixc_edge_reader.getVarValue("xover_height_cor")))
                self.load_tide_sol1 = np.concatenate((self.load_tide_sol1, pixc_edge_reader.getVarValue("load_tide_sol1")))
                self.load_tide_sol2 = np.concatenate((self.load_tide_sol2, pixc_edge_reader.getVarValue("load_tide_sol2")))
                self.pole_tide = np.concatenate((self.pole_tide, pixc_edge_reader.getVarValue("pole_tide")))
                self.solid_earth_tide = np.concatenate((self.solid_earth_tide, pixc_edge_reader.getVarValue("solid_earth_tide")))
                self.geoid = np.concatenate((self.geoid, pixc_edge_reader.getVarValue("geoid")))
                self.surface_type_flag = np.concatenate((self.surface_type_flag, pixc_edge_reader.getVarValue("surface_type_flag")))
                self.pixc_qual = np.concatenate((self.pixc_qual, pixc_edge_reader.getVarValue("pixc_qual")))
    
                # 5.3 - Info of the nadir point associated to the PixC
                self.nadir_time = np.concatenate((self.nadir_time, pixc_edge_reader.getVarValue("nadir_time")))
                self.nadir_time_tai = np.concatenate((self.nadir_time, pixc_edge_reader.getVarValue("nadir_time_tai")))
                self.nadir_longitude = np.concatenate((self.nadir_longitude, pixc_edge_reader.getVarValue("nadir_longitude")))
                self.nadir_latitude = np.concatenate((self.nadir_latitude, pixc_edge_reader.getVarValue("nadir_latitude")))
                self.nadir_x = np.concatenate((self.nadir_x, pixc_edge_reader.getVarValue("nadir_x")))
                self.nadir_y = np.concatenate((self.nadir_y, pixc_edge_reader.getVarValue("nadir_y")))
                self.nadir_z = np.concatenate((self.nadir_z, pixc_edge_reader.getVarValue("nadir_z")))
                self.nadir_vx = np.concatenate((self.nadir_vx, pixc_edge_reader.getVarValue("nadir_vx")))
                self.nadir_vy = np.concatenate((self.nadir_vy, pixc_edge_reader.getVarValue("nadir_vy")))
                self.nadir_vz = np.concatenate((self.nadir_vz, pixc_edge_reader.getVarValue("nadir_vz")))
                self.nadir_sc_event_flag = np.concatenate((self.nadir_sc_event_flag, pixc_edge_reader.getVarValue("nadir_sc_event_flag")))
                self.nadir_tvp_qual = np.concatenate((self.nadir_tvp_qual, pixc_edge_reader.getVarValue("nadir_tvp_qual")))
    
                # 5.4 - Update metadata from PixC info
                self.pixc_metadata["start_time"] = min(self.pixc_metadata["start_time"], str(pixc_edge_reader.getAttValue("start_time")))
                self.pixc_metadata["stop_time"] = max(self.pixc_metadata["stop_time"], str(pixc_edge_reader.getAttValue("stop_time")))
    
    
            # 6 - Close file
            pixc_edge_reader.close()

            retour = out_nb_pix, out_tile_number, out_nb_pix_azimuth

        return retour
        
    # ----------------------------------------

    def swath_global_relabeling(self):
        """
        This function gives new labels to entities gathered at tile edges.
        """
        logger = logging.getLogger(self.__class__.__name__)
        # 1 - Deal with all edges

        for i_edge, tile_num in enumerate(self.tile_num[:-1]):  # Loop over tile edges

            logger.info(" ***** Processing edge of tiles %d and %d *****" % (tile_num, tile_num + 1))

            # 1.1 - Get indices of pixels processed at the current edge
            tile_idx1 = np.where(self.tile_index == tile_num)[0]
            tile_idx2 = np.where(self.tile_index == tile_num + 1)[0]

            # 1.2 - If one tile does not have pixel to process, continue to next iteration
            if (tile_idx1.size == 0) or (tile_idx2.size == 0):
                logger.info("")
                continue
            logger.info(" %d pixels of tile %d will be match with %d pixels of tile %d" % (tile_idx1.size, tile_num, tile_idx2.size, tile_num + 1))

            # 1.3 - New labels for pixels at the edge of current edge
            new_labels_subset = self.gather_edge_entities(tile_idx1, tile_idx2)
            
            # 1.4 - Link old labels to new ones
            self.label_matching(tile_idx1, tile_idx2, new_labels_subset)
            
            logger.info("")

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

    # ----------------------------------------

    def gather_edge_entities(self, in_tile_idx1, in_tile_idx2):
        """
        This function gives new labels for pixels at current tile edge.
            1. Pixels within a buffer around tile edge are selected.
            2. A water mask is computed
            3. The mask is labeled
            
        :param in_tile_idx1: indices of pixels edge of tile 1
        :type in_tile_idx1: 1D-array of int
        :param in_tile_idx2: indices of pixels edge of tile 2
        :type in_tile_idx2: 1D-array of int
        
        :return: out_new_labels_subset = new labels given to pixels of edge entities
        :rtype: 1D-array of int
        """
        logger = logging.getLogger(self.__class__.__name__)

        # 1 - Pixels at top / bottom of tile are loaded
        rg, az = self.select_edge_pixels(in_tile_idx1, in_tile_idx2)

        # 2 - Equivalent matrix size in azimuth and range
        # 2.1 - Equivalent matrix size in azimuth and range
        nb_pix_range = max(rg) + 1
        nb_pix_azimuth = max(az) + 1
        # 2.2 - Compute water mask over the subset
        water_mask = my_tools.computeBinMat(nb_pix_range, nb_pix_azimuth, rg, az)

        # 3 - Label entities over the subset of PixC at the edge tile
        sep_entities, nb_obj = my_tools.labelRegion(water_mask)

        # 4 - Convert into 1D-array
        out_new_labels_subset = my_tools.convert2dMatIn1dVec(rg, az, sep_entities)

        logger.info("%d separate entities located at the edge tile" % np.unique(out_new_labels_subset).size)

        return out_new_labels_subset

    # ----------------------------------------

    def select_edge_pixels(self, in_tile_idx1, in_tile_idx2):
        """
        This function selects pixels at top and bottom of tile 1 and 2
        
        :param in_tile_idx1: indices of pixels edge of tile 1
        :type in_tile_idx1: 1D-array of int
        :param in_tile_idx2: indices of pixels edge of tile 2
        :type in_tile_idx2: 1D-array of int
        
        :return: out_rg = range indices of pixels at the edge
        :rtype: 1D-array of int
        :return: out_az = azimuth indices of pixels at the edge
        :rtype: 1D-array of int
        """

        # 1 - Distinguish top / bottom edge
        rg1, az1 = self.get_edge_pixels(in_tile_idx1, "top")
        rg2, az2 = self.get_edge_pixels(in_tile_idx2, "bottom")

        # 2 - Concatenate pixels range in a numpy array
        rg = np.concatenate((rg1, rg2))

        # 3 - Concatenate pixels azimuth in a numpy array
        az = np.concatenate((az1, az2 + max(az1)))


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
        This function matches labels computed in LakeTile_edge file with labels computed in the gatherEdgeEntities() function.
        Pixels belonging to the same entity but cut by the tiling process are gathered with a single new label.
            
        :param in_tile_idx1: indices of pixels in tile 1
        :type in_tile_idx1: 1D-array of int
        :param in_tile_idx2: indices of pixels in tile 2
        :type in_tile_idx2: 1D-array of int
        :param in_new_labels_subset: new labels computed in gatherEdgeEntities()
        :type in_new_labels_subset: 1D-array of int
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Matching old labels from LakeTile_edge with new labels")

        # 1 - Get old labels computed in LakeTile_edge
        old_labels_subset1, old_labels_subset2 = self.select_edge_labels(in_tile_idx1, in_tile_idx2)

        # 2 - Get new labels
        # NB: the first part of in_new_labels_subset contains labels of tile 1, the end of the array contains labels of tile 2
        in_new_labels_subset1 = in_new_labels_subset[:old_labels_subset1.size]
        in_new_labels_subset2 = in_new_labels_subset[old_labels_subset1.size:]

        # the structure of in_new_labels_subset{1-2} and old_labels_subset{1-2} is the exactly the same, it contains the labels of all pixels within the subset defined the azimuth buffer

        in_new_labels_subset_unique = np.unique(in_new_labels_subset)

        # correspondance contains a list of tuples. Each tuple will correspond to a new entity and will contains the old labels of tiles 1 and 2
        correspondance = []

        # Matching the labels of the subset of tile 1 and 2
        for new_l in in_new_labels_subset_unique:
            # get old labels (lake tile labels) of tiles 1 et 2
            old_l1 = np.unique(old_labels_subset1[np.where(in_new_labels_subset1 == new_l)])
            old_l2 = np.unique(old_labels_subset2[np.where(in_new_labels_subset2 == new_l)])

            # Most current case : at one label of tile1 corresponds one label of tile 2
            if old_l1.size == 1 and old_l2.size == 1:
                # appending to correspondance a tuple with old label 1 and old label 2
                correspondance.append((old_l1[0], old_l2[0]))

            # At one label of tile 1 do not corresponds a label of tile 2. The case happens when the lake is entierly located at the boundary of on tile bot does not cross the border.
            elif old_l2.size == 0:
                for idx1 in np.arange(old_l1.size):
                    correspondance.append((None, old_l1[idx1]))

            # At one label of tile 2 do not corresponds a label of tile 1.
            elif old_l1.size == 0:
                for idx2 in np.arange(old_l2.size):
                    correspondance.append((None, old_l2[idx2]))

            # Case that rarely occurs : the lake meanders along the border between two tiles, in this case severals labels from tile1 matches with several labels of tile2.
            else:
                for idx1 in np.arange(old_l1.size):
                    for idx2 in np.arange(old_l2.size):
                        correspondance.append((old_l1[idx1], old_l2[idx2]))

        # To give an explicit example, labels of tile1 of replaced by letters belonging to ['a', ..] and labels of tile2 are replaces by lettres belonging to [ 'r', ...]
        # correspondance : [(a, r), (a, s), (b, s), (b, t), (c, u), (d, u)]
        # labels a, b and r, s, t belong to the same entity, labels c, d, u belong to a separate entity.

        # The fonction lib.matchLabels returns reorganised correspondance labels :
        # label_matched : [([a, b], [r, s, t]), ([c, d], [u])]
        label_matched = match_labels(correspondance)

        logger.debug(" %d labels of first tile matched with %d labels of second tile into %d entities" % (np.unique(old_labels_subset1).size, np.unique(old_labels_subset2).size, len(label_matched)))

        # need new labels for global lake_sp processings
        unique_label = np.arange(len(label_matched)).astype('int') + max(self.labels) + 1
        # local labels from half edge tiles 1 and 2 are moved into new labels global only for Lake_sp processings
        # Ex : label_matched : [([a, b],[r,s,t]),([c,d],[u])]
        # The for loop iterates over entities at current tile edge
        for idx, label in enumerate(label_matched):
            # l contains the tuple corresponding to the idx^th entity
            # Ex : l : ([a, b], [r, s, t])

            # new_labels[idx] contains a label specefic for lake_SP processings
            new_label = unique_label[idx]

            for old_l1 in label[0]:
                # At one label of tile 2 do not corresponds a label of tile 1. The case happens when the lake is entierly located at the boundary of on tile bot does not cross the border.
                # In this case, old_l1 is setted to None, then, it is not processed.
                if old_l1:
                    # Get label of entity already computed in the case of a lake covering more than two tiles
                    labels_concerned = np.unique(self.labels[in_tile_idx1][np.where(np.logical_and(self.edge_label[in_tile_idx1] == old_l1, self.edge_loc[in_tile_idx1] == 2))])
                    # Deleting label = 0 as those label are not already computed
                    labels_concerned = np.delete(labels_concerned, np.where(labels_concerned == 0))
                    # For each already processed and more than two tile lake, update the global labels
                    for label_to_relabel in labels_concerned:
                        self.labels[np.where(self.labels == label_to_relabel)] = new_label

                    # Set global label
                    self.labels[in_tile_idx1[np.where(self.edge_label[in_tile_idx1] == old_l1)]] = new_label

            for old_l2 in label[1]:
                # At one label of tile 1 do not corresponds a label of tile 2. The case happens when the lake is entierly located at the boundary of on tile bot does not cross the border.
                # In this case, old_l2 is setted to None, then, it is not processed.
                if old_l2:
                    # Get label of entity already computed in the case of a lake covering more than two tiles
                    labels_concerned = (np.unique(self.labels[in_tile_idx2][np.where(np.logical_and(self.edge_label[in_tile_idx2] == old_l2, self.edge_loc[in_tile_idx2] == 2))]))
                    # Deleting label = 0 as those label are not already computed
                    labels_concerned = np.delete(labels_concerned, np.where(labels_concerned == 0))

                    # For each already processed and more than two tile lake, update the global labels
                    for label_to_relabel in labels_concerned:
                        self.labels[np.where(self.labels == label_to_relabel)] = new_label

                    # Set global label
                    self.labels[in_tile_idx2[np.where(self.edge_label[in_tile_idx2] == old_l2)]] = new_label

        nb_edge_entities = np.unique(np.concatenate((self.labels[in_tile_idx1], self.labels[in_tile_idx2]))).size

        logger.debug(" %d separate entities are located at tile edge" % nb_edge_entities)

    # ----------------------------------------

    def select_edge_labels(self, in_tile_idx1, in_tile_idx2):
        """
        This function selects old labels from LakeTile_edge at top and bottom of tile 1 and 2

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
        This function returns the LakeTile_edge labels of pixels within the buffer zone
        
        :param in_tile_idx: indices of pixels at the edge of tile
        :type in_tile_idx: 1D array of int
        :param in_edge_loc_str: edge location = "top" or "bottom"
        :type in_edge_loc_str: string
        
        :return: LakeTile_edge labels of pixels within the buffer zone
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
            lake_azimuth_idx[np.where(lake_tile_idx > tile)] = lake_azimuth_idx[np.where(lake_tile_idx > tile)] + max(self.azimuth_index[in_indices][np.where(lake_tile_idx == tile)]) + 1

        return lake_azimuth_idx

    # ----------------------------------------

    def get_majority_pixels_tile_ref(self, in_label):
        """
        This fuction returns the tile reference of the tile containing the larger number of pixels with the given label.
            
        :param in_label : labels of lake to process
        :type in_label: int
        
        :return: tile reference
        :rtype: string
        """
        
        # 1 - Get unique values and counts of tiles for pixels with label in_label
        # unique, counts = np.unique(self.tile_idx[np.where(self.labels == in_label)], return_counts=True)
        unique, counts = np.unique(self.tile_index[np.where(self.labels == in_label)], return_counts = True)

        # 2 - Get tile ref corresponding to the max number of pixels
        OUT_tile_max_pix = str(unique[np.where(counts == max(counts))][0]) + self.swath_side

        return OUT_tile_max_pix
        
    # ----------------------------------------

    def get_lake_tile_label(self, in_new_label):
        """
        This function is designed to retrieve old labels of PGE_LakeTile. The given new label corresponds to a global label, corresponding to several old labels.
        The old label involving the largest number of pixels is return.
            
        :param in_new_label: global new label
        :type in_new_label: int
        
        :return: LakeTile_edge label involving the largest number of pixels
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



#######################################


def match_labels(in_liste):
    """
    This function reorganise labels in order to group labels by entities
        
    :param IN_liste: ex : [(a, r), (a, s), (b, s), (b, t), (c, u), (d, u)]. Labels a, b and r, s, t belong to the same entity, labels c, d, u belong to a separate entity.
    
    :return: labels gathered by entities
             ex : [(set([a, b]), set([r, s, t])), (set([c, d]), set([u]))] <=> [([a, b], [r, s, t]), ([c, d], [u])]
    """
    return group_by_second(group_by_first(in_liste))


def group_by_first(in_liste):
    """
    This function take a list of tuples. The list is grouped by the first element of tuple and returned as a dictionary.
        
    :param IN_liste: la list of tuple. Ex : [ (a, r), (b, r),  (c, s), (c, t)]
    
    :return out_dico : ex : {a: set([r]), b: set([r]), c : set([t]}
    """

    out_dico = {}
    
    for (ind_i, ind_j) in in_liste:
        if ind_i not in out_dico:
            out_dico[ind_i] = set()
        out_dico[ind_i].add(ind_j)

    return out_dico


def group_by_second(in_dict):
    """
    This function take a dictionary. The dictionary is grouped by second elements and returned as a list of tuple of set.
        
    :param in_dict: result of group by first function: {a: set([r]), b: set([r]), c : set([t]}
    
    :return: a list of tuple of set. ex : [(set([a, b]), set([r])), (set([c]), set([s, t]))]
    """
    
    # Init
    out_results = []
    
    for key, value in list(in_dict.items()):
        
        has_been_found = False
        
        for ind, result in enumerate(out_results):
            if checklist_inter(value, result):
                has_been_found = True
                out_results[ind] = merge(result, key, value)

        if not has_been_found:
            out_results.append((set([key]), value))
            
    return out_results


def checklist_inter(in_value, in_result):
    """
    Check if an element of value is contained in result
        
    :param in_value: a set of elements ex : set([r])
    :param in_result: a tuple of set ex : (set([a]), set([r]))
    
    :return: set([r])
    """
    return in_result[1].intersection(in_value)


def merge(in_result, in_key, in_value):
    """
    Merge value with element of tuple key in result.
        
    :param in_result: a tuple of set ex : (set([a]), set([r]))
    :param in_key: key of first element of tuple. ex : b
    :param in_value: set to merge ex : set([r])
    
    :return: merged tuple of set. ex :set([a, b]), set([r])
    """
    return (in_result[0].union(set([in_key])), in_result[1].union(in_value))
