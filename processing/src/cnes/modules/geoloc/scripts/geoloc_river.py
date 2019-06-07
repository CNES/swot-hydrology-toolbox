#!/usr/bin/env python
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
'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''

from cnes.modules.geoloc.lib.interface import Sensor, RiverTile, PixcvecRiver
import cnes.modules.geoloc.lib.geoloc as geoloc

from cnes.common.lib.my_variables import GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE

import numpy as np
import argparse
import os
import os.path
from collections import defaultdict
from numpy import polyfit
import shutil
from netCDF4 import Dataset
from scipy import interpolate
import cnes.common.service_error as service_error
import logging

# used to filter out fill value and clearly wrong heights, dead sea is about -410 (for now)
MIN_HEIGHT_VALUE = -800


class GeolocRiver(object):
    """
        class GeolocRiver
    """
    def __init__(self, pixc, pixcvec, sensor, rivertile):
        """
        Constructor of GeolocRiver
          :param pixc: pixel cloud object
          :type pixc: class PixelCloud
          :param pixcvec: PIXCVecRiver object
          :type pixcvec: class PixcvecRiver
          :param sensor: sensor object
          :type sensor: class Sensor
          :param rivertile: rivertile object
          :type class RiverTile

        Variables of the object:
        - pixc / pixel cloud from which to compute lake products
        - pixcvec / pixel cloud complementary from which to compute lake products
        - sensor / associated data for along-track pixels
        - rivertile / associated data for each river pixel of pixel cloud


        """
        # 1 - Initiate logging service
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("GeolocRiver initialization")

        # 2 - Initialization of the GeolocRiver variables
        self.pixc = pixc
        self.pixcvec = pixcvec
        self.sensor = sensor
        self.rivertile = rivertile
        # Build a dict of (reach id, node id) to node index in the node list
        # Used to find node index (in the node file) corresponding to the (reach index, nodex index) tuple from the pixel cloud
        self.node_reach_dict = {}
        for rivertile_index in range(self.rivertile.node_indx.size):
            rivertile_reach_indx = self.rivertile.reach_indx[rivertile_index]
            rivertile_node_indx = self.rivertile.node_indx[rivertile_index]
            self.node_reach_dict[(rivertile_reach_indx, rivertile_node_indx)] = rivertile_index

        # Build a dict of reach_id -> list of nodes ids
        # Used to interpolate node heights per reach
        self.reach_dict = defaultdict(list)
        for rivertile_index in range(self.rivertile.node_indx.size):
            rivertile_reach_indx = self.rivertile.reach_indx[rivertile_index]
            rivertile_node_indx = self.rivertile.node_indx[rivertile_index]
            self.reach_dict[rivertile_reach_indx].append(rivertile_node_indx)

    def fit_node_heights(self):
        """Linear fit on node heights per reach"""

        # dict of reach id -> linear fit coefficients
        reach_fit = {}

        # For each reach, compute the fit
        for reach_id in self.reach_dict.keys():

            # x variable are the node ids. This assumes that nodes are equidistant
            # to be improved using the along reach coordinate of the node instead
            x = self.reach_dict[reach_id]

            # y variable is height
            # get the node heights for the current reach from the reach dict
            y = np.array([self.rivertile.h_n_ave[self.node_reach_dict[(reach_id, node_indx)]] for node_indx in x])

            # filter bad data in y (nan, -inf, fill values)
            good = (y >= MIN_HEIGHT_VALUE)
            x = np.array(x)[good]
            y = np.array(y)[good]

            # Linear fitting (if enough points)
            if np.sum(good) >= 2:
                a, b = polyfit(x, y, 1)
                reach_fit[reach_id] = (a, b)

            else:
                reach_fit[reach_id] = None

        # Reconstruct node heights using the computed coefficients
        for rivertile_index in range(self.rivertile.node_indx.size):
            rivertile_reach_indx = self.rivertile.reach_indx[rivertile_index]
            rivertile_node_indx = self.rivertile.node_indx[rivertile_index]
            if reach_fit[rivertile_reach_indx] is not None:
                a, b = reach_fit[rivertile_reach_indx]
                self.rivertile.h_n_ave_fit[rivertile_index] = a * rivertile_node_indx + b

    def estimate_pixvec_height_for_geoloc(self, interpolate_pixc_between_nodes):
        """
        Compute new heights for each pixcvec point,
        either the associated node height (interpolate_pixc_between_nodes == False),
        or interpolating between the associated node and the next/previous closest (interpolate_pixc_between_nodes == True)

            :param interpolate_pixc_between_nodes: active or deactivate the interpolation
            :type interpolate_pixc_between_nodes: boolean
        """

        self.new_height = np.copy(self.pixcvec.height)

        unfound_keys = 0

        # For each pixcvec point
        for point_index in range(self.pixcvec.height.size):
            reach_index = self.pixcvec.reach_index[point_index]
            node_index = self.pixcvec.node_index[point_index]

            # Check if the reach and node of the point exist in the node file
            key = (reach_index, node_index)
            if key not in self.node_reach_dict:
                unfound_keys += 1
            else:
                rivertile_index = self.node_reach_dict[key]

                # If possible, interpolate new point height from the associated node and the next/previous closest
                if interpolate_pixc_between_nodes and (reach_index, node_index - 1) in self.node_reach_dict and \
                                (reach_index, node_index + 1) in self.node_reach_dict:
                    self.interpolate_pixc_height(reach_index,node_index,point_index,rivertile_index)

                # Simple case: just use the node height to estimate improved_height
                else:
                    new_height = self.rivertile.h_n_ave[rivertile_index]
                    # avoid wrong data
                    if new_height >= MIN_HEIGHT_VALUE:
                        self.new_height[point_index] = new_height

        if unfound_keys > 0:
            #  - Initiate logging service
            logger = logging.getLogger(self.__class__.__name__)
            logger.warning("Warning: {} points' (reach, node) were not found in the node file".format(unfound_keys))

    def interpolate_pixc_height(self,reach_index, node_index,point_index, rivertile_index):


        rivertile_index_less_one = self.node_reach_dict[(reach_index, node_index - 1)]
        rivertile_index_plus_one = self.node_reach_dict[(reach_index, node_index + 1)]

        # 's' means curvilinear coordinate (distance along reach)
        s_pixc = self.pixcvec.along_reach[point_index]
        s_node = self.rivertile.s[rivertile_index]
        s_node_previous = self.rivertile.s[rivertile_index_less_one]
        s_node_after = self.rivertile.s[rivertile_index_plus_one]

        # logic to handle if point is before/after the associated node
        if s_node < s_pixc:
            d1 = s_pixc - s_node
            d2 = s_node_after - s_pixc
            h1 = self.rivertile.h_n_ave_fit[rivertile_index]
            h2 = self.rivertile.h_n_ave_fit[rivertile_index_plus_one]
        else:
            d1 = s_pixc - s_node_previous
            d2 = s_node - s_pixc
            h1 = self.rivertile.h_n_ave_fit[rivertile_index_less_one]
            h2 = self.rivertile.h_n_ave_fit[rivertile_index]

        # Linear interpolation of height
        if h1 >= MIN_HEIGHT_VALUE and h2 >= MIN_HEIGHT_VALUE:
            self.new_height[point_index] = h1 + (h2 - h1) / (d1 + d2) * d1


    def apply_improved_geoloc(self, method='taylor'):
        """
        Compute the new lat, lon, height using the new heights

        :param method: using method to compute the new lat, lon, height
        :type method: string
        """
        if method == 'taylor':
            self.taylor_improved_geoloc()
        else:
            message = "the methode " + str(method) + " is undefined"
            raise service_error.ParameterError("apply_improved_geoloc", message)

    def taylor_improved_geoloc(self):
        """
        Improve the height of noisy point (in object sensor)
        """

        nb_pix = self.pixcvec.height.size
        # Convert geodetic coordinates (lat, lon, height) to cartesian coordinates (x, y, z)
        x, y, z = geoloc.convert_llh2ecef(self.pixcvec.latitude, self.pixcvec.longitude, self.pixcvec.height, GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE)

        # Get position of associated along-track pixels (in cartesian coordinates)
        nadir_x = self.sensor.nadir_x
        nadir_y = self.sensor.nadir_y
        nadir_z = self.sensor.nadir_z

        # Get velocity of associated along-track pixels (in cartesian coordinates)
        nadir_vx = self.sensor.nadir_vx
        nadir_vy = self.sensor.nadir_vy
        nadir_vz = self.sensor.nadir_vz

        # Get distance from satellite to target point
        ri = self.pixcvec.near_range + self.pixcvec.range_idx * self.pixcvec.range_spacing

        # Init output vectors
        self.out_lat_corr = np.zeros(nb_pix)  # Improved latitudes
        self.out_lon_corr = np.zeros(nb_pix)  # Improved longitudes
        self.out_height_corr = np.zeros(nb_pix)  # Improved heights

        # need to remap illumnation time to nearest sensor index
        # TODO replace this by a call to a get_sensor_index or equivalent function
        # that either interpolates the sensor or does something more efficient
        f = interpolate.interp1d(self.sensor.time, range(len(self.sensor.time)))
        sensor_s = (np.rint(f(self.pixc.illumination_time))).astype(int)

        # Loop over each pixel (could be vectorized)
        # vectorisation
        h_noisy = self.pixcvec.height
        nadir_x_vect = np.zeros(nb_pix)
        nadir_y_vect = np.zeros(nb_pix)
        nadir_z_vect = np.zeros(nb_pix)
        nadir_vx_vect = np.zeros(nb_pix)
        nadir_vy_vect = np.zeros(nb_pix)
        nadir_vz_vect = np.zeros(nb_pix)

        for i in np.arange(nb_pix):
            ind_sensor = sensor_s[i]
            nadir_x_vect[i] = nadir_x[ind_sensor]
            nadir_y_vect[i] = nadir_y[ind_sensor]
            nadir_z_vect[i] = nadir_z[ind_sensor]
            nadir_vx_vect[i] = nadir_vx[ind_sensor]
            nadir_vy_vect[i] = nadir_vy[ind_sensor]
            nadir_vz_vect[i] = nadir_vz[ind_sensor]

        # improuve height with vectorised pixel
        p_final, p_final_llh, h_mu, (iter_grad, nfev_minimize_scalar) = geoloc.pointcloud_height_geoloc_vect(np.transpose(np.array([x, y, z])),
                                                                                                             h_noisy,
                                                                                                             np.transpose(np.array(
                                                                                                                 [nadir_x_vect, nadir_y_vect,
                                                                                                                  nadir_z_vect])),
                                                                                                             np.transpose(np.array(
                                                                                                                 [nadir_vx_vect, nadir_vy_vect,
                                                                                                                  nadir_vz_vect])),
                                                                                                             ri, self.new_height,
                                                                                                             recompute_doppler=True,
                                                                                                             recompute_range=True, verbose=False,
                                                                                                             max_iter_grad=1, height_goal=1.e-3)

        self.out_lat_corr, self.out_lon_corr, self.out_height_corr = np.transpose(p_final_llh)


def geoloc_river(pixc, pixcvec, sensor, rivertile, fit_heights_per_reach, interpolate_pixc_between_nodes, method='taylor'):
    """
    Improved river geolocation

      :param pixc: pixel cloud object
      :type pixc: class PixelCloud
      :param pixcvec: PIXCVecRiver object
      :type pixcvec: class PixcvecRiver
      :param sensor: sensor object
      :type sensor: class Sensor
      :param rivertile: rivertile object
      :type rivertile:class RiverTile
      :param fit_heights_per_reach: active or deactivate the linear fit on node heights per reach
      :type fit_heights_per_reach: boolean
      :param interpolate_pixc_between_nodes: active or deactivate the interpolation
      :type interpolate_pixc_between_nodes: boolean
      :param method: using method to compute the new latitude, longitude and height
      :type method: string
    """
    geolocriver = GeolocRiver(pixc, pixcvec, sensor, rivertile)
    # Do the improved river geolocation
    # 1 - Initiate logging service
    logger = logging.getLogger()
    logger.info("Improved geolocation")
    if fit_heights_per_reach:
        geolocriver.fit_node_heights()

    geolocriver.estimate_pixvec_height_for_geoloc(interpolate_pixc_between_nodes)
    geolocriver.apply_improved_geoloc(method=method)

    return geolocriver.out_lat_corr, geolocriver.out_lon_corr, geolocriver.out_height_corr
