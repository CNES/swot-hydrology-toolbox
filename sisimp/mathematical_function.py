#!/usr/bin/env python
"""
module mathematical_function.py

module author : Capgemini

Copyright (c) 2018 CNES. All rights reserved
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from math import sin

import os
import re

import lib.my_api as my_api

from lib.my_variables import RAD2DEG, DEG2RAD, GEN_APPROX_RAD_EARTH

def calc_delta_h(IN_angles, IN_noise_height, IN_height_bias_std):
    """
    Calculate the delta h values and add noise

    :param IN_angles : the angles
    :type IN_angles: 1D-array of int

    :param IN_noise_height : 
    :type IN_noise_height :

    :param IN_height_bias_std :
    :type IN_height_bias_std :

    :return OUT_noisy_h: the noisy height values
    :rtype OUT_noisy_h: 1D-array of float
    """

    if (IN_angles.size != 0) and (np.max(IN_angles*RAD2DEG) > np.max(IN_noise_height[:, 0])):
        my_api.printInfo("One or more incidence angles are greater than the max value defined in the noise file ! Values higher than {0} degrees will be set to      the maximum value defined in the file.".format(np.max(IN_noise_height[:, 0])))

    OUT_noisy_h = 0
    if (IN_noise_height[:, 1] < 1.e-5).any() and not IN_height_bias_std < 1.e-5:  # Case noise file as one or more zeros
        OUT_noisy_h = np.random.normal(0, IN_height_bias_std) + np.interp(IN_angles*RAD2DEG, IN_noise_height[:, 0], IN_noise_height[:, 1])
    elif not (IN_noise_height[:, 1] < 1.e-5).any() and IN_height_bias_std < 1.e-5:  # Case height bias equals zero
        OUT_noisy_h = np.random.normal(0, np.interp(IN_angles*RAD2DEG, IN_noise_height[:, 0], IN_noise_height[:, 1]))
    elif (IN_noise_height[:, 1] < 1.e-5).any() and IN_height_bias_std < 1.e-5:  # Case both are equals to zero
        OUT_noisy_h = np.interp(IN_angles*RAD2DEG, IN_noise_height[:, 0], IN_noise_height[:, 1])
    else:  # Case none are equals to zero
        OUT_noisy_h = np.random.normal(0, IN_height_bias_std) + np.random.normal(0, np.interp(IN_angles*RAD2DEG, IN_noise_height[:, 0], IN_noise_height[:     , 1]))

    return OUT_noisy_h


def calc_delta_jitter(IN_orbit_heading, IN_lat, IN_orbit_jitter):
    """
    Calculate the jitter

    :param IN_orbit_heading: the orbit headings
    :type IN_orbit_heading: 1D-array of float

    :param IN_lat: the orbit latitudes
    :type IN_lat: 1D-array of float

    :param IN_orbit_jitter :
    :type IN_orbit_jitter :

    :return: the values of delta jitter
    :rtype: 1D-array of float
    """

    # Random jitter of +/- random_jitter
    orbit_jitter = [0, int(np.random.normal(0, IN_orbit_jitter))][IN_orbit_jitter != 0]
    if orbit_jitter > IN_orbit_jitter:
        orbit_jitter = IN_orbit_jitter
    elif orbit_jitter < -IN_orbit_jitter:
        orbit_jitter = -IN_orbit_jitter

    my_api.printInfo("orbit_jitter = %.6f" % orbit_jitter)

    return (np.cos(IN_orbit_heading) * orbit_jitter) / (GEN_APPROX_RAD_EARTH * np.cos(np.mean(IN_lat)*DEG2RAD))


def calc_delta_sensor(IN_delta_h, IN_orbit_altitudes, IN_y):
    """
    Calculate the noise of the sensor.

    :param IN_delta_h: noise in heights
    :type IN_delta_h: 1D-array of float
    :param IN_y: the cross track distance
    :type IN_y: 1D-array of float
    :param IN_orbit_altitudes: the orbit altitudes
    :type IN_orbit_altitudes: 1D-array of float

    :return: the sensor noise
    :rtype: 1D-array of float
    """
    return IN_delta_h * IN_orbit_altitudes / IN_y


def lonlat_from_azy(IN_az, IN_y, IN_lat_init, IN_lon_init, IN_heading_init, IN_unit="rad"):
    """
    Convert coordinates from azimuth-y to lon-lat for a given track

    :param IN_az: azimuth coordinate of given points
    :type IN_az: 1D-array of float
    :param IN_y: crosstrack distance of given points
    :type IN_y: 1D-array of float
    :param IN_unit: "rad" (default) ou "deg" to output coordinates in radians or degrees
    :type IN_unit: string

    :return: OUT_lon = longitude of points
    :rtype: OUT_lon = 1D-array of float
    :return: OUT_lat = latitude of points
    :rtype: OUT_lat = 1D-array of float
    """
    # Differential coordinates in km
    du = np.cos(linear_extrap(IN_az, np.arange(len(IN_lat_init)), IN_heading_init)) * IN_y
    dv = -np.sin(linear_extrap(IN_az, np.arange(len(IN_lat_init)), IN_heading_init)) * IN_y

    # Spherical coordinates
    OUT_lat = linear_extrap(IN_az, np.arange(len(IN_lat_init)), IN_lat_init) + dv / GEN_APPROX_RAD_EARTH
    OUT_lon = linear_extrap(OUT_lat, IN_lat_init, IN_lon_init) + du / (GEN_APPROX_RAD_EARTH * np.cos(linear_extrap(IN_az, np.arange(len(IN_lat_init)), IN_lat_init)))

    if IN_unit == "deg":
        return OUT_lon*RAD2DEG, OUT_lat*RAD2DEG  # Output in degrees
    return OUT_lon, OUT_lat  # Output in radians


def linear_extrap(IN_x, IN_xp, IN_yp):
   if IN_xp[0] < IN_xp[-1]:
       OUT_y = np.interp(IN_x, IN_xp, IN_yp)
       OUT_y = np.where(IN_x<IN_xp[0], IN_yp[0]+(IN_x-IN_xp[0])*(IN_yp[0]-IN_yp[1])/(IN_xp[0]-IN_xp[1]), OUT_y)
       OUT_y = np.where(IN_x>IN_xp[-1], IN_yp[-1]+(IN_x-IN_xp[-1])*(IN_yp[-1]-IN_yp[-2])/(IN_xp[-1]-IN_xp[-2]), OUT_y)
   else:
       OUT_y = np.interp(IN_x, IN_xp[::-1], IN_yp[::-1])
       OUT_y = np.where(IN_x>IN_xp[0], IN_yp[0], OUT_y)
       OUT_y = np.where(IN_x<IN_xp[-1], IN_yp[-1], OUT_y)
   return OUT_y
