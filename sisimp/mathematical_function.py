#!/usr/bin/env python
"""
module mathematical_function.py

module author : Capgemini

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from math import sin

import os
import re

import lib.my_api as my_api

from lib.my_variables import RAD2DEG, DEG2RAD, GEN_APPROX_RAD_EARTH


def calc_delta_h(IN_angles, IN_noise_height, IN_height_bias_std, seed=None):
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
    
    np.random.seed(seed)

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
        OUT_noisy_h = np.random.normal(0, IN_height_bias_std) + np.random.normal(0, np.interp(IN_angles*RAD2DEG, IN_noise_height[:, 0], IN_noise_height[:, 1]))

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


def lonlat_from_azy(IN_az, IN_ri, IN_attributes, IN_swath, h=0, IN_unit="rad"):
    """
    Convert coordinates from azimuth-y to lon-lat for a given track

    :param IN_az: azimuth coordinate of given points
    :type IN_az: 1D-array of float
    :param IN_ri: crosstrack distance of given points
    :type IN_ri: 1D-array of float
    :param IN_attributes:
    :type IN_attributes:
    :param IN_swath:
    :type IN_swath:
    :param h:
    :type h:
    :param IN_unit: "rad" (default) ou "deg" to output coordinates in radians or degrees
    :type IN_unit: string

    :return: OUT_lon = longitude of points
    :rtype: OUT_lon = 1D-array of float
    :return: OUT_lat = latitude of points
    :rtype: OUT_lat = 1D-array of float
    """

    lat0 = IN_attributes.lat[IN_az]
    phi0 = IN_attributes.lon[IN_az]
    psi0 = IN_attributes.heading_init[IN_az]
    IN_Alt = IN_attributes.alt[IN_az]
    theta0 = np.pi/2 - lat0
    
    costheta_0 = IN_attributes.costheta_init[IN_az]
    sintheta_0 = IN_attributes.sintheta_init[IN_az]
    cosphi_0 = IN_attributes.cosphi_init[IN_az]
    sinphi_0 = IN_attributes.sinphi_init[IN_az]
    cospsi_0 = IN_attributes.cospsi_init[IN_az]
    sinpsi_0 = IN_attributes.sinpsi_init[IN_az]
    
    if IN_swath == 'Left':
        sign = -1
    else:
        sign = 1.

    mu = sign*np.arccos((IN_ri**2-((GEN_APPROX_RAD_EARTH+h)**2+(GEN_APPROX_RAD_EARTH+IN_Alt)**2))/(-2*(GEN_APPROX_RAD_EARTH + h)*(GEN_APPROX_RAD_EARTH+IN_Alt)))

    Cx = (np.cos(mu)*sintheta_0*cosphi_0 + np.sin(mu)*(sinpsi_0*costheta_0*cosphi_0-cospsi_0*sinphi_0))
    Cy = (np.cos(mu)*sintheta_0*sinphi_0 + np.sin(mu)*(sinpsi_0*costheta_0*sinphi_0+cospsi_0*cosphi_0))
    Cz = (np.cos(mu)*costheta_0 - np.sin(mu)*sinpsi_0*sintheta_0)
    
    xp_yp_zp = np.where(np.logical_and(Cx > 0, np.logical_and(Cy > 0, Cz > 0)))
    xp_yp_zm = np.where(np.logical_and(Cx > 0, np.logical_and(Cy > 0, Cz < 0)))
    xp_ym_zp = np.where(np.logical_and(Cx > 0, np.logical_and(Cy < 0, Cz > 0)))
    xm_yp_zp = np.where(np.logical_and(Cx < 0, np.logical_and(Cy > 0, Cz > 0)))
    xp_ym_zm = np.where(np.logical_and(Cx > 0, np.logical_and(Cy < 0, Cz < 0)))
    xm_ym_zp = np.where(np.logical_and(Cx < 0, np.logical_and(Cy < 0, Cz > 0)))
    xm_yp_zm = np.where(np.logical_and(Cx < 0, np.logical_and(Cy > 0, Cz < 0)))
    xm_ym_zm = np.where(np.logical_and(Cx < 0, np.logical_and(Cy < 0, Cz < 0)))

    OUT_lat = np.zeros(len(mu), float)
    OUT_lon = np.zeros(len(mu), float)
     
    OUT_lat[xp_yp_zp] = np.pi/2. - np.arccos(Cz[xp_yp_zp])
    OUT_lon[xp_yp_zp] = np.arctan(Cy[xp_yp_zp]/Cx[xp_yp_zp])

    OUT_lat[xp_yp_zm] = np.pi/2. - np.arccos(Cz[xp_yp_zm])
    OUT_lon[xp_yp_zm] = np.arctan(Cy[xp_yp_zm]/Cx[xp_yp_zm]) 
    
    OUT_lat[xp_ym_zp] = np.pi/2. - np.arccos(Cz[xp_ym_zp])
    OUT_lon[xp_ym_zp] = np.arctan(Cy[xp_ym_zp]/Cx[xp_ym_zp])
    
    OUT_lat[xm_yp_zp] = np.pi/2. - np.arccos(Cz[xm_yp_zp])
    OUT_lon[xm_yp_zp] = np.arctan(Cy[xm_yp_zp]/Cx[xm_yp_zp]) + np.pi
    
    OUT_lat[xp_ym_zm] = np.pi/2. - np.arccos(Cz[xp_ym_zm])
    OUT_lon[xp_ym_zm] = np.arctan(Cy[xp_ym_zm]/Cx[xp_ym_zm])
    
    OUT_lat[xm_ym_zp] = np.pi/2. - np.arccos(Cz[xm_ym_zp])
    OUT_lon[xm_ym_zp] = np.arctan(Cy[xm_ym_zp]/Cx[xm_ym_zp]) - np.pi
    
    OUT_lat[xm_yp_zm] = np.pi/2. - np.arccos(Cz[xm_yp_zm])
    OUT_lon[xm_yp_zm] = np.arctan(Cy[xm_yp_zm]/Cx[xm_yp_zm]) + np.pi
    
    OUT_lat[xm_ym_zm] = np.pi/2. - np.arccos(Cz[xm_ym_zm])
    OUT_lon[xm_ym_zm] = np.arctan(Cy[xm_ym_zm]/Cx[xm_ym_zm]) - np.pi/2.

    if IN_unit == "deg":
        return OUT_lon*RAD2DEG, OUT_lat*RAD2DEG  # Output in degrees
    return OUT_lon, OUT_lat  # Output in radians


def lonlat_from_azy_old(IN_az, IN_y, IN_lat_init, IN_lon_init, IN_heading_init, IN_unit="rad"):
    """
    Convert coordinates from azimuth-y to lon-lat for a given track

    :param IN_az: azimuth coordinate of given points
    :type IN_az: 1D-array of float
    :param IN_y: crosstrack distance of given points
    :type IN_y: 1D-array of float
    :param IN_lat_init:
    :type IN_lat_init:
    :param IN_lon_init:
    :type IN_lon_init:
    :param IN_heading_init:
    :type IN_heading_init:
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
        OUT_y = np.where(IN_x < IN_xp[0], IN_yp[0]+(IN_x-IN_xp[0])*(IN_yp[0]-IN_yp[1])/(IN_xp[0]-IN_xp[1]), OUT_y)
        OUT_y = np.where(IN_x > IN_xp[-1], IN_yp[-1]+(IN_x-IN_xp[-1])*(IN_yp[-1]-IN_yp[-2])/(IN_xp[-1]-IN_xp[-2]), OUT_y)
    else:
        OUT_y = np.interp(IN_x, IN_xp[::-1], IN_yp[::-1])
        OUT_y = np.where(IN_x > IN_xp[0], IN_yp[0], OUT_y)
        OUT_y = np.where(IN_x < IN_xp[-1], IN_yp[-1], OUT_y)
    return OUT_y
