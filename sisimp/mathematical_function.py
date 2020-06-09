#!/usr/bin/env python
"""
module mathematical_function.py

module author : Capgemini

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""

from __future__ import absolute_import, division, print_function, unicode_literals
from scipy import ndimage

import numpy as np
from math import sin

import os
import re
import sys

import lib.my_api as my_api
import pygeodesy

from lib.my_variables import RAD2DEG, DEG2RAD, GEN_APPROX_RAD_EARTH


def calc_delta_h(IN_water_pixels, IN_angles_water, IN_angles, IN_noise_height, IN_height_bias_std, sensor_wavelength, baseline, near_range, seed=None):
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

    ## Conversion temporaire du fichier erreur de hauteur en fichier erreur de phase
    IN_noise_phase = IN_noise_height[:,1]*2*np.pi/(sensor_wavelength*near_range*np.sin(IN_noise_height[:,0]*DEG2RAD)/baseline)
    ## Conversion temporaire du fichier erreur de hauteur en fichier erreur de phase



    if (IN_angles.size != 0) and (np.max(IN_angles*RAD2DEG) > np.max(IN_noise_height[:, 0])):
        my_api.printInfo("One or more incidence angles are greater than the max value defined in the noise file ! Values higher than {0} degrees will be set to      the maximum value defined in the file.".format(np.max(IN_noise_height[:, 0])))
        stdv = np.interp(IN_angles*RAD2DEG, IN_noise_height[:, 0], IN_noise_phase)
        stdv[np.isnan(stdv)]= 0.
        stdv_2d = np.interp(IN_angles_water*RAD2DEG, IN_noise_height[:,0], IN_noise_phase)
        stdv_2d[np.isnan(stdv_2d)]=0.
        OUT_noisy_phi = np.random.normal(0, stdv_2d, (len(IN_water_pixels), len(IN_water_pixels[0])))   
        OUT_noisy_phi[np.where(IN_angles*RAD2DEG) > np.max(IN_noise_height[:, 1])] = IN_noise_height[-1, 1]
        
    if (IN_noise_height[:, 1] < 1.e-5).any() :  # Case noise file as one or more zeros
        OUT_noisy_phi = np.zeros((len(IN_water_pixels), len(IN_water_pixels[0])))
        stdv = np.zeros((len(IN_angles)))
    else:
        stdv = np.interp(IN_angles*RAD2DEG, IN_noise_height[:, 0], IN_noise_phase)
        stdv[np.isnan(stdv)]= 0.
        stdv_2d = np.interp(IN_angles_water*RAD2DEG, IN_noise_height[:,0], IN_noise_phase)
        stdv_2d[np.isnan(stdv_2d)]=0.
        OUT_noisy_phi = np.random.normal(0, stdv_2d, (len(IN_water_pixels), len(IN_water_pixels[0])))

        
    h_amb = sensor_wavelength*near_range*np.sin(IN_angles_water)/baseline
    OUT_noisy_h = np.angle(np.exp(1j*OUT_noisy_phi)) * h_amb/(2*np.pi)
    phase_noise_std = stdv
    
    if IN_height_bias_std != 0.:
        OUT_noisy_h += np.random.normal(0, IN_height_bias_std) 
    
    return OUT_noisy_h, phase_noise_std, sensor_wavelength*near_range*np.sin(IN_angles)/baseline/2/np.pi

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

    my_api.printInfo("[mathematical_function] [calc_delta_jitter] orbit_jitter = %.6f" % orbit_jitter)

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


    cosmu = (IN_ri**2-((GEN_APPROX_RAD_EARTH+h)**2+(GEN_APPROX_RAD_EARTH+IN_Alt)**2))/(-2*(GEN_APPROX_RAD_EARTH + h)*(GEN_APPROX_RAD_EARTH+IN_Alt))
    
    cosmu[np.where(cosmu >1)] = 1.
    cosmu[np.where(cosmu <-1)] = -1.
    
    mu = sign*np.arccos(cosmu)


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

def calc_sigma0(IN_water_pixels, classification_tab,  water_flag, water_land_flag, darkwater_flag, \
                sigma_water_in_db = 10., sigma_darkwater_in_db = 0., sigma_land_in_db = 0., az_ne_look = 4.):
    """
    Compute sigma0 for different pixels (land, water, darkwater) coordinates from azimuth-y to lon-lat for a given track

    :param IN_water_pixels: az/rg array of pixel cloud pixels
    :type IN_az: 2D-array of float
    
    :param classification_tab: classification of givens points (1&2 land, 3&4 water, 24 darkwater)
    :type IN_az: 1D-array of float
    
    :return: sigma0 = sigma0 of points
    :rtype: sigma0 = 1D-array of float 
    """
    
    #Convert sig0 in linear
    sig0_water = 10.**(sigma_water_in_db/10.)
    sig0_darkwater = 10.**(sigma_darkwater_in_db/10.)
    sig0_land = 10.**(sigma_land_in_db/10.)
    
    sig0_water_tab = np.random.normal(sig0_water, sig0_water/np.sqrt(az_ne_look), (len(IN_water_pixels), len(IN_water_pixels[0])))
    sig0_darkwater_tab = np.random.normal(sig0_darkwater, sig0_darkwater/np.sqrt(az_ne_look), (len(IN_water_pixels), len(IN_water_pixels[0])))
    sig0_land_tab = np.random.normal(sig0_land, sig0_land/np.sqrt(az_ne_look), (len(IN_water_pixels), len(IN_water_pixels[0])))
        
    sig0_water_tab_conv = ndimage.convolve(sig0_water_tab, np.array([[1/9, 1/9, 1/9], [1/9, 1/9, 1/9], [1/9, 1/9, 1/9]]))
    sig0_darkwater_tab_conv = ndimage.convolve(sig0_darkwater_tab, np.array([[1/9, 1/9, 1/9], [1/9, 1/9, 1/9], [1/9, 1/9, 1/9]]))
    sig0_land_tab_conv = ndimage.convolve(sig0_land_tab, np.array([[1/9, 1/9, 1/9], [1/9, 1/9, 1/9], [1/9, 1/9, 1/9]]))
    
    sig0_tab = np.copy(sig0_land_tab)
    ind=np.where(IN_water_pixels!=0)
    sig0_tab_flat = sig0_tab[ind]
    
    ind_dw = np.where(classification_tab == darkwater_flag) 
    ind_water = np.where(classification_tab == water_flag) 
    ind_land_water = np.where(classification_tab == water_land_flag) 
    
    sig0_tab_flat[ind_dw] = sig0_darkwater_tab_conv[ind][ind_dw]
    sig0_tab_flat[ind_water] = sig0_water_tab_conv[ind][ind_water]
    sig0_tab_flat[ind_land_water] = sig0_water_tab_conv[ind][ind_land_water]

    return sig0_tab_flat

def calc_geoid(lat, lon, geoid_file):
    
    geoid = pygeodesy.GeoidPGM(geoid_file)
    h_geoid= geoid.height(lat, lon)
    h_geoid += (1.0 + 0.3)*np.sqrt(5/(4*np.pi))*(-0.31466)*(1.5*np.sin(lat*DEG2RAD)*np.sin(lat*DEG2RAD)-0.5)
    
    return h_geoid
