# -*- coding: utf8 -*-
'''
.. module lib.py
    :synopsis: Library with generic functions
    Created on 03/02/2017

.. module author: Claire POTTIER - CNES DSO/SI/TR

Copyright (c) 2017 CNES. All rights reserved.
'''

import numpy as np
import mahotas as mh
from subprocess import check_call

# Global variables
GEN_RAD_EARTH_EQ = 6378137.0 # Radius (in meters) of the Earth model (WGS84 ellipsoid) at the equator
GEN_RAD_EARTH_POLE = 6356752.31425 # Radius (in meters) of the Earth model to the pole 
GEN_APPROX_RAD_EARTH = (2*GEN_RAD_EARTH_EQ + GEN_RAD_EARTH_POLE)/3 # Radius (in meters) of the sphere equivalent to ellipsoid


def computeBinMat(IN_sizeX, IN_sizeY, IN_X, IN_Y):
    '''
    Creates a 2D binary matrix from Y and Y 1D vectors
    i.e. for each ind: mat[IN_X[ind],IN_Y[ind]] = 1
        
    :param IN_sizeX: number of pixels in X dimension
    :type IN_sizeX: int
    :param IN_sizeY: number of pixels in Y dimension
    :type IN_sizeY: int
    :param IN_X: X indices of "1" pixels
    :type IN_X: 1D vector of int
    :param IN_Y: Y indices of "1" pixels
    :type IN_Y: 1D vector of int
        
    :return: 2D matrix with "1" for each (IN_X_i, IN_Y_i) and 0 elsewhere
    :rtype: 2D binary matrix of int 0/1
    '''
    print("[lib] == computeBinMat ==")
    
    # 0 - Deal with exceptions
    # 0.1 - Input vectors size must be the same
    if IN_X.size != IN_Y.size:
        raise ValueError("computeBinMat(IN_X, IN_Y) : IN_X and IN_Y must be the same size ; currently : IN_X = %d and IN_Y = %d" % ( IN_X.size , IN_Y.size ))
    else:
        nb_pts = IN_X.size
    print("> Nb pixels to deal with = %d" % ( nb_pts ))
    # 0.2 - max(X) < IN_sizeX
    if np.max(IN_X) >= IN_sizeX:
        raise ValueError("computeBinMat(IN_X, IN_Y) : elements of IN_X must be less than IN_sizeX")
    # 0.3 - max(X) < IN_sizeX
    if np.max(IN_Y) >= IN_sizeY:
        raise ValueError("computeBinMat(IN_X, IN_Y) : elements of IN_Y must be less than IN_sizeY")
    
    # 1 - Init output binary image
    OUT_binIm = np.zeros((IN_sizeY, IN_sizeX))
    print("> Binary matrix size = (X=%d , Y=%d)" % ( IN_sizeX , IN_sizeY ))
        
    # 2 - Put 1 for every pixels defined by the input vectors
    for ind in range(nb_pts): 
        OUT_binIm[IN_Y[ind], IN_X[ind]] = 1
    
    return OUT_binIm
        
#######################################
    
def labelRegion(IN_binMat):
    '''
    Identifies all separate regions in a 2D binary matrix
    
    :param IN_binMat: 2D binary matrix
    :type IN_binMat: 2D binary matrix of int 0/1
    
    :return: a 2D matrix, in which each pixel has a unique value ranging from 0 to N (number of identified regions)
    :rtype: 2D matrix of int 
    :return: the number of objects
    :rtype: int
    '''
    print("[lib] == labelRegion ==")
    
    OUT_mLabelRegions, OUT_nbObj = mh.label(IN_binMat)
    print("> Number of detected objects = %d" % ( OUT_nbObj )) 
    
    return OUT_mLabelRegions, OUT_nbObj
    
#######################################

def convert2dMatIn1dVec(IN_X, IN_Y, IN_mat):
    '''
    Converts a 2D matrix [ IN_X[ind], IN_Y[ind] ] in a 1D vector [ind]
    
    :param IN_mat: 2D matrix
    :type IN_mat: 2D matrix of values (int, float, ...)
    
    :return: a 1D vector, each value at "ind" being the value of IN_mat[ IN_X[ind], IN_Y[ind] ]
    :rtype: 2D matrix of values (int, float, ...)
    '''
    print("[lib] == convert2dMatIn1dVec ==")
    
    # Number of pixels in IN_X
    nb_pix = IN_X.size
    # Init output vector (same size of input IN_X and IN_Y)   
    OUT_vector = np.zeros(nb_pix)
    # Fill the output vector
    for indp in range(nb_pix):
        OUT_vector[indp] = IN_mat[IN_Y[indp], IN_X[indp]]
    
    return OUT_vector
    
#######################################

def cptDigits(IN_real):
    '''
    Compute the position of the first significant figure of a real number (ex: 110 => 2 ; 0.000156 => -4)
        
    :param IN_real: real number
    :type IN_real: float
    
    :return the position of the first significant figure of a real number
    :rtype: int 
    '''
    
    # Compteur de decimales
    if abs(IN_real) > 1:
        OUT_cpt = 1
    else:
        OUT_cpt = 0
    
    # Cas des nombres >= 10
    while abs(IN_real) >= 10:
        IN_real /= 10.0
        OUT_cpt += 1
    
    # Cas des nombres <= 1
    while abs(IN_real) <= 1:
        IN_real *= 10.0
        OUT_cpt -= 1
    
    return OUT_cpt
    
#######################################

def deg2rad(IN_deg):
    '''
    Convert angles in degrees to radians
    
    :param IN_deg: angles to convert
    :type IN_deg: scalar or 1D-array
    
    :return angles in radians
    :rtype: same as input
    '''
    return IN_deg * np.pi / 180.

def rad2deg(IN_rad):
    '''
    Convert angles in radians to degrees
    
    :param IN_rad: angles to convert
    :type IN_rad: scalar or 1D-array
    
    :return angles in degrees
    :rtype: same as input
    '''
    return IN_rad * 180. / np.pi
    
#######################################

def computeAz(IN_lon, IN_lat, IN_vNadirLon, IN_vNadirLat):
    '''
    Compute the azimuth index associated to the point with coordinates P(IN_lon, IN_lat) with respect to the nadir track with coordinates (IN_vNadirLon, IN_vNadirLat) for each point.
    NB: the method used is to minimize the scalar product between P and the points of the nadir track.
    
    :param IN_lon: longitude in degrees east
    :type IN_lon: float
    :param IN_lat: latitude in degrees north
    :type IN_lat: float
    :param IN_vNadirLon: longitude of each nadir point in degrees east
    :type IN_vNadirLon: 1D-array
    :param IN_vNadirLat: latitude of each nadir point in degrees north
    :type IN_vNadirLat: 1D-array
    
    :return OUT_idx_min: azimuth index corresponding to input point (IN_lon, IN_lat)
    :type: int
    '''
    
    # 1 - Init variables
    nb_nadir = IN_vNadirLon.size
    list_scalar = np.zeros(nb_nadir)
    
    # 2 - Compute scalar product for each nadir point
    # 2.1 - 1st point
    list_scalar[0] = (IN_vNadirLon[0]-IN_lon)*(IN_vNadirLon[1]-IN_vNadirLon[0]) + (IN_vNadirLat[0]-IN_lat)*(IN_vNadirLat[1]-IN_vNadirLat[0])
    # 2.2 - Middle points
    for ind in np.arange(1,nb_nadir-1):
        list_scalar[ind] = (IN_vNadirLon[ind]-IN_lon)*(IN_vNadirLon[ind+1]-IN_vNadirLon[ind-1]) + (IN_vNadirLat[ind]-IN_lat)*(IN_vNadirLat[ind+1]-IN_vNadirLat[ind-1])
    # 2.3 - Last point
    list_scalar[-1] = (IN_vNadirLon[-1]-IN_lon)*(IN_vNadirLon[-1]-IN_vNadirLon[-2]) + (IN_vNadirLat[-1]-IN_lat)*(IN_vNadirLat[-1]-IN_vNadirLat[-2])
    
    # 3 - Find min scalar
    OUT_idx_min = np.argmin(np.absolute(list_scalar))
    
    # 4 - Return corresponding azimuth index
    return OUT_idx_min
    
#######################################

def computeDist(IN_lon1, IN_lat1, IN_lon2, IN_lat2):
    '''
    Compute Euler distance between 2 points of the plan: P1(lon1, lat1) and P2(lon2, lat2)
    
    :param IN_lon1: longitude of P1 in degrees east
    :type IN_lon1: float
    :param IN_lat1: latitude of P1 in degrees north
    :type IN_lat1: float
    :param IN_lon2: longitude of P2 in degrees east
    :type IN_lon2: float
    :param IN_lat2: latitude of P2 in degrees north
    :type IN_lat2: float
    
    :return the spherical distance (in meters) between P1 and P2
    :type: float
    '''
    
    # 1 - Convert degrees to radians
    lon1 = deg2rad(IN_lon1)
    lat1 = deg2rad(IN_lat1)
    lon2 = deg2rad(IN_lon2)
    lat2 = deg2rad(IN_lat2)
    
    # 1 - Distance angulaire en radians
    s12 = np.arccos(np.sin(lat1)*np.sin(lat2) + np.cos(lat1)*np.cos(lat2)*np.cos(lon2-lon1))
    
    return s12 * GEN_APPROX_RAD_EARTH
    
#######################################

def llh2xyz(IN_lon, IN_lat, IN_height, IN_flag_rad=True):
    '''
    Convert geographic coordinates (longitude, latitude, height) to cartesian coordinates (x, y, z) 
        
    :param IN_lon: longitude in degrees east
    :type IN_lon: scalar or 1D-array
    :param IN_lat: latitude in degrees north
    :type IN_lat: scalar or 1D-array
    :param IN_height: height in meters
    :type IN_height: scalar or 1D-array
    :param IN_flag_rad: =False if IN_lan and IN_lat are in degrees, =True if they are in radians (default)
    :type IN_flag_rad: boolean
    
    :return OUT_x, OUT_y, OUT_z: cartesian coordinates
    :type: scalars or 1D-array
    
    .. author: Damien DESROCHES & Alejandro BOHE - CNES DSO/SI/TR
    '''
    
    # 1 - Convert geodetic coordinates in radians
    if ( IN_flag_rad ):
        lat = IN_lat
        lon = IN_lon
    else:
        lat = deg2rad(IN_lat)
        lon = deg2rad(IN_lon)
    
    # 2 - Compute earth excentricity
    e = np.sqrt(GEN_RAD_EARTH_EQ ** 2 - GEN_RAD_EARTH_POLE ** 2) / GEN_RAD_EARTH_EQ

    # 3 - Compute earth radius for latitude lat
    Rn = GEN_RAD_EARTH_EQ / np.sqrt(1 - e**2 * np.sin(lat)**2 )
    
    # 4 - Compute cartesian coordinates
    OUT_x = (Rn + IN_height) * np.cos(lat) * np.cos(lon)
    OUT_y = (Rn + IN_height) * np.cos(lat) * np.sin(lon)
    OUT_z = (Rn * (1.-e**2) + IN_height) * np.sin(lat)

    return OUT_x, OUT_y, OUT_z   

def xyz2llh(IN_x, IN_y, IN_z):
    '''
    Convert cartesian coordinates (x, y, z) to geographic coordinates (lat, lon, height) 
    
    :param IN_x: coordinate along x-axis
    :type IN_x: scalar or 1D-array
    :param IN_y: coordinate along y-axis
    :type IN_y: scalar or 1D-array
    :param IN_z: coordinate along z-axis
    :type IN_z: scalar or 1D-array
    
    :return OUT_long, OUT_lat, OUT_height: geographic coordinates (longitude, latitude, height) 
    :rtype: same as input (scalar or 1D-array)
    '''
    
    # 3.0 - Compute ellipsoide variables
    e = np.sqrt(GEN_RAD_EARTH_EQ ** 2 - GEN_RAD_EARTH_POLE ** 2) / GEN_RAD_EARTH_EQ # Ellipsoid excentricity
    f = 1 - np.sqrt(1 - e**2) # Flattening factor
    r = np.sqrt(IN_x**2 + IN_y**2 + IN_z**2) # Distance between center of earth and point
    d = np.sqrt(IN_x**2 + IN_y**2) # Distance between center of earth and point with latitude 0
    
    # 3.1 - Compute longitudes
    OUT_lon = np.arctan2(IN_y, IN_x)
    
    # 3.2 - Compute latitudes
    mu = np.arctan(IN_z/d * ( (1-f) + e**2 * GEN_RAD_EARTH_EQ/r))
    OUT_lat = np.arctan( (IN_z*(1-f)+e**2 * GEN_RAD_EARTH_EQ * np.sin(mu)**3) / ((1-f)*(d-e**2 * GEN_RAD_EARTH_EQ * np.cos(mu)**3)) )
    
    # 3.3 - Compute heights
    OUT_height = (d*np.cos(OUT_lat)) + (IN_z*np.sin(OUT_lat)) - GEN_RAD_EARTH_EQ*np.sqrt(1-(e**2)*np.sin(OUT_lat)**2)
    
    return rad2deg(OUT_lon), rad2deg(OUT_lat), OUT_height

def compter(liste):
    compte = {}.fromkeys(set(liste),0)
    for valeur in liste:
        compte[valeur] += 1
    return compte 
    
def compter_area(liste, tab_area):
    
    compte = {}.fromkeys(set(liste),0)
    for i, valeur in enumerate(liste):
        compte[valeur] += tab_area[i]
    return compte
    
    
def gdal_check_call(args):
    #~ gdal_cmd = " ".join(args)
    gdal_cmd = args
    cmd = "module purge; module load gdal_hdf5; {}".format(gdal_cmd)
    print(cmd)
    check_call(cmd, shell=True)
    
#######################################

if __name__ == '__main__':
    
    #print xyz2llh(5289084.19826832, 27760.13308768, 4983805.89485445)
    print(llh2xyz(0.007955, 43.525029, 9.240165, False))
