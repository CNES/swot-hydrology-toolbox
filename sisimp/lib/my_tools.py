# -*- coding: utf-8 -*-
'''
.. module my_tools.py
    :synopsis: library with generic functions
    03/02/2017 - Version 0.10 = creation 
    02/02/2018 - Version 0.11 = merging of SAM and stand-alone specific code (use of my_api library)
                                + add computeMean_2sigma function
    03/06/2018 - Version 0.12 = add alpha_shape function
    
***** RESTE A FAIRE *****
- faire une version alpha_shape / concavHull commune avec Damien et Emmanuelle (utilisation sur le floodplain DEM)

.. module author: Claire POTTIER - CNES DSO/SI/TR

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


'''

import os
import numpy as np
from scipy.ndimage.measurements import label
import math
from scipy.spatial import Delaunay
from shapely.ops import cascaded_union
import shapely.geometry as geometry

import lib.my_api as my_api
from lib.my_variables import GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE, GEN_APPROX_RAD_EARTH
        

def testFile(IN_file, IN_extent=None):
    """
    Test if full path in input is an existing file and, optionnally, a file in the awaited format
    
    :param IN_file: input full path
    :type IN_file: string
    :param IN_extent: awaited format file (optionnal)
    :type IN_extent: string
    """
    
    if os.path.exists(IN_file):
        if not os.path.isfile(IN_file):
            my_api.exitWithError("ERROR = %s is not a file" % IN_file)
        if IN_extent and not IN_file.endswith(IN_extent):
            my_api.exitWithError("ERROR = %s must have the %s extent" % (IN_file, IN_extent))
    else:
        my_api.exitWithError("ERROR = %s doesn't exist" % IN_file)
        
        
def testDir(IN_dir):
    '''
    Test if input path is an existing directory
    
    :param IN_dir: input path
    :type IN_dir: string
    '''
    
    if ( os.path.exists(IN_dir) ):
        if ( os.path.isdir(IN_dir) == False ):
            my_api.exitWithError("ERROR = %s is not a directory" % ( IN_dir ))
    else:
        my_api.exitWithError("ERROR = %s doesn't exist" % ( IN_dir ))
        
#######################################
    
def convert_to_m180_180(in_long):
    """
    Convert longitudes from [0;360[ to [-180;180[
    
    :param in_long: longitudes to convert
    :type in_long: float or 1D-array of float
    
    :return: out_long = converted longitude
    :rtype: same as input = float or 1D-array of float
    """
    
    out_long = in_long
    ind = np.where(in_long > 180.0)
    if ind is not None:
        out_long[ind] -= 360.
        
    return out_long
    

def convert_to_0_360(in_long):
    """
    Convert longitudes from [-180;180[ to [0;360[ 
    
    :param in_long: longitudes to convert
    :type in_long: float or 1D-array of float
    
    :return: out_long = converted longitude
    :rtype: same as input = float or 1D-array of float
    """
    
    out_long = in_long
    ind = np.where(in_long < 0.0)
    if ind is not None:
        out_long[ind] += 360.
        
    return out_long

#######################################

def compute_interferogram(plus_y_ant_x, plus_y_ant_y, plus_y_ant_z, minus_y_ant_x, minus_y_ant_y, minus_y_ant_z, tile_x, tile_y, tile_z):
    dist_plus = np.sqrt((plus_y_ant_x - tile_x)**2 \
    +(plus_y_ant_y - tile_y)**2 \
    +(plus_y_ant_z - tile_z)**2)

    dist_minus = np.sqrt((minus_y_ant_x - tile_x)**2 \
    +(minus_y_ant_y - tile_y)**2 \
    +(minus_y_ant_z - tile_z)**2)

    phase_ref = -2*np.pi/0.008385803020979*(dist_plus - dist_minus)
    interferogram = np.exp(1.j*phase_ref)

    return [interferogram.real, interferogram.imag]
        
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
    my_api.printDebug("[my_tools] == computeBinMat ==")
    
    # 0 - Deal with exceptions
    # 0.1 - Input vectors size must be the same
    if IN_X.size != IN_Y.size:
        raise ValueError("computeBinMat(IN_X, IN_Y) : IN_X and IN_Y must be the same size ; currently : IN_X = %d and IN_Y = %d" % ( IN_X.size , IN_Y.size ))
    else:
        nb_pts = IN_X.size
    my_api.printDebug("> Nb pixels to deal with = %d" % ( nb_pts ))
    # 0.2 - max(X) < IN_sizeX
    if np.max(IN_X) >= IN_sizeX:
        raise ValueError("computeBinMat(IN_X, IN_Y) : elements of IN_X must be less than IN_sizeX")
    # 0.3 - max(X) < IN_sizeX
    if np.max(IN_Y) >= IN_sizeY:
        raise ValueError("computeBinMat(IN_X, IN_Y) : elements of IN_Y must be less than IN_sizeY")
    
    # 1 - Init output binary image
    OUT_binIm = np.zeros((IN_sizeY, IN_sizeX))
    my_api.printDebug("> Binary matrix size = (X=%d , Y=%d)" % ( IN_sizeX , IN_sizeY ))
        
    # 2 - Put 1 for every pixels defined by the input vectors
    for ind in range(nb_pts): 
        OUT_binIm[IN_Y[ind], IN_X[ind]] = 1
    
    return OUT_binIm
        
#######################################
    
def labelRegion(IN_binMat):
    '''
    Identifies all separate regions in a 2D binary matrix
    
    :param IN_api: SAM api to print log
    :type IN_api: SAM api
    :param IN_binMat: 2D binary matrix
    :type IN_binMat: 2D binary matrix of int 0/1
    
    :return: a 2D matrix, in which each pixel has a unique value ranging from 0 to N (number of identified regions)
    :rtype: 2D matrix of int 
    :return: the number of objects
    :rtype: int
    '''
    my_api.printDebug("[my_tools] == labelRegion ==")
    
    OUT_mLabelRegions, OUT_nbObj = label(IN_binMat)
    my_api.printDebug("> Number of detected objects = %d" % ( OUT_nbObj ) )
    
    return OUT_mLabelRegions, OUT_nbObj
    
#######################################

def convert2dMatIn1dVec(IN_X, IN_Y, IN_mat):
    '''
    Converts a 2D matrix [ IN_X[ind], IN_Y[ind] ] in a 1D vector [ind]
    
    :param IN_api: SAM api to print log
    :type IN_api: SAM api
    :param IN_mat: 2D matrix
    :type IN_mat: 2D matrix of values (int, float, ...)
    
    :return: a 1D vector, each value at "ind" being the value of IN_mat[ IN_X[ind], IN_Y[ind] ]
    :rtype: 2D matrix of values (int, float, ...)
    '''
    my_api.printDebug("[my_tools] == convert2dMatIn1dVec ==")
    
    # Number of pixels in IN_X
    nb_pix = IN_X.size
    # Init output vector (same size of input IN_X and IN_Y)   
    OUT_vector = np.zeros(nb_pix)
    # Fill the output vector
    for indp in range(nb_pix):
        OUT_vector[indp] = IN_mat[IN_Y[indp], IN_X[indp]]
    
    return OUT_vector
    
#######################################

def alpha_shape(IN_coords, IN_alpha):
    '''
    Compute the alpha shape (concave hull) of a set of points.
    
    :param IN_coords: set of points coordinates
    :type IN_coords: 2D-array of size (nb_pixels, 2=lon/lat)
    :param IN_alpha: alpha value to influence the gooeyness of the border. Smaller numbers don't fall inward as much as larger numbers. Too large, and you lose everything!
    :type IN_alpha: int
    
    :return Shapely.MultiPolygons which is the hull of the input set of points
    
    Retrieved from http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
    '''
    
    # Number of points
    nb_pts = IN_coords.shape[0]
    
    # 0 - Particular case of nb points <= 3: return convex hull
    if ( nb_pts <= 3 ):
        # 0.1 - Init Shapely set of points
        points = []
        # 0.2 - Aggregate coordinates to the set of points
        for indp in np.arange(nb_pts):
            points.append(geometry.point.Point(IN_coords[indp,0],IN_coords[indp,1],0))
        # 0.3 - Return convex hull
        return geometry.MultiPoint(list(points)).convex_hull

    # 1 - Compute Delaunay triangulation
    tri = Delaunay(IN_coords) # tri = specific object with attributes (cf. https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.Delaunay.html)
    
    # 2 - Loop over triangles
    list_triangle = []
    for ia, ib, ic in tri.vertices: # ia, ib, ic = indices of corner points of the triangle
        pa = IN_coords[ia]
        pb = IN_coords[ib]
        pc = IN_coords[ic]
        # Lengths of sides of triangle
        a = math.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
        b = math.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
        c = math.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)
        # Semiperimeter of triangle
        s = (a + b + c)/2.0
        # Area of triangle by Heron's formula
        area = math.sqrt(s*(s-a)*(s-b)*(s-c))
        circum_r = 0
        if ( area != 0 ):
            circum_r = a*b*c/(4.0*area)
        # Here's the radius filter
        if circum_r < 1.0/IN_alpha:
            gg = [(pa[0], pa[1]),(pb[0], pb[1]),(pc[0], pc[1])]
            tt = geometry.Polygon(gg)
            list_triangle.append(tt)    
    
    return cascaded_union(list_triangle)
    
#######################################

def cptDigits(IN_real):
    '''
    Compute the position of the first significant figure of a real number (ex: 110 => 3 ; 0.000156 => -4)
        
    :param IN_real: real number
    :type IN_real: float
    
    :return the position of the first significant figure of a real number
    :rtype: int 
    '''
    
    # Cas particulier de 0
    if ( IN_real == 0 ):
        OUT_cpt = 1
    
    elif ( abs(IN_real) >= 10 ): # Cas des nombres >= 10
        OUT_cpt = 1 
        while ( abs(IN_real) >= 10 ):
            IN_real /= 10.0
            OUT_cpt += 1
            
    elif ( abs(IN_real) < 1 ): # Cas des nombres < 1
        OUT_cpt = 0
        while ( abs(IN_real) < 1 ):
            IN_real *= 10.0
            OUT_cpt -= 1
            
    else:
        OUT_cpt = 1
        
    return OUT_cpt
    
#######################################

def convertSec2Time(IN_time, txt_format=1):
    '''
    Convert a date_time in seconds towards a string with the format specified as in put parameter
    
    :param IN_time: date_time in seconds
    :type IN_time: float
    :param txt_format(optionnal): 1="hh:mm:ss" 2="[hh]h [mm]min [ss]s" 3="hh:mm:ss.ms" 4="DD hh:mm:ss"
    :type txt_format: int (defaut = 1)
    
    :return: la duree correspondant a IN_time au format choisi
    :rtype: str
    '''
    
    # Calcul des secondes
    tmp_MM = math.floor(IN_time / 60.0)
    SS = IN_time - 60.0 * tmp_MM
        
    # Calcul des minutes
    HH = math.floor(tmp_MM / 60.0)
    MM = tmp_MM - 60.0 * HH
        
    # Calcul des heures si > 24h
    DD = math.floor(HH / 24.0)
    HH -= 24.0 * DD
    if ( txt_format == 1 ):
        return "%02d:%02d:%02d" % ( HH , MM , SS )
    elif ( txt_format == 2 ):
        return "%02dh%02dmin%02ds" % ( HH , MM , SS )
    elif ( txt_format == 3 ):
        return "%02d:%02d:%02.3f" % ( HH , MM , SS )
    elif ( txt_format == 4 ):
        return "Day %02d %02d:%02d:%02d" % ( DD+1 , HH , MM , SS ) # DD+1 car le 1e jour est le jour 1
    
#######################################

def computeMean_2sigma(IN_vVal):
    '''
    Compute the mean of the input values, after remove of non-sense values, i.e. below or above the median +/- 2*standard deviation
    
    :param IN_vVal: vector of values for which the mean is computed
    :type IN_vVal: 1D-array of float
        
    :return the mean value
    :rtype: float
    '''
    
    # 1 - Compute statistical values over this vector
    # 1.1 - Median
    med = np.median(IN_vVal)
    # 1.2 - Standard deviation
    std = np.std(IN_vVal)
    # 1.3 - Nb values
    nb_val = IN_vVal.size
        
    # 2 - Remove indices with value out of 2 sigma
    # 2.1 - Lower than mean - 2 sigma
    idx_lower = np.where(IN_vVal < med-2*std)[0]
    vIndices = np.delete(np.arange(nb_val), idx_lower)
    # 2.2 - Higher than mean + 2 sigma
    idx_higher = np.where(IN_vVal[vIndices] > med+2*std)[0]
    vIndices = np.delete(vIndices, idx_higher)
    
    my_api.printDebug("[my_tools/computeMean_2sigma] %d / %d pixels used for mean computation" % ( vIndices.size , nb_val ))
        
    # 3 - Return the mean of cleaned vector
    return np.mean(IN_vVal[vIndices])
    
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

#######################################

def coords_from_labels(matrix):
    nb_col, nb_line = matrix.shape
    labels_coords = {}

    processing = np.round(np.linspace(0, nb_col, 11), 0)
    for i in range(nb_col):
        if i in processing :
            my_api.printInfo("[my_tools] [coords_from_labels] Processing %d %%" %( int(100 * (i+1)/nb_col) ) )
        for j in range(nb_line):
            val = matrix[i, j]
            labels_coords.setdefault(val, []).append((i,j))

    return labels_coords


if __name__ == '__main__':
    
    a = np.arange(10)
    a[1] = -1000
    a[9] = 1000
    print(a)
    print(computeMean_2sigma(a))
