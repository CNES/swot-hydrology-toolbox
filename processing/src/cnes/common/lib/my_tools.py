# -*- coding: utf8 -*-
"""
.. module:: my_tools.py
    :synopsis: library with generic functions
    Created on 03/02/2017

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR
    
.. todo:: faire une version alpha_shape / concavHull commune avec Damien et Emmanuelle (utilisation sur le floodplain DEM)

Copyright (c) 2017 CNES. All rights reserved.
"""

import math
import numpy as np
import os
from osgeo import osr
from scipy.ndimage.measurements import label
from skimage.morphology import square
from sklearn.cluster import KMeans
import logging

import skimage
if skimage.__version__ >= "0.11":
    from skimage.filters import median as median_filter
else:
    from skimage.filter.rank import median as median_filter

import cnes.common.lib.my_api as my_api
from cnes.common.lib.my_variables import GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE, GEN_APPROX_RAD_EARTH
import cnes.common.lib_lake.locnes_variables as my_var
import cnes.common.serviceError as serviceError


def testFile(IN_file, IN_extent=None):
    """
    Test if full path in input is an existing file and, optionnally, a file in the awaited format
    
    :param IN_file: input full path
    :type IN_file: string
    :param IN_extent: awaited format file (optionnal)
    :type IN_extent: string
    """
    logger = logging.getLogger("my_tools")
    
    if os.path.exists(IN_file):
        if not os.path.isfile(IN_file):
            message = "ERROR = %s is not a file" % IN_file
            raise serviceError.ProcessingError(message, logger)
        if IN_extent and not IN_file.endswith(IN_extent):
            message = "ERROR = %s must have the %s extent" % (IN_file, IN_extent)
            raise serviceError.ProcessingError(message, logger)
    else:
        message = "ERROR = %s doesn't exist" % IN_file
        raise serviceError.ProcessingError(message, logger)


def testDir(IN_dir):
    """
    Test if input path is an existing directory
    
    :param IN_dir: input path
    :type IN_dir: string
    """
    logger = logging.getLogger("my_tools")
    
    if os.path.exists(IN_dir):
        if not os.path.isdir(IN_dir):
            message = "ERROR = %s is not a directory" % IN_dir
            raise serviceError.ProcessingError(message, logger)
    else:
        message = "ERROR = %s doesn't exist" % IN_dir
        raise serviceError.ProcessingError(message, logger)


#######################################


def deg2rad(IN_deg):
    """
    Convert angles in degrees to radians
    
    :param IN_deg: angles to convert
    :type IN_deg: scalar or 1D-array
    
    :return angles in radians
    :rtype: same as input
    """
    return IN_deg * np.pi / 180.


def rad2deg(IN_rad):
    """
    Convert angles in radians to degrees
    
    :param IN_rad: angles to convert
    :type IN_rad: scalar or 1D-array
    
    :return angles in degrees
    :rtype: same as input
    """
    return IN_rad * 180. / np.pi


#######################################


def computeBinMat(IN_sizeX, IN_sizeY, IN_X, IN_Y):
    """
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
    """
    logger = logging.getLogger("my_tools")
    logger.debug("[my_tools] == computeBinMat ==")
    
    # 0 - Deal with exceptions
    # 0.1 - Input vectors size must be the same
    if IN_X.size != IN_Y.size:
        message = "computeBinMat(IN_X, IN_Y) : IN_X and IN_Y must be the same size ; currently : IN_X = %d and IN_Y = %d" % (IN_X.size, IN_Y.size)
        raise serviceError.ProcessingError(message, logger)
    else:
        nb_pts = IN_X.size
    logger.debug("> Nb pixels to deal with = %d" % nb_pts)
    # 0.2 - max(X) < IN_sizeX
    if np.max(IN_X) >= IN_sizeX:
        message = "computeBinMat(IN_X, IN_Y) : elements of IN_X must be less than IN_sizeX"
        raise serviceError.ProcessingError(message, logger)
    # 0.3 - max(X) < IN_sizeX
    if np.max(IN_Y) >= IN_sizeY:
        message = "computeBinMat(IN_X, IN_Y) : elements of IN_Y must be less than IN_sizeY"
        raise serviceError.ProcessingError(message, logger)
    # 1 - Init output binary image
    OUT_binIm = np.zeros((IN_sizeY, IN_sizeX))
    logger.debug("> Binary matrix size = (X=%d , Y=%d)" % (IN_sizeX, IN_sizeY))
        
    # 2 - Put 1 for every pixels defined by the input vectors
    for ind in range(nb_pts): 
        OUT_binIm[IN_Y[ind], IN_X[ind]] = 1
    
    return OUT_binIm


#######################################


def labelRegion(IN_binMat):
    """
    Identifies all separate regions in a 2D binary matrix

    :param IN_binMat: 2D binary matrix
    :type IN_binMat: 2D binary matrix of int 0/1
    
    :return: a 2D matrix, in which each pixel has a unique value ranging from 0 to N (number of identified regions)
    :rtype: 2D matrix of int 
    :return: the number of objects
    :rtype: int
    """
    logger = logging.getLogger("my_tools")
    logger.debug("[my_tools] == labelRegion ==")
    
    OUT_mLabelRegions, OUT_nbObj = label(IN_binMat)
    message = "> Number of detected objects = %d" % OUT_nbObj
    logger.debug(message)
    
    return OUT_mLabelRegions, OUT_nbObj


def relabelLakeUsingSegmentationHeigth(IN_X, IN_Y, IN_height):
    """
    This function main interest is to determine the number of lakes inside a subset of PixC in radar geometry.
    In most of the cases, only one lake is located in the given subset, but in some case of different lakes can be gather inside one because of radar geometric distortions or in the case a of a dam.
    Steps :
        1. Creates a 2D height matrix from X and Y 1D vectors
        2. Unsupervised heigth classification to determine number a classes to get an STD over each classe lower than STD_HEIGHT_MAX
        3. Smooth result using a 2D median filter
        4. Return labels

    :param IN_X: X indices of "1" pixels
    :type IN_X: 1D vector of int
    :param IN_Y: Y indices of "1" pixels
    :type IN_Y: 1D vector of int
    :param IN_height: height of pixels
    :type IN_height: 1D vector of float

    :return: labels recomputed over 1 of several lakes
    :rtype: 1D vector of int
    """

    logger = logging.getLogger("my_tools")
    logger.debug("[my_tools] == relabelLakeUsingSegmentationHeigth ==")

    logger.debug("[my_tools] Building heigth matrix")

    # 0 - Deal with exceptions
    # 0.1 - Input vectors size must be the same
    if (IN_X.size != IN_Y.size) or (IN_X.size != IN_height.size):
        raise ValueError("relabelLakeUsingSegmentationHeigth(IN_X, IN_Y, IN_height) : IN_X and IN_Y must be the same size ; currently : IN_X = %d and IN_Y = %d" % (IN_X.size, IN_Y.size))
    else:
        nb_pts = IN_X.size

    # 1 - Compute height matrix
    # 1.1 - Init heigth image
    heightImg = np.zeros((np.max(IN_Y) + 1, np.max(IN_X) + 1))
    message = "> Height matrix size = (X=%d , Y=%d)" % (np.max(IN_Y.size) + 1, np.max(IN_X.size) + 1)
    logger.debug(message)

    # 1.2 - Put height for every pixels defined by the input vectors
    for ind in range(nb_pts):
        heightImg[IN_Y[ind], IN_X[ind]] = IN_height[ind]

    logger.debug("[my_tools] K-means processing")

    # 2 - Unsupervised clustering to determine number of classes
    std_heigth = np.std(IN_height)
    nb_classes = 1

    # 2.1 - Cas with only one class
    if std_heigth < my_var.STD_HEIGHT_MAX:  # If only one lake is in the given pixc subset
        return np.ones(nb_pts)  # Return one unique label for all pixc

    # 2.2 - Init a k-means classifier
    kmeans_classif = KMeans()
    while std_heigth > my_var.STD_HEIGHT_MAX:

        nb_classes += 1

        # 2.3 - Cluster height over nb_classes classes
        kmeans_classif = KMeans(n_clusters=nb_classes)
        kmeans_classif.fit(IN_height.reshape(-1, 1))  # Turn line into column for IN_height

        # 2.3 - Compute heigth std inside each class
        std_heigth = np.max([np.std(IN_height[np.where(kmeans_classif.labels_ == curLabel)]) for curLabel in np.unique(kmeans_classif.labels_)])

        # If number of classes upper than 10 => stop iteration
        if nb_classes > 10:
            break

    message = "[my_tools] NB classes : %d, max std : %f " % (nb_classes, std_heigth)
    logger.debug(message)

    # 3 - Format output vector

    # 3.1 - Build a labeled matrix
    labeledImg = np.zeros(heightImg.shape).astype('int')
    for ind in range(nb_pts):
        labeledImg[IN_Y[ind], IN_X[ind]] = int(kmeans_classif.labels_[ind])

    # 3.2 - Median filter on a 3x3 window to smooth output
    labeledImg_filted = median_filter(labeledImg.astype('uint8'), square(2).astype('uint8'))

    # 3.3 - Init and fill output label array
    OUT_labels = np.zeros(IN_height.shape)
    for ind in range(nb_pts):
        OUT_labels[ind] = labeledImg_filted[IN_Y[ind], IN_X[ind]] + 1

    return OUT_labels


#######################################


def convert2dMatIn1dVec(IN_X, IN_Y, IN_mat):
    """
    Converts a 2D matrix [ IN_X[ind], IN_Y[ind] ] in a 1D vector [ind]
    
    :param IN_X: column indices of points of IN_mat to return
    :type IN_X: 1D array
    :param IN_Y: row indices of points of IN_mat to return
    :type IN_Y: 1D array
    :param IN_mat: 2D array
    :type IN_mat: 2D array of values (int, float, ...)
    
    :return: a 1D vector, each value at "ind" being the value of IN_mat[ IN_X[ind], IN_Y[ind] ]
    :rtype: 2D matrix of values (int, float, ...)
    """
    logger = logging.getLogger("my_tools")
    logger.debug("[my_tools] == convert2dMatIn1dVec ==")
    
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
    """
    Compute the position of the first significant figure of a real number (ex: 110 => 3 ; 0.000156 => -4)
        
    :param IN_real: real number
    :type IN_real: float
    
    :return the position of the first significant figure of a real number
    :rtype: int 
    """
    
    # Cas particulier de 0
    if IN_real == 0:
        OUT_cpt = 1
    
    elif abs(IN_real) >= 10:  # Cas des nombres >= 10
        OUT_cpt = 1 
        while abs(IN_real) >= 10:
            IN_real /= 10.0
            OUT_cpt += 1
            
    elif abs(IN_real) < 1:  # Cas des nombres < 1
        OUT_cpt = 0
        while abs(IN_real) < 1:
            IN_real *= 10.0
            OUT_cpt -= 1
            
    else:
        OUT_cpt = 1
        
    return OUT_cpt


#######################################


def convertSec2Time(IN_time, txt_format=1):
    """
    Convert a date_time in seconds towards a string with the format specified as in put parameter
    
    :param IN_time: date_time in seconds
    :type IN_time: float
    :param txt_format: [optionnal] 1="hh:mm:ss" 2="[hh]h [mm]min [ss]s" 3="hh:mm:ss.ms" 4="Day DD hh:mm:ss"
    :type txt_format: int (defaut = 1)
    
    :return: la duree correspondant a IN_time au format choisi
    :rtype: str
    """
    
    # Calcul des secondes
    tmp_MM = math.floor(IN_time / 60.0)
    SS = IN_time - 60.0 * tmp_MM
        
    # Calcul des minutes
    HH = math.floor(tmp_MM / 60.0)
    MM = tmp_MM - 60.0 * HH
        
    # Calcul des heures si > 24h
    DD = math.floor(HH / 24.0)
    HH -= 24.0 * DD
    if txt_format == 1:
        return "%02d:%02d:%02d" % (HH, MM, SS)
    elif txt_format == 2:
        return "%02dh%02dmin%02ds" % (HH, MM, SS)
    elif txt_format == 3:
        return "%02d:%02d:%02.3f" % (HH, MM, SS)
    elif txt_format == 4:
        return "Day %02d %02d:%02d:%02d" % (DD+1, HH, MM, SS)  # DD+1 car le 1e jour est le jour 1


#######################################


def computeMean_2sigma(IN_vVal):
    """
    Compute the mean of the input values, after remove of non-sense values, i.e. below or above the median +/- 2*standard deviation
    
    :param IN_vVal: vector of values for which the mean is computed
    :type IN_vVal: 1D-array of float
        
    :return the mean value
    :rtype: float
    """
    
    logger = logging.getLogger("my_tools")
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
    
    message = "[my_tools/computeMean_2sigma] %d / %d pixels used for mean computation" % (vIndices.size, nb_val)
    logger.debug(message)
        
    # 3 - Return the mean of cleaned vector
    return np.mean(IN_vVal[vIndices])


#######################################


def computeAz(IN_lon, IN_lat, IN_vNadirLon, IN_vNadirLat):
    """
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
    """
    
    # 1 - Init variables
    nb_nadir = IN_vNadirLon.size
    list_scalar = np.zeros(nb_nadir)
    
    # 2 - Compute scalar product for each nadir point
    # 2.1 - 1st point
    list_scalar[0] = (IN_vNadirLon[0]-IN_lon)*(IN_vNadirLon[1]-IN_vNadirLon[0]) + (IN_vNadirLat[0]-IN_lat)*(IN_vNadirLat[1]-IN_vNadirLat[0])
    # 2.2 - Middle points
    for ind in np.arange(1, nb_nadir-1):
        list_scalar[ind] = (IN_vNadirLon[ind]-IN_lon)*(IN_vNadirLon[ind+1]-IN_vNadirLon[ind-1]) + (IN_vNadirLat[ind]-IN_lat)*(IN_vNadirLat[ind+1]-IN_vNadirLat[ind-1])
    # 2.3 - Last point
    list_scalar[-1] = (IN_vNadirLon[-1]-IN_lon)*(IN_vNadirLon[-1]-IN_vNadirLon[-2]) + (IN_vNadirLat[-1]-IN_lat)*(IN_vNadirLat[-1]-IN_vNadirLat[-2])
    
    # 3 - Find min scalar
    OUT_idx_min = np.argmin(np.absolute(list_scalar))
    
    # 4 - Return corresponding azimuth index
    return OUT_idx_min


#######################################


def computeDist(IN_lon1, IN_lat1, IN_lon2, IN_lat2):
    """
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
    """
    
    # 1 - Convert degrees to radians
    lon1 = deg2rad(IN_lon1)
    lat1 = deg2rad(IN_lat1)
    lon2 = deg2rad(IN_lon2)
    lat2 = deg2rad(IN_lat2)
    
    # 1 - Distance angulaire en radians
    s12 = np.arccos(np.sin(lat1)*np.sin(lat2) + np.cos(lat1)*np.cos(lat2)*np.cos(lon2-lon1))
    
    return s12 * GEN_APPROX_RAD_EARTH


#######################################


def getArea(IN_polygon, IN_centroid):
    """
    This function projects IN_polygon from geographic coordinate into centroid corresponding UTM zone.
    Once projected, the area of the polygon is computed and returned
    
    :param IN_polygon: polygon
    :type IN_polygon: OGR geometry
    :param IN_centroid: centroid of IN_polygon
    :type IN_centroid: tuple of coordinates
    
    :return: area of IN_polygon
    :rtype: float
    """
    
    # 0 - Computation of EPSG code corresponding UTM zone
    epsg_code = "32"  # Initialization

    centroid_lon = IN_centroid[0]
    centroid_lat = IN_centroid[1]

    if centroid_lat > 0:
        epsg_code += "6"
    else:
        epsg_code += "7"

    epsg_code += str(int(1 + (centroid_lon + 180)//6))

    # 1 - Projection of IN_polygon into UTM
    src_source = osr.SpatialReference()
    src_source.ImportFromEPSG(4326)

    src_target = osr.SpatialReference()
    src_target.ImportFromEPSG(int(epsg_code))

    transform = osr.CoordinateTransformation(src_source, src_target)

    IN_polygon.Transform(transform)

    # 2 - Compute and return area
    return IN_polygon.GetArea()


#######################################


def llh2xyz(IN_lon, IN_lat, IN_height, IN_flag_rad=True):
    """
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
    """
    
    # 1 - Convert geodetic coordinates in radians
    if IN_flag_rad:
        lat = IN_lat
        lon = IN_lon
    else:
        lat = deg2rad(IN_lat)
        lon = deg2rad(IN_lon)
    
    # 2 - Compute earth excentricity
    e = np.sqrt(GEN_RAD_EARTH_EQ ** 2 - GEN_RAD_EARTH_POLE ** 2) / GEN_RAD_EARTH_EQ

    # 3 - Compute earth radius for latitude lat
    Rn = GEN_RAD_EARTH_EQ / np.sqrt(1 - e**2 * np.sin(lat)**2)
    
    # 4 - Compute cartesian coordinates
    OUT_x = (Rn + IN_height) * np.cos(lat) * np.cos(lon)
    OUT_y = (Rn + IN_height) * np.cos(lat) * np.sin(lon)
    OUT_z = (Rn * (1.-e**2) + IN_height) * np.sin(lat)

    return OUT_x, OUT_y, OUT_z   


def xyz2llh(IN_x, IN_y, IN_z):
    """
    Convert cartesian coordinates (x, y, z) to geographic coordinates (lat, lon, height) 
    
    :param IN_x: coordinate along x-axis
    :type IN_x: scalar or 1D-array
    :param IN_y: coordinate along y-axis
    :type IN_y: scalar or 1D-array
    :param IN_z: coordinate along z-axis
    :type IN_z: scalar or 1D-array
    
    :return OUT_long, OUT_lat, OUT_height: geographic coordinates (longitude, latitude, height) 
    :rtype: same as input (scalar or 1D-array)
    """
    
    # 1 - Compute ellipsoide variables
    e = np.sqrt(GEN_RAD_EARTH_EQ ** 2 - GEN_RAD_EARTH_POLE ** 2) / GEN_RAD_EARTH_EQ  # Ellipsoid excentricity
    f = 1 - np.sqrt(1 - e**2)  # Flattening factor
    r = np.sqrt(IN_x**2 + IN_y**2 + IN_z**2)  # Distance between center of earth and point
    d = np.sqrt(IN_x**2 + IN_y**2)  # Distance between center of earth and point with latitude 0
    
    # 2 - Compute longitudes
    OUT_lon = np.arctan2(IN_y, IN_x)
    
    # 3 - Compute latitudes
    mu = np.arctan(IN_z/d * ((1-f) + e**2 * GEN_RAD_EARTH_EQ/r))
    OUT_lat = np.arctan((IN_z*(1-f)+e**2 * GEN_RAD_EARTH_EQ * np.sin(mu)**3) / ((1-f)*(d-e**2 * GEN_RAD_EARTH_EQ * np.cos(mu)**3)))
    
    # 4 - Compute heights
    OUT_height = (d*np.cos(OUT_lat)) + (IN_z*np.sin(OUT_lat)) - GEN_RAD_EARTH_EQ*np.sqrt(1-(e**2)*np.sin(OUT_lat)**2)
    
    return rad2deg(OUT_lon), rad2deg(OUT_lat), OUT_height


#######################################


if __name__ == '__main__':
    
    a = np.arange(10)
    a[1] = -1000
    a[9] = 1000
    print(a)
    print(computeMean_2sigma(a))
