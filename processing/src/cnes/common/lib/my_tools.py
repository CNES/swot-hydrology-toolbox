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

from cnes.common.lib.my_variables import GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE, GEN_APPROX_RAD_EARTH
import cnes.common.lib_lake.locnes_variables as my_var
import cnes.common.serviceError as serviceError


def testFile(in_file, IN_extent=None):
    """
    Test if full path in input is an existing file and, optionnally, a file in the awaited format
    
    :param in_file: input full path
    :type in_file: string
    :param IN_extent: awaited format file (optionnal)
    :type IN_extent: string
    """
    logger = logging.getLogger("my_tools")
    
    if os.path.exists(in_file):
        if not os.path.isfile(in_file):
            message = "ERROR = %s is not a file" % in_file
            raise serviceError.ProcessingError(message, logger)
        if IN_extent and not in_file.endswith(IN_extent):
            message = "ERROR = %s must have the %s extent" % (in_file, IN_extent)
            raise serviceError.ProcessingError(message, logger)
    else:
        message = "ERROR = %s doesn't exist" % in_file
        raise serviceError.ProcessingError(message, logger)


def testDir(in_dir):
    """
    Test if input path is an existing directory
    
    :param in_dir: input path
    :type in_dir: string
    """
    logger = logging.getLogger("my_tools")
    
    if os.path.exists(in_dir):
        if not os.path.isdir(in_dir):
            message = "ERROR = %s is not a directory" % in_dir
            raise serviceError.ProcessingError(message, logger)
    else:
        message = "ERROR = %s doesn't exist" % in_dir
        raise serviceError.ProcessingError(message, logger)


#######################################


def deg2rad(in_deg):
    """
    Convert angles in degrees to radians
    
    :param in_deg: angles to convert
    :type in_deg: scalar or 1D-array
    
    :return angles in radians
    :rtype: same as input
    """
    return in_deg * np.pi / 180.


def rad2deg(in_rad):
    """
    Convert angles in radians to degrees
    
    :param in_rad: angles to convert
    :type in_rad: scalar or 1D-array
    
    :return angles in degrees
    :rtype: same as input
    """
    return in_rad * 180. / np.pi


#######################################


def computeBinMat(in_size_x, in_size_y, in_x, in_y):
    """
    Creates a 2D binary matrix from Y and Y 1D vectors
    i.e. for each ind: mat[in_x[ind],in_y[ind]] = 1
        
    :param in_size_x: number of pixels in X dimension
    :type in_size_x: int
    :param in_size_y: number of pixels in Y dimension
    :type in_size_y: int
    :param in_x: X indices of "1" pixels
    :type in_x: 1D vector of int
    :param in_y: Y indices of "1" pixels
    :type in_y: 1D vector of int
        
    :return: 2D matrix with "1" for each (in_x_i, in_y_i) and 0 elsewhere
    :rtype: 2D binary matrix of int 0/1
    """
    logger = logging.getLogger("my_tools")
    logger.debug("[my_tools] == computeBinMat ==")
    
    # 0 - Deal with exceptions
    # 0.1 - Input vectors size must be the same
    if in_x.size != in_y.size:
        message = "computeBinMat(in_x, in_y) : in_x and in_y must be the same size ; currently : in_x = %d and in_y = %d" % (in_x.size, in_y.size)
        raise serviceError.ProcessingError(message, logger)
    else:
        nb_pts = in_x.size
    logger.debug("> Nb pixels to deal with = %d" % nb_pts)
    # 0.2 - max(X) < in_size_x
    if np.max(in_x) >= in_size_x:
        message = "computeBinMat(in_x, in_y) : elements of in_x must be less than in_size_x"
        raise serviceError.ProcessingError(message, logger)
    # 0.3 - max(X) < in_size_x
    if np.max(in_y) >= in_size_y:
        message = "computeBinMat(in_x, in_y) : elements of in_y must be less than in_size_y"
        raise serviceError.ProcessingError(message, logger)
    # 1 - Init output binary image
    out_bin_im = np.zeros((in_size_y, in_size_x))
    logger.debug("> Binary matrix size = (X=%d , Y=%d)" % (in_size_x, in_size_y))
        
    # 2 - Put 1 for every pixels defined by the input vectors
    for ind in range(nb_pts): 
        out_bin_im[in_y[ind], in_x[ind]] = 1
    
    return out_bin_im


#######################################


def labelRegion(in_bin_mat):
    """
    Identifies all separate regions in a 2D binary matrix

    :param in_bin_mat: 2D binary matrix
    :type in_bin_mat: 2D binary matrix of int 0/1
    
    :return: a 2D matrix, in which each pixel has a unique value ranging from 0 to N (number of identified regions)
    :rtype: 2D matrix of int 
    :return: the number of objects
    :rtype: int
    """
    logger = logging.getLogger("my_tools")
    logger.debug("[my_tools] == labelRegion ==")
    
    out_m_label_regions, out_nb_obj = label(in_bin_mat)
    message = "> Number of detected objects = %d" % out_nb_obj
    logger.debug(message)
    
    return out_m_label_regions, out_nb_obj


def relabelLakeUsingSegmentationHeigth(in_x, in_y, in_height):
    """
    This function main interest is to determine the number of lakes inside a subset of PixC in radar geometry.
    In most of the cases, only one lake is located in the given subset, but in some case of different lakes can be gather inside one because of radar geometric distortions or in the case a of a dam.
    Steps :
        1. Creates a 2D height matrix from X and Y 1D vectors
        2. Unsupervised heigth classification to determine number a classes to get an STD over each classe lower than STD_HEIGHT_MAX
        3. Smooth result using a 2D median filter
        4. Return labels

    :param in_x: X indices of "1" pixels
    :type in_x: 1D vector of int
    :param in_y: Y indices of "1" pixels
    :type in_y: 1D vector of int
    :param in_height: height of pixels
    :type in_height: 1D vector of float

    :return: labels recomputed over 1 of several lakes
    :rtype: 1D vector of int
    """

    logger = logging.getLogger("my_tools")
    logger.debug("[my_tools] == relabelLakeUsingSegmentationHeigth ==")

    logger.debug("[my_tools] Building heigth matrix")

    # 0 - Deal with exceptions
    # 0.1 - Input vectors size must be the same
    if (in_x.size != in_y.size) or (in_x.size != in_height.size):
        raise ValueError("relabelLakeUsingSegmentationHeigth(in_x, in_y, in_height) : in_x and in_y must be the same size ; currently : in_x = %d and in_y = %d" % (in_x.size, in_y.size))
    else:
        nb_pts = in_x.size

    # 1 - Compute height matrix
    # 1.1 - Init heigth image
    height_img = np.zeros((np.max(in_y) + 1, np.max(in_x) + 1))
    message = "> Height matrix size = (X=%d , Y=%d)" % (np.max(in_y.size) + 1, np.max(in_x.size) + 1)
    logger.debug(message)

    # 1.2 - Put height for every pixels defined by the input vectors
    for ind in range(nb_pts):
        height_img[in_y[ind], in_x[ind]] = in_height[ind]

    logger.debug("[my_tools] K-means processing")

    # 2 - Unsupervised clustering to determine number of classes
    std_heigth = np.std(in_height)
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
        kmeans_classif.fit(in_height.reshape(-1, 1))  # Turn line into column for in_height

        # 2.3 - Compute heigth std inside each class
        std_heigth = np.max([np.std(in_height[np.where(kmeans_classif.labels_ == curLabel)]) for curLabel in np.unique(kmeans_classif.labels_)])

        # If number of classes upper than 10 => stop iteration
        if nb_classes > 10:
            break

    message = "[my_tools] NB classes : %d, max std : %f " % (nb_classes, std_heigth)
    logger.debug(message)

    # 3 - Format output vector

    # 3.1 - Build a labeled matrix
    labeled_img = np.zeros(height_img.shape).astype('int')
    for ind in range(nb_pts):
        labeled_img[in_y[ind], in_x[ind]] = int(kmeans_classif.labels_[ind])

    # 3.2 - Median filter on a 3x3 window to smooth output
    labeled_img_filted = median_filter(labeled_img.astype('uint8'), square(2).astype('uint8'))

    # 3.3 - Init and fill output label array
    out_labels = np.zeros(in_height.shape)
    for ind in range(nb_pts):
        out_labels[ind] = labeled_img_filted[in_y[ind], in_x[ind]] + 1

    return out_labels


#######################################


def convert2dMatIn1dVec(in_x, in_y, in_mat):
    """
    Converts a 2D matrix [ in_x[ind], in_y[ind] ] in a 1D vector [ind]
    
    :param in_x: column indices of points of in_mat to return
    :type in_x: 1D array
    :param in_y: row indices of points of in_mat to return
    :type in_y: 1D array
    :param in_mat: 2D array
    :type in_mat: 2D array of values (int, float, ...)
    
    :return: a 1D vector, each value at "ind" being the value of in_mat[ in_x[ind], in_y[ind] ]
    :rtype: 2D matrix of values (int, float, ...)
    """
    logger = logging.getLogger("my_tools")
    logger.debug("[my_tools] == convert2dMatIn1dVec ==")
    
    # Number of pixels in in_x
    nb_pix = in_x.size
    # Init output vector (same size of input in_x and in_y)   
    out_vector = np.zeros(nb_pix)
    # Fill the output vector
    for indp in range(nb_pix):
        out_vector[indp] = in_mat[in_y[indp], in_x[indp]]
    
    return out_vector


#######################################


def cptDigits(in_real):
    """
    Compute the position of the first significant figure of a real number (ex: 110 => 3 ; 0.000156 => -4)
        
    :param in_real: real number
    :type in_real: float
    
    :return the position of the first significant figure of a real number
    :rtype: int 
    """
    
    # Cas particulier de 0
    if in_real == 0:
        out_cpt = 1
    
    elif abs(in_real) >= 10:  # Cas des nombres >= 10
        out_cpt = 1 
        while abs(in_real) >= 10:
            in_real /= 10.0
            out_cpt += 1
            
    elif abs(in_real) < 1:  # Cas des nombres < 1
        out_cpt = 0
        while abs(in_real) < 1:
            in_real *= 10.0
            out_cpt -= 1
            
    else:
        out_cpt = 1
        
    return out_cpt


#######################################


def convertSec2Time(in_time, txt_format=1):
    """
    Convert a date_time in seconds towards a string with the format specified as in put parameter
    
    :param in_time: date_time in seconds
    :type in_time: float
    :param txt_format: [optionnal] 1="hh:mm:ss" 2="[hh]h [mm]min [ss]s" 3="hh:mm:ss.ms" 4="Day DD hh:mm:ss"
    :type txt_format: int (defaut = 1)
    
    :return: la duree correspondant a in_time au format choisi
    :rtype: str
    """
    
    # Calcul des secondes
    tmp_mm = math.floor(in_time / 60.0)
    SS = in_time - 60.0 * tmp_mm
        
    # Calcul des minutes
    HH = math.floor(tmp_mm / 60.0)
    MM = tmp_mm - 60.0 * HH
        
    # Calcul des heures si > 24h
    DD = math.floor(HH / 24.0)
    HH -= 24.0 * DD
    if txt_format == 1:
        retour = "%02d:%02d:%02d" % (HH, MM, SS)
    elif txt_format == 2:
        retour = "%02dh%02dmin%02ds" % (HH, MM, SS)
    elif txt_format == 3:
        retour = "%02d:%02d:%02.3f" % (HH, MM, SS)
    elif txt_format == 4:
        retour = "Day %02d %02d:%02d:%02d" % (DD+1, HH, MM, SS)  # DD+1 car le 1e jour est le jour 1

    return retour

#######################################


def computeMean_2sigma(in_v_val):
    """
    Compute the mean of the input values, after remove of non-sense values, i.e. below or above the median +/- 2*standard deviation
    
    :param in_v_val: vector of values for which the mean is computed
    :type in_v_val: 1D-array of float
        
    :return the mean value
    :rtype: float
    """
    
    logger = logging.getLogger("my_tools")
    # 1 - Compute statistical values over this vector
    # 1.1 - Median
    med = np.median(in_v_val)
    # 1.2 - Standard deviation
    std = np.std(in_v_val)
    # 1.3 - Nb values
    nb_val = in_v_val.size
    
    if (len(in_v_val.shape)) == 2:
        in_v_val = in_v_val[0]
    # 2 - Remove indices with value out of 2 sigma
    # 2.1 - Lower than mean - 2 sigma
    idx_lower = np.where(in_v_val < med-2*std)[0]
    v_indices = np.delete(np.arange(nb_val), idx_lower)
    # 2.2 - Higher than mean + 2 sigma

    idx_higher = np.where(in_v_val[v_indices] > med+2*std)[0]
    v_indices = np.delete(v_indices, idx_higher)
    
    message = "[my_tools/computeMean_2sigma] %d / %d pixels used for mean computation" % (v_indices.size, nb_val)
    logger.debug(message)
        
    # 3 - Return the mean of cleaned vector
    return np.mean(in_v_val[v_indices])


#######################################


def computeAz(in_lon, in_lat, in_v_nadir_lon, in_v_nadir_lat):
    """
    Compute the azimuth index associated to the point with coordinates P(in_lon, in_lat) with respect to the nadir track with coordinates (in_v_nadir_lon, in_v_nadir_lat) for each point.
    NB: the method used is to minimize the scalar product between P and the points of the nadir track.
    
    :param in_lon: longitude in degrees east
    :type in_lon: float
    :param in_lat: latitude in degrees north
    :type in_lat: float
    :param in_v_nadir_lon: longitude of each nadir point in degrees east
    :type in_v_nadir_lon: 1D-array
    :param in_v_nadir_lat: latitude of each nadir point in degrees north
    :type in_v_nadir_lat: 1D-array
    
    :return out_idx_min: azimuth index corresponding to input point (in_lon, in_lat)
    :type: int
    """
    
    # 1 - Init variables
    nb_nadir = in_v_nadir_lon.size
    list_scalar = np.zeros(nb_nadir)
    
    # 2 - Compute scalar product for each nadir point
    # 2.1 - 1st point
    list_scalar[0] = (in_v_nadir_lon[0]-in_lon)*(in_v_nadir_lon[1]-in_v_nadir_lon[0]) + (in_v_nadir_lat[0]-in_lat)*(in_v_nadir_lat[1]-in_v_nadir_lat[0])
    # 2.2 - Middle points
    for ind in np.arange(1, nb_nadir-1):
        list_scalar[ind] = (in_v_nadir_lon[ind]-in_lon)*(in_v_nadir_lon[ind+1]-in_v_nadir_lon[ind-1]) + (in_v_nadir_lat[ind]-in_lat)*(in_v_nadir_lat[ind+1]-in_v_nadir_lat[ind-1])
    # 2.3 - Last point
    list_scalar[-1] = (in_v_nadir_lon[-1]-in_lon)*(in_v_nadir_lon[-1]-in_v_nadir_lon[-2]) + (in_v_nadir_lat[-1]-in_lat)*(in_v_nadir_lat[-1]-in_v_nadir_lat[-2])
    
    # 3 - Find min scalar
    out_idx_min = np.argmin(np.absolute(list_scalar))
    
    # 4 - Return corresponding azimuth index
    return out_idx_min


#######################################


def computeDist(in_lon1, in_lat1, in_lon2, in_lat2):
    """
    Compute Euler distance between 2 points of the plan: P1(lon1, lat1) and P2(lon2, lat2)
    
    :param in_lon1: longitude of P1 in degrees east
    :type in_lon1: float
    :param in_lat1: latitude of P1 in degrees north
    :type in_lat1: float
    :param in_lon2: longitude of P2 in degrees east
    :type in_lon2: float
    :param in_lat2: latitude of P2 in degrees north
    :type in_lat2: float
    
    :return the spherical distance (in meters) between P1 and P2
    :type: float
    """
    
    # 1 - Convert degrees to radians
    lon1 = deg2rad(in_lon1)
    lat1 = deg2rad(in_lat1)
    lon2 = deg2rad(in_lon2)
    lat2 = deg2rad(in_lat2)
    
    # 1 - Distance angulaire en radians
    s12 = np.arccos(np.sin(lat1)*np.sin(lat2) + np.cos(lat1)*np.cos(lat2)*np.cos(lon2-lon1))
    
    return s12 * GEN_APPROX_RAD_EARTH


#######################################


def getArea(in_polygon, in_centroid):
    """
    This function projects in_polygon from geographic coordinate into centroid corresponding UTM zone.
    Once projected, the area of the polygon is computed and returned
    
    :param in_polygon: polygon
    :type in_polygon: OGR geometry
    :param in_centroid: centroid of in_polygon
    :type in_centroid: tuple of coordinates
    
    :return: area of in_polygon
    :rtype: float
    """
    
    # 0 - Computation of EPSG code corresponding UTM zone
    epsg_code = "32"  # Initialization

    centroid_lon = in_centroid[0]
    centroid_lat = in_centroid[1]

    if centroid_lat > 0:
        epsg_code += "6"
    else:
        epsg_code += "7"

    epsg_code += str(int(1 + (centroid_lon + 180)//6))

    # 1 - Projection of in_polygon into UTM
    src_source = osr.SpatialReference()
    src_source.ImportFromEPSG(4326)

    src_target = osr.SpatialReference()
    src_target.ImportFromEPSG(int(epsg_code))

    transform = osr.CoordinateTransformation(src_source, src_target)

    in_polygon.Transform(transform)

    # 2 - Compute and return area
    return in_polygon.GetArea()


#######################################


def llh2xyz(in_lon, in_lat, in_height, IN_flag_rad=True):
    """
    Convert geographic coordinates (longitude, latitude, height) to cartesian coordinates (x, y, z) 
        
    :param in_lon: longitude in degrees east
    :type in_lon: scalar or 1D-array
    :param in_lat: latitude in degrees north
    :type in_lat: scalar or 1D-array
    :param in_height: height in meters
    :type in_height: scalar or 1D-array
    :param IN_flag_rad: =False if IN_lan and in_lat are in degrees, =True if they are in radians (default)
    :type IN_flag_rad: boolean
    
    :return out_x, out_y, out_z: cartesian coordinates
    :type: scalars or 1D-array
    
    .. author: Damien DESROCHES & Alejandro BOHE - CNES DSO/SI/TR
    """
    
    # 1 - Convert geodetic coordinates in radians
    if IN_flag_rad:
        lat = in_lat
        lon = in_lon
    else:
        lat = deg2rad(in_lat)
        lon = deg2rad(in_lon)
    
    # 2 - Compute earth excentricity
    e = np.sqrt(GEN_RAD_EARTH_EQ ** 2 - GEN_RAD_EARTH_POLE ** 2) / GEN_RAD_EARTH_EQ

    # 3 - Compute earth radius for latitude lat
    Rn = GEN_RAD_EARTH_EQ / np.sqrt(1 - e**2 * np.sin(lat)**2)
    
    # 4 - Compute cartesian coordinates
    out_x = (Rn + in_height) * np.cos(lat) * np.cos(lon)
    out_y = (Rn + in_height) * np.cos(lat) * np.sin(lon)
    out_z = (Rn * (1.-e**2) + in_height) * np.sin(lat)

    return out_x, out_y, out_z   


def xyz2llh(in_x, in_y, in_z):
    """
    Convert cartesian coordinates (x, y, z) to geographic coordinates (lat, lon, height) 
    
    :param in_x: coordinate along x-axis
    :type in_x: scalar or 1D-array
    :param in_y: coordinate along y-axis
    :type in_y: scalar or 1D-array
    :param in_z: coordinate along z-axis
    :type in_z: scalar or 1D-array
    
    :return out_long, out_lat, out_height: geographic coordinates (longitude, latitude, height) 
    :rtype: same as input (scalar or 1D-array)
    """
    
    # 1 - Compute ellipsoide variables
    e = np.sqrt(GEN_RAD_EARTH_EQ ** 2 - GEN_RAD_EARTH_POLE ** 2) / GEN_RAD_EARTH_EQ  # Ellipsoid excentricity
    f = 1 - np.sqrt(1 - e**2)  # Flattening factor
    r = np.sqrt(in_x**2 + in_y**2 + in_z**2)  # Distance between center of earth and point
    d = np.sqrt(in_x**2 + in_y**2)  # Distance between center of earth and point with latitude 0
    
    # 2 - Compute longitudes
    out_lon = np.arctan2(in_y, in_x)
    
    # 3 - Compute latitudes
    mu = np.arctan(in_z/d * ((1-f) + e**2 * GEN_RAD_EARTH_EQ/r))
    out_lat = np.arctan((in_z*(1-f)+e**2 * GEN_RAD_EARTH_EQ * np.sin(mu)**3) / ((1-f)*(d-e**2 * GEN_RAD_EARTH_EQ * np.cos(mu)**3)))
    
    # 4 - Compute heights
    out_height = (d*np.cos(out_lat)) + (in_z*np.sin(out_lat)) - GEN_RAD_EARTH_EQ*np.sqrt(1-(e**2)*np.sin(out_lat)**2)
    
    return rad2deg(out_lon), rad2deg(out_lat), out_height


#######################################


if __name__ == '__main__':
    
    a = np.arange(10)
    a[1] = -1000
    a[9] = 1000
    print(a)
    print(computeMean_2sigma(a))
