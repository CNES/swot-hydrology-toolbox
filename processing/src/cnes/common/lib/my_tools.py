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
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: my_tools.py
    :synopsis: library with generic functions
     Created on 2017/03/02

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""

import datetime
from functools import partial
import logging
import math
import numpy as np
import os
from osgeo import osr, ogr
import pyproj
from scipy.ndimage.measurements import label
from shapely.ops import transform
from skimage.morphology import square
from sklearn.cluster import KMeans

import skimage
if skimage.__version__ >= "0.11":
    from skimage.filters import median as median_filter
else:
    from skimage.filter.rank import median as median_filter
from scipy.spatial import distance

import cnes.common.service_error as service_error

import cnes.common.lib.my_variables as my_var


def test_file(in_file, in_extent=None):
    """
    Test if full path in input is an existing file and, optionnally, a file in the awaited format

    :param in_file: input full path
    :type in_file: string
    :param in_extent: awaited format file (optionnal)
    :type in_extent: string
    """
    logger = logging.getLogger("my_tools")

    if os.path.exists(in_file):
        if not os.path.isfile(in_file):
            message = "ERROR = %s is not a file" % in_file
            raise service_error.ProcessingError(message, logger)
        if in_extent and not in_file.endswith(in_extent):
            message = "ERROR = %s must have the %s extent" % (in_file, in_extent)
            raise service_error.ProcessingError(message, logger)
    else:
        message = "ERROR = %s doesn't exist" % in_file
        raise service_error.ProcessingError(message, logger)


def test_list_of_files(in_list_file, in_extent=None):
    """
    Test if list of full paths in input is an existing file and, optionnally, a file in the awaited format

    :param in_file: list of input full path separed with ";"
    :type in_file: string
    :param in_extent: awaited format file (optionnal)
    :type in_extent: string
    """

    for file in in_list_file:
        test_file(file, in_extent)


def test_dir(in_dir):
    """
    Test if input path is an existing directory

    :param in_dir: input path
    :type in_dir: string
    """
    logger = logging.getLogger("my_tools")

    if os.path.exists(in_dir):
        if not os.path.isdir(in_dir):
            message = "ERROR = %s is not a directory" % in_dir
            raise service_error.ProcessingError(message, logger)
    else:
        message = "ERROR = %s doesn't exist" % in_dir
        raise service_error.ProcessingError(message, logger)


#######################################


def test_key(in_dict, in_key):
    """
    Test existence of in_key as a key of input dictionary in_dict
    
    :param in_dict: input dictionary
    :type in_dict: dict
    :param in_key: key to test presence in dictionary
    :type in_key: str
    
    :return: output_value = in_key if it's is a key of input dictionary in_dict; = None if not
    :rtype: str
    """
    output_value = in_key
    if in_key not in in_dict.keys():
        output_value = None
    return output_value


def get_value(in_dict, in_key):
    """
    Return value of in_key if in_key is a key of input dictionary in_dict; else return None
    
    :param in_dict: input dictionary
    :type in_dict: dict
    :param in_key: key to get value of
    :type in_key: str
    
    :return: output_value = value related to in_key if it's is a key of input dictionary in_dict; = None if not
    :rtype: str
    """
    output_value = None
    if (in_key is not None) and (in_key in in_dict.keys()):
        output_value = in_dict[in_key]
    return output_value


#######################################


def deg2rad(in_deg):
    """
    Convert angles in degrees to radians

    :param in_deg: angles to convert
    :type in_deg: scalar or 1D-array

    :return: angles in radians
    :rtype: same as input
    """
    return in_deg * np.pi / 180.


def rad2deg(in_rad):
    """
    Convert angles in radians to degrees

    :param in_rad: angles to convert
    :type in_rad: scalar or 1D-array

    :return: angles in degrees
    :rtype: same as input
    """
    return in_rad * 180. / np.pi


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
    if np.iterable(in_long) :
        ind = np.where(in_long > 180.0)[0]
        if len(ind) > 0:
            out_long[ind] -= 360.
    else :
        if in_long > 180.0:
            out_long -= 360
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
    if np.iterable(in_long):
        ind = np.where(in_long < 0.0)
        if ind is not None:
            out_long[ind] += 360.
    else :
        if in_long < 0.0:
            out_long += 360.
    return out_long


#######################################


def swot_timeformat(in_datetime, in_format=0):
    """
    Convert time into appropriate string format

    :param in_datetime: time value to write as a string
    :type in_datetime: datetime.datetime
    :param in_format: string format option 0 (default)="YYYY-MM-DD hh:mm:ss.ssssssZ" 1="YYYY-MM-DD hh:mm:ss"
                                            2="YYYY-MM-DDThh:mm:ss.ssssssZ" 3="YYYY-MM-DDThh:mm:ss"
    :type in_format: int

    :return: out_datetime: time value written as a string
    :rtype: string
    """

    out_datetime = ""

    if in_format == 1:
        out_datetime = in_datetime.strftime("%Y-%m-%d %H:%M:%S")
    elif in_format == 2:
        out_datetime = "%sZ" % in_datetime.strftime("%Y-%m-%dT%H:%M:%S.%f")
    elif in_format == 3:
        out_datetime = in_datetime.strftime("%Y-%m-%dT%H:%M:%S")
    else:
        out_datetime = "%sZ" % in_datetime.strftime("%Y-%m-%d %H:%M:%S.%f")

    return out_datetime


def convert_sec_2_time(in_time, txt_format=1):
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


def convert_utc_to_str(in_utc_time):
    """
    Convert UTC time to appropriate string format

    :param in_utc_time: date_time in seconds from 01/01/2000 00:00:00
    :type in_utc_time: float

    :return: UTC time from 01/01/2000 00:00:00 as a string
    :rtype: string
    """
    return swot_timeformat(datetime.datetime(2000, 1, 1) + datetime.timedelta(seconds=in_utc_time), in_format=3)


#######################################


def convert_fillvalue(in_data, in_flag="nc2shp"):
    """
    Convert NetCDF fill value in data into shapefile fill value (or reverse)

    :param in_data: data in which fill values will be replaced
    :type in_data: scalar or numpy.array
    :param in_flag: flag to specify conversion "nc2shp"(default)=NetCDF to Shapefile "shp2nc"=Shapefile to NetCDF
    :type in_flag: string

    :return: out_data = data in which fill values have been replaced
    :rtype: same as in_data in input
    """

    # Output vector init
    out_data = in_data

    if np.isscalar(in_data):

        # PROCESS FOR SCALAR

        # 1 - Retrieve data type
        type_name = type(in_data)

        # 2 - Get associated fill values
        fv_nc = my_var.FV_NETCDF[type_name]  # For NetCDF files
        fv_shp = my_var.FV_SHP[type_name]  # For Shapefiles

        # 3 - Make conversion
        if in_flag == "nc2shp":
            if in_data == fv_nc:
                out_data = fv_shp
        elif in_flag == "shp2nc":
            if in_data == fv_shp:
                out_data = fv_nc

    else:

        # PROCESS FOR ARRAYS

        # 1 - Retrieve data type
        type_name = in_data.dtype.name
        if type_name.startswith("str"):
            type_name = "str"

        # 2 - Get associated fill values
        fv_nc = my_var.FV_NETCDF[type_name]  # For NetCDF files
        fv_shp = my_var.FV_SHP[type_name]  # For Shapefiles

        # 3 - Make conversion
        if in_flag == "nc2shp":
            nan_idx = np.where(in_data == fv_nc)
            # Otherwise, _FillValues don't work
            if (type_name == 'uint8') or (type_name == 'uint16') or (type_name == 'int8') or (type_name == 'int16'): 
                out_data = in_data.astype('int32')
            out_data[nan_idx] = fv_shp
        elif in_flag == "shp2nc":
            nan_idx = np.where(in_data == fv_shp)
            out_data[nan_idx] = fv_nc

    return out_data


#######################################


def compute_bin_mat(in_size_x, in_size_y, in_x, in_y, verbose=True):
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
    :param verbose: print log debug flag
    :type verbose: bool

    :return: 2D matrix with "1" for each (in_x_i, in_y_i) and 0 elsewhere
    :rtype: 2D binary matrix of int 0/1
    """
    logger = logging.getLogger("my_tools")
    if verbose:
        logger.debug("- start -")

    # 0 - Deal with exceptions
    # 0.1 - Input vectors size must be the same
    if in_x.size != in_y.size:
        message = "compute_bin_mat(in_x, in_y) : in_x and in_y must be the same size ; currently : in_x = %d and in_y = %d" % (in_x.size, in_y.size)
        raise service_error.ProcessingError(message, logger)
    else:
        nb_pts = in_x.size
    if verbose:
        logger.debug("> Nb pixels to deal with = %d" % nb_pts)
    # 0.2 - max(X) < in_size_x
    if np.max(in_x) >= in_size_x:
        message = "compute_bin_mat(in_x, in_y) : elements of in_x must be less than in_size_x"
        raise service_error.ProcessingError(message, logger)
    # 0.3 - max(X) < in_size_x
    if np.max(in_y) >= in_size_y:
        message = "compute_bin_mat(in_x, in_y) : elements of in_y must be less than in_size_y"
        raise service_error.ProcessingError(message, logger)
    # 1 - Init output binary image
    out_bin_im = np.zeros((in_size_y, in_size_x))
    if verbose:
        logger.debug("> Binary matrix size = (X=%d , Y=%d)" % (in_size_x, in_size_y))

    # 2 - Put 1 for every pixels defined by the input vectors
    for ind in range(nb_pts):
        out_bin_im[in_y[ind], in_x[ind]] = 1

    return out_bin_im


#######################################


def label_region(in_bin_mat):
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
    logger.debug("- start -")

    out_m_label_regions, out_nb_obj = label(in_bin_mat)
    message = "> Number of detected objects = %d" % out_nb_obj
    logger.debug(message)

    return out_m_label_regions, out_nb_obj


def relabel_lake_using_segmentation_heigth(in_x, in_y, in_height, in_std_height_max):
    """
    This function main interest is to determine the number of lakes inside a subset of PixC in radar geometry.
    In most of the cases, only one lake is located in the given subset, but in some case of different lakes can
    be gather inside one because of radar geometric distortions or in the case a of a dam.

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
    :param in_std_height_max: maximal standard deviation of height inside a lake
    :type in_std_height_max: float

    :return: labels recomputed over 1 of several lakes
    :rtype: 1D vector of int
    """
    # Get instance of service config file
    logger = logging.getLogger("my_tools")
    logger.debug("- start -")

    # 0 - Deal with exceptions
    # 0.1 - Input vectors size must be the same
    if (in_x.size != in_y.size) or (in_x.size != in_height.size):
        raise ValueError("relabel_lake_using_segmentation_heigth(in_x, in_y, in_height) : in_x and in_y must be the same size ;" + \
                         "currently : in_x = %d and in_y = %d" % (in_x.size, in_y.size))
    else:
        nb_pts = in_x.size

    # 1 - Compute height matrix
    # 1.1 - Init heigth image
    logger.debug("Building heigth matrix")
    height_img = np.zeros((np.max(in_y) + 1, np.max(in_x) + 1))
    message = "> Height matrix size = (X=%d , Y=%d)" % (np.max(in_y.size) + 1, np.max(in_x.size) + 1)
    logger.debug(message)

    # 1.2 - Put height for every pixels defined by the input vectors
    for ind in range(nb_pts):
        height_img[in_y[ind], in_x[ind]] = in_height[ind]

    logger.debug("K-means processing")

    # 2 - Unsupervised clustering to determine number of classes
    std_heigth = np.std(in_height)
    nb_classes = 1

    # 2.1 - Cas with only one class
    if std_heigth < in_std_height_max:  # If only one lake is in the given pixc subset
        retour = np.ones(nb_pts)  # Return one unique label for all pixc
    else:
        # 2.2 - Init a k-means classifier
        kmeans_classif = KMeans()
        while std_heigth > in_std_height_max:

            nb_classes += 1

            # 2.3 - Cluster height over nb_classes classes
            kmeans_classif = KMeans(n_clusters=nb_classes)
            kmeans_classif.fit(in_height.reshape(-1, 1))  # Turn line into column for in_height

            # 2.3 - Compute heigth std inside each class
            std_heigth = np.max([np.std(in_height[np.where(kmeans_classif.labels_ == curLabel)]) for curLabel in np.unique(kmeans_classif.labels_)])

            # If number of classes upper than 10 => stop iteration
            if nb_classes > 10:
                break

        message = "NB classes : %d, max std : %f " % (nb_classes, std_heigth)
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

        retour = out_labels

    return retour


#######################################


def convert_2d_mat_in_1d_vec(in_x, in_y, in_mat):
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
    logger.debug("- start -")

    # Number of pixels in in_x
    nb_pix = in_x.size
    # Init output vector (same size of input in_x and in_y)
    out_vector = np.zeros(nb_pix)
    # Fill the output vector
    for indp in range(nb_pix):
        out_vector[indp] = in_mat[in_y[indp], in_x[indp]]

    return out_vector


#######################################


def compute_mean_2sigma(in_v_val, in_nan=None):
    """
    Compute the mean of the input values, after remove of non-sense values, i.e. below or above the median +/- 2*standard deviation
    If set, remove NaN values from the computation

    :param in_v_val: vector of values for which the mean is computed
    :type in_v_val: 1D-array of float
    :param in_nan: value to consider as NaN (default=None)
    :type in_nan: depends on the initial vector

    :return: the mean value (None if no finite value)
    :rtype: float
    """
    logger = logging.getLogger("my_tools")

    retour = None  # Init retour to None
    flag_ok = True  # Flag to compute the sum
    v_val = in_v_val  # Vector on which compute the sum
    
    # 0 - Consider only non-NaN values
    if in_nan is not None:
        not_nan_idx = np.where(v_val < in_nan)[0]
        if len(not_nan_idx) == 0:
            flag_ok = False
            retour = None
        else:
            v_val = in_v_val[not_nan_idx]

    if flag_ok:
        
        # 1 - Compute statistical values over this vector
        # 1.1 - Median
        med = np.median(v_val)
        # 1.2 - Standard deviation
        std = np.std(v_val)
        # 1.3 - Nb values
        nb_val = v_val.size

        # 2 - Remove indices with value out of 2 sigma

        # 2.1 - Lower than mean - 2 sigma
        idx_lower = np.where(v_val < med-2*std)[0]
        v_indices = np.delete(np.arange(nb_val), idx_lower)
        # 2.2 - Higher than mean + 2 sigma
        idx_higher = np.where(v_val[v_indices] > med+2*std)[0]
        v_indices = np.delete(v_indices, idx_higher)

        message = "%d / %d pixels used for mean computation" % (v_indices.size, nb_val)
        logger.debug(message)

        # 3 - Return the mean of clean vector
        retour = np.mean(v_val[v_indices])

    return retour


def compute_std(in_v_val, in_nan=None):
    """
    Compute the standard deviation of the input values
    If set, remove NaN values from the computation

    :param in_v_val: vector of values for which the standard deviation is computed
    :type in_v_val: 1D-array of float
    :param in_nan: value to consider as NaN (default=None)
    :type in_nan: depends on the initial vector

    :return: the mean value (None if no finite value)
    :rtype: float
    """

    retour = None  # Init retour to None
    flag_ok = True  # Flag to compute the sum
    v_val = in_v_val  # Vector on which compute the sum
    
    # Consider only non-NaN values
    if in_nan is not None:
        not_nan_idx = np.where(v_val < in_nan)[0]
        if len(not_nan_idx) == 0:
            flag_ok = False
            retour = None
        else:
            v_val = in_v_val[not_nan_idx]

    # Return the standard deviation of the clean vector elements
    if flag_ok:
        retour = np.std(v_val)
        
    return retour


def compute_sum(in_v_val, in_nan=None):
    """
    Compute the sum of the input values
    If set, remove NaN values from the computation

    :param in_v_val: vector of values for which the sum is computed
    :type in_v_val: 1D-array of float
    :param in_nan: value to consider as NaN (default=None)
    :type in_nan: depends on the initial vector

    :return: the mean value (None if no finite value)
    :rtype: float
    """
    
    retour = None  # Init retour to None
    flag_ok = True  # Flag to compute the sum
    v_val = in_v_val  # Vector on which compute the sum
    
    # Consider only non-NaN values
    if in_nan is not None:
        not_nan_idx = np.where(v_val < in_nan)[0]
        if len(not_nan_idx) == 0:
            flag_ok = False
            retour = None
        else:
            v_val = in_v_val[not_nan_idx]

    # Return the sum of the clean vector elements
    if flag_ok:
        retour = np.sum(v_val)
        
    return retour


#######################################


def compute_az(in_lon, in_lat, in_v_nadir_lon, in_v_nadir_lat):
    """
    Compute the azimuth index associated to the point with coordinates P(in_lon, in_lat) with respect
    to the nadir track with coordinates (in_v_nadir_lon, in_v_nadir_lat) for each point.
    NB: the method used is to minimize the scalar product between P and the points of the nadir track.

    :param in_lon: longitude in degrees east
    :type in_lon: float
    :param in_lat: latitude in degrees north
    :type in_lat: float
    :param in_v_nadir_lon: longitude of each nadir point in degrees east
    :type in_v_nadir_lon: 1D-array
    :param in_v_nadir_lat: latitude of each nadir point in degrees north
    :type in_v_nadir_lat: 1D-array

    :return: out_idx_min azimuth index corresponding to input point (in_lon, in_lat)
    :type: int
    """

    # 1 - Stack data
    in_coords = np.vstack((np.array(in_lon), np.array(in_lat)))
    nadir_coords = np.vstack((np.array(in_v_nadir_lon), np.array(in_v_nadir_lat)))

    # 2 - Compute distances between point (in_lon, in_lat) and vector (in_v_nadir_lon, in_v_nadir_lat)
    distances = (distance.cdist(in_coords.transpose(), nadir_coords.transpose())).transpose()
    # get indice of minimum distance
    out_idx_min = np.argmin(distances)

    # 3 - Return corresponding azimuth index
    return out_idx_min


#######################################


def compute_dist(in_lon1, in_lat1, in_lon2, in_lat2):
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

    :return: the spherical distance (in meters) between P1 and P2
    :type: float
    """

    # 1 - Convert degrees to radians
    lon1 = deg2rad(in_lon1)
    lat1 = deg2rad(in_lat1)
    lon2 = deg2rad(in_lon2)
    lat2 = deg2rad(in_lat2)

    # 1 - Distance angulaire en radians
    s12 = np.arccos(np.sin(lat1)*np.sin(lat2) + np.cos(lat1)*np.cos(lat2)*np.cos(lon2-lon1))

    return s12 * my_var.GEN_APPROX_RAD_EARTH


#######################################


def get_utm_epsg_code(lon, lat):
    """
    get EPSG code
    :param lon: longitude of P1 in degrees east
    :type lon: float
    :param lat: latitude of P1 in degrees north
    :type lat: float

    :return: ESPG code
    :type: int
    """
    utm_band = str((math.floor((lon + 180 )/6)%60)+1)
    if len(utm_band) == 1:
        utm_band = '0'+utm_band
    if lat >= 0:
        epsg = '326' + utm_band
    else:
        epsg = '326' + utm_band
    return epsg


def get_utm_coords(in_lon, in_lat):
    """
    get UTM coordinates
    :param in_lon: longitude of P1 in degrees east
    :type in_lon: float
    :param in_lat: latitude of P1 in degrees north
    :type in_lat: float

    :return: UTM coordinates
    :type: float
    """
    lat_mean = np.mean(in_lat)
    lon_mean = np.mean(in_lon)

    # code commenté x_c, y_c, zone_number, zone_lettre = utm.from_latlon(lat_mean, lon_mean)
    latlon = pyproj.Proj(init="epsg:4326")
    epsg=get_utm_epsg_code(lon_mean, lat_mean)
    utm_proj = pyproj.Proj(init='epsg:' + epsg)
    X, Y = pyproj.transform(latlon, utm_proj, in_lon, in_lat)
    return (X,Y, )


def get_area(in_polygon, centroid=None):
    """
    This function projects in_polygon from geographic coordinate into centroid corresponding UTM zone.
    Once projected, the area of the polygon is computed and returned

    :param in_polygon: polygon
    :type in_polygon: OGR geometry
    :param centroid: centroid of in_polygon (optional)
    :type centroid: tuple of coordinates

    :return: area of in_polygon
    :rtype: float
    """
    
    # 0 - Copy polygon
    tmp_poly = in_polygon.Clone()

    # 1 - Computation of EPSG code corresponding UTM zone
    epsg_code = "32"  # Initialization
    # Retrieve centroid coordinates
    if centroid is None:
        poly_centroid = tmp_poly.Centroid().GetPoint(0)
        centroid_lon = poly_centroid[0]
        centroid_lat = poly_centroid[1]
    else:
        centroid_lon = centroid[0]
        centroid_lat = centroid[1]

    epsg_code = get_utm_epsg_code(centroid_lon, centroid_lat)

    # 1 - Projection of tmp_poly into UTM
    src_source = osr.SpatialReference()
    src_source.ImportFromEPSG(4326)

    src_target = osr.SpatialReference()
    src_target.ImportFromEPSG(int(epsg_code))

    transform = osr.CoordinateTransformation(src_source, src_target)

    tmp_poly.Transform(transform)

    # 2 - Compute and return area
    return tmp_poly.GetArea()


#######################################


def get_utm_coords_from_lonlat(in_long, in_lat, utm_epsg_code=None):
    """
    This function projects coordinates from geographic coordinate into UTM zone corresponding to the mean of in_lon, in, lat.

    :param in_long: longitudes of points
    :type in_long: 1D-array of float
    :param in_lat: latitudes of points
    :type in_lat: 1D-array of float

    :return: X coordinates, Y coordintes, UTM_epsg_code
    :rtype: tuple of (1D-array of float, 1D-array of float, int)
    """

    if utm_epsg_code == None :
        lat_mean = np.mean(in_lat)
        lon_mean = np.mean(in_long)
        utm_epsg_code = get_utm_epsg_code(lon_mean, lat_mean)

    latlon_proj = pyproj.Proj(init="epsg:4326")
    utm_proj = pyproj.Proj(init='epsg:' + utm_epsg_code)

    X, Y = pyproj.transform(latlon_proj, utm_proj, in_long, in_lat)

    return X, Y, utm_epsg_code


def get_lon_lat_polygon_from_utm(in_utm_poly, in_utm_epsg_code):
    """
    This function projects polygon from UTM projection given in in_utm_epsg_code into lon lat coordinates.

    :param in_utm_poly: utm polygon
    :type in_utm_poly: shapely Polygon or Multipolygon geometry
    :param in_utm_epsg_code: code of utm zone
    :type in_utm_epsg_code: int

    :return: polygon in geographical lon lat coordinates
    :rtype: shapely Polygon or Multipolygon geometry
    """
    utm_proj = pyproj.Proj(init='epsg:' + in_utm_epsg_code)
    latlon_proj_wrap = pyproj.Proj(init="epsg:4326")
    projection_wm_func = partial(pyproj.transform, utm_proj, latlon_proj_wrap)

    lonlat_poly = transform(projection_wm_func, in_utm_poly)

    return lonlat_poly


#######################################


def load_pixels_to_mem_layer(in_lon, in_lat):
    """
    Write pixels given their longitude and latitude coordinates to a memory layer
    
    :param in_lon: longitudes of pixels
    :type in_lon: 1D-array of float
    :param in_lat: latitudes of pixels
    :type in_lat: 1D-array of float
    
    :return out_data_source: data source of output layer
    :rtype out_data_source: OGRdata_source
    :return out_layer: output layer
    :rtype out_layer: OGRlayer
    """
    
    # 1 - Create memory layer
    mem_driver = ogr.GetDriverByName('MEMORY')  # Memory driver
    # Open the memory datasource with write access
    out_data_source = mem_driver.CreateDataSource('memData')
    # Set spatial projection
    srs = ogr.osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    # Create output layer
    out_layer = out_data_source.CreateLayer(str('point'), srs=srs, geom_type=ogr.wkbPoint)

    # 2 - Add pixels to layer
    for (lon, lat) in zip(in_lon, in_lat):
        
        # 2.1 - Convert current (lon,lat) to point geometry
        p = ogr.Geometry(ogr.wkbPoint)
        p.AddPoint(lon, lat)
        
        # 2.2 - Init point geometry wrt layer definition and set it to point
        point = ogr.Feature(out_layer.GetLayerDefn())
        point.SetGeometry(p)
        
        # 2.3 - Add point to layer
        out_layer.CreateFeature(point)

    # 3 - Return layer and data source
    return out_layer, out_data_source


def get_layer_fields_name_and_type(in_layer):
    """
    Return name and type of each field of a layer
    
    :param in_layer: layer
    :type in_layer: OGRlayer
    
    :return: out_list_field_type = list of type of each field
    :rtype: dict
    """
    
    # 0 - Init output variable = list of field type
    out_list_type = {}  # Init list of fields
    
    # 1 - Retrieve layer definition
    layer_defn = in_layer.GetLayerDefn()
    
    # 2 - Store field type for each field name
    for att in range(layer_defn.GetFieldCount()) :
        # 2.1 - Get attribute name
        field_name = layer_defn.GetFieldDefn(att).GetName()
        # 2.2 - Get attribute type
        field_type = layer_defn.GetFieldDefn(att).GetFieldTypeName(layer_defn.GetFieldDefn(att).GetType())
        # 2.3 - Save in dictionary 
        out_list_type[field_name] = field_type
        
    return out_list_type
