# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:3.1.0:DM:#91:2021/05/21:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: my_segmentation.py
    :synopsis: handle height-segmentation computation
     Created on 2021/01/14

.. moduleauthor:: Cécile CAZALS - CS

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""

import logging
import numpy as np
from numpy.lib.stride_tricks import as_strided

from scipy.cluster import hierarchy
from scipy.ndimage import median_filter
import scipy.ndimage as ndimage

from skimage.morphology import square
from skimage import filters
from skimage.segmentation import watershed
from skimage import segmentation
from skimage.feature import peak_local_max
import skimage
if skimage.__version__ >= "0.11":
    from skimage.filters import median as median_filter
else:
    from skimage.filter.rank import median as median_filter

from sklearn.cluster import KMeans, MeanShift, estimate_bandwidth

import cnes.common.lib.my_tools as my_tools
import cnes.common.service_error as service_error


def relabel_lake_using_segmentation_heigth(in_x, in_y, in_height, in_pix_area, min_size, height_segmentation_method):
    """
    This function main interest is to determine the number of lakes inside a subset of PixC in radar geometry.
    In most of the cases, only one lake is located in the given subset, but in some case of different lakes can
    be gather inside one because of radar geometric distortions or in the case a of a dam.
    Each resulting lake should have a area over min_size, otherwise, it is not split.

    Steps :

     1. Creates a 2D height matrix from X and Y 1D vectors
     2. Segmentation method define by height_segmentation_method
     3. Smooth result using a 2D median filter to filter sublakes with area < min_size
     4. Return labels

    :param in_x: X indices of pixels
    :type in_x: 1D vector of int
    :param in_y: Y indices of pixels
    :type in_y: 1D vector of int
    :param in_height: height of pixels
    :type in_height: 1D vector of float
    :param in_pix_area: area of pixels
    :type in_pix_area: 1D vector of float
    :param min_size: minimum size for object
    :type min_size: float
    :param height_segmentation_method: segmentation method
    :type height_segmentation_method: int

    :return: labels recomputed over 1 of several lakes
    :rtype: 1D vector of int
    """
    # Get instance of service config file
    logger = logging.getLogger("my_tools")

    nb_pts = in_x.size
    area = np.sum(in_pix_area)

    logger.debug("Number of pixels : %d ; area %d ; Min - Max height: %f - %f" %(nb_pts, int(area), np.min(in_height), np.max(in_height)))

    # split lake if contains more than 5 pixels or area >= 2*min_size
    compute_labels = True

    # if area < 2*min_size or contains less than 5 pixels, lake will not be split.
    if area < 2 * min_size or nb_pts < 5 :
        out_labels = np.ones((nb_pts), dtype='int')
        compute_labels = False

    elif height_segmentation_method == 1 :
        logger.debug("Lake segmentation method = Felzenszwalb")
        height_img = get_2d_from_1d(in_height, in_x, in_y)
        labeled_img = segmentation.felzenszwalb(height_img)
        out_labels = get_1d_from_2d(labeled_img, in_x, in_y)

    elif height_segmentation_method == 2 :
        logger.debug("Lake segmentation method = SLIC")
        height_img = get_2d_from_1d(in_height, in_x, in_y)
        labeled_img = segmentation.slic(height_img, n_segments=20)
        out_labels = get_1d_from_2d(labeled_img, in_x, in_y)

    elif height_segmentation_method == 3 :
        logger.debug("Lake segmentation method = Unsupervised thresholding")
        height_img = get_2d_from_1d(in_height, in_x, in_y)
        # compute threshold
        text_threshold = filters.threshold_local(height_img,block_size=51, offset=10)
        # apply threshold
        labeled_img = height_img > text_threshold
        out_labels = get_1d_from_2d(labeled_img, in_x, in_y)

    elif height_segmentation_method == 4 :
        logger.debug("Lake segmentation method = Watershed segmentation method")
        height_img = get_2d_from_1d(in_height, in_x, in_y)
        nb_line, nb_col = height_img.shape
        k = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
        # smmoth images
        height_img_smoothed = ndimage.convolve(height_img, k, mode='constant', cval=0.0)
        distance = ndimage.distance_transform_edt(height_img_smoothed)
        min_size = np.min([nb_line, nb_col])
        # compute local maximums
        if min_size < 5 :
            local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((min_size, min_size)))
        else :
            local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((5,5)))
        # if image contains more than 1 maximum local
        if np.sum(local_maxi) >= 2:
            local_max_dist = distance[np.where(local_maxi==True)]
            flat = local_max_dist.flatten()
            flat.sort()
            second_max_val = flat[-2]
            local_maxi = distance * local_maxi  >= second_max_val
            # Watershed segmentation
            markers = ndimage.label(local_maxi)[0]
            distance = ndimage.distance_transform_edt(height_img)
            labeled_img = watershed(-distance, markers, mask = height_img)
            out_labels = get_1d_from_2d(labeled_img, in_x, in_y)
        else :
            out_labels = np.zeros((nb_pts))

    elif height_segmentation_method == 5:
        logger.debug("Lake segmentation method = k-means clustering method")
        std_heigth = np.std(in_height)
        nb_classes = 1

        min_count = 100
        nb_pix_min_by_class = 30

        labels = None

        in_std_height_max = 1
        out_labels = np.ones(in_height.shape)
        while std_heigth > in_std_height_max and min_count > nb_pix_min_by_class:
            nb_classes += 1
            # Cluster height over nb_classes classes
            kmeans_classif = KMeans(n_clusters=nb_classes)
            kmeans_classif.fit(in_height.reshape(-1, 1))  # Turn line into column for in_height

            out_labels = kmeans_classif.labels_
            # Compute heigth std inside each class
            std_heigth = np.max([np.std(in_height[np.where(labels == curLabel)]) for curLabel in np.unique(labels)])

            unique_labels, counts = np.unique(labels, return_counts=True)
            min_count = np.min(counts)
            # If number of classes upper than 10 => stop iteration
            if nb_classes > 10:
                break
        logger.debug("NB classes : %d, max std : %f " % (nb_classes, std_heigth))

    elif height_segmentation_method == 6:
        logger.debug("Lake segmentation method = hierarchical clustering")
        # buid matrix
        thresh = 1.5
        mat = np.zeros((in_height.size, 2))
        mat[:, 0] = in_height
        mat[:, 1] = in_height
        tmp_labels = hierarchy.fclusterdata(mat, thresh, criterion="distance")
        out_labels = split_labels_by_region(in_x, in_y, tmp_labels)

    elif height_segmentation_method == 7:
        logger.debug("Lake segmentation method = Otsu thresholding")
        # Compute first threshold
        out_labels, split_lake_in_two = get_otsu_labelling(in_height, in_pix_area, min_size)
        # If lake slitted, try to split again separate parts
        if split_lake_in_two :
            # try to split first and second sub-lakes
            for i in [1, 2]:
                # Compute threshold for sub lakes
                out_labels_tmp, split_lake_in_more_than_two = get_otsu_labelling(in_height[np.where(out_labels == i)], in_pix_area, min_size)
                # if sub-lakes need to be split, build new label to avoid that 1 labels corresponds to 2 separate lakes
                if split_lake_in_more_than_two :
                    out_labels[np.where(out_labels == i)] = out_labels_tmp + 2
                    out_labels_tmp2 = out_labels
                    for i, label in enumerate(np.unique(out_labels)):
                        out_labels_tmp2[np.where(out_labels == label)] = i + 1
                    out_labels = out_labels_tmp2

    elif height_segmentation_method == 8:
        logger.debug("Lake segmentation method = MeanShift")

        height_img = get_2d_from_1d(in_height, in_x, in_y)

        nb_line, nb_col = height_img.shape
        logger.debug("nb_lin, nb_col : %d %d" %(nb_line, nb_col))

        # The following bandwidth can be automatically detected using
        data = np.vstack([in_height, in_x, in_y]).transpose()
        bandwidth = estimate_bandwidth(data, quantile=0.1, n_samples=500)

        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(data)
        out_labels = split_labels_by_region(in_x, in_y, ms.labels_)

    else :
        raise service_error.ProcessingError("Lake segmentation method %d unknown" % height_segmentation_method, logger)

    if(compute_labels):
        # filter labels to keep only labels with area > min_area
        out_labels = filter_image_by_region_size(in_x, in_y, out_labels, in_pix_area, min_size)
        # reorder labels to have labels from 1 to N in decreasing number of pixels.
        out_labels = reorder_labels(out_labels)

        labels, counts = np.unique(out_labels, return_counts=True)
        if len(labels) > 1 :
            logger.debug("Current lake contains %d sublakes" %(len(labels)))
            for i, label in enumerate(labels):
                std = np.std(in_height[np.where(out_labels == label)])
                mean = np.mean(in_height[np.where(out_labels == label)])
                logger.debug("Sublakes %d contains %d pixels with mean_height %f and std_height %f" % (label, counts[i], mean, std))

        # # =================================================================================================================
        # # TO REMOVE WHEN TESTS ARE OK
        # if np.unique(out_labels).size > 1:
        #
        #     fig_path = "/work/ALT/swot/swothr/users/cazalsc/Issue_70/Illustrations/figure_%d_pix.png" % (nb_pts)
        #     i=1
        #
        #     while os.path.exists(fig_path):
        #         fig_path.replace('.png', '_%d.png' %i)
        #         i+=1
        #     plot_height_and_segementation(in_height, out_labels, in_x, in_y, fig_path, cMap, norm)
        # # =================================================================================================================

    return out_labels


def get_1d_from_2d(in_img, in_x, in_y):
    """
    This function returns 1D array with values if pixels whith coordinates in_x and in_y in 2D matrix in_img.
    :param in_img: image
    :type in_img: 2D array of int or float
    :param in_x: pixels x coordinates
    :type in_x: 1D array of int
    :param in_y: pixels y coordinates
    :type in_y: 1D array of int

    :return: values of pixels in matrix
    :type: 1D array of int or float
    """

    out_vect = np.zeros(in_x.shape, dtype=in_img.dtype)
    nb_pts = in_x.size

    for ind in range(nb_pts):
        out_vect[ind] = in_img[in_y[ind], in_x[ind]]
    return out_vect


def get_2d_from_1d(in_values, in_x, in_y):
    """
    This function returns 2D array with values of pixels whith coordinates in_x and in_y filled from in_values.
    :param in_values: values of pixels
    :type in_values: 1D array of int or float
    :param in_x: pixels x coordinates
    :type in_x: 1D array of int
    :param in_y: pixels y coordinates
    :type in_y: 1D array of int

    :return: Matrix with in_values set
    :type: 2D array of int or float
    """

    if (in_x.size != in_y.size) or (in_x.size != in_values.size):
        raise ValueError("in_values, in_x and in_y must be the same size but are : in_values = %d in_x = %d and in_y = %d" \
                         % (in_values.size, in_x.size, in_y.size))
    nb_pts = in_x.size

    out_img = np.zeros((np.max(in_y) + 1, np.max(in_x) + 1), dtype=in_values.dtype)
    for ind in range(nb_pts):
        out_img[in_y[ind], in_x[ind]] = in_values[ind]

    return out_img


def split_labels_by_region(in_x, in_y, in_labels):
    """
    This function checks that each unique label of in_labels corresponds to one single separate entity.
    It it is not the case, a new label is affected to separate entities
    :param in_x: X indices of "1" pixels
    :type in_x: 1D vector of int
    :param in_y: Y indices of "1" pixels
    :type in_y: 1D vector of int
    :param in_labels: labels of pixels
    :type in_labels: 1D vector of int

    :return out_labels: new labels
    :type out_labels: 1D vector of int
    """
    out_labels = np.zeros((in_labels.shape))

    for i, label in enumerate(np.unique(in_labels)):
        idx = np.where(in_labels == label)
        mat = get_2d_from_1d(in_labels[idx], in_x[idx], in_y[idx])
        new_mat = my_tools.label_region(mat)[0]
        new_labels = get_1d_from_2d(new_mat, in_x[idx], in_y[idx])
        out_labels[idx] = np.max(out_labels) + new_labels

    return out_labels


def get_otsu_labelling(in_height, in_pix_area, min_size):
    """
    This function returns labels and split_lake_flag. A Otsu threshold is applied to in_height, if threshold is different
    form mean +/- 2std of each resulting classes and that each classe is bigger than min_size, the in_height is slitted
    in 2 classes. Labels are carried in out_labels, split_lake_flag indicates if in_height have been slitted.

    :param in_height: height of pixels
    :type in_height: 1D vector of float
    :param in_pix_area: area of pixels
    :type in_pix_area: 1D vector of float
    :param min_size: minimum objet size
    :type min_size: float

    :return out_labels: labels of pixels
    :type in_x: 1D vector of int
    :return split_lake_flag: indicates if in_height have been splitted
    :type split_lake_flag: boolean
    """
    # Check that more than one value can be found in in_height
    if len(set(in_height)) < 2:
        split_lake_flag = False
        labels = np.ones((in_height.size))
    else :
        # compute Otsu threshold
        threshold = filters.threshold_otsu(in_height)
        # compute corresponding labels
        labels = (in_height > threshold).astype('int') + 1
        labels_unique, counts = np.unique(labels, return_counts=True)
        # Check if a class contains is bigger enough to become a new lake
        if (counts * np.mean(in_pix_area) < min_size).any() or in_height.size <= 5 :
            split_lake_flag = False
        else :
            split_lake_flag = True
            #  Check if threshold is different from mean+/-2std for each class to split lake
            for i, label in enumerate(labels_unique):
                height_sub = in_height[np.where(labels == label)]
                std = np.std(height_sub)
                mean = np.mean(height_sub)
                if np.logical_and(threshold > mean - 2 * std, threshold < mean + 2 * std):
                    split_lake_flag = False

    if split_lake_flag == True:
        out_labels = labels
    else :
        out_labels = np.ones((in_height.size))

    return out_labels, split_lake_flag


def sliding_window(arr, window_size):
    """
    Construct a sliding window view of the array

    Code from Stack Over Flow : https://stackoverflow.com/a/11000193
    """
    arr = np.asarray(arr)
    window_size = int(window_size)
    if arr.ndim != 2:
        raise ValueError("need 2-D input")
    if not (window_size > 0):
        raise ValueError("need a positive window size")
    shape = (arr.shape[0] - window_size + 1,
             arr.shape[1] - window_size + 1,
             window_size, window_size)
    if shape[0] <= 0:
        shape = (1, shape[1], arr.shape[0], shape[3])
    if shape[1] <= 0:
        shape = (shape[0], 1, shape[2], arr.shape[1])
    strides = (arr.shape[1]*arr.itemsize, arr.itemsize,
               arr.shape[1]*arr.itemsize, arr.itemsize)
    return as_strided(arr, shape=shape, strides=strides)


def cell_neighbors(arr, i, j, d=1):
    """
    Return d-th neighbors of cell (i, j)

    Code from Stack Over Flow : https://stackoverflow.com/a/11000193
    """
    w = sliding_window(arr, 2*d+1)

    ix = np.clip(i - d, 0, w.shape[0]-1)
    jx = np.clip(j - d, 0, w.shape[1]-1)

    i0 = max(0, i - d - ix)
    j0 = max(0, j - d - jx)
    i1 = w.shape[2] - max(0, d - i + ix)
    j1 = w.shape[3] - max(0, d - j + jx)

    return w[ix, jx][i0:i1,j0:j1].ravel()


def get_sub_img(in_labels_img, in_x, in_y):
    """
    This fuction returns a subset of in_labels_img, including pixels with coordinates in_x and in_y and their neighbors.
    :param in_labels_img: image
    :type in_labels_img: 2D matrix of float or int
    :param in_x: pixels x coordinates
    :type in_x: 1D array of int
    :param in_y: pixels y coordinates
    :type in_y: 1D array of int

    :return: subset image
    :type: 2D matrix of float or int
    """
    size_x, size_y = in_labels_img.shape

    min_x = np.min(in_x) - 1
    if min_x <0:
        min_x = 0

    max_x = np.max(in_x) + 2
    if max_x>size_x:
        max_x = size_x

    min_y = np.min(in_y) - 1
    if min_y <0:
        min_y = 0

    max_y = np.max(in_y) + 2
    if max_y>size_y:
        max_y = size_y

    return in_labels_img[min_x:max_x, min_y:max_y]


def replace_label_with_maj_neighbor_value(in_labels_img, in_label):
    """
    This function returns in_image with label in_label replaced by the majority value in the neighborhood of pixels with label in_label.

    :param in_labels_img: image of labels
    :type in_labels_img: 2d matrix of int
    :param in_label: label to replace
    :type in_label: int

    :return: in_labels_img modified
    :type: 2d matrix of int
    """
    # Get instance of service config file
    logger = logging.getLogger("my_segmentation")


    idx_x, idx_y = np.where(in_labels_img == in_label)

    sub_img = get_sub_img(in_labels_img, idx_x, idx_y)

    sub_img_dilat = skimage.morphology.binary_dilation((sub_img == in_label))
    cell_neighbors = np.logical_xor(sub_img == in_label, sub_img_dilat)
    neighbors_values = sub_img[cell_neighbors == True]
    neighbors_values = neighbors_values[neighbors_values!=0]
    if neighbors_values.size > 1:
        maj_val = np.bincount(neighbors_values).argmax()
    else :
        maj_val = neighbors_values

    if maj_val.size == 0:
        logger.warning("Error in my_segementation : lake is not splitted")
        maj_val = 1
    in_labels_img[(idx_x, idx_y)] = maj_val

    return in_labels_img


def filter_image_by_region_size(in_x, in_y, in_labels, in_pix_area, min_size):
    """
    Replace labels of in_labels with area < min_size by majority value in the neighborhood.
    This function deals with small entities that cannot ben considered as a lake

    :param in_x: pixels x coordinates
    :type in_x: 1D array of int
    :param in_y: pixels y coordinates
    :type in_y: 1D array of int
    :param in_labels: label of pixels
    :type in_y: 1D array of int
    :param in_pix_area: pixels area
    :type in_pix_area: 1D array of float
    :param min_size: minimum size of object
    :type in_y: float

    :return: new labels
    :type: 1D array of int
    """
    # get 2S labels and pixel area
    labels_img = get_2d_from_1d(in_labels, in_x, in_y).astype('int')
    pix_area_img = get_2d_from_1d(in_pix_area, in_x, in_y)

    # build a new label image
    relabels_img = np.zeros(labels_img.shape, dtype='int')
    max_label = 1

    # for each label, compute binary matrix to get separate entities with label l, and associate a new label l to each separate entity.
    # ex: [[1, 0, 0], [0, 1, 1]]  -> [[1, 0, 0], [0, 2, 2]]
    for l in np.unique(in_labels):
        new_labels = my_tools.label(labels_img==l)[0]
        relabels_img[new_labels>0] += new_labels[new_labels>0]+max_label
        max_label = np.max(relabels_img)

    new_labels, counts = np.unique(relabels_img, return_counts=True)

    # List labels that with too small area and set major_value_exists to true if a label is bigger enough to be considered as a lake.
    label_too_small = []
    major_value_exists = False
    for nl in new_labels :
        if nl == 0:
            continue
        idx = np.where(relabels_img == nl)
        area = np.sum(pix_area_img[idx])

        if area < min_size or pix_area_img[idx].size <= 5 :
            label_too_small.append(nl)
        else :
            major_value_exists = True

    #  if a label is bigger enough, replace labels too small by neighborhood major value, otherwise set all labels to 1.
    if major_value_exists :
        for l in label_too_small:
            relabels_img = replace_label_with_maj_neighbor_value(relabels_img, l)
    else :
        relabels_img = np.ones(relabels_img.shape)

    out_labels = get_1d_from_2d(relabels_img, in_x, in_y)

    return out_labels


def reorder_labels(in_labels):
    """
    Reorders labels of in_labels from O to N, with a decreasing number of pixels.
    :param in_labels: label of pixels
    :type in_y: 1D array of int

    :return: new labels
    :type: 1D array of int
    """
    out_labels = np.zeros(in_labels.shape, dtype='int')
    unique, counts = np.unique(in_labels, return_counts=True)

    for i, label in enumerate(unique[np.argsort(-counts)]):
        out_labels[in_labels==label] = i+1
        
    return out_labels






# #======================================================================================================================
# # TO REMOVE WHEN TESTS ARE OK
# import os
# import matplotlib.pyplot as plt
#
# from matplotlib.colors import ListedColormap
# from matplotlib.colors import BoundaryNorm
#
# color_array = np.load("/work/ALT/swot/swothr/users/cazalsc/Issue_70/color_plot.npz")["color_array"]
# color_list = color_array
# cMap = ListedColormap(color_list)
#
# bounds = np.arange(-0.5, len(color_list) + 0.5, 1)
# norm = BoundaryNorm(bounds, cMap.N)
#
# def plot_height(in_height, in_x, in_y):
#
#     height_img = get_2d_from_1d(in_height, in_x, in_y)
#     height_img[np.where(height_img == 0)] = np.nan
#
#     plt.imshow(height_img)
#     plt.colorbar()
#
#     plt.show()
#
#     # Clear the current axes.
#     plt.cla()
#     # Clear the current figure.
#     plt.clf()
#     # Closes all the figure windows.
#     plt.close('all')
#
# def plot_height_and_segementation(in_height, in_labels, in_x, in_y, fig_path, cMap, norm):
#
#     height_img = get_2d_from_1d(in_height, in_x, in_y)
#     height_img[np.where(height_img == 0)] = np.nan
#     # np.ma.masked_where(y > 0.7, y)
#     fig, ax = plt.subplots(1, 2, figsize=(10,5))
#     tmp = ax[0].imshow(height_img)
#
#     out_labels_img = get_2d_from_1d(in_labels, in_x, in_y)
#     ax[1].imshow(out_labels_img, cmap=cMap, interpolation='nearest', norm=norm)
#
#
#     fig.colorbar(tmp, ax=ax[0])
#     if fig_path:
#         plt.savefig(fig_path)
#     else :
#         plt.show()
#
#     # Clear the current axes.
#     plt.cla()
#     # Clear the current figure.
#     plt.clf()
#     # Closes all the figure windows.
#     plt.close('all')
# #======================================================================================================================


