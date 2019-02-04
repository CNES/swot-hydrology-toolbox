# -*- coding: utf8 -*-
"""
.. module proc_pixc.py
    :synopsis: Deals with SWOT pixel cloud product
    02/28/2017 - Creation

.. module author: Claire POTTIER - CNES DSO/SI/TR
    
.. todo: sortir les pixels oceans de la liste des pixels à traiter
.. todo: contenu exact du PIXC à màj qd décidé
.. todo: entités qui contiennent pixels rivière et autres

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from __future__ import absolute_import, division, print_function, unicode_literals

import ast
import datetime
import numpy as np
import os
from osgeo import ogr, osr
from scipy import interpolate
import logging

import cnes.common.lib.my_basins as my_basins
import cnes.common.lib.my_netcdf_file as my_nc
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib_lake.locnes_variables as my_var


class PixelCloud(object):

    def __init__(self, IN_pixc_file, IN_idx_reject, IN_use_fractional_inundation=None):
        """
        Constructor: retrieve needed data from pixel cloud
        
        :param IN_pixc_file: full path of L2_HR_PIXC file
        :type IN_pixc_file: string
        :param IN_idx_reject: list of indices to reject before all processing
        :type IN_idx_reject: 1D-array of int
        :param IN_use_fractional_inundation:
                For which classes should the inundation fraction be used?
                The default is to assume that interior pixels are 100% water,
                But to use both land and water edge pixels partially by using the fractional inundation kwd.
        :type IN_use_fractional_inundation: bool list, default None

        Variables of the object:
        - From group pixc in L2_HR_PIXC file
            nb_pix_range / int: number of pixels in range dimension (= global attribute named nr_pixels in L2_HR_PIXC)
            nb_pix_azimuth / int: number of pixels in azimuth dimension (= global attribute named nr_lines in L2_HR_PIXC)
            cycle_num / int: cycle number (= global attribute named cycle_number in L2_HR_PIXC)
            pass_num / int: pass number (= global attribute named pass_number in L2_HR_PIXC)
            tile_ref / int: tile reference (= global attribute named tile_ref in L2_HR_PIXC)
            origin_classif / 1D-array of int: classification value of pixels (= variable named classification in L2_HR_PIXC)
            range_idx / 1D-array of int: range indices of water pixels (= variable named range_index in L2_HR_PIXC)
            azimuth_idx / 1D-array of int: azimuth indices of water pixels (= variable named azimuth_index in L2_HR_PIXC)
            pixel_area / 1D-array of int: area of water pixels (= variable named pixel_area in L2_HR_PIXC)
            longitude / 1D-array of float: longitude of water pixels (= variable named longitude in L2_HR_PIXC)
            latitude / 1D-array of float: latitude of water pixels (= variable named latitude in L2_HR_PIXC)
            height / 1D-array of float: height of water pixels (= variable named height in L2_HR_PIXC)
            crosstrack / 1D-array of float: cross-track distance from nadir to center of water pixels (= variable named crosstrack in L2_HR_PIXC)
        - From group tvp in L2_HR_PIXC file
            nadir_time / 1D-array of float: observation time of each nadir pixel (= variable named time in L2_HR_PIXC file)
            nadir_longitude / 1D-array of float: longitude of each nadir pixel (= variable named longitude in L2_HR_PIXC file)
            nadir_latitude / 1D-array of float: latitude of each nadir pixel (= variable named latitude in L2_HR_PIXC file)   
            nadir_[x|y|z] / 1D-array of float: [x|y|z] cartesian coordinates of each nadir pixel (= variables named [x|y|z] in L2_HR_PIXC file)
            nadir_[vx|vy|vz] / 1D-array of float: velocity vector of each nadir pixel in cartesian coordinates (= variables named velocity_unit_[x|y|z] in L2_HR_PIXC file)
            near_range / 1D-array of float: near range distance for each nadir point (= variable named near_range in L2_HR_PIXC file)
        - From processing
            nb_water_pix / int: number of water pixels
            pixc_vec_river / PixelCloudVecRiver object: associed pixel cloud river complementary file
            selected_idx / 1D-array of int: indices from original 1D-arrays of not rejected pixels with classification indices listed in IN_classif
            nb_selected / int: number of selected pixels (=selected_idx.size)
            tile_poly / ogr.Polygon: polygon of the PixC tile
            continent / string: continent covered by the tile (if global var CONTINENT_FILE exists)
            labels / 1D-array of int: labelled regions associated to each PixC water pixel; pixels of this vector correspond one-to-one to the pixels of data from L2_HR_PIXC and L2_HR_PIXC_VEC_RIVER
            nb_obj / int : number of separate entities in the PixC tile
            labels_inside / 1D-array of int: label of objects entirely inside the tile
            nb_obj_inside / int : number of these objects
            labels_edge / 1D-array of int: label of objects at the top/bottom edges of the tile
            nb_obj_edge / int : number of these objects 
            edge_idx / 1D-array of int: indices of pixels contained in objects at top/bottom edges
            edge_label / 1D-array of int: object label for each pixel contained in objects at top/bottom edges
            edge_loc / 1D-array of int: object edge location (0=bottom 1=top 2=both) for each pixel contained in objects at top/bottom edges
            nb_edge_pix / int: number of pixels contained in objects at top/bottom edges
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("L2_HR_PIXC file = %s" % IN_pixc_file)
        TMP_print = "Keeping pixels with classification flags ="
        if my_var.FLAG_WATER != "":
            TMP_print += " %s (WATER)" % my_var.FLAG_WATER
        if my_var.FLAG_DARK != "":
            TMP_print += " %s (DARK)" % my_var.FLAG_DARK
        if my_var.FLAG_LAYOVER != "":
            TMP_print += " %s (LAYOVER)" % my_var.FLAG_LAYOVER
        logger.info(TMP_print)

        # 1 - Retrieve needed information from pixel cloud file
        pixc_reader = my_nc.myNcReader(IN_pixc_file)
        pixc_group = pixc_reader.content.groups['pixel_cloud']
        sensor_group = pixc_reader.content.groups['tvp']
        
                # 1.1 - Number of pixels in range dimension
        self.nb_pix_range = int(pixc_reader.getAttValue("interferogram_size_range", IN_group=pixc_group))
        # 1.2 - Number of pixels in azimuth dimension
        self.nb_pix_azimuth = int(pixc_reader.getAttValue("interferogram_size_azimuth", IN_group=pixc_group))
        # 1.3 - Cycle number
        self.cycle_num = int(pixc_reader.getAttValue("cycle_number"))
        # 1.4 - Pass number
        self.pass_num = int(pixc_reader.getAttValue("pass_number"))
        # 1.5 - Tile reference
        self.tile_ref = pixc_reader.getAttValue("tile_name")
        # 1.6 - Classification flag
        self.origin_classif = pixc_reader.getVarValue("classification", IN_group=pixc_group)
        # 1.7 - Range indices of water pixels
        self.origin_range_idx = pixc_reader.getVarValue("range_index", IN_group=pixc_group)
        # 1.8 - Azimuth indices of water pixels
        self.origin_azimuth_idx = pixc_reader.getVarValue("azimuth_index", IN_group=pixc_group)
        # 1.9 - Sensor azimuth position of water pixels
        origin_illumination_time = pixc_reader.getVarValue("illumination_time", IN_group=pixc_group)
                        
        # 1.10 - Continent overflown
        # Create polygon of tile from global attributes
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(float(pixc_reader.getAttValue("inner_first_longitude")), float(pixc_reader.getAttValue("inner_first_latitude")))
        ring.AddPoint(float(pixc_reader.getAttValue("outer_first_longitude")), float(pixc_reader.getAttValue("outer_first_latitude")))
        ring.AddPoint(float(pixc_reader.getAttValue("outer_last_longitude")), float(pixc_reader.getAttValue("outer_last_latitude")))
        ring.AddPoint(float(pixc_reader.getAttValue("inner_last_longitude")), float(pixc_reader.getAttValue("inner_last_latitude")))
        ring.AddPoint(float(pixc_reader.getAttValue("inner_first_longitude")), float(pixc_reader.getAttValue("inner_first_latitude")))
        self.tile_poly = ogr.Geometry(ogr.wkbPolygon)
        self.tile_poly.AddGeometry(ring)
        self.continent = my_basins.link_poly_to_continent(self.tile_poly)
                
        # 2 - Get indices corresponding to input classification flags
        # 2.1 - Flag rejected pixels in a specific way
        TMP_classif = self.origin_classif  # Init temporary classif vector
        if IN_idx_reject is None:
            logger.info("No pixel to be rejected")
        else:
            logger.info("%d pixels to be rejected" % IN_idx_reject.size)
            TMP_classif[IN_idx_reject] = 100  # Dummy value to ease selection hereafter
        # 2.2 - Build list of classification flags to keep
        TMP_flags = ""
        if my_var.FLAG_WATER != "":
            TMP_flags += my_var.FLAG_WATER.replace('"','')  # Water flags
        if my_var.FLAG_DARK != "":
            if TMP_flags != "":
                TMP_flags += ";"
            TMP_flags += my_var.FLAG_DARK.replace('"','')  # Dark water flags
        if my_var.FLAG_LAYOVER != "":
            if TMP_flags != "":
                TMP_flags += ";"
            TMP_flags += my_var.FLAG_LAYOVER.replace('"','')  # Layover flags
        # 2.3 - Get list of classification flags
        list_classif_flags = TMP_flags.split(";")
        # 2.4 - Get list of selected indices
        self.selected_idx = None  # Init wanted indices vector
        for classif_flag in list_classif_flags:
            vInd = np.where(TMP_classif == int(classif_flag))[0]
            logger.info("%d pixels with classification flag = %d" % (vInd.size, int(classif_flag)))
            if vInd.size != 0:
                if self.selected_idx is None:
                    self.selected_idx = vInd
                else:
                    self.selected_idx = np.concatenate((self.selected_idx, vInd))
        if self.selected_idx is None:
            self.nb_selected = 0
        else:
            self.nb_selected = self.selected_idx.size
        logger.info("=> %d pixels to keep" % self.nb_selected)

        # 3 - Keep PixC data only for selected pixels
        if self.nb_selected != 0:
            
            # 3.1 - In PixC group
            
            # Classification flags
            self.classif = self.origin_classif[self.selected_idx]
            # Range indices of water pixels
            self.range_idx = self.origin_range_idx[self.selected_idx]
            # Number of water pixels
            self.nb_water_pix = self.range_idx.size
            # Range indices of water pixels
            self.azimuth_idx = self.origin_azimuth_idx[self.selected_idx]
            # Sensor azimuth position of water pixels
            illumination_time = origin_illumination_time[self.selected_idx]
            # Pixel area
            self.pixel_area = pixc_reader.getVarValue("pixel_area", IN_group=pixc_group)[self.selected_idx]
            # Longitude
            self.longitude = pixc_reader.getVarValue("longitude", IN_group=pixc_group)[self.selected_idx]
            # Latitude
            self.latitude = pixc_reader.getVarValue("latitude", IN_group=pixc_group)[self.selected_idx]
            # Height
            self.height = pixc_reader.getVarValue("height", IN_group=pixc_group)[self.selected_idx]
            # Cross-track distance
            self.crosstrack = pixc_reader.getVarValue("cross_track", IN_group=pixc_group)[self.selected_idx]

            # Inundation fraction
            # The inundation fraction is, in theory, a number between 0 and 1 estimating the fraction of the pixel covered with water.
            # In practice, because of noise, it can be outside the bounds and even be negative!
            # It will produced an ensemble unbiased estimate if the class mean cross sections are known.
            if IN_use_fractional_inundation is None:

                logger.info("Fractional inundation not used")
                # All water pixels are supposed inundated
                self.fractional_inundation = None
                self.inundated_area = self.pixel_area

            else:
                logger.info("Fractional inundation used!")

                # Get continuous classification (water fraction)
                self.fractional_inundation = pixc_reader.getVarValue("continuous_classification", IN_group=pixc_group)[self.selected_idx]

                # Get classification at the selected indexes
                classif_idx = self.origin_classif[self.selected_idx]

                # Get flags to know for which classes should the inundation fraction be used
                use_fractional_inundation = IN_use_fractional_inundation.split(";")

                # Initialize inundated area
                self.inundated_area = self.pixel_area

                # Update it when use_fractional_inundation for class i is True
                for i, k in enumerate(list_classif_flags):
                    if ast.literal_eval(use_fractional_inundation[i]):
                        logger.info("[use_fractional_inundation[{}]: {} for class: {}".format(i, use_fractional_inundation[i], k))
                        index = np.where(classif_idx == int(k))
                        self.inundated_area[index] = self.pixel_area[index] * self.fractional_inundation[index]

            # 3.2 - In TVP group
            
            # Interpolate nadir_time wrt illumination time
            TMP_nadir_time = pixc_reader.getVarValue("time", IN_group=sensor_group)  # Read nadir_time values
            f = interpolate.interp1d(TMP_nadir_time, range(len(TMP_nadir_time)))  # Interpolator
            nadir_idx = (np.rint(f(illumination_time))).astype(int)  # Link between PixC and nadir pixels
            
            # Nadir time
            self.nadir_time = TMP_nadir_time[nadir_idx]
            # Nadir longitude
            self.nadir_longitude = pixc_reader.getVarValue("longitude", IN_group=sensor_group)[nadir_idx]
            # Nadir latitude
            self.nadir_latitude = pixc_reader.getVarValue("latitude", IN_group=sensor_group)[nadir_idx]
            # Nadir cartesian coordinates
            self.nadir_x = pixc_reader.getVarValue("x", IN_group=sensor_group)[nadir_idx]
            self.nadir_y = pixc_reader.getVarValue("y", IN_group=sensor_group)[nadir_idx]
            self.nadir_z = pixc_reader.getVarValue("z", IN_group=sensor_group)[nadir_idx]
            # Nadir velocity in cartesian coordinates
            self.nadir_vx = pixc_reader.getVarValue("vx", IN_group=sensor_group)[nadir_idx]
            self.nadir_vy = pixc_reader.getVarValue("vy", IN_group=sensor_group)[nadir_idx]
            self.nadir_vz = pixc_reader.getVarValue("vz", IN_group=sensor_group)[nadir_idx]
            # Near range distance at each nadir point 
            self.near_range = pixc_reader.getAttValue("near_range")
                    
        # 4.8 - Close file
        pixc_reader.close()

        # 5 - Other initialization
        self.labels = None  # Vector of entity labels associated to each pixel
        self.nb_obj = None  # Number of separate entities
        self.labels_inside = None  # Labels of entities entirely inside the tile
        self.nb_obj_inside = None  # Number of entities inside the tile
        self.labels_at_top_edge = None  # Labels of entities at the top edge of the tile
        self.nb_obj_at_top_edge = None  # Number of entities at the top edge of the tile
        self.labels_at_bottom_edge = None  # Labels of entities at the bottom edge of the tile
        self.nb_obj_at_bottom_edge = None  # Number of entities at the bottom edge of the tile
        self.labels_at_both_edges = None  # Labels of entities at the top and bottom edges of the tile
        self.nb_obj_at_both_edges = None  # Number of entities at the top and bottom edges of the tile
        self.edge_idx = None  # Indices of pixels contained in objects at top/bottom edges
        self.edge_label = None  # Object label for each pixel contained in objects at top/bottom edges
        self.edge_loc = None  # Object edge location (0=bottom 1=top 2=both) for each pixel contained in objects at top/bottom edges
        self.nb_edge_pix = 0  # Number of pixels contained in objects at top/bottom edges

    # ----------------------------------------

    def computeWaterMask(self):
        """
        Create the water mask (i.e. a 2D binary matrix) in radar geometry,
        from the pixel cloud (1D-array layers of azimuth_index, range_index, classification and continuous classification)
        
        :return: water mask in radar geometry, i.e. a 2D matrix with "1" for each (IN_X_i, IN_Y_i)  having classification=0 and 0 elsewhere
        :rtype: 2D binary matrix of int 0/1
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")

        return my_tools.computeBinMat(self.nb_pix_range, self.nb_pix_azimuth, self.range_idx, self.azimuth_idx)

    def computeSeparateEntities(self):
        """
        Identify all separate entities in the water mask
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")

        # 1 - Create the water mask
        water_mask = self.computeWaterMask()

        # 2 - Identify all separate entities in a 2D binary mask
        sepEntities, self.nb_obj = my_tools.labelRegion(water_mask)

        # 3 - Convert 2D labelled mask in 1D-array layer of the same size of the L2_HR_PIXC layers
        self.labels = my_tools.convert2dMatIn1dVec(self.range_idx, self.azimuth_idx, sepEntities)
        self.labels = self.labels.astype(int)  # Conversion from float to integer

        # 4 - For each label : check if only one lake is in each label and relabels if necessary
        labels_tmp = np.zeros(self.labels.shape)

        for label in np.unique(self.labels):
            idx = np.where(self.labels == label)

            min_rg = min(self.range_idx[idx])
            min_az = min(self.azimuth_idx[idx])

            relabel_obj = my_tools.relabelLakeUsingSegmentationHeigth(self.range_idx[idx] - min_rg,
                                                                      self.azimuth_idx[idx] - min_az,
                                                                      self.height[idx])

            labels_tmp[self.labels == label] = np.max(labels_tmp) + relabel_obj

        self.labels = labels_tmp
        self.nb_obj = np.unique(self.labels).size

    def computeObjInsideTile(self):
        """
        Separate labels of lakes and new objects entirely inside the tile, from labels of objects at top or bottom of the tile
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")

        # 1 - Identify objects at azimuth = 0
        idx_at_az_0 = np.where(self.azimuth_idx == 0)[0]
        labels_at_az_0 = np.unique(self.labels[idx_at_az_0])

        # 2 - Identify objects at azimuth = max
        idx_at_az_max = np.where(self.azimuth_idx == self.nb_pix_azimuth - 1)[0]
        labels_at_az_max = np.unique(self.labels[idx_at_az_max])

        # 3 - Identify objects at azimuth = 0 and azimuth = max
        self.labels_at_both_edges = np.intersect1d(labels_at_az_0, labels_at_az_max)
        self.nb_obj_at_both_edges = self.labels_at_both_edges.size
        logger.info("> %d labels at bottom (az=0) AND top (az=%d) of the tile" % (self.nb_obj_at_both_edges, self.nb_pix_azimuth))
        for ind in np.arange(self.nb_obj_at_both_edges):
            logger.debug("%d" % self.labels_at_both_edges[ind])

        # 4 - Identify labels...
        # 4.1 - Only at azimuth = 0
        self.labels_at_bottom_edge = np.setdiff1d(labels_at_az_0, self.labels_at_both_edges)
        self.nb_obj_at_bottom_edge = self.labels_at_bottom_edge.size
        logger.info("> %d labels at bottom of the tile (az=0)" % self.nb_obj_at_bottom_edge)
        for ind in np.arange(self.nb_obj_at_bottom_edge):
            logger.debug("%d" % self.labels_at_bottom_edge[ind])
        # 4.2 - Only at azimuth = max
        self.labels_at_top_edge = np.setdiff1d(labels_at_az_max, self.labels_at_both_edges)
        self.nb_obj_at_top_edge = self.labels_at_top_edge.size
        logger.info("> %d labels at top of the tile (az=%d)" % (self.nb_obj_at_top_edge, self.nb_pix_azimuth))
        for ind in np.arange(self.nb_obj_at_top_edge):
            logger.debug("%d" % self.labels_at_top_edge[ind])

        # 5 - Get labels of objects entirely inside the tile
        self.labels_inside = np.arange(self.nb_obj) + 1  # Initialisation
        if self.nb_obj_at_bottom_edge != 0:
            self.labels_inside = np.setdiff1d(self.labels_inside, self.labels_at_bottom_edge)  # Delete labels at bottom
        if self.nb_obj_at_top_edge != 0:
            self.labels_inside = np.setdiff1d(self.labels_inside, self.labels_at_top_edge)  # Delete labels at top
        if self.nb_obj_at_both_edges != 0:
            self.labels_inside = np.setdiff1d(self.labels_inside,
                                              self.labels_at_both_edges)  # Delete labels at top and bottom
        self.nb_obj_inside = self.labels_inside.size
        logger.info("> %d objects entirely inside the tile" % self.nb_obj_inside)

    # ----------------------------------------

    def write_edge_file(self, IN_out_file, noval, compress=False):
        """
        Save pixels related to objects at top/bottom edge of the PixC tile in a NetCDF file.
        These 1D-arrays indicate pixel index, associated label, location (top/bottom/both edges) and needed PixC variables.
        If there is no pixel at tile edge, the file is empty but exists.
        
        :param IN_out_file: output file full path
        :type IN_out_file: string
        :param noval: No data value
        :type noval: float
        :param compress: parameter the define to compress or not the file
        :type compress: boolean
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Output L2_HR_LAKE_TILE edge file = %s" % IN_out_file)

        # 1 - Open file in writing mode
        data = my_nc.myNcWriter(IN_out_file)

        if noval is None:
            noval = -999000000.

        # 2 - Write file
        if self.nb_selected == 0:
            logger.info("NO selected PixC => empty edge file generated")
            data.add_dimension('record', 0)

        else:

            if (self.nb_obj_at_top_edge + self.nb_obj_at_bottom_edge + self.nb_obj_at_both_edges) == 0:
                logger.info(
                    "No object at top/bottom edge of the PixC tile => empty edge file generated")
                data.add_dimension('record', 0)

            else:

                # 2.1 - Get edge pixels indices and their associated label
                flag_start = False
                # 2.1.1 - Fill with bottom edge objects
                if self.nb_obj_at_bottom_edge != 0:
                    flag_start = True
                    self.edge_idx = np.where(self.labels == self.labels_at_bottom_edge[0])[
                        0]  # Get PixC pixels related to bottom edge object
                    self.edge_label = np.ones(self.edge_idx.size) * self.labels_at_bottom_edge[
                        0]  # Associated label vector
                    self.edge_loc = np.zeros(self.edge_idx.size)  # Associated location: 0=bottom 1=top 2=both
                    for indl in np.arange(1, self.nb_obj_at_bottom_edge):  # Loop over the other bottom edge objects
                        tmp_idx = np.where(self.labels == self.labels_at_bottom_edge[indl])[
                            0]  # Get pixels related to edge object
                        self.edge_idx = np.append(self.edge_idx, tmp_idx)
                        self.edge_label = np.append(self.edge_label,
                                                    np.ones(tmp_idx.size) * self.labels_at_bottom_edge[indl])
                        self.edge_loc = np.append(self.edge_loc, np.zeros(tmp_idx.size))
                # 2.1.2 - Fill with top edge objects
                if self.nb_obj_at_top_edge != 0:
                    if flag_start is False:
                        flag_start = True
                        idx_start = 1
                        self.edge_idx = np.where(self.labels == self.labels_at_top_edge[0])[
                            0]  # Get PixC pixels related to top edge object
                        self.edge_label = np.ones(self.edge_idx.size) * self.labels_at_top_edge[
                            0]  # Associated label vector
                        self.edge_loc = np.zeros(self.edge_idx.size) + 1  # Associated location: 0=bottom 1=top 2=both
                    else:
                        idx_start = 0
                    for indl in np.arange(idx_start, self.nb_obj_at_top_edge):  # Loop over the other top edge objects
                        tmp_idx = np.where(self.labels == self.labels_at_top_edge[indl])[
                            0]  # Get pixels related to edge object
                        self.edge_idx = np.append(self.edge_idx, tmp_idx)
                        self.edge_label = np.append(self.edge_label,
                                                    np.ones(tmp_idx.size) * self.labels_at_top_edge[indl])
                        self.edge_loc = np.append(self.edge_loc, np.zeros(tmp_idx.size) + 1)
                # 2.1.3 - Fill with bottom and top edges objects
                if self.nb_obj_at_both_edges != 0:
                    if flag_start is False:
                        idx_start = 1
                        self.edge_idx = np.where(self.labels == self.labels_at_both_edges[0])[
                            0]  # Get PixC pixels related to both edges object
                        self.edge_label = np.ones(self.edge_idx.size) * self.labels_at_both_edges[
                            0]  # Associated label vector
                        self.edge_loc = np.zeros(self.edge_idx.size) + 2  # Associated location: 0=bottom 1=top 2=both
                    else:
                        idx_start = 0
                    for indl in np.arange(idx_start, self.nb_obj_at_both_edges):  # Loop over the other both edges objects
                        tmp_idx = np.where(self.labels == self.labels_at_both_edges[indl])[0]  # Get pixels related to edge object
                        self.edge_idx = np.append(self.edge_idx, tmp_idx)
                        self.edge_label = np.append(self.edge_label,
                                                    np.ones(tmp_idx.size) * self.labels_at_both_edges[indl])
                        self.edge_loc = np.append(self.edge_loc, np.zeros(tmp_idx.size) + 2)
                # 2.1.4 - Number of edge pixels
                self.nb_edge_pix = self.edge_idx.size

                # 2.2 - Write output NetCdf file

                data.add_dimension('record', self.nb_edge_pix)

                data.add_variable('edge_idx', np.int32, 'record', IN_fill_value=np.int(noval), IN_compress=compress)
                data.fill_variable('edge_idx', self.selected_idx[self.edge_idx])
                data.add_variable('edge_label', np.int32, 'record', IN_fill_value=np.int(noval), IN_compress=compress)
                data.fill_variable('edge_label', self.edge_label)
                data.add_variable('edge_loc', np.int32, 'record', IN_fill_value=np.int(noval), IN_compress=compress)
                data.fill_variable('edge_loc', self.edge_loc)

                data.add_variable('azimuth_index', np.int32, 'record', IN_fill_value=np.int(noval), IN_compress=compress)
                data.fill_variable('azimuth_index', self.azimuth_idx[self.edge_idx])
                data.add_variable('range_index', np.int32, 'record', IN_fill_value=np.int(noval), IN_compress=compress)
                data.fill_variable('range_index', self.range_idx[self.edge_idx])

                data.add_variable('classif', np.int32, 'record', IN_fill_value=np.int(noval), IN_compress=compress)
                data.fill_variable('classif', self.origin_classif[self.edge_idx])
                data.add_variable('pixel_area', np.float64, 'record', IN_fill_value=np.float(noval), IN_compress=compress)
                data.fill_variable('pixel_area', self.pixel_area[self.edge_idx])

                data.add_variable('latitude', np.double, 'record', IN_fill_value=noval, IN_compress=compress)
                data.add_variable_attribute('latitude', 'units', 'degrees_north')
                data.fill_variable('latitude', self.latitude[self.edge_idx])
                data.add_variable('longitude', np.double, 'record', IN_fill_value=noval, IN_compress=compress)
                data.add_variable_attribute('longitude', 'units', 'degrees_east')
                data.fill_variable('longitude', self.longitude[self.edge_idx])
                data.add_variable('height', np.float64, 'record', IN_fill_value=np.float(noval), IN_compress=compress)
                data.add_variable_attribute('height', 'units', 'm')
                data.fill_variable('height', self.height[self.edge_idx])
                data.add_variable('crosstrack', np.float64, 'record', IN_fill_value=np.float(noval), IN_compress=compress)
                data.fill_variable('crosstrack', self.crosstrack[self.edge_idx])

                data.add_variable('nadir_time', np.float64, 'record', IN_fill_value=np.float(noval), IN_compress=compress)
                data.fill_variable('nadir_time', self.nadir_time[self.edge_idx])
                data.add_variable('nadir_lon', np.double, 'record', IN_fill_value=noval, IN_compress=compress)
                data.add_variable_attribute('nadir_lon', 'units', 'degrees_east')
                data.fill_variable('nadir_lon', self.nadir_longitude[self.edge_idx])
                data.add_variable('nadir_lat', np.double, 'record', IN_fill_value=noval, IN_compress=compress)
                data.add_variable_attribute('nadir_lat', 'units', 'degrees_north')
                data.fill_variable('nadir_lat', self.nadir_latitude[self.edge_idx])
                data.add_variable('nadir_x', np.float64, 'record', IN_fill_value=np.float(noval), IN_compress=compress)
                data.fill_variable('nadir_x', self.nadir_x[self.edge_idx])
                data.add_variable('nadir_y', np.float64, 'record', IN_fill_value=np.float(noval), IN_compress=compress)
                data.fill_variable('nadir_y', self.nadir_y[self.edge_idx])
                data.add_variable('nadir_z', np.float64, 'record', IN_fill_value=np.float(noval), IN_compress=compress)
                data.fill_variable('nadir_z', self.nadir_z[self.edge_idx])
                data.add_variable('nadir_vx', np.float64, 'record', IN_fill_value=np.float(noval), IN_compress=compress)
                data.fill_variable('nadir_vx', self.nadir_vx[self.edge_idx])
                data.add_variable('nadir_vy', np.float64, 'record', IN_fill_value=np.float(noval), IN_compress=compress)
                data.fill_variable('nadir_vy', self.nadir_vy[self.edge_idx])
                data.add_variable('nadir_vz', np.float64, 'record', IN_fill_value=np.float(noval), IN_compress=compress)
                data.fill_variable('nadir_vz', self.nadir_vz[self.edge_idx])
                
        data.add_global_attribute('nadir_near_range', self.near_range)
        data.add_global_attribute('producer', my_var.PRODUCER)
        data.add_global_attribute('creation_date', str(datetime.datetime.now()))
        data.add_global_attribute('cycle_number', self.cycle_num)
        data.add_global_attribute('pass_number', self.pass_num)
        data.add_global_attribute('tile_ref', self.tile_ref)
        if my_var.CONTINENT_FILE is not None:
            data.add_global_attribute('continent', self.continent)
        data.add_global_attribute('nr_pixels', self.nb_pix_range)
        data.add_global_attribute('nr_lines', self.nb_pix_azimuth)
        data.add_global_attribute('lake_db', my_var.LAKE_DB)
        if my_var.CONTINENT_FILE is not None:
            data.add_global_attribute('continent_file', my_var.CONTINENT_FILE)
        data.add_global_attribute('flag_water', my_var.FLAG_WATER)
        data.add_global_attribute('flag_dark', my_var.FLAG_DARK)
        data.add_global_attribute('flag_layover', my_var.FLAG_LAYOVER)
        data.add_global_attribute('min_size', my_var.MIN_SIZE)
        TMP_geoloc = 0
        if my_var.IMP_GEOLOC:
            TMP_geoloc = 1
        data.add_global_attribute('imp_geoloc', TMP_geoloc)
        data.add_global_attribute('hull_method', my_var.HULL_METHOD)
        data.add_global_attribute('std_height_max', my_var.STD_HEIGHT_MAX)
        data.add_global_attribute('biglake_model', my_var.BIGLAKE_MODEL)
        data.add_global_attribute('biglake_min_size', my_var.BIGLAKE_MIN_SIZE)
        data.add_global_attribute('biglake_grid_spacing', my_var.BIGLAKE_GRID_SPACING)
        data.add_global_attribute('biglake_grid_res', my_var.BIGLAKE_GRID_RES)

        # 3 - Close file
        data.close()

    def write_edge_file_asShp(self, IN_filename):
        """
        Write PixC subset related to edge (top/bottom) objects as a shapefile

        :param IN_filename: full path of the output file
        :type IN_filename: string
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Output L2_HR_LAKE_TILE edge shapefile = %s" % IN_filename)

        # 1 - Initialisation du fichier de sortie
        # 1.1 - Driver
        shpDriver = ogr.GetDriverByName(str("ESRI Shapefile"))
        # 1.2 - Creation du fichier
        if os.path.exists(IN_filename):
            shpDriver.DeleteDataSource(IN_filename)
        outDataSource = shpDriver.CreateDataSource(IN_filename)
        # 1.3 - Creation de la couche
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)  # WGS84
        outLayer = outDataSource.CreateLayer(str(str(os.path.basename(IN_filename)).replace('.shp', '')), srs, geom_type=ogr.wkbPoint)
        # 1.4 - Creation des attributs
        outLayer.CreateField(ogr.FieldDefn(str('EDGE_IDX'), ogr.OFTInteger))
        outLayer.CreateField(ogr.FieldDefn(str('EDGE_LAB'), ogr.OFTInteger))
        outLayer.CreateField(ogr.FieldDefn(str('EDGE_LOC'), ogr.OFTInteger))
        outLayer.CreateField(ogr.FieldDefn(str('AZ_INDEX'), ogr.OFTInteger))
        outLayer.CreateField(ogr.FieldDefn(str('R_INDEX'), ogr.OFTInteger))
        tmpField = ogr.FieldDefn(str('HEIGHT'), ogr.OFTReal)
        tmpField.SetWidth(10)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('NADIR_T'), ogr.OFTReal)
        tmpField.SetWidth(10)
        tmpField.SetPrecision(3)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('NAD_LON'), ogr.OFTReal)
        tmpField.SetWidth(10)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        tmpField = ogr.FieldDefn(str('NAD_LAT'), ogr.OFTReal)
        tmpField.SetWidth(10)
        tmpField.SetPrecision(6)
        outLayer.CreateField(tmpField)
        outLayerDefn = outLayer.GetLayerDefn()

        # 2 - On traite point par point
        for indp in range(self.nb_edge_pix):
            tmp_idx = self.edge_idx[indp]
            # 2.1 - On cree l'objet dans le format de la couche de sortie
            outFeature = ogr.Feature(outLayerDefn)
            # 2.2 - On lui assigne le point
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(self.longitude[tmp_idx], self.latitude[tmp_idx])
            outFeature.SetGeometry(point)
            # 2.3 - On lui assigne les attributs
            outFeature.SetField(str('EDGE_IDX'), int(self.selected_idx[tmp_idx]))
            outFeature.SetField(str('EDGE_LAB'), int(self.edge_label[indp]))
            outFeature.SetField(str('EDGE_LOC'), int(self.edge_loc[indp]))
            outFeature.SetField(str('AZ_INDEX'), int(self.azimuth_idx[tmp_idx]))
            outFeature.SetField(str('R_INDEX'), int(self.range_idx[tmp_idx]))
            outFeature.SetField(str('HEIGHT'), float(self.height[tmp_idx]))
            outFeature.SetField(str('NADIR_T'), float(self.nadir_time[tmp_idx]))
            outFeature.SetField(str('NAD_LON'), float(self.nadir_longitude[tmp_idx]))
            outFeature.SetField(str('NAD_LAT'), float(self.nadir_latitude[tmp_idx]))
            # 2.4 - On ajoute l'objet dans la couche de sortie
            outLayer.CreateFeature(outFeature)

        # 3 - Destroy the data sources to free resources
        outDataSource.Destroy()
