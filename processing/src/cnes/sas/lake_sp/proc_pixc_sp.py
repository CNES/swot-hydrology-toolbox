# -*- coding: utf8 -*-
"""
.. module:: proc_pixc_sp.py
    :synopsis: Deals with subset of pixel cloud, for objects located at top or bottom edge of a tile; i.e. gather pixels involved in edge lake product retrieved from all tiles of L2_HR_LakeTile_edge files
    Created on 27/09/2017

.. moduleauthor:: Cécile Cazals - CS

Copyright (c) 2017 CNES. All rights reserved.
"""
from __future__ import absolute_import, division, print_function, unicode_literals 

import numpy as np

import cnes.common.lib.my_api as my_api
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_netcdf_file as my_nc


class PixC_Edge_SP(object):

    def __init__(self, IN_ascending, IN_lake_tile_edge_path_list):
        """
        This class is designed to process all L2_HR_LakeTile edge files of one swath. LakeTile edge files contain pixel cloud information (from L2_HR_PIXC) for pixels involved in lakes located at the top/bottom edges of tiles.
        The LakeTile edge file contains only PixC variables needed information, but also additional information like:
            - edge_loc field, the edge location of the lake : top of the tile, bottom of the tile or both top and bottom (0=bottom 1=top 2=both)
            - edge_label contains the object labels retrieved from PGE_LakeTile labeling process.
            - edge_idx contains the L2_HR_PIXC initial pixels indices. This information is needed to update the improved geoloc and tag fields of LakeTile pixcvec files.
        This class processes all LakeTile edge files of a single swath to gather all entities at the edge of tile into lake entities.
        
        :param IN_ascending: orbit orientation (True = ascending, False = descending)
        :type IN_ascending: boolean
        :param IN_lake_tile_edge_path_list: list of LakeTile edge files full path
        :type IN_lake_tile_edge_path_list: list of string

        Variables of the object:
            - From L2_HR_LakeTile edge file:
            
                Global attributes :
                    nb_pix_range / int: number of pixels in range dimension (= global attribute named nb_pix_range in LakeTile_edge)
                    nb_pix_azimuth / int: number of pixels in azimuth dimension (= global attribute named nb_pix_azimuth in LakeTile_edge)
                    cycle_num / int: cycle number (= global attribute named cycle_number in LakeTile_edge)
                    pass_num / int: pass number (= global attribute named pass_number in LakeTile_edge)

                Variables :
                    range_idx / 1D-array of int: range indices of water pixels (= variable named range_index in  LakeTile_edge)
                    azimuth_idx / 1D-array of int: azimuth indices of water pixels (= variable named azimuth_index in LakeTile_edge)
                    pixel_area / 1D-array of int: area of water pixels (= variable named pixel_area in LakeTile_edge)
                    height / 1D-array of float: height of water pixels (= variable named height_medium in L2_HR_PIXC_main and LakeTile_edge)
                    crosstrack / 1D-array of float: cross-track distance from nadir to center of water pixels (= variable named crosstrack_medium in LakeTile_edge)
                    nadir_time / 1D-array of float: observation time of each nadir pixel (= variable named time in LakeTile_edge)
                    nadir_longitude / 1D-array of float: longitude of each nadir pixel (= variable named nadir_lon in LakeTile_edge)
                    nadir_lattitude / 1D-array of float: latitude of each nadir pixel (= variable named latitude in LakeTile_edge)
                    nadir_[x|y|z] / 1D-array of float: [x|y|z] cartesian coordinates of each nadir pixel (= variables named nadir_[x|y|z] in LakeTile_edge)
                    nadir_[vx|vy|vz] / 1D-array of float: velocity vector of each nadir pixel in cartesian coordinates (= variables named nadir_[vx|vy|vz] in LakeTile_edge)
                    nadir_alt / 1D-array of float: satellite altitude at each nadir point (= variable named nadir_alt in LakeTile_edge)
                    near_range / 1D-array of float: near range distance for each nadir point (= variable named nadir_near_range in LakeTile_edge)
                    latitude / 1D-array of float: latitude of water pixels (= variable named latitude_medium in L2_HR_PIXC_main and LakeTile_edge)
                    longitude / 1D-array of float: longitude of water pixels (= variable named longitude_medium in L2_HR_PIXC_main and LakeTile_edge)

        - From process:
            tile_ref / list of string : list of tile references to process ex: ['42N-R', '43N-R', '44N-R', '45N-R']
            nb_tiles / int : number of tiles to process
            nb_pixels / int : number of pixels to process
            tile_idx / tuple of ndarrays : indices of tile for each pixel
            ascending / bool : orbit orientation (True = ascending, False=descending)
            labels / 1D-array of int : arrays of new labels
        """
        my_api.printInfo("[PixelCloudSP] == INIT ==")

        # List of tile references to process ex: ['42N-R', '43N-R', '44N-R', '45N-R']
        self.tile_ref = []

        # Number of pixels to process
        self.nb_pixels = 0

        # Tile reference of each pixel
        self.tile_idx = []

        # Cycle number
        self.cycle_num = 0
        # Orbit number
        self.pass_num = 0
        # Continent 
        self.continent = None
        # Number of pixel in range
        self.nb_pix_range = 0

        # Orbit orientation (True = ascending, False=descending)
        self.ascending = IN_ascending

        # Merging of all PIXC_edge info into 1D vectors
        for lake_tile_edge_path in IN_lake_tile_edge_path_list:  # For each LakeTile edge file
            
            my_api.printInfo("[PixelCloudSP] Loading L2_HR_LakeTile edge file = %s" % lake_tile_edge_path)
            
            # Load data
            nb_pix_loaded, tile_ref = self.loadData(lake_tile_edge_path)

            my_api.printInfo("[PixelCloudSP] --> %d pixels loaded" % nb_pix_loaded)

            # Update nb_pixel
            self.nb_pixels += nb_pix_loaded

            # check if tile is neighboring the previus tile
            if len(self.tile_ref) >= 1:
                current_tile_number = int(tile_ref[:-3])
                previus_tile_number = int(self.tile_ref[-1][:-3])
                if current_tile_number != previus_tile_number + 1:
                    # if current tile is not adjacent to previous tile, add an empty tile to tile_ref
                    self.tile_ref.append("")

            # self.tile_idx must be set before incrementation of self.tile_ref to start by 0
            self.tile_idx += [len(self.tile_ref)] * nb_pix_loaded
            self.tile_ref.append(tile_ref)

        # Convert list to numpy array
        self.tile_idx = np.array(self.tile_idx)
        my_api.printInfo("[PixelCloudSP] %d PixC loaded for current swath" % self.nb_pixels)

        # Init labels to 0
        self.labels = np.zeros(self.nb_pixels)

    def loadData(self, IN_lake_tile_edge_filename):
        """
        This function loads NetCDF information.
        
        :param IN_lake_tile_edge_filename: full path of the NetCDF file to load
        :type IN_lake_tile_edge_filename: string
        
        :return OUT_nb_pix: number of pixel loaded
        :rtype OUT_nb_pix: integer
        :return OUT_tile_ref: reference of loaded tile
        :rtype OUT_tile_ref: string
        """

        # 1 - Open input NetCDF file in reading mode
        pixc_edge_reader = my_nc.myNcReader(IN_lake_tile_edge_filename)

        # 2 - Get tile reference
        OUT_tile_ref = pixc_edge_reader.getAttValue("tile_ref")

        # 3 - Initialization of object variables if not already done
        if not self.tile_ref:
            
            # 3.1 - Get global attributes
            self.cycle_num = pixc_edge_reader.getAttValue("cycle_number")
            self.pass_num = pixc_edge_reader.getAttValue("pass_number")
            try:
                self.continent = pixc_edge_reader.getAttValue("continent")
            except:
                pass
            self.nb_pix_range = pixc_edge_reader.getAttValue("nr_pixels")
            self.nb_pix_azimuth = pixc_edge_reader.getAttValue("nr_lines")

            # 3.2 - Initialize variables to numpy arrays
            
            # 3.2.1 - Edge objects info
            self.edge_loc = np.array(()).astype('int')  # Localization (top/bottom/both)
            self.edge_label = np.array(()).astype('int')  # Label in tile
            self.edge_idx = np.array(()).astype('int')  # Index in original L2_HR_PIXC product
            
            # 3.2.2 - Variables from L2_HR_PIXC product
            self.azimuth_idx = np.array(()).astype('int')  # Azimuth indices
            self.range_idx = np.array(()).astype('int')  # Range indices
            self.classif= np.array(())  # Classification
            self.pixel_area = np.array(())  # Pixel area
            self.latitude = np.array(())  # Latitude
            self.longitude = np.array(())  # Longitude
            self.height = np.array(())  # Height
            self.crosstrack = np.array(())  # Cross-track distance
            
            # 3.2.3 - Info of the nadir point associated to the PixC
            self.nadir_time = np.array(())  # Time
            self.nadir_longitude = np.array(())  # Longitude
            self.nadir_latitude = np.array(())  # Latitude
            self.nadir_x = np.array(())  # X cartesian coordinate
            self.nadir_y = np.array(())  # Y cartesian coordinate
            self.nadir_z = np.array(())  # Z cartesian coordinate
            self.nadir_vx = np.array(())  # Velocity in X coordinate
            self.nadir_vy = np.array(())  # Velocity in Y coordinate
            self.nadir_vz = np.array(())  # Velocity in Z coordinate
            self.near_range = np.array(())  # Near range distance, ie distance from satellite to the ground pixel the nearest of the nadir

        # 4 - Get number of pixels  
        OUT_nb_pix = pixc_edge_reader.getDimValue('record')

        # 5 - Update vectors if there are pixels
        if OUT_nb_pix > 0:
            
            # 5.1 - Edge objects info
            self.edge_label = np.concatenate((self.edge_label, pixc_edge_reader.getVarValue("edge_label")))
            self.edge_idx = np.concatenate((self.edge_idx, pixc_edge_reader.getVarValue("edge_idx")))
            self.edge_loc = np.concatenate((self.edge_loc, pixc_edge_reader.getVarValue("edge_loc")))
            
            # 5.2 - Variables from L2_HR_PIXC product
            self.azimuth_idx = np.concatenate((self.azimuth_idx, pixc_edge_reader.getVarValue("azimuth_index")))
            self.range_idx = np.concatenate((self.range_idx, pixc_edge_reader.getVarValue("range_index")))
            self.pixel_area = np.concatenate((self.pixel_area, pixc_edge_reader.getVarValue("pixel_area")))
            self.classif = np.concatenate((self.classif, pixc_edge_reader.getVarValue("classif")))
            self.latitude = np.concatenate((self.latitude, pixc_edge_reader.getVarValue("latitude")))
            self.longitude = np.concatenate((self.longitude, pixc_edge_reader.getVarValue("longitude")))
            self.height = np.concatenate((self.height, pixc_edge_reader.getVarValue("height")))
            self.crosstrack = np.concatenate((self.crosstrack, pixc_edge_reader.getVarValue("crosstrack")))
            
            # 5.3 - Info of the nadir point associated to the PixC
            self.nadir_time = np.concatenate((self.nadir_time, pixc_edge_reader.getVarValue("nadir_time")))
            self.nadir_longitude = np.concatenate((self.nadir_longitude, pixc_edge_reader.getVarValue("nadir_lon")))
            self.nadir_latitude = np.concatenate((self.nadir_latitude, pixc_edge_reader.getVarValue("nadir_lat")))
            self.nadir_x = np.concatenate((self.nadir_x, pixc_edge_reader.getVarValue("nadir_x")))
            self.nadir_y = np.concatenate((self.nadir_y, pixc_edge_reader.getVarValue("nadir_y")))
            self.nadir_z = np.concatenate((self.nadir_z, pixc_edge_reader.getVarValue("nadir_z")))
            self.nadir_vx = np.concatenate((self.nadir_vx, pixc_edge_reader.getVarValue("nadir_vx")))
            self.nadir_vy = np.concatenate((self.nadir_vy, pixc_edge_reader.getVarValue("nadir_vy")))
            self.nadir_vz = np.concatenate((self.nadir_vz, pixc_edge_reader.getVarValue("nadir_vz")))
            self.near_range = np.concatenate((self.near_range, pixc_edge_reader.getVarValue("nadir_near_range")))

        # 6 - Close file
        pixc_edge_reader.close()
        
        return OUT_nb_pix, OUT_tile_ref
        
    # ----------------------------------------

    def edgeGlobalRelabeling(self):
        """
        This function gives new labels to entities gathered at tile edges.
        """

        # 1 - Deal with all edges
        for i_edge in range(len(self.tile_ref)-1):  # Loop over tile edges
            
            my_api.printInfo("[PixelCloudSP] ***** Processing edge of tiles %s and %s *****" % (self.tile_ref[i_edge], self.tile_ref[i_edge+1]))

            # 1.1 - Get indices of pixels processed at the current edge
            tile_idx1 = np.where(self.tile_idx == i_edge)[0]
            tile_idx2 = np.where(self.tile_idx == i_edge+1)[0]

            # 1.2 - If one tile does not have pixel to process, continue to next iteration
            if (tile_idx1.size == 0) or (tile_idx2.size == 0):
                my_api.printInfo("")
                continue

            # 1.3 - New labels for pixels at the edge of current edge
            new_labels_subset = self.gatherEdgeEntities(tile_idx1, tile_idx2)
            
            # 1.4 - Link old labels to new ones
            self.labelMatching(tile_idx1, tile_idx2, new_labels_subset)
            
            my_api.printInfo("")

        # 2 - Deal with edge at beginning or end of pass
        if (self.labels == 0).any():
            
            # 2.1 - Get tiles with unprocessed labels (still initialized to zero)
            tiles_to_process = np.unique(self.tile_idx[np.where(self.labels == 0)])

            for tile in tiles_to_process:

                # 2.2 - Get associated tile index
                tile_idx = np.where(self.tile_idx == tile)[0]

                # 2.3 - Get old labels from PGE_LakeTile
                old_labels = np.unique(self.edge_label[tile_idx][np.where(self.labels[tile_idx] == 0)])

                # 2.4 - Compute a global new label
                new_labels = np.arange(old_labels.size) + np.max(self.labels) + 1

                # 2.5 - Update global labels
                for idx in range(old_labels.size):
                    self.labels[tile_idx[np.where(self.edge_label[tile_idx] == old_labels[idx])]] = new_labels[idx]

    # ----------------------------------------

    def gatherEdgeEntities(self, IN_tile_idx1, IN_tile_idx2):
        """
        This function gives new labels for pixels at current tile edge.
            1. Pixels within a buffer around tile edge are selected.
            2. A water mask is computed
            3. The mask is labeled
            
        :param IN_tile_idx1: indices of pixels edge of tile 1
        :type IN_tile_idx1: 1D-array of int
        :param IN_tile_idx2: indices of pixels edge of tile 2
        :type IN_tile_idx2: 1D-array of int
        
        :return: OUT_new_labels_subset = new labels given to pixels of edge entities
        :rtype: 1D-array of int
        """
        my_api.printDebug("[PixelCloudSP] == gatherEdgeEntities ==")

        # 1 - Pixels at top / bottom of tile are loaded
        rg, az = self.selectEdgePixels(IN_tile_idx1, IN_tile_idx2)

        # 2 - Equivalent matrix size in azimuth and range
        # 2.1 - Equivalent matrix size in azimuth and range
        nb_pix_range = max(rg) + 1
        nb_pix_azimuth = max(az) + 1
        # 2.2 - Compute water mask over the subset
        waterMask = my_tools.computeBinMat(nb_pix_range, nb_pix_azimuth, rg, az)

        # 3 - Label entities over the subset of PixC at the edge tile
        sepEntities, nb_obj = my_tools.labelRegion(waterMask)

        # 4 - Convert into 1D-array
        OUT_new_labels_subset = my_tools.convert2dMatIn1dVec(rg, az, sepEntities)

        my_api.printInfo("[PixelCloudSP] %d separate entities located at the edge tile" % np.unique(OUT_new_labels_subset).size)

        return OUT_new_labels_subset

    def selectEdgePixels(self, IN_tile_idx1, IN_tile_idx2):
        """
        This function selects pixels at top and bottom of tile 1 and 2
        
        :param IN_tile_idx1: indices of pixels edge of tile 1
        :type IN_tile_idx1: 1D-array of int
        :param IN_tile_idx2: indices of pixels edge of tile 2
        :type IN_tile_idx2: 1D-array of int
        
        :return: OUT_rg = range indices of pixels at the edge
        :rtype: 1D-array of int
        :return: OUT_az = azimuth indices of pixels at the edge
        :rtype: 1D-array of int
        """
        my_api.printDebug("[PixelCloudSP] == selectEdgePixels ==")

        # 1 - Distinguish top / bottom edge considering the pass orientation (ascending vs descending)
        # For ascending passes, top = North edge and bottom = South edge. For descending passes, it's the opposite.
        if self.ascending:
            rg1, az1 = self.getEdgePixels(IN_tile_idx1, "top")
            rg2, az2 = self.getEdgePixels(IN_tile_idx2, "bottom")
        else:
            rg1, az1 = self.getEdgePixels(IN_tile_idx1, "bottom")
            rg2, az2 = self.getEdgePixels(IN_tile_idx2, "top")

        # 2 - Concatenate pixels range in a numpy array
        rg = np.concatenate((rg1, rg2))

        # 3 - Concatenate pixels azimuth in a numpy array
        if self.ascending:
            # if ascending : tile1 have higher azimuth indices and tile 2 have lower azimuth indices
            az = np.concatenate((az1, az2 + max(az1)))
        else:
            az = np.concatenate((az1 + max(az2), az2))

        # 4 - Reduce range and azimuth values in order to reduce the size of generated water mask
        OUT_az = az - min(az)
        OUT_rg = rg - min(rg)

        # Return reduced range and azimuth indices
        return OUT_rg, OUT_az

    def getEdgePixels(self, IN_tile_idx, IN_edge_loc_str):
        """
        The function returns range and azimuth of pixels located within a buffer around the tile edge specified in edge_loc_str.
            
        :param IN_tile_idx: indices of pixels edge of tile
        :type IN_tile_idx: 1D array of int
        :param IN_edge_loc_str: edge location = "top" or "bottom"
        :type IN_edge_loc_str: string
        
        :return: range and azimuth indices of pixels
        :rtype: 1D array of int
        """

        idx_edge_buf = None
        
        # Associated location
        if IN_edge_loc_str == "bottom":
            # bottom of tile, get indices of all pixels with azimuth zero
            idx_edge_buf = np.where(self.azimuth_idx[IN_tile_idx] == 0)[0]
            
        elif IN_edge_loc_str == "top":
            # top of tile, get indices of all pixels with maximal azimuth
            az_max = max(self.azimuth_idx[IN_tile_idx])
            idx_edge_buf = np.where(self.azimuth_idx[IN_tile_idx] == az_max)[0]

        else:
            my_api.exitWithError("IN_edge_loc_str input variable has to be 'top' or 'bottom'")

        return self.range_idx[IN_tile_idx][idx_edge_buf], self.azimuth_idx[IN_tile_idx][idx_edge_buf]

    # ----------------------------------------

    def labelMatching(self, IN_tile_idx1, IN_tile_idx2, IN_new_labels_subset):
        """
        This function matches labels computed in LakeTile_edge file with labels computed in the gatherEdgeEntities() function.
        Pixels belonging to the same entity but cut by the tiling process are gathered with a single new label.
            
        :param IN_tile_idx1: indices of pixels in tile 1
        :type IN_tile_idx1: 1D-array of int
        :param IN_tile_idx2: indices of pixels in tile 2
        :type IN_tile_idx2: 1D-array of int
        :param IN_new_labels_subset: new labels computed in gatherEdgeEntities()
        :type IN_new_labels_subset: 1D-array of int
        """
        my_api.printDebug("[PixelCloudSP] Matching old labels from LakeTile_edge with new labels")

        # 1 - Get old labels computed in LakeTile_edge
        old_labels_subset1, old_labels_subset2 = self.selectEdgeLabels(IN_tile_idx1, IN_tile_idx2)

        # 2 - Get new labels
        # NB: the first part of IN_new_labels_subset contains labels of tile 1, the end of the array contains labels of tile 2
        IN_new_labels_subset1 = IN_new_labels_subset[:old_labels_subset1.size]
        IN_new_labels_subset2 = IN_new_labels_subset[old_labels_subset1.size:]

        # the structure of IN_new_labels_subset{1-2} and old_labels_subset{1-2} is the exactly the same, it contains the labels of all pixels within the subset defined the azimuth buffer

        IN_new_labels_subset_unique = np.unique(IN_new_labels_subset)

        # correspondance contains a list of tuples. Each tuple will correspond to a new entity and will contains the old labels of tiles 1 and 2
        correspondance = []

        # Matching the labels of the subset of tile 1 and 2
        for new_l in IN_new_labels_subset_unique:
            # get old labels (lake tile labels) of tiles 1 et 2
            old_l1 = np.unique(old_labels_subset1[np.where(IN_new_labels_subset1 == new_l)])
            old_l2 = np.unique(old_labels_subset2[np.where(IN_new_labels_subset2 == new_l)])

            # Most current case : at one label of tile1 corresponds one label of tile 2
            if old_l1.size == 1 and old_l2.size == 1:
                # appending to correspondance a tuple with old label 1 and old label 2
                correspondance.append((old_l1[0], old_l2[0]))

            # At one label of tile 1 do not corresponds a label of tile 2. The case happens when the lake is entierly located at the boundary of on tile bot does not cross the border.
            elif old_l2.size == 0:
                for idx1 in np.arange(old_l1.size):
                    correspondance.append((None, old_l1[idx1]))

            # At one label of tile 2 do not corresponds a label of tile 1.
            elif old_l1.size == 0:
                for idx2 in np.arange(old_l2.size):
                    correspondance.append((None, old_l2[idx2]))

            # Case that rarely occurs : the lake meanders along the border between two tiles, in this case severals labels from tile1 matches with several labels of tile2.
            else:
                for idx1 in np.arange(old_l1.size):
                    for idx2 in np.arange(old_l2.size):
                        correspondance.append((old_l1[idx1], old_l2[idx2]))

        # To give an explicit example, labels of tile1 of replaced by letters belonging to ['a', ..] and labels of tile2 are replaces by lettres belonging to [ 'r', ...]
        # correspondance : [(a, r), (a, s), (b, s), (b, t), (c, u), (d, u)]
        # labels a, b and r, s, t belong to the same entity, labels c, d, u belong to a separate entity.

        # The fonction lib.matchLabels returns reorganised correspondance labels :
        # label_matched : [([a, b], [r, s, t]), ([c, d], [u])]
        label_matched = matchLabels(correspondance)

        my_api.printDebug("[PixelCloudSP] > %d labels of first tile matched with %d labels of second tile into %d entities" % (np.unique(old_labels_subset1).size, np.unique(old_labels_subset2).size, len(label_matched)))

        # need new labels for global lake_sp processings
        unique_label = np.arange(len(label_matched)).astype('int') + max(self.labels) + 1
        # local labels from half edge tiles 1 and 2 are moved into new labels global only for Lake_sp processings
        # Ex : label_matched : [([a, b],[r,s,t]),([c,d],[u])]
        # The for loop iterates over entities at current tile edge
        for idx, label in enumerate(label_matched):
            # l contains the tuple corresponding to the idx^th entity
            # Ex : l : ([a, b], [r, s, t])

            # new_labels[idx] contains a label specefic for lake_SP processings
            new_label = unique_label[idx]

            for old_l1 in label[0]:
                # At one label of tile 2 do not corresponds a label of tile 1. The case happens when the lake is entierly located at the boundary of on tile bot does not cross the border.
                # In this case, old_l1 is setted to None, then, it is not processed.
                if old_l1:
                    # Get label of entity already computed in the case of a lake covering more than two tiles
                    labels_concerned = np.unique(self.labels[IN_tile_idx1][np.where(np.logical_and(self.edge_label[IN_tile_idx1] == old_l1, self.edge_loc[IN_tile_idx1] == 2))])
                    # Deleting label = 0 as those label are not already computed
                    labels_concerned = np.delete(labels_concerned, np.where(labels_concerned == 0))
                    # For each already processed and more than two tile lake, update the global labels
                    for label_to_relabel in labels_concerned:
                        self.labels[np.where(self.labels == label_to_relabel)] = new_label

                    # Set global label
                    self.labels[IN_tile_idx1[np.where(self.edge_label[IN_tile_idx1] == old_l1)]] = new_label

            for old_l2 in label[1]:
                # At one label of tile 1 do not corresponds a label of tile 2. The case happens when the lake is entierly located at the boundary of on tile bot does not cross the border.
                # In this case, old_l2 is setted to None, then, it is not processed.
                if old_l2:
                    # Get label of entity already computed in the case of a lake covering more than two tiles
                    labels_concerned = (np.unique(self.labels[IN_tile_idx2][np.where(np.logical_and(self.edge_label[IN_tile_idx2] == old_l2, self.edge_loc[IN_tile_idx2] == 2))]))
                    # Deleting label = 0 as those label are not already computed
                    labels_concerned = np.delete(labels_concerned, np.where(labels_concerned == 0))

                    # For each already processed and more than two tile lake, update the global labels
                    for label_to_relabel in labels_concerned:
                        self.labels[np.where(self.labels == label_to_relabel)] = new_label

                    # Set global label
                    self.labels[IN_tile_idx2[np.where(self.edge_label[IN_tile_idx2] == old_l2)]] = new_label

        nb_edge_entities = np.unique(np.concatenate((self.labels[IN_tile_idx1], self.labels[IN_tile_idx2]))).size

        my_api.printDebug("[PixelCloudSP] > %d separate entities are located at tile edge" % nb_edge_entities)

    def selectEdgeLabels(self, IN_tile_idx1, IN_tile_idx2):
        """
        This function selects old labels from LakeTile_edge at top and bottom of tile 1 and 2

        :param IN_tile_idx1: indices of edge pixels in tile 1
        :type IN_tile_idx1: 1D-array of int
        :param IN_tile_idx2: indices of edge pixels in tile 2
        :type IN_tile_idx2: 1D-array of int
        
        :return: labels of edge pixels at edge 1 and 2
        :rtype: 1D-array of int
        """
        
        # In ascending case, the edge is located at the top of tile1 and at the bottom of tile2. In descending case, it's the opposite.
        if self.ascending:
            label1 = self.getEdgeLabels(IN_tile_idx1, "top")
            label2 = self.getEdgeLabels(IN_tile_idx2, "bottom")
        else:
            label1 = self.getEdgeLabels(IN_tile_idx1, "bottom")
            label2 = self.getEdgeLabels(IN_tile_idx2, "top")

        return label1, label2

    def getEdgeLabels(self, IN_tile_idx, IN_edge_loc_str):
        """
        This function returns the LakeTile_edge labels of pixels within the buffer zone
        
        :param IN_tile_idx: indices of pixels at the edge of tile
        :type IN_tile_idx: 1D array of int
        :param IN_edge_loc_str: edge location = "top" or "bottom"
        :type IN_edge_loc_str: string
        
        :return: LakeTile_edge labels of pixels within the buffer zone
        :rtype: 1D-array of int
        """

        idx_edge_buf = None
        
        # Associated location
        if IN_edge_loc_str == "bottom":
            # Bottom of tile: get indices of all pixels with azimuth zero
            idx_edge_buf = np.where(self.azimuth_idx[IN_tile_idx] == 0)[0]
            
        elif IN_edge_loc_str == "top":
            # Top of tile, get indices of all pixels with maximal azimuth
            az_max = max(self.azimuth_idx[IN_tile_idx])
            idx_edge_buf = np.where(self.azimuth_idx[IN_tile_idx] == az_max)[0]

        else:
            my_api.exitWithError("IN_edge_loc_str input variable has to be 'top' or 'bottom'")

        return self.edge_label[IN_tile_idx[idx_edge_buf]]
        
    # ----------------------------------------

    def getAzimuthOfLake(self, IN_indices):
        """
            This function returns a re-computed azimuth index in order to have a continous azimuth along a lake at the edge of tiles

        :param IN_indices: indices of pixels of a lake
        :type IN_indices: 1D-array of int

        :return: recomputed azimuth_idx of the lake
        :rtype: 1D-array of int
        """

        tiles = np.unique(self.tile_idx[IN_indices])
        lake_tile_idx = self.tile_idx[IN_indices]
        lake_azimuth_idx = self.azimuth_idx[IN_indices]

        if self.pass_num % 2 == 0:  # ascending pass
            for tile in tiles:
                lake_azimuth_idx[np.where(lake_tile_idx > tile)] = lake_azimuth_idx[np.where(lake_tile_idx > tile)] + max(self.azimuth_idx[IN_indices][np.where(lake_tile_idx == tile)])+1
        else:
            for tile in tiles[::-1]:  # descending pass
                lake_azimuth_idx[np.where(lake_tile_idx > tile)] = lake_azimuth_idx[np.where(lake_tile_idx > tile)] + max(self.azimuth_idx[IN_indices][np.where(lake_tile_idx == tile)])+1

        return lake_azimuth_idx

    # ----------------------------------------

    def getMajorityPixelsTileRef(self, IN_label):
        """
        This fuction returns the tile reference of the tile containing the larger number of pixels with the given label.
            
        :param IN_label : labels of lake to process
        :type IN_label: int
        
        :return: tile reference
        :rtype: string
        """
        
        # 1 - Get unique values and counts of tiles for pixels with label IN_label
        # unique, counts = np.unique(self.tile_idx[np.where(self.labels == IN_label)], return_counts=True)
        unique, idx = np.unique(self.tile_idx[np.where(self.labels == IN_label)], return_inverse=True)
        counts = np.bincount(idx)

        # 2 - Get tile ref corresponding to the max number of pixels
        OUT_tile_max_pix = self.tile_ref[unique[np.where(counts == max(counts))][0]]

        return OUT_tile_max_pix
        
    # ----------------------------------------

    def getLakeTileLabel(self, IN_new_label):
        """
        This function is designed to retrieve old labels of PGE_LakeTile. The given new label corresponds to a global label, corresponding to several old labels.
        The old label involving the largest number of pixels is return.
            
        :param IN_new_label: global new label
        :type IN_new_label: int
        
        :return: LakeTile_edge label involving the largest number of pixels
        :rtype: string
        """

        # 1 - Get tiles concerned by the current new label
        tiles_concerned = np.unique(self.tile_idx[np.where(self.labels == IN_new_label)])

        nb_max_pix = 0
        OUT_final_label = 0

        for tile in tiles_concerned:

            # 2 - Get indices of IN_new_label pixels in tile
            label_idx = np.where(self.labels == IN_new_label)

            # 3 - Get lake_tile label value and number of pixels
            # unique, count = np.unique(self.edge_label[tile_idx], return_counts=True)
            unique, idx = np.unique(self.edge_label[label_idx], return_inverse=True)
            count = np.bincount(idx)
            np_max_pix_tmp = np.max(count)
            label_max_pix_tmp = unique[np.where(count == np_max_pix_tmp)][0]

            if np_max_pix_tmp > nb_max_pix:
                nb_max_pix = np_max_pix_tmp
                OUT_final_label = label_max_pix_tmp

        # 4 - Returns the lake tile label involving the largest number of pixels
        return str(OUT_final_label)


#######################################


def matchLabels(IN_liste):
    """
    This function reorganise labels in order to group labels by entities
        
    :param IN_liste: ex : [(a, r), (a, s), (b, s), (b, t), (c, u), (d, u)]. Labels a, b and r, s, t belong to the same entity, labels c, d, u belong to a separate entity.
    
    :return: labels gathered by entities
             ex : [(set([a, b]), set([r, s, t])), (set([c, d]), set([u]))] <=> [([a, b], [r, s, t]), ([c, d], [u])]
    """
    return groupBySecond(groupByFirst(IN_liste))


def groupByFirst(IN_liste):
    """
    This function take a list of tuples. The list is grouped by the first element of tuple and returned as a dictionary.
        
    :param IN_liste: la list of tuple. Ex : [ (a, r), (b, r),  (c, s), (c, t)]
    
    :return OUT_dico : ex : {a: set([r]), b: set([r]), c : set([t]}
    """

    OUT_dico = {}
    
    for (ind_i, ind_j) in IN_liste:
        if ind_i not in OUT_dico:
            OUT_dico[ind_i] = set()
        OUT_dico[ind_i].add(ind_j)

    return OUT_dico


def groupBySecond(IN_dict):
    """
    This function take a dictionary. The dictionary is grouped by second elements and returned as a list of tuple of set.
        
    :param IN_dict: result of group by first function: {a: set([r]), b: set([r]), c : set([t]}
    
    :return: a list of tuple of set. ex : [(set([a, b]), set([r])), (set([c]), set([s, t]))]
    """
    
    # Init
    OUT_results = []
    
    for key, value in list(IN_dict.items()):
        
        hasBeenFound = False
        
        for ind, result in enumerate(OUT_results):
            if checklistInter(value, result):
                hasBeenFound = True
                OUT_results[ind] = merge(result, key, value)

        if not hasBeenFound:
            OUT_results.append((set([key]), value))
            
    return OUT_results


def checklistInter(IN_value, IN_result):
    """
    Check if an element of value is contained in result
        
    :param IN_value: a set of elements ex : set([r])
    :param IN_result: a tuple of set ex : (set([a]), set([r]))
    
    :return: set([r])
    """
    return IN_result[1].intersection(IN_value)


def merge(IN_result, IN_key, IN_value):
    """
    Merge value with element of tuple key in result.
        
    :param IN_result: a tuple of set ex : (set([a]), set([r]))
    :param IN_key: key of first element of tuple. ex : b
    :param IN_value: set to merge ex : set([r])
    
    :return: merged tuple of set. ex :set([a, b]), set([r])
    """
    return (IN_result[0].union(set([IN_key])), IN_result[1].union(IN_value))
