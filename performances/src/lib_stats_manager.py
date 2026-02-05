#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2020/09/11: initial version
# FIN-HISTORIQUE
# ======================================================
"""
.. module lib_stats_manager.py
    :synopsis: Class for handling statistics computation and visualization

.. moduleauthor: CNES (L. Labat-Allee)

Copyright (c) 2020 CNES. All rights reserved.
"""
import logging
import os
import matplotlib as mpl

if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
from src.lib_relative_errors import RelativeErrors, plot_nb_obj_as_a_function_of_water_body_size
import pandas as pd

print_debug = False


######################################


class StatsManager(object):
    """
    Main class for managing statistics (rivers and lakes)
    """

    def __init__(self,
                 lines,
                 metrix_name='area',
                 threshold=None,
                 tags_dict=None,
                 begin_slices_list=None,
                 end_slices_list=None,
                 water_body=None,
                 which_perfos='req',
                 which_area='total',
                 outlier_lim=None,
                 save_to_disk=None
                 ):
        """
        Initialize statistics manager

        Args:
            lines (str): jon text file lines
            metrix_name (str): matrix name to consider: 'area' or 'wse' (only for lakes)
            threshold (float): to filter abnormal values
            which_perfos (str): to filter perfos over lake size / river width: 'req' (requirements) or 'goal'
            which_area (str): area type to be compared among 'detected' and 'total' (detected + dark water)
            outlier_lim (list): outlier limits to clip
        """
        logger = logging.getLogger(self.__class__.__name__)

        # Diff data stored into a json structure
        self._lines = lines

        # Exit if no line left
        if len(self._lines) == 0:
            logger.error("No info left to compute statistics")
            sys.exit(0)

        self._metrix_name = metrix_name
        self._threshold = threshold
        self._water_body = water_body
        self._tags_dict = tags_dict
        self._begin_slices_list = begin_slices_list
        self._end_slices_list = end_slices_list
        self._which_perfos = which_perfos
        self._which_area = which_area
        self._outlier_lim = outlier_lim
        if save_to_disk is not None:
            if os.path.exists(save_to_disk):
                self._dirname = save_to_disk
        self._save_to_disk = save_to_disk

    # -------------------------------------------

    def compute_statistics(self):
        """
        Compute statistics and visualize errors distribution
        """
        logger = logging.getLogger(self.__class__.__name__)

        # Initialize list of tags to extract
        tags = self._tags_dict[self._metrix_name]
        tags_values_list = []

        # Loop on lines, read all tag values
        for t, tag in enumerate(tags):
            # print("tags[{0}]: {1}".format(t, tag))
            tag_values_list = self.extract_values(tag)
            tags_values_list.append(tag_values_list)

        # to check nan values
        is_nan = 0

        # to check that there are values into one slice at least
        is_one_slice_filled = 0

        # Need to check that all list are filled
        for tv, tag_values in enumerate(tags_values_list):
            if len(tag_values) != 0:
                is_nan = 1

        # Initialization
        water_body_size = None
        xrange = None
        fig_list = []

        if is_nan != 0:

            try:
                # river width or lake size
                water_body_size = tags_values_list[0]
                if self._water_body == 'lake size':
                    # lake surface area are in km2 => convert in m2
                    water_body_size = [i * 1e6 for i in water_body_size]
                xrange = tags_values_list[1]
                # relative_error (for AREA) or absolute_error (for WSE)
                errors = np.asarray(tags_values_list[2])
            except:
                logger.error("Error when retrieving tag''s values with metrix from {0}".format(self._metrix_name))

            # Type of errors
            if self._metrix_name == 'wse':
                error_type = "Absolute"
                unit = "m"
            else:
                error_type = "Relative"
                unit = "%"
                
            #  Initialize RelativeErrors object
            relative_errors = RelativeErrors(errors, metrix_name=self._metrix_name, water_bodies_size=water_body_size,
                                             xrange=xrange)
            # Compute relative error(s) one sigma
            relative_errors.compute_error_one_sigma(self._threshold, self._begin_slices_list,
                                                    self._end_slices_list, self._water_body, self._which_perfos)

            # Check that one slice is filled
            for s, error in enumerate(relative_errors._errors_one_sigma[2:]):
                logger.info("{0} error one sigma, slice {1}: {2:.3f} {3}".format(error_type,
                                                                                 relative_errors._slices_title_list[
                                                                                     s],
                                                                                 relative_errors._errors_one_sigma[
                                                                                     s + 1], unit))
                if error != 0:
                    is_one_slice_filled = 1

            # Bar plot of relative errors one sigma as a function of water body size (river width or lake surface)
            if is_one_slice_filled != 0:
                fig = relative_errors.plot_nb_obs_as_a_function_of_water_body_size(self._water_body)
                fig_list.append(fig)
                fig = relative_errors.plot_one_sigma_error_distribution(self._water_body, self._which_perfos,
                                                                        self._which_area)
                fig_list.append(fig)

            # Statistics on rivers or lakes
            # For perf assessment based on RiverObs (using reaches), plot relative errors distribution along swath's width
            fig = relative_errors.plot_along_swath(self._water_body, self._which_area,
                                                   outlier_lim=self._outlier_lim)
            fig_list.append(fig)
            fig = relative_errors.plot_along_swath_popular(self._threshold, self._which_area)
            fig_list.append(fig)

            # Plot number of lakes
            df_filtered = pd.DataFrame(self._lines)
            print("Number of PLD lakes: {}".format(df_filtered['lake_id'].str.split(':').str[1].nunique()))
            print("Number of matchups: {}".format(len(df_filtered['lake_id'])))

            df_filtered["lake_id_subset"] = df_filtered["lake_id"].str.extract(r":(\d+)$")[0].astype("Int64")
            df_data_unique = df_filtered.drop_duplicates(subset=["lake_id_subset"])

            fig = plot_nb_obj_as_a_function_of_water_body_size(
                np.array(df_data_unique['Reference_surface_area'].astype(float)) * 1e6,
                "lake size", self._begin_slices_list, self._end_slices_list)
            fig_list.append(fig)

            if self._save_to_disk is None:
                # Show figures
                plt.show()
            else:
                # Save figures to the disk
                new_names = ["nb_diff_per_size", "1s_error_per_size", "error_per_xtrack", "med_error_per_xtrack",
                             "nb_lakes_per_size"]
                for i, (fig, new_name) in enumerate(zip(fig_list, new_names)):
                    fig_name = f"{i+1}_{new_name}.png"
                    fig.savefig(os.path.join(self._dirname, fig_name))
                    print(f"Figure saved: {fig_name} in {self._dirname}")
                    plt.close(fig)

        else:
            logger.info("No data found. Error one sigma is not computed.")

    def extract_values(self, tag):
        """
        Extract tag values, for all lines of json

        Args:
            tag (str): tag

        Returns:
            (list): list of values
        """
        logger = logging.getLogger(self.__class__.__name__)

        list_values = []

        # Loop on lines, read relative error values
        for i, line in enumerate(self._lines):
            value = line[tag]
            if value != "NA":
                extract = np.float64(value)
                list_values.append(extract)
            else:
                logger.warning("Value for <%s> is set to NA in line %d" % (tag, i))

        return list_values
