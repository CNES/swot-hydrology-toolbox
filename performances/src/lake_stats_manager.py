#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Class for handling lakes statistics computation and visualization

Copyright (c) 2020 CNES.
"""
import os
import numpy as np
from src.lib_stats_manager import StatsManager
from src.lib_diff_data import DiffData
from src.lib_relative_errors import plot_nb_obj_as_a_function_of_water_body_size
import pandas as pd

# Define slices to plot errors one sigma evolution depending on lake size
begin_slices_list_4goal = [10000, 62500, 250000, 1000000, 4000000]  # m2
end_slices_list_4goal = [62500, 250000, 1000000, 4000000, "+"]  # m2
begin_slices_list_4req = [62500, 250000, 1000000, 4000000]  # m2
end_slices_list_4req = [250000, 1000000, 4000000, "+"]  # m2
# begin_slices_list = [0.01, 0.0625, 0.25, 1, 4]  # km2
# end_slices_list = [0.0625, 0.25, 1, 4, 25]  # km2


def get_tags():

    """
    Get lake tags list
    :return: dict
    """

    # Define tags to find wanted info in the lake metrics file

    # Define tags list when performance assessment is based on surface area:
    # PLD reference lake size, polygon-centroid position, relative error
    area_tags = ('Reference_surface_area', 'XLocMean', 'Relative_error')

    # Define tags list when performance assessment is based on WSE
    # PLD reference lake size, polygon-centroid position, absolute error
    wse_tags = ('Reference_surface_area', 'XLocMean', 'Absolute_error')

    # Define a dictionary containing begin and end tags
    lake_tags_dict = {'reference_surface_area': area_tags,
                      'reference_wse': wse_tags
                      }

    # Fields of a line of metrics file
    fields = ['cycle', 'pass', 'lake_id', 'XLocMean', 'Reference_surface_area', 'Relative_error',
              'Absolute_error', 'dark_frac', 'wse_levelled']

    return lake_tags_dict, fields


#######################################


class LakeStatsManager(StatsManager):
    """
    Class for managing statistics on lakes
    """

    def __init__(self, working_dir,
                 metrix_name='reference_surface_area',
                 which_perfos="req",
                 threshold=None,
                 which_area="total",
                 outlier_lim=None,
                 gdf=None,
                 levelled=True):
        """
        Initialize statistics manager
        
        :param working_dir: working directory
        :type working_dir: string
        :param metrix_name: metrics name to consider, among ('reference_surface_area', 'reference_wse')
        :type metrix_name: string
        :param which_perfos: to filter over lake size: 'goal'=lakes bigger than 100mx100m; 'req'=lakes bigger than 250mx250m (default)
        :type which_perfos: string
        :param threshold: to filter abnormal values
        :type threshold: float
        :param which_area: area type to be compared among 'detected' and 'total' (detected + dark water)
        :type which_area: string
        :param outlier_lim: outlier limits to clip (area_limit, wse_limit)
        :type outlier_lim: list
        :param gdf: if DiffData comes from Hydro database, read gdf
        :type gdf: GeoDataFrame
        :param levelled: if levelled = True, only keep levelled data in stats, if False keep every data but remove mean(wse) of a given lake_id for each absolute_wse value 
        :type gdf: GeoDataFrame
        """

        # Init slices depending on which_perfos
        begin_slices_list = begin_slices_list_4req
        end_slices_list = end_slices_list_4req
        if which_perfos == 'goal':
            begin_slices_list = begin_slices_list_4goal
            end_slices_list = end_slices_list_4goal

        # Create output directory if it doesn't exist
        output_dir = os.path.join(working_dir, which_perfos)
        if not os.path.exists(output_dir):
            print("[WARNING] %s directory doesn't exist => create it" % output_dir)
            os.makedirs(output_dir)

        # Update surface area type
        lake_tags_dict, fields = get_tags()

        # Get diff data
        diff_data = DiffData()
        # print("Get data from hydro database")
        diff_data.from_lakesp_prior_diff_table(gdf, fields, which_area, metrix_name)

        # To filter
        if metrix_name == 'reference_wse':
            # to keep only nivelled stations
            if levelled:
                diff_data.filter_levelled()
            else:
                diff_data.convert_absolute_error_into_relative()

        # Initialize the base class
        StatsManager.__init__(self,
                              lines=diff_data.json,
                              metrix_name=metrix_name,
                              threshold=threshold,
                              tags_dict=lake_tags_dict,
                              begin_slices_list=begin_slices_list,
                              end_slices_list=end_slices_list,
                              water_body='lake size',
                              which_perfos=which_perfos,
                              which_area=which_area,
                              filter_by_scene_tile='all',
                              outlier_lim=outlier_lim,
                              save_to_disk=output_dir)

        # Plot number of lakes
        df_filtered = pd.DataFrame(diff_data.json)

        print("Number of PLD lakes: {}".format(df_filtered['lake_id'].str.split(':').str[1].nunique()))
        print("Number of matchups: {}".format(len(df_filtered['lake_id'])))

        df_filtered["lake_id_subset"] = df_filtered["lake_id"].str.extract(r":(\d+)$")[0].astype("Int64")
        df_data_unique = df_filtered.drop_duplicates(subset=["lake_id_subset"])

        # Save to figure. Example: fig_surface-area_scene-all_swath-nominal_5.png
        fig = plot_nb_obj_as_a_function_of_water_body_size(np.array(df_data_unique['Reference_surface_area'].astype(float)) * 1e6,
                                                             "lake size", begin_slices_list, end_slices_list)
        fig_name = 'fig_{}_scene-{}_swath-{}_4'.format(metrix_name.replace("reference_", "").replace("_", "-"),
                                                       'all', "nominal")
        fig.savefig(os.path.join(output_dir, fig_name))
        print("Number of lakes figure is saved to the disk! Here: {0}".format(os.path.join(output_dir, fig_name)))
        fig.clf()