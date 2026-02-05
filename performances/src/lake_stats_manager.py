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


def get_tags(which_area="total"):
    """
    Get lake tags list
    :param which_area: area type to be compared among 'detected' and 'total' (detected + dark water)
    :type which_area: string
    :return: dict
    """

    area = which_area.capitalize()  # 'Detected' or 'Total'

    # Define tags to find wanted info in the lake metrics file

    # Define tags list when performance assessment is based on surface area:
    # Detected/Total surface area, reference surface area, reference lake size, polygon-centroid position
    area_tags = ('{}_surface_area'.format(area), 'Reference_surface_area', 'Reference_surface_area', 'XLocMean',
                 'Relative_error')

    # Define tags list when performance assessment is based on WSE
    # Detected WSE, reference WSE, reference lake size, polygon-centroid position
    wse_tags = ('Detected_WSE', 'Reference_WSE', 'Reference_surface_area', 'XLocMean', 'Absolute_error')

    # Define a dictionary containing begin and end tags
    lake_tags_dict = {'reference_surface_area': area_tags,
                      'reference_surface_area_from_db': area_tags,
                      'reference_wse': wse_tags,
                      'reference_wse_from_db': wse_tags}

    # Fields of a line of metrics file
    fields = ['scene', 'cycle', 'pass', 'tile', 'case', 'lake_id_compared', 'lake_id_ref', 'XLocMean',
              '{}_surface_area'.format(area), 'Reference_surface_area', 'Relative_error',
              'Detected_WSE', 'Reference_WSE', 'Absolute_error', 'dark_frac', 'wse_levelled']

    return lake_tags_dict, fields


#######################################


class LakeStatsManager(StatsManager):
    """
    Class for managing statistics on lakes
    """

    def __init__(self, working_dir,
                 metrix_name='reference_surface_area',
                 subdataset=None,
                 subdataset_filter=None,
                 which_perfos="req",
                 threshold=None,
                 by_lake_size=None,
                 along_swath=None,
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
        :param subdataset: to filter over a single scene (value = scene_number),
                                      or  a single tile (value = scene_number;<pppp>_<tttt>)
                                      or wrt a listing of tile quality
        :type subdataset: string
        :param subdataset_filter: to filter over tile quality if subdataset is a listing file; ex: 0-1, 0-1-3
        :type subdataset_filter: string
        :param which_perfos: to filter over lake size: 'goal'=lakes bigger than 100mx100m; 'req'=lakes bigger than 250mx250m (default)
        :type which_perfos: string
        :param threshold: to filter abnormal values
        :type threshold: float
        :param by_lake_size: to plot errors one sigma distribution as a function of lake size
        :type by_lake_size: string
        :param along_swath: to plot errors one sigma distribution along the swath width
        :type along_swath: string
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
        lake_tags_dict, fields = get_tags(which_area)

        # Get diff data
        diff_data = DiffData()
        # print("Get data from hydro database")
        diff_data.from_lakesp_prior_diff_table(gdf, fields, which_area, metrix_name)

        # To filter
        if metrix_name == 'reference_wse_from_db':
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
                              by_water_body_size=by_lake_size,
                              water_body='lake size',
                              along_swath=along_swath,
                              which_perfos=which_perfos,
                              which_area=which_area,
                              filter_by_scene_tile='all',
                              outlier_lim=outlier_lim,
                              save_to_disk=output_dir)

        # Plot number of lakes
        df_filtered = pd.DataFrame(diff_data.json)

        print("Number of PLD lakes: {}".format(df_filtered['lake_id_ref'].str.split(':').str[1].nunique()))
        print("Number of matchups: {}".format(len(df_filtered['lake_id_ref'])))

        df_filtered["lake_id"] = df_filtered["lake_id_ref"].str.extract(r":(\d+)$")[0].astype("Int64")
        df_data_unique = df_filtered.drop_duplicates(subset=["lake_id"])

        # Save to figure. Example: fig_surface-area_scene-all_swath-nominal_5.png
        fig = plot_nb_obj_as_a_function_of_water_body_size(np.array(df_data_unique['Reference_surface_area'].astype(float)) * 1e6,
                                                             "lake size", begin_slices_list, end_slices_list)
        fig_name = 'fig_{}_scene-{}_swath-{}_4'.format(metrix_name.replace("reference_", "").replace("_from_db","").replace("_", "-"),
                                                       filter_by_scene_tile, "nominal")
        fig.savefig(os.path.join(output_dir, fig_name))
        print("Number of lakes figure is saved to the disk! Here: {0}".format(os.path.join(output_dir, fig_name)))
        fig.clf()