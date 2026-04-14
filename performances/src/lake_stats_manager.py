#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 This file is part of the SWOT Hydrology Toolbox
 Copyright (c) 2026 CNES
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY,
 read LICENSE.txt for further details.
"""

from src.lib_stats_manager import StatsManager
from src.lib_diff_data import DiffData

# Define slices to plot errors one sigma evolution depending on lake size
# GOAL
begin_slices_list_4goal = [10000, 62500, 250000, 1000000, 4000000]  # m2
end_slices_list_4goal = [62500, 250000, 1000000, 4000000, "+"]  # m2
# REQ (for requirement)
begin_slices_list_4req = [62500, 250000, 1000000, 4000000]  # m2
end_slices_list_4req = [250000, 1000000, 4000000, "+"]  # m2


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
    lake_tags_dict = {'area_total': area_tags, 'area_geom': area_tags, 'wse': wse_tags}

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
                 metrix_name='area_total',
                 which_perfos="req",
                 threshold=None,
                 outlier_lim=None,
                 gdf=None,
                 levelled=True):
        """
        Initialize statistics manager
        
        :param working_dir: working directory
        :type working_dir: string
        :param metrix_name: metrics name to consider, among 'wse' or 'area_total' or 'area_geom'
        :type metrix_name: string
        :param which_perfos: to filter over lake size: 'goal'=lakes bigger than 100mx100m
                                                        'req'=lakes bigger than 250mx250m (default)
        :type which_perfos: string
        :param threshold: to filter abnormal values
        :type threshold: float
        :param outlier_lim: outlier limits to clip (area_limit, wse_limit)
        :type outlier_lim: list
        :param gdf: if DiffData comes from Hydro database, read gdf
        :type gdf: GeoDataFrame
        :param levelled: if levelled = True, only keep levelled insitu data in stats
                         if levelled = False, keep every data but remove mean(wse) of a given lake_id
                         for each absolute_wse value
        :type gdf: GeoDataFrame
        """

        # Init slices depending on which_perfos
        begin_slices_list = begin_slices_list_4req
        end_slices_list = end_slices_list_4req
        if which_perfos == 'goal':
            begin_slices_list = begin_slices_list_4goal
            end_slices_list = end_slices_list_4goal

        # Get tags
        lake_tags_dict, fields = get_tags()

        # Initialize diff data
        diff_data = DiffData()
        diff_data.from_lakesp_prior_diff_table(gdf, fields, metrix_name)

        # Filter levelled stations
        if metrix_name == 'wse':
            if levelled:
                diff_data.filter_levelled()
            else:
                diff_data.convert_absolute_error_into_relative()

        # Initialize statistics manager
        StatsManager.__init__(self,
                              lines=diff_data.json,
                              metrix_name=metrix_name,
                              threshold=threshold,
                              tags_dict=lake_tags_dict,
                              begin_slices_list=begin_slices_list,
                              end_slices_list=end_slices_list,
                              water_body='lake size',
                              which_perfos=which_perfos,
                              outlier_lim=outlier_lim,
                              save_to_disk=working_dir)