# -*- coding: utf-8 -*-
"""
.. module:: stats_lakesp.py
    :synopsis: compute performances over LakeSP product
    Created on Mon Feb 9 08:47:18 2026

.. moduleauthor: Claire POTTIER (CNES DTN/TPI/TR), Lucie LABAT-ALLEE (CS-GROUP)
"""
import argparse
import getpass
import os

import logging
LOGGER = logging.getLogger(__name__)

from src.lake_stats_manager import LakeStatsManager


def plot_stats(working_dir, df_data, which_stats="wse", which_perfos="goal"):
    """

    :param working_dir: output working directory to store output figures
    :param df_data: diff dataframe
    :param which_stats: 'wse' or 'area_total' or 'area_geom'
    :param which_perfos: to filter over lake size: 'goal'=lakes bigger than 100mx100m; 'req'=lakes bigger than 250mx250m
    :return:
    """
    
    # 0 - Create working directory if it doesn't exist
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    
    # 2.1 - Init LakeStatsManager object
    try:
        stats_manager = LakeStatsManager(working_dir,
                                         metrix_name=which_stats,
                                         which_perfos=which_perfos,
                                         gdf=df_data)
    except ValueError:
        LOGGER.error("CANNOT compute stats")
        return

    # 2.2 - Compute statistics
    stats_manager.compute_statistics()

    
#################################################
    

if __name__ == "__main__":

    # TODO: à revoir
    parser = argparse.ArgumentParser(description="Compute perfos per source")
    parser.add_argument('--working_dir', help='full path of working directory', nargs='?', type=str, default="/work/EXPERT_CENTER/mwec/workspace/HR/%s/perfos_lakesp" % getpass.getuser())
    parser.add_argument('--which_stats', help='variable for perfos among wse (default) and area_total', choices=['wse', 'area_total', 'area_geom'], nargs='?', type=str, default="wse")
    parser.add_argument('--input_csv', help='input CSV with diff with prefarea', nargs='?', type=str)
    args = vars(parser.parse_args())
    
    LOGGER.info("===== RUN perfos_lakesp =====")
    LOGGER.info("> working directory = %s" % args["working_dir"])
    LOGGER.info("> variable on which to compute perfos = %s" % args["which_stats"])
    LOGGER.info("===== END OF perfos_lakesp =====")
    