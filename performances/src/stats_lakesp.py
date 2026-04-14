# -*- coding: utf-8 -*-
"""
 This file is part of the SWOT Hydrology Toolbox
 Copyright (c) 2026 CNES
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY,
 read LICENSE.txt for further details.
"""

import argparse
import pandas as pd
import os
from datetime import date

import logging
LOGGER = logging.getLogger(__name__)

from src.lake_stats_manager import LakeStatsManager


def plot_stats(working_dir, df_data, which_stats="wse", which_perfos="goal"):
    """
    Compute statistics and save errors distribution figures
    :param working_dir: output working directory to store output figures
    :param df_data: diff dataframe
    :param which_stats: 'wse' or 'area_total' or 'area_geom'
    :param which_perfos: to filter over lake size: 'goal'=lakes bigger than 100mx100m; 'req'=lakes bigger than 250mx250m
    :return:
    """
    
    # 0 - Create working directory if it doesn't exist
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    
    # 1 - Init LakeStatsManager object
    try:
        stats_manager = LakeStatsManager(working_dir,
                                         metrix_name=which_stats,
                                         which_perfos=which_perfos,
                                         gdf=df_data)
    except ValueError:
        LOGGER.error("CANNOT compute stats")
        return

    # 2 - Compute statistics
    stats_manager.compute_statistics()

    
#################################################
    
def main():
    parser = argparse.ArgumentParser(description="Compute statistics for LakeSP product")
    parser.add_argument('--working_dir', help='working directory', nargs='?', type=str,
                        default="/work/EXPERT_CENTER/mwec/workspace/HR/labatl/scripts/swot-hydrology-toolbox"
                                "/performances/examples/output")
    parser.add_argument('--which_stats', help='which variable need to be evaluate: wse or area_total',
                        choices=['wse', 'area_total', 'area_geom'], nargs='?', type=str, default="wse")
    parser.add_argument('--which_perfos', help='compare the statistics with the scientific specifications: '
                                               'requirement or goal',
                        choices=['goal', 'req'], nargs='?', type=str, default="req")
    parser.add_argument('--csv', help='input file containing diff data', nargs='?', type=str)
    args = vars(parser.parse_args())

    LOGGER.info("===== stats_lakesp = BEGIN =====")

    LOGGER.info("> Working directory = %s" % args["working_dir"])
    LOGGER.info("> Compute the statistics on %s variable" % args["which_stats"].upper())
    LOGGER.info("> Compare the statistics with the scientific specifications : %s" % args["which_perfos"].upper())

    # Get csv diff dataset
    df_diff = pd.read_csv(args["csv"], sep=";")

    # Get today date
    today = date.today()
    cur_date = today.strftime("%Y%m%d")

    # Output directory to save statistics images to the disk
    output_dir = os.path.join(args["working_dir"], "lakes_perfos_{}_{}".format(args["which_stats"], cur_date))

    # Run statistics
    plot_stats(output_dir, df_diff, which_stats=args["which_stats"], which_perfos=args["which_perfos"])

    LOGGER.info("===== stats_lakesp = END =====")


if __name__ == "__main__":

    main()


    