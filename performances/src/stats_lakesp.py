# -*- coding: utf-8 -*-
"""
.. module:: stats_lakesp.py
    :synopsis: compute performances over LakeSP product
    Created on Fri May 24 08:47:18 2024

.. moduleauthor: Claire POTTIER (CNES DTN/TPI/TR), Lucie LABAT-ALLEE (CS-GROUP)
"""
import argparse
import getpass
import os
import shutil

import logging
LOGGER = logging.getLogger(__name__)

from src.lake_stats_manager import LakeStatsManager

#################################################


def plot_stats(working_dir, df_data, which_stats="wse", which_perfos="goal", flag_levelled=False):
    """

    :param working_dir:
    :param df_data:
    :param which_stats:
    :param which_perfos:
    :param flag_levelled: if True, only keep levelled data for stats / if False, use every data and compute
    relative error (remove mean) - option not used for SURFACE AREA perfos
    :return:
    """
    
    # 0 - Create working directory if it doesn't exist
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    
    # 1 - Plot histo with nb of lakes wrt lake size
    if (which_stats == "wse") and flag_levelled:
        df_data_tmp = df_data[df_data["wse_levelled"]]
        df_data_uniq = df_data_tmp.drop_duplicates(subset=["lake_id"])
    else:
        df_data_uniq = df_data.drop_duplicates(subset=["lake_id"])

    # 2.1 - Build input parameters
    if which_stats == "surface-area":
        metrix_name = "reference_surface_area"
        which_area = "total"  # Area to be compared with: 'total', 'detected' or 'geom' for stats on AREA

    else:
        metrix_name = "reference_wse"
        which_area = "prior"  # Area to be compared with: 'prior' for stats on WSE

    flag_by_lake_size = "yes"  # to plot errors one sigma distribution as a function of lake size
    threshold = None  # Option to filter abnormal values. Ex: 1
    flag_along_swath = "yes"  # to plot errors one sigma distribution along the swath width
    outlier_lim = None  # Outlier limits to clip on graphs (area_limit (%), wse_limit (m)). Ex: (200, 1)

    # 2.2 - Init LakeStatsManager object
    try:
        stats_manager = LakeStatsManager(working_dir,
                                         metrix_name=metrix_name,
                                         which_perfos=which_perfos,
                                         threshold=threshold,
                                         by_lake_size=flag_by_lake_size,
                                         along_swath=flag_along_swath,
                                         which_area=which_area,
                                         outlier_lim=outlier_lim,
                                         levelled=flag_levelled,
                                         gdf=df_data)
    except ValueError:
        LOGGER.error("CANNOT compute stats")
        return
    
    # 2.3 - Compute statistics
    stats_manager.compute_statistics()
    
    # 3 - Move and changes filenames
    for ind, new_name in enumerate(["nb_diff_per_size", "1s_error_per_size", "error_per_xtrack", "med_error_per_xtrack", "nb_lakes_per_size"]):
        os.replace(os.path.join(working_dir, which_perfos, "fig_%s_scene-all_swath-nominal_%d.png" % (which_stats, ind)),
                   os.path.join(working_dir, "%d_%s.png" % (ind+1, new_name)))
    shutil.rmtree(os.path.join(working_dir, which_perfos))
    
    
#################################################
    

if __name__ == "__main__":

    # TODO: à revoir
    parser = argparse.ArgumentParser(description="Compute perfos per source")
    parser.add_argument('--working_dir', help='full path of working directory', nargs='?', type=str, default="/work/EXPERT_CENTER/mwec/workspace/HR/%s/perfos_lakesp" % getpass.getuser())
    parser.add_argument('--which_stats', help='variable for perfors among wse (default) and area_total', choices=['wse', 'surface_area'], nargs='?', type=str, default="wse")
    parser.add_argument('--input_csv', help='input CSV with diff with prefarea', nargs='?', type=str)
    args = vars(parser.parse_args())
    
    LOGGER.info("===== RUN perfos_lakesp =====")
    LOGGER.info("> working directory = %s" % args["working_dir"])
    LOGGER.info("> variable on which to compute perfos = %s" % args["which_stats"])
    LOGGER.info("===== END OF perfos_lakesp =====")
    