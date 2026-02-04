# -*- coding: utf-8 -*-
"""
.. module:: stats_lakesp.py
    :synopsis: compute performances over LakeSP product
    Created on Fri May 24 08:47:18 2024

.. moduleauthor: Claire POTTIER (CNES DTN/TPI/TR), Lucie LABAT-ALLEE (CS-GROUP)
"""
import argparse
import getpass
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

import logging
LOGGER = logging.getLogger(__name__)

from src.lake_stats_manager import LakeStatsManager

#################################################


def plot_nb_lakes_as_a_function_of_water_body_size(water_bodies_size_tab, which_perfos="goal"):
    
    # 1.1 - Prepare data and X-axis title
    if which_perfos == "goal":
        begin_slices_list = [10000, 62500, 250000, 1000000, 4000000]  # m2
        end_slices_list = [62500, 250000, 1000000, 4000000, "+"]  # m2
    else:
        begin_slices_list = [62500, 250000, 1000000, 4000000]  # m2
        end_slices_list = [250000, 1000000, 4000000, "+"]  # m2
    nb_slices = len(begin_slices_list)
    _slices_title_list = []
    tab_nb_items = np.zeros(nb_slices)
    for s, cur_slice in enumerate(begin_slices_list):
        
        # Build list of indices corresponding to slices
        index_sup = (water_bodies_size_tab > int(cur_slice))
        if (s == nb_slices-1) and (end_slices_list[s] == "+"):
            water_bodies_size_selection = index_sup
        else:
            # Select water_bodies_size_tab index where size value is <= end slice
            index_inf = (water_bodies_size_tab <= int(end_slices_list[s]))
            # Select water_bodies_size_tab index where both conditions are gathered
            water_bodies_size_selection = np.logical_and(index_sup, index_inf)
        tab_nb_items[s] = len(np.where(water_bodies_size_selection)[0])
        
        # Build text of water body size (for plotting X labels)
        if (s == nb_slices-1) and (end_slices_list[s] == "+"):
            slice_title = '[>{0}{1}]'.format(int(np.sqrt(cur_slice)),'\u00B2')  # for example: >1000\u00B2
        else:
            slice_title = '[{0}{1}-{2:d}{3}]'.format(int(np.sqrt(cur_slice)),'\u00B2', int(np.sqrt(end_slices_list[s])),'\u00B2')  # for example: [50\u00B2-100\u00B2]
        # Add slice title to the list
        _slices_title_list.append(slice_title)
    nb_total = np.sum(tab_nb_items)
    if nb_total == 0:
        return None
    
    # 1.2 - Compute related array of x indexes
    ind_x = range(nb_slices)
    
    # 2 - Init figure
    fig = plt.figure(figsize=(10, 4))
    plt.ylim(0, np.max(tab_nb_items)*1.1)
    
    # 3 - Plot histogram of nb of items vs water body size slice
    barlist = plt.bar(ind_x, tab_nb_items)
    plt.xticks(ind_x, _slices_title_list)
    # Define a different color for each water bodies size
    if len(_slices_title_list) == 4:
        colorbar = {0: 'darkviolet', 1: 'violet', 2: 'darkorange', 3: 'gold'}
    else:
        colorbar = {0: 'darkblue', 1: 'darkviolet', 2: 'violet', 3: 'darkorange', 4: 'gold', 5: 'darkkhaki'}
    for ind in ind_x:
        barlist[ind].set_color(colorbar[ind])
    
    # 4 - Add title and labels
    # 4.1 - Title
    title = '%d lakes with matching in situ data' % nb_total
    plt.title(title, fontsize=20)
    # 4.2 - Label X
    water_body = 'lake size'
    water_body_unit = 'm'+'\u00B2'
    xlabel='{0} ({1})'.format(water_body.capitalize(), water_body_unit)
    plt.xlabel(xlabel, fontsize=14)
    # 4.3 - Label Y
    plt.ylabel('Number of lakes', fontsize=14)
    # 4.4 - Nb and associated percentage
    for ind, value in enumerate(tab_nb_items):
        #print_str = "{0:d} (={1:.0f}%)".format(int(value), round(value/nb_total*100.0))
        print_str = "{0:.0f}%".format(round(value/nb_total*100.0))
        if max(tab_nb_items) <= 100:
            coeff = 0.1
        else:
            coeff = 1
        plt.text(ind-0.1, value+coeff, print_str, color=colorbar[ind], fontweight='bold')
        
    return fig
    

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
    fig = plot_nb_lakes_as_a_function_of_water_body_size(np.array(df_data_uniq["p_ref_area"])*1e6, which_perfos=which_perfos)
    # Save as png
    if fig is not None:
        fig_name = '0_nb_lakes_per_size.png'
        fig.savefig(os.path.join(working_dir, fig_name))
        fig.clf()
    
    # 2.1 - Build input parameters
    if which_stats == "surface-area":
        metrix_name = "reference_surface_area_from_db"
        which_area = "total"  # Area to be compared with: 'total', 'detected' or 'geom' for stats on AREA

    else:
        metrix_name = "reference_wse_from_db"
        which_area = "prior"  # Area to be compared with: 'prior' for stats on WSE

    dark_frac_filter = None  # Option to filter over fractional area of dark water. Ex: 0.9, 0.8, etc.
    basin_id_filter = None  # Option to filter by basin_id. Ex: dark_frac_filter='783'
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
                                         dark_frac_filter=dark_frac_filter,
                                         outlier_lim=outlier_lim,
                                         basin_id_filter=basin_id_filter,
                                         levelled=flag_levelled,
                                         gdf=df_data)
    except ValueError:
        LOGGER.error("CANNOT compute stats")
        return
    
    # 2.3 - Compute statistics
    stats_manager.compute_statistics()
    
    # 3 - Move and changes filenames
    for ind, new_name in enumerate(["nb_diff_per_size", "1s_error_per_size", "error_per_xtrack", "med_error_per_xtrack"]):
        os.replace(os.path.join(working_dir, which_perfos, "fig_%s_scene-all_swath-nominal_%d.png" % (which_stats, ind)), os.path.join(working_dir, "%d_%s.png" % (ind+1, new_name)))
    shutil.rmtree(os.path.join(working_dir, which_perfos))
    
    
#################################################
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Compute perfos per source")
    parser.add_argument('--working_dir', help='full path of working directory', nargs='?', type=str, default="/work/EXPERT_CENTER/mwec/workspace/HR/%s/perfos_lakesp" % getpass.getuser())
    parser.add_argument('--which_stats', help='variable for perfors among wse (default) and area_total', choices=['wse', 'surface_area'], nargs='?', type=str, default="wse")
    parser.add_argument('--input_csv', help='input CSV with diff with prefarea', nargs='?', type=str)
    args = vars(parser.parse_args())
    
    LOGGER.info("===== RUN perfos_lakesp =====")
    LOGGER.info("> working directory = %s" % args["working_dir"])
    LOGGER.info("> variable on which to compute perfos = %s" % args["which_stats"])

    # MAJ
    compute_perfos_from_csv_file(args["working_dir"], args["input_csv"])
    
    LOGGER.info("===== END OF perfos_lakesp =====")
    