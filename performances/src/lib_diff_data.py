#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2023/10/19: initial version
# FIN-HISTORIQUE
# ======================================================
"""
.. module lib_diff_data.py
    :synopsis: Class for getting DiffData, from:
    - metrics.txt file for SWOT simulated data
    - hydro database, for SWOT real data
    Convert input data to a common format: json

.. moduleauthor: CNES (L. Labat-Allee)

Copyright (c) 2023 CNES. All rights reserved.
"""
import logging
FV_REAL = -999999999999

class DiffData(object):
    """
    Class DiffData
    """

    def __init__(self):
        """
        Constructor
        """

        self.json = []

    def from_lakesp_prior_diff_table(self, gdf_diff, fields, which_area, metrix_name):
        """
        Get data from hydro database (table LakeSP_Prior_diff)
        :param gdf_diff: geodataframe
        :param fields: line fields
        :param which_area: 'total', 'detected' or 'geom'
        :param metrix_name: 'reference_surface_area_from_db' or 'reference_wse_from_db'
        :return: list of dictionary (json format)
        """

        # Read GeoDataFrame
        diff_observations = gdf_diff

        # Reading line by line from GeoDataFrame
        for line in diff_observations.index:
            # Create a dictionary for all variables
            # Note that Total_surface_area is set to Reference_surface_area (not used to compute relative error)
            # Note that Detected_WSE and Reference_WSE are not used to compute absolute error
            if metrix_name == 'reference_surface_area_from_db':
                relative_error = str(diff_observations.loc[line, "area_{}".format(which_area)])
                absolute_error = "None"
                wse_levelled = False  # parameter is not used for surface area stats
            else:
                # reference_wse_from_db
                relative_error = "None"
                absolute_error = str(diff_observations.loc[line, "wse"])
                wse_levelled = bool(diff_observations.loc[line, "wse_levelled"])

            dict_fields = {fields[0]: "None",  # scene
                           fields[1]: str(diff_observations.loc[line, "cycle_number"]),  # cycle
                           fields[2]: str(diff_observations.loc[line, "pass_number"]),  # pass
                           fields[3]: "None",  # tile
                           fields[4]: "WC",  # case
                           # fields[5]: str(line[0]),  # lake_id_compared (index 0 is feature_id) API v1.6
                           # fields[6]: str(line[0]),  # lake_id_ref
                           fields[5]: str(diff_observations.loc[line, "_fuid"]),  # lake_id_compared API v1.7
                           fields[6]: str(diff_observations.loc[line, "_fuid"]),  # lake_id_ref
                           fields[7]: str(int(diff_observations.loc[line, "xtrk_dist"] * 1e-3)),  # XLocMean in km
                           fields[8]: str(diff_observations.loc[line, "p_ref_area"]),  # Total_surface_area
                           fields[9]: str(diff_observations.loc[line, "p_ref_area"]),  # Reference_surface_area
                           fields[10]: relative_error,  # Relative_error
                           fields[11]: 999,  # Detected_WSE
                           fields[12]: 999,  # Reference_WSE
                           fields[13]: absolute_error,  # Absolute_error
                           fields[14]: str(diff_observations.loc[line, "dark_frac"]),  # dark_frac
                           fields[15]: wse_levelled  # wse_levelled
                           }
            self.json.append(dict_fields)

    def from_riversp_reach_diff_table(self, diff_observations, fields, which_area='total'):
        """
        Get data from hydro database (RiverSP_reach_diff table)
        :param diff_observations: dataframe
        :param fields: line fields
        :param which_area: 'total', 'detected' or 'geom'
        :return: list of dictionary (json format)
        """

        # Reading line by line from DataFrame
        for line in range(len(diff_observations)):
            # Create a dictionary for all variables
            # Note that Total|Detected|Reference_surface_area are not used to compute relative error
            dict_fields = {fields[0]: "None",  # scene
                           fields[1]: str(diff_observations.loc[line, "cycle_number"]),  # cycle
                           fields[2]: str(diff_observations.loc[line, "pass_number"]),  # pass
                           fields[3]: "None",  # tile
                           fields[4]: "WC",  # case
                           fields[5]: str(diff_observations.loc[line, "_fuid"]),  # reach
                           fields[6]: str(int(diff_observations.loc[line, "xtrk_dist"] * 1e-3)),  # XLocMean in km
                           fields[7]: 999,  # length (not used)
                           fields[8]: str(diff_observations.loc[line, "p_width"]),  # p_width
                           fields[9]: 999,  # Total|Detected_surface_area (not used)
                           fields[10]: 999,  # Reference_surface_area (not used)
                           fields[11]: str(diff_observations.loc[line, "area_{}".format(which_area)]),  # Relative_error
                           fields[12]: str(diff_observations.loc[line, "dark_frac"])  # dark_frac
                           }
            self.json.append(dict_fields)

    def filter_levelled(self):
        """
        Filter data that are not nivelled        
        """
        logger = logging.getLogger(self.__class__.__name__)
        nb_lines = len(self.json)
        self.json = list(filter(lambda x: bool(x['wse_levelled']) == True, self.json))
        print("Keep %d / %d lines after filter over filter_nivelled" % (len(self.json), nb_lines))

    def convert_absolute_error_into_relative(self):
        # to be optimized
        """
        Compute relative error by removing height mean value  
        """
        logger = logging.getLogger(self.__class__.__name__)    
        dict_mean = {}
        for (indl, cur_line) in enumerate(self.json):
            if cur_line["lake_id_ref"] not in dict_mean :
                dict_mean[cur_line["lake_id_ref"]]= [float(cur_line["Absolute_error"]), 1]
            else:
                dict_mean[cur_line["lake_id_ref"]] = [dict_mean[cur_line["lake_id_ref"]][0] + float(cur_line["Absolute_error"]), dict_mean[cur_line["lake_id_ref"]][1] + 1]
            
        # compute mean height for a given serie
        for key in dict_mean.keys():   
            dict_mean[key] = dict_mean[key][0]/float(dict_mean[key][1])
        
        # for each absolute error wse, remove mean height (name remains absolute error, that is not really good)
        for (indl, cur_line) in enumerate(self.json):
            cur_line["Absolute_error"] = str(float(cur_line["Absolute_error"]) - dict_mean[cur_line["lake_id_ref"]])
            

if __name__ == "__main__":
    # Exampple of input metrics file
    stat_file = "/work/EXPERT_CENTER/mwec/workspace/HR/labatl/scripts/perfos/output/lakes/nominal-nominal_wrt_truth_20231019/stats/metrics_on_nominal-nominal_wrt_truth_20231019_tot_Prior.txt"

    # Example of fields (for lake statistics)
    fields = ['scene', 'cycle', 'pass', 'tile', 'case', 'lake_id_compared', 'lake_id_ref', 'XLocMean',
              'Total_surface_area', 'Reference_surface_area', 'Relative_error',
              'Detected_WSE', 'Reference_WSE', 'Absolute_error', 'dark_frac']

    diff_data = DiffData(fields)

    diff_data.from_file(stat_file)
