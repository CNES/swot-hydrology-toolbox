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
import pandas as pd
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

    def from_file(self, csv_path, fields, metrix_name):
        """
        Get data from a CSV file and store it in a JSON format.

        :param csv_path: Path to the CSV file
        :param fields: List of field names for the output dictionary
        :param metrix_name: 'wse' or 'area_total' or 'area_geom'
        :return: list of dictionary (json format)
        """

        # Read CSV file into a DataFrame
        df = pd.read_csv(csv_path)

        # Iterate over each line in the DataFrame
        for line in df.index:
            # Create a dictionary for all variables
            if metrix_name.startswith("area_"):
                relative_error = str(df.loc[line, f"{metrix_name}"])
                absolute_error = "None"
                wse_levelled = False  # parameter is not used for surface area stats
            else:
                # wse
                relative_error = "None"
                absolute_error = str(df.loc[line, "wse"])
                wse_levelled = bool(df.loc[line, "wse_levelled"])

            dict_fields = {
                fields[0]: str(df.loc[line, "cycle_number"]),  # cycle
                fields[1]: str(df.loc[line, "pass_number"]),  # pass
                fields[2]: str(df.loc[line, "_fuid"]),  # lake_id
                fields[3]: str(int(df.loc[line, "xtrk_dist"] * 1e-3)),  # XLocMean in km
                fields[4]: str(df.loc[line, "p_ref_area"]),  # Reference_surface_area
                fields[5]: relative_error,  # Relative_error
                fields[6]: absolute_error,  # Absolute_error
                fields[7]: str(df.loc[line, "dark_frac"]),  # dark_frac
                fields[8]: wse_levelled  # wse_levelled
            }
            self.json.append(dict_fields)

    def from_lakesp_prior_diff_table(self, gdf_diff, fields, metrix_name):
        """
        Get data from hydro database (table LakeSP_Prior_diff)
        :param gdf_diff: geodataframe
        :param fields: line fields
        :param metrix_name: 'wse' or 'area_total' or 'area_geom'
        :return: list of dictionary (json format)
        """

        # Read GeoDataFrame
        diff_observations = gdf_diff

        # Reading line by line from GeoDataFrame
        for line in diff_observations.index:
            # Create a dictionary for all variables
            if metrix_name.startswith("area_"):
                relative_error = str(diff_observations.loc[line, "{}".format(metrix_name)])
                absolute_error = "None"
                wse_levelled = False  # parameter is not used for surface area stats
            else:
                # wse
                relative_error = "None"
                absolute_error = str(diff_observations.loc[line, "wse"])
                wse_levelled = bool(diff_observations.loc[line, "wse_levelled"])

            dict_fields = {fields[0]: str(diff_observations.loc[line, "cycle_number"]),  # cycle
                           fields[1]: str(diff_observations.loc[line, "pass_number"]),  # pass
                           fields[2]: str(diff_observations.loc[line, "_fuid"]),  # lake_id
                           fields[3]: str(int(diff_observations.loc[line, "xtrk_dist"] * 1e-3)),  # XLocMean in km
                           fields[4]: str(diff_observations.loc[line, "p_ref_area"]),  # Reference_surface_area
                           fields[5]: relative_error,  # Relative_error
                           fields[6]: absolute_error,  # Absolute_error
                           fields[7]: str(diff_observations.loc[line, "dark_frac"]),  # dark_frac
                           fields[8]: wse_levelled  # wse_levelled
                           }
            self.json.append(dict_fields)

    def filter_levelled(self):
        """
        Filter data that are not nivelled        
        """
        logger = logging.getLogger(self.__class__.__name__)
        nb_lines = len(self.json)
        self.json = list(filter(lambda x: bool(x['wse_levelled']) == True, self.json))
        print("Keep %d / %d lines after filter over nivelled insitu data" % (len(self.json), nb_lines))

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

    # Example of csv loading file
    stat_file = "/work/EXPERT_CENTER/mwec/workspace/HR/labatl/scripts/swot-hydrology-toolbox/performances/examples/" \
                "input/diff_area.csv"

    # Example of fields (for lake statistics)
    fields = ['cycle', 'pass', 'lake_id', 'XLocMean', 'Reference_surface_area', 'Relative_error',
              'Absolute_error', 'dark_frac', 'wse_levelled']

    diff_data = DiffData(fields)

    diff_data.from_file(stat_file, fields, 'total', 'area')
