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
import json
import os
import datetime

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

    def from_file(self, stat_file, fields):
        """
        Get data from metrics file .txt
        :param stat_file: text file (output of run/run_lakes_estimation.py)
        :param fields: line fields
        :return: list of dictionary (json format)
        """

        logger = logging.getLogger(self.__class__.__name__)

        # Read txt file and store data into a json structure
        with open(stat_file) as fh:

            for line in fh:

                # reading line by line from the text file
                description = line.strip().split(' - ', 15)

                # for output see below
                # print(description)

                # loop variable
                i = 0
                # intermediate dictionary
                dict_fields = {}
                while i < len(fields):
                    # creating dictionary for variable
                    value = description[i].split(None, 3)
                    dict_fields[fields[i]] = value[1]
                    i += 1

                # logger.debug('dict_fields[{}]: {}'.format(i, dict_fields))

                self.json.append(dict_fields)

        print(len(self.json))

        # creating json file
        cur_date = datetime.datetime.now().strftime("%Y%m%d")
        username = os.environ.get("USERNAME") 
        tmp_dir = "/work/scratch/data/{}/lakes_perfos_{}/stats/".format(username, cur_date)
        if not os.path.exists(tmp_dir):
            print("[WARNING] %s directory doesn't exist => create it" % tmp_dir)
            os.makedirs(tmp_dir)
        out_file = open(os.path.join(tmp_dir, "metrics_from_file.json"), "w")
        json.dump(self.json, out_file, indent=4)
        out_file.close()

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
            
    def filter_wse_fill_value(self):
        """
        Filter data where Detected WSE = -999999999999 (fill_value)
        """
        logger = logging.getLogger(self.__class__.__name__)

        nb_lines = len(self.json)
        try:
            self.json = list(filter(lambda x: float(x['Detected_WSE']) == FV_REAL, self.json))
            logger.info("Keep %d / %d lines after filter WSE fill value: %d" % (len(self.json), nb_lines,
                                                                                FV_REAL))
        except:
            pass

    def filter_scenes_tiles(self, subdataset=None, subdataset_filter=None):
        """
        Filter data by:
        - specific scene IF subdataset represents a scene number or name
        - specific tile IF subdataset = <scene_number>;<pppp>_<tttt>
        - listing of scenes and tiles IF subdataset is an existing file; in this case,
        filter (scene, tile) pairs with quality > subdataset_filter
        Available quality values =
        0 : No particular problems (especially favorable case)
        1 : Moderate phenomenology-related difficulties
        2 : Major phenomenology-related difficulties (severe layover, only near/far range...)
        3 : Simulation/truth problems with limited impact
        4 : Simulation/truth problems with major impact
        :param subdataset: to filter over a single scene (value = scene_number),
        a single tile (value = scene_number;<pppp>_<tttt>) or wrt a listing of tile quality or by basin_id
        :type subdataset: str
        :param subdataset_filter: to filter over tile quality if subdataset is a listing file; ex: 0-1, 0-1-3
        :type subdataset_filter: str
        """
        logger = logging.getLogger(self.__class__.__name__)

        # By default, no filter
        filter_by_scene_tile = "all"

        # 0 - Filter by scene and tile
        if subdataset is not None:
            # 1 - Build dict listing tiles to process for all specified scenes
            list_scene_tiles = dict()
            if os.path.isfile(subdataset):
                list_filter = subdataset_filter.split("-")
                cpt_keep = 0
                cpt_remove = 0
                reader = open(subdataset, "r")
                for cur_line in reader:
                    tmp_split = cur_line.split(";")
                    scene_number = tmp_split[0]
                    tile_ref = tmp_split[1]
                    quality = tmp_split[2][0]
                    if quality in list_filter:
                        logger.debug("Keep in listing / scene = %s and tile = %s" % (scene_number, tile_ref))
                        cpt_keep += 1
                        if scene_number in list_scene_tiles.keys():
                            list_scene_tiles[scene_number] += ";%s" % tile_ref
                        else:
                            list_scene_tiles[scene_number] = "%s" % tile_ref
                    else:
                        logger.debug("Remove from listing / scene = %s and tile = %s" % (scene_number, tile_ref))
                        cpt_remove += 1
                reader.close()
                filter_by_scene_tile = "listing-filter-%s" % subdataset_filter
                logger.info("Keep %d / %d (scene, tile) pairs" % (cpt_keep, cpt_keep + cpt_remove))

            elif ";" in subdataset:
                scene_number = subdataset.split(";")[0]
                tile_number = subdataset.split(";")[1]
                list_scene_tiles[scene_number] = tile_number
                filter_by_scene_tile = "%s_tile-%s" % (str(scene_number), str(tile_number))

            else:
                list_scene_tiles[subdataset] = "all"
                filter_by_scene_tile = str(subdataset)

            # 2 - Keep wanted tiles
            tmp_json = []  # Init list of indices of lines to keep

            # 2.1 - Loop over all lines

            for (indl, cur_line) in enumerate(self.json):
                scene_number = cur_line['scene']
                if scene_number in list_scene_tiles.keys():
                    if list_scene_tiles[scene_number] != "all":
                        tile_number = "{}_{}".format(cur_line['pass'], cur_line['tile'])
                        if tile_number not in list_scene_tiles[scene_number]:
                            logger.debug(
                                "Remove line with scene number %s and tile %s = %s" % (scene_number, tile_number, cur_line))
                        else:
                            logger.debug(
                                "KEEP line with scene number %s and tile %s = %s" % (scene_number, tile_number, cur_line))
                            tmp_json.append(cur_line)
                    else:
                        logger.debug("KEEP line with scene number %s = %s" % (scene_number, cur_line))
                        tmp_json.append(cur_line)
                else:
                    logger.debug("Remove line with scene number %s = %s" % (scene_number, cur_line))

            logger.info("Keep %d / %d lines after filter over scenes and tiles " % (len(tmp_json), len(self.json)))

            # 2.2 - Keep only selected lines
            self.json = tmp_json

        return filter_by_scene_tile




if __name__ == "__main__":
    # Exampple of input metrics file
    stat_file = "/work/EXPERT_CENTER/mwec/workspace/HR/labatl/scripts/perfos/output/lakes/nominal-nominal_wrt_truth_20231019/stats/metrics_on_nominal-nominal_wrt_truth_20231019_tot_Prior.txt"

    # Example of fields (for lake statistics)
    fields = ['scene', 'cycle', 'pass', 'tile', 'case', 'lake_id_compared', 'lake_id_ref', 'XLocMean',
              'Total_surface_area', 'Reference_surface_area', 'Relative_error',
              'Detected_WSE', 'Reference_WSE', 'Absolute_error', 'dark_frac']

    diff_data = DiffData(fields)

    diff_data.from_file(stat_file)
