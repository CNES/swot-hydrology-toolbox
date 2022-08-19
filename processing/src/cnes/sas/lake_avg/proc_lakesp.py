# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2021/10/27:Version initiale
# VERSION:2.0.0:DM:#91:2022/05/05:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: proc_lakeavg.py
    :synopsis: Deals with LakeAvg shapefile product
     Created on 2021/04/23

.. moduleauthor: Claire POTTIER (CNES DSO/SI/TR)

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import numpy as np
from osgeo import ogr

import cnes.common.service_config_file as service_config_file

import cnes.common.lib_lake.locnes_filenames as locnes_filenames


class LakeSPProduct(object):
    """
    class LakeSPProduct
    Manage LakeSP products
    """
    
    def __init__(self, in_basin_code):
        """
        Constructor

        :param in_basin_code: level-3 basin code (ie 3 digits = CBB)
        :type in_basin_code: int

        Variables of the object:
            - cfg / service_config_file.cfg: instance of LOCNES configuration file
            - obj_lake_db / lake_db.lakeDb_shp or lake_db.lakeDb_sqlite: lake database
            - content / LakeAvgProduct: container of the LakeAvg product
        """
        # Get instance of service config file
        self.cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # 1 - Init variables
        self.basin_code = in_basin_code
        self.list_att = ["time", "time_tai", "time_str", \
                         "wse", "wse_u", \
                         "area_total", "area_tot_u", \
                         "ds1_l", "ds1_l_u", "ds1_q", "ds1_q_u", \
                         "ds2_l", "ds2_l_u", "ds2_q", "ds2_q_u", \
                         "partial_f", "geoid_hght"]
        
        # 2 - Initialize LakeSP dictionnary to store product
        self.lakesp_archive = dict()
            
    def set_from_lakesp_files(self, in_lakesp_files):
        """
        Set variables from LakeSP shapefiles
        
        :param in_lakesp_files: list of full path of LakeSP_Prior shapefiles 
        :type in_lakesp_files: list
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        for cur_file in in_lakesp_files:
            logger.debug("INPUT file = %s" % cur_file)
            
            # 1 - Open shapefile in read-only access
            shp_driver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Shapefile driver
            lakesp_ds = shp_driver.Open(cur_file, 0)

            # 2 - Get the LakeSP layer
            lakesp_layer = lakesp_ds.GetLayer()
            nb_features_all = lakesp_layer.GetFeatureCount()
            
            # 3.1 - Filter wrt basin_code
            lakesp_layer.SetAttributeFilter("lake_id LIKE '{}%'".format(str(self.basin_code)))
            nb_features_begin = lakesp_layer.GetFeatureCount()
            logger.debug("> %d / %d PLD lakes located in basin_code %s" % (nb_features_begin, nb_features_all, str(self.basin_code)))
            lakesp_layer.SetAttributeFilter(None)
            # 3.2 - Filter wrt valid time, wse, and area_total
            lakesp_layer.SetAttributeFilter("lake_id LIKE '{}%' AND time > 0 AND wse > -1e10 AND area_total > 0".format(str(self.basin_code)))
            nb_features = lakesp_layer.GetFeatureCount()
            
            if nb_features  == 0:
                logger.debug("> No valid feature with basin_code = %s" % str(self.basin_code))
                
            else:
                logger.debug("> %d / %d features in basin_code %s are VALID" % (nb_features, nb_features_begin, str(self.basin_code)))
                
                # 4 - Store info
                
                # 4.1 - Retrieve pass number
                tmp_dict = locnes_filenames.get_info_from_filename(cur_file, "LakeSP")
                
                for cur_lake in lakesp_layer:
                    
                    # 4.2 - Retrieve lake_id
                    lake_id = cur_lake.GetField(str("lake_id"))
                        
                    # 4.3 - Create dict for lake_id if it doesn't exist
                    if lake_id not in self.lakesp_archive.keys():
                        self.lakesp_archive[lake_id] = dict()
                        self.lakesp_archive[lake_id]["pass"] = list()
                        self.lakesp_archive[lake_id]["geom"] = list()
                        for cur_att in self.list_att:
                            self.lakesp_archive[lake_id][cur_att] = list()
                    
                    # 4.4 - Update the lists
                    # TODO: delete if loop when LakeSP Issue #265 is corrected
                    # Issue #265 = [LakeSP] Certains lake_id sont en doublon dans le fichier _Prior
                    if tmp_dict["pass"] not in self.lakesp_archive[lake_id]["pass"]:
                        self.lakesp_archive[lake_id]["pass"].append(tmp_dict["pass"])
                        cur_geom = cur_lake.GetGeometryRef()
                        self.lakesp_archive[lake_id]["geom"].append(cur_geom.Clone())
                        for cur_att_name in self.list_att:
                            cur_att_value = cur_lake.GetField(str(cur_att_name))
                            if (cur_att_value == float) and (cur_att_value < -1e10):
                                cur_att_value = np.nan
                            self.lakesp_archive[lake_id][cur_att_name].append(cur_att_value)
                    else:
                        logger.warning("ISSUE in LakeSP product: more than 1 feature for PLD lake %s in the current LakeSP product" % lake_id)

            # 5 - Close shapefile
            lakesp_ds.Destroy()
                        
