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
import os
from osgeo import ogr

import cnes.common.service_config_file as service_config_file

import cnes.common.lib.my_shp_file as my_shp
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_variables as my_var
import cnes.common.lib_lake.lake_db as lake_db
import cnes.common.lib_lake.locnes_products_shapefile as shp_file


class LakeAvgProduct(object):
    """
    class LakeAvgProduct
    Manage LakeAvg products
    """
    
    def __init__(self, in_obj_lakedb, in_obj_lakesp, in_layer_name):
        """
        Constructor

        :param in_obj_lakedb: lake database
        :type in_obj_lakedb: lake_db.lakeDb_shp or lake_db.lakeDb_sqlite
        :param in_obj_lakesp: needed data retrieved from LakeSP products related to this LakeAvg product
        :type in_obj_lakesp: proc_lakesp.LakeSPProduct
        :param in_layer_name: name for LakeAvg layer
        :type in_layer_name: string

        Variables of the object:
            - cfg / service_config_file.cfg: instance of LOCNES configuration file
            - obj_lakedb / lake_db.lakeDb_shp or lake_db.lakeDb_sqlite: lake database
            - content / LakeAvgProduct: container of the LakeAvg product
        """
        # Get instance of service config file
        self.cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # 1 - Init variable
        # Lake database object
        self.obj_lakedb = in_obj_lakedb
        # LakeSP input files
        self.obj_lakesp = in_obj_lakesp
        
        # 2 - Initialize lake product content
        self.content = shp_file.LakeAvgProduct(in_layer_name)        
        
    def free_memory(self):
        """
        Destroy memory layers
        """
        self.content.free()

    # ----------------------------------------

    def compute_lake_features(self):
        """
        Computes LakeAvg features for features retrieved in related LakeSP products.
        """
        logger = logging.getLogger(self.__class__.__name__)

        for cur_lakeid, cur_lakesp_archive in self.obj_lakesp.lakesp_archive.items():
            logger.debug("===== Deal with PLD lake %s =====" % cur_lakeid)

            if self.obj_lakedb.lake_layer is not None:
                # 1 - Remove from the list of PLD located over the area
                try:
                    self.obj_lakedb.list_lakeid.remove(cur_lakeid)
                except:
                    msg = "STRANGE... PLD lake with lake_id=%s is not in the list => removed from LakeAvg product" % cur_lakeid
                    logger.warning(msg)
                    continue
            
            # 2.1 - Create prior lake object
            obj_plake = lake_db.PriorLake(self.obj_lakedb, cur_lakeid)
            # 2.2 - Format PLD attributes (i.e. convert None to required _FillValue if needed)
            pld_attributes = obj_plake.format_attributes()
            
            # 3 - Compute averaged geometry and attributes of PLD feature
            avg_geom, avg_attributes = self.form_avg_feature(cur_lakesp_archive)
            
            # 4 - Storage change
            if obj_plake.ok_to_compute_stocc:
                # 4.1 - Compute averaged storage change with direct approach
                direct_storage_change_values = self.compute_direct_storage_change(obj_plake, avg_attributes)  
                # 4.2 - Compute averaged storage change with incremental approach
                incremental_storage_change_values = self.compute_incremental_storage_change(obj_plake, avg_attributes)                   
                # 4.3 - Add prior feature to _Prior layer
                self.content.add_feature(avg_geom, {**avg_attributes, **pld_attributes, **direct_storage_change_values,\
                                                    **incremental_storage_change_values})
            else:
                logger.debug("Not able to compute averaged storage change")
                # 4.3 - Add prior feature to _Prior layer
                self.content.add_feature(avg_geom, {**avg_attributes, **pld_attributes})
            
        # 5 - Deal with PLD lakes which have been observed by SWOT during the cycle (if asked)
        if self.cfg.getboolean("CONFIG_PARAMS", "ADD_ALL"):
            self.add_pld_features_not_observed()
        
    # ----------------------------------------
    # Functions specific to lake feature (i.e.
    # geometry + attributes) computation
    # ----------------------------------------
    
    def form_avg_feature(self, in_lakesp_archive):
        """
        Process LakeSP features related to current PLD lake
        to build the LakeAvg feature boundary and associated attributes
        
        :param in_lakesp_archive: LakeSP archive related to the current PLD lake
        :type in_lakesp_archive: dict
        
        :return: out_geom = geometry of LakeAvg feature
        :rtype: out_geom = OGRPolygon
        :return: out_attributes = lake attributes (None for each item which could not be computed)
        :rtype: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Compute _hmin/_hmed/_hmax attributes
        out_attributes_hxxx = self.compute_hxxx_attributes(in_lakesp_archive)
        
        # 2 - Compute _avg attributes, except area and storage change
        out_attributes_avg = self.compute_avg_attributes(in_lakesp_archive)
        logger.debug("PLD lake observed {} times during this cycle".format(out_attributes_avg["npass"]))
        logger.debug("{} FULL passes = {}".format(out_attributes_avg["npass_full"], out_attributes_avg["pass_full"]))
        logger.debug("{} PARTIAL passes = {}".format(out_attributes_avg["npass_part"], out_attributes_avg["pass_part"]))
        
        # 3 - Build feature boundary and set area
        out_geom, out_attributes_geom = self.build_avg_boundary(in_lakesp_archive, out_attributes_avg)
        
        return out_geom, {**out_attributes_hxxx, **out_attributes_avg, **out_attributes_geom}
    
    def compute_hxxx_attributes(self, in_lakesp_archive):
        """
        Compute _hmin/_hmed/_hmax attributes related to current PLD lake
        
        :param in_lakesp_archive: LakeSP archive related to the current PLD lake
        :type in_lakesp_archive: dict

        :return: out_attributes = LakeAvg attributes (None for each item which could not be computed)
        :rtype: dict
        """
        
        nb_pass = len(in_lakesp_archive["pass"])
        
        if nb_pass == 1:
            out_attributes = compute_avg_values("med", in_lakesp_archive)
            
        elif nb_pass == 2:
            att_min = compute_avg_values("min", in_lakesp_archive)
            att_max = compute_avg_values("max", in_lakesp_archive)
            out_attributes = {**att_min, **att_max}
            
        else:
            att_min = compute_avg_values("min", in_lakesp_archive)
            att_med = compute_avg_values("med", in_lakesp_archive)
            att_max = compute_avg_values("max", in_lakesp_archive)
            out_attributes = {**att_min, **att_med, **att_max}
            
        return out_attributes
        
    def compute_avg_attributes(self, in_lakesp_archive):
        """
        Compute pass and averaged attributes, except area and storage change
        
        :param in_lakesp_archive: LakeSP archive related to the current PLD lake
        :type in_lakesp_archive: dict

        :return: out_attributes = LakeAvg attributes (None for each item which could not be computed)
        :rtype: dict
        """
        
        # Init output dictionary
        out_attributes = dict()
        
        # 1 - Attributes related to pass
        
        # 1.1 - Compute pass_full and pass_part lists
        pass_full = None
        pass_part = None
        for num_pass, partial_f in zip(in_lakesp_archive["pass"], in_lakesp_archive["partial_f"]):
            if partial_f == 0:
                # pass_full
                if pass_full is None:
                    pass_full  = "%s" % num_pass
                else:
                    pass_full += ";%s" % num_pass
            else:
                # pass_part
                if pass_part is None:
                    pass_part = "%s" % num_pass
                else:
                    pass_part += ";%s" % num_pass
                
        # 1.2 - Set related attributes
        # Related to pass_full
        out_attributes["pass_full"] = pass_full
        if pass_full is None:
            out_attributes["npass_full"] = 0
        else:
            out_attributes["npass_full"] = len(pass_full.split(";"))
        # Related to pass_part
        out_attributes["pass_part"] = pass_part
        if pass_part is None:
            out_attributes["npass_part"] = 0
        else:
            out_attributes["npass_part"] = len(pass_part.split(";"))
            
        # 1.3 - Number of passes
        out_attributes["npass"] = out_attributes["npass_full"] + out_attributes["npass_part"]
        
        # 1.4 - Deduce partial_f
        out_attributes["partial_f"] = 0
        if out_attributes["npass_full"] == 0:
            out_attributes["partial_f"] = 1
            
        # 2 - Attributes related to time
        out_attributes["t_avg"] = my_tools.value_or_none(np.nanmean(in_lakesp_archive["time"]))
        out_attributes["t_tai_avg"] = my_tools.value_or_none(np.nanmean(in_lakesp_archive["time_tai"]))
        if out_attributes["t_avg"] is None:
            out_attributes["t_str_avg"] = None
        else:
            out_attributes["t_str_avg"] = my_tools.convert_utc_to_str(out_attributes["t_avg"])  # Time in UTC as a string
            
        # 3 - Averaged wse
        out_attributes["wse_avg"] = my_tools.value_or_none(np.nanmean(in_lakesp_archive["wse"]))
        out_attributes["wse_avg_u"] = my_tools.value_or_none(np.nanmean(in_lakesp_archive["wse_u"]))

        # 4 - Compute other values
        # 4.1 - quality_f
        # 4.2 - geoid_hght
        out_attributes["geoid_hght"] = np.nanmean(in_lakesp_archive["geoid_hght"])
        if not np.isfinite(out_attributes["geoid_hght"]):
            out_attributes["geoid_hght"] = None
            
        return out_attributes
    
    def build_avg_boundary(self, in_lakesp_archive, in_attributes_avg):
        """
        Build LakeAvg feature boundary from list of LakeSP geometries given in input
        IF at least one full observation exists: the averaged geometry corresponds to the LakeSP geometry whose wse is closest to wse_avg
        ELSE: average geometry = Union of all partial polygons
        
        :param in_lakesp_archive: LakeSP archive related to the current PLD lake
        :type in_lakesp_archive: dict
        :param in_attributes_avg: average attributes of the current PLD lake (except area and storage change)
        :type in_attributes_avg: dict
        
        :return: out_geom = geometry of LakeAvg feature
        :rtype: out_geom = OGRPolygon
        :return: out_attributes = LakeAvg attributes (None for each item which could not be computed)
        :rtype: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Init output dictionary
        out_attributes = dict()
        out_attributes["area_avg"] = None
        out_attributes["area_avg_u"] = None
        
        if in_attributes_avg["npass_full"] == 0:
            
            # 1 - Compute geometry
            out_geom = None
            for geom in in_lakesp_archive["geom"]:
                if out_geom is None:
                    out_geom = geom
                else:
                    out_geom = out_geom.Union(geom)
                    
            # 2 - Compute related attributes
            out_attributes["area_avg"] = my_tools.get_area(out_geom) / 10**6
            out_attributes["area_avg_u"] = my_tools.value_or_none(np.nanmean(in_lakesp_archive["area_tot_u"]))
                    
        else:
            
            # 0 - Retrieve index of the LakeSP feature with the wse closest to wse_avg
            ind_mean = np.argmin(abs(in_lakesp_archive["wse"] - in_attributes_avg["wse_avg"]))
            logger.debug("Pass {} with wse={:.3f} is the closest of mean wse={:.3f}".format(in_lakesp_archive["pass"][ind_mean], in_lakesp_archive["wse"][ind_mean], in_attributes_avg["wse_avg"]))
            
            # 1 - Retrieve and store associated geometry
            out_geom = in_lakesp_archive["geom"][ind_mean]
            
            # 2 - Retrieve and store associated area attributes
            out_attributes["area_avg"] = in_lakesp_archive["area_total"][ind_mean]
            out_attributes["area_avg_u"] = in_lakesp_archive["area_tot_u"][ind_mean]
        
        return out_geom.Clone(), out_attributes
    
    # ------------------------------------
    # Functions specific to PLD attributes
    # ------------------------------------
    
    def add_pld_features_not_observed(self):
        """
        Add PLD lakes which were not observed
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        nb_missing = len(self.obj_lakedb.list_lakeid)
        
        if nb_missing == 0:
            logger.debug("ALL PLD lakes have been observed")
            
        else:
            logger.debug("%d PLD lakes have NOT been observed" % nb_missing)
            
            for cur_lakeid in self.obj_lakedb.list_lakeid:
                
                # 1.1 - Create prior lake object
                obj_plake = lake_db.PriorLake(self.obj_lakedb, cur_lakeid)
                # 1.2 - Format PLD attributes (i.e. convert None to required _FillValue if needed)
                pld_attributes = obj_plake.format_attributes()
                
                # 2 - Add lake_id
                pld_attributes["lake_id"] = cur_lakeid
                
                # 3 - Add prior feature to _Prior layer
                self.content.add_feature(None, pld_attributes)

    # ------------------------------------------------
    # Functions specific to storage change computation
    # ------------------------------------------------
    
    def compute_direct_storage_change(self, in_obj_plake, in_avg_attributes):
        """
        Compute storage change from cycle-averaged WSE and area for the PLD lake
        with DIRECT approach
        
        :param in_obj_plake: PLD lake object
        :type in_obj_plake: lake_db.PriorLake
        :return: in_prior_attributes = cycle-averaged attributes related to PLD lake
        :rtype: in_prior_attributes = dict
        
        :return: out_storage_values = storage change values related to current PLD lake
        :rtype: out_storage_values = dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # 0 - Init output
        out_storage_values = dict()
        out_storage_values["ds1_l_avg"] = None
        out_storage_values["ds1l_avg_u"] = None
        out_storage_values["ds1q_avg"] = None
        out_storage_values["ds1q_avg_u"] = None
        
        # 1 - Build dictionary going in input of storage change function
        list_obs = dict()
        list_obs["pld"] = dict()
        list_obs["pld"]["area"] = in_avg_attributes["area_avg"]
        list_obs["pld"]["area_u"] = in_avg_attributes["area_avg_u"]
        list_obs["pld"]["wse"] = in_avg_attributes["wse_avg"]
        list_obs["pld"]["wse_u"] = in_avg_attributes["wse_avg_u"]
        # Test if extreme values
        for key in ["area", "wse"]:
            if (list_obs["pld"][key] is None) or (not np.isfinite(list_obs["pld"][key])) or (list_obs["pld"][key] < -9e11):
                logger.debug("Lake has {}={} => storage change not computed".format(key, list_obs["pld"][key]))
                list_obs = None
                break
        if list_obs is not None:
            for key in ["area_u", "wse_u"]:
                if (list_obs["pld"][key] is None) or (not np.isfinite(list_obs["pld"][key])) or (list_obs["pld"][key] < -9e11):
                    key2 = key.replace("_u", "")
                    logger.debug("Lake has {}={} => {} set to {}={}".format(key, list_obs["pld"][key], key, key2, list_obs["pld"][key2]))
                    list_obs["pld"][key] = list_obs["pld"][key2]
            
        # 2 - Compute storage change for PLD lake
        out_storage_values["ds1_l_avg"], out_storage_values["ds1l_avg_u"],\
        out_storage_values["ds1q_avg"], out_storage_values["ds1q_avg_u"] = in_obj_plake.run_stocc(list_obs)
        
        return out_storage_values
    
    def compute_incremental_storage_change(self, in_obj_plake, in_avg_attributes):
        """
        Compute storage change from cycle-averaged WSE and area for the PLD lake
        with INCREMENTAL approach
        
        :param in_obj_plake: PLD lake object
        :type in_obj_plake: lake_db.PriorLake
        :return: in_prior_attributes = cycle-averaged attributes related to PLD lake
        :rtype: in_prior_attributes = dict
        
        :return: out_storage_values = storage change values related to current PLD lake
        :rtype: out_storage_values = dict
        """
        
        # 0 - Init output
        out_storage_values = dict()
        out_storage_values["ds2_l_avg"] = my_var.FV_REAL
        out_storage_values["ds2l_avg_u"] = my_var.FV_REAL
        out_storage_values["ds2_q_avg"] = my_var.FV_REAL
        out_storage_values["ds2q_avg_u"] = my_var.FV_REAL
            
        return out_storage_values

    # ----------------------------------------
    
    def write_tmp_file(self, in_filename):
        """
        Write LakeAvg memory layer as a shapefile
        
        :param in_filename: full path of the output file
        :type in_filename: string
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # 1 - Write them in the output shapefile
        my_shp.write_mem_layer_as_shp(self.content.layer, in_filename)
    

#######################################
        
    
def compute_avg_values(in_type, in_archive):
    """
    Return min/med/max values of the archive in input
    
    :param in_type: the selected operation to apply on the archive; =min/med/max
    :type in_type: string
    :param in_archive: archive of observations
    :type in_archive: dict
    
    :return: out_attributes = lake attributes (None for each item which could not be computed)
    :rtype: dict
    """
    logger = logging.getLogger("proc_lakeavg")
    #logger.debug("== compute_avg_values ; type = %s ==" % in_type)
    
    # Init output dictionary
    out_attributes = dict()
    
    # List of attributes retrieved from LakeSP_Prior file
    list_att_sp = ["time", "time_tai", "time_str", \
                   "wse", "wse_u", \
                   "area_total", "area_tot_u", \
                   "ds1_l", "ds1_l_u", "ds1_q", "ds1_q_u", \
                   "ds2_l", "ds2_l_u", "ds2_q", "ds2_q_u", \
                   "partial_f"]
    # List of corresponding attributes in LakeAvg product
    list_att_avg = ["t", "t_tai", "t_str", \
                    "wse", "wse_u", \
                    "area", "are_u", \
                    "ds1_l", "ds1l_u", "ds1_q", "ds1q_u", \
                    "ds2_l", "ds2l_u", "ds2_q", "ds2q_u", \
                    "partf"]
        
    # 1 - Retrieve indice of min/med/max wse
    if in_type == "min":
        ind = np.argmin(in_archive["wse"])
    elif in_type == "med":
        med = np.median(in_archive["wse"])
        ind = np.argmin(abs(in_archive["wse"] - med))
    elif in_type == "max":
        ind = np.argmax(in_archive["wse"])
    else:
        logger.error("in_type = %s => UNKNOWN, should be min/med/max")
        
    # 2 - Compute attributes
    for att_sp, att_avg in zip(list_att_sp, list_att_avg):
        
        # 2.1 - Form key
        if att_avg.endswith("_u"):
            if att_avg.startswith("ds"):
                key = "%sh%s_u" % (att_avg[:-2], in_type)
            else:
                key = "%s_h%s_u" % (att_avg[:-2], in_type)
        else:
            key = "%s_h%s" % (att_avg, in_type)
            
        # 2.2 - Retrieve value
        out_attributes[key] = in_archive[att_sp][ind]
        
    return out_attributes
    

#######################################        
        
    
def write_file(in_list_tmp_lakeavg_shp, in_filename, in_proc_metadata):
    """
    Write the combined LakeAvg shapefile
    
    :param in_list_tmp_lakeavg_shp: list of temporary LakeAvg shapefiles to combine
    :type in_list_tmp_lakeavg_shp: set
    :param in_filename: full path of the output file
    :type in_filename: string
    :param in_proc_metadata: processing metadata
    :type in_proc_metadata: dict
    """
    cfg = service_config_file.get_instance()
    logger = logging.getLogger("proc_lakeavg")
    logger.debug("== write_file ==")
    
    if len(in_list_tmp_lakeavg_shp) == 0:
        logger.debug("No temporary shapefile generated => EMPTY LakeAvg product to generate")
        
        # 1 - Init empty LakeAvg layer
        obj_lakeavg = LakeAvgProduct(None, None, os.path.splitext(os.path.basename(in_filename))[0])
        
        # 2 - Write empty layer to output file
        obj_lakeavg.write_tmp_file(in_filename)
        
    else:

        # 1 - Merge shapefiles
        flag_add_nogeom_feat = True
        if not cfg.getboolean("CONFIG_PARAMS", "ADD_ALL"):
            flag_add_nogeom_feat = False
        my_shp.merge_shp(in_filename, in_list_tmp_lakeavg_shp,
                                                      flag_add_nogeom_feat=flag_add_nogeom_feat)

        # 2 - Delete temporary shapefiles if asked
        if cfg.getboolean("CONFIG_PARAMS", "DEL_TMP_SHP"):
            shp_driver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Driver for shapefiles
            for cur_tmp_file in in_list_tmp_lakeavg_shp:
                logger.debug("Delete temporary shapefile = %s" % os.path.basename(cur_tmp_file))
                shp_driver.DeleteDataSource(cur_tmp_file)
        else:
            logger.debug("Keep temporary sub-basin LakeAvg shapefiles")
    
    # 4 - Write XML metadatafile for shapefile
    logger.debug("Writing associated metadata file = %s.xml" % in_filename)
    tmp_lakeavg = shp_file.LakeAvgProduct(os.path.splitext(os.path.basename(in_filename))[0], 
                                          flag_create=False)
    tmp_lakeavg.update_and_write_metadata("%s.xml" % in_filename, 
                                          in_proc_metadata=in_proc_metadata)
    
