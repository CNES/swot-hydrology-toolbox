# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2019/05/17:version initiale.
# VERSION:2.0.0:DM:#91:2020/07/03:Poursuite industrialisation
# VERSION:3.0.0:DM:#91:2021/03/12:Poursuite industrialisation
# VERSION:3.1.0:DM:#91:2021/05/21:Poursuite industrialisation
# VERSION:3.2.0:DM:#91:2021/10/27:Poursuite industrialisation
# VERSION:4.0.0:DM:#91:2022/05/05:Poursuite industrialisation
# VERSION:5.0.0:DM:#91:2022/12/08:Poursuite industrialisation
# VERSION:5.0.1:DM:#91:2023/02/24:Poursuite industrialisation
# VERSION:5.1.0:DM:#91:2023/03/31:Poursuite industrialisation
# VERSION:5.2.0:DM:#91:2023/05/23:Poursuite industrialisation
# VERSION:5.2.1:DM:#91:2023/07/24:Poursuite industrialisation
# VERSION:5.3.0:DM:#91:2023/10/04:Poursuite industrialisation
# VERSION:5.3.1:DM:#91:2024/01/08:Poursuite industrialisation
# VERSION:6.1.0:DM:#91:2024/09/20:Poursuite industrialisation
# VERSION:6.2.0:DM:#91:2024/11/29:Poursuite industrialisation
# VERSION:6.2.2:DM:#822:2025/06/30:Nouveau patch Lake_TILE_V6.2.1
# VERSION:6.2.3:FA:#839:2025/11/10:[LAKE_TILE] L2_HR_LAKETILE fails on C032P429T305R with message ValueError: list.remove(x): x not in list
# VERSION:6.2.6:DM:#900:2026/07/06: [PGE] [LAKES] Ajout d'une interface niveau moyen des lacs
# VERSION:6.2.6:DM:#900:2026/07/06:SSALTO-565 [PGE] [LAKES] Ajout d'une interface coeff NN pour calcul ice_dyn_f
# VERSION:6.2.7:DM:#906:2026/07/28:SSALTO-872 [PGE][LAKES] Gestion des erreurs/warning GDAL
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: proc_lake.py
    :synopsis: Deals with LakeTile and LakeSP shapefile products
     Created on 2017/02/28

.. moduleauthor: Claire POTTIER (CNES DSO/SI/TR) and Cécile CAZALS (CS)

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import absolute_import, division, print_function, unicode_literals
import json
import logging
from lxml import etree as ET
import numpy as np
import os
from osgeo import ogr
from shapely import geometry, prepared, vectorized

import cnes.common.service_config_file as service_config_file
import cnes.common.service_error as service_error

from cnes.modules.geoloc.scripts.biglake_model import BigLakeModel

import cnes.common.lib.my_hull as my_hull
import cnes.common.lib.my_shp_file as my_shp
import cnes.common.lib.my_timer as my_timer
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_variables as my_var
import cnes.common.lib_lake.lake_db as lake_db
import cnes.common.lib_lake.locnes_filenames as locnes_filenames
import cnes.common.lib_lake.locnes_products_shapefile as shp_file
import cnes.common.lib_lake.proc_pixc_vec as proc_pixc_vec

import udes.pixel_ice as pixel_ice


class LakeProduct(object):
    """
    class LakeProduct
    Manage LakeTile and LakeSP products
    """
    
    def __init__(self, in_product_type, in_obj_pixc, in_obj_pixcvec, in_obj_lake_db, in_obj_mls, in_layer_name, in_swath_side=None):
        """
        Constructor

        :param in_product_type: type of product among "SP"=LakeSP and "TILE"=LakeTile
        :type in_product_type: string
        :param in_obj_pixc: pixel cloud from which to compute lake products
        :type in_obj_pixc: proc_pixc.PixelCloud or proc_pixc_sp.PixelCloudSP object
        :param in_obj_pixcvec: pixel cloud complementary file from which to compute lake products
        :type in_obj_pixcvec: proc_pixc_vec.PixelCloudVec or proc_pixc_vec_sp.PixelCloudVecSP object
        :param in_obj_lake_db: lake database
        :type in_obj_lake_db: lake_db.lakeDb_shp or lake_db.lakeDb_sqlite
        :param in_obj_mls: mean lake surface for large lakes
        :type in_obj_mls: proc_mls.MeanLakeSurface object
        :param in_layer_name: name for lake product layer
        :type in_layer_name: string
        :param in_swath_side: swath side
        :type in_swath_side: string

        Variables of the object:
            - cfg / service_config_file.cfg: instance of LOCNES configuration file
            - product_type / string: type of product among "SP"=LakeSP and "TILE"=LakeTile
            - obj_pixc / proc_pixc.PixelCloud or proc_pixc_sp.PixelCloud: pixel cloud from which to compute lake products
            - obj_pixcvec / proc_pixc_vec.PixelCloudVec or proc_pixc_vec_sp.PixelCloudVec: extra info for pixel cloud
            - obj_lake_db / lake_db.lakeDb_shp or lake_db.lakeDb_sqlite: lake database
            - content_obs / LakeSPShpProduct: container of the lake "obs" product
            - content_prior / LakeSPShpProduct: container of the lake "prior" product
            - content_unassigned / LakeSPShpProduct: container of the lake "unassigned" product
            - lakeid_uniq / set: list of uniq prior identifiers linked to all observed objects
            - compare_stats / dict: store parameters for global comparison between _Obs and _Prior lake products
            - compare_stats_param / list: list of parameters to compare between _Obs and _Prior lake products
        """
        # Get instance of service config file
        self.cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        
        # If ADD_ALL in lake_tile_param.cfg: retrieved as a string => getboolean needed to convert as a boolean
        # If not: init as a boolean => must use get and not getboolean
        try:
            flag_add_all = self.cfg.getboolean("CONFIG_PARAMS", "ADD_ALL")
        except:
            flag_add_all = self.cfg.get("CONFIG_PARAMS", "ADD_ALL")
        
        # 1 - Init variables
        # Product type
        if (in_product_type != "TILE") and (in_product_type != "SP"):
            message = "ERROR = product type is %s ; should be SP or TILE" % in_product_type
            raise service_error.ProcessingError(message, logger)
        else:
            self.product_type = in_product_type
        # Pixel cloud object
        self.obj_pixc = in_obj_pixc
        # Pixel cloud complementary file object
        self.obj_pixcvec = in_obj_pixcvec
        # Lake database object
        self.obj_lake_db = in_obj_lake_db
        # Mean Lake Surface object
        self.obj_mls = in_obj_mls
        # Swath side to be used by proc_pixc_vec.run_hconstr_geoloc
        if in_swath_side is None:
            # Case of LakeTile
            try:
                self.swath_side = in_obj_pixc.pixc_metadata["swath_side"]
            except ValueError:
                logger.error("Global attribute swath_side is not in obj_pixc object.", exc_info=True)
                raise
        else:
            # Case of LakeSP
            self.swath_side = in_swath_side
        
        # 2 - Initialize lake product contents
        if self.product_type == "TILE":
            
            self.content_obs = shp_file.LakeTileObsProduct(in_layer_name)  # Obs-oriented file type
            self.content_prior = shp_file.LakeTilePriorProduct(in_layer_name)  # PLD-oriented file type
            self.content_unassigned = shp_file.LakeTileUnassignedProduct(in_layer_name)  # Obs-oriented file type for unassigned water bodies
            
            # 2.1 - Find continents and basins associated to tile
            continent_id_list, continent_code_list, basin_code_list = self.obj_lake_db.link_poly_to_continent_and_basin(self.obj_pixc.tile_poly)
            self.obj_pixc.pixc_metadata["continent_id"] = continent_id_list
            self.obj_pixc.pixc_metadata["continent_code"] = continent_code_list
            self.obj_pixc.pixc_metadata["basin_code"] = basin_code_list
            self.obj_pixcvec.pixcvec_metadata["continent_id"] = continent_id_list
            self.obj_pixcvec.pixcvec_metadata["continent_code"] = continent_code_list
            
        else:
            self.content_obs = shp_file.LakeSPObsProduct(in_layer_name)  # Obs-oriented file type
            self.content_prior = shp_file.LakeSPPriorProduct(in_layer_name)  # PLD-oriented file type
            # Init with all prior lakes over the granule if requested
            if flag_add_all:
                self.obj_lake_db.init_prior_layer(self.content_prior, in_date=self.obj_pixc.date)
            self.content_unassigned = shp_file.LakeSPUnassignedProduct(in_layer_name)  # Obs-oriented file type for unassigned water bodies
            
        # 3 - Other variables
        self.lakeid_uniq = set()  # List of uniq prior identifiers linked to all observed objects
        # Dictionnary to store parameters for global comparison between _Obs and _Prior lake products
        self.compare_stats = {}
        self.compare_stats_params = ["area_total", "area_detct"]
        for key in ["obs", "prior"]:
            self.compare_stats[key] = {}
            for param in self.compare_stats_params:
                self.compare_stats[key][param] = 0.0
        
    def free_memory(self):
        """
        Destroy memory layers
        """
        self.content_obs.free()  # Obs-oriented layer
        self.content_prior.free()  # PLD-oriented layer
        self.content_unassigned.free()  # Obs-oriented file type for unassigned water bodies

    # ----------------------------------------

    def manage_param_geoloc(self):
        """
        Check geoloc parameters
        Log info
        """
        logger = logging.getLogger(self.__class__.__name__)

        # Retrieve needed configuration parameters
        hull_method = self.cfg.getfloat("CONFIG_PARAMS", "HULL_METHOD")
        # If IMP_GEOLOC in lake_tile_param.cfg: retrieved as a string => getboolean needed to convert as a boolean
        # If not: init as a boolean => must use get and not getboolean
        try:
            flag_geoloc = self.cfg.getboolean('CONFIG_PARAMS', 'IMP_GEOLOC')
        except:
            flag_geoloc = self.cfg.get('CONFIG_PARAMS', 'IMP_GEOLOC')
       
        # Print geoloc
        if flag_geoloc:
            logger.debug("Improve geoloc = YES")
        else:
            logger.debug("Improve geoloc = NO")
        # Print hull_method
        # TODO : les égalités de réel ne sont pas autorisé il faut changer ces tests
        if hull_method == 0:
            logger.debug("Use of CONVEX HULL for lake boundaries")
        elif hull_method == 1.0:
            logger.debug("Use of CONCAVE HULL based on Delaunay triangulation with CGAL library for lake boundaries")
        elif hull_method == 1.1:
            logger.debug("Use of CONCAVE HULL based on Delaunay triangulation")
        elif hull_method == 1.2:
            logger.debug("Use of CONCAVE HULL based on Delaunay triangulation with varying alpha param for lake boundaries")
        elif hull_method == 2.0:
            logger.debug("Hull computation method : radar vectorization with reordered lon lat (default)")
        elif hull_method == 2.1:
            logger.debug("Hull computation method : radar vectorization check segments building polygon")
        else:
            message = "HULL_METHOD values unkown (%d); should be 0(convex) 1(concave) 2(radar vect)" % hull_method
            raise service_error.ProcessingError(message, logger)

    def compute_lake_features(self, in_list_labels, labels_biglakes=[], lakeid_biglakes=[]):
        """
        Computes lake features for pixels for which label is in in_list_labels.

         - NB: This processing is limited to water bodies being of a minimum size defined by MIN_SIZE.
         - NB2: Improved geolocation is computed for all water regions, not only those > MIN_SIZE

        :param in_list_labels: list of labels to process
        :type in_list_labels: 1D-array of int
        :param labels_biglakes: list of labels related to big lakes which will have a specific processing ("SP" processing only)
        :type labels_biglakes: list of int
        :param lakeid_biglakes: list of PLD lake_id related to big lakes which will have a specific processing ("SP" processing only)
        :type lakeid_biglakes: list of int
        """
        logger = logging.getLogger(self.__class__.__name__)

        # Check geoloc parameter
        self.manage_param_geoloc()

        # Get min_size parameter
        min_size = self.cfg.getfloat('CONFIG_PARAMS', 'MIN_SIZE')
        # Print min_size 
        logger.debug("Minimum size for lakes = %0.2f km2" % min_size)
        # Get maximum tolerated difference between improved height and initial height
        max_diff_height = self.cfg.getfloat("CONFIG_PARAMS", "MAX_DIFF_HEIGHT_GEOLOC")
  
        ########################################################
        # Compute geometry and attributes of observed features #
        # except those related to PLD (ie lake_id, overlap,    #
        # p_ attributes and storage change)                    #
        ########################################################
        
        # Init variables
        cpt_too_small = 0  # Counter of too small objects
        cpt_obj = 0  # Counter of processed objects
        cpt_id = 1  # Incremental counter for obs_id
        
        timer = my_timer.Timer()
        timer.start()
        
        for indl, label in enumerate(in_list_labels):  # Loop on inside tile objects
            logger.debug("")
            
            # Flag to identify big lakes
            flag_biglake = False
            if label in labels_biglakes:
                flag_biglake = True

            # ==================================================
            # 1 - Get pixels indices associated to current label
            # ==================================================
            pix_index = np.where(self.obj_pixc.labels == label)[0]
            obj_nb_pix = pix_index.size
            if obj_nb_pix == 0:
                logger.warning("[STRANGE...] label %s corresponds to 0 pixel..." % label)
                continue
            
            mean_long = np.mean(self.obj_pixc.longitude[pix_index])
            mean_lat = np.mean(self.obj_pixc.latitude[pix_index])
            logger.debug(\
			"===== compute_product %d over %d / label = %d / nb pixels = %d / long = %.6f - lat = %.6f =====" \
                        % (indl+1, len(in_list_labels), label, obj_nb_pix, mean_long, mean_lat))
            
            # ==========================================================================
            # 2 - Remove pixels outside the specified cross-track distance range, if any
            #     and compute partial flag
            # ==========================================================================
            
            # 2.1 - Test if feature is partial in the across-track direction
            pix_index, flag_partial = self.obj_pixc.remove_pixels_outside_window(pix_index)
            
            # 2.2 - Test if feature is partial in the along-track direction
            flag_partial2 = self.obj_pixc.compute_azimuth_partial_flag(pix_index)
            
            # 2.3 - Merge both
            if flag_partial2 == 1:
                flag_partial = 1
            
            # ==========================================================
            # 3 - Subset of self.obj_pixc.classif_dict for current label
            # ==========================================================
            classif = dict()
            for key in self.obj_pixc.classif_dict.keys():
                classif[key] = self.obj_pixc.classif_dict[key][pix_index]
            # Add mandatory pixels if there are not enough pixels to compute wse
            nb_pixc_wse = np.sum(classif["4wse"])
            if nb_pixc_wse < self.cfg.getint("CONFIG_PARAMS", "MIN_NB_PIXC"):
                logger.debug("WARNING NOT ENOUGH valid 4wse pixels to compute wse (%d instead of %d) => try with mandatory pixels" % (nb_pixc_wse, self.cfg.getint("CONFIG_PARAMS", "MIN_NB_PIXC")))
                classif["4wse"] = classif["mandatory"] & classif["not_specular_ringing"]
                nb_pixc_wse = np.sum(classif["4wse"])
                if nb_pixc_wse < self.cfg.getint("CONFIG_PARAMS", "MIN_NB_PIXC"):
                    logger.debug("WARNING Still NOT ENOUGH valid mandatory pixels to compute wse (%d) => water feature removed from product" % nb_pixc_wse)
                    continue
                logger.debug("NOW %s valid mandatory pixels to compute wse => continue processing of water feature" % nb_pixc_wse)

            # ==========================================
            # 4 - Compute object median and mean heights
            # ==========================================
            
            # Compute median height using open water pixels without specular ringing (weighting using uncertainties)
            # Use of external method common with RiverObs to compute median/mean/std
            obj_median_height = self.obj_pixc.compute_height(pix_index[classif["4wse"]], method="median")
            obj_mean_height = self.obj_pixc.compute_height(pix_index[classif["4wse"]])
            if not np.isfinite(obj_median_height):
                logger.debug("WARNING median height for the water feature = NaN => height-constrained geolocation cannot be computed => water feature removed from product")
                continue

            # ============================
            # 5 - Compute object raw sizes
            # ============================
            
            # 5.1 - Total area (=water + dark water)
            obj_area_total = np.nansum(self.obj_pixc.inundated_area[pix_index]) 
            # In m2
            obj_area_total /= 10**6  # Conversion in km2
            
            # 5.2 - Mandatory water area 
            obj_area_mandatory = np.nansum(self.obj_pixc.inundated_area[pix_index[classif["mandatory"]]])
            # In m2
            obj_area_mandatory /= 10**6  # Conversion in km2

            # ====================================================
            # 6 - Compute lake feature 
            #     if interior water area of object is large enough 
            #     and all pixels are not aligned
            # ====================================================
            is_aligned = my_tools.are_pixels_aligned(self.obj_pixc.range_index[pix_index], self.obj_pixc.azimuth_index[pix_index])

            if obj_area_mandatory >= min_size and not is_aligned:
                logger.debug("> Water region PROCESSED: area of mandatory pixels (raw total area) = %.2f km2 (%.2f km2)"\
                             % (obj_area_mandatory, obj_area_total))
            
                # 6.1 - Compute height-constrained geolocation
                imp_lon, imp_lat, imp_height = self.compute_hconstr_geoloc(pix_index, classif, obj_median_height, obj_mean_height, obj_area_total)
                
                # 6.2 - Remove pixels with high diff between improved height and initial height
                diff_height = abs(self.obj_pixc.height[pix_index] - imp_height)
                idx_nok = np.where(diff_height > max_diff_height)
                if len(idx_nok) > 0 :
                    pix_index = np.delete(pix_index, idx_nok)
                    obj_nb_pix = pix_index.size
                    imp_lon = np.delete(imp_lon, idx_nok)
                    imp_lat = np.delete(imp_lat, idx_nok)
                    imp_height = np.delete(imp_height, idx_nok)
                    classif = dict()
                    for key in self.obj_pixc.classif_dict.keys():
                        classif[key] = self.obj_pixc.classif_dict[key][pix_index]
                        
                if not flag_biglake:
                    # 6.3 - Update PIXCVec with height-constrained geolocation infos
                    self.update_pixcvec_with_hconstr_geoloc(pix_index, imp_lon, imp_lat, imp_height)
                    
                # 6.4 - Select valid PIXC wrt value of longitude, latitude or height
                #       which must be finite
                valid_flag, not_nan_index = select_valid_pixels(obj_nb_pix, imp_lon, imp_lat, imp_height)
                # No need to continue if there is no valid PIXC
                if valid_flag == 0:  
                    continue

                # 6.5 - Prepare needed arrays of valid PIXC
                if valid_flag == 2:
                    valid_pix_index = pix_index
                    valid_classif = classif
                    valid_lon = imp_lon
                    valid_lat = imp_lat
                else:
                    valid_pix_index = pix_index[not_nan_index]
                    valid_classif = dict()
                    for key in self.obj_pixc.classif_dict.keys():
                        valid_classif[key] = self.obj_pixc.classif_dict[key][valid_pix_index]
                    valid_lon = imp_lon[not_nan_index]
                    valid_lat = imp_lat[not_nan_index]
                        
                if flag_biglake:
                        
                    # 6.6 - Add patch to _Prior layer
                    self.add_prior_patch(label, cpt_id, lakeid_biglakes,
                                         valid_pix_index, valid_classif, 
                                         valid_lon, valid_lat)
                
                else:
                    
                    # 6.6 - Add observed feature
                    obs_id, pixcvec_lakeid = self.add_obs_feature(label, cpt_id, 
                                                                  valid_pix_index, valid_classif, 
                                                                  valid_lon, valid_lat, 
                                                                  flag_partial=flag_partial)
                    
                    if obs_id is not None:
                        # 6.7 - Update PIXCVec obs_id and lake_id attributes
                        self.update_pixcvec_with_ids(valid_pix_index, obs_id, pixcvec_lakeid)
                        cpt_id += 1
                        # 6.8 - Update PIXCVec with ice_dyn_f
                        self.update_pixcvec_with_ice_dyn_f(valid_pix_index)

                # 6.9 - Increase counter of processed objects
                cpt_obj += 1

            else:
                logger.debug("> Water region NOT PROCESSED: interior water area = %.3f km2 < %.2f km2 (total area = %.3f km2 with %d pixels)" 
                                % (obj_area_mandatory, min_size, obj_area_total, obj_nb_pix))
                cpt_too_small += 1  # Increase counter of too small objects

        logger.debug("")
        logger.debug(timer.stop(message="Computation time for observed features: "))
        logger.debug("> %d objects processed" % cpt_obj)
        logger.debug("> %d objects not processed because too small" % cpt_too_small)
        logger.debug("")

        ##################################################################
        # Compute _p attributes and storage change for observed features #
        # + compute geometry and attributes of PLD features              #
        ##################################################################
        
        timer.start()
        
        nb_prior = len(self.lakeid_uniq)
        
        if nb_prior == 0:
            logger.debug("NO observed feature linked to a PLD lake")
            
        else:
            # Save _Fillvalue in NetCDF and shapefile for ice clim flag attribute
            fv_shp_icef = my_var.FV_SHP[my_var.PLD_FIELD_LAKE_TYPE["ice_clim_flag"]]
            fv_nc_icef = my_var.FV_NETCDF[my_var.PLD_FIELD_LAKE_TYPE["ice_clim_flag"]]
            
            # Deal with PLD lakes observed by SWOT
            logger.debug("==================")
            logger.debug("Compute PLD features")
            logger.debug("==================")
            logger.debug("")
            logger.debug("%d PLD lakes linked to observed lakes" % nb_prior)

            for i, cur_lakeid in enumerate(self.lakeid_uniq):
                logger.debug("")
                logger.debug("===== Deal with PLD lake %s (%d/%d) =====" % (cur_lakeid, i+1, len(self.lakeid_uniq)))
        
                # 7.1.1 - Create prior lake object
                obj_plake = lake_db.PriorLake(self.obj_lake_db, cur_lakeid)
                # 7.1.2 - Format PLD attributes (i.e. convert None to required _FillValue if needed)
                pld_attributes = obj_plake.format_attributes()
                
                # 7.2 - Update p_ attributes of observed features strongly connected to this PLD lake
                # 7.2.1 - Select them
                self.content_obs.layer.SetAttributeFilter("lake_id LIKE '{}%'".format(cur_lakeid))
                nb_obslake = self.content_obs.layer.GetFeatureCount()
                logger.debug("{} observed lake(s) are strongly connected to this PLD lake".format(nb_obslake))
                # 7.2.2 - Set p_ attributes from PLD infos to all observed lakes having this PLD lake as main overlap
                if nb_obslake > 0:
                    self.set_pld_attributes(pld_attributes)
                # 7.2.3 - Reinit layer attribute filter
                self.content_obs.layer.SetAttributeFilter(None)
                
                # 7.3 - Select all observed lakes overlapping the PLD lake
                self.content_obs.layer.SetAttributeFilter("lake_id LIKE '%{}%'".format(cur_lakeid))
                nb_obslake = self.content_obs.layer.GetFeatureCount()
                logger.debug("{} observed lake(s) are connected to this PLD lake".format(nb_obslake))
                
                if nb_obslake == 0:
                    logger.warning("[STRANGE...] PLD lake listed in _Obs file corresponds to NO obs lake...")
                    
                    if not self.obj_lake_db.flag_prior_init:
                        # 8.4 - Add prior feature to _Prior layer
                        self.content_prior.add_feature(None, pld_attributes)
                    
                else:
        
                    # 7.4 - Retrieve PIXC indices corresponding to prior feature
                    # NB: use of selected_index to convert PIXCVec indices to PIXC indices reference
                    if self.product_type == "SP":
                        pixc_index = np.where(self.obj_pixcvec.lake_id == obj_plake.lake_id.encode())[0]
                    else:
                        pixc_index = np.where(self.obj_pixcvec.lake_id[self.obj_pixc.selected_index] == obj_plake.lake_id.encode())[0]
                        
                    if len(pixc_index) == 0:
                        logger.warning("WARNING - No pixel is actually related to the PLD lake")
                        
                        if not self.obj_lake_db.flag_prior_init:
                            # 7.5 - Add prior feature to _Prior layer
                            self.content_prior.add_feature(None, pld_attributes)
                            
                    else:
                        
                        # 7.6 - Subset of self.obj_pixc.classif_dict for prior feature
                        classif = dict()
                        for key in self.obj_pixc.classif_dict.keys():
                            classif[key] = self.obj_pixc.classif_dict[key][pixc_index]
                            
                        # 7.7 - Add prior feature
                        self.add_prior_feature(obj_plake, pixc_index, classif, 
                                               pld_attributes)
                        
                        # 7.8 - Set PIXCVec used_for_lake_wse attribute
                        self.update_pixcvec_with_used_for_lake_wse(pixc_index, classif["4wse"])
                        
                        # 7.9 - Update PIXCVec ice_clim_f attribute
                        if pld_attributes["ice_clim_f"] == fv_shp_icef:
                            # Need to convert _FillValue from shapefile to NetCDF format
                            self.update_pixcvec_with_ice_clim_f(pixc_index, fv_nc_icef)
                        else:
                            self.update_pixcvec_with_ice_clim_f(pixc_index, pld_attributes["ice_clim_f"])
                
                # 7.9 - Reinit layer attribute filter
                self.content_obs.layer.SetAttributeFilter(None)
                
        logger.debug(timer.stop(message="Computation time for prior features: "))

    # ----------------------------------------
    # Fonctions specific to height-constrained 
    # geolocation process 
    # ----------------------------------------
    
    def compute_hconstr_geoloc(self, in_pix_index, in_classif, in_median_height, in_mean_height, in_area_total):
        """
        Prepare and run the height-constrained geolocation algorithm over the feature
            defined by the PIXC of indices definied in input.
        
        :param in_pix_index: indices of PIXC defining the feature
        :type in_pix_index: 1D-array of int
        :param in_classif: dictionary listing gatherings of pixels (water is used here)
        :type in_classif: dict of 1D-array of boolean
        :param in_median_height: median height of input PIXC
        :type in_median_height: float
        :param in_mean_height: mean height of input PIXC
        :type in_mean_height: float
        :param in_area_total: total area of input PIXC
        :type in_area_total: float
        
        :return: out_lon = improved longitudes
        :rtype: out_lon = 1D-array of float
        :return: out_lat = improved latitudes
        :rtype: out_lat = 1D-array of float
        :return: out_height = improved heights
        :rtype: out_height = 1D-array of float
        """
        
        # If IMP_GEOLOC in lake_tile_param.cfg: retrieved as a string => getboolean needed to convert as a boolean
        # If not: init as a boolean => must use get and not getboolean
        try:
            flag_geoloc = self.cfg.getboolean('CONFIG_PARAMS', 'IMP_GEOLOC')
        except:
            flag_geoloc = self.cfg.get('CONFIG_PARAMS', 'IMP_GEOLOC')
        
        if flag_geoloc:

            # 2a - Fit lake height model depending on lake size
            height_model = self.compute_height_model(in_pix_index, in_classif, in_median_height, in_mean_height, in_area_total)

            # 2b - Compute imp geolocation 
            out_lon, out_lat, out_height, p_final = proc_pixc_vec.run_hconstr_geoloc(self.product_type, self.obj_pixc, in_pix_index, height_model, self.swath_side)
            
            # 2c - Compute flattened interferogram in order to compute coherence during uncertainties estimation
            self.obj_pixc.compute_interferogram_flattened(in_pix_index, p_final)

        else:
            out_lon = self.obj_pixc.longitude[in_pix_index]
            out_lat = self.obj_pixc.latitude[in_pix_index]
            out_height = self.obj_pixc.height[in_pix_index]
                
        return out_lon, out_lat, out_height
    
    def compute_height_model(self, in_pix_index, in_classif, in_median_height, in_mean_height, in_area_total):
        """
        Compute the height model related to the feature defined by the PIXC
        of indices defined in input. This height model depends on the feature size.
        
        :param in_pix_index: indices of PIXC defining the feature
        :type in_pix_index: 1D-array of int
        :param in_classif: dictionary listing gatherings of pixels (water is used here)
        :type in_classif: dict of 1D-array of boolean
        :param in_median_height: median height of input PIXC
        :type in_median_height: float
        :param in_mean_height: mean height of input PIXC
        :type in_mean_height: float
        :param in_area_total: total area (ie water + dark) of input PIXC
        :type in_area_total: float
        
        :return: out_height_model = height model specific to the feature 
        :rtype: out_height_model = 1D-array of float
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Retrieve needed config parameters
        biglake_min_size = self.cfg.getfloat("CONFIG_PARAMS", "BIGLAKE_MIN_SIZE")
        biglake_model = self.cfg.get("CONFIG_PARAMS", "BIGLAKE_MODEL")
        biglake_grid_spacing = self.cfg.getint("CONFIG_PARAMS", "BIGLAKE_GRID_SPACING")
        biglake_grid_res = self.cfg.getint("CONFIG_PARAMS", "BIGLAKE_GRID_RES")
        
        # Compute height model depending on the feature size
        if (biglake_model != 'no') and (in_area_total >= biglake_min_size):

            biglakemodel = BigLakeModel(biglake_model)

            logger.debug("Using {} biglake model for improved geolocation (lake total area = {} km2)".format(biglake_model, in_area_total))

            if biglake_model == 'grid':
                out_height_model = biglakemodel.fit_biglake_model(self.obj_pixc,
                                                                  in_pix_index,
                                                                  grid_spacing=biglake_grid_spacing,
                                                                  grid_resolution=biglake_grid_res)

            elif biglake_model == 'polynomial':                                 
                out_height_model = biglakemodel.fit_biglake_model_polyfit(self.obj_pixc, in_pix_index, in_classif)
                                                 
            else:
                logger.debug("No height model defined, assume Median Height model")
                out_height_model = np.full(self.obj_pixc.height[in_pix_index].shape, in_median_height)

        else:
            logger.debug("Using lake median height = {} m for improved geolocation (lake total area = {} km2, lake average height = {} m)".format(in_median_height, in_area_total, in_mean_height))
            out_height_model = np.full(self.obj_pixc.height[in_pix_index].shape, in_median_height)
            
        # Return
        return out_height_model

    # ----------------------------------------
    # Functions specific to lake feature (i.e.
    # geometry + attributes) computation
    # ----------------------------------------
    
    def add_obs_feature(self, in_label, in_number_for_id, in_pixc_index, in_classif, in_lon, in_lat, flag_partial=0):
        """
        Process valid PIXC related to current feature
        to build the observed feature boundary and associated attributes
        and add the observed feature to the obs-oriented layer
        
        :param in_label: label of the feature in the scene
        :type in_label: int
        :param in_number_for_id: number for obs identifier
        :type in_number_for_id: int
        :param in_pixc_index: indices of the PIXC related to the feature
        :type in_pixc_index: 1D-array of int
        :param in_classif: dictionary listing pixels selected or not for different processes (hull/wse/area computation)
        :type in_classif: dict of 1D-array of boolean
        :param in_lon: improved longitudes vector for PIXC of the feature
        :type in_lon: 1D-array of float
        :param in_lat: improved latitudes vector for PIXC of the feature
        :type in_lat: 1D-array of float
        :param flag_partial: =1 if water body is partially observed; =0 otherwise (default)
        :type flag_partial: boolean
        
        :return: out_obs_id = obs_id identifier of the feature
        :rtype: out_obs_id = string
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Deal with object number = {}".format(in_label))
        
        # 1 - Build feature boundary
        feature_geom = self.build_obs_boundary(in_pixc_index[in_classif["4hull"]], in_lon[in_classif["4hull"]], in_lat[in_classif["4hull"]])
        
        # 2.1 - Compute common attributes
        feature_attributes = self.compute_common_attributes(feature_geom, in_pixc_index, in_classif=in_classif, flag_partial=flag_partial)
        # 2.2 - Compute ice_dyn_f attribute
        feature_attributes["ice_dyn_f"] = pixel_ice.compute_ice_dyn_f_4lake(self.obj_pixc.ice_dyn_f[in_pixc_index], self.obj_pixc.sig0[in_pixc_index])
        
        # 3 - Link to PLD
        lake_id, overlap, out_pixcvec_lakeid = self.link_obs_geom_to_pld_geom(feature_geom, in_lon, in_lat)
        feature_attributes["lake_id"] = lake_id
        if lake_id != "no_data":
            self.lakeid_uniq = self.lakeid_uniq.union(set(lake_id.split(";")))
        feature_attributes["overlap"] = overlap
        feature_attributes["n_overlap"] = len(overlap.split(";"))
        feature_attributes["reach_id"] = None
        
        # 4.1 - Compute 3 digits of basin code (CBB)
        if (self.product_type == "TILE") and (int(in_label) in self.obj_pixc.labels_coastal_zones):
            lake_basin_code = my_var.BASIN_4_OCEAN
        elif lake_id != "no_data":
            lake_basin_code = lake_id[0:3]
        else:
            lake_basin_code = str(self.obj_lake_db.link_poly_to_basin(feature_geom)[0]).rjust(3, str('0'))
        # 4.2 - Form obs_id
        if self.product_type == "TILE":  # "TILE" case: add current tile number and swath info
            tile_ref = str(self.obj_pixc.pixc_metadata["tile_number"]).rjust(3, str('0')) + str(self.obj_pixc.pixc_metadata["swath_side"])
        else:  # "SP" case: add main tile info
            tile_ref = self.obj_pixc.get_majority_pixels_tile_ref(in_pixc_index, int(in_label))
        # 4.3 - Compute obs number
        nb_digits = self.cfg.getint('CONFIG_PARAMS', 'NB_DIGITS')  # Retrieve nb_digits from config file
        if self.product_type == "TILE":  # "TILE" case: only add label
            obs_number = str(in_number_for_id).rjust(nb_digits, str('0'))
        else:  # "SP" case: add main tile info
            obs_number = str(self.obj_pixc.get_lake_tile_label(in_label)).rjust(nb_digits, str('0'))
        feature_attributes["obs_id"] = "%s%s%s" % (lake_basin_code, tile_ref, obs_number)
        
        # 5 - Add feature to layer, if it exists
        logger.debug("obs_id = %s / lake_id = %s" % (str(feature_attributes["obs_id"]), lake_id))
        if feature_geom is not None:
            if lake_id == "no_data":
                logger.debug("Feature added to _Unassigned layer")
                self.content_unassigned.add_feature(feature_geom, feature_attributes)
            else:
                logger.debug("Feature added to _Obs layer")
                self.content_obs.add_feature(feature_geom, feature_attributes)
            out_obs_id = feature_attributes["obs_id"]
        else:
            logger.warning("Feature not added to layer because geometry is None")
            out_obs_id = None
            
        return out_obs_id, out_pixcvec_lakeid
    
    def build_obs_boundary(self, in_pixc_index, in_lon, in_lat):
        """
        Build observed feature boundary from PIXC defined by the input coordinates
        
        :param in_pixc_index: indices of the PIXC selected to compute the boundary of the feature (classif "4hull")
        :type in_pixc_index: 1D-array of int
        :param in_lon: improved longitudes vector for these PIXC
        :type in_lon: 1D-array of float
        :param in_lat: improved latitudes vector for these PIXC
        :type in_lat: 1D-array of float
        
        :return: out_geom = boundary of the feature
        :rtype: out_geom = OGRPolygon
        """
        
        if self.product_type == 'SP':
            out_geom = my_hull.compute_lake_boundaries(in_lon,
                                                       in_lat,
                                                       self.obj_pixc.get_range_of_lake(in_pixc_index),
                                                       self.obj_pixc.azimuth_index_in_full_pass[in_pixc_index],
                                                       self.obj_pixc.nb_pix_range, self.obj_pixc.swath_side, self.obj_pixc.pass_num)
            
        else:
            out_geom = my_hull.compute_lake_boundaries(in_lon,
                                                       in_lat,
                                                       self.obj_pixc.range_index[in_pixc_index],
                                                       self.obj_pixc.azimuth_index[in_pixc_index],
                                                       self.obj_pixc.nb_pix_range, self.obj_pixc.pixc_metadata["swath_side"],
                                                       self.obj_pixc.pixc_metadata["pass_number"])


        return out_geom
    
    def add_prior_feature(self, in_obj_plake, in_pixc_index, in_classif, in_pld_infos):
        """
        Process valid PIXC related to input PLD lake
        to build the prior feature boundary and associated attributes
        compute storage change
        and add the prior feature to the PLD-oriented layer
        
        :param in_obj_plake: currently processed PLD lake
        :type in_obj_plake: lake_db.PriorLake
        :param in_pixc_index: indices of the PIXC related to the feature
        :type in_pixc_index: 1D-array of int
        :param in_classif: dictionary listing pixels selected or not for different processes (hull/wse/area computation)
        :type in_classif: dict of 1D-array of boolean
        :param in_pld_infos: values of available p_ attributes
        :type in_pld_infos: dict
        :param flag_partial: =1 if water body is partially observed; =0 otherwise (default)
        :type flag_partial: boolean
        """
        logger = logging.getLogger(self.__class__.__name__)
        #logger.debug("Deal with PLD lake = {}".format(in_obj_plake.lake_id))
        logger.debug("- start -")
        
        # 1 - Compute observed geometry and common attributes of PLD feature
        prior_geom, prior_attributes = self.form_prior_feature(in_obj_plake, in_pixc_index, in_classif)
        
        if in_obj_plake.ok_to_compute_stocc and (prior_attributes["partial_f"] == 0):
            
            # 2 - Compute direct storage change for this PLD lake and all observed lakes overlapping it
            direct_storage_change_values = self.compute_direct_storage_change(in_obj_plake, in_pixc_index, in_classif, prior_attributes)
            
            # 3 - Compute incremental storage change for this PLD lake and all observed lakes overlapping it
            incremental_storage_change_values = self.compute_incremental_storage_change(in_obj_plake, in_pixc_index, prior_attributes)
            
            # 4 - Add/update feature to/of _Prior layer
            if self.obj_lake_db.flag_prior_init:
                self.content_prior.update_feature(in_obj_plake.lake_id, prior_geom, 
                                                  {**prior_attributes, \
                                                   **direct_storage_change_values, \
                                                   **incremental_storage_change_values})
            else:
                self.content_prior.add_feature(prior_geom, {**prior_attributes, \
                                                            **in_pld_infos, \
                                                            **direct_storage_change_values, \
                                                            **incremental_storage_change_values})
        else:
            
            msg = "Not able to compute storage change"
            if not in_obj_plake.ok_to_compute_stocc:
                msg += "; parameters not available in lake DB"
            if prior_attributes["partial_f"] == 1:
                msg += "; lake is partially observed"
            logger.debug(msg)
            
            # 4 - Add/update feature to/of _Prior layer
            if self.obj_lake_db.flag_prior_init:
                self.content_prior.update_feature(in_obj_plake.lake_id, prior_geom, prior_attributes)
            else:
                self.content_prior.add_feature(prior_geom, {**prior_attributes, **in_pld_infos})
    
    def form_prior_feature(self, in_obj_plake, in_pixc_index, in_classif):
        """
        Create and initialize prior feature 
        
        :param in_obj_plake: PLD lake object
        :type in_obj_plake: lake_db.PriorLake
        :param in_pixc_index: indices of PIXC related to PLD lake
        :type in_pixc_index: Numpy 1D-array of int
        :param in_classif: dictionary listing pixels selected or not for different processes (hull/wse/area computation)
        :type in_classif: dict of 1D-array of boolean
        
        :return: out_geom = geometry of prior feature
        :rtype: out_geom = OGRPolygon
        :return: out_attributes = attributes of prior feature
        :rtype: out_attributes = dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Deal with PLD lake = {}".format(in_obj_plake.lake_id))
        
        # 1 - Build feature boundary
        selected_pixc = in_pixc_index[in_classif["4hull"]]
        if self.cfg.getboolean("CONFIG_PARAMS", "DEL_1PIX_IN_AZ_4HULL"):
            selected_pixc = np.array(self.obj_pixc.compute_lakeshore_erosion(in_pixc_index[in_classif["4hull"]]))
        out_geom, list_obs_id, list_overlap = self.build_prior_boundary(in_obj_plake, selected_pixc)
                
        # 2 - Compute partial flag
        nb_obs_features = self.content_obs.layer.GetFeatureCount()
        if nb_obs_features == 1:
            feature = self.content_obs.layer.GetNextFeature()
            flag_partial = int(feature.GetField("partial_f"))
        else:
            list_flags = list()
            for feature in self.content_obs.layer:
                list_flags.append(int(feature.GetField("partial_f")))
            flag_partial = max(list_flags)
        self.content_obs.layer.ResetReading()
        
        # 3.1 - Compute common attributes
        out_attributes = self.compute_common_attributes(out_geom, in_pixc_index, in_classif=in_classif, flag_partial=flag_partial)
        
        # Save for comparison
        for param in self.compare_stats_params:
            if (param in out_attributes.keys()) and (out_attributes[param] is not None):
                self.compare_stats["prior"][param] += out_attributes[param]
        # 3.2 - Compute ice_dyn_f attribute
        if ("ice_start" in in_obj_plake.pld_info.keys()) and ("ice_end" in in_obj_plake.pld_info.keys()):
            out_attributes["ice_dyn_f"] = pixel_ice.compute_ice_dyn_f_4prior(self.obj_pixc.ice_dyn_f[in_pixc_index], 
                                                                              self.obj_pixc.sig0[in_pixc_index],
                                                                              in_obj_plake.pld_info["ice_start"], 
                                                                              in_obj_plake.pld_info["ice_end"],
                                                                              out_attributes["time_str"],
                                                                              time_window=self.cfg.getint("CONFIG_PARAMS", "WINDOW_4_ICE_CLIM"))
        else:
            logger.debug("ice_start and/or ice_stop not in PLD => DO NOT use prior ice_clim_f")
            out_attributes["ice_dyn_f"] = pixel_ice.compute_ice_dyn_f_4lake(self.obj_pixc.ice_dyn_f[in_pixc_index], self.obj_pixc.sig0[in_pixc_index])
        # 3.3 - Compute prior attributes
        prior_attributes = self.compute_prior_attributes(in_obj_plake, in_pixc_index, out_attributes["area_total"], in_classif=in_classif)
        out_attributes = {**out_attributes, **prior_attributes}
        
        # 4 - Add identifiers and overlap values
        out_attributes["lake_id"] = in_obj_plake.lake_id
        if list_obs_id:
            out_attributes["obs_id"] = ';'.join(list_obs_id)
        else:
            out_attributes["obs_id"] = "no_data"
        if list_overlap:
            out_attributes["overlap"] = ';'.join(list_overlap)
            out_attributes["n_overlap"] = len(list_overlap)
        else:
            out_attributes["overlap"] = "no_data"
            out_attributes["n_overlap"] = 0
            
        # 5 - Update qual_f_b bitflag with diff_from_pld_area_suspect bit
        if out_attributes["area_total"] is not None:
            ratio = out_attributes["area_total"] / in_obj_plake.pld_info["ref_area"]
            if (ratio < 1./3.) or (ratio > 3):
                out_attributes["qual_f_b"] += 64
        
        return out_geom, out_attributes
    
    def build_prior_boundary(self, in_obj_plake, in_pixc_index):
        """
        Build prior feature boundary from PIXCVec pixels with the same lake_id
        
        :param in_obj_plake: PLD lake object
        :type in_obj_plake: lake_db.PriorLake
        :param in_pixc_index: list of indices of the PIXC selected to compute the boundary of the feature (classif "4hull")
        :type in_pixc_index: 1D-array of int
        
        :return: out_geom = geometry of prior feature
        :rtype: out_geom = OGRPolygon
        :return: out_list_obs_id = list of observed features intersecting the PLD lake
        :rtype: out_list_obs_id = list
        :return: out_list_overlap = list of fractions of PLD feature covered by each observed lake
        :rtype: out_list_overlap = list
        """
        logger = logging.getLogger(self.__class__.__name__)
        #logger.debug("Deal with PLD lake = {}".format(in_obj_plake.lake_id))
        
        out_geom = None
        out_list_obs_id = []
        out_list_overlap = []
        
        # Case with 1 obs <-> 1 or N PLD lake(s)
        nb_obs_inter = self.content_obs.layer.GetFeatureCount()
        if nb_obs_inter == 1:
            
            # Retrieve associated obs feature
            cur_feature = self.content_obs.layer.GetNextFeature()
            out_list_obs_id.append(cur_feature.GetField("obs_id"))
            logger.debug("obs_id = %s / lake_id = %s" % (out_list_obs_id[0], cur_feature.GetField("lake_id")))
            
            # Case with 1 obs <-> N PLD lakes
            # PLD lake geometry is part of the obs geometry
            logger.debug("Case with 1 obs <-> N PLD lakes => keep part of OBS geometry")

            tmp_pixcvec_index = self.get_pixcvec_index_from_pixc_index(in_pixc_index)
            tmp_geom = self.build_obs_boundary(in_pixc_index,
                                                self.obj_pixcvec.longitude_vectorproc[tmp_pixcvec_index],
                                                self.obj_pixcvec.latitude_vectorproc[tmp_pixcvec_index])

            # Compute overlaping area
            area_pld = my_tools.get_area(in_obj_plake.geom)
            try:
                geom_inter = tmp_geom.Intersection(in_obj_plake.geom)
            except Exception as e:
                logger.warning(f"PLD lake contains an invalid geometry at lake_id : {cur_feature.GetField("lake_id")}, causing the following error : {e}")
                geom_inter = tmp_geom.Intersection(in_obj_plake.geom.Buffer(0))
            if geom_inter is not None:
                area_inter = my_tools.get_area(geom_inter)
                out_list_overlap.append(str(round(area_inter/area_pld*100.)))
                logger.debug("PLD lake area = {} m2 - PLD/obs intersection area = {} m2 - overlap = {}%".format(area_pld, area_inter, \
                             out_list_overlap[0]))
        
        # Case with N obs <-> 1 or N PLD lake(s)
        else:
            
            logger.debug("%d obs features correspond to PLD lake" % nb_obs_inter)
            
            # Init output geometry
            tmp_geom = ogr.Geometry(ogr.wkbMultiPolygon)

            tmp_list_obs_id = []
            tmp_list_overlap = []
            
            for cur_feature in self.content_obs.layer:
                
                # Retrieve corresponding obs_id:
                cur_obs_id = cur_feature.GetField("obs_id")
                tmp_list_obs_id.append(cur_obs_id)
                logger.debug("Deal with obs_id %s " % cur_obs_id)

                # Complex case with N obs <-> N PLD lakes
                logger.debug("Complex case with obs <-> N PLD lakes => keep part of OBS geometry")
            
                in_pixcvec_index = self.get_pixcvec_index_from_pixc_index(in_pixc_index)
                obsid_over_in_lakeid = self.obj_pixcvec.obs_id[in_pixcvec_index]
                tmp_pixc_index = in_pixc_index[np.where(obsid_over_in_lakeid == cur_obs_id.encode())]
                if len(tmp_pixc_index) == 0:
                    logger.debug("WARNING Strange... NO pixel finally associated to obsid %s" % cur_obs_id)
                    continue
                tmp_pixcvec_index = self.get_pixcvec_index_from_pixc_index(tmp_pixc_index)
                obs_poly = self.build_obs_boundary(tmp_pixc_index,
                                                    self.obj_pixcvec.longitude_vectorproc[tmp_pixcvec_index],
                                                    self.obj_pixcvec.latitude_vectorproc[tmp_pixcvec_index])
                    
                # Compute overlaping area
                area_pld = my_tools.get_area(in_obj_plake.geom)
                try:
                    geom_inter = obs_poly.Intersection(in_obj_plake.geom)
                except Exception as e:
                    logger.warning(f"PLD lake contains an invalid geometry at lake_id : {cur_feature.GetField("lake_id")}, causing the following error : {e}")
                    geom_inter = obs_poly.Intersection(in_obj_plake.geom.Buffer(0))

                if geom_inter is not None:
                    area_inter = my_tools.get_area(geom_inter)
                    tmp_overlap = str(round(area_inter/area_pld*100.))
                    tmp_list_overlap.append(tmp_overlap)
                    logger.debug("PLD lake area = {} m2 - PLD/obs intersection area = {} m2 - overlap = {}%".format(area_pld, area_inter, \
                                 tmp_overlap))
                else :
                    logger.warning("PLD lakeid %s geometry does not intersect observed geometry %s " %(str(in_obj_plake.lake_id), \
                                                                                                       str(tmp_list_obs_id)))
                # Add current geometry to output geometry
                tmp_geom = tmp_geom.Union(obs_poly)

            # Sort obs_id and overlap fractions by decreasing area intersection
            sorted_idx = sorted(range(len(tmp_list_overlap)), key=lambda k: tmp_list_overlap[k], reverse=True)
            out_list_obs_id = [tmp_list_obs_id[idx] for idx in sorted_idx]
            out_list_overlap = [tmp_list_overlap[idx] for idx in sorted_idx]

        if tmp_geom is not None:
            out_geom = tmp_geom.Clone()
        
        return out_geom, out_list_obs_id, out_list_overlap
    
    def add_prior_patch(self, in_label, in_number_for_id, lakeid_biglakes, in_pixc_index, in_classif, in_lon, in_lat):
        """
        Specific processing for patches for big lakes in LakeSP processing:
            Convert all _Obs features to a single _Prior feature and add it to the _Prior layer
            
        :param in_label: label of the feature in the scene
        :type in_label: int
        :param in_number_for_id: number for obs identifier
        :type in_number_for_id: int
        :param lakeid_biglakes: list of PLD lake_id related to big lakes which will have a specific processing ("SP" processing only)
        :type lakeid_biglakes: list of int
        :param in_pixc_index: indices of the PIXC related to the feature
        :type in_pixc_index: 1D-array of int
        :param in_classif: dictionary listing pixels selected or not for different processes (hull/wse/area computation)
        :type in_classif: dict of 1D-array of boolean
        :param in_lon: improved longitudes vector for PIXC of the feature
        :type in_lon: 1D-array of float
        :param in_lat: improved latitudes vector for PIXC of the feature
        :type in_lat: 1D-array of float
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Dealing with a PATCH for big lakes")
        
        # 1 - Build feature boundary
        feature_geom = self.build_obs_boundary(in_pixc_index[in_classif["4hull"]], in_lon[in_classif["4hull"]], in_lat[in_classif["4hull"]])
        if feature_geom is None:
            logger.warning("Feature not added to _Prior layer because geometry is None")
            return
        
        # 2 - Attributes of feature
        feature_attributes = dict()
        
        # 2.1 - Link to PLD
        lake_id, overlap, out_pixcvec_lakeid = self.link_obs_geom_to_pld_geom(feature_geom, in_lon, in_lat)
        if lake_id not in lakeid_biglakes:
            return
        feature_attributes["lake_id"] = lake_id
        
        # 2.2 - Compute 3 digits of basin code (CBB)
        if (self.product_type == "TILE") and (int(in_label) in self.obj_pixc.labels_coastal_zones):
            lake_basin_code = my_var.BASIN_4_OCEAN
        elif lake_id != "no_data":
            lake_basin_code = lake_id[0:3]
        else:
            lake_basin_code = str(self.obj_lake_db.link_poly_to_basin(feature_geom)[0]).rjust(3, str('0'))
        # 4.2 - Form obs_id
        if self.product_type == "TILE":  # "TILE" case: add current tile number and swath info
            tile_ref = str(self.obj_pixc.pixc_metadata["tile_number"]).rjust(3, str('0')) + str(self.obj_pixc.pixc_metadata["swath_side"])
        else:  # "SP" case: add main tile info
            tile_ref = self.obj_pixc.get_majority_pixels_tile_ref(in_pixc_index, int(in_label))
        # 4.3 - Compute obs number
        nb_digits = self.cfg.getint('CONFIG_PARAMS', 'NB_DIGITS')  # Retrieve nb_digits from config file
        if self.product_type == "TILE":  # "TILE" case: only add label
            obs_number = str(in_number_for_id).rjust(nb_digits, str('0'))
        else:  # "SP" case: add main tile info
            obs_number = str(self.obj_pixc.get_lake_tile_label(in_label)).rjust(nb_digits, str('0'))
        feature_attributes["obs_id"] = "%s%s%s" % (lake_basin_code, tile_ref, obs_number)
        
        # 2.3 - Set time_str to indicate that the feature is a patch for big lake in _Prior
        feature_attributes["time_str"] = "patch"
        
        # 3 - Add feature to layer, if it exists
        logger.debug("lake_id = %s / obs_id = %s" % (lake_id, str(feature_attributes["obs_id"])))
        self.content_prior.add_feature(feature_geom, feature_attributes)
        logger.debug("Patch added to _Prior layer")
                
    def compute_common_attributes(self, in_geom, in_pixc_index, in_classif=None, flag_partial=0):
        """
        Computes common attributes from PIXC related to the current feature.
        This subset of pixel cloud includes pixels for which self.obj_pixc.labels=in_label

        :param in_geom: geometry of the current feature
        :type in_geom: ogr.MultiPolygon
        :param in_pixc_index: list of indices of PIXC defining the current feature
        :type in_pixc_index: 1D-array of int
        :param in_classif: dictionary listing pixels selected or not for different processes (hull/wse/area computation)
        :type in_classif: dict of 1D-array of boolean
        :param flag_partial: =1 if water body is partially observed; =0 otherwise (default)
        :type flag_partial: boolean

        :return: out_attributes = lake attributes (None for each item which could not be computed)
        :rtype: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # Init output dictionary
        out_attributes = dict()

        # ========
        # 1 - Time
        # ========
        
        # Mean datetime of observation
        if np.isfinite(self.obj_pixc.nadir_time[in_pixc_index]).any():
            tmp_time = np.nanmean(self.obj_pixc.nadir_time[in_pixc_index])
        else:
            tmp_time = None
        out_attributes["time"] = my_tools.value_or_none(tmp_time)  # UTC time

        if np.isfinite(self.obj_pixc.nadir_time_tai[in_pixc_index]).any():
            tmp_time_tai = np.nanmean(self.obj_pixc.nadir_time_tai[in_pixc_index])
        else:
            tmp_time_tai = None
        out_attributes["time_tai"] = my_tools.value_or_none(tmp_time_tai)  # TAI time
        if out_attributes["time"] is None:
            out_attributes["time_str"] = None
        else:
            out_attributes["time_str"] = my_tools.convert_utc_to_str(out_attributes["time"])  # Time in UTC as a string
            
        # Related to nearest xovers
        if np.isfinite(self.obj_pixc.prev_xover[in_pixc_index]).any():
            tmp_prev_xover = np.nanmean(self.obj_pixc.prev_xover[in_pixc_index])
        else:
            tmp_prev_xover = None
        out_attributes["t_prv_xovr"] = my_tools.value_or_none(tmp_prev_xover)  # Time from previous xover
        
        if np.isfinite(self.obj_pixc.next_xover[in_pixc_index]).any():
            tmp_next_xover = np.nanmean(self.obj_pixc.next_xover[in_pixc_index])
        else:
            tmp_next_xover = None
        out_attributes["t_nxt_xovr"] = my_tools.value_or_none(tmp_next_xover)  # Time to next xover
        
        # =================================
        # 2 - Measured hydrology parameters
        # =================================
        
        # Area and uncertainties
        area_total, area_total_unc, area_detct, area_detct_unc = self.obj_pixc.compute_area_with_uncertainties(in_pixc_index[in_classif["4area"]],\
                                                                  method='composite')
        area_method = self.cfg.get("CONFIG_PARAMS", "AREA_METHOD")
        if area_method == "polygon":
            out_attributes["area_total"] = my_tools.get_area(in_geom) / 10**6
        else:
            out_attributes["area_total"] = my_tools.value_or_none(area_total / 10**6)
        out_attributes["area_tot_u"] = my_tools.value_or_none(area_total_unc / 10**6)
        out_attributes["area_detct"] = my_tools.value_or_none(area_detct / 10**6)
        out_attributes["area_det_u"] = my_tools.value_or_none(area_detct_unc / 10**6)
        out_attributes["n_area"] = np.sum(in_classif["4area"])
        
        # Water surface elevation and uncertainty
        out_attributes["n_wse"] = np.sum(in_classif["4wse"])
        mean_height, height_std, height_unc = self.obj_pixc.compute_height_with_uncertainties(in_pixc_index[in_classif["4wse"]])
        out_attributes["wse"] = my_tools.value_or_none(mean_height)
        out_attributes["wse_u"] = my_tools.value_or_none(height_std)
        out_attributes["wse_r_u"] = my_tools.value_or_none(height_unc)
        out_attributes["wse_std"] = my_tools.value_or_none(my_tools.compute_std_2sigma(self.obj_pixc.corrected_height\
                                        [in_pixc_index[in_classif["4wse"]]], name="wse_std"))
        
        # Metric of layover effect = layover area
        out_attributes["layovr_val"] = my_tools.value_or_none(my_tools.compute_mean_2sigma(
                                       self.obj_pixc.layover_impact[in_pixc_index[in_classif["4wse"]]], name="layovr_val"))
        
        # Median distance from PIXC to the satellite ground track
        out_attributes["xtrk_dist"] = my_tools.value_or_none(my_tools.compute_mean_2sigma(self.obj_pixc.cross_track[in_pixc_index]))
        
        # ============================
        # 3 - KaRIn Sigma0 Information
        # ============================
        
        # sigma0 and uncertainties
        sig0, sig0_std, sig0_u = self.obj_pixc.compute_sig0_with_uncertainties(in_pixc_index[in_classif["4wse"]])
        out_attributes["sig0"] = my_tools.value_or_none(sig0)
        out_attributes["sig0_u"] = my_tools.value_or_none(sig0_u)
        out_attributes["sig0_std"] = my_tools.value_or_none(sig0_std)
        
        # polarization of sigma0
        if self.product_type == "SP":
            out_attributes["pol"] = "+".join(list(set(self.obj_pixc.polarization)))
        else:
            out_attributes["pol"] = self.obj_pixc.polarization
        
        # ==========================
        # 4 - Geophysical references
        # ==========================
        
        # Geoid model height
        out_attributes["geoid_hght"] = my_tools.value_or_none(self.obj_pixc.compute_geophysical_ref('geoid_hght', in_pixc_index[in_classif["4wse"]]))
        # Earth tide
        out_attributes["solid_tide"] = my_tools.value_or_none(self.obj_pixc.compute_geophysical_ref('solid_tide', in_pixc_index[in_classif["4wse"]]))
        # Pole tide
        out_attributes["pole_tide"] = my_tools.value_or_none(self.obj_pixc.compute_geophysical_ref('pole_tide', in_pixc_index[in_classif["4wse"]]))
        # Load tide
        out_attributes["load_tidef"] = my_tools.value_or_none(self.obj_pixc.compute_geophysical_ref('load_tidef', in_pixc_index[in_classif["4wse"]]))
        out_attributes["load_tideg"] = my_tools.value_or_none(self.obj_pixc.compute_geophysical_ref('load_tideg', in_pixc_index[in_classif["4wse"]]))
        
        # =================================
        # 5 - Geophysical range corrections
        # =================================
        
        # Dry tropo corr
        out_attributes["dry_trop_c"] = my_tools.value_or_none(my_tools.compute_mean_2sigma(self.obj_pixc.model_dry_tropo_cor[in_pixc_index], \
                                                              name="dry_trop_c"))
        # Wet tropo corr
        out_attributes["wet_trop_c"] = my_tools.value_or_none(my_tools.compute_mean_2sigma(self.obj_pixc.model_wet_tropo_cor[in_pixc_index], \
                                                              name="wet_trop_c"))
        # Iono corr
        out_attributes["iono_c"] = my_tools.value_or_none(my_tools.compute_mean_2sigma(self.obj_pixc.iono_cor_gim_ka[in_pixc_index], name="iono_c"))
        
        # ==========================
        # 6 - Instrument corrections
        # ==========================
        
        # KaRIn correction from crossover cal processing evaluated for lake
        out_attributes["xovr_cal_c"] = my_tools.value_or_none(my_tools.compute_mean_2sigma(self.obj_pixc.height_cor_xover[in_pixc_index],\
                                                              name="xovr_cal_c"))
        
        # ======================
        # 7 - Quality indicators
        # ======================
        
        # Fractional area of pixels used to compute wse
        if out_attributes["n_wse"] > 0:
            wse_area = self.obj_pixc.compute_area(in_pixc_index[in_classif["4wse"]], method='composite')
            out_attributes["wse_frac"] = my_tools.value_or_none(wse_area / area_total)
        else:
            out_attributes["wse_frac"] = 0.0
        
        # Fractional area of water-near-land pixels
        if in_classif["water_edge"].any():
            edge_area = self.obj_pixc.compute_area(in_pixc_index[in_classif["water_edge"]], method='composite')
            out_attributes["edge_frac"] = my_tools.value_or_none(edge_area / area_total)
        else:
            out_attributes["edge_frac"] = 0.0
        
        # Fractional area of dark water
        out_attributes["dark_frac"] = (area_total - area_detct) / area_total
        if in_classif is not None:
            if in_classif["dark"] is None:
                out_attributes["dark_frac"] = 0.0
                
        # Fractional area of low-coherence pixels
        if in_classif["low_coherence"].any():
            low_area = self.obj_pixc.compute_area(in_pixc_index[in_classif["low_coherence"]], method='composite')
            out_attributes["low_frac"] = my_tools.value_or_none(low_area / area_total)
        else:
            out_attributes["low_frac"] = 0.0
                
        # Fractional area of pixels with specular ringing
        if in_classif["specular_ringing"].any():
            spec_area = self.obj_pixc.compute_area(in_pixc_index[in_classif["specular_ringing"]], method='composite')
            out_attributes["spec_frac"] = my_tools.value_or_none(spec_area / area_total)
        else:
            out_attributes["spec_frac"] = 0.0
        
        # Partial flag
        out_attributes["partial_f"] = flag_partial
        
        # Quality of the cross-over calibration
        out_attributes["xovr_cal_q"] = 0
        # Depends on geolocation_qual variable 
        # bit 23 = xovercal_missing => xovr_cal_q = 2
        # bit 6 = xovercal_suspect => xovr_cal_q = 1
        for flag_mask, out_value in zip(['00000000100000000000000000000000', '00000000000000000000000001000000'], [2,1]):
            mask_nok = np.uint32(int(flag_mask, 2))
            tmp_comp = np.bitwise_and(self.obj_pixc.geolocation_qual[in_pixc_index], mask_nok)
            nb_nok = np.sum(tmp_comp > 0)
            if nb_nok > 0:
                out_attributes["xovr_cal_q"] = out_value
                break
        
        # Summary and bitwise quality flags
        out_attributes["quality_f"], out_attributes["qual_f_b"] = self.compute_quality_flags(in_pixc_index, 
                                                                                              out_attributes["wse_std"],
                                                                                              out_attributes["xovr_cal_c"],
                                                                                              out_attributes["xovr_cal_q"])
        
        return out_attributes
                
    def compute_quality_flags(self, in_pixc_index, in_wse_std, in_xovr_cal_c, in_xovr_cal_q):
        """
        Computes quality flags from PIXC related to the current feature.
        This subset of pixel cloud includes pixels for which self.obj_pixc.labels=in_label

        :param in_pixc_index: list of indices of PIXC defining the current feature
        :type in_pixc_index: 1D-array of int
        :param in_wse_std: standard deviation of the water surface elevation of all the pixels composing the lake
        :type in_wse_std: float
        :param in_xovr_cal_c: cross-over calibration correction
        :type in_xovr_cal_c: float
        :param in_xovr_cal_q: flag that indicates the quality of the cross-over calibration
        :type in_xovr_cal_q: float

        :return: out_quality_f = summary quality flag
        :rtype: out_quality_f = int
        :return: out_qual_f_b = quality bitflag
        :rtype: out_qual_f_b = int
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # Number of pixels of the current feature
        nb_pixc = len(in_pixc_index)
        
        # 1 - Populate quality bitflag
        out_qual_f_b = 0
        #­ 1.1 - classification_qual_suspect
        idx1 = np.where(self.obj_pixc.classification_qual[in_pixc_index] > 0)[0]
        idx2 = np.where(self.obj_pixc.classification_qual[in_pixc_index] < 65536)[0]
        idx = np.intersect1d(idx1, idx2)
        if len(idx)/nb_pixc > 0.25:
            out_qual_f_b += 1
        # 1.2 - geolocation_qual_suspect
        idx1 = np.where(self.obj_pixc.geolocation_qual[in_pixc_index] >= 4)[0]
        idx2 = np.where(self.obj_pixc.geolocation_qual[in_pixc_index] < 65536)[0]
        idx = np.intersect1d(idx1, idx2)
        if len(idx)/nb_pixc > 0.25:
            out_qual_f_b += 2
        # 1.3 - few_open_water_suspect
        idx = np.where(self.obj_pixc.classif[in_pixc_index] == 4)[0]
        if len(idx)/nb_pixc < 0.5:
            out_qual_f_b += 8
        # 1.4 - wse_std_suspect
        if (in_wse_std is not None) and (in_wse_std > 3.0):
            out_qual_f_b += 32
        # 1.5 - diff_from_pld_area_suspect
        # Only for LakeTile/SP_Prior layer
        # 1.6 - xovr_cal_suspect
        if in_xovr_cal_q == 1:
            out_qual_f_b += 256
        # 1.7 - no_prior_suspect
        tmp_comp = np.bitwise_and(self.obj_pixc.classification_qual[in_pixc_index], 4)
        idx = np.argwhere(tmp_comp > 0)
        if len(idx)/nb_pixc > 0.25:
            out_qual_f_b += 1024
        # 1.8 - water_false_detection_rate_suspect
        tmp_comp = np.bitwise_and(self.obj_pixc.classification_qual[in_pixc_index], 16)
        idx = np.argwhere(tmp_comp > 0)
        if len(idx)/nb_pixc > 0.25:
            out_qual_f_b += 2048
        # 1.9 - classification_qual_degraded
        idx = np.where(self.obj_pixc.classification_qual[in_pixc_index] >= 65536)[0]
        if len(idx)/nb_pixc > 0.25:
            out_qual_f_b += 32768
        # 1.10 - geolocation_qual_degraded
        idx = np.where(self.obj_pixc.geolocation_qual[in_pixc_index] >= 65536)[0]
        if len(idx)/nb_pixc > 0.25:
            out_qual_f_b += 65536	
        # 1.11 - low_coh_degraded and low_coh_bad
        idx1 = np.where(self.obj_pixc.classif[in_pixc_index] == 6)[0]
        idx2 = np.where(self.obj_pixc.classif[in_pixc_index] == 7)[0]
        ratio = (len(idx1)+len(idx2))/nb_pixc
        if ratio > 0.75:
            out_qual_f_b += 8388608
        elif ratio > 0.5:
            out_qual_f_b += 262144
        # 1.12 - specular_ringing_degraded and specular_ringing_bad
        tmp_comp = np.bitwise_and(self.obj_pixc.classification_qual[in_pixc_index], 524288)
        idx = np.argwhere(tmp_comp > 0)
        if len(idx)/nb_pixc > 0.5:
            out_qual_f_b += 67108864
        elif len(idx)/nb_pixc > 0.25:
            out_qual_f_b += 524288
        # 1.13 - xovr_cal_bad
        if in_xovr_cal_q == 2:
            out_qual_f_b += 16777216
        
        # 2 - Deduce summary bitflag
        if out_qual_f_b == 0:
            out_quality_f = 0
        elif out_qual_f_b < 16384:
            out_quality_f = 1
        elif out_qual_f_b < 4194304:
            out_quality_f = 2
        else:
            out_quality_f = 3
        
        return out_quality_f, out_qual_f_b
                
    def compute_prior_attributes(self, in_obj_plake, in_pixc_index, in_feat_area, in_classif=None):
        """
        Computes attributes specific to _Prior shapefile, from the PIXC entirely inside the PLD polygon.

        :param in_obj_plake: PLD lake object
        :type in_obj_plake: lake_db.PriorLake
        :param in_pixc_index: list of indices of PIXC defining the current feature
        :type in_pixc_index: 1D-array of int
        :param in_feat_area: total area of the observed feature (area_total, in km2)
        :type in_feat_area: float
        :param in_classif: dictionary listing pixels selected or not for different processes (hull/wse/area computation)
        :type in_classif: dict of 1D-array of boolean

        :return: out_attributes = lake prior attributes (None for each item which could not be computed)
        :rtype: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # Init output dictionary
        out_attributes = dict()
        
        # ===================================================
        # 1 - Compute subset of pixels inside the PLD polygon
        # ===================================================
        
        # 1.1.1 - Convert PLD polygon from OGR.Polygon into Shapely.Polygon
        geojson = json.loads(in_obj_plake.geom.ExportToJson())
        pld_geom_shapely = geometry.shape(geojson)
        # 1.1.2 - Prepare PLD polygon to optimize testing
        pld_geom_prepared = prepared.prep(pld_geom_shapely)
        
        # 1.2 - Retrieve PIXCVec coordinates
        pixcvec_index = self.get_pixcvec_index_from_pixc_index(in_pixc_index)  # Indices in PIXCVec ref
        pixcvec_long = self.obj_pixcvec.longitude_vectorproc[pixcvec_index]  # Indices between 0 and N-1 (N=nb elements of in_pixc_index)
        pixcvec_lat = self.obj_pixcvec.latitude_vectorproc[pixcvec_index]  # Indices between 0 and N-1
        
        # 1.3 - High level filtering over PLD polygon bounding box
        xmin, xmax, ymin, ymax = in_obj_plake.geom.GetEnvelope()
        bbox_mask = ((pixcvec_long >= xmin) & (pixcvec_long <= xmax) & (pixcvec_lat >= ymin) & (pixcvec_lat <= ymax))
        sub_pixcvec_bbox = np.flatnonzero(bbox_mask)  # Indices among 0 and N-1
        logger.debug("> %d / %d pixels after filtering over PLD polygon BOUNDING BOX" % (len(sub_pixcvec_bbox), len(in_pixc_index)))
        
        # 1.4 - Retrieve indices of pixels inside the PLD polygone
        sub_pixcvec_bbox_long = pixcvec_long[sub_pixcvec_bbox]
        sub_pixcvec_bbox_lat = pixcvec_lat[sub_pixcvec_bbox]
        mask = vectorized.contains(pld_geom_prepared, sub_pixcvec_bbox_long, sub_pixcvec_bbox_lat)   # bool array
        inside_idx = sub_pixcvec_bbox[mask]  # Indices among 0 and N-1
        sub_pixcvec_pld = np.asarray(inside_idx, dtype=np.int64)  # Indices among 0 and N-1
        sub_pixc_pld = in_pixc_index[sub_pixcvec_pld]  # Indices in PIXC ref
        nb_pixc_pld = len(sub_pixc_pld)
        logger.debug("> %d / %d pixels INSIDE the PLD polygon" % (nb_pixc_pld, len(in_pixc_index)))
        
        # 1.5 - Determine pixels inside the PLD polygon to be used to compute wse2
        sub_classif_pld_4wse = in_classif["4wse"][sub_pixcvec_pld]  # Array of True/False, same length as sub_pixcvec_pld
        sub_pixc_pld_4wse = sub_pixc_pld[sub_classif_pld_4wse]  # Array of indices in PIXC ref, same length as nb of True in sub_classif_pld_4wse
        nb_pixc_pld_4wse = len(sub_pixc_pld_4wse)
        logger.debug("> %d / %d pixels kept to compute CONSTRAINED WSE (wrt %d for observed feature)" % (nb_pixc_pld_4wse, nb_pixc_pld, len(in_pixc_index[in_classif["4wse"]])))
        # Highlight in PIXCVec which pixels will be used to compute wse2
        self.obj_pixcvec.used_for_lake_wse2[pixcvec_index[sub_pixcvec_pld]] = np.array(sub_classif_pld_4wse, dtype=int)
        
        # ======================
        # 2 - Quality indicators
        # ======================
        
        # Fractional area of pixels within the PLD polygon
        # NB: output of compute area in m2 whereas in_feat_area in km2
        pld_area = 0.0
        if nb_pixc_pld > 0:
            pld_area = self.obj_pixc.compute_area(sub_pixc_pld, method='composite') / 10**6
        if in_feat_area is not None :
            out_attributes["pld_frac"] = pld_area / in_feat_area
        
        # Fractional area of pixels used to compute wse2
        # NB: output of compute area in m2 whereas in_feat_area in km2
        wse_area = 0.0
        if nb_pixc_pld_4wse > 0:
            wse_area = self.obj_pixc.compute_area(sub_pixc_pld_4wse, method='composite') / 10**6
        if in_feat_area is not None :
            out_attributes["wse2_frac"] = wse_area / in_feat_area
        
        # =========================================================
        # 3 - Extreme case when no PIXC 4wse inside the PLD polygon
        # =========================================================
            
        # Number of PIXC available to compute wse2
        out_attributes["n_wse2"] = nb_pixc_pld_4wse
        
        if nb_pixc_pld_4wse == 0:
            logger.debug("WARNING Still NO valid mandatory pixel to compute wse2 (%d) => DON'T COMPUTE additional prior attributes" % nb_pixc_pld_4wse)
            return out_attributes
        
       # ==========================
        # 4 - Geophysical references
        # ==========================
        
        # Difference between SWOT-observed mean lake surface and the EGM2008 geoid model
        # Interpolate MLS over pixels inside the PLD polygon
        pixcvec_mls = self.obj_mls.interpolate(in_obj_plake.lake_id, pixcvec_long[sub_pixcvec_pld], pixcvec_lat[sub_pixcvec_pld])  # Same length as sub_*_pld arrays
        # Compute its difference with geoid model (array of 0 for small lakes which haven't MLS)
        if pixcvec_mls is None:
            pixcvec_delta_mls = np.zeros(nb_pixc_pld, dtype=float)
        else:
            pixcvec_delta_mls = self.obj_pixc.geoid[sub_pixc_pld] - pixcvec_mls  # Same length as sub_*_pld arrays
            
        # Update PIXCVec pixels with delta_mls
        self.obj_pixcvec.delta_mean_lake_surface[pixcvec_index[sub_pixcvec_pld]] = pixcvec_delta_mls
            
        # Compute the lake-averaged value
        out_attributes["delta_mls"] = my_tools.value_or_none(self.obj_pixc.compute_geophysical_ref('delta_mls', sub_pixc_pld_4wse, delta_mls=pixcvec_delta_mls[sub_classif_pld_4wse]))
        
        # =================================
        # 5 - Measured hydrology parameters
        # =================================
        
        # Constrained water surface elevation and uncertainties
        mean_height, height_std, height_unc = self.obj_pixc.compute_height_with_uncertainties(sub_pixc_pld_4wse, delta_mls=pixcvec_delta_mls[sub_classif_pld_4wse])
        out_attributes["wse2"] = my_tools.value_or_none(mean_height)
        out_attributes["wse2_u"] = my_tools.value_or_none(height_std)
        out_attributes["wse2_r_u"] = my_tools.value_or_none(height_unc)
        out_attributes["wse2_std"] = my_tools.value_or_none(my_tools.compute_std_2sigma(self.obj_pixc.corrected_height[sub_pixc_pld_4wse], name="wse_std"))
        
        # ============================
        # 6 - KaRIn Sigma0 Information
        # ============================
        
        # Constrained sigma0 and uncertainties
        sig0, sig0_std, sig0_u = self.obj_pixc.compute_sig0_with_uncertainties(sub_pixc_pld_4wse)
        out_attributes["sig0_2"] = my_tools.value_or_none(sig0)
        out_attributes["sig0_2_u"] = my_tools.value_or_none(sig0_u)
        out_attributes["sig0_2_std"] = my_tools.value_or_none(sig0_std)
            
        return out_attributes

    # -----------------------------------
    # Functions specific to link with PLD
    # -----------------------------------
    
    def link_obs_geom_to_pld_geom(self, in_geom, in_lon, in_lat):
        """
        Link observed geometry to PLD lakes geometries
        and retrieve corresponding information
        
        :param in_geom: observed geometry
        :type in_geom: OGRPolygon
        :param in_lon: improved longitudes vector for PIXC of the feature
        :type in_lon: 1D-array of float
        :param in_lat: improved latitudes vector for PIXC of the feature
        :type in_lat: 1D-array of float
        
        :return: out_lake_id = list of identifiers of PLD lakes intersecting the observed lake
        :rtype: out_lake_id = string
        :return: out_overlap = list of %age of PLD lake overlapping the observed lake
        :rtype: out_overlap = string
        :return: out_pixcvec_lakeid = lake identifiers from PLD to which each PIXC corresponds one-to-one
        :rtype: out_pixcvec_lakeid = list of string
        """
        
        # 1 - Intersect observed geometry with polygons of PLD
        list_prior, list_fraction, out_pixcvec_lakeid = self.obj_lake_db.link_to_db(in_geom, in_lon, in_lat)

        # 3 - Form lake_id attribute
        if list_prior:
            out_lake_id = ';'.join(list_prior)
        else:
            out_lake_id = "no_data"
            
        # 4 - Form overlap attribute
        if list_fraction:
            out_overlap = ';'.join(list_fraction)
        else:
            out_overlap = "no_data"
            
        return out_lake_id, out_overlap, out_pixcvec_lakeid
    
    def set_pld_attributes(self, in_pld_infos):
        """
        Set p_ attributes of all obs features linked to current PLD lake (i.e. currently selected in the _Obs layer)
        to prior values of current PLD lake 
        
        :param in_pld_infos: values of available p_ attributes
        :type in_pld_infos: dict
        """
        # log logger = logging.getLogger(self.__class__.__name__)
        # log logger.debug("- start -")
        
        for obs_lake in self.content_obs.layer:
            
            # Loop only over needed attributes from PLD
            for key in my_var.PLD_FIELD_TO_KEEP_IN_OBS:
                
                # Retrieve current value in obs
                tmp_value = obs_lake.GetField(str(key))
                
                if in_pld_infos[key] is not None:  # No need to do anything if no value in PLD
                    
                    # Process depending on attribute
                    if key.startswith("ice_"):  # Keep max value
                        value_to_write = my_tools.value_or_none(np.nanmax([int(tmp_value), int(in_pld_infos[key])]))
                    elif key == "p_res_id":
                        if int(in_pld_infos[key]) != 0:
                            value_to_write = in_pld_infos[key]
                        else:
                            value_to_write = tmp_value
                    else:  # Concat strings
                        if tmp_value == my_var.FV_STRING_SHP:
                            value_to_write = in_pld_infos[key]
                        else:
                            value_to_write = "%s;%s" % (tmp_value, in_pld_infos[key])
                            
                    # Update obs attribute
                    obs_lake.SetField(str(key), str(value_to_write))
                    
            # Rewrite obs feature with updated attributes
            self.content_obs.layer.SetFeature(obs_lake)

    # ------------------------------------------------
    # Functions specific to storage change computation
    # ------------------------------------------------
    
    def compute_direct_storage_change(self, in_obj_plake, in_pixc_index, in_classif, in_prior_attributes):
        """
        Head function for storage change computation with the DIRECT approach
        Run specific child function depending on the STOCC_INPUT config parameter
        
        :param in_obj_plake: PLD lake object
        :type in_obj_plake: lake_db.PriorLake
        :param in_pixc_index: indices of PIXC related to PLD lake
        :type in_pixc_index: Numpy 1D-array of int
        :param in_classif: dictionary listing pixels selected or not for different processes (hull/wse/area computation)
        :type in_classif: dict of 1D-array of boolean
        :return: in_prior_attributes = attributes of prior feature related to PLD lake
        :rtype: in_prior_attributes = dict
        
        :return: out_storage_values = storage change values related to current PLD lake
        :rtype: out_storage_values = dict
        """
        
        # 1 - Retrieve config parameter for choice of input data for storage change computation
        stocc_input = self.cfg.get("CONFIG_PARAMS", "STOCC_INPUT")
        
        # 2 - Run specific function depending on its value
        if stocc_input == "obs":
            out_storage_values = self.compute_direct_storage_change_obs(in_obj_plake, in_pixc_index, in_classif, in_prior_attributes)
        else:
            out_storage_values = self.compute_direct_storage_change_pld(in_obj_plake, in_prior_attributes)
            
        return out_storage_values
    
    def compute_direct_storage_change_pld(self, in_obj_plake, in_prior_attributes):
        """
        Compute storage change precisely from WSE and area average over the PLD lake
        With this method, storage change is not set to related observed features
        
        :param in_obj_plake: PLD lake object
        :type in_obj_plake: lake_db.PriorLake
        :return: in_prior_attributes = attributes of prior feature related to PLD lake
        :rtype: in_prior_attributes = dict
        
        :return: out_storage_values = storage change values related to current PLD lake
        :rtype: out_storage_values = dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # 0 - Init output
        out_storage_values = dict()
        out_storage_values["ds1_l"] = my_var.FV_REAL
        out_storage_values["ds1_l_u"] = my_var.FV_REAL
        out_storage_values["ds1_q"] = my_var.FV_REAL
        out_storage_values["ds1_q_u"] = my_var.FV_REAL
        
        # 1 - Build dictionary going in input of storage change function
        list_obs = dict()
        list_obs["pld"] = dict()
        list_obs["pld"]["area"] = in_prior_attributes["area_total"]
        list_obs["pld"]["area_u"] = in_prior_attributes["area_tot_u"]
        list_obs["pld"]["wse"] = in_prior_attributes["wse"]
        list_obs["pld"]["wse_u"] = in_prior_attributes["wse_u"]
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
        out_storage_values["ds1_l"], out_storage_values["ds1_l_u"], out_storage_values["ds1_q"], out_storage_values["ds1_q_u"] = \
            in_obj_plake.run_stocc(list_obs)
        
        return out_storage_values

    def compute_direct_storage_change_obs(self, in_obj_plake, in_pixc_index, in_classif, in_prior_attributes):
        """
        Compute storage change precisely from WSE and area of all observed features linked the the PLD lake
        Set storage change values for all observed features linked to current PLD lake
        NB: these observed features have been selected in the _Obs layer before the use of this function

        Envisionned cases:
            - 1 prior lake <=> 1 observed lake
            - 1 prior lake <=> 2 or more observed lakes
            - 1 observed lake <=> 2 or more prior lakes
            - mixed lakes...
        
        :param in_obj_plake: PLD lake object
        :type in_obj_plake: lake_db.PriorLake
        :param in_pixc_index: indices of PIXC related to PLD lake
        :type in_pixc_index: Numpy 1D-array of int
        :param in_classif: dictionary listing pixels selected or not for different processes (hull/wse/area computation)
        :type in_classif: dict of 1D-array of boolean
        :return: in_prior_attributes = attributes of prior feature related to PLD lake
        :rtype: in_prior_attributes = dict
        
        :return: out_storage_values = storage change values related to current PLD lake
        :rtype: out_storage_values = dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 0 - Init output
        out_storage_values = dict()
        out_storage_values["ds1_l"] = my_var.FV_REAL
        out_storage_values["ds1_l_u"] = my_var.FV_REAL
        out_storage_values["ds1_q"] = my_var.FV_REAL
        out_storage_values["ds1_q_u"] = my_var.FV_REAL
        
        # Nb of observed lakes related to in_obj_plake
        nb_obs_lake = self.content_obs.layer.GetFeatureCount()
        logger.debug("Case PLD lake -> %d obs lakes" % nb_obs_lake)
        
        # 1 - Build dictionary going in input of storage change function
        
        # 1.0 - Compute prior overlap coefficients
        prior_overlap = dict()
        tmp_prior_obsid = in_prior_attributes["obs_id"].split(";")
        tmp_prior_overlap = np.array(in_prior_attributes["overlap"].split(";"), dtype="float")
        sum_tmp_prior_overlap = np.sum(tmp_prior_overlap)
        for obs_id, overlap in zip(tmp_prior_obsid, tmp_prior_overlap):
            prior_overlap[obs_id] = overlap / sum_tmp_prior_overlap
            
        list_obs = dict()
        for obs_lake in self.content_obs.layer:

            # 1.1 - Retrieve lake feature ID and init dedicated dict
            obs_id = obs_lake.GetField(str("obs_id"))
            list_obs[obs_id] = dict()
            
            # 1.2 - Compute related wse and area values
            if ";" in obs_lake.GetField(str("lake_id")):
                # Case the observed lake is related to 2 or more prior lakes
                if nb_obs_lake == 1:
                    # Case 1 obs lake <-> N PLD lakes
                    logger.debug("> sub-case obs lake -> N PLD lakes")
                    list_obs[obs_id]["area"] = in_prior_attributes["area_total"]
                    list_obs[obs_id]["area_u"] = in_prior_attributes["area_tot_u"]
                    list_obs[obs_id]["wse"] = in_prior_attributes["wse"]
                    list_obs[obs_id]["wse_u"] = in_prior_attributes["wse_u"]
                else:
                    # Case multi obs <-> PLD lakes associations
                    logger.debug("> sub-case multi obs <-> PLD lakes associations")
                    list_obs[obs_id]["area"], list_obs[obs_id]["area_u"], tmp1, tmp2 = \
                        self.obj_pixc.compute_area_with_uncertainties(in_pixc_index[in_classif["4area"]], flag_all=False)
                    list_obs[obs_id]["wse"], list_obs[obs_id]["wse_u"], tmp1 = \
                        self.obj_pixc.compute_height_with_uncertainties(in_pixc_index[in_classif["4wse"]])
                    
                # Retrieve storage change previous values
                list_obs[obs_id]["lake_id"] = obs_lake.GetField(str("lake_id")).split(";")
                list_obs[obs_id]["overlap"] = np.array(obs_lake.GetField(str("overlap")).split(";"), dtype="float")
                list_obs[obs_id]["ds1_l"] = float(obs_lake.GetField(str("ds1_l")))
                list_obs[obs_id]["ds1_l_u"] = float(obs_lake.GetField(str("ds1_l_u")))
                list_obs[obs_id]["ds1_q"] = float(obs_lake.GetField(str("ds1_q")))
                list_obs[obs_id]["ds1_q_u"] = float(obs_lake.GetField(str("ds1_q_u")))
                if list_obs[obs_id]["ds1_l"] < -9e11:
                    list_obs[obs_id]["ds1_l"] = 0.0
                    list_obs[obs_id]["ds1_l_u"] = 0.0
                    list_obs[obs_id]["ds1_q"] = 0.0
                    list_obs[obs_id]["ds1_q_u"] = 0.0
                
            else:
                # Case the observed lake is related to only 1 prior lake
                logger.debug("> sub-case obs lake -> 1 PLD lake")
                list_obs[obs_id]["area"] = float(obs_lake.GetField(str("area_total")))
                list_obs[obs_id]["area_u"] = float(obs_lake.GetField(str("area_tot_u")))
                list_obs[obs_id]["wse"] = float(obs_lake.GetField(str("wse")))
                list_obs[obs_id]["wse_u"] = float(obs_lake.GetField(str("wse_u")))
            
            # 1.3 - Compute alpha coefficient from prior overlap attribute
            if nb_obs_lake == 1:
                list_obs[obs_id]["alpha"] = 1.0
            else:
                list_obs[obs_id]["alpha"] = prior_overlap[obs_id]
                
            # 1.4 - Test if extreme values
            flag_ok = True
            for key in ["area", "wse"]:
                if (list_obs[obs_id][key] is None) or (not np.isfinite(list_obs[obs_id][key])) or (list_obs[obs_id][key] < -9e11):
                    logger.debug("Lake with obs_id={} has {}={} => removed from storage change computation".format(obs_id, key,\
                                                                                                                   list_obs[obs_id][key]))
                    del list_obs[obs_id]
                    flag_ok = False
                    break
            if flag_ok:
                for key in ["area_u", "wse_u"]:
                    if (list_obs[obs_id][key] is None) or (not np.isfinite(list_obs[obs_id][key])) or (list_obs[obs_id][key] < -9e11):
                        key2 = key.replace("_u", "")
                        logger.debug("Lake with obs_id={} has {}={} => {} set to {}={}".format(obs_id, key, list_obs[obs_id][key], key, \
                                                                                                               key2, list_obs[obs_id][key2]))
                        list_obs[obs_id][key] = list_obs[obs_id][key2]
            
        # Reset reading
        self.content_obs.layer.ResetReading()
            
        # 2 - Compute storage change for PLD lake
        out_storage_values["ds1_l"], out_storage_values["ds1_l_u"], out_storage_values["ds1_q"], out_storage_values["ds1_q_u"] = \
            in_obj_plake.run_stocc(list_obs)
                
        # 3 - Compute storage change for observed features and set values
        for obs_lake in self.content_obs.layer:
            
            # 3.1 - Retrieve lake feature ID and init dedicated dict
            obs_id = obs_lake.GetField(str("obs_id"))
            
            # 3.2 - Set fields
            if "overlap" in list_obs[obs_id].keys():
                # Case the observed lake is related to 2 or more prior lakes
                for indi, lakeid in enumerate(list_obs[obs_id]["lake_id"]):
                    if lakeid == in_obj_plake.lake_id:
                        coeff = list_obs[obs_id]["overlap"][indi] / np.sum(list_obs[obs_id]["overlap"])
                        break
                obs_lake.SetField(str("ds1_l"), list_obs[obs_id]["ds1_l"] + out_storage_values["ds1_l"]*list_obs[obs_id]["alpha"])
                obs_lake.SetField(str("ds1_l_u"), list_obs[obs_id]["ds1_l_u"] + out_storage_values["ds1_l_u"]*coeff*list_obs[obs_id]["alpha"])
                obs_lake.SetField(str("ds1_q"), list_obs[obs_id]["ds1_q"] + out_storage_values["ds1_q"]*list_obs[obs_id]["alpha"])
                obs_lake.SetField(str("ds1_q_u"), list_obs[obs_id]["ds1_q_u"] + out_storage_values["ds1_q_u"]*coeff*list_obs[obs_id]["alpha"])
            else:
                # Case the observed lake is related to only 1 prior lake
                obs_lake.SetField(str("ds1_l"), out_storage_values["ds1_l"]*list_obs[obs_id]["alpha"])
                obs_lake.SetField(str("ds1_l_u"), out_storage_values["ds1_l_u"]*list_obs[obs_id]["alpha"])
                obs_lake.SetField(str("ds1_q"), out_storage_values["ds1_q"]*list_obs[obs_id]["alpha"])
                obs_lake.SetField(str("ds1_q_u"), out_storage_values["ds1_q_u"]*list_obs[obs_id]["alpha"])
            
            # 3.3 - Rewrite feature with storage change values
            self.content_obs.layer.SetFeature(obs_lake)
            
        # Reset reading
        self.content_obs.layer.ResetReading()
        
        return out_storage_values
    
    def compute_incremental_storage_change(self, in_obj_plake, in_pixc_index, in_prior_attributes):
        """
        Compute storage change with the INCREMENTAL approach
        
        :param in_obj_plake: PLD lake object
        :type in_obj_plake: lake_db.PriorLake
        :param in_pixc_index: indices of PIXC related to PLD lake
        :type in_pixc_index: Numpy 1D-array of int
        :return: in_prior_attributes = attributes of prior feature related to PLD lake
        :rtype: in_prior_attributes = dict
        
        :return: out_storage_values = storage change values related to current PLD lake
        :rtype: out_storage_values = dict
        """
        
        # 0 - Init output
        out_storage_values = dict()
        out_storage_values["ds2_l"] = my_var.FV_REAL
        out_storage_values["ds2_l_u"] = my_var.FV_REAL
        out_storage_values["ds2_q"] = my_var.FV_REAL
        out_storage_values["ds2_q_u"] = my_var.FV_REAL
            
        return out_storage_values

    # ----------------------------------------
    # Functions specific to PIXCVec update
    # ----------------------------------------
    
    def get_pixcvec_index_from_pixc_index(self, in_pix_index):
        """
        Compute PIXCVec indices corresponding PIXC indices depending on product type
        
        :param in_pix_index: indices of PIXC 
        :type in_pix_index: 1D-array of int
        
        :return: out_pixcvec_indices = corresponding indices of PIXCVec
        :type out_pixcvec_indices: 1D-array of int
        """
        
        if self.product_type == "SP":
            out_pixcvec_indices = in_pix_index
        else:
            out_pixcvec_indices = self.obj_pixc.selected_index[in_pix_index]
            
        return out_pixcvec_indices
    
    def update_pixcvec_with_hconstr_geoloc(self, in_pix_index, in_lon, in_lat, in_height):
        """
        Save improved values in PIXCVec object. Indices depends on the product type
        
        :param in_pix_index: indices of PIXC defining the feature
        :type in_pix_index: 1D-array of int
        :param in_lon: improved longitudes to save
        :type in_lon: 1D-array of float
        :param in_lat: improved latitudes to save
        :type in_lat: 1D-array of float
        :param in_height: improved heights to save
        :type in_height: 1D-array of float
        """
        
        # 1 - Compute PIXCVec indices corresponding PIXC indices
        pixcvec_index = self.get_pixcvec_index_from_pixc_index(in_pix_index)
            
        # 2 - Save the improved values
        self.obj_pixcvec.longitude_vectorproc[pixcvec_index] = in_lon
        self.obj_pixcvec.latitude_vectorproc[pixcvec_index] = in_lat
        self.obj_pixcvec.height_vectorproc[pixcvec_index] = in_height
    
    def update_pixcvec_with_used_for_lake_wse(self, in_pix_index, in_classif_4wse):
        """
        Save improved values in PIXCVec object. Indices depends on the product type
        
        :param in_pix_index: indices of PIXC defining the feature
        :type in_pix_index: 1D-array of int
        :param in_classif_4wse: =True if the pixel is used to compute wse, =False otherwise
        :type in_classif_4wse: 1D-array of boolean
        """
        
        # 1 - Compute PIXCVec indices corresponding PIXC indices
        pixcvec_index = self.get_pixcvec_index_from_pixc_index(in_pix_index)
            
        # 2 - Save the flag values
        self.obj_pixcvec.used_for_lake_wse[pixcvec_index] = np.array(in_classif_4wse, dtype=int)
    
    def update_pixcvec_with_ids(self, in_pix_index, in_obs_id, in_pixcvec_lakeid):
        """
        Update PIXCVec lake_id and obs_id for current lake product
        
        :param in_pix_index: list of indices of PIXC belonging to the same feature
        :type in_pix_index: 1D-array of int
        :param in_obs_id: obs_id identifier of this feature
        :type in_obs_id: string
        :param in_pixcvec_lakeid: list of indices of PIXC belonging to the same feature
        :type in_pixcvec_lakeid: 1D-array of int
        """
        
        # 1 - Compute PIXCVec indices corresponding PIXC indices
        pixcvec_index = self.get_pixcvec_index_from_pixc_index(in_pix_index)
            
        # 2 - Update identifiers
        self.obj_pixcvec.obs_id[pixcvec_index] = in_obs_id
        self.obj_pixcvec.lake_id[pixcvec_index] = in_pixcvec_lakeid
    
    def update_pixcvec_with_ice_clim_f(self, in_pix_index, in_ice_clim_f):
        """
        Update PIXCVec ice_clim_f for current lake product
        
        :param in_pix_index: list of indices of PIXC belonging to the same feature
        :type in_pix_index: 1D-array of int
        :param in_ice_clim_f: ice clim flag of this feature
        :type in_ice_clim_f: string
        """
        
        # 1 - Compute PIXCVec indices corresponding PIXC indices
        pixcvec_index = self.get_pixcvec_index_from_pixc_index(in_pix_index)
            
        # 2 - Update ice clim flag
        self.obj_pixcvec.ice_clim_f[pixcvec_index] = in_ice_clim_f
        
    def update_pixcvec_with_ice_dyn_f(self, in_pix_index):
        """
        Update PIXCVec ice_dyn_f for current lake product
        
        :param in_pix_index: list of indices of PIXC belonging to the same feature
        :type in_pix_index: 1D-array of int
        """
        
        # 1 - Compute PIXCVec indices corresponding PIXC indices
        pixcvec_index = self.get_pixcvec_index_from_pixc_index(in_pix_index)
            
        # 2 - Update ice dyn flag
        self.obj_pixcvec.ice_dyn_f[pixcvec_index] = self.obj_pixc.ice_dyn_f[in_pix_index]

    # ----------------------------------------
    # Function dedicated to high-level stats
    # ----------------------------------------
    
    def compute_obs_stats(self):
        """
        Compute sum of values listed in self.compare_stats_params for the _Obs shapefile
        """
    
        # 1 - Compute and store stats
        for lake in self.content_obs.layer:
            for param in self.compare_stats_params:
                tmp_value = float(lake.GetField(param))
                if tmp_value > my_var.FV_REAL:
                    self.compare_stats["obs"][param] += tmp_value
                
        # 2 - Reinitialize reader pointer
        self.content_obs.layer.ResetReading()

    # ----------------------------------------
    # Writing functions
    # ----------------------------------------
    
    def write_obs_file(self, in_filename, in_proc_metadata):
        """
        Write the observation-oriented file, i.e.
        observed water features related to at least one PLD lake
        
        :param in_filename: full path of the output file
        :type in_filename: string
        :param in_proc_metadata: processing metadata
        :type in_proc_metadata: dict
        """
        # log logger = logging.getLogger(self.__class__.__name__)
        # log logger.debug("- start -")
        
        # 1 - Write them in the output shapefile
        my_shp.write_mem_layer_as_shp(self.content_obs.layer, in_filename)
        
        # 2 - Estimate time_coverage_start and time_coverage_end
        self.content_obs.layer.ResetReading()
        time_utc_values = []
        for feature in self.content_obs.layer:
            time_utc_feat = feature.GetField("time")
            if time_utc_feat > 0.0:
                time_utc_values.append(time_utc_feat)
        self.content_obs.layer.ResetReading()
        time_str_dict = dict()
        if len(time_utc_values) > 0:
            time_str_dict["time_coverage_start"] = my_tools.convert_utc_to_str(min(time_utc_values), in_format=2)
            time_str_dict["time_coverage_end"] = my_tools.convert_utc_to_str(max(time_utc_values), in_format=2)
        else:
            time_str_dict["time_coverage_start"] = "None"
            time_str_dict["time_coverage_end"] = "None"
        
        # 3 - Write XML metadatafile for shapefile
        #logger.debug("Writing associated metadata file = %s.xml" % in_filename)
        if self.obj_pixc is None:
            self.content_obs.update_and_write_metadata("%s.xml" % in_filename, 
                                                       in_proc_metadata={**in_proc_metadata, **time_str_dict})
        else:
            self.content_obs.update_and_write_metadata("%s.xml" % in_filename, 
                                                         in_inprod_metadata=self.obj_pixc.pixc_metadata,
                                                         in_proc_metadata={**in_proc_metadata, **time_str_dict})
    
    def write_prior_file(self, in_filename, in_proc_metadata):
        """
        Write the PLD-oriented file, i.e. PLD lakes overflown by SWOT
        
        :param in_filename: full path of the output file
        :type in_filename: string
        :param in_proc_metadata: processing metadata
        :type in_proc_metadata: dict
        """
        # log logger = logging.getLogger(self.__class__.__name__)
        # log logger.debug("- start -")
        
        # 1 - Write them in the output shapefile
        my_shp.write_mem_layer_as_shp(self.content_prior.layer, in_filename)
        
        # 2 - Estimate time_coverage_start and time_coverage_end
        self.content_prior.layer.ResetReading()
        time_utc_values = []
        for feature in self.content_prior.layer:
            time_utc_feat = feature.GetField("time")
            if time_utc_feat > 0.0:
                time_utc_values.append(time_utc_feat)
        self.content_prior.layer.ResetReading()
        time_str_dict = dict()
        if len(time_utc_values) > 0:
            time_str_dict["time_coverage_start"] = my_tools.convert_utc_to_str(min(time_utc_values), in_format=2)
            time_str_dict["time_coverage_end"] = my_tools.convert_utc_to_str(max(time_utc_values), in_format=2)
        else:
            time_str_dict["time_coverage_start"] = "None"
            time_str_dict["time_coverage_end"] = "None"
        
        # 3 - Write XML metadatafile for shapefile
        #logger.debug("Writing associated metadata file = %s.xml" % in_filename)
        if self.obj_pixc is None:
            self.content_prior.update_and_write_metadata("%s.xml" % in_filename, 
                                                         in_proc_metadata={**in_proc_metadata, **time_str_dict})
        else:
            self.content_prior.update_and_write_metadata("%s.xml" % in_filename, 
                                                     in_inprod_metadata=self.obj_pixc.pixc_metadata,
                                                     in_proc_metadata={**in_proc_metadata, **time_str_dict})
    
    def write_unknown_file(self, in_filename, in_proc_metadata):
        """
        Write the file containing water features unassigned to any prior features (ie neither PRD nor PLD)
        
        :param in_filename: full path of the output file
        :type in_filename: string
        :param in_proc_metadata: processing metadata
        :type in_proc_metadata: dict
        """
        # log logger = logging.getLogger(self.__class__.__name__)
        # log logger.debug("- start -")
        
        # 1 - Write them in the output shapefile
        my_shp.write_mem_layer_as_shp(self.content_unassigned.layer, in_filename)
            
        # 2 - Estimate time_coverage_start and time_coverage_end
        self.content_unassigned.layer.ResetReading()
        time_utc_values = []
        for feature in self.content_unassigned.layer:
            time_utc_feat = feature.GetField("time")
            if time_utc_feat > 0.0:
                time_utc_values.append(time_utc_feat)
        self.content_unassigned.layer.ResetReading()
        time_str_dict = dict()
        if len(time_utc_values) > 0:
            time_str_dict["time_coverage_start"] = my_tools.convert_utc_to_str(min(time_utc_values), in_format=2)
            time_str_dict["time_coverage_end"] = my_tools.convert_utc_to_str(max(time_utc_values), in_format=2)
        else:
            time_str_dict["time_coverage_start"] = "None"
            time_str_dict["time_coverage_end"] = "None"
        
        # 3 - Write XML metadatafile for shapefile
        #logger.debug("Writing associated metadata file = %s.xml" % in_filename)
        if self.obj_pixc is None:
            self.content_unassigned.update_and_write_metadata("%s.xml" % in_filename,
                                                              in_proc_metadata={**in_proc_metadata, **time_str_dict})
        else:
            self.content_unassigned.update_and_write_metadata("%s.xml" % in_filename, 
                                                              in_inprod_metadata=self.obj_pixc.pixc_metadata,
                                                              in_proc_metadata={**in_proc_metadata, **time_str_dict})
    

#######################################
        
        
def write_lakesp_file(in_list_tmp_lakesp_shp, in_list_laketile_shp, 
                      in_filename, in_proc_metadata, 
                      continent_id=None, flag_del_tmp_shp=True, flag_add_nogeom_feat=False):
    """
    Write combined LakeSP_Obs or _Prior or _Unassigned shapefile
    
    :param in_list_tmp_lakesp_shp: list of temporary LakeSP shapefiles to combine
    :type in_list_tmp_lakesp_shp: list
    :param in_list_laketile_shp: list of LakeTile shapefiles to combine
    :type in_list_laketile_shp: list
    :param in_filename: full path of the output file
    :type in_filename: string
    :param in_proc_metadata: processing metadata
    :type in_proc_metadata: dict
    :param continent_id: 2-letter identifier of the processed continent
    :type continent_id: string
    :param flag_del_tmp_shp: =True to delete temporary shapefiles (default), =False otherwise)
    :type flag_del_tmp_shp: boolean
    :param flag_add_nogeom_feat: =True to add features with no geometry in the output layer, =False otherwise (default)
    :type flag_add_nogeom_feat: boolean
    """
    logger = logging.getLogger("proc_lake")
    logger.debug("== write_lakesp_file ==")

    # 1 - Merge shapefiles LakeSP swath R and L and LakeTile_shp
    in_list_tmp_shp = in_list_tmp_lakesp_shp + in_list_laketile_shp
    my_shp.merge_shp(in_filename, in_list_tmp_shp, in_continent_id=continent_id, flag_add_nogeom_feat=flag_add_nogeom_feat)
    # Open data source and layer
    data_source_sp = my_shp.shp_driver.Open(in_filename, 1)  # Open in writing mode
    layer_sp = data_source_sp.GetLayer()  # Get the layer

    # 2 - Estimate time_coverage_start and time_coverage_end
    time_utc_values = []
    for feature in layer_sp:
        time_utc_feat = feature.GetField("time")
        if time_utc_feat > 0.0:
            time_utc_values.append(time_utc_feat)
    layer_sp.ResetReading()
    time_str_dict = dict()
    if len(time_utc_values) > 0:
        time_str_dict["time_coverage_start"] = my_tools.convert_utc_to_str(min(time_utc_values), in_format=2)
        time_str_dict["time_coverage_end"] = my_tools.convert_utc_to_str(max(time_utc_values), in_format=2)
    else:
        time_str_dict["time_coverage_start"] = "None"
        time_str_dict["time_coverage_end"] = "None"
    
    # 3 - Write XML metadatafile for shapefile
    logger.debug("Writing associated metadata file = %s.xml" % in_filename)
    # 3.1 - Build dictionary for metadata update
    updated_metadata = update_lakesp_metadata(in_list_tmp_lakesp_shp)
    # 3.2 - Init LakeSP product; type is retrieve from in_filename pattern
    if os.path.basename(in_filename).startswith(locnes_filenames.LAKE_SP_PREFIX["obs"]):
        tmp_lakesp = shp_file.LakeSPObsProduct(None, flag_create=False)
    elif os.path.basename(in_filename).startswith(locnes_filenames.LAKE_SP_PREFIX["prior"]):
        tmp_lakesp = shp_file.LakeSPPriorProduct(None, flag_create=False)
        # For _Prior file only: merge R and L features corresponding to the same PLD lake
        tmp_lakesp.merge_duplicate_features(layer_sp)
    elif os.path.basename(in_filename).startswith(locnes_filenames.LAKE_SP_PREFIX["unknown"]):
        tmp_lakesp = shp_file.LakeSPUnassignedProduct(None, flag_create=False)
    else:
        message = "Filename pattern is UNKNOWN"
        raise service_error.ProcessingError(message, logger)
    # 3.3 - Update and write metadatafile
    tmp_lakesp.update_and_write_metadata("%s.xml" % in_filename, 
                                          in_proc_metadata={**in_proc_metadata, **time_str_dict, **updated_metadata})
    
    # 4 - Delete temporary shapefiles if asked
    if flag_del_tmp_shp:
        shp_driver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Driver for shapefiles
        for cur_tmp_file in in_list_tmp_lakesp_shp:
            logger.debug("Delete temporary shapefile = %s" % os.path.basename(cur_tmp_file))
            shp_driver.DeleteDataSource(cur_tmp_file)
            os.remove("%s.xml" % cur_tmp_file)  # Also remove associated .shp.xml file
    else:
        logger.debug("Keep temporary swath LakeSP shapefiles")

    # 5 - Close merged shapefile
    data_source_sp.Destroy()
    
def update_lakesp_metadata(in_list_shp):
    """
    Compute matadata values wrt metadata from multiple .shp.xml files
    
    :param in_list_shp: list of input shapefiles
    :type in_list_shp: list of string
    
    :return: out_metadata = updated metadata
    :rtype: out_metadata = dict
    """
    logger = logging.getLogger("proc_lake")
    logger.debug("== update_lakesp_metadata ==")
    
    # 0 - Init variables
    # 0.1 - List of attributes to update
    list_attributes = ["time_granule_start", "time_granule_end", \
                       "geospatial_lon_min", "geospatial_lon_max", "geospatial_lat_min", "geospatial_lat_max", \
                       "left_first_longitude", "left_first_latitude", "left_last_longitude", "left_last_latitude", \
                       "right_first_longitude", "right_first_latitude", "right_last_longitude", "right_last_latitude"]
    # 0.2 - Output metadata
    out_metadata = dict()
    for key in list_attributes:
        if key.startswith("time"):
            out_metadata[key] = "None"
        else:
            out_metadata[key] = -9999.0
    
    for cur_shp in in_list_shp:
        
        # 1 - Load related .shp.xml file
        metadata = ET.parse(cur_shp + ".xml")
        
        # 2 - Update output metadata wrt currently loaded metadata
        for key in list_attributes:
            
            if key.startswith("time"):
                # Retrieve key value in current metadata
                cur_value = metadata.xpath("//swot_product/global_metadata/%s" % key)[0].text
                
                # Update output value for key 
                if (cur_value is not None) and (cur_value != "None"):
                    if out_metadata[key] == "None":
                        out_metadata[key] = cur_value
                    elif key == "time_coverage_start":
                        if cur_value < out_metadata[key]:
                            out_metadata[key] = cur_value
                    else:
                        if cur_value > out_metadata[key]:
                            out_metadata[key] = cur_value
                
            else:
                # Retrieve key value in current metadata
                cur_value = float(metadata.xpath("//swot_product/global_metadata/%s" % key)[0].text)
                
                # Update output value for key 
                if cur_value != -9999.0:
                    if out_metadata[key] == -9999.0:
                        out_metadata[key] = cur_value
                    elif ("min" in key) or ("first" in key):
                        if cur_value < out_metadata[key]:
                            out_metadata[key] = cur_value
                    else:
                        if cur_value > out_metadata[key]:
                            out_metadata[key] = cur_value
                        
    return out_metadata
    

#######################################

    
def select_valid_pixels(in_nb_pixels, in_lon, in_lat, in_height):
    """
    Test potential issue in output of height-constrained geolocation process,
    i.e. select valid PIXC = PIXC without NaNs for longitude or latitude or height
    
    :param in_nb_pixels: nb of PIXC
    :type in_nb_pixels: int
    :param in_lon: improved longitudes to test
    :type in_lon: 1D-array of float
    :param in_lat: improved latitudes to test
    :type in_lat: 1D-array of float
    :param in_height: improved heights to test
    :type in_height: 1D-array of float
    
    :return: out_valid_flag =0 if there is no valid PIXC; =2 if all PIXC are valid; =1 otherwise
    :rtype: out_valid_flag = int
    :return: out_valid_index = indices of valid PIXC when out_valid_flag == 1; =None otherwise
    :rtype: out_valid_index = 1D-array of float
    """
    logger = logging.getLogger("proc_lake.py")
    
    # 1 - Init output values
    out_valid_flag = 2
    out_valid_index = None
    
    # 2 - Retrieve indices of PIXC with either longitude, or latitude or heights being NaN
    tmp_nan_index1 = np.where(np.isnan(in_lon))[0]
    tmp_nan_index2 = np.where(np.isnan(in_lat))[0]
    tmp_nan_index3 = np.where(np.isnan(in_height))[0]
    tmp_nan_index = np.concatenate((tmp_nan_index1, tmp_nan_index2, tmp_nan_index3))
    tmp_nan_index_clean = np.unique(tmp_nan_index)
    nb_nan = len(tmp_nan_index_clean)
    
    # 3 - Compute valid indices
    if nb_nan == in_nb_pixels:
        logger.debug("WARNING !!! All the pixels have NaN for improved geolocation => feature not computed")
        out_valid_flag = 0
    elif nb_nan != 0:
        logger.debug("WARNING !!! %d pixels have NaN for improved geolocation => removed from computation" % nb_nan)
        out_valid_flag = 1
        tmp_valid_index1 = np.where(np.isfinite(in_lon))[0]
        tmp_valid_index2 = np.where(np.isfinite(in_lat))[0]
        tmp_valid_index3 = np.where(np.isfinite(in_height))[0]
        tmp_valid_index = np.concatenate((tmp_valid_index1, tmp_valid_index2, tmp_valid_index3))
        out_valid_index = np.unique(tmp_valid_index)
        
    return out_valid_flag, out_valid_index
    
