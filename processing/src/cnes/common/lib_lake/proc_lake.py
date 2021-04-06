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

import logging
import numpy as np
import os
from osgeo import ogr

import cnes.common.service_config_file as service_config_file
import cnes.common.service_error as service_error

from cnes.modules.geoloc.scripts.biglake_model import BigLakeModel

import cnes.common.lib.my_hull as my_hull
import cnes.common.lib.my_shp_file as my_shp
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_variables as my_var
import cnes.common.lib_lake.lake_db as lake_db
import cnes.common.lib_lake.locnes_products_shapefile as shp_file
import cnes.common.lib_lake.proc_pixc_vec as proc_pixc_vec


class LakeProduct(object):
    """
    class LakeProduct
    Manage LakeTile and LakeSP products
    """
    
    def __init__(self, in_product_type, in_obj_pixc, in_obj_pixc_vec, in_obj_lake_db, in_layer_name):
        """
        Constructor

        :param in_product_type: type of product among "SP"=LakeSP and "TILE"=LakeTile
        :type in_product_type: string
        :param in_obj_pixc: pixel cloud from which to compute lake products
        :type in_obj_pixc: proc_pixc.PixelCloud or proc_pixc_sp.PixelCloudSP object
        :param in_obj_pixc_vec: pixel cloud complementary file from which to compute lake products
        :type in_obj_pixc_vec: proc_pixc_vec.PixelCloudVec or proc_pixc_vec_sp.PixelCloudVecSP object
        :param in_obj_lake_db: lake database
        :type in_obj_lake_db: lake_db.lakeDb_shp or lake_db.lakeDb_sqlite
        :param in_layer_name: name for lake product layer
        :type in_layer_name: string

        Variables of the object:
            - cfg / service_config_file.cfg: instance of LOCNES configuration file
            - product_type / string: type of product among "SP"=LakeSP and "TILE"=LakeTile
            - obj_pixc / proc_pixc.PixelCloud or proc_pixc_sp.PixelCloud: pixel cloud from which to compute lake products
            - obj_pixc_vec / proc_pixc_vec.PixelCloudVec or proc_pixc_vec_sp.PixelCloudVec: extra info for pixel cloud
            - obj_lake_db / lake_db.lakeDb_shp or lake_db.lakeDb_sqlite: lake database
            - content_obs / LakeSPShpProduct: container of the lake "obs" product
            - content_prior / LakeSPShpProduct: container of the lake "prior" product
            - lakeid_uniq / set: list of uniq prior identifiers linked to all observed objects
            - compare_stats / dict: store parameters for global comparison between _Obs and _Prior lake products
            - compare_stats_param / list: list of parameters to compare between _Obs and _Prior lake products
        """
        # Get instance of service config file
        self.cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
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
        self.obj_pixc_vec = in_obj_pixc_vec
        # Lake database object
        self.obj_lake_db = in_obj_lake_db
        
        # 2 - Initialize lake product contents
        if self.product_type == "TILE":
            self.content_obs = shp_file.LakeTileObsProduct(in_layer_name)  # Obs-oriented file type
            self.content_prior = shp_file.LakeTilePriorProduct(in_layer_name)  # PLD-oriented file type
        else:
            self.content_obs = shp_file.LakeSPObsProduct(in_layer_name)  # Obs-oriented file type
            self.content_prior = shp_file.LakeSPPriorProduct(in_layer_name)  # PLD-oriented file type
        
        # 3 - Other variables
        self.lakeid_uniq = set()  # List of uniq prior identifiers linked to all observed objects
        # Dictionnary to store parameters for global comparison between _Obs and _Prior lake products
        self.compare_stats = {}
        self.compare_stats_params = ["delta_s_l", "delta_s_q", "area_total", "area_detct"]
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

    # ----------------------------------------

    def compute_lake_features(self, in_list_labels):
        """
        Computes lake features for pixels for which label is in in_list_labels.

         - NB: This processing is limited to water bodies being of a minimum size defined by MIN_SIZE.
         - NB2: Improved geolocation is computed for all entities, not only those > MIN_SIZE

        :param in_list_labels: list of labels to process
        :type in_list_labels: 1D-array of int
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 0 - Retrieve needed configuration parameters
        min_size = self.cfg.getfloat('CONFIG_PARAMS', 'MIN_SIZE')
        hull_method = self.cfg.getfloat("CONFIG_PARAMS", "HULL_METHOD")
        flag_geoloc = self.cfg.getboolean('CONFIG_PARAMS', 'IMP_GEOLOC')
        nb_digits = self.cfg.getint('CONFIG_PARAMS', 'NB_DIGITS')
        
        # Print min_size
        logger.info("Minimum size for lakes = %0.2f km2" % min_size)
        # Print geoloc
        if flag_geoloc:
            logger.info("Improve geoloc = YES")
        else:
            logger.info("Improve geoloc = NO")
        # Print hull_method
        # TODO : les égalités de réel ne sont pas autorisé il faut changer ces tests
        if hull_method == 0:
            logger.info("Use of CONVEX HULL for lake boundaries")
        elif hull_method == 1.0:
            logger.info("Use of CONCAVE HULL based on Delaunay triangulation with CGAL library for lake boundaries")
        elif hull_method == 1.1:
            logger.info("Use of CONCAVE HULL based on Delaunay triangulation")
        elif hull_method == 1.2:
            logger.info("Use of CONCAVE HULL based on Delaunay triangulation with varying alpha param for lake boundaries")
        elif hull_method == 2:
            logger.info("Use of RADAR VECTORISATION for lake boundaries")
        else:
            message = "HULL_METHOD values unkown (%d); should be 0(convex) 1(concave) 2(radar vect)" % hull_method
            raise service_error.ProcessingError(message, logger)

        ########################################################
        # Compute geometry and attributes of observed features #
        # except those related to PLD (ie lake_id, overlap,    #
        # p_ attributes and storage change)                    #
        ########################################################
        
        # Init variables
        cpt_too_small = 0  # Counter of too small objects
        cpt_obj = 1  # Counter of processed objects
        
        for indl, label in enumerate(in_list_labels):  # Loop on inside tile objects

            # ======================================================
            # 1 - Get pixels indices for associated to current label
            # ======================================================
            pix_index = np.where(self.obj_pixc.labels == label)[0]
            obj_nb_pix = pix_index.size
            if obj_nb_pix == 0:
                logger.warning("[STRANGE...] label %s corresponds to 0 pixel..." % label)
                continue
            
            # ===============================================
            # 2 - Compute categories wrt classification flags
            # ===============================================
            classif = self.sort_pixels_wrt_classif_flags(pix_index)

            # ========================
            # 3 - Compute object sizes
            # ========================
            
            # 3.1 - Total area (=water + dark water)
            obj_area_total = np.nansum(self.obj_pixc.pixel_area[pix_index[select_water_dark_pixels(classif, in_flag_water=True, in_flag_dark=True)]]) 
            # In m2
            obj_area_total /= 10**6  # Conversion in km2
            
            # 3.2 - Detected area (=only water)
            obj_area_detected = np.nansum(self.obj_pixc.pixel_area[pix_index[select_water_dark_pixels(classif, in_flag_water=True, in_flag_dark=False)]])
            # In m2
            obj_area_detected /= 10**6  # Conversion in km2
            
            logger.info("")
            logger.info("===== compute_product %d over %d / label = %d / nb pixels = %d / detected area (total area) = %.2f km2 (%.2f km2) =====" \
                        % (indl+1, len(in_list_labels), label, obj_nb_pix, obj_area_detected, obj_area_total))
            
            # ================================
            # 4 - Compute improved geolocation
            # ================================
            
            # 4.1 - Compute improved geolocation
            imp_lon, imp_lat, imp_height = self.improve_geoloc(pix_index, classif, obj_area_total)
            
            # 4.2 - Update PIXCVec with improved geolocation infos
            self.update_pixcvec_with_imp_geoloc(pix_index, imp_lon, imp_lat, imp_height)
            
            # 4.3 - Select valid PIXC wrt value of longitude, latitude or height
            #       which must be finite
            valid_flag, not_nan_index = select_valid_pixels(obj_nb_pix, imp_lon, imp_lat, imp_height)
            # No need to continue if there is no valid PIXC
            if valid_flag == 0:  
                #self.add_pld_features_not_observed()  # Add not observed PLD lakes to _Prior layer
                continue

            # ====================================================
            # 5 - Compute lake feature 
            #     if detected water area of object is large enough 
            #     and all pixels are not aligned
            # ====================================================
            is_aligned = my_tools.are_pixels_aligned(self.obj_pixc.range_index[pix_index], self.obj_pixc.azimuth_index[pix_index])

            if obj_area_detected >= min_size and not is_aligned:

                # 5.1 - Compute obs number
                if self.product_type == "TILE":  # "TILE" case: only add cpt_obj
                    obs_number = str(cpt_obj).rjust(nb_digits, str('0'))
                else :  # "SP" case: add main tile info
                    obs_number = str(self.obj_pixc.get_lake_tile_label(label)).rjust(nb_digits, str('0'))
                    
                # 5.2 - Prepare needed arrays of valid PIXC
                if valid_flag == 2:
                    valid_pix_index = pix_index
                    valid_classif = classif
                    valid_lon = imp_lon
                    valid_lat = imp_lat
                else:
                    valid_pix_index = pix_index[not_nan_index]
                    valid_classif = self.sort_pixels_wrt_classif_flags(pix_index[not_nan_index])
                    valid_lon = imp_lon[not_nan_index]
                    valid_lat = imp_lat[not_nan_index]
                    
                # 5.3 - Add observed feature
                obs_id, pixcvec_lakeid = self.add_obs_feature(obs_number, valid_pix_index, valid_classif, valid_lon, valid_lat)
                
                # 5.4 - Update PIXCVec obs_id and lake_id attributes
                if obs_id is not None:
                    self.update_pixcvec_with_ids(valid_pix_index, obs_id, pixcvec_lakeid)

                # 5.5 - Increase counter of processed objects
                cpt_obj += 1

            else:
                logger.info("> Detected area of object %d is too small (= %.2f km2; total area = %.2f km2 \
                             with %d pixels)"%(label, obj_area_detected, obj_area_total, obj_nb_pix))
                cpt_too_small += 1  # Increase counter of too small objects

        logger.info("> %d objects not processed because too small" % cpt_too_small)

        ##################################################################
        # Compute _p attributes and storage change for observed features #
        # + compute geometry and attributes of PLD features              #
        ##################################################################
        
        nb_prior = len(self.lakeid_uniq)
        
        if nb_prior == 0:
            logger.info("NO observed feature linked to a PLD lake")
            
        else:
            # Deal with PLD lakes observed by SWOT
            logger.info("==================")
            logger.info("Compute PLD features")
            logger.info("==================")
            logger.info("")

            logger.info("%d PLD lakes linked to observed lakes" % nb_prior)

            for cur_lakeid in self.lakeid_uniq:
                logger.info("===== Deal with PLD lake %s =====" %cur_lakeid)

                # 6.0 - Remove from the list of PLD located over the area
                try:
                    self.obj_lake_db.list_lakeid.remove(cur_lakeid)
                except:  # To handle obs lakes at the edge of the tile, intersecting PLD lakes not covered by the tile
                    msg = "PLD lake with lake_id=%s is not covered by the PIXC tile BUT intersects an observed feature" % cur_lakeid
                    logger.warning(msg)
        
                # 6.1.1 - Create prior lake object
                obj_plake = lake_db.PriorLake(self.obj_lake_db, cur_lakeid)
                # 6.1.2 - Format PLD attributes (i.e. convert None to required _FillValue if needed)
                pld_attributes = obj_plake.format_attributes()
                
                # 6.2 - Update p_ attributes of observed features strongly connected to this PLD lake
                # 6.2.1 - Select them
                self.content_obs.layer.SetAttributeFilter("lake_id LIKE '{}%'".format(cur_lakeid))
                nb_obslake = self.content_obs.layer.GetFeatureCount()
                logger.debug("{} observed lake(s) are strongly connected to this PLD lake".format(nb_obslake))
                # 6.2.2 - Set p_ attributes from PLD infos to all observed lakes having this PLD lake as main overlap
                if nb_obslake > 0:
                    self.set_pld_attributes(pld_attributes)
                # 6.2.3 - Reinit layer attribute filter
                self.content_obs.layer.SetAttributeFilter(None)
                
                # 6.3 - Select all observed lakes overlapping the PLD lake
                self.content_obs.layer.SetAttributeFilter("lake_id LIKE '%{}%'".format(cur_lakeid))
                nb_obslake = self.content_obs.layer.GetFeatureCount()
                logger.debug("{} observed lake(s) are connected to this PLD lake".format(nb_obslake))
                
                if nb_obslake == 0:
                    logger.warning("[STRANGE...] PLD lake listed in _Obs file corresponds to NO obs lake...")
                    
                    # 6.6 - Add prior feature to _Prior layer
                    self.content_prior.add_feature(None, pld_attributes)
                    
                else:
        
                    # 6.4 - Retrieve PIXCVec indices corresponding to prior feature
                    # NB: use of selected_index to convert PIXCVec indices to PIXC indices reference
                    if self.product_type == "SP":
                        pixc_index = np.where(self.obj_pixc_vec.lake_id == obj_plake.lake_id.encode())[0]
                    else:
                        pixc_index = np.where(self.obj_pixc_vec.lake_id[self.obj_pixc.selected_index] == obj_plake.lake_id.encode())[0]
                
                    # 6.5 - Compute observed geometry and common attributes of PLD feature
                    prior_geom, prior_attributes = self.form_prior_feature(obj_plake, pixc_index)
                    
                    if obj_plake.ok_to_compute_stocc:
                        # 6.6 - Compute storage change for this PLD lake and all observed lakes overlapping it
                        storage_change_values = self.compute_storage_change(obj_plake, pixc_index, prior_attributes)                     
                        # 6.7 - Add prior feature to _Prior layer
                        self.content_prior.add_feature(prior_geom, {**prior_attributes, **pld_attributes, **storage_change_values})
                    else:
                        logger.debug("Not able to compute storage change")
                        # 6.7 - Add prior feature to _Prior layer
                        self.content_prior.add_feature(prior_geom, {**prior_attributes, **pld_attributes})
                
                # 6.8 - Reinit layer attribute filter
                self.content_obs.layer.SetAttributeFilter(None)
                
        # 6.8 - Deal with PLD lakes which should have been observed by SWOT
        if self.product_type == "TILE":
            self.add_pld_features_not_observed()

    # ----------------------------------------
    # Fonctions specific to height-constrained 
    # geolocation process 
    # ----------------------------------------
    
    def improve_geoloc(self, in_pix_index, in_classif, in_area_total):
        """
        Prepare and run the height-constrained geolocation algorithm over the feature
            defined by the PIXC of indices definied in input.
        
        :param in_pix_index: indices of PIXC defining the feature
        :type in_pix_index: 1D-array of int
        :param in_classif: classification of input PIXC
        :type in_classif: 1D-array of int
        :param in_area_total: total area (ie water + dark) of input PIXC
        :type in_area_total: float
        
        :return: out_lon = improved longitudes
        :rtype: out_lon = 1D-array of float
        :return: out_lat = improved latitudes
        :rtype: out_lat = 1D-array of float
        :return: out_height = improved heights
        :rtype: out_height = 1D-array of float
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # 1 - Compute mean height using both water and dark water pixel (weighting using uncertainties)
        # Use of external method common with RiverObs to compute mean/std
        mean_height = self.obj_pixc.compute_height(in_pix_index)
        
        # 2 - Improve geolocation
        if self.cfg.getboolean('CONFIG_PARAMS', 'IMP_GEOLOC'):

            # 2a - Fit lake height model depending on lake size
            height_model = self.compute_height_model(in_pix_index, in_classif, mean_height, in_area_total)

            # 2b - Compute imp geolocation 
            out_lon, out_lat, out_height, p_final = proc_pixc_vec.compute_imp_geoloc(self.product_type, self.obj_pixc, in_pix_index, height_model)
            
            # 2c - Compute flattened interferogram in order to compute coherence during uncertainties estimation
            self.obj_pixc.compute_interferogram_flattened(in_pix_index, p_final)

        else:
            out_lon = self.obj_pixc.longitude[in_pix_index]
            out_lat = self.obj_pixc.latitude[in_pix_index]
            out_height = self.obj_pixc.height[in_pix_index]
                
        return out_lon, out_lat, out_height
    
    def compute_height_model(self, in_pix_index, in_classif, in_mean_height, in_area_total):
        """
        Compute the height model related to the feature defined by the PIXC
        of indices defined in input. This height model depends on the feature size.
        
        :param in_pix_index: indices of PIXC defining the feature
        :type in_pix_index: 1D-array of int
        :param in_classif: classification of input PIXC
        :type in_classif: 1D-array of int
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

            logger.info("Using {} biglake model for improved geolocation (lake total area = {} km2)".format(biglake_model, in_area_total))

            if biglake_model == 'grid':
                out_height_model = biglakemodel.fit_biglake_model(self.obj_pixc,
                                                                  in_pix_index,
                                                                  grid_spacing=biglake_grid_spacing,
                                                                  grid_resolution=biglake_grid_res)

            elif biglake_model == 'polynomial':                                 
                out_height_model = biglakemodel.fit_biglake_model_polyfit(self.obj_pixc, in_pix_index, in_classif)
                                                 
            else:
                logger.debug("No height model defined, assume Mean Height model")
                out_height_model = np.full(self.obj_pixc.height[in_pix_index].shape, in_mean_height)

        else:
            logger.info("Using lake average height = {} m for improved geolocation (lake total area = {} km2)".format(in_mean_height, in_area_total))
            out_height_model = np.full(self.obj_pixc.height[in_pix_index].shape, in_mean_height)
            
        # Return
        return out_height_model

    # ----------------------------------------
    # Functions specific to lake feature (i.e.
    # geometry + attributes) computation
    # ----------------------------------------
    
    def add_obs_feature(self, in_number, in_indices, in_classif, in_lon, in_lat):
        """
        Process valid PIXC related to current feature
        to build the observed feature boundary and associated attributes
        and add the observed feature to the obs-oriented layer
        
        :param in_number: number of the feature in the scene
        :type in_number: str
        :param in_indices: list of indices of the PIXC related to the feature
        :type in_indices: 1D-array of int
        :param in_classif_dict: dictionary of indices of pixels of in_indices corresponding to categories "water" and "dark"
        :type in_classif_dict: dict (output of self.sort_pixels_wrt_classif_flags)
        :param in_lon: improved longitudes vector for PIXC of the feature
        :type in_lon: 1D-array of float
        :param in_lat: improved latitudes vector for PIXC of the feature
        :type in_lat: 1D-array of float
        
        :return: out_obs_id = obs_id identifier of the feature
        :rtype: out_obs_id = string
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Deal with object number = {}".format(in_number))
        
        # 1 - Build feature boundary
        feature_geom = self.build_obs_boundary(in_indices, in_lon, in_lat)
        
        # 2 - Compute common attributes
        feature_attributes = self.compute_common_attributes(in_indices, in_classif_dict=in_classif)
        
        # 3 - Link to PLD
        lake_id, overlap, out_pixcvec_lakeid = self.link_obs_geom_to_pld_geom(feature_geom, in_lon, in_lat)
        feature_attributes["lake_id"] = lake_id
        feature_attributes["overlap"] = overlap
        
        # 4.1 - Compute 3 digits of basin identifier (CBB)
        if lake_id != "no_data":
            lake_basin_id = lake_id[0:3]
        else:
            lake_basin_id = str(self.obj_lake_db.link_poly_to_basin(feature_geom)[0]).rjust(3, str('0'))
        # 4.2 - Form obs_id
        if self.product_type == "TILE":  # "TILE" case: only add cpt_obj
            tile_ref = str(self.obj_pixc.pixc_metadata["tile_number"]).rjust(3, str('0')) + str(self.obj_pixc.pixc_metadata["swath_side"])
        else:  # "SP" case: add main tile info
            tile_ref = self.obj_pixc.get_majority_pixels_tile_ref(int(in_number))
        feature_attributes["obs_id"] = "%s%s%s" % (lake_basin_id, tile_ref, in_number)
        
        # 5 - Add feature to layer, if it exists
        logger.debug("obs_id = %s / lake_id = %s" % (str(feature_attributes["obs_id"]), lake_id))
        if feature_geom is not None:
            self.content_obs.add_feature(feature_geom, feature_attributes)
            out_obs_id = feature_attributes["obs_id"]
        else:
            logger.warning("Feature not added to layer because geometry is None")
            out_obs_id = None
            
        return out_obs_id, out_pixcvec_lakeid
    
    def build_obs_boundary(self, in_indices, in_lon, in_lat):
        """
        Build observed feature boundary from PIXC defined by the input coordinates
        
        :param in_indices: list of indices of the PIXC related to the feature
        :type in_indices: 1D-array of int
        :param in_lon: improved longitudes vector for PIXC of the feature
        :type in_lon: 1D-array of float
        :param in_lat: improved latitudes vector for PIXC of the feature
        :type in_lat: 1D-array of float
        
        :return: out_geom = boundary of the feature
        :rtype: out_geom = OGRPolygon
        """
        
        if self.product_type == 'SP':
            out_geom = my_hull.compute_lake_boundaries(in_lon,
                                                       in_lat,
                                                       self.obj_pixc.get_range_of_lake(in_indices),
                                                       self.obj_pixc.get_azimuth_of_lake(in_indices),
                                                       self.obj_pixc.nb_pix_range)
            
        else:
            out_geom = my_hull.compute_lake_boundaries(in_lon,
                                                       in_lat,
                                                       self.obj_pixc.range_index[in_indices],
                                                       self.obj_pixc.azimuth_index[in_indices],
                                                       self.obj_pixc.nb_pix_range)
            
        return out_geom
    
    def form_prior_feature(self, in_obj_plake, in_pixc_index):
        """
        Create and initialize prior feature 
        
        :param in_obj_plake: PLD lake object
        :type in_obj_plake: lake_db.PriorLake
        :param in_pixc_index: indices of PIXC related to PLD lake
        :type in_pixc_index: Numpy 1D-array of int
        
        :return: out_geom = geometry of prior feature
        :rtype: out_geom = OGRPolygon
        :return: out_attributes = attributes of prior feature
        :rtype: out_attributes = dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Deal with PLD lake = {}".format(in_obj_plake.lake_id))
        
        # 2 - Compute common attributes
        out_attributes = self.compute_common_attributes(in_pixc_index)
        # Save for comparison
        for param in self.compare_stats_params:
            if (param in out_attributes.keys()) and (out_attributes[param] is not None):
                self.compare_stats["prior"][param] += out_attributes[param]
        
        # 3 - Build feature boundary
        out_geom, list_obs_id, list_overlap = self.build_prior_boundary(in_obj_plake, in_pixc_index)
        
        # 4 - Add identifiers and overlap values
        out_attributes["lake_id"] = in_obj_plake.lake_id
        if list_obs_id:
            out_attributes["obs_id"] = ';'.join(list_obs_id)
        else:
            out_attributes["obs_id"] = "no_data"
        if list_obs_id:
            out_attributes["overlap"] = ';'.join(list_overlap)
        else:
            out_attributes["overlap"] = "no_data"
        
        return out_geom, out_attributes
    
    def build_prior_boundary(self, in_obj_plake, in_pixc_index):
        """
        Build prior feature boundary by intersecting the influence area of the PLD lake
        with the observed features selected in _Obs layer.
        
        :param in_obj_plake: PLD lake object
        :type in_obj_plake: lake_db.PriorLake
        :param in_pixc_index: list of indices of the PIXC related to the PLD feature
        :type in_pixc_index: 1D-array of int
        
        :return: out_geom = geometry of prior feature
        :rtype: out_geom = OGRPolygon
        :return: out_list_obs_id = list of observed features intersecting the PLD lake
        :rtype: out_list_obs_id = list
        :return: out_list_overlap = list of fractions of PLD feature covered by each observed lake
        :rtype: out_list_overlap = list
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Deal with PLD lake = {}".format(in_obj_plake.lake_id))
        
        out_geom = None
        out_list_obs_id = []
        out_list_overlap = []
        
        # Case with 1 obs <-> 1 or N PLD lake(s)
        nb_obs_inter = self.content_obs.layer.GetFeatureCount()
        if nb_obs_inter == 1:
            
            # Retrieve associated obs feature
            cur_feature = self.content_obs.layer.GetNextFeature()
            cur_geom = cur_feature.GetGeometryRef()
            out_list_obs_id.append(cur_feature.GetField("obs_id"))
            logger.debug("obs_id = %s / lake_id = %s" % (out_list_obs_id[0], cur_feature.GetField("lake_id")))
            
            # Simple case with 1 obs <-> 1 PLD lake
            if ";" not in cur_feature.GetField("lake_id"):
                logger.debug("Simple case with 1 obs <-> 1 PLD lake => keep OBS geometry")
                tmp_geom = cur_geom
            
            # Complex case with 1 obs <-> N PLD lakes
            # PLD lake geometry is part of the obs geometry
            else:
                logger.debug("Complex case with 1 obs <-> N PLD lakes => keep part of OBS geometry")
                
                # Get the influence area polygon
                influence_area_poly = self.obj_lake_db.get_influence_area_poly(in_obj_plake.lake_id)
                
                if influence_area_poly is None:
                    tmp_pixcvec_index = self.get_pixcvec_index_from_pixc_index(in_pixc_index)
                    tmp_geom = self.build_obs_boundary(in_pixc_index, 
                                                       self.obj_pixc_vec.longitude_vectorproc[tmp_pixcvec_index],
                                                       self.obj_pixc_vec.latitude_vectorproc[tmp_pixcvec_index])
                    
                else:
                    tmp_geom = cur_geom.Intersection(influence_area_poly)
        
            # Compute overlaping area
            area_pld = my_tools.get_area(in_obj_plake.geom)
            geom_inter = tmp_geom.Intersection(in_obj_plake.geom)
            if geom_inter is not None:
                area_inter = my_tools.get_area(geom_inter)
                out_list_overlap.append(str(round(area_inter/area_pld*100.)))
                logger.debug("PLD lake area = {} m2 - PLD/obs intersection area = {} m2 - overlap = {}%".format(area_pld, area_inter, out_list_overlap[0]))
        
        # Case with N obs <-> 1 PLD lake
        else:
            
            logger.debug("%d obs features correspond to PLD lake" % nb_obs_inter)
            
            # Init output geometry
            tmp_geom = ogr.Geometry(ogr.wkbMultiPolygon)
            
            # Get the influence area polygon
            influence_area_poly = self.obj_lake_db.get_influence_area_poly(in_obj_plake.lake_id)
            
            tmp_list_obs_id = []
            tmp_list_overlap = []
            
            for cur_feature in self.content_obs.layer:
                
                cur_geom = cur_feature.GetGeometryRef()
                
                # Retrieve corresponding obs_id:
                cur_obs_id = cur_feature.GetField("obs_id")
                tmp_list_obs_id.append(cur_obs_id)
                logger.debug("Deal with obsid %s " % cur_obs_id)
                # Simple case with N obs <-> 1 PLD lake
                if ";" not in cur_feature.GetField("lake_id"):
                    logger.debug("Simple case with obs <-> 1 PLD lake => keep OBS geometry")
                    obs_poly = cur_geom.Clone()
                    
                # Complex case with N obs <-> N PLD lakes
                else:
                    logger.debug("Complex case with obs <-> N PLD lakes => keep part of OBS geometry")
                
                    if influence_area_poly is None:
                        in_pixcvec_index = self.get_pixcvec_index_from_pixc_index(in_pixc_index)
                        obsid_over_in_lakeid = self.obj_pixc_vec.obs_id[in_pixcvec_index]
                        tmp_pixc_index = in_pixc_index[np.where(obsid_over_in_lakeid == cur_obs_id.encode())]
                        tmp_pixcvec_index = self.get_pixcvec_index_from_pixc_index(tmp_pixc_index)
                        obs_poly = self.build_obs_boundary(tmp_pixc_index, 
                                                           self.obj_pixc_vec.longitude_vectorproc[tmp_pixcvec_index],
                                                           self.obj_pixc_vec.latitude_vectorproc[tmp_pixcvec_index])
                        
                    else:
                        obs_poly = cur_geom.Intersection(influence_area_poly)
                
                # Compute overlaping area
                area_pld = my_tools.get_area(in_obj_plake.geom)
                geom_inter = obs_poly.Intersection(in_obj_plake.geom)
                if geom_inter is not None:
                    area_inter = my_tools.get_area(geom_inter)
                    tmp_overlap = str(round(area_inter/area_pld*100.))
                    tmp_list_overlap.append(tmp_overlap)
                    logger.debug("PLD lake area = {} m2 - PLD/obs intersection area = {} m2 - overlap = {}%".format(area_pld, area_inter, tmp_overlap))
                else :
                    logger.warning("PLD lakeid %s geometrie do not intersects observed geometry %s " %(str(in_obj_plake.lake_id), str(tmp_list_obs_id)))
                # Add current geometry to output geometry
                tmp_geom = tmp_geom.Union(obs_poly)

            # Sort obs_id and overlap fractions by decreasing area intersection
            sorted_idx = sorted(range(len(tmp_list_overlap)), key=lambda k: tmp_list_overlap[k], reverse=True)
            out_list_obs_id = [tmp_list_obs_id[idx] for idx in sorted_idx]
            out_list_overlap = [tmp_list_overlap[idx] for idx in sorted_idx]

        if tmp_geom is not None:
            out_geom = tmp_geom.Clone()
        
        return out_geom, out_list_obs_id, out_list_overlap
                
    def compute_common_attributes(self, in_pixc_index, in_classif_dict=None):
        """
        Computes common attributes from PIXC related to the current feature.
        This subset of pixel cloud include pixels for which self.obj_pixc.labels=in_label

        :param in_pixc_index: list of indices of PIXC defining the current feature
        :type in_pixc_index: 1D-array of int
        :param in_classif_dict: dictionary of indices of pixels of in_indices corresponding to categories "water" and "dark"
        :type in_classif_dict: dict (output of self.sort_pixels_wrt_classif_flags)

        :return: out_attributes = lake attributes (None for each item which could not be computed)
        :rtype: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # Init output dictionary
        out_attributes = dict()

        # ==================================
        # 1 - Median datetime of observation
        # ==================================
        
        tmp_time = np.nanmedian(self.obj_pixc.nadir_time[in_pixc_index])
        out_attributes["time"] = value_or_none(tmp_time)  # UTC time
        out_attributes["time_tai"] = value_or_none(np.nanmedian(self.obj_pixc.nadir_time_tai[in_pixc_index]))  # TAI time
        if out_attributes["time"] is None:
            out_attributes["time_str"] = None
        else:
            out_attributes["time_str"] = my_tools.convert_utc_to_str(out_attributes["time"])  # Time in UTC as a string
        
        # =================================
        # 2 - Measured hydrology parameters
        # =================================
        
        # Area and uncertainty
        area_total, area_total_unc, area_detct, area_detct_unc = self.obj_pixc.compute_area_with_uncertainties(in_pixc_index, method='composite')
        out_attributes["area_total"] = value_or_none(area_total / 10**6)
        out_attributes["area_tot_u"] = value_or_none(area_total_unc / 10**6)
        out_attributes["area_detct"] = value_or_none(area_detct / 10**6)
        out_attributes["area_det_u"] = value_or_none(area_detct_unc / 10**6)
        
        # Water surface elevation and uncertainty
        mean_height, height_std, height_unc = self.obj_pixc.compute_height_with_uncertainties(in_pixc_index)
        out_attributes["wse"] = value_or_none(mean_height)
        out_attributes["wse_r_u"] = value_or_none(height_unc)
        out_attributes["wse_u"] = out_attributes["wse_r_u"]
        out_attributes["wse_std"] = value_or_none(height_std)
        
        # Metric of layover effect = layover area
        out_attributes["layovr_val"] = None
        
        # Median distance from PIXC to the satellite ground track
        out_attributes["xtrk_dist"] = np.nanmedian(self.obj_pixc.cross_track[in_pixc_index])
        
        # ======================
        # 3 - Quality indicators
        # ======================
        
        # Summary quality indicator
        out_attributes["quality_f"] = None
        
        # Fractional area of dark water
        out_attributes["dark_frac"] = (area_total - area_detct) / area_total * 100.0
        if in_classif_dict is not None:
            if in_classif_dict["dark"] is None:
                out_attributes["dark_frac"] = 0.0
            
        # Ice cover flags
        out_attributes["ice_clim_f"] = None
        out_attributes["ice_dyn_f"] = None
        
        # Partial flag: =1 if the lake is partially covered by the swath, 0 otherwise
        range_index = self.obj_pixc.range_index[in_pixc_index]
        nr_edge_pix = np.where(range_index == 0)[0]
        fr_edge_pix = np.where(range_index == self.obj_pixc.nb_pix_range-1)[0]
        # Specific processing for 1st and last tiles
        az_edge_pix = np.array([])
        if self.product_type == "SP":
            # Get pixels that belong to the first or last azimuth line
            az_edge_pix = np.where(self.obj_pixc.is_boundary_pix[in_pixc_index])[0]
        if nr_edge_pix.size+fr_edge_pix.size+az_edge_pix.size == 0:
            out_attributes["partial_f"] = 0
        else:
            out_attributes["partial_f"] = 1
        
        # ==========================
        # 4 - Geophysical references
        # ==========================
        
        # Geoid model height
        out_attributes["geoid_hght"] = value_or_none(my_tools.compute_mean_2sigma(self.obj_pixc.geoid[in_pixc_index]))
        # Earth tide
        out_attributes["solid_tide"] = value_or_none(my_tools.compute_mean_2sigma(self.obj_pixc.solid_earth_tide[in_pixc_index]))
        # Pole tide
        out_attributes["pole_tide"] = value_or_none(my_tools.compute_mean_2sigma(self.obj_pixc.pole_tide[in_pixc_index]))
        # Load tide
        out_attributes["load_tidef"] = value_or_none(my_tools.compute_mean_2sigma(self.obj_pixc.load_tide_fes[in_pixc_index]))
        out_attributes["load_tideg"] = value_or_none(my_tools.compute_mean_2sigma(self.obj_pixc.load_tide_got[in_pixc_index]))
        
        # =================================
        # 5 - Geophysical range corrections
        # =================================
        
        # Dry tropo corr
        out_attributes["dry_trop_c"] = value_or_none(my_tools.compute_mean_2sigma(self.obj_pixc.model_dry_tropo_cor[in_pixc_index]))
        # Wet tropo corr
        out_attributes["wet_trop_c"] = value_or_none(my_tools.compute_mean_2sigma(self.obj_pixc.model_wet_tropo_cor[in_pixc_index]))
        # Iono corr
        out_attributes["iono_c"] = value_or_none(my_tools.compute_mean_2sigma(self.obj_pixc.iono_cor_gim_ka[in_pixc_index]))
        
        # ==========================
        # 6 - Instrument corrections
        # ==========================
        
        # KaRIn correction from crossover cal processing evaluated for lake
        out_attributes["xovr_cal_c"] = value_or_none(my_tools.compute_mean_2sigma(self.obj_pixc.height_cor_xover[in_pixc_index]))
        
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
        # Save IDs of intersecting PLD lake 
        if list_prior:
            for cur_lake_id in list_prior:
                self.lakeid_uniq.add(cur_lake_id)

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
        Set p_ attributes to prior values to current PLD lake 
        and for all obs features linked to current PLD lake (i.e. currently selected in the _Obs layer)
        
        :param in_pld_infos: values of p_ attributes to populate
        :type in_pld_infos: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        for obs_lake in self.content_obs.layer:
            # Set attributes related to PLD lake
            for key, value in in_pld_infos.items():
                obs_lake.SetField(str(key), value)
            # Rewrite obs feature with updated attributes
            self.content_obs.layer.SetFeature(obs_lake)
    
    def add_pld_features_not_observed(self):
        """
        Add PLD lakes which were not observed to _Prior layer
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        nb_missing = len(self.obj_lake_db.list_lakeid)
        
        if nb_missing == 0:
            logger.info("ALL PLD lakes have been observed")
            
        else:
            logger.info("%d PLD lakes have NOT been observed" % nb_missing)
            
            for cur_lakeid in self.obj_lake_db.list_lakeid:
                
                # 1.1 - Create prior lake object
                obj_plake = lake_db.PriorLake(self.obj_lake_db, cur_lakeid)
                # 1.2 - Format PLD attributes (i.e. convert None to required _FillValue if needed)
                pld_attributes = obj_plake.format_attributes()
                
                # 2 - Add lake_id
                pld_attributes["lake_id"] = cur_lakeid
                
                # 3 - Add prior feature to _Prior layer
                self.content_prior.add_feature(None, pld_attributes)

    # ------------------------------------------------
    # Functions specific to storage change computation
    # ------------------------------------------------

    def compute_storage_change(self, in_obj_plake, in_pixc_index, in_prior_attributes):
        """
        Set storage change values for all observed features linked to current PLD lake
        and return them
        NB: these observed objects have been selected in the _Obs layer before the use of this function

        Envisionned cases:
            - 1 prior lake <=> 1 observed lake
            - 1 prior lake <=> 2 or more observed lakes
            - 1 observed lake <=> 2 or more prior lakes
            - mixed lakes...
        
        :param in_obj_plake: PLD lake object
        :type in_obj_plake: lake_db.PriorLake
        :param in_pixc_index: indices of PIXC related to PLD lake
        :type in_pixc_index: Numpy 1D-array of int
        :return: in_prior_attributes = attributes of prior feature related to PLD lake
        :rtype: in_prior_attributes = dict
        
        :return: out_storage_values = storage change values related to current PLD lake
        :rtype: out_storage_values = dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 0 - Init output
        out_storage_values = dict()
        out_storage_values["delta_s_l"] = my_var.FV_REAL
        out_storage_values["ds_l_u"] = my_var.FV_REAL
        out_storage_values["delta_s_q"] = my_var.FV_REAL
        out_storage_values["ds_q_u"] = my_var.FV_REAL
        
        # Nb of observed lakes related to in_obj_plake
        nb_obs_lake = self.content_obs.layer.GetFeatureCount()
        logger.debug("Case PLD lake -> %d obs lakes" % nb_obs_lake)
        
        # 1 - Build dictionary going in intput of storage change function
        
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
                    list_obs[obs_id]["area"], list_obs[obs_id]["area_u"], tmp1, tmp2 = self.obj_pixc.compute_area_with_uncertainties(in_pixc_index, flag_all=False)
                    list_obs[obs_id]["wse"], list_obs[obs_id]["wse_u"], tmp1 = self.obj_pixc.compute_height_with_uncertainties(in_pixc_index)
                    
                # Retrieve storage change previous values
                list_obs[obs_id]["lake_id"] = obs_lake.GetField(str("lake_id")).split(";")
                list_obs[obs_id]["overlap"] = np.array(obs_lake.GetField(str("overlap")).split(";"), dtype="float")
                list_obs[obs_id]["delta_s_l"] = float(obs_lake.GetField(str("delta_s_l")))
                list_obs[obs_id]["ds_l_u"] = float(obs_lake.GetField(str("ds_l_u")))
                list_obs[obs_id]["delta_s_q"] = float(obs_lake.GetField(str("delta_s_q")))
                list_obs[obs_id]["ds_q_u"] = float(obs_lake.GetField(str("ds_q_u")))
                if list_obs[obs_id]["delta_s_l"] < -9e11:
                    list_obs[obs_id]["delta_s_l"] = 0.0
                    list_obs[obs_id]["ds_l_u"] = 0.0
                    list_obs[obs_id]["delta_s_q"] = 0.0
                    list_obs[obs_id]["ds_q_u"] = 0.0
                
            else:
                # Case the observed lake is related to only 1 prior lake
                logger.debug("> sub-case obs lake -> 1 PLD lake")
                list_obs[obs_id]["area"] = float(obs_lake.GetField(str("area_total")))
                list_obs[obs_id]["area_u"] = float(obs_lake.GetField(str("area_tot_u")))
                list_obs[obs_id]["wse"] = float(obs_lake.GetField(str("wse")))
                list_obs[obs_id]["wse_u"] = float(obs_lake.GetField(str("wse_u")))
            
            # 1.4 - Compute alpha coefficient from prior overlap attribute
            if nb_obs_lake == 1:
                list_obs[obs_id]["alpha"] = 1.0
            else:
                list_obs[obs_id]["alpha"] = prior_overlap[obs_id]
            
        # Reset reading
        self.content_obs.layer.ResetReading()
            
        # 2 - Compute storage change for PLD lake
        out_storage_values["delta_s_l"], out_storage_values["ds_l_u"], out_storage_values["delta_s_q"], out_storage_values["ds_q_u"] = in_obj_plake.run_stocc(list_obs) 
                
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
                obs_lake.SetField(str("delta_s_l"), list_obs[obs_id]["delta_s_l"] + out_storage_values["delta_s_l"]*list_obs[obs_id]["alpha"])
                obs_lake.SetField(str("ds_l_u"), list_obs[obs_id]["ds_l_u"] + out_storage_values["ds_l_u"]*coeff*list_obs[obs_id]["alpha"])
                obs_lake.SetField(str("delta_s_q"), list_obs[obs_id]["delta_s_q"] + out_storage_values["delta_s_q"]*list_obs[obs_id]["alpha"])
                obs_lake.SetField(str("ds_q_u"), list_obs[obs_id]["ds_q_u"] + out_storage_values["ds_q_u"]*coeff*list_obs[obs_id]["alpha"])
            else:
                # Case the observed lake is related to only 1 prior lake
                obs_lake.SetField(str("delta_s_l"), out_storage_values["delta_s_l"]*list_obs[obs_id]["alpha"])
                obs_lake.SetField(str("ds_l_u"), out_storage_values["ds_l_u"]*list_obs[obs_id]["alpha"])
                obs_lake.SetField(str("delta_s_q"), out_storage_values["delta_s_q"]*list_obs[obs_id]["alpha"])
                obs_lake.SetField(str("ds_q_u"), out_storage_values["ds_q_u"]*list_obs[obs_id]["alpha"])
            
            # 3.3 - Rewrite feature with storage change values
            self.content_obs.layer.SetFeature(obs_lake)
            
        # Reset reading
        self.content_obs.layer.ResetReading()
            
        # 4 - Add values to statistics dictionary
        if out_storage_values["delta_s_l"] < my_var.FV_REAL:
            self.compare_stats["obs"]["delta_s_l"] += out_storage_values["delta_s_l"]
        if out_storage_values["delta_s_q"] < my_var.FV_REAL:
            self.compare_stats["obs"]["delta_s_q"] += out_storage_values["delta_s_q"]
        if out_storage_values["delta_s_l"] < my_var.FV_REAL:
            self.compare_stats["prior"]["delta_s_l"] += out_storage_values["delta_s_l"]
        if out_storage_values["delta_s_q"] < my_var.FV_REAL:
            self.compare_stats["prior"]["delta_s_q"] += out_storage_values["delta_s_q"]
        
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
    
    def update_pixcvec_with_imp_geoloc(self, in_pix_index, in_lon, in_lat, in_height):
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
        self.obj_pixc_vec.longitude_vectorproc[pixcvec_index] = in_lon
        self.obj_pixc_vec.latitude_vectorproc[pixcvec_index] = in_lat
        self.obj_pixc_vec.height_vectorproc[pixcvec_index] = in_height
    
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
        self.obj_pixc_vec.obs_id[pixcvec_index] = in_obs_id
        self.obj_pixc_vec.lake_id[pixcvec_index] = in_pixcvec_lakeid

    # ----------------------------------------

    def sort_pixels_wrt_classif_flags(self, in_ind):
        """
        Sort the subset of PixC with indices in_ind wrt classification flags

        :param in_ind: indices of subset of PixC
        :type in_ind: 1D-array of int
        
        :return: dictionary of indices of pixels corresponding to categories "water" and "dark"
        :rtype: dict
        """

        # 0 - Init output dictionary
        out_dict = {}
        out_dict["water"] = None
        out_dict["dark"] = None

        # 1 - Get subset of PixC corresponding to input indices
        tmp_classif = self.obj_pixc.classif[in_ind]

        # 2 - Deal with water flags
        flag_water = self.cfg.get("CONFIG_PARAMS", "FLAG_WATER")
        flag_dark = self.cfg.get("CONFIG_PARAMS", "FLAG_DARK")
        list_classif_flags = flag_water.replace('"','').split(";")
        for classif_flag in list_classif_flags:
            v_ind = np.where(tmp_classif == int(classif_flag))[0]
            if v_ind.size != 0:
                if out_dict["water"] is None:
                    out_dict["water"] = v_ind
                else:
                    out_dict["water"] = np.concatenate((out_dict["water"], v_ind))

        # 3 - Deal with dark water flags
        list_classif_flags = flag_dark.replace('"','').split(";")
        for classif_flag in list_classif_flags:
            v_ind = np.where(tmp_classif == int(classif_flag))[0]
            if v_ind.size != 0:
                if out_dict["dark"] is None:
                    out_dict["dark"] = v_ind
                else:
                    out_dict["dark"] = np.concatenate((out_dict["dark"], v_ind))

        return out_dict

    # ----------------------------------------
    # Function dedicated to high-level stats
    # ----------------------------------------
    
    def compute_obs_stats(self):
        """
        Compute sum of values listed in self.compare_stats_params for the _Obs shapefile
        """
    
        # 1 - Select water features
        self.content_obs.layer.SetAttributeFilter("lake_id != 'no_data'")
        
        # 2 - Compute and store stats
        for lake in self.content_obs.layer:
            for param in self.compare_stats_params:
                tmp_value = float(lake.GetField(param))
                if tmp_value > my_var.FV_REAL:
                    self.compare_stats["obs"][param] += tmp_value
                
        # 3 - Reinitialize reader pointer
        self.content_obs.layer.ResetReading()
        
        # 4 - Remove filter over memory layer
        self.content_obs.layer.SetAttributeFilter(None)
    

#######################################


class LakeTileProduct(LakeProduct):
    """
    class LakeTileProduct
    Manage LakeTile products, as a child class of main LakeProduct class
    """
    def __init__(self, in_obj_pixc, in_obj_pixc_vec, in_obj_lake_db, in_layer_name):
        """
        Constructor
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Init LakeProduct class
        super().__init__("TILE", in_obj_pixc, in_obj_pixc_vec, in_obj_lake_db, in_layer_name)

        # 2 - Find continent associated to tile
        continent = self.obj_lake_db.link_poly_to_continent(self.obj_pixc.tile_poly)
        self.obj_pixc.pixc_metadata["continent"] = continent
        self.obj_pixc_vec.pixcvec_metadata["continent"] = continent

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
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Select water features
        self.content_obs.layer.SetAttributeFilter("lake_id != 'no_data'")
        
        # 2 - Write them in the output shapefile
        my_shp.write_mem_layer_as_shp(self.content_obs.layer, in_filename)
        
        # 3 - Write XML metadatafile for shapefile
        logger.debug("Writing associated metadata file = %s.xml" % in_filename)
        self.content_obs.update_and_write_metadata("%s.xml" % in_filename, 
                                                     in_inprod_metadata=self.obj_pixc.pixc_metadata,
                                                     in_proc_metadata=in_proc_metadata)
        
        # 4 - Remove filter over memory layer
        self.content_obs.layer.SetAttributeFilter(None)
    
    def write_prior_file(self, in_filename, in_proc_metadata):
        """
        Write the PLD-oriented file, i.e. PLD lakes overflown by SWOT
        
        :param in_filename: full path of the output file
        :type in_filename: string
        :param in_proc_metadata: processing metadata
        :type in_proc_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Write them in the output shapefile
        my_shp.write_mem_layer_as_shp(self.content_prior.layer, in_filename)
        
        # 2 - Write XML metadatafile for shapefile
        logger.debug("Writing associated metadata file = %s.xml" % in_filename)
        self.content_prior.update_and_write_metadata("%s.xml" % in_filename, 
                                                     in_inprod_metadata=self.obj_pixc.pixc_metadata,
                                                     in_proc_metadata=in_proc_metadata)
    
    def write_unknown_file(self, in_filename, in_proc_metadata):
        """
        Write the file containing water features unassigned to any prior features (ie neither PRD nor PLD)
        
        :param in_filename: full path of the output file
        :type in_filename: string
        :param in_proc_metadata: processing metadata
        :type in_proc_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Init layer
        tmp_content_unknown = shp_file.LakeTileUnassignedProduct(os.path.basename(in_filename), in_filename=in_filename)
        
        # 2 - Select water features
        self.content_obs.layer.SetAttributeFilter("lake_id = 'no_data'")
        
        # 3 - Write them in the output shapefile
        for cur_feature in self.content_obs.layer:
            cur_att = {}
            for att_name in tmp_content_unknown.attribute_metadata.keys():
                cur_att[att_name] = cur_feature.GetField(str(att_name))
            tmp_content_unknown.add_feature(cur_feature.GetGeometryRef(), cur_att)
        
        # 4 - Write XML metadatafile for shapefile
        logger.debug("Writing associated metadata file = %s.xml" % in_filename)
        tmp_content_unknown.update_and_write_metadata("%s.xml" % in_filename, 
                                                      in_inprod_metadata=self.obj_pixc.pixc_metadata,
                                                      in_proc_metadata=in_proc_metadata)
        
        # 5 - Close shapefile
        tmp_content_unknown.free()
        
        # 6 - Remove filter over memory layer
        self.content_obs.layer.SetAttributeFilter(None)
    

#######################################


class LakeSPProduct(object):
    """
    class LakeSPProduct
    Manage LakeSP products, as a class composed of 2 LakeProduct objects, one for each swath
    """
    def __init__(self, in_obj_pixc_sp, in_obj_pixcvec_sp, in_obj_lake_db, in_layer_name_root, in_continent):
        """
        Constructor

        :param in_obj_pixc_sp: pixel cloud from which to compute lake products
        :type in_obj_pixc: proc_pixc_sp.PixelCloudSP
        :param in_obj_pixcvec_sp: pixel cloud complementary file from which to compute lake products
        :type in_obj_pixcvec_sp: proc_pixc_vec_sp.PixelCloudVecSP
        :param in_obj_lake_db: lake database
        :type in_obj_lake_db: lake_db.lakeDb_shp or lake_db.lakeDb_sqlite
        :param in_layer_name_root: root for the name of right and left lake product layers
        :type in_layer_name_root: string
        :param in_continent: continent code
        :type in_continent: string

        Variables of the object:
            - swath_r / LakeProduct: lake product of right swath
            - swath_l / LakeProduct: lake product of left swath
            - continent / string: continent covered by the LakeSP product
            
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
    
        # 1 - Init right swath
        logger.info("Init lake product object for swath R")
        self.swath_r = LakeProduct("SP",
                                   in_obj_pixc_sp.pixc_edge_r,
                                   in_obj_pixcvec_sp.pixcvec_r,
                                   in_obj_lake_db,
                                   in_layer_name_root + "_R")
    
        # 2 - Init left swath
        logger.info("Init lake product object for swath L")
        self.swath_l = LakeProduct("SP",
                                   in_obj_pixc_sp.pixc_edge_l,
                                   in_obj_pixcvec_sp.pixcvec_l,
                                   in_obj_lake_db,
                                   in_layer_name_root + "_L")
        
        # 3 - Others
        self.continent = in_continent  # Continent
        
    def free_memory(self):
        """
        Destroy memory layer
        """
        self.swath_r.content_obs.free()
        self.swath_r.content_prior.free()
        self.swath_l.content_obs.free()
        self.swath_l.content_prior.free()

    # ----------------------------------------
    
    def write_obs_file(self, in_filename, in_pixc_metadata, in_proc_metadata, in_list_laketile_obs_files):
        """
        Write the observation-oriented file, i.e.
        observed water features related to at least one PLD lake
        
        :param in_filename: full path of the output file
        :type in_filename: string
        :param in_pixc_metadata: metadata retrieved from LakeTile_edge
        :type in_pixc_metadata: dict
        :param in_proc_metadata: processing metadata
        :type in_proc_metadata: dict
        :param in_list_laketile_obs_files: list of LakeTile_Obs shapefiles
        :type in_list_laketile_obs_files: list
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Select water features for each swath
        self.swath_r.content_obs.layer.SetAttributeFilter("lake_id != 'no_data'")
        self.swath_l.content_obs.layer.SetAttributeFilter("lake_id != 'no_data'")
        
        # 2 - Merge right and left swath SP layers into 1 single layer
        dataSource_sp1, layer_sp = my_shp.merge_2_layers(self.swath_r.content_obs.layer, self.swath_l.content_obs.layer, self.continent)
        
        # 3 - Merge obtained layer with LakeTile_Obs shapefiles
        dataSource_sp2, layer_sp = my_shp.merge_mem_layer_with_shp(in_list_laketile_obs_files, layer_sp, self.continent)
        
        # 4 - Write the obtained layer in the output shapefile
        my_shp.write_mem_layer_as_shp(layer_sp, in_filename)
        
        # 5 - Write XML metadatafile for shapefile
        logger.debug("Writing associated metadata file = %s.xml" % in_filename)
        self.swath_r.content_obs.update_and_write_metadata("%s.xml" % in_filename, 
                                                           in_inprod_metadata=in_pixc_metadata,
                                                           in_proc_metadata=in_proc_metadata)
        
        # 6 - Remove filter over memory layer
        self.swath_r.content_obs.layer.SetAttributeFilter(None)
        self.swath_l.content_obs.layer.SetAttributeFilter(None)

        # 7 - Close temporary dataSources
        dataSource_sp1.Destroy()
        dataSource_sp2.Destroy()
    
    def write_prior_file(self, in_filename, in_pixc_metadata, in_proc_metadata, in_list_laketile_prior_files):
        """
        Write the PLD-oriented file, i.e.
        observed water features related to at least one PLD lake
        
        :param in_filename: full path of the output file
        :type in_filename: string
        :param in_pixc_metadata: metadata retrieved from LakeTile_edge
        :type in_pixc_metadata: dict
        :param in_proc_metadata: processing metadata
        :type in_proc_metadata: dict
        :param in_list_laketile_obs_files: list of LakeTile_Obs shapefiles
        :type in_list_laketile_obs_files: list
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Merge right and left swath SP layers into 1 single layer
        dataSource_sp1, layer_sp = my_shp.merge_2_layers(self.swath_r.content_prior.layer, self.swath_l.content_prior.layer, self.continent)
        
        # 2 - Merge obtained layer with LakeTile_Prior shapefiles
        dataSource_sp2, layer_sp = my_shp.merge_mem_layer_with_shp(in_list_laketile_prior_files, layer_sp, self.continent)
        
        # 3 - Write the obtained layer in the output shapefile
        my_shp.write_mem_layer_as_shp(layer_sp, in_filename)
        
        # 4 - Write XML metadatafile for shapefile
        logger.debug("Writing associated metadata file = %s.xml" % in_filename)
        self.swath_r.content_prior.update_and_write_metadata("%s.xml" % in_filename,
                                                             in_inprod_metadata=in_pixc_metadata,
                                                             in_proc_metadata=in_proc_metadata)

        # 5 - Close temporary dataSources
        dataSource_sp1.Destroy()
        dataSource_sp2.Destroy()
    
    def write_unknown_file(self, in_filename, in_pixc_metadata, in_proc_metadata, in_list_laketile_unknown_files):
        """
        Write the file containing water features unassigned to any prior features (ie neither PRD nor PLD)
        
        :param in_filename: full path of the output file
        :type in_filename: string
        :param in_pixc_metadata: metadata retrieved from LakeTile_edge
        :type in_pixc_metadata: dict
        :param in_proc_metadata: processing metadata
        :type in_proc_metadata: dict
        :param in_list_laketile_unknown_files: list of LakeTile_Unassigned shapefiles
        :type in_list_laketile_unknown_files: list
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Init layer
        tmp_content_unknown = shp_file.LakeSPUnassignedProduct(os.path.basename(in_filename))
        
        # 2 - Select water features for each swath
        self.swath_r.content_obs.layer.SetAttributeFilter("lake_id = 'no_data'")
        self.swath_l.content_obs.layer.SetAttributeFilter("lake_id = 'no_data'")
        
        # 3 - Write them in the output shapefile
        # 3.1 - Right swath
        for cur_feature in self.swath_r.content_obs.layer:
            cur_att = {}
            for att_name in tmp_content_unknown.attribute_metadata.keys():
                cur_att[att_name] = cur_feature.GetField(str(att_name))
            tmp_content_unknown.add_feature(cur_feature.GetGeometryRef(), cur_att)
        # 3.2 - Left swath
        for cur_feature in self.swath_l.content_obs.layer:
            cur_att = {}
            for att_name in tmp_content_unknown.attribute_metadata.keys():
                cur_att[att_name] = cur_feature.GetField(str(att_name))
            tmp_content_unknown.add_feature(cur_feature.GetGeometryRef(), cur_att)
        
        # 4 - Merge obtained layer with LakeTile_Unassigned shapefiles
        dataSource_sp2, layer_sp = my_shp.merge_mem_layer_with_shp(in_list_laketile_unknown_files, tmp_content_unknown.layer, self.continent)
        
        # 5 - Write the obtained layer in the output shapefile
        my_shp.write_mem_layer_as_shp(layer_sp, in_filename)
        
        # 6 - Write XML metadatafile for shapefile
        logger.debug("Writing associated metadata file = %s.xml" % in_filename)
        tmp_content_unknown.update_and_write_metadata("%s.xml" % in_filename, 
                                                      in_inprod_metadata=in_pixc_metadata,
                                                      in_proc_metadata=in_proc_metadata)
        
        # 7 - Remove filter over memory layer
        self.swath_r.content_obs.layer.SetAttributeFilter(None)
        self.swath_l.content_obs.layer.SetAttributeFilter(None)
        
        # 6 - Close temporary dataSources
        tmp_content_unknown.free()
        dataSource_sp2.Destroy()


#######################################


def select_water_dark_pixels(in_classif_dict, in_flag_water=False, in_flag_dark=False):
    """
    Merge vectors of indices of classification dictionary wrt to kind of flags wanted
    
    :param in_classif_dict: dictionary of indices of pixels corresponding to categories "water" and "dark"
    :type in_classif_dict: dict (output of self.sort_pixels_wrt_classif_flags)
    :param in_flag_water: =True if water flags selected; =False otherwise (default)
    :type in_flag_water: boolean
    :param in_flag_dark: =True if dark water flags selected; =False otherwise (default)
    :type in_flag_dark: boolean

    :return: list of indices of selected pixels
    :rtype: 1D-array of int
    """

    out_ind = []

    if in_flag_water and (in_classif_dict["water"] is not None):
        out_ind = in_classif_dict["water"]

    if in_flag_dark and (in_classif_dict["dark"] is not None):
        if len(out_ind) == 0:
            out_ind = in_classif_dict["dark"]
        else:
            out_ind = np.concatenate((out_ind, in_classif_dict["dark"]))
            
    return out_ind

    
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
    logger = logging.getLogger("proc_lake.py/select_valid_pixels")
    
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
        logger.warning("!!! All the pixels have NaN for improved geolocation => feature not computed")
        out_valid_flag = 0
    elif nb_nan != 0:
        logger.warning("!!! %d pixels have NaN for improved geolocation => removed from computation" % nb_nan)
        out_valid_flag = 1
        tmp_valid_index1 = np.where(np.isfinite(in_lon))[0]
        tmp_valid_index2 = np.where(np.isfinite(in_lat))[0]
        tmp_valid_index3 = np.where(np.isfinite(in_height))[0]
        tmp_valid_index = np.concatenate((tmp_valid_index1, tmp_valid_index2, tmp_valid_index3))
        out_valid_index = np.unique(tmp_valid_index)
        
    return out_valid_flag, out_valid_index
    

def value_or_none(in_value):
    """
    Test if in_value is finite or not (case of NaN); if not, return None
    
    :param in_value: value to test
    :type in_value: float
    
    :return: out_value = in_value or None if NaN
    :rtype: out_value = float
    """
    
    out_value = in_value
    
    if isinstance(out_value, float):
        if not np.isfinite(out_value):
            out_value = None
            
    return out_value
