# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2019/05/17:version initiale.
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

import cnes.common.service_config_file as service_config_file
import cnes.common.service_error as service_error

from cnes.modules.geoloc.scripts.biglake_model import BigLakeModel

import cnes.common.lib.my_hull as my_hull
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_variables as my_var
import cnes.common.lib_lake.locnes_products_shapefile as shp_file
import cnes.common.lib_lake.proc_pixc_vec as proc_pixc_vec
import cnes.common.lib_lake.storage_change as storage_change


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
            - obj_pixc / proc_pixc.PixelCloud or proc_pixc_sp.PixelCloud: pixel cloud from which to compute lake products
            - obj_pixc_vec / proc_pixc_vec.PixelCloudVec or proc_pixc_vec_sp.PixelCloudVec: extra info for pixel cloud
            - obj_lake_db / lake_db.lakeDb_shp or lake_db.lakeDb_sqlite: lake database
            - shp_mem_layer / LakeTileShp_product: shapefile memory layer of the lake product
            - uniq_prior_id / set: list of uniq prior identifiers linked to observed objects
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
            self.type = in_product_type
        # Pixel cloud object
        self.obj_pixc = in_obj_pixc
        # Pixel cloud complementary file object
        self.obj_pixc_vec = in_obj_pixc_vec
        # Lake database object
        self.obj_lake_db = in_obj_lake_db
        
        # 2 - Initialize lake product memory layer
        self.shp_mem_layer = shp_file.LakeSPShpProduct(self.type, in_layer_name)
        
        # 3 - Other variables
        self.uniq_prior_id = set()  # List of uniq prior identifiers linked to all observed objects

        # 4 - ONLY FOR LakeTile: find continent associated to tile
        if self.type == "TILE":
            continent = self.obj_lake_db.link_poly_to_continent(self.obj_pixc.tile_poly)
            self.obj_pixc.pixc_metadata["continent"] = continent
            self.obj_pixc_vec.pixcvec_metadata["continent"] = continent

    # ----------------------------------------

    def compute_lake_products(self, in_list_labels):
        """
        Computes lake data products for pixels for which label is in in_list_labels.
        These products are stored in the shapefile defined by self.shpOut.

         - NB: This processing is limited to water bodies being of a minimum size defined by MIN_SIZE.
         - NB2: Improved geolocation is computed for all entities

        :param in_list_labels: list of labels to process
        :type in_list_labels: 1D-array of int
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Minimum size for lakes = %0.1f km2" % self.cfg.getfloat('CONFIG_PARAMS', 'MIN_SIZE'))
        if self.cfg.getboolean('CONFIG_PARAMS', 'IMP_GEOLOC'):
            logger.info("Improve geoloc = YES")
        else:
            logger.info("Improve geoloc = NO")
        hull_method = self.cfg.getfloat("CONFIG_PARAMS", "HULL_METHOD")
        # TODO : les égalités de réel ne sont pas autorisé il faut changer ces tests

        if hull_method == 0:
            logger.info("Use of CONVEX HULL for lake boundaries")
        elif hull_method == 1.0:
            logger.info("Use of CONCAVE HULL based on Delaunay triangulation with CGAL library for lake boundaries")
        elif hull_method == 1.1:
            logger.info("Use of CONCAVE HULL based on Delaunay triangulation with varying alpha param for lake boundaries")
        elif hull_method == 2:
            logger.info("Use of RADAR VECTORISATION for lake boundaries")
        else:
            message = "HULL_METHOD values unkown (%d); should be 0(convex) 1(concave) 2(radar vect)" % hull_method
            raise service_error.ProcessingError(message, logger)

        # 0 - Init variables
        cpt_too_small = 0  # Counter of too small objects
        cpt_obj = 1  # Counter of processed objects

        biglake_min_size = self.cfg.getfloat("CONFIG_PARAMS", "BIGLAKE_MIN_SIZE")
        biglake_model = self.cfg.get("CONFIG_PARAMS", "BIGLAKE_MODEL")
        biglake_grid_spacing = self.cfg.getint("CONFIG_PARAMS", "BIGLAKE_GRID_SPACING")
        biglake_grid_res = self.cfg.getint("CONFIG_PARAMS", "BIGLAKE_GRID_RES")
        min_size = self.cfg.getfloat('CONFIG_PARAMS', 'MIN_SIZE')
        nb_digits = self.cfg.getint('ID', 'NB_DIGITS')
        
        # Compute all attributes except those related to PLD (ie p_ attributes and storage change)
        for i, label in enumerate(in_list_labels):  # Loop on inside tile objects

            # 1 - Get pixels indices for associated to current label
            pix_index = np.where(self.obj_pixc.labels == label)[0]
            obj_nb_pix = pix_index.size
            if obj_nb_pix == 0:
                logger.warning("[STRANGE...] label %s corresponds to 0 pixel..." % label)
                continue
            
            # 2 - Compute categories wrt classification flags
            classif = self.sort_pixels_wrt_classif_flags(pix_index)
            
            # 3 - Compute object size = detected area
            obj_size = np.sum(self.obj_pixc.pixel_area[pix_index[select_water_dark_pixels(classif, in_flag_water=True, in_flag_dark=True)]])  # In m2
            obj_size /= 10**6  # Conversion in km2
            
            logger.info("")
            logger.info("===== compute_product %d over %d / label = %d / nb pixels = %d / size = %.2f km2 =====" \
                        % (i, len(in_list_labels), label, obj_nb_pix, obj_size))
                        
            # 4 - Compute mean height using both water and dark water pixel (weighting using uncertainties)
            # Modification using JPL aggreagtion method to compute mean/std
            mean_height, weight = my_tools.compute_height_with_uncertainties(self.obj_pixc, pix_index[classif["dark_and_water"]], height = 'only_height')

            
            # 5 - Compute improved geolocation if wanted
            if self.cfg.getboolean('CONFIG_PARAMS', 'IMP_GEOLOC'):

                # 5a - Fit lake height model depending on lake size

                if (biglake_model != 'no') and (obj_size >= biglake_min_size):

                    biglakemodel = BigLakeModel(biglake_model)
                    height_model = biglakemodel.height_model

                    logger.info("Using {} biglake model for improved geolocation (lake size {} km2)".format(height_model, obj_size))

                    if height_model == 'grid':
                        height_model = biglakemodel.fit_biglake_model(self.obj_pixc,
                                                                      pix_index,
                                                                      grid_spacing=biglake_grid_spacing,
                                                                      grid_resolution=biglake_grid_res)


                    elif height_model == 'polynomial':                                 
                        height_model = biglakemodel.fit_biglake_model_polyfit(self.obj_pixc, pix_index, classif)
                                                         
                    else:
                        logger.debug("No height model defined, assume Mean Height model")
                        height_model = np.full(self.obj_pixc.height[pix_index].shape, mean_height)

                else:
                    logger.debug("Using lake average height {} m for improved geolocation (lake size {} km2)".format(mean_height, obj_size))
                    height_model = np.full(self.obj_pixc.height[pix_index].shape, mean_height)

                # 5b - Compute imp geolocation 
                imp_lon, imp_lat, imp_height, p_final = proc_pixc_vec.compute_imp_geoloc(self.type, self.obj_pixc, pix_index, height_model)
                
                # 5c - Compute flattened interferogram in order to compute coherence during uncertainties estimation
                self.obj_pixc.interferogram_flattened[pix_index] = my_tools.compute_interferogram_flatten(self.obj_pixc, pix_index, p_final)
          
            else:
                imp_lon = self.obj_pixc.longitude[pix_index]
                imp_lat = self.obj_pixc.latitude[pix_index]
                imp_height = self.obj_pixc.height[pix_index]


            # Save improved values in obj_pixc_vec
            if self.type == "SP":  # Indices of obj_pixc_vec change depending on product type
                tmp_index = pix_index
            else:
                tmp_index = self.obj_pixc.selected_index[pix_index]
            self.obj_pixc_vec.longitude_vectorproc[tmp_index] = imp_lon
            self.obj_pixc_vec.latitude_vectorproc[tmp_index] = imp_lat
            self.obj_pixc_vec.height_vectorproc[tmp_index] = imp_height

            # Compute lake product if object area large enough
            if obj_size >= min_size:

                # 6 - Compute lakeobs number
                if self.type == "TILE":  # "TILE" case: only add cpt_obj
                    lakeobs_num = str(cpt_obj).rjust(nb_digits, str('0'))
                else :  # "SP" case: add main tile info
                    lakeobs_num = str(self.obj_pixc.get_lake_tile_label(label)).rjust(nb_digits, str('0'))

                # 7 - Compute lake object (geometry and attributes)
                # 7.1 - Test potential issue in output of impGeoloc
                nan_index = np.where(np.isnan(imp_lon))[0]
                nb_nan = len(nan_index)
                # 7.2 - Process on non NaN points 
                feature_geom = None
                if nb_nan == 0:
                    feature_geom, attributes = self.compute_product(lakeobs_num, pix_index, classif, obj_size, mean_height, imp_lon, imp_lat)
                else:
                    if nb_nan == obj_nb_pix:
                        logger.warning("!!! All the pixels have NaN for improved geolocation => object not computed")
                    else:
                        logger.warning("!!! %d pixels have NaN for improved geolocation => removed from computation" % nb_nan)
                        not_nan_index = np.where(np.isfinite(imp_lon))[0]
                        not_nan_classif = self.sort_pixels_wrt_classif_flags(pix_index[not_nan_index])
                        feature_geom, attributes = self.compute_product(lakeobs_num, pix_index[not_nan_index], not_nan_classif,
                                                                        obj_size, mean_height, imp_lon[not_nan_index], imp_lat[not_nan_index])
                logger.debug("obs_id: " + str(attributes["obs_id"]))
                try:
                    logger.debug("lake_id: " + str(attributes["lake_id"]))
                except:
                    logger.debug("lake_id: --------")

                # 8 - Add feature to layer, if it exists
                if feature_geom is not None:
                    self.shp_mem_layer.add_feature(feature_geom, attributes)

                # 9 - Increase counter of processed objects
                cpt_obj += 1

            else:
                logger.info("> Object %d too small (%d pixels = %.2f km2)" % (label, obj_nb_pix, obj_size))
                cpt_too_small += 1  # Increase counter of too small objects

        logger.info("> %d objects not processed because too small" % cpt_too_small)

        # 10 - Compute storage change
        nb_prior = len(self.uniq_prior_id)
        if nb_prior == 0:
            logger.info("NO object linked to a priori lake => NO storage change computed")
            
        else:
            logger.info("%d prior lakes linked to observed lakes" % nb_prior)
            
            for p_id in self.uniq_prior_id:
                logger.info("> Deal with prior lake %s" % p_id)
                
                # 10.1 - Get prior infos
                p_name, p_grand, p_max_wse, p_max_area, p_ref_date, p_ref_ds, p_storage = self.obj_lake_db.get_prior_values(p_id)  
                
                # 10.2 - Set p_ attributes from PLD infos to all observed lakes having this PLD lake as main overlap
                self.shp_mem_layer.layer.SetAttributeFilter("lake_id LIKE '{}%'".format(p_id))
                logger.debug(". {} observed lake(s) are stongly connected to this PLD lake".format(self.shp_mem_layer.layer.GetFeatureCount()))
                self.set_prior_info(p_name, p_grand, p_max_wse, p_max_area, p_ref_date, p_ref_ds, p_storage)
                self.shp_mem_layer.layer.SetAttributeFilter(None)  # Reinit layer attribute filter           
                
                # 10.3 - Compute storage change for all observed lakes overlapping this PLD lake
                logger.debug(". Compute storage change")
                self.shp_mem_layer.layer.SetAttributeFilter("lake_id LIKE '%{}%'".format(p_id))
                self.set_storage_change(p_max_wse, p_max_area, p_ref_ds) 
                self.shp_mem_layer.layer.SetAttributeFilter(None)  # Reinit layer attribute filter
                
    def compute_product(self, in_lakeobs_num, in_indices, in_classif_dict, in_size, in_mean_height, in_imp_lon, in_imp_lat):
        """
        Computes lake product from a subset of pixel cloud, i.e. pixels for which self.obj_pixc.labels=in_label

        :param in_lakeobs_num: lake number
        :type in_lakeobs_num: str
        :param in_indices: list of indices of pixels with a in_id label
        :type in_indices: 1D-array of int
        :param in_classif_dict: dictionary of indices of pixels of in_indices corresponding to categories "water" and "dark"
        :type in_classif_dict: dict (output of self.sort_pixels_wrt_classif_flags)
        :param in_size: size of input object
        :type in_size: float
        :param in_mean_height: mean height over the object
        :type in_mean_height: float
        :param in_imp_lon: improved longitudes vector for pixels of the object
        :type in_imp_lon: 1D-array of float
        :param in_imp_lat: improved latitudes vector for pixels of the object
        :type in_imp_lat: 1D-array of float

        :return: out_geom = lake geometry
        :rtype: OGRPolygon
        :return: out_attributes = lake attributes
        :rtype: dict
        """

        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Deal with object number = {}".format(in_lakeobs_num))
            
        # 1 - Compute geometry
        # 1.1 - Compute the lake boundaries
        if self.type == 'SP':
            out_geom = my_hull.compute_lake_boundaries(in_imp_lon,
                                                       in_imp_lat,
                                                       self.obj_pixc.get_range_of_lake(in_indices),
                                                       self.obj_pixc.get_azimuth_of_lake(in_indices),
                                                       self.obj_pixc.nb_pix_range)
        else:
            out_geom = my_hull.compute_lake_boundaries(in_imp_lon,
                                                       in_imp_lat,
                                                       self.obj_pixc.range_index[in_indices],
                                                       self.obj_pixc.azimuth_index[in_indices],
                                                       self.obj_pixc.nb_pix_range)
        # 1.2 - Get centroid
        poly_centroid = out_geom.Centroid().GetPoint(0)
        # 1.3 - Get crosstrack distance and observation time (UTC and TAI) of the centroid
        centroid_ct_dist, centroid_time, centroid_time_tai = self.compute_ct_time(poly_centroid)
        # Set crosstrack sign
        if self.obj_pixc.cross_track[in_indices[0]] < 0:
            centroid_ct_dist = -centroid_ct_dist
            
        # 2 - Link to lake a priori database
        list_prior, list_fraction, pixc_vec_tag = self.obj_lake_db.link_to_db(out_geom, in_imp_lon, in_imp_lat)
        # Save IDs of intersecting PLD lake 
        if list_prior:
            for tmp_lake_id in list_prior:
                self.uniq_prior_id.add(tmp_lake_id)
        
        # 3 - Compute lake attributes
        out_attributes = dict()  # Init dictionary 
        
        # 3.1 - Identifiers

        # Compute 3 digits of basin identifier (CBB)
        if list_prior:
            lake_basin_id = list_prior[0][0:3]
        else:
            lake_basin_id = str(self.obj_lake_db.link_poly_to_basin(out_geom)[0]).rjust(3, str('0'))

        # 3.1.1 - Compute obs_id
        if self.type == "TILE":  # "TILE" case: only add cpt_obj
            tile_ref = str(self.obj_pixc.pixc_metadata["tile_number"]).rjust(3, str('0')) + str(self.obj_pixc.pixc_metadata["swath_side"])
        else:  # "SP" case: add main tile info
            tile_ref = self.obj_pixc.get_majority_pixels_tile_ref(int(in_lakeobs_num))
        out_attributes["obs_id"] = "%s%s%s" % (lake_basin_id, tile_ref, in_lakeobs_num)
        
        # 3.1.2 - Compute lake_id
        if list_prior:
            out_attributes["lake_id"] = ';'.join(list_prior)
        else:
            out_attributes["lake_id"] = "no_data"
            
        # 3.1.3 - Compute overlap
        if list_fraction:
            out_attributes["overlap"] = ';'.join(list_fraction)
        else:
            out_attributes["overlap"] = "no_data"

        # Update PIXCVec lake_id and obs_id for current lake product
        if self.type == "SP":  # Indices of obj_pixc_vec change depending on product type
            tmp_index = in_indices
        else:
            tmp_index = self.obj_pixc.selected_index[in_indices]
        for i, ind in enumerate(tmp_index):
            self.obj_pixc_vec.obs_id[ind] = out_attributes["obs_id"]  
            self.obj_pixc_vec.lake_id[ind] = pixc_vec_tag[i]

        # 3.2 - Mean datetime of observation
        out_attributes["time"] = centroid_time  # UTC time
        out_attributes["time_tai"] = centroid_time_tai  # TAI time
        out_attributes["time_str"] = my_tools.convert_utc_to_str(centroid_time)  # Time in UTC as a string
        
        # 3.3 - Measured hydrology parameters
        if in_classif_dict["dark_and_water"] is None:
            selected_indices = in_indices
        else:
            selected_indices = in_indices[in_classif_dict["dark_and_water"]]  # Indices of water (and dark) pixels
        # Water surface elevation uncertainties 
        mean_height, height_std_out, height_uncert_out, lat_uncert_out, lon_uncert_out = my_tools.compute_height_with_uncertainties(self.obj_pixc, selected_indices, height = 'corrected')
        
        out_attributes["wse"] = mean_height
        out_attributes["wse_r_u"] = height_uncert_out
        out_attributes["wse_std"] = height_std_out
# ~ =======
        # ~ # Mean water surface elevation over the lake and uncertainty
        # ~ selected_wse = selected_height

        # ~ # Use only values and not fill values
        # ~ valid_geoid = np.where(selected_geoid_hght != my_var.FV_NETCDF[str(selected_geoid_hght.dtype)])
        # ~ selected_wse[valid_geoid] -= selected_geoid_hght[valid_geoid]
# ~ >>>>>>> e94af6fe1ab64d6b285283276400c4d05f600361

        # TO BE COMPLETED, USING SYSTEMATIC ERRORS
        out_attributes["wse_u"] = my_var.FV_REAL
        # TO BE COMPLETED, USING SYSTEMATIC ERRORS

        
        # Sigma0 uncertainties and aggregation (not in product for the moment)
        rdr_sig0, rdr_sig0_std, rdr_sig0_u = my_tools.sig0_with_uncert(self.obj_pixc, selected_indices)
        
        # Total water area and uncertainty
        area, area_unc, area_det, area_det_unc = my_tools.area_with_uncert(self.obj_pixc, selected_indices, method='composite')
        out_attributes["area_total"] = area 
        out_attributes["area_detct"] = area_det
        out_attributes["area_tot_u"] = area_unc
        out_attributes["area_det_u"] = area_det_unc

# ~ =======
        # ~ out_attributes["area_total"] = in_size
        # ~ # Uncertainty
        # ~ polygon_area = my_tools.get_area(out_geom.Clone(), centroid=poly_centroid) / 10**6
        # ~ tmp_area_u = in_size - polygon_area
        # ~ out_attributes["area_tot_u"] = my_var.FV_REAL
        # ~ if tmp_area_u <= 1e12:  # If larger than field width ; **** to be modified later ****
            # ~ out_attributes["area_tot_u"] = tmp_area_u
            
        # ~ # Area of detected water pixels and uncertainty
        # ~ tmp_area_water = my_tools.compute_sum(self.obj_pixc.pixel_area[selected_indices]) / 10**6
        # ~ out_attributes["area_detct"] = tmp_area_water
        # ~ # Uncertainty
        # ~ out_attributes["area_det_u"] = my_var.FV_REAL
# ~ >>>>>>> e94af6fe1ab64d6b285283276400c4d05f600361
        
        # Metric of layover effect = layover area
        out_attributes["layovr_val"] = my_var.FV_REAL
        
        # Average distance from polygon centroid to the satellite ground track
        out_attributes["xtrk_dist"] = centroid_ct_dist
        
        # 3.5 - Quality indicators
        
        # Summary quality indicator
        out_attributes["quality_f"] = 0
        
        # Fractional area of dark water
        if in_classif_dict["dark"] is None:
            out_attributes["dark_frac"] = 0
        else:
            out_attributes["dark_frac"] = my_tools.compute_sum(self.obj_pixc.pixel_area[in_indices[in_classif_dict["dark"]]]) / 10**6 / in_size
            
        # Ice cover flags
        out_attributes["ice_clim_f"] = 255
        out_attributes["ice_dyn_f"] = 255
        
        # Partial flag: =1 if the lake is partially covered by the swath, 0 otherwise
        range_index = self.obj_pixc.range_index[in_indices]
        nr_edge_pix = np.where(range_index == 0)[0]
        fr_edge_pix = np.where(range_index == self.obj_pixc.nb_pix_range-1)[0]
        # Specific processing for 1st and last tiles
        az_edge_pix = np.array([])
        if self.type == "SP":
            # get pixels that belong to the fist or last azimuth line
            az_edge_pix = np.where(self.obj_pixc.is_boundary_pix[in_indices])[0]
        if nr_edge_pix.size+fr_edge_pix.size+az_edge_pix.size == 0:
            out_attributes["partial_f"] = 0
        else:
            out_attributes["partial_f"] = 1
            
        # 3.6 - Geophysical references
        # Geoid model height
        out_attributes["geoid_hght"] = my_tools.compute_mean_2sigma(self.obj_pixc.geoid[selected_indices], in_nan=my_var.FV_FLOAT)
        # Earth tide
        out_attributes["solid_tide"] = my_tools.compute_mean_2sigma(self.obj_pixc.solid_earth_tide[selected_indices], in_nan=my_var.FV_FLOAT)
        # Pole tide
        out_attributes["pole_tide"] = my_tools.compute_mean_2sigma(self.obj_pixc.pole_tide[selected_indices], in_nan=my_var.FV_FLOAT)
        # Load tide
        out_attributes["load_tide1"] = my_tools.compute_mean_2sigma(self.obj_pixc.load_tide_sol1[selected_indices], in_nan=my_var.FV_FLOAT)
        out_attributes["load_tide2"] = my_tools.compute_mean_2sigma(self.obj_pixc.load_tide_sol2[selected_indices], in_nan=my_var.FV_FLOAT)
        
        # 2.7 - Geophysical range corrections
        # Dry tropo corr
        out_attributes["dry_trop_c"] = my_tools.compute_mean_2sigma(self.obj_pixc.model_dry_tropo_cor[in_indices], in_nan=my_var.FV_FLOAT)
        # Wet tropo corr
        out_attributes["wet_trop_c"] = my_tools.compute_mean_2sigma(self.obj_pixc.model_wet_tropo_cor[in_indices], in_nan=my_var.FV_FLOAT)
        # Iono corr
        out_attributes["iono_c"] = my_tools.compute_mean_2sigma(self.obj_pixc.iono_cor_gim_ka[in_indices], in_nan=my_var.FV_FLOAT)
        
        # 2.8 - Instrument corrections
        # KaRIn correction from crossover cal processing evaluated for lake
        out_attributes["xovr_cal_c"] = my_tools.compute_mean_2sigma(self.obj_pixc.xover_height_cor[in_indices], in_nan=my_var.FV_FLOAT)
        
        return out_geom, out_attributes

    # ----------------------------------------
    
    def set_prior_info(self, in_name, in_grand, in_max_wse, in_max_area, in_ref_date, in_ref_ds, in_storage):
        """
        Set p_ attributes to prior values for all objects linked to current PLD lake
        
        :param in_name: name of the prior lake
        :type in_name: string
        :param in_grand: GRanD identifier 
        :type in_grand: int
        :param in_max_wse: PLD lake maximum water surface elevation
        :type in_max_wse: float
        :param in_max_area: PLD lake maximum area
        :type in_max_area: float
        :param in_ref_date: reference date for storage change computation
        :type in_ref_date: string
        :param in_ref_ds: reference storage change for storage change computation
        :type in_ref_ds: float
        :param in_storage: PLD lake maximum storage value
        :type in_storage: float
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        for obs_lake in self.shp_mem_layer.layer:
            
            # 1 - Set PLD name
            if in_name is not None:
                obs_lake.SetField(str("p_name"), in_name)
            else:
                obs_lake.SetField(str("p_name"), my_var.FV_STRING_SHP)
                    
            # 2 - Set GRanD identifier
            if in_grand is not None:
                obs_lake.SetField(str("p_grand_id"), in_grand)
            else:
                obs_lake.SetField(str("p_grand_id"), my_var.FV_INT9_SHP)

            # 3 - Set max water surface elevation
            if in_max_wse is not None:
                obs_lake.SetField(str("p_max_wse"), in_max_wse)
            else:
                obs_lake.SetField(str("p_max_wse"), my_var.FV_REAL)

            # 4 - Set max area
            if in_max_area is not None:
                obs_lake.SetField(str("p_max_area"), in_max_area)
            else:
                obs_lake.SetField(str("p_max_area"), my_var.FV_REAL)

            # 5 - Set reference date
            if in_ref_date is not None:
                obs_lake.SetField(str("p_ref_date"), in_ref_date)
            else:
                obs_lake.SetField(str("p_ref_date"), my_var.FV_STRING_SHP)

            # 6 - Set refence storage change
            if in_ref_ds is not None:
                obs_lake.SetField(str("p_ref_ds"), in_ref_ds)
            else:
                obs_lake.SetField(str("p_ref_ds"), my_var.FV_REAL)
                    
            # 5 - Update maximum water storage value
            if in_storage is not None:
                obs_lake.SetField(str("p_storage"), in_storage)
            else:
                obs_lake.SetField(str("p_storage"), my_var.FV_REAL)
        
            # 6 - Rewrite feature with updated attributes
            self.shp_mem_layer.layer.SetFeature(obs_lake)

    def set_storage_change(self, in_max_wse, in_max_area, in_ref_ds):
        """
        Set storage change value for all objects linked to current prior lake
        
        :param in_max_height: maximum water surface elevation
        :type in_max_height: float
        :param in_max_area: maximum area
        :type in_max_area: float
        :param in_ref_ds: reference storage change
        :type in_ref_ds: float

        Envisionned cases:
            - 1 prior lake <=> 1 observed lake
            - 1 prior lake <=> 2 or more observed lakes
            - 1 observed lake <=> 2 or more prior lakes
            - mixed lakes...
        """

        # 1 - Number of observed objects linked to this a priori lake
        nb_obs_lake = self.shp_mem_layer.layer.GetFeatureCount()
            
        # 2 - Process wrt to case
        if nb_obs_lake == 1:

            # Get lake feature and values
            obs_lake = self.shp_mem_layer.layer.GetNextFeature()
            obs_height = obs_lake.GetField(str("wse"))
            obs_area = obs_lake.GetField(str("area_total"))

            if ";" in obs_lake.GetField(str("lake_id")):  # Case 1 observed lake <=> 2 or more prior lakes
                pass

            else:  # Case 1 prior lake <=> 1 observed lake

                # Compute linear storage change
                stoc_val, stoc_u = storage_change.stocc_linear(obs_height, obs_area, in_max_wse, in_max_area)
                # Fill associated field
                if stoc_val is not None:
                    if in_ref_ds is not None:
                        obs_lake.SetField(str("delta_s_L"), stoc_val - in_ref_ds)
                    else:
                        obs_lake.SetField(str("delta_s_L"), stoc_val)
                if stoc_u is not None:
                    obs_lake.SetField(str("ds_L_u"), stoc_u)

                # Compute quadratic storage change
                stoc_val, stoc_u = storage_change.stocc_quadratic(obs_height, obs_area, in_max_wse, in_max_area)
                # Fill associated field
                if stoc_val is not None:
                    if in_ref_ds is not None:
                        obs_lake.SetField(str("delta_s_Q"), stoc_val - in_ref_ds)
                    else:
                        obs_lake.SetField(str("delta_s_Q"), stoc_val)
                if stoc_u is not None:
                    obs_lake.SetField(str("ds_Q_u"), stoc_u)

            # Rewrite feature with storage change values
            self.shp_mem_layer.layer.SetFeature(obs_lake)
                
        else:  # Case 1 prior lake <=> 2 or more observed lakes
            pass

    # ----------------------------------------

    def sort_pixels_wrt_classif_flags(self, in_ind):
        """
        Sort the subset of PixC with indices in_ind wrt classification flags

        :param in_ind: indices of subset of PixC
        :type in_ind: 1D-array of int
        
        :return: dictionary of indices of pixels corresponding to categories "water" and "dark" and "dark_and_water"
        :rtype: dict
        """

        # 0 - Init output dictionary
        out_dict = {}
        out_dict["water"] = None
        out_dict["dark"] = None
        out_dict["dark_and_water"] = None

        # 1 - Get subset of PixC corresponding to input indices
        tmp_classif = self.obj_pixc.classif[in_ind]

        # 2 - Deal with water flags
        flag_water = self.cfg.get("CONFIG_PARAMS", "FLAG_WATER")
        flag_dark = self.cfg.get("CONFIG_PARAMS", "FLAG_DARK")
        list_classif_flags_water = flag_water.replace('"','').split(";")
        for classif_flag in list_classif_flags_water:
            v_ind = np.where(tmp_classif == int(classif_flag))[0]
            if v_ind.size != 0:
                if out_dict["water"] is None:
                    out_dict["water"] = v_ind
                else:
                    out_dict["water"] = np.concatenate((out_dict["water"], v_ind))
        # 3 - Deal with dark water flags
        list_classif_flags_dark = flag_dark.replace('"','').split(";")
        for classif_flag in list_classif_flags_dark:
            v_ind = np.where(tmp_classif == int(classif_flag))[0]
            if v_ind.size != 0:
                if out_dict["dark"] is None:
                    out_dict["dark"] = v_ind
                else:
                    out_dict["dark"] = np.concatenate((out_dict["dark"], v_ind))

        list_classif_flags = list_classif_flags_dark + list_classif_flags_water
        for classif_flag in list_classif_flags:
            v_ind = np.where(tmp_classif == int(classif_flag))[0]
            if v_ind.size != 0:
                if out_dict["dark_and_water"] is None:
                    out_dict["dark_and_water"] = v_ind
                else:
                    out_dict["dark_and_water"] = np.concatenate((out_dict["dark_and_water"], v_ind))
                    
        return out_dict

    # ----------------------------------------

    def compute_ct_time(self, in_point):
        """
        Compute cross-track distance and observation time for the in_point among PixC around this point.

        :param in_point: point coordinates
        :type in_point: OGRPoint
        :return: out_ct_dist cross-track distance of the nearest PixC pixel (ie self.obj_pixc.crosstrack_medium)
        :rtype: float
        :return: out_time = UTC time (ie self.obj_pixc.nadir_time) of the nadir point of the same azimuth index as the nearest PixC pixel
        :rtype: float
        :return: out_time_tai = TAI time (ie self.obj_pixc.nadir_time_tai) of the nadir point of the same azimuth index as the nearest PixC pixel
        :rtype: float
        """

        # 1 - Get coordinates of the point
        centroid_lon = in_point[0]
        centroid_lat = in_point[1]

        # 2 - Compute associated azimuth index
        centroid_az = my_tools.compute_az(centroid_lon, centroid_lat, self.obj_pixc.nadir_longitude, self.obj_pixc.nadir_latitude)

        # 3 - Get crosstrack distance
        out_ct_dist = my_tools.compute_dist(centroid_lon, centroid_lat, self.obj_pixc.nadir_longitude[centroid_az],
                                            self.obj_pixc.nadir_latitude[centroid_az])

        # 4 - Get observation time
        out_time = self.obj_pixc.nadir_time[centroid_az]
        out_time_tai = self.obj_pixc.nadir_time[centroid_az]
        
        return out_ct_dist, out_time, out_time_tai
    

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

