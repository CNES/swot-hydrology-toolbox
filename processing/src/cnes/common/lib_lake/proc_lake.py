# -*- coding: utf8 -*-
"""
.. module:: proc_lake.py
    :synopsis: Deals with LakeTile and LakeSP shapefile products
     Created on 2017/02/28

.. moduleauthor: Claire POTTIER (CNES DSO/SI/TR) and Cécile CAZALS (CS)

.. todo:: revoir la date
.. todo:: revoir le lien de màj du PIXC_VEC tag si plrs lacs associés à l'objet ou si aucun lac (on prend l'identifiant de la tuile ?)
.. todo:: revoir le calcul des erreurs
.. todo:: améliorer le calcul de la hauteur moyenne
.. todo:: voir gestion des shapefiles avec Emmanuelle

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import logging

from cnes.modules.geoloc.scripts.biglake_model import BigLakeModel

import cnes.common.lib.my_hull as my_hull
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_variables as my_var2
import cnes.common.lib_lake.locnes_products_shapefile as shp_file
import cnes.common.lib_lake.proc_pixc_vec as proc_pixc_vec
import cnes.common.lib_lake.storage_change as storage_change
import cnes.common.service_config_file as service_config_file
import cnes.common.service_error as service_error


class LakeProduct(object):
    """
        class LakeProduct
    """
    def __init__(self, in_product_type, in_obj_pixc, in_obj_pixc_vec, in_obj_lake_db, in_layer_name, in_id_prefix=""):
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
        :param in_id_prefix: prefix for the lake identifier (default="")
        :type in_id_prefix: string

        Variables of the object:
            - obj_pixc / proc_pixc.PixelCloud or proc_pixc_sp.PixelCloud: pixel cloud from which to compute lake products
            - obj_pixc_vec / proc_pixc_vec.PixelCloudVec or proc_pixc_vec_sp.PixelCloudVec: extra info for pixel cloud
            - obj_lake_db / lake_db.lakeDb_shp or lake_db.lakeDb_sqlite: lake database
            - id_prefix / string: prefix for LAKE_ID
            - shp_mem_layer / LakeTileShp_product: shapefile memory layer of the lake product
            - id_prefix / string: prefix for LAKE_ID
            - uniq_prior_id / set: list of uniq prior identifiers linked to observed objects
        """
        # Get instance of service config file
        self.cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # Init variables
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

        # Prefix for lake identifier
        self.id_prefix = in_id_prefix  
        
        # Initialize lake product layer
        if self.type == "TILE":
            self.shp_mem_layer = shp_file.LakeTileShp_product(in_layer_name, )
        
        # Other variables
        self.uniq_prior_id = set()  # List of uniq prior identifiers linked to observed objects

        # Other variables
#        self.data_source = None  # Data source of the product shapefile
#        self.layer = None  # Layer of the product shapefile
#        self.layer_defn = None  # Layer definition of the product shapefile

    # ----------------------------------------

    def computeLakeProducts(self, in_list_labels):
        """
        Computes lake data products for pixels for which label is in in_list_labels.
        These products are stored in the shapefile defined by self.shpOut.

         - NB: This processing is limited to water bodies being of a minimum size defined by MIN_SIZE.
         - NB2: Improved geolocation is computed for all entities

        :param in_list_labels: list of labels to process
        :type in_list_labels: 1D-array of int
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Minimum size for lakes = %0.1f ha" % self.cfg.getfloat('CONFIG_PARAMS', 'MIN_SIZE'))
        if self.cfg.getboolean('CONFIG_PARAMS', 'IMP_GEOLOC'):
            logger.info("Improve geoloc = YES")
        else:
            logger.info("Improve geoloc = NO")
        hull_method = self.cfg.getfloat("CONFIG_PARAMS", "HULL_METHOD")
        # TODO : les égalités de réel ne sont pas autorisé changer il faut changer ces tests

        if hull_method == 0:
            logger.info("Use of CONVEX HULL for lake boundaries")
        elif hull_method == 1.0:
            logger.info("Use of CONCAVE HULL based on Delaunay triangulation with CGAL library for lake boundaries")
        elif hull_method == 1.1:
            logger.info("Use of CONCAVE HULL based on Delaunay triangulation with varying alpha param for lake boundaries")
        elif hull_method == 1.2:
            logger.info("Use of CONCAVE HULL based on basic Delaunay triangulation for lake boundaries")
        elif hull_method == 2:
            logger.info("Use of RADAR VECTORISATION for lake boundaries")
        else:
            message = "HULL_METHOD values unkown (%d); should be 0(convex) 1(concave with CGAL) 2(radar vect)" % hull_method
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

        for label in in_list_labels:  # Loop on inside tile objects

            # 1 - Get pixels indices for associated to current label
            pix_index = np.where(self.obj_pixc.labels == label)[0]
            obj_nb_pix = pix_index.size
            if obj_nb_pix == 0:
                logger.warning("[STRANGE...] label %s corresponds to 0 pixel..." % label)
                continue
            
            # 2 - Compute categories wrt classification flags
            classif = self.sortPixelsWrtClassifFlags(pix_index)

            # 3 - Compute object size = detected area
            obj_size = np.sum(self.obj_pixc.pixel_area[pix_index[selectWaterDarkPixels(classif, in_flag_water=True, in_flag_dark=True)]])
            
            logger.info("")
            logger.info("===== compute_product / label = %d / nb pixels = %d / size = %.2f m2 =====" % (label, obj_nb_pix, obj_size))
                        
            # 4 - Compute mean height ONLY over water pixels except if there is only dark water
            # TODO: to improve later
            if classif["water"] is None:
                mean_height = my_tools.compute_mean_2sigma(self.obj_pixc.height[pix_index[classif["dark"]]], in_nan=my_var2.FV_FLOAT)
            else:
                mean_height = my_tools.compute_mean_2sigma(self.obj_pixc.height[pix_index[classif["water"]]], in_nan=my_var2.FV_FLOAT)
            # == END-TODO ==
            
            # 5 - Compute improved geolocation if wanted
            if self.cfg.getboolean('CONFIG_PARAMS', 'IMP_GEOLOC'):

                # 5a - Fit lake height model depending on lake size

                if (biglake_model != 'no') and (obj_size >= biglake_min_size):

                    biglakemodel = BigLakeModel(biglake_model)
                    height_model = biglakemodel.height_model

                    logger.info("Using {} biglake model for improved geolocation (lake size {} ha)".format(height_model, obj_size))

                    if height_model == 'grid':
                        height_model = biglakemodel.fit_biglake_model(self.obj_pixc,
                                                                      pix_index,
                                                                      grid_spacing=biglake_grid_spacing,
                                                                      grid_resolution=biglake_grid_res,
                                                                      plot=False)

                    elif height_model == 'polynomial':                                 
                        height_model = biglakemodel.fit_biglake_model_polyfit(self.obj_pixc, pix_index, classif)
                                                         
                    else:
                        logger.debug("No height model defined, assume Mean Height model")
                        height_model = np.full(self.obj_pixc.height[pix_index].shape, mean_height)

                else:
                    logger.debug("Using lake average height {} m for improved geolocation (lake size {} m2)".format(mean_height, obj_size))
                    height_model = np.full(self.obj_pixc.height[pix_index].shape, mean_height)

                # 5b - Compute imp geolocation 
                imp_lon, imp_lat, imp_height = proc_pixc_vec.computeImpGeoloc(self.type, self.obj_pixc, pix_index, height_model)

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

                # 6 - Compute lake identifier
                if self.type == "TILE":  # "TILE" case: only add cpt_obj
                    lake_id = "%s%s" % (self.id_prefix, str(cpt_obj).rjust(nb_digits, str('0')))
                elif self.type == "SP":  # "SP" case: add main tile info
                    lake_id = "%s%s_%s" % (self.id_prefix, self.obj_pixc.getMajorityPixelsTileRef(label), str(self.obj_pixc.getLakeTileLabel(label)).rjust(nb_digits, str('0')))

                # 7 - Compute lake object (geometry and attributes)
                # 7.1 - Test potential issue in output of impGeoloc
                nan_index = np.where(np.isnan(imp_lon))[0]
                nb_nan = len(nan_index)
                # 7.2 - Process on non NaN points 
                feature_geom = None
                if nb_nan == 0:
                    feature_geom, attributes = self.compute_product(lake_id, pix_index, classif, obj_size, mean_height, imp_lon, imp_lat)
                else:
                    if nb_nan == obj_nb_pix:
                        logger.warning("!!! All the pixels have NaN for improved geolocation => object not computed")
                    else:
                        logger.warning("!!! %d pixels have NaN for improved geolocation => removed from computation" % nb_nan)
                        not_nan_index = np.where(np.isfinite(imp_lon))[0]
                        not_nan_classif = self.sortPixelsWrtClassifFlags(pix_index[not_nan_index])
                        feature_geom, attributes = self.compute_product(lake_id, pix_index[not_nan_index], not_nan_classif, obj_size, mean_height, imp_lon[not_nan_index], imp_lat[not_nan_index])

                # 8 - Add feature to layer, if it exists
                if feature_geom is not None:
                    self.shp_mem_layer.add_feature(feature_geom, attributes)
                
                # 9 - Increase counter of processed objects
                cpt_obj += 1

            else:
                logger.info("> Object %d too small (%d pixels = %.2f m2)" % (label, obj_nb_pix, obj_size))
                cpt_too_small += 1  # Increase counter of too small objects

        logger.info("> %d objects not processed because too small" % cpt_too_small)

        # 10 - Compute storage change
        nb_linked = len(self.uniq_prior_id)
        if nb_linked == 0:
            logger.info("NO object linked to a priori lake => NO storage change computed")
        else:
            logger.info("%d objects linked to a priori lake" % nb_linked)
            logger.info("=> Compute storage change")
            self.computeStorageChange()

    def compute_product(self, in_lake_id, in_indices, in_classif_dict, in_size, in_mean_height, in_imp_lon, in_imp_lat):
        """
        Computes lake product from a subset of pixel cloud, i.e. pixels for which self.obj_pixc.labels=in_label

        :param in_lake_id: identifier for the lake
        :type in_lake_id: string
        :param in_indices: list of indices of pixels with a in_id label
        :type in_indices: 1D-array of int
        :param in_classif_dict: dictionary of indices of pixels of in_indices corresponding to categories "water" and "dark"
        :type in_classif_dict: dict (output of self.sortPixelsWrtClassifFlags)
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
            
        # 1 - Compute geometry
        # 1.1 - Find if lake crosses Greenwich meridian
        # If so, longitudes converted in -180/180 (even in shapefile)
        min_long = min(in_imp_lon)
        max_long = max(in_imp_lon)
        if (min_long<180.0) and (max_long>180.0):
            logger.info("Lake %s crosses Greenwich meridian" % in_lake_id)
            geom_long = my_tools.convert_to_m180_180(in_imp_lon)
            self.obj_pixc_vec.greenwich_idx.append(in_indices)
        else:
            geom_long = in_imp_lon
        # 1.2 - Compute the lake boundaries
        if self.type == 'SP':
            out_geom = my_hull.compute_lake_boundaries(geom_long,
                                                       in_imp_lat,
                                                       self.obj_pixc.range_index[in_indices],
                                                       self.obj_pixc.getAzimuthOfLake(in_indices),
                                                       self.obj_pixc.nb_pix_range)
        else:
            out_geom = my_hull.compute_lake_boundaries(geom_long,
                                                       in_imp_lat,
                                                       self.obj_pixc.range_index[in_indices],
                                                       self.obj_pixc.azimuth_index[in_indices],
                                                       self.obj_pixc.nb_pix_range)
        # 1.3 - Get centroid
        poly_centroid = out_geom.Centroid().GetPoint(0)
        # 1.4 - Get crosstrack distance and observation time (UTC and TAI) of the centroid
        centroid_ct_dist, centroid_time, centroid_time_tai = self.computeCtTime(poly_centroid)
        # Set crosstrack sign
        if self.obj_pixc.cross_track[in_indices[0]] < 0:
            centroid_ct_dist = -centroid_ct_dist
        
        # 2 - Update attributes
        
        out_attributes = dict()
        
        # 2.1 - Lake identifier
        out_attributes["obslake_id"] = in_lake_id
        
        # 2.2 - Link to a priori database, if specified
        list_prior = None
        pixc_vec_tag = None
        if self.obj_lake_db is not None:
            list_prior, pixc_vec_tag = self.obj_lake_db.link_to_db(out_geom, geom_long, in_imp_lat)

        if self.type == "SP":  # Indices of obj_pixc_vec change depending on product type
            tmp_index = in_indices
        else:
            tmp_index = self.obj_pixc.selected_index[in_indices]

        if list_prior is None:  # PIXCVec_tag = the id of the lake within the tile
            # Update only PIXCVec_tag
            for ind in tmp_index:
                if self.obj_pixc_vec.other_tag[ind] == "":
                    self.obj_pixc_vec.other_tag[ind] = in_lake_id
                else:
                    self.obj_pixc_vec.other_tag[ind] += ";" + in_lake_id
        else:  
            # Update PIXCVec_tag
            for indpixc, ind in enumerate(tmp_index):
                if self.obj_pixc_vec.lake_tag[ind] == "":   
                    self.obj_pixc_vec.lake_tag[ind] = pixc_vec_tag[indpixc]
                else:
                    self.obj_pixc_vec.taglake_[ind] += ";" + pixc_vec_tag[indpixc]
                out_attributes["prior_id"] = list_prior
            
            # Handle prior_id
            if type(list_prior) == str:
                out_attributes["prior_id"] = list_prior  # Update SHP_prior_id
                self.uniq_prior_id.add(str(list_prior))  # Update list of uniq values of prior IDs
            else:
                out_attributes["prior_id"] = ';'.join(list_prior)   # Update SHP_prior_id
                for ind, p_id in enumerate(list_prior):
                    self.uniq_prior_id.add(str(p_id))  # Update list of uniq values of prior IDs

        # 2.3 - Mean date of observation
        out_attributes["time"] = centroid_time  # UTC time
        out_attributes["time_tai"] = centroid_time_tai  # TAI time
        out_attributes["time_str"] = my_tools.convert_utc_to_str(centroid_time)
        
        # 2.4 - Mean height over the lake and uncertainty
        out_attributes["height"] = in_mean_height
        out_attributes["height_u"] = my_var2.FV_REAL
        
        # 2.5 - Height standard deviation (only for big lakes)
        if in_size >= self.cfg.getfloat("CONFIG_PARAMS", "BIGLAKE_MIN_SIZE"):
            out_attributes["height_std"] = my_tools.compute_std(self.obj_pixc_vec.height_vectorproc[in_indices[in_classif_dict["water"]]], in_nan=my_var2.FV_FLOAT)
            
        # 2.6 - Area of detected water pixels and uncertainty
        tmp_area_water = my_tools.compute_sum(self.obj_pixc.pixel_area[in_indices[in_classif_dict["water"]]])
        out_attributes["area_detct"] = tmp_area_water
        # Uncertainty
        out_attributes["area_det_u"] = my_var2.FV_REAL
        
        # 2.7 - Total water area and uncertainty
        out_attributes["area_total"] = in_size
        # Uncertainty
        polygon_area = my_tools.getArea(out_geom.Clone(), poly_centroid)
        tmp_area_u = in_size - polygon_area
        out_attributes["area_tot_u"] = my_var2.FV_REAL
        if tmp_area_u <= 1e12:  # If larger than field width ; **** to be modified later ****
            out_attributes["area_tot_u"] = tmp_area_u
        
        # 2.8 - Area of pixels used to compute height
        out_attributes["area_of_ht"] = tmp_area_water
        
        # 2.9 - Metric of layover effect = layover area
        #out_attributes["layovr_val"] = my_var2.FV_REAL
        
        # 2.10 - Average distance from polygon centroid to the satellite ground track
        out_attributes["xtrk_dist"] = centroid_ct_dist
        
        # 2.11 - Dark water flag
        if in_classif_dict["dark"] is not None:
            out_attributes["dark_f"] = 1
        else:
            out_attributes["dark_f"] = 0
            
        # 2.12 - Ice flag
        #out_attributes["frozen_f"] = my_var2.FV_REAL
        
        # 2.13 - Layover flag
        #out_attributes["layover_f"] = my_var2.FV_REAL
            
        # 2.14 - Quality indicator
        #out_attributes["quality_f", my_var2.FV_REAL)
        
        # 2.15 - Partial flag: =1 if the lake is partially covered by the swath, 0 otherwise
        range_index = self.obj_pixc.range_index[in_indices]
        nr_edge_pix = np.where(range_index == 0)[0]
        fr_edge_pix = np.where(range_index == self.obj_pixc.nb_pix_range-1)[0]
        # Specific processing for 1st and last tiles
        nb_az_edge_pix = 0
        if self.type == "SP":
            if (self.obj_pixc.tile_index[in_indices] == 0).any() or (self.obj_pixc.tile_index[in_indices] == np.max(self.obj_pixc.tile_index)).any():
                azimuth_index = self.obj_pixc.azimuth_index[in_indices]
                if self.obj_pixc.ascending:
                    az_edge_pix = np.where(azimuth_index == 0)[0]
                else:
                    az_edge_pix = np.where(azimuth_index == self.obj_pixc.nb_pix_azimuth - 1)[0]
                nb_az_edge_pix = az_edge_pix.size
        if nr_edge_pix.size+fr_edge_pix.size+nb_az_edge_pix == 0:
            out_attributes["partial_f"] = 0
        else:
            out_attributes["partial_f"] = 1
            
        # 2.16 - Quality of cross-over calibrations
        #out_attributes["xovr_cal_f"] = my_var2.FV_REAL
        
        # 2.17 - Geoid model height
        out_attributes["geoid_hght"] = my_tools.compute_mean_2sigma(self.obj_pixc.geoid, in_nan=my_var2.FV_FLOAT)
        # 2.18 - Earth tide
        out_attributes["earth_tide"] = my_tools.compute_mean_2sigma(self.obj_pixc.solid_earth_tide, in_nan=my_var2.FV_FLOAT)
        # 2.19 - Pole tide
        out_attributes["pole_tide"] = my_tools.compute_mean_2sigma(self.obj_pixc.pole_tide, in_nan=my_var2.FV_FLOAT)
        # 2.20 - Load tide
        out_attributes["load_tide1"] = my_tools.compute_mean_2sigma(self.obj_pixc.load_tide_sol1, in_nan=my_var2.FV_FLOAT)
        out_attributes["load_tide2"] = my_tools.compute_mean_2sigma(self.obj_pixc.load_tide_sol2, in_nan=my_var2.FV_FLOAT)
        
        # 2.21 - Dry tropo corr
        out_attributes["dry_trop_c"] = my_tools.compute_mean_2sigma(self.obj_pixc.model_dry_tropo_cor, in_nan=my_var2.FV_FLOAT)
        # 2.22 - Wet tropo corr
        out_attributes["wet_trop_c"] = my_tools.compute_mean_2sigma(self.obj_pixc.model_wet_tropo_cor, in_nan=my_var2.FV_FLOAT)
        # 2.23 - Iono corr
        out_attributes["iono_c"] = my_tools.compute_mean_2sigma(self.obj_pixc.iono_cor_gim_ka, in_nan=my_var2.FV_FLOAT)
        
        # 2.24 - KaRIn measured backscatter averaged for lake
        out_attributes["sig0"] = my_tools.compute_mean_2sigma(self.obj_pixc.sig0, in_nan=my_var2.FV_FLOAT)
        # 2.25 - KaRIn measured backscatter uncertainty for lake 
        #out_attributes["sig0_u"] = my_var2.FV_REAL
        # 2.26 - KaRin instrument sigma0 calibration 
        #out_attributes["sig0_cal"] = my_var2.FV_REAL
        # 2.27 - sigma0 atmospheric correction within the swath from model data 
        #out_attributes["sig0_atm_c"] = my_var2.FV_REAL
        # 2.28 - KaRIn correction from crossover cal processing evaluated for lake 
        out_attributes["xovr_cal_c"] = my_tools.compute_mean_2sigma(self.obj_pixc.xover_height_cor, in_nan=my_var2.FV_FLOAT)
        
        # 2.29 - Height correction from KaRIn orientation (attitude) determination
        #out_attributes["kar_att_c"] = my_var2.FV_REAL
        # 2.30 - Overall instrument system height bias
        #out_attributes["h_bias_c"] = my_var2.FV_REAL
        # 2.31 - KaRIn to s/c CG correction to height
        #out_attributes["sys_cg_c"] = my_var2.FV_REAL
        # 2.32 - Corrections on height deduced from instrument internal calibrations if applicable 
        #out_attributes["intr_cal_c"] = my_var2.FV_REAL
        
        return out_geom, out_attributes
    
    def computeStorageChange(self):
        """
        Compute storage change for all objects linked to lakes in a priori database

        Envisionned cases:
            - 1 prior lake <=> 1 observed lake
            - 1 prior lake <=> 2 or more observed lakes
            - 1 observed lake <=> 2 or more prior lakes
            - mixed lakes...
        """
        logger = logging.getLogger(self.__class__.__name__)

        for p_id in self.uniq_prior_id:

            logger.debug("Deal with prior lake %s" % p_id)

            # 1 - Get reference height and area
            ref_height, ref_area = self.obj_lake_db.get_ref_values(p_id)

            # 2 - Get observed lakes linked to this a priori lake
            self.shp_mem_layer.layer.SetAttributeFilter("prior_id LIKE '%{}%'".format(p_id))
            
            # 3 - Number of observed objects linked to this a priori lake
            nb_obs_lake = self.shp_mem_layer.layer.GetFeatureCount()
            
            # 4 - Process wrt to case
            if nb_obs_lake == 1:

                # Get lake feature and values
                obs_lake = self.shp_mem_layer.layer.GetNextFeature()
                obs_height = obs_lake.GetField(str("height"))
                obs_area = obs_lake.GetField(str("area_total"))

                if ";" in obs_lake.GetField(str("prior_id")):  # Case 1 observed lake <=> 2 or more prior lakes
                    pass

                else:  # Case 1 prior lake <=> 1 observed lake

                    # Compute linear storage change
                    stoc_val, stoc_u = storage_change.STOCC_linear(obs_height, obs_area, ref_height, ref_area)
                    # Fill associated field
                    if stoc_val is not None:
                        obs_lake.SetField(str("delta_s_L"), stoc_val)
                    if stoc_u is not None:
                        obs_lake.SetField(str("ds_L_u"), stoc_u)

                    # Compute quadratic storage change
                    stoc_val, stoc_u = storage_change.STOCC_quadratic(obs_height, obs_area, ref_height, ref_area)
                    # Fill associated field
                    if stoc_val is not None:
                        obs_lake.SetField(str("delta_s_Q"), stoc_val)
                    if stoc_u is not None:
                        obs_lake.SetField(str("ds_Q_u"), stoc_u)

                # Rewrite feature with storage change values
                self.shp_mem_layer.layer.SetFeature(obs_lake)
                
            else:  # Case 1 prior lake <=> 2 or more observed lakes
                pass
        
            self.shp_mem_layer.layer.SetAttributeFilter(None)

    # ----------------------------------------

    def sortPixelsWrtClassifFlags(self, in_ind):
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

    def computeCtTime(self, in_point):
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
        centroid_az = my_tools.computeAz(centroid_lon, centroid_lat, self.obj_pixc.nadir_longitude, self.obj_pixc.nadir_latitude)

        # 3 - Get crosstrack distance
        out_ct_dist = my_tools.computeDist(centroid_lon, centroid_lat, self.obj_pixc.nadir_longitude[centroid_az], self.obj_pixc.nadir_latitude[centroid_az])

        # 4 - Get observation time
        out_time = self.obj_pixc.nadir_time[centroid_az]
        out_time_tai = self.obj_pixc.nadir_time[centroid_az]
        
        return out_ct_dist, out_time, out_time_tai

#######################################


def selectWaterDarkPixels(in_classif_dict, in_flag_water=False, in_flag_dark=False):
    """
    Merge vectors of indices of classification dictionary wrt to kind of flags wanted
    
    :param in_classif_dict: dictionary of indices of pixels corresponding to categories "water" and "dark"
    :type in_classif_dict: dict (output of self.sortPixelsWrtClassifFlags)
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

