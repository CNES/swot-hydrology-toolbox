# -*- coding: utf8 -*-
"""
.. module:: proc_lake.py
    :synopsis: Deals with LakeTile and LakeSP shapefile products
    Created on 02/28/2017

.. moduleauthor: Claire POTTIER (CNES DSO/SI/TR) and Cécile CAZALS (CS)
    
.. todo:: revoir la date
.. todo:: revoir le lien de màj du PIXC_VEC tag si plrs lacs associés à l'objet ou si aucun lac (on prend l'identifiant de la tuile ?)
.. todo:: revoir le calcul des erreurs
.. todo:: revoir les identifiants (ID_LAKE - PRIOR_ID - PIXC_VEC tag)
.. todo:: améliorer le calcul de la hauteur moyenne
.. todo:: voir gestion des shapefiles avec Emmanuelle
.. todo:: lien avec la BD: utiliser pyspatialite plutôt que OGR ?
.. todo:: lors du lien avec la BD, màj champs DELTA_S, DS_UNC (avec calculs Jeff) et AREA_EST (plus tard)
.. todo:: pk le prior id est en NULL italique + foncé alors que les autres champs non remplis sont en NULL droit et pâle ?

Copyright (c) 2017 CNES. All rights reserved.
"""
from __future__ import absolute_import, division, print_function, unicode_literals 

import datetime
from lxml import etree as ET
import math
import numpy as np
from osgeo import osr, ogr
import logging

from cnes.modules.geoloc.scripts.biglake_model import BigLakeModel

import cnes.common.lib.my_api as my_api
import cnes.common.lib.my_hull as my_hull
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib_lake.locnes_variables as my_var
import cnes.common.lib_lake.proc_pixc_vec as proc_pixc_vec


class LakeProduct(object):
    """
        class LakeProduct
    """
    def __init__(self, in_product_type, in_obj_pixc, in_obj_pixc_vec, in_obj_lake_db, in_layer_name, IN_id_prefix=""):
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
        - type / string: 
        - obj_pixc / proc_pixc.PixelCloud or proc_pixc_sp.PixelCloud: pixel cloud from which to compute lake products
        - obj_pixc_vec / proc_pixc_vec.PixelCloudVec or proc_pixc_vec_sp.PixelCloudVec: extra info for pixel cloud
        - obj_lake_db / lake_db.lakeDb_shp or lake_db.lakeDb_sqlite: lake database
        - id_prefix / string: prefix for LAKE_ID
        - dataSource / OGR_data_source: data source of the product shapefile
        - layer / OGRLayer: layer of the product shapefile
        - layer_defn / OGR_layer_definition: layer definition of the product shapefile
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("[LakeProduct] == INIT ==")
        
        # Init variables
        # Product type
        if (in_product_type != "TILE") and (in_product_type != "SP"):
            my_api.exitWithError("[LakeProduct] ERROR = product type is %s ; should be SP or TILE" % in_product_type)
        else:
            self.type = in_product_type
        # Pixel cloud object
        self.obj_pixc = in_obj_pixc  
        # Pixel cloud complementary file object
        self.obj_pixc_vec = in_obj_pixc_vec 
        # Lake database object
        self.obj_lake_db = in_obj_lake_db
        
        # Prefix for lake identifier
        self.id_prefix = IN_id_prefix        
        
        # Other variables
        self.dataSource = None  # Data source of the product shapefile
        self.layer = None  # Layer of the product shapefile
        self.layer_defn = None  # Layer definition of the product shapefile
        
        # Initialize lake product layer
        self.initProduct(in_layer_name)

    # ----------------------------------------
    
    def initProduct(self, in_layer_name):
        """
        Initialize lake product memory layer and attributes creation
        
        :param in_layer_name: name for lake product layer        
        :type in_layer_name: string
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("[LakeProduct] == initProduct ==")
        
        mem_driver = ogr.GetDriverByName(str('MEMORY'))  # Driver for memory layers

        # 1 - Create memory layer
        logger.info('> Creating memory layer')
        self.dataSource = mem_driver.CreateDataSource('memData')
        
        # 2 - Create layer
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)  # WGS84
        self.layer = self.dataSource.CreateLayer(str(in_layer_name), srs, geom_type=ogr.wkbMultiPolygon)
        
        # 3 - Create attributes
        # 3.1 - Object identifier
        self.layer.CreateField(ogr.FieldDefn(str('obslake_id'), ogr.OFTString))
        # 3.2 - List of lakes in the a priori DB related to the object
        self.layer.CreateField(ogr.FieldDefn(str('prior_id'), ogr.OFTString))
        # 3.3 - Mean date of observation
        self.layer.CreateField(ogr.FieldDefn(str('time_day'), ogr.OFTInteger))  # Time in UTC days
        self.layer.CreateField(ogr.FieldDefn(str('time_sec'), ogr.OFTInteger))  # Time in the day in UTC seconds
        self.layer.CreateField(ogr.FieldDefn(str('time_str'), ogr.OFTString))  # Time in UTC as a string
        # 3.4 - Mean height over the lake and uncertainty
        self.addField_real('height', 14, 2)
        self.addField_real('height_u', 14, 2)
        # 3.5 - Height standard deviation (only for big lakes)
        self.addField_real('height_std', 14, 2)
        # 3.6 - Area of detected water pixels and uncertainty
        self.addField_real('area_detct', 14, 2)
        self.addField_real('area_det_u', 14, 2)
        # 3.7 - Total water area and uncertainty
        self.addField_real('area_total', 14, 2)
        self.addField_real('area_tot_u', 14, 2)
        # 3.8 - Area of pixels used to compute height
        self.addField_real('area_of_ht', 14, 2)
        # 3.9 - Metric of layover effect
        self.addField_real('layovr_val', 14, 2)
        # 3.10 - Average distance from polygon centroid to the satellite ground track
        self.addField_real('xtrk_dist', 13, 3)
        # 3.11 - Storage change and uncertainty
        self.addField_real('delta_s', 13, 3)
        self.addField_real('ds_u', 13, 3)
        # 3.12 - Dark water flag
        self.layer.CreateField(ogr.FieldDefn(str('f_dark'), ogr.OFTInteger))
        # 3.13 - Ice flag
        self.layer.CreateField(ogr.FieldDefn(str('f_ice'), ogr.OFTInteger))
        # 3.14 - Layover flag
        self.layer.CreateField(ogr.FieldDefn(str('f_layover'), ogr.OFTInteger))
        # 3.15 - Quality indicator
        self.layer.CreateField(ogr.FieldDefn(str('f_quality'), ogr.OFTInteger))
        # 3.16 - Partial flag: =1 if the lake is partially covered by the swath, 0 otherwise
        self.layer.CreateField(ogr.FieldDefn(str('f_partial'), ogr.OFTInteger))
        # 3.17 - Quality of cross-over calibrations
        self.layer.CreateField(ogr.FieldDefn(str('f_xovr_cal'), ogr.OFTInteger))
        # 3.18 - Geoid model height
        self.addField_real('geoid_hght', 13, 3)
        # 3.19 - Earth tide
        self.addField_real('earth_tide', 13, 3)
        # 3.20 - Pole tide
        self.addField_real('pole_tide', 13, 3)
        # 3.21 - Earth tide
        self.addField_real('load_tide', 13, 3)
        # 3.22 - Dry tropo corr
        self.addField_real('c_dry_trop', 13, 3)
        # 3.23 - Wet tropo corr
        self.addField_real('c_wet_trop', 13, 3)
        # 3.24 - Iono corr
        self.addField_real('c_iono', 13, 3)
        # 3.25 - KaRIn measured backscatter averaged for lake
        self.addField_real('rdr_sigma0', 13, 3)
        # 3.26 - KaRIn measured backscatter uncertainty for lake 
        self.addField_real('rdr_sig0_u', 13, 3)
        # 3.27 - KaRin instrument sigma0 calibration 
        self.addField_real('sigma0_cal', 13, 3)
        # 3.28 - sigma0 atmospheric correction within the swath from model data 
        self.addField_real('c_sig0_atm', 13, 3)
        # 3.29 - KaRIn correction from crossover cal processing evaluated for lake 
        self.addField_real('c_xovr_cal', 13, 3)
        # 3.30 - Height correction from KaRIn orientation (attitude) determination
        self.addField_real('c_kar_att', 13, 3)
        # 3.31 - Overall instrument system height bias
        self.addField_real('c_h_bias', 13, 3)
        # 3.32 - KaRIn to s/c CG correction to height
        self.addField_real('c_sys_cg', 13, 3)
        # 3.33 - Corrections on height deduced from instrument internal calibrations if applicable
        self.addField_real('c_intr_cal', 13, 3)
        
        # 4 - Get layer definition
        self.layer_defn = self.layer.GetLayerDefn()
        
    def addField_real(self, in_name, in_width, in_precision):
        """
        Add a real field to current layer
        
        :param in_name: name of the field
        :type in_name: string
        :param in_width: size of the field (in_width + in_precision = 16)
        :type in_width: int
        :param in_precision: precision of the field (in_width + in_precision = 16)
        :type in_precision: int
        """
        tmp_field = ogr.FieldDefn(str(in_name), ogr.OFTReal)  # Init field
        tmp_field.SetWidth(in_width)  # Set width
        tmp_field.SetPrecision(in_precision)  # Set precision
        self.layer.CreateField(tmp_field)  # Add field to the current layer
    
    def writeMetadataFile(self, in_filename):
        """
        Write the metadata file associated to the lake shapefile
        
        :param in_filename: full path for the metadata file
        :type in_filename: string
        """
        
        if self.type == "TILE":
            root = ET.Element("LakeTile_shp")
        elif self.type == "SP":
            root = ET.Element("LakeSP_shp")
        else:
            root = ET.Element("root")

        creation = ET.SubElement(root, "creation")
        ET.SubElement(creation, "producer").text = str(my_var.PRODUCER)
        ET.SubElement(creation, "date").text = str(datetime.datetime.now())

        tile = ET.SubElement(root, "tile_info")
        ET.SubElement(tile, "cycle_num").text = "%03d" % self.obj_pixc.cycle_num
        ET.SubElement(tile, "pass_num").text = "%03d" % self.obj_pixc.pass_num
        if self.type == "TILE":
            ET.SubElement(tile, "tile_ref").text = str(self.obj_pixc.tile_ref)
        if my_var.CONTINENT_FILE is not None:
            ET.SubElement(tile, "continent").text = str(self.obj_pixc.continent)
        
        config = ET.SubElement(root, "config_params")
        ET.SubElement(config, "lake_db").text = str(my_var.LAKE_DB)
        ET.SubElement(config, "lake_db_id").text = str(my_var.LAKE_DB_ID)
        if my_var.CONTINENT_FILE is not None:
            ET.SubElement(config, "continent_file").text = str(my_var.CONTINENT_FILE)
        ET.SubElement(config, "flag_water").text = str(my_var.FLAG_WATER)
        ET.SubElement(config, "flag_dark").text = str(my_var.FLAG_DARK)
        ET.SubElement(config, "flag_layover").text = str(my_var.FLAG_LAYOVER)
        ET.SubElement(config, "min_size").text = str(my_var.MIN_SIZE)
        ET.SubElement(config, "imp_geoloc").text = str(my_var.IMP_GEOLOC)
        ET.SubElement(config, "hull_method").text = str(my_var.HULL_METHOD)
        ET.SubElement(config, "std_height_max").text = str(my_var.STD_HEIGHT_MAX)
        ET.SubElement(config, "biglake_model").text = str(my_var.BIGLAKE_MODEL)
        ET.SubElement(config, "biglake_min_size").text = str(my_var.BIGLAKE_MIN_SIZE)
        ET.SubElement(config, "biglake_grid_spacing").text = str(my_var.BIGLAKE_GRID_SPACING)
        ET.SubElement(config, "biglake_grid_res").text = str(my_var.BIGLAKE_GRID_RES)

        tree = ET.ElementTree(root)
        tree.write(in_filename, pretty_print=True, xml_declaration=True, encoding='utf-8')

    # ----------------------------------------
    
    def computeLakeProducts(self, in_list_labels):
        """
        Computes lake data products for pixels for which label is in in_list_labels.
        These products are stored in the shapefile defined by self.shpOut.
        NB: This processing is limited to water bodies being of a minimum size defined by MIN_SIZE. 
        NB2: Improved geolocation is computed for all entities
        
        :param in_list_labels: list of labels to process
        :type in_list_labels: 1D-array of int
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("[LakeProduct] == computeLakeProducts ==")
        logger.info("[LakeProduct] Minimum size for lakes = %0.1f ha" % my_var.MIN_SIZE)
        if my_var.IMP_GEOLOC:
            logger.info("[LakeProduct] Improve geoloc = YES")
        else:
            logger.info("[LakeProduct] Improve geoloc = NO")
        if my_var.HULL_METHOD == 0:
            logger.info("[LakeProduct] Use of CONVEX HULL for lake boundaries")
        elif math.floor(my_var.HULL_METHOD) == 1:
            logger.info("[LakeProduct] Use of CONCAVE HULL for lake boundaries")
        elif my_var.HULL_METHOD == 2:
            logger.info("[LakeProduct] Use of CONCAVE HULL RADAR VECT for lake boundaries")
        else:
            my_api.exitWithError("[LakeProduct] HULL_METHOD values unkown (%d); should be 0(convex) 1(concav) 2(radar vect)" % my_var.HULL_METHOD)
                    
        # 0 - Init variables
        cpt_too_small = 0  # Counter of too small objects
        cpt_obj = 1  # Counter of processed objects
        
        for label in in_list_labels:  # Loop on inside tile objects
            
            # 1 - Get pixels indices for associated to current label
            pix_idx = np.where(self.obj_pixc.labels == label)[0]
            obj_nb_pix = pix_idx.size
            
            # 2 - Compute categories wrt classification flags
            classif = self.sortPixelsWrtClassifFlags(pix_idx)

            # 3 - Compute object size = detected area (in m2 => convert in ha)
            obj_size = np.sum(self.obj_pixc.pixel_area[pix_idx[selectWaterDarkLayoverPixels(classif, in_flag_water=True, in_flag_dark=True)]])*1e-4
            
            logger.info("")
            logger.info("[LakeProduct] ===== computeProduct / label = %d / nb pixels = %d / size = %.2f ha =====" % (label, obj_nb_pix, obj_size))
            
            # 4 - Compute mean height ONLY over water pixels
            mean_height = my_tools.computeMean_2sigma(self.obj_pixc.height[pix_idx[classif["water"]]])
            
            # 5 - Compute improved geolocation if wanted
            if my_var.IMP_GEOLOC:
            
                # 5a - Fit lake height model depending on lake size
                if (my_var.BIGLAKE_MIN_SIZE is not None) and (obj_size >= my_var.BIGLAKE_MIN_SIZE):
                    
                    biglakemodel = BigLakeModel(my_var.BIGLAKE_MODEL)
                    height_model = biglakemodel.height_model

                    logger.info("[LakeProduct] Using {} biglake model for improved geolocation (lake size {} ha)".format(height_model, obj_size))
                    
                    if height_model == 'grid':
                        height_model = biglakemodel.fit_biglake_model(self.obj_pixc,
                                                                      pix_idx,
                                                                      grid_spacing=my_var.BIGLAKE_GRID_SPACING,
                                                                      grid_resolution=my_var.BIGLAKE_GRID_RES,
                                                                      plot=False)
                        
                    elif height_model == 'polynomial':                                 
                        height_model = biglakemodel.fit_biglake_model_polyfit(self.obj_pixc, pix_idx)
                                                         
                    else:
                        my_api.printDebug("No height model defined, assume Mean Height model")
                        height_model = np.full(self.obj_pixc.height[pix_idx].shape, mean_height)

                else:
                    my_api.printDebug("[LakeProduct] Using lake average height {} m for improved geolocation (lake size {} ha)".format(mean_height, obj_size))
                    height_model = np.full(self.obj_pixc.height[pix_idx].shape, mean_height)

                # 5b - Compute imp geolocation 
                imp_lon, imp_lat, imp_height = proc_pixc_vec.computeImpGeoloc(self.type, self.obj_pixc, pix_idx, height_model)

            else:
                imp_lon = self.obj_pixc.longitude[pix_idx]
                imp_lat = self.obj_pixc.latitude[pix_idx]
                imp_height = self.obj_pixc.height[pix_idx]

            # Save improved values in obj_pixc_vec
            if self.type == "SP":  # Indices of obj_pixc_vec change depending on product type
                tmp_idx = pix_idx
            else:
                tmp_idx = self.obj_pixc.selected_idx[pix_idx]
            self.obj_pixc_vec.longitude_vectorproc[tmp_idx] = imp_lon
            self.obj_pixc_vec.latitude_vectorproc[tmp_idx] = imp_lat
            self.obj_pixc_vec.height_vectorproc[tmp_idx] = imp_height

            # Compute lake product if object area large enough
            if obj_size >= my_var.MIN_SIZE:
                
                # 6 - Compute lake identifier
                if self.type == "TILE":  # "TILE" case: only add cpt_obj
                    lake_id = "%s%s" % (self.id_prefix, str(cpt_obj).rjust(my_var.NB_DIGITS, str('0')))
                elif self.type == "SP":  # "SP" case: add main tile info
                    lake_id = "%s%s_%s" % (self.id_prefix, self.obj_pixc.getMajorityPixelsTileRef(label), str(self.obj_pixc.getLakeTileLabel(label)).rjust(my_var.NB_DIGITS, str('0')))
                
                # 7 - Compute lake object (geometry and attributes)
                feature = self.computeProduct(lake_id, pix_idx, classif, obj_size, mean_height, imp_lon, imp_lat)
                
                # 8 - Add feature to layer
                self.layer.CreateFeature(feature)
                
                # 9 - Increase counter of processed objects
                cpt_obj += 1
                
            else:
                logger.info("> Object %d too small (%d pixels = %.2f ha)" % (label, obj_nb_pix, obj_size))
                cpt_too_small += 1  # Increase counter of too small objects
        
        logger.info("> %d objects not processed because too small" % cpt_too_small)
    
    def computeProduct(self, in_lake_id, in_indices, in_classif_dict, in_size, in_mean_height, in_imp_lon, in_imp_lat):
        """
        Computes lake product from a subset of pixel cloud, i.e. pixels for which self.obj_pixc.labels=in_label
        
        :param in_lake_id: identifier for the lake
        :type in_lake_id: string
        :param in_indices: list of indices of pixels with a in_id label
        :type in_indices: 1D-array of int
        :param in_classif_dict: dictionary of indices of pixels of in_indices corresponding to categories "water" "dark" and "layover"
        :type in_classif_dict: dict (output of self.sortPixelsWrtClassifFlags)
        :param in_size: size of input object
        :type in_size: float
        :param in_mean_height: mean height over the object
        :type in_mean_height: float
        :param in_imp_lon: improved longitudes vector for pixels of the object
        :type in_imp_lon: 1D-array of float
        :param in_imp_lat: improved latitudes vector for pixels of the object
        :type in_imp_lat: 1D-array of float
        
        :return: lake product
        :rtype: OGRFeature
        """
        
        # 1 - Create the object
        out_feature = ogr.Feature(self.layer_defn)
            
        # 2 - Compute geometry
        # 2.1 - Compute the lake boundaries
        if self.type == 'SP':
            lake_hull = my_hull.compute_lake_boundaries(in_imp_lon,
                                                       in_imp_lat,
                                                       self.obj_pixc.range_idx[in_indices],
                                                       self.obj_pixc.getAzimuthOfLake(in_indices),
                                                       self.obj_pixc.nb_pix_range)
        else:
            lake_hull = my_hull.compute_lake_boundaries(in_imp_lon,
                                                       in_imp_lat,
                                                       self.obj_pixc.range_idx[in_indices],
                                                       self.obj_pixc.azimuth_idx[in_indices],
                                                       self.obj_pixc.nb_pix_range)

        # 2.2 - Add polygon geometry
        out_feature.SetGeometry(lake_hull)
        # 2.3 - Get centroid
        poly_centroid = lake_hull.Centroid().GetPoint(0)
        # 2.4 - Get crosstrack distance and observation time of the centroid
        centroid_ct_dist, centroid_time = self.computeCtTime(poly_centroid)
        # Set crosstrack sign
        if self.obj_pixc.crosstrack[in_indices[0]] < 0:
            centroid_ct_dist = -centroid_ct_dist
        
        # 3 - Update attributes
        
        # 3.1 - Lake identifier
        out_feature.SetField(str("obslake_id"), str(in_lake_id))
        
        # 3.2 - Link to a priori database, if specified
        list_prior = None
        pixc_vec_tag = None
        if self.obj_lake_db is not None:
            list_prior, pixc_vec_tag = self.obj_lake_db.linkToDb(lake_hull, in_imp_lon, in_imp_lat)

        if self.type == "SP":  # Indices of obj_pixc_vec change depending on product type
            tmp_idx = in_indices
        else:
            tmp_idx = self.obj_pixc.selected_idx[in_indices]

        if list_prior is None:  # PIXCVec_tag = the id of the lake within the tile
            for ind in tmp_idx:
                if self.obj_pixc_vec.tag[ind] == "":
                    self.obj_pixc_vec.tag[ind] = in_lake_id
                else:
                    self.obj_pixc_vec.tag[ind] += ";" + in_lake_id
        else:  
            for ind in tmp_idx:
                if self.obj_pixc_vec.tag[ind] == "":
                    self.obj_pixc_vec.tag[ind] = pixc_vec_tag
                else:
                    self.obj_pixc_vec.tag[ind] += ";" + pixc_vec_tag
            out_feature.SetField(str("prior_id"), str(list_prior))
            
        # 3.3 - Mean date of observation
        #out_feature.SetField(str("time_day"), -9999)  # Time in UTC days
        #out_feature.SetField(str("time_sec"), -9999)  # Time in the day in UTC seconds
        out_feature.SetField(str("time_str"), my_tools.convertSec2Time(centroid_time))
        
        # 3.4 - Mean height over the lake and uncertainty
        out_feature.SetField(str("height"), float(in_mean_height))
        #out_feature.SetField(str("height_u"), -9999)
        
        # 3.5 - Height standard deviation (only for big lakes)
        if in_size >= my_var.BIGLAKE_MIN_SIZE:
            out_feature.SetField(str("height_std"), np.std(self.obj_pixc_vec.height_vectorproc[in_indices[in_classif_dict["water"]]]))
            
        # 3.6 - Area of detected water pixels and uncertainty (pixel_area in m2 => convert in ha)
        tmp_area_water = np.sum(self.obj_pixc.pixel_area[in_indices[in_classif_dict["water"]]])*1e-4
        out_feature.SetField(str("area_detct"), tmp_area_water)
        # Uncertainty
        #out_feature.SetField(str("area_det_u"), -9999)
        # 3.7 - Total water area and uncertainty
        out_feature.SetField(str("area_total"), float(in_size))
        # Uncertainty
        polygon_area = my_tools.getArea(lake_hull.Clone(), poly_centroid)*1e-4  # Area in m2 => convert in ha
        tmp_area_u = float(in_size - polygon_area)
        if tmp_area_u <= 1e8:  # If larger than field width ; **** to be modified later ****
            out_feature.SetField(str("area_tot_u"), tmp_area_u)
        
        # 3.8 - Area of pixels used to compute height
        out_feature.SetField(str("area_of_ht"), tmp_area_water)
        
        # 3.9 - Metric of layover effect = layover area
        if in_classif_dict["layover"] is not None:
            out_feature.SetField(str("layovr_val"), np.sum(self.obj_pixc.pixel_area[in_indices[in_classif_dict["layover"]]])*1e-4)
        else:
            out_feature.SetField(str("layovr_val"), 0)
        
        # 3.10 - Average distance from polygon centroid to the satellite ground track
        out_feature.SetField(str("xtrk_dist"), float(centroid_ct_dist))
        
        # 3.11 - Storage change and uncertainty
        #out_feature.SetField("delta_s", -9999.0)
        #out_feature.SetField("ds_u", -9999.0)
        
        # 3.12 - Dark water flag
        if in_classif_dict["dark"] is not None:
            out_feature.SetField(str("f_dark"), 1)
        else:
            out_feature.SetField(str("f_dark"), 0)
        # 3.13 - Ice flag
        #out_feature.SetField(str("f_ice"), int(np.where(self.obj_pixc.ice_flag[in_indices] == 1)[0].size))
        # 3.14 - Layover flag
        if in_classif_dict["layover"] is not None:
            out_feature.SetField(str("f_layover"), 1)
        else:
            out_feature.SetField(str("f_layover"), 0)
        # 3.15 - Quality indicator
        #out_feature.SetField("f_quality", -9999.0)
        
        # 3.16 - Partial flag: =1 if the lake is partially covered by the swath, 0 otherwise
        range_idx = self.obj_pixc.range_idx[in_indices]
        nr_edge_pix = np.where(range_idx == 0)[0]
        fr_edge_pix = np.where(range_idx == self.obj_pixc.nb_pix_range-1)[0]
        # Specific processing for 1st and last tiles
        nb_az_edge_pix = 0
        if self.type == "SP":
            if (self.obj_pixc.tile_idx[in_indices] == 0).any() or (self.obj_pixc.tile_idx[in_indices] == np.max(self.obj_pixc.tile_idx)).any():
                azimuth_idx = self.obj_pixc.azimuth_idx[in_indices]
                if self.obj_pixc.ascending:
                    az_edge_pix = np.where(azimuth_idx == 0)[0]
                else:
                    az_edge_pix = np.where(azimuth_idx == self.obj_pixc.nb_pix_azimuth - 1)[0]
                nb_az_edge_pix = az_edge_pix.size
        if nr_edge_pix.size+fr_edge_pix.size+nb_az_edge_pix == 0:
            out_feature.SetField(str("f_partial"), 0)
        else:
            out_feature.SetField(str("f_partial"), 1)
            
        # 3.17 - Quality of cross-over calibrations
        #out_feature.SetField("f_xovr_cal", -9999.0)
        
        # 3.18 - Geoid model height
        #out_feature.SetField("geoid_hght", -9999.0)
        # 3.19 - Earth tide
        #out_feature.SetField("earth_tide", -9999.0)
        # 3.20 - Pole tide
        #out_feature.SetField("pole_tide", -9999.0)
        # 3.21 - Earth tide
        #out_feature.SetField("load_tide", -9999.0)
        
        # 3.22 - Dry tropo corr
        #out_feature.SetField("c_dry_trop", -9999.0)
        # 3.23 - Wet tropo corr
        #out_feature.SetField("c_wet_trop", -9999.0)
        # 3.24 - Iono corr
        #out_feature.SetField("c_iono", -9999.0)
        
        # 3.25 - KaRIn measured backscatter averaged for lake
        #out_feature.SetField("rdr_sigma0", -9999.0)
        # 3.26 - KaRIn measured backscatter uncertainty for lake 
        #out_feature.SetField("rdr_sig0_u", -9999.0)
        # 3.27 - KaRin instrument sigma0 calibration 
        #out_feature.SetField("sigma0_cal", -9999.0)
        # 3.28 - sigma0 atmospheric correction within the swath from model data 
        #out_feature.SetField("c_sig0_atm", -9999.0)
        # 3.29 - KaRIn correction from crossover cal processing evaluated for lake 
        #out_feature.SetField("c_xovr_cal", -9999.0)
        
        # 3.30 - Height correction from KaRIn orientation (attitude) determination
        #out_feature.SetField("c_kar_att", -9999.0)
        # 3.31 - Overall instrument system height bias
        #out_feature.SetField("c_h_bias", -9999.0)
        # 3.32 - KaRIn to s/c CG correction to height
        #out_feature.SetField("c_sys_cg", -9999.0)
        # 3.33 - Corrections on height deduced from instrument internal calibrations if applicable 
        #out_feature.SetField("c_intr_cal", -9999.0)
        
        return out_feature

    # ----------------------------------------
    
    def sortPixelsWrtClassifFlags(self, in_ind):
        """
        Sort the subset of PixC with indices in_ind wrt classification flags
        
        :param in_ind: indices of subset of PixC
        :type in_ind: 1D-array of int
        
        :return: dictionary of indices of pixels corresponding to categories "water" "dark" and "layover"
        :rtype: dict
        """
        
        # 0 - Init output dictionary
        out_dict = {}
        out_dict["water"] = None
        out_dict["dark"] = None
        out_dict["layover"] = None
        
        # 1 - Get subset of PixC corresponding to input indices
        tmp_classif = self.obj_pixc.classif[in_ind]
        
        # 2 - Deal with water flags
        list_classif_flags = my_var.FLAG_WATER.replace('"','').split(";")
        for classif_flag in list_classif_flags:
            v_ind = np.where(tmp_classif == int(classif_flag))[0]
            if v_ind.size != 0:
                if out_dict["water"] is None:
                    out_dict["water"] = v_ind
                else:
                    out_dict["water"] = np.concatenate((out_dict["water"], v_ind))
        
        # 3 - Deal with dark water flags
        list_classif_flags = my_var.FLAG_DARK.replace('"','').split(";")
        for classif_flag in list_classif_flags:
            v_ind = np.where(tmp_classif == int(classif_flag))[0]
            if v_ind.size != 0:
                if out_dict["dark"] is None:
                    out_dict["dark"] = v_ind
                else:
                    out_dict["dark"] = np.concatenate((out_dict["dark"], v_ind))
        
        # 4 - Deal with layover flags
        list_classif_flags = my_var.FLAG_LAYOVER.replace('"','').split(";")
        for classif_flag in list_classif_flags:
            v_ind = np.where(tmp_classif == int(classif_flag))[0]
            if v_ind.size != 0:
                if out_dict["layover"] is None:
                    out_dict["layover"] = v_ind
                else:
                    out_dict["layover"] = np.concatenate((out_dict["layover"], v_ind))
                    
        return out_dict

    # ----------------------------------------
    
    def computeCtTime(self, in_point):
        """
        Compute cross-track distance and observation time for the in_point among PixC around this point.
        
        :param in_point: point coordinates
        :type in_point: OGRPoint
        
        :return: out_ct_dist = cross-track distance of the nearest PixC pixel (ie self.obj_pixc.crosstrack_medium)
        :rtype: float
        :return out_time = observation time (ie self.obj_pixc.nadir_time) of the nadir point of the same azimuth index as the nearest PixC pixel
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
        
        return out_ct_dist, out_time


#######################################
            
    
def selectWaterDarkLayoverPixels(in_classif_dict, in_flag_water=False, in_flag_dark=False, in_flag_layover=False):
    """
    Merge vectors of indices of classification dictionary wrt to kind of flags wanted
    
    :param in_classif_dict: dictionary of indices of pixels corresponding to categories "water" "dark" and "layover"
    :type in_classif_dict: dict (output of self.sortPixelsWrtClassifFlags)
    :param in_flag_water: =True if water flags selected; =False otherwise (default)
    :type in_flag_water: boolean
    :param in_flag_dark: =True if dark water flags selected; =False otherwise (default)
    :type in_flag_dark: boolean
    :param in_flag_layover: =True if layover flags selected; =False otherwise (default)
    :type in_flag_layover: boolean
    
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
        
    if in_flag_layover and (in_classif_dict["layover"] is not None):
        if len(out_ind) == 0:
            out_ind = in_classif_dict["layover"]
        else:
            out_ind = np.concatenate((out_ind, in_classif_dict["layover"]))
            
    return out_ind
    
