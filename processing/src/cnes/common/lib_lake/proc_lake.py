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

    def __init__(self, IN_productType, IN_objPixc, IN_objPixcVec, IN_objLakeDb, IN_layer_name, IN_id_prefix=""):
        """
        Constructor
        
        :param IN_productType: type of product among "SP"=LakeSP and "TILE"=LakeTile
        :type IN_productType: string
        :param IN_objPixc: pixel cloud from which to compute lake products
        :type IN_objPixc: proc_pixc.PixelCloud or proc_pixc_sp.PixelCloudSP object
        :param IN_objPixcVec: pixel cloud complementary file from which to compute lake products
        :type IN_objPixcVec: proc_pixc_vec.PixelCloudVec or proc_pixc_vec_sp.PixelCloudVecSP object
        :param IN_objLakeDb: lake database
        :type IN_objLakeDb: lake_db.lakeDb_shp or lake_db.lakeDb_sqlite
        :param IN_layer_name: name for lake product layer        
        :type IN_layer_name: string
        :param IN_id_prefix: prefix for the lake identifier (default="")
        :type IN_id_prefix: string

        Variables of the object:
        - type / string: 
        - objPixc / proc_pixc.PixelCloud or proc_pixc_sp.PixelCloud: pixel cloud from which to compute lake products
        - objPixcVec / proc_pixc_vec.PixelCloudVec or proc_pixc_vec_sp.PixelCloudVec: extra info for pixel cloud
        - objLakeDb / lake_db.lakeDb_shp or lake_db.lakeDb_sqlite: lake database
        - id_prefix / string: prefix for LAKE_ID
        - dataSource / OGR_dataSource: data source of the product shapefile
        - layer / OGRLayer: layer of the product shapefile
        - layerDefn / OGR_layer_definition: layer definition of the product shapefile
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("[LakeProduct] == INIT ==")
        
        # Init variables
        # Product type
        if (IN_productType != "TILE") and (IN_productType != "SP"):
            my_api.exitWithError("[LakeProduct] ERROR = product type is %s ; should be SP or TILE" % IN_productType)
        else:
            self.type = IN_productType
        # Pixel cloud object
        self.objPixc = IN_objPixc  
        # Pixel cloud complementary file object
        self.objPixcVec = IN_objPixcVec 
        # Lake database object
        self.objLakeDb = IN_objLakeDb
        
        # Prefix for lake identifier
        self.id_prefix = IN_id_prefix        
        
        # Other variables
        self.dataSource = None  # Data source of the product shapefile
        self.layer = None  # Layer of the product shapefile
        self.layerDefn = None  # Layer definition of the product shapefile
        
        # Initialize lake product layer
        self.initProduct(IN_layer_name)

    # ----------------------------------------
    
    def initProduct(self, IN_layer_name):
        """
        Initialize lake product memory layer and attributes creation
        
        :param IN_layer_name: name for lake product layer        
        :type IN_layer_name: string
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("[LakeProduct] == initProduct ==")
        
        memDriver = ogr.GetDriverByName(str('MEMORY'))  # Driver for memory layers

        # 1 - Create memory layer
        logger.info('> Creating memory layer')
        self.dataSource = memDriver.CreateDataSource('memData')
        
        # 2 - Create layer
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)  # WGS84
        self.layer = self.dataSource.CreateLayer(str(IN_layer_name), srs, geom_type=ogr.wkbMultiPolygon)
        
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
        self.layerDefn = self.layer.GetLayerDefn()
        
    def addField_real(self, IN_name, IN_width, IN_precision):
        """
        Add a real field to current layer
        
        :param IN_name: name of the field
        :type IN_name: string
        :param IN_width: size of the field (IN_width + IN_precision = 16)
        :type IN_width: int
        :param IN_precision: precision of the field (IN_width + IN_precision = 16)
        :type IN_precision: int
        """
        tmpField = ogr.FieldDefn(str(IN_name), ogr.OFTReal)  # Init field
        tmpField.SetWidth(IN_width)  # Set width
        tmpField.SetPrecision(IN_precision)  # Set precision
        self.layer.CreateField(tmpField)  # Add field to the current layer
    
    def writeMetadataFile(self, IN_filename):
        """
        Write the metadata file associated to the lake shapefile
        
        :param IN_filename: full path for the metadata file
        :type IN_filename: string
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
        ET.SubElement(tile, "cycle_num").text = "%03d" % self.objPixc.cycle_num
        ET.SubElement(tile, "pass_num").text = "%03d" % self.objPixc.pass_num
        if self.type == "TILE":
            ET.SubElement(tile, "tile_ref").text = str(self.objPixc.tile_ref)
        if my_var.CONTINENT_FILE is not None:
            ET.SubElement(tile, "continent").text = str(self.objPixc.continent)
        
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
        tree.write(IN_filename, pretty_print=True, xml_declaration=True, encoding='utf-8')

    # ----------------------------------------
    
    def computeLakeProducts(self, IN_listLabels):
        """
        Computes lake data products for pixels for which label is in IN_listLabels.
        These products are stored in the shapefile defined by self.shpOut.
        NB: This processing is limited to water bodies being of a minimum size defined by MIN_SIZE. 
        NB2: Improved geolocation is computed for all entities
        
        :param IN_listLabels: list of labels to process
        :type IN_listLabels: 1D-array of int
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
        cpt_tooSmall = 0  # Counter of too small objects
        cpt_obj = 1  # Counter of processed objects
        
        for label in IN_listLabels:  # Loop on inside tile objects
            
            # 1 - Get pixels indices for associated to current label
            pix_idx = np.where(self.objPixc.labels == label)[0]
            obj_nb_pix = pix_idx.size
            
            # 2 - Compute categories wrt classification flags
            classif = self.sortPixelsWrtClassifFlags(pix_idx)

            # 3 - Compute object size = detected area (in m2 => convert in ha)
            obj_size = np.sum(self.objPixc.pixel_area[pix_idx[selectWaterDarkLayoverPixels(classif, IN_flag_water=True, IN_flag_dark=True)]])*1e-4
            
            logger.info("")
            logger.info("[LakeProduct] ===== computeProduct / label = %d / nb pixels = %d / size = %.2f ha =====" % (label, obj_nb_pix, obj_size))
            
            # 4 - Compute mean height ONLY over water pixels
            mean_height = my_tools.computeMean_2sigma(self.objPixc.height[pix_idx[classif["water"]]])
            
            # 5 - Compute improved geolocation if wanted
            if my_var.IMP_GEOLOC:
            
                # 5a - Fit lake height model depending on lake size
                if (my_var.BIGLAKE_MIN_SIZE is not None) and (obj_size >= my_var.BIGLAKE_MIN_SIZE):
                    
                    biglakemodel = BigLakeModel(my_var.BIGLAKE_MODEL)
                    height_model = biglakemodel.height_model

                    logger.info("[LakeProduct] Using {} biglake model for improved geolocation (lake size {} ha)".format(height_model, obj_size))
                    
                    if height_model == 'grid':
                        height_model = biglakemodel.fit_biglake_model(self.objPixc,
                                                                      pix_idx,
                                                                      grid_spacing=my_var.BIGLAKE_GRID_SPACING,
                                                                      grid_resolution=my_var.BIGLAKE_GRID_RES,
                                                                      plot=False)
                        
                    elif height_model == 'polynomial':                                 
                        height_model = biglakemodel.fit_biglake_model_polyfit(self.objPixc, pix_idx)
                                                         
                    else:
                        my_api.printDebug("No height model defined, assume Mean Height model")
                        height_model = np.full(self.objPixc.height[pix_idx].shape, mean_height)

                else:
                    my_api.printDebug("[LakeProduct] Using lake average height {} m for improved geolocation (lake size {} ha)".format(mean_height, obj_size))
                    height_model = np.full(self.objPixc.height[pix_idx].shape, mean_height)

                # 5b - Compute imp geolocation 
                imp_lon, imp_lat, imp_height = proc_pixc_vec.computeImpGeoloc(self.type, self.objPixc, pix_idx, height_model)

            else:
                imp_lon = self.objPixc.longitude[pix_idx]
                imp_lat = self.objPixc.latitude[pix_idx]
                imp_height = self.objPixc.height[pix_idx]

            # Save improved values in objPixcVec
            if self.type == "SP":  # Indices of objPixcVec change depending on product type
                TMP_idx = pix_idx
            else:
                TMP_idx = self.objPixc.selected_idx[pix_idx]
            self.objPixcVec.longitude_vectorproc[TMP_idx] = imp_lon
            self.objPixcVec.latitude_vectorproc[TMP_idx] = imp_lat
            self.objPixcVec.height_vectorproc[TMP_idx] = imp_height

            # Compute lake product if object area large enough
            if obj_size >= my_var.MIN_SIZE:
                
                # 6 - Compute lake identifier
                if self.type == "TILE":  # "TILE" case: only add cpt_obj
                    lake_id = "%s%s" % (self.id_prefix, str(cpt_obj).rjust(my_var.NB_DIGITS, str('0')))
                elif self.type == "SP":  # "SP" case: add main tile info
                    lake_id = "%s%s_%s" % (self.id_prefix, self.objPixc.getMajorityPixelsTileRef(label), str(self.objPixc.getLakeTileLabel(label)).rjust(my_var.NB_DIGITS, str('0')))
                
                # 7 - Compute lake object (geometry and attributes)
                feature = self.computeProduct(lake_id, pix_idx, classif, obj_size, mean_height, imp_lon, imp_lat)
                
                # 8 - Add feature to layer
                self.layer.CreateFeature(feature)
                
                # 9 - Increase counter of processed objects
                cpt_obj += 1
                
            else:
                logger.info("> Object %d too small (%d pixels = %.2f ha)" % (label, obj_nb_pix, obj_size))
                cpt_tooSmall += 1  # Increase counter of too small objects
        
        logger.info("> %d objects not processed because too small" % cpt_tooSmall)
    
    def computeProduct(self, IN_lake_id, IN_indices, IN_classif_dict, IN_size, IN_meanHeight, IN_imp_lon, IN_imp_lat):
        """
        Computes lake product from a subset of pixel cloud, i.e. pixels for which self.objPixc.labels=IN_label
        
        :param IN_lake_id: identifier for the lake
        :type IN_lake_id: string
        :param IN_indices: list of indices of pixels with a IN_id label
        :type IN_indices: 1D-array of int
        :param IN_classif_dict: dictionary of indices of pixels of IN_indices corresponding to categories "water" "dark" and "layover"
        :type IN_classif_dict: dict (output of self.sortPixelsWrtClassifFlags)
        :param IN_size: size of input object
        :type IN_size: float
        :param IN_meanHeight: mean height over the object
        :type IN_meanHeight: float
        :param IN_imp_lon: improved longitudes vector for pixels of the object
        :type IN_imp_lon: 1D-array of float
        :param IN_imp_lat: improved latitudes vector for pixels of the object
        :type IN_imp_lat: 1D-array of float
        
        :return: lake product
        :rtype: OGRFeature
        """
        
        # 1 - Create the object
        OUT_feature = ogr.Feature(self.layerDefn)
            
        # 2 - Compute geometry
        # 2.1 - Compute the lake boundaries
        if self.type == 'SP':
            lakeHull = my_hull.computeLakeBoundaries(IN_imp_lon, 
                                                     IN_imp_lat, 
                                                     self.objPixc.range_idx[IN_indices], 
                                                     self.objPixc.getAzimuthOfLake(IN_indices), 
                                                     self.objPixc.nb_pix_range)
        else:
            lakeHull = my_hull.computeLakeBoundaries(IN_imp_lon,
                                                     IN_imp_lat, 
                                                     self.objPixc.range_idx[IN_indices], 
                                                     self.objPixc.azimuth_idx[IN_indices], 
                                                     self.objPixc.nb_pix_range)

        # 2.2 - Add polygon geometry
        OUT_feature.SetGeometry(lakeHull)
        # 2.3 - Get centroid
        poly_centroid = lakeHull.Centroid().GetPoint(0)
        # 2.4 - Get crosstrack distance and observation time of the centroid
        centroid_ct_dist, centroid_time = self.computeCtTime(poly_centroid)
        # Set crosstrack sign
        if self.objPixc.crosstrack[IN_indices[0]] < 0:
            centroid_ct_dist = -centroid_ct_dist
        
        # 3 - Update attributes
        
        # 3.1 - Lake identifier
        OUT_feature.SetField(str("obslake_id"), str(IN_lake_id))
        
        # 3.2 - Link to a priori database, if specified
        list_prior = None
        pixc_vec_tag = None
        if self.objLakeDb is not None:
            list_prior, pixc_vec_tag = self.objLakeDb.linkToDb(lakeHull, IN_imp_lon, IN_imp_lat)

        if self.type == "SP":  # Indices of objPixcVec change depending on product type
            TMP_idx = IN_indices
        else:
            TMP_idx = self.objPixc.selected_idx[IN_indices]

        if list_prior is None:  # PIXCVec_tag = the id of the lake within the tile
            for ind in TMP_idx:
                if self.objPixcVec.tag[ind] == "":
                    self.objPixcVec.tag[ind] = IN_lake_id
                else:
                    self.objPixcVec.tag[ind] += ";" + IN_lake_id
        else:  
            for ind in TMP_idx:
                if self.objPixcVec.tag[ind] == "":
                    self.objPixcVec.tag[ind] = pixc_vec_tag
                else:
                    self.objPixcVec.tag[ind] += ";" + pixc_vec_tag
            OUT_feature.SetField(str("prior_id"), str(list_prior))
            
        # 3.3 - Mean date of observation
        #OUT_feature.SetField(str("time_day"), -9999)  # Time in UTC days
        #OUT_feature.SetField(str("time_sec"), -9999)  # Time in the day in UTC seconds
        OUT_feature.SetField(str("time_str"), my_tools.convertSec2Time(centroid_time))
        
        # 3.4 - Mean height over the lake and uncertainty
        OUT_feature.SetField(str("height"), float(IN_meanHeight))
        #OUT_feature.SetField(str("height_u"), -9999)
        
        # 3.5 - Height standard deviation (only for big lakes)
        if IN_size >= my_var.BIGLAKE_MIN_SIZE:
            OUT_feature.SetField(str("height_std"), np.std(self.objPixcVec.height_vectorproc[IN_indices[IN_classif_dict["water"]]]))
            
        # 3.6 - Area of detected water pixels and uncertainty (pixel_area in m2 => convert in ha)
        TMP_area_water = np.sum(self.objPixc.pixel_area[IN_indices[IN_classif_dict["water"]]])*1e-4
        OUT_feature.SetField(str("area_detct"), TMP_area_water)
        # Uncertainty
        #OUT_feature.SetField(str("area_det_u"), -9999)
        # 3.7 - Total water area and uncertainty
        OUT_feature.SetField(str("area_total"), float(IN_size))
        # Uncertainty
        polygon_area = my_tools.getArea(lakeHull.Clone(), poly_centroid)*1e-4  # Area in m2 => convert in ha
        TMP_area_u = float(IN_size - polygon_area)
        if TMP_area_u <= 1e8:  # If larger than field width ; **** to be modified later ****
            OUT_feature.SetField(str("area_tot_u"), TMP_area_u)
        
        # 3.8 - Area of pixels used to compute height
        OUT_feature.SetField(str("area_of_ht"), TMP_area_water)
        
        # 3.9 - Metric of layover effect = layover area
        if IN_classif_dict["layover"] is not None:
            OUT_feature.SetField(str("layovr_val"), np.sum(self.objPixc.pixel_area[IN_indices[IN_classif_dict["layover"]]])*1e-4)
        else:
            OUT_feature.SetField(str("layovr_val"), 0)
        
        # 3.10 - Average distance from polygon centroid to the satellite ground track
        OUT_feature.SetField(str("xtrk_dist"), float(centroid_ct_dist))
        
        # 3.11 - Storage change and uncertainty
        #OUT_feature.SetField("delta_s", -9999.0)
        #OUT_feature.SetField("ds_u", -9999.0)
        
        # 3.12 - Dark water flag
        if IN_classif_dict["dark"] is not None:
            OUT_feature.SetField(str("f_dark"), 1)
        else:
            OUT_feature.SetField(str("f_dark"), 0)
        # 3.13 - Ice flag
        #OUT_feature.SetField(str("f_ice"), int(np.where(self.objPixc.ice_flag[IN_indices] == 1)[0].size))
        # 3.14 - Layover flag
        if IN_classif_dict["layover"] is not None:
            OUT_feature.SetField(str("f_layover"), 1)
        else:
            OUT_feature.SetField(str("f_layover"), 0)
        # 3.15 - Quality indicator
        #OUT_feature.SetField("f_quality", -9999.0)
        
        # 3.16 - Partial flag: =1 if the lake is partially covered by the swath, 0 otherwise
        range_idx = self.objPixc.range_idx[IN_indices]
        nr_edge_pix = np.where(range_idx == 0)[0]
        fr_edge_pix = np.where(range_idx == self.objPixc.nb_pix_range-1)[0]
        # Specific processing for 1st and last tiles
        nb_az_edge_pix = 0
        if self.type == "SP":
            if (self.objPixc.tile_idx[IN_indices] == 0).any() or (self.objPixc.tile_idx[IN_indices] == np.max(self.objPixc.tile_idx)).any():
                azimuth_idx = self.objPixc.azimuth_idx[IN_indices]
                if self.objPixc.ascending:
                    az_edge_pix = np.where(azimuth_idx == 0)[0]
                else:
                    az_edge_pix = np.where(azimuth_idx == self.objPixc.nb_pix_azimuth - 1)[0]
                nb_az_edge_pix = az_edge_pix.size
        if nr_edge_pix.size+fr_edge_pix.size+nb_az_edge_pix == 0:
            OUT_feature.SetField(str("f_partial"), 0)
        else:
            OUT_feature.SetField(str("f_partial"), 1)
            
        # 3.17 - Quality of cross-over calibrations
        #OUT_feature.SetField("f_xovr_cal", -9999.0)
        
        # 3.18 - Geoid model height
        #OUT_feature.SetField("geoid_hght", -9999.0)
        # 3.19 - Earth tide
        #OUT_feature.SetField("earth_tide", -9999.0)
        # 3.20 - Pole tide
        #OUT_feature.SetField("pole_tide", -9999.0)
        # 3.21 - Earth tide
        #OUT_feature.SetField("load_tide", -9999.0)
        
        # 3.22 - Dry tropo corr
        #OUT_feature.SetField("c_dry_trop", -9999.0)
        # 3.23 - Wet tropo corr
        #OUT_feature.SetField("c_wet_trop", -9999.0)
        # 3.24 - Iono corr
        #OUT_feature.SetField("c_iono", -9999.0)
        
        # 3.25 - KaRIn measured backscatter averaged for lake
        #OUT_feature.SetField("rdr_sigma0", -9999.0)
        # 3.26 - KaRIn measured backscatter uncertainty for lake 
        #OUT_feature.SetField("rdr_sig0_u", -9999.0)
        # 3.27 - KaRin instrument sigma0 calibration 
        #OUT_feature.SetField("sigma0_cal", -9999.0)
        # 3.28 - sigma0 atmospheric correction within the swath from model data 
        #OUT_feature.SetField("c_sig0_atm", -9999.0)
        # 3.29 - KaRIn correction from crossover cal processing evaluated for lake 
        #OUT_feature.SetField("c_xovr_cal", -9999.0)
        
        # 3.30 - Height correction from KaRIn orientation (attitude) determination
        #OUT_feature.SetField("c_kar_att", -9999.0)
        # 3.31 - Overall instrument system height bias
        #OUT_feature.SetField("c_h_bias", -9999.0)
        # 3.32 - KaRIn to s/c CG correction to height
        #OUT_feature.SetField("c_sys_cg", -9999.0)
        # 3.33 - Corrections on height deduced from instrument internal calibrations if applicable 
        #OUT_feature.SetField("c_intr_cal", -9999.0)
        
        return OUT_feature

    # ----------------------------------------
    
    def sortPixelsWrtClassifFlags(self, IN_ind):
        """
        Sort the subset of PixC with indices IN_ind wrt classification flags
        
        :param IN_ind: indices of subset of PixC
        :type IN_ind: 1D-array of int
        
        :return: dictionary of indices of pixels corresponding to categories "water" "dark" and "layover"
        :rtype: dict
        """
        
        # 0 - Init output dictionary
        OUT_dict = {}
        OUT_dict["water"] = None
        OUT_dict["dark"] = None
        OUT_dict["layover"] = None
        
        # 1 - Get subset of PixC corresponding to input indices
        TMP_classif = self.objPixc.classif[IN_ind]
        
        # 2 - Deal with water flags
        list_classif_flags = my_var.FLAG_WATER.replace('"','').split(";")
        for classif_flag in list_classif_flags:
            vInd = np.where(TMP_classif == int(classif_flag))[0]
            if vInd.size != 0:
                if OUT_dict["water"] is None:
                    OUT_dict["water"] = vInd
                else:
                    OUT_dict["water"] = np.concatenate((OUT_dict["water"], vInd))
        
        # 3 - Deal with dark water flags
        list_classif_flags = my_var.FLAG_DARK.replace('"','').split(";")
        for classif_flag in list_classif_flags:
            vInd = np.where(TMP_classif == int(classif_flag))[0]
            if vInd.size != 0:
                if OUT_dict["dark"] is None:
                    OUT_dict["dark"] = vInd
                else:
                    OUT_dict["dark"] = np.concatenate((OUT_dict["dark"], vInd))
        
        # 4 - Deal with layover flags
        list_classif_flags = my_var.FLAG_LAYOVER.replace('"','').split(";")
        for classif_flag in list_classif_flags:
            vInd = np.where(TMP_classif == int(classif_flag))[0]
            if vInd.size != 0:
                if OUT_dict["layover"] is None:
                    OUT_dict["layover"] = vInd
                else:
                    OUT_dict["layover"] = np.concatenate((OUT_dict["layover"], vInd))
                    
        return OUT_dict

    # ----------------------------------------
    
    def computeCtTime(self, IN_point):
        """
        Compute cross-track distance and observation time for the IN_point among PixC around this point.
        
        :param IN_point: point coordinates
        :type IN_point: OGRPoint
        
        :return: OUT_ct_dist = cross-track distance of the nearest PixC pixel (ie self.objPixc.crosstrack_medium)
        :rtype: float
        :return OUT_time = observation time (ie self.objPixc.nadir_time) of the nadir point of the same azimuth index as the nearest PixC pixel
        :rtype: float
        """
        
        # 1 - Get coordinates of the point
        centroid_lon = IN_point[0]
        centroid_lat = IN_point[1]
        
        # 2 - Compute associated azimuth index
        centroid_az = my_tools.computeAz(centroid_lon, centroid_lat, self.objPixc.nadir_longitude, self.objPixc.nadir_latitude) 
        
        # 3 - Get crosstrack distance
        OUT_ct_dist = my_tools.computeDist(centroid_lon, centroid_lat, self.objPixc.nadir_longitude[centroid_az], self.objPixc.nadir_latitude[centroid_az])
        
        # 4 - Get observation time
        OUT_time = self.objPixc.nadir_time[centroid_az]
        
        return OUT_ct_dist, OUT_time


#######################################
            
    
def selectWaterDarkLayoverPixels(IN_classif_dict, IN_flag_water=False, IN_flag_dark=False, IN_flag_layover=False):
    """
    Merge vectors of indices of classification dictionary wrt to kind of flags wanted
    
    :param IN_classif_dict: dictionary of indices of pixels corresponding to categories "water" "dark" and "layover"
    :type IN_classif_dict: dict (output of self.sortPixelsWrtClassifFlags)
    :param IN_flag_water: =True if water flags selected; =False otherwise (default)
    :type IN_flag_water: boolean
    :param IN_flag_dark: =True if dark water flags selected; =False otherwise (default)
    :type IN_flag_dark: boolean
    :param IN_flag_layover: =True if layover flags selected; =False otherwise (default)
    :type IN_flag_layover: boolean
    
    :return: list of indices of selected pixels
    :rtype: 1D-array of int
    """
    
    OUT_ind = []
    
    if IN_flag_water and (IN_classif_dict["water"] is not None):
        OUT_ind = IN_classif_dict["water"]
    
    if IN_flag_dark and (IN_classif_dict["dark"] is not None):
        if len(OUT_ind) == 0:
            OUT_ind = IN_classif_dict["dark"]
        else:
            OUT_ind = np.concatenate((OUT_ind, IN_classif_dict["dark"]))
        
    if IN_flag_layover and (IN_classif_dict["layover"] is not None):
        if len(OUT_ind) == 0:
            OUT_ind = IN_classif_dict["layover"]
        else:
            OUT_ind = np.concatenate((OUT_ind, IN_classif_dict["layover"]))
            
    return OUT_ind
    
