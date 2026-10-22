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
# VERSION:5.2.1:DM:#91:2023/07/24:Poursuite industrialisation
# VERSION:5.3.0:DM:#91:2023/10/04:Poursuite industrialisation
# VERSION:6.1.0:DM:#91:2024/09/20:Poursuite industrialisation
# VERSION:6.2.0:DM:#91:2024/11/29:Poursuite industrialisation
# VERSION:6.2.6:DM:#900:2026/07/06: [PGE] [LAKES] Ajout d'une interface niveau moyen des lacs
# VERSION:6.2.6:DM:#900:2026/07/06:SSALTO-565 [PGE] [LAKES] Ajout d'une interface coeff NN pour calcul ice_dyn_f
# VERSION:6.2.7:DM:#906:2026/07/28:SSALTO-872 [PGE][LAKES] Gestion des erreurs/warning GDAL
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: locnes_products_shapefile.py
    :synopsis: Deals with SWOT shapefile products
     Created on 02/26/2019

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
"""
from collections import OrderedDict
import logging
import os
import numpy as np
from lxml import etree as ET
from osgeo import ogr, osr

import cnes.common.service_config_file as service_config_file
import cnes.common.service_error as service_error

import cnes.common.lib.my_variables as my_var
import cnes.common.lib_lake.locnes_product_general as my_prod
import cnes.common.lib.my_tools as my_tools


class ShapefileProduct(my_prod.LocnesProduct):
    """
    Deal with SWOT Shapefile products: 
        - LakeTile_Obs, LakeTile_Prior, LakeTile_Unassigned
        - LakeSP_Obs, LakeSP_Prior, LakeSP_Unassigned
        - LakeAvg
    """
    
    def __init__(self, in_xml_file, in_layer_name, in_filename=None, 
                 in_inprod_metadata=None, in_proc_metadata=None,
                 flag_create=True):
        """
        Constructor: set the general values
        
        :param in_xml_file: XML file to populate the variables of the object
        :type in_xml_file: string
        :param in_layer_name: name for shapefile layer        
        :type in_layer_name: string
        :param in_filename: full path of the shapefile for the layer; if None, create a memory layer (default)
        :param in_filename: string     
        :param in_inprod_metadata: metadata specific to input data product
        :type in_inprod_metadata: dict
        :param in_proc_metadata: metadata specific to current processing
        :type in_proc_metadata: dict
        :param flag_create: =True (default) to create a shapefile or a memory layer, =False otherwise
        :type flag_create: boolean
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # 1 - Inheritance of LocnesProduct
        super().__init__(in_xml_file,
                         in_inprod_metadata=in_inprod_metadata,
                         in_proc_metadata=in_proc_metadata)
        
        # 2 - Shapefile attributes
        self.data_source = None  # Data source
        self.layer = None  # Data layer
        self.layer_defn = None  # Layer definition
        self.xml_tree = None  # XML tree
            
        # 3 - Create layer
        if flag_create:
            self.create_layer(in_layer_name, ogr.wkbMultiPolygon, in_filename=in_filename)
        
    def free(self):
        """
        Destroy memory layer
        """
        if self.data_source is not None:
            self.data_source.Destroy()

    #----------------------------------------
    
    def set_from_xml(self, in_xml_file):
        """
        Populate the variables of the object from an XML file
        This XML file has been generated from the Product Description Document with the tool
        tools/pdd2xml/run_pdd2xml.py
        
        :param in_xml_file: XML file to populate the variables of the object
        :type in_xml_file: string
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("Read content from XML file = %s" % in_xml_file)
        
        # 1 - Load the XML file
        xml_reader = ET.parse(in_xml_file)
        root = xml_reader.getroot()
        
        # 2 - Get the file type
        for element_niv1 in root:
            if element_niv1.tag == "science":
                self.file_type = element_niv1.get("uid")
        # Infer primary key
        self.primary_key = "obs_id"
        if self.file_type.endswith("_prior") or self.file_type.endswith("_lakeavg"):
            self.primary_key = "lake_id"
                        
        # 3 - Scan global metadata
        level_global_metadata = root[0][0]
        for element in level_global_metadata:
            cur_metadata = element.get("name")
            
            # Test existence (should not already be in the list)
            if cur_metadata in self.global_metadata:
                logger.warning("Metadata %s is repeated in the XML %s" % (cur_metadata, in_xml_file))
                
            else:  # Store information in a dictionary
                
                # Init
                self.global_metadata[cur_metadata] = OrderedDict()
                # Type attributes
                self.global_metadata[cur_metadata]["type"] = str(element.tag)
                # Description
                annotation = element[0]
                self.global_metadata[cur_metadata]["description"] = annotation.get("description")
                # Init value
                if cur_metadata in ["Conventions", "title", "short_name", "product_file_id", "platform", "references",\
                                    "reference_document", "pge_name"]:
                    self.global_metadata[cur_metadata]["value"] = self.global_metadata[cur_metadata]["description"]
                else:
                    if self.global_metadata[cur_metadata]["type"] == "string":
                        self.global_metadata[cur_metadata]["value"] = ""
                    else:
                        self.global_metadata[cur_metadata]["value"] = -9999
            
        # 4 - Scan attribute metadata
        level_attribute_metadata = root[0][1]
        for element in level_attribute_metadata:
            cur_attribute = element.get("name")
        
            # On initialise le dictionnaire associé s'il n'existe pas
            if cur_attribute in self.attribute_metadata.keys():
                logger.warning("Variable %s is repeated in the XML %s" % (cur_attribute, in_xml_file))
            else:
                self.attribute_metadata[cur_attribute] = OrderedDict()
                
            # Type
            self.attribute_metadata[cur_attribute]["type"] = str(element.tag)
            if "width" in element.keys():
                self.attribute_metadata[cur_attribute]["width"] = int(element.get("width"))
            if "precision" in element.keys():
                self.attribute_metadata[cur_attribute]["precision"] = int(element.get("precision"))
            
            # Other
            annotation = element[0]
            for key, value in annotation.items():
                self.attribute_metadata[cur_attribute][key] = value
                
    #----------------------------------------
    
    def create_layer(self, in_layer_name, in_geom_type, in_spatial_ref=4326, in_filename=None):
        """
        Create product layer and attributes
        
        :param in_layer_name: name for product layer        
        :type in_layer_name: string
        :param in_geom_type: type of geometry
        :type in_geom_type: OGRGeometry
        :param in_spatial_ref: name of the spatial reference (default=4326=WGS84)
        :type in_spatial_ref: int   
        :param in_filename: full path of the shapefile for the layer; if None, create a memory layer (default)
        :param in_filename: string     
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Create data source
        # 1.1 - Memory layer
        if in_filename is None:
            logger.debug('> Creating memory layer')
            driver = ogr.GetDriverByName(str('MEMORY'))  # Driver for memory container
            self.data_source = driver.CreateDataSource('memData')
        else:
            logger.debug('> Creating shapefile layer')
            driver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Driver for shapefile 
            if os.path.exists(in_filename):
                logger.warning("Output shapefile %s already exists => delete file" % in_filename)
                driver.DeleteDataSource(in_filename)
            else:
                logger.debug("Output shapefile = %s" % in_filename)
            self.data_source = driver.CreateDataSource(in_filename)
        
        # 2 - Create layer
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(in_spatial_ref)  # WGS84
        self.layer = self.data_source.CreateLayer(str(in_layer_name), srs, geom_type=in_geom_type)
        
        # 3 - Create attributes
        self.create_attributes()
        
        # 4 - Get layer definition
        self.layer_defn = self.layer.GetLayerDefn()
        
    def create_attributes(self):
        """
        Create layer attributes depending on their type
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        for key, value in self.attribute_metadata.items():
            
            if value['type'] in ["int", "int4", "int9", "float", "text"]:
                tmp_field = ogr.FieldDefn(str(key), my_var.FORMAT_OGR[value['type']])  # Init field
                #logger.debug("> Create attribute %s of type %s" % (key, my_var.FORMAT_OGR_STR[value['type']]))
                
                # Set width
                if value['type'] == "int4":
                    tmp_field.SetWidth(4)
                elif value['type'] in ["int", "int9"]:
                    tmp_field.SetWidth(9)
                elif "width" in value.keys():
                    tmp_field.SetWidth(value['width'])
                
                # Set precision
                if "precision" in value.keys():
                    tmp_field.SetPrecision(value['precision'])
                
            else:
                logger.error("Format %s unknown for shapefile => attribute %s CANNOT be created" % (value["type"], key))
                raise
            
            self.layer.CreateField(tmp_field)

    #----------------------------------------
    
    def add_feature(self, in_geom, in_attributes):
        """
        Add feature with geometry given by in_geom and attribute values given by in_attributes
        
        :param in_geom: geometry to add to feature
        :type in_geom: OGRGeometry
        :param in_attributes: list of attributes and their value 
        :type in_attributes: dict
        """
        
        # 1 - Create the object
        feature = ogr.Feature(self.layer_defn)
        
        # 2 - Add geometry
        feature.SetGeometry(in_geom)    
        
        # 3 - Add attribute values
        feature = self.fill_feature_from_dict(feature, in_attributes)
            
        # 4 - Add feature
        self.layer.CreateFeature(feature)
        
        # 5 - Destroy the feature to free resources
        feature.Destroy()

    def fill_feature_from_dict(self, in_feature, in_dict):
        """
        Fill feature from dictionary information. Attributes and type are loaded from self.attribute_metadata.
    
        :param in_feature: empty feature
        :type in_feature: ogr.Feature
        :param in_dict: attributes information to store in feature.
        :type in_dict: dict
    
        :return: in_feature = feature filled
        :rtype: in_feature = ogr.Feature
        """
                
        att_keys = in_dict.keys()  # List of attributes to value
        
        for att_name, attr_carac in self.attribute_metadata.items():  # Loop over the whole list of the existing attributes
            value_to_write = None  # Init value to write
            if att_name in att_keys:  # If in the list of attributes to modify
                if in_dict[att_name] is None:
                    value_to_write = my_var.FV_SHP[attr_carac["type"]]  # Fill to _FillValue if None
                else:
                    if in_dict[att_name] != my_var.FV_SHP[attr_carac["type"]]:
                        in_dict[att_name] = self.check_validity_value(in_dict[att_name], att_name)
                    value_to_write = convert_wrt_type(in_dict[att_name], attr_carac["type"])  # Fill with appropriate type
            else:
                value_to_write = my_var.FV_SHP[attr_carac["type"]]  # Fill to _FillValue if not in list
    
            in_feature.SetField(str(att_name), value_to_write)  # Fill the field with the wanted value
    
        return in_feature

    def fill_dict_from_feature(self, in_feature):
        """
        Fill dictionary from in_feature information. Attributes and type are loaded from self.attribute_metadata.
    
        :param in_feature: feature with information
        :type in_feature: ogr.Feature
    
        :return: out_dict = fill dictionary
        :rtype: out_dict = depends on the wanted type
        """
        
        out_dict = {}
        for key, value in self.attribute_metadata.items():
            val = in_feature.GetField(key)
            if val == my_var.FV_SHP[self.attribute_metadata[key]['type']]:
                out_dict[key] = None
            else:
                out_dict[key] = val
    
        return out_dict

    #----------------------------------------

    def check_validity_value(self, in_value, att_name):
        """
        Check values validity, replace value by fill value if under valid_min ou over valid_max.

        :param in_value: value
        :type in_value: int, string, char.
        :param in_att_name: name of attribute value
        :type in_att_name: string

        :return: out_value = valid values
        :rtype: int, string, char.
        """

        # 0 - Init output array
        out_value = in_value

        # 1 - Test if a value valid_min exists in attribute_metadata
        if "valid_min" in self.attribute_metadata[att_name] and self.attribute_metadata[att_name]["valid_min"] != "TBD":
            # Replace values under valid_min by no data value
            if in_value < float(self.attribute_metadata[att_name]["valid_min"]):
                out_value = my_var.FV_SHP[self.attribute_metadata[att_name]["type"]]
                logging.debug("Value %s is under min value for attribute %s" % (str(in_value), att_name))

        # 2 - Test if a value valid_max exists in attribute_metadata
        if "valid_max" in self.attribute_metadata[att_name] and self.attribute_metadata[att_name]["valid_max"] != "TBD":
            # Replace values over valid_max by no data value
            if in_value > float(self.attribute_metadata[att_name]["valid_max"]):
                out_value = my_var.FV_SHP[self.attribute_metadata[att_name]["type"]]
                logging.debug("Value %s is over max value for attribute %s" % (str(in_value), att_name))

        return out_value
    #----------------------------------------

    def update_feature(self, in_id, in_geom=None, in_attributes=None):
        """
        Update feature with input identifier with geometry given by in_geom and attribute values given by in_attributes
        
        :param in_id: identifier of the feature to update
        :type in_id: string
        :param in_geom: geometry to set to feature (may be None)
        :type in_geom: OGRGeometry
        :param in_attributes: list of attributes and their value (may be None)
        :type in_attributes: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # 1 - Retrieve the feature
        self.layer.SetAttributeFilter("%s = '%s'" % (self.primary_key, in_id))
        
        nb_features = self.layer.GetFeatureCount()
        if nb_features == 0:
            logger.debug("WARNING There is no feature init with lake_id = %s" % in_id)
            
        elif nb_features > 1:
            logger.debug("WARNING There is %d features init with lake_id = %s" % (nb_features, in_id))
        
        else:
            feature = self.layer.GetNextFeature()
            
            # 2 - Update geometry
            if in_geom is not None:
                feature.SetGeometry(in_geom)    
            
            # 3 - Update attribute values
            if in_attributes is not None:
                feature = self.update_fields(feature, in_attributes)

        # 4 - Add feature
        self.layer.SetFeature(feature)
        
        # 5 - Reset filter
        self.layer.SetAttributeFilter(None)
        
        # 6 - Destroy the feature to free resources
        feature.Destroy()
        
    def update_fields(self, in_feature, in_dict):
        """
        Update fields of input feature, listed in input dictionary. Attributes and type are loaded from self.attribute_metadata.
    
        :param in_feature: feature to update
        :type in_feature: ogr.Feature
        :param in_dict: attributes fields to update.
        :type in_dict: dict
    
        :return: in_feature = updated feature
        :rtype: in_feature = ogr.Feature
        """
        
        att_keys = self.attribute_metadata.keys()  # List of feature attributes
        
        for att_name in in_dict.keys():
            if att_name in att_keys:  # If in the list of attributes
                if in_dict[att_name] is not None: 
                    # Convert the value to the appropriate format
                    if in_dict[att_name] != my_var.FV_SHP[self.attribute_metadata[att_name]["type"]]:
                        in_dict[att_name] = self.check_validity_value(in_dict[att_name], att_name)
                    value_to_write = convert_wrt_type(in_dict[att_name],
                                                      self.attribute_metadata[att_name]["type"])  
                    # Update the field with the wanted value
                    in_feature.SetField(str(att_name), value_to_write)
    
        return in_feature

    #----------------------------------------
        
    def update_and_write_metadata(self, in_filename, in_inprod_metadata=None, in_proc_metadata=None):
        """
        Write .shp.xml file
        
        :param in_filename: full path for the metadata file
        :type in_filename: string   
        :param in_inprod_metadata: metadata specific to input data product
        :type in_inprod_metadata: dict
        :param in_proc_metadata: metadata specific to processing
        :type in_proc_metadata: dict
        """
        #log logger = logging.getLogger(self.__class__.__name__)
        #log logger.debug("- start -")
        
        # 1.1 - Update metadata retrieved from input data product
        if in_inprod_metadata is not None:
            self.set_metadata_val(in_inprod_metadata)
        # 1.2 - Update processing metadata
        if in_proc_metadata is not None:
            self.set_metadata_val(in_proc_metadata)
            
        # 2 - Update file path if necessary
        self.fullpath_or_basename()

        # 3 - Build the XML tree
        self.build_xml()
        
        # 4 - Write the .shp.xml file
        self.xml_tree.write(in_filename, pretty_print=True, xml_declaration=True, encoding='utf-8')
    
    def build_xml(self):
        """
        Write the metadata file associated to the lake shapefile
        
        :param in_filename: full path for the metadata file
        :type in_filename: string
        """
        
        # 1 - Define root element
        root = ET.Element("swot_product")
        
        # 2 - Add global metadata
        sub1 = ET.SubElement(root, "global_metadata")
        for name, dict_att in self.global_metadata.items():
            ET.SubElement(sub1, str(name)).text = str(dict_att["value"])
        
        # 3 - Add attributes metadata
        sub2 = ET.SubElement(root, "attribute_metadata")
        for name, dict_att in self.attribute_metadata.items():  # 1st sub-level
            sub_att = ET.SubElement(sub2, str(name))
            for key, value in dict_att.items():  # 2nd sub-level
                if key not in ["width", "precision"]:
                    ET.SubElement(sub_att, str(key)).text = str(value)
                
        # 4 - Add config params metadata in case of LakeTile product
        if "laketile" in self.file_type:
            sub3 = ET.SubElement(root, "processing_parameters")
            for key, value in self.config_metadata.items():
                ET.SubElement(sub3, str(key)).text = str(value)

        # 5 - Save the XML tree
        self.xml_tree = ET.ElementTree(root)


#######################################
        

class LakeTileShpProduct(ShapefileProduct):
    """
    Deal with LakeTile shapefiles
    """
    
    def __init__(self, in_xml_file, in_layer_name, in_filename=None, in_pixc_metadata=None, in_proc_metadata=None):
        """
        Constructor; specific to LakeTile shapefiles
        
        :param in_xml_file: XML file to populate the variables of the object
        :type in_xml_file: string
        :param in_layer_name: name for shapefile layer        
        :type in_layer_name: string
        :param in_filename: full path of the shapefile for the layer; if None, create a memory layer (default)
        :param in_filename: string     
        :param in_pixc_metadata: metadata specific to L2_HR_PIXC product
        :type in_pixc_metadata: dict
        :param in_proc_metadata: metadata specific to LakeTile processing
        :type in_proc_metadata: dict
        """
        # Get instance of service config file
        cfg = service_config_file.get_instance()
        
        # 1 - Inheritance of ShapefileProduct
        super().__init__(in_xml_file,
                         in_layer_name, 
                         in_filename=in_filename, 
                         in_inprod_metadata=in_pixc_metadata,
                         in_proc_metadata=in_proc_metadata)
        
        # 2 - Init configuration parameters metadata
        self.config_metadata = OrderedDict()
        self.config_metadata["lake_db_id"] = str(cfg.get('DATABASES', 'LAKE_DB_ID'))
        self.config_metadata["classif_list"] = cfg.get('CONFIG_PARAMS', 'CLASSIF_LIST')
        self.config_metadata["classif_mandatory"] = cfg.get('CONFIG_PARAMS', 'CLASSIF_MANDATORY')
        self.config_metadata["classif_water"] = cfg.get('CONFIG_PARAMS', 'CLASSIF_WATER')
        self.config_metadata["classif_4hull"] = cfg.get('CONFIG_PARAMS', 'CLASSIF_4HULL')
        self.config_metadata["del_1pix_in_az_4hull"] = cfg.get('CONFIG_PARAMS', 'DEL_1PIX_IN_AZ_4HULL')
        self.config_metadata["classif_4wse"] = cfg.get('CONFIG_PARAMS', 'CLASSIF_4WSE')
        self.config_metadata["classif_4area"] = cfg.get('CONFIG_PARAMS', 'CLASSIF_4AREA')
        self.config_metadata["del_1pix_in_az_4area"] = cfg.get('CONFIG_PARAMS', 'DEL_1PIX_IN_AZ_4AREA')
        self.config_metadata["min_xtrack"] = cfg.get('CONFIG_PARAMS', 'MIN_XTRACK')
        self.config_metadata["max_xtrack"] = cfg.get('CONFIG_PARAMS', 'MAX_XTRACK')
        self.config_metadata["min_nb_pixc"] = cfg.get('CONFIG_PARAMS', 'MIN_NB_PIXC')
        self.config_metadata["min_size"] = cfg.get('CONFIG_PARAMS', 'MIN_SIZE')
        self.config_metadata["max_nb_tiles_full_az"] = cfg.get('CONFIG_PARAMS', 'MAX_NB_TILES_FULL_AZ')
        self.config_metadata["assignment_layer"] = cfg.get('CONFIG_PARAMS', 'ASSIGNMENT_LAYER')
        self.config_metadata["min_overlap"] = cfg.get('CONFIG_PARAMS', 'MIN_OVERLAP')
        self.config_metadata["area_method"] = cfg.get('CONFIG_PARAMS', 'AREA_METHOD')
        self.config_metadata["segmentation_method"] = cfg.get('CONFIG_PARAMS', 'SEGMENTATION_METHOD')
        self.config_metadata["imp_geoloc"] = cfg.get('CONFIG_PARAMS', 'IMP_GEOLOC')
        self.config_metadata["max_diff_height_geoloc"] = cfg.get('CONFIG_PARAMS', 'MAX_DIFF_HEIGHT_GEOLOC')
        self.config_metadata["hull_method"] = cfg.get('CONFIG_PARAMS', 'HULL_METHOD')
        self.config_metadata["nb_pix_max_delauney"] = cfg.get('CONFIG_PARAMS', 'NB_PIX_MAX_DELAUNEY')
        self.config_metadata["biglake_model"] = cfg.get('CONFIG_PARAMS', 'BIGLAKE_MODEL')
        self.config_metadata["biglake_min_size"] = cfg.get('CONFIG_PARAMS', 'BIGLAKE_MIN_SIZE')
        self.config_metadata["biglake_grid_spacing"] = cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_SPACING')
        self.config_metadata["biglake_grid_res"] = cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_RES')
        self.config_metadata["stocc_input"] = cfg.get('CONFIG_PARAMS', 'STOCC_INPUT')
        self.config_metadata["window_4_ice_clim"] = cfg.get('CONFIG_PARAMS', 'WINDOW_4_ICE_CLIM')
        self.config_metadata["add_all"] = cfg.get('CONFIG_PARAMS', 'ADD_ALL')
        self.config_metadata["nb_digits"] = cfg.get('CONFIG_PARAMS', 'NB_DIGITS')


class LakeTileObsProduct(LakeTileShpProduct):
    """
    Deal with LakeTile_Obs file, the obs-oriented shapefile
    """
    
    def __init__(self, in_layer_name, in_filename=None, in_pixc_metadata=None, in_proc_metadata=None):
        """
        Constructor; specific to LakeTile_Obs shapefile
        
        :param in_layer_name: name for shapefile layer        
        :type in_layer_name: string
        :param in_filename: full path of the shapefile for the layer; if None, create a memory layer (default)
        :param in_filename: string     
        :param in_pixc_metadata: metadata specific to L2_HR_PIXC product
        :type in_pixc_metadata: dict
        :param in_proc_metadata: metadata specific to LakeTile processing
        :type in_proc_metadata: dict
        """
        
        # 1 - Inheritance of LakeTileShpProduct
        super().__init__(os.path.join(os.path.dirname( __file__ ), "xml/laketile_obs.xml"),
                          in_layer_name, 
                          in_filename=in_filename, 
                          in_pixc_metadata=in_pixc_metadata,
                          in_proc_metadata=in_proc_metadata)


class LakeTilePriorProduct(LakeTileShpProduct):
    """
    Deal with LakeTile_Prior file, the PLD-oriented shapefile
    """
    
    def __init__(self, in_layer_name, in_filename=None, in_pixc_metadata=None, in_proc_metadata=None):
        """
        Constructor; specific to LakeTile_Prior shapefile
        
        :param in_layer_name: name for shapefile layer        
        :type in_layer_name: string
        :param in_filename: full path of the shapefile for the layer; if None, create a memory layer (default)
        :param in_filename: string     
        :param in_pixc_metadata: metadata specific to L2_HR_PIXC product
        :type in_pixc_metadata: dict
        :param in_proc_metadata: metadata specific to LakeTile processing
        :type in_proc_metadata: dict
        """
        
        # 1 - Inheritance of LakeTileShpProduct
        super().__init__(os.path.join(os.path.dirname( __file__ ), "xml/laketile_prior.xml"),
                          in_layer_name, 
                          in_filename=in_filename, 
                          in_pixc_metadata=in_pixc_metadata,
                          in_proc_metadata=in_proc_metadata)


class LakeTileUnassignedProduct(LakeTileShpProduct):
    """
    Deal with LakeTile_Unassigned file, the shapefile dedicated to unassigned water features
    """
    
    def __init__(self, in_layer_name, in_filename=None, in_pixc_metadata=None, in_proc_metadata=None):
        """
        Constructor; specific to LakeTile_Unassigned shapefile
        
        :param in_layer_name: name for shapefile layer        
        :type in_layer_name: string
        :param in_filename: full path of the shapefile for the layer; if None, create a memory layer (default)
        :param in_filename: string     
        :param in_pixc_metadata: metadata specific to L2_HR_PIXC product
        :type in_pixc_metadata: dict
        :param in_proc_metadata: metadata specific to LakeTile processing
        :type in_proc_metadata: dict
        """
        
        # 1 - Inheritance of LakeTileShpProduct
        super().__init__(os.path.join(os.path.dirname( __file__ ), "xml/laketile_unassigned.xml"),
                          in_layer_name, 
                          in_filename=in_filename, 
                          in_pixc_metadata=in_pixc_metadata,
                          in_proc_metadata=in_proc_metadata)


#######################################


class LakeSPObsProduct(ShapefileProduct):
    """
    Deal with LakeSP_Obs file, the obs-oriented shapefile
    """
    
    def __init__(self, in_layer_name, in_filename=None, in_laketile_metadata=None, in_proc_metadata=None, flag_create=True):
        """
        Constructor; specific to LakeSP_Obs shapefile
        
        :param in_layer_name: name for shapefile layer        
        :type in_layer_name: string
        :param in_filename: full path of the shapefile for the layer; if None, create a memory layer (default)
        :param in_filename: string     
        :param in_laketile_metadata: metadata specific to L2_HR_LakeTile product
        :type in_laketile_metadata: dict
        :param in_proc_metadata: metadata specific to LakeSP processing
        :type in_proc_metadata: dict
        :param flag_create: =True (default) to create a shapefile or a memory layer, =False otherwise
        :type flag_create: boolean
        """
        
        # 1 - Inheritance of ShapefileProduct
        super().__init__(os.path.join(os.path.dirname( __file__ ), "xml/lakesp_obs.xml"),
                          in_layer_name, 
                          in_filename=in_filename, 
                          in_inprod_metadata=in_laketile_metadata,
                          in_proc_metadata=in_proc_metadata,
                          flag_create=flag_create)


class LakeSPPriorProduct(ShapefileProduct):
    """
    Deal with LakeSP_Prior file, the PLD-oriented shapefile
    """
    
    def __init__(self, in_layer_name, in_filename=None, in_laketile_metadata=None, in_proc_metadata=None, flag_create=True):
        """
        Constructor; specific to LakeSP_PLD shapefile
        
        :param in_layer_name: name for shapefile layer        
        :type in_layer_name: string
        :param in_filename: full path of the shapefile for the layer; if None, create a memory layer (default)
        :param in_filename: string     
        :param in_laketile_metadata: metadata specific to L2_HR_LakeTile product
        :type in_laketile_metadata: dict
        :param in_proc_metadata: metadata specific to LakeSP processing
        :type in_proc_metadata: dict
        :param flag_create: =True (default) to create a shapefile or a memory layer, =False otherwise
        :type flag_create: boolean
        """
        
        # 1 - Inheritance of ShapefileProduct
        super().__init__(os.path.join(os.path.dirname( __file__ ), "xml/lakesp_prior.xml"),
                          in_layer_name, 
                          in_filename=in_filename, 
                          in_inprod_metadata=in_laketile_metadata,
                          in_proc_metadata=in_proc_metadata,
                          flag_create=flag_create)

    #----------------------------------------

    def merge_duplicate_features(self, lake_sp_lyr):
        """
        A same PLD lake can be observed in swath R and L. In this case, 
        the lake is duplicated in LakeSP_Prior file. A merge
        between the 2 observations should be conducted.

        :param lake_sp_lyr: lake_sp_prior layer
        :type lake_sp_lyr: OGR.Layer
        """
        logger = logging.getLogger(self.__class__.__name__)

        # 1 - Get lake_id of duplicated lakes
        # 1.1 - Retrieve lake_id
        lake_id_list = []
        for feature in lake_sp_lyr:
            lake_id_list.append(feature.GetField("lake_id"))
        lake_sp_lyr.ResetReading()
        # 1.2 - List duplicate lake_id
        tmp_duplicates = [lake_id for lake_id in lake_id_list if lake_id_list.count(lake_id) > 1]
        duplicates = set(tmp_duplicates)
        logger.debug("%d lakes are in swath R and L, should be merged" % len(duplicates))
        # 1.3 - Gather duplicate feature FID
        duplicate_lake_id = []
        duplicate_lake_id_fid = []
        for feature in lake_sp_lyr:
            lake_id = feature.GetField("lake_id")
            if lake_id in duplicates:
                lake_fid = feature.GetFID()
                if lake_id not in duplicate_lake_id:
                    duplicate_lake_id.append(lake_id)
                    duplicate_lake_id_fid.append([lake_fid])
                else:
                    duplicate_lake_id_fid[duplicate_lake_id.index(lake_id)].append(lake_fid)
        logger.debug("List of lake_id to merge: %s" % (";".join(duplicate_lake_id)))

        # 2 - Merge duplicated lakes
        fid_to_delete = []
        # 2.1 - Merge duplicate features into single features
        for lake_id, lake_fids in zip(duplicate_lake_id, duplicate_lake_id_fid):
            fid_to_delete += self.merge_lakes(lake_id, lake_fids, lake_sp_lyr)
        # 2.2 - Delete original features
        fid_to_delete.sort(reverse=True)
        for fid in fid_to_delete:
            lake_sp_lyr.DeleteFeature(fid)

    def merge_lakes(self, lake_id, lake_fids, lake_sp_lyr):
        """
        PLD lake lake_id have been observed on swath R and L, corresponding 
        features are merged into a single feature and deleted.

        :param lake_id: lake_id of lake to merge
        :type lake_id: str
        :param lake_fids: list of FID to merge
        :type lake_fids: list of int
        :param lake_sp_lyr: lake_sp_prior layer
        :type lake_sp_lyr: OGR.Layer
        """
        #log logger = logging.getLogger(self.__class__.__name__)
        #log logger.debug("Merge lake %s " % lake_id)
        
        #log if len(lake_fids) > 2:
            #log logger.debug("WARNING !!! LakeSP %s contains more than 2 lakes => all features will be merged" % lake_id)
            
        fid_to_delete = []
        while len(lake_fids) >= 2:
            
            # 1 - Retrieve first 2 features and their obs_id
            lake_fid_1 = lake_fids[0]
            lake_fid_2 = lake_fids[1]
            feat_lake_1 = lake_sp_lyr.GetFeature(lake_fid_1)
            feat_lake_2 = lake_sp_lyr.GetFeature(lake_fid_2)
            obs_id_1 = feat_lake_1.GetField("obs_id")
            obs_id_2 = feat_lake_2.GetField("obs_id")
            time_str_1 = feat_lake_1.GetField("time_str")
            time_str_2 = feat_lake_2.GetField("time_str")

            if obs_id_1 == 'no_data' and obs_id_2 == 'no_data':
                #logger.debug("Both features empty")
                fid_to_delete.append(lake_fid_2)
                lake_fids.pop(1)
                
            elif obs_id_1 == 'no_data':
                #logger.debug("Feature 1 empty")
                fid_to_delete.append(lake_fid_1)
                lake_fids.pop(0)
                
            elif obs_id_2 == 'no_data':
                #logger.debug("Feature 2 empty")
                fid_to_delete.append(lake_fid_2)
                lake_fids.pop(1)
                
            elif (time_str_1 == "patch") or (time_str_2 == "patch"):
                fid_to_delete.append(lake_fid_1)
                fid_to_delete.append(lake_fid_2)
                
                # 2 - Merge geometries
                geom_merged = self.merge_geometries(feat_lake_1, feat_lake_2)
                
                # 3 - Compute merged feature from Feature which is not a patch
                if time_str_1 == "patch":
                    #logger.debug("Feature 1 = patch for big lakes")
                    feat_lake_merged = feat_lake_2.Clone()
                else:
                    #logger.debug("Feature 2 = patch for big lakes")
                    feat_lake_merged = feat_lake_1.Clone()
                feat_lake_merged.SetGeometry(geom_merged)
                
                # 4 - Add merged feature to layer
                lake_sp_lyr.CreateFeature(feat_lake_merged)
                
                # 5 - Add new FID to list
                lake_fids.append(feat_lake_merged.GetFID())
                
                # 6 - Remove old FIDs
                lake_fids.pop(0)
                lake_fids.pop(0)
                
            else:
                #logger.debug("Both features filled")
                fid_to_delete.append(lake_fid_1)
                fid_to_delete.append(lake_fid_2)
                
                # 2 - Merge features
                feat_lake_merged = self.merge_features(feat_lake_1, feat_lake_2)
                
                # 3 - Add merged feature to layer
                lake_sp_lyr.CreateFeature(feat_lake_merged)
                
                # 4 - Add new FID to list
                lake_fids.append(feat_lake_merged.GetFID())
                
                # 5 - Remove old FIDs
                lake_fids.pop(0)
                lake_fids.pop(0)

        return fid_to_delete
    
    def merge_geometries(self, feat_lake_1, feat_lake_2):
        """
        Merge geometres of features 1 and 2 into a single one.

        :param feat_lake_1: lake feature 1
        :type feat_lake_1: ogr.Feature
        :param feat_lake_2: lake feature 2
        :type feat_lake_2: ogr.Feature

        :return: geom_merged = merged geometry
        :rtype: geom_merged = ogr.Feature
        """
        
        if feat_lake_1.GetGeometryRef():
            geom_lake_1 = feat_lake_1.GetGeometryRef()
        else:
            geom_lake_1 = ogr.Geometry(ogr.wkbPolygon)
        if feat_lake_2.GetGeometryRef():
            geom_lake_2 = feat_lake_2.GetGeometryRef()
        else:
            geom_lake_2 = ogr.Geometry(ogr.wkbPolygon)
            
        geom_merged = geom_lake_1.Union(geom_lake_2)
        
        return geom_merged

    def merge_features(self, feat_lake_1, feat_lake_2):
        """
        Merge feature 1 and feature 2 into a single one.

        :param feat_lake_1: lake feature 1
        :type feat_lake_1: ogr.Feature
        :param feat_lake_2: lake feature 2
        :type feat_lake_2: ogr.Feature

        :return: feat_lake_merged = merged lake feature
        :rtype: feat_lake_merged = ogr.Feature
        """
        
        # 1 - Merge geometries
        geom_merged = self.merge_geometries(feat_lake_1, feat_lake_2)

        # 2 - Merge attributes
        
        # 2.1 - Retrieve attributes
        lake_1_dict = self.fill_dict_from_feature(feat_lake_1)
        lake_2_dict = self.fill_dict_from_feature(feat_lake_2)

        # 2.2 - Init ratio between 2 observations
        if lake_1_dict['area_total'] and lake_2_dict['area_total']:
            ratio = lake_1_dict['area_total'] / (lake_1_dict['area_total'] + lake_2_dict['area_total'])
        elif lake_1_dict['area_total']:
            ratio = 1
        elif lake_2_dict['area_total']:
            ratio = 0
        else:
            ratio = None

        # 2.3 - Compute merge 
        lake_merged_dict = {}
        
        lake_merged_dict["lake_id"] = lake_1_dict["lake_id"]
        
        if lake_1_dict["reach_id"]:
            lake_merged_dict["reach_id"] = lake_1_dict["reach_id"]
        elif lake_2_dict["reach_id"]:
            lake_merged_dict["reach_id"] = lake_2_dict["reach_id"]
        else:
            lake_merged_dict["reach_id"] = None

        lake_obs_list = []
        overlap_list = []
        n_overlap = 0
        if lake_1_dict["obs_id"]:
            lake_obs_list += lake_1_dict["obs_id"].split(";")
            overlap_list += [elem for elem in lake_1_dict["overlap"].split(";")]
            n_overlap += int(lake_1_dict["n_overlap"])

        if lake_2_dict["obs_id"]:
            lake_obs_list += lake_2_dict["obs_id"].split(";")
            overlap_list += [elem for elem in lake_2_dict["overlap"].split(";")]
            n_overlap += int(lake_2_dict["n_overlap"])

        while "" in overlap_list:
            overlap_list.remove("")
        overlap_list = [int(elem) for elem in overlap_list if elem != '']

        idx_sorted = np.argsort(overlap_list)[::-1]

        lake_obs_list = [str(lake_obs_list[idx]) for i, idx in enumerate(idx_sorted) if idx < len(lake_obs_list)]
        overlap_list = [str(overlap_list[idx]) for i, idx in enumerate(idx_sorted) if idx < len(overlap_list)]

        if lake_obs_list:
            lake_merged_dict["obs_id"] = ";".join(lake_obs_list)[:254] # truncated to shapefile max str length: 254
            lake_merged_dict["overlap"] = ";".join(overlap_list)[:254] # truncated to shapefile max str length: 254
            lake_merged_dict["n_overlap"] = str(n_overlap)
        else:
            lake_merged_dict["obs_id"] = None
            lake_merged_dict["overlap"] = None
            lake_merged_dict["n_overlap"] = None

        lake_merged_dict["time"] = attributes_fusion(lake_1_dict["time"], lake_2_dict["time"], 'ratio', ratio)
        lake_merged_dict["time_tai"] = attributes_fusion(lake_1_dict["time_tai"], lake_2_dict["time_tai"], 'ratio', ratio)
        lake_merged_dict['time_str'] = my_tools.convert_utc_to_str(lake_merged_dict["time"])
        lake_merged_dict['t_prv_xovr'] = attributes_fusion(lake_1_dict["t_prv_xovr"], lake_2_dict["t_prv_xovr"], 'ratio', ratio)
        lake_merged_dict['t_nxt_xovr'] = attributes_fusion(lake_1_dict["t_nxt_xovr"], lake_2_dict["t_nxt_xovr"], 'ratio', ratio)

        lake_merged_dict["wse"] = attributes_fusion(lake_1_dict["wse"], lake_2_dict["wse"], 'ratio', ratio)
        lake_merged_dict["wse_u"] = attributes_fusion(lake_1_dict["wse_u"], lake_2_dict["wse_u"], 'ratio', ratio)
        lake_merged_dict["wse_r_u"] = attributes_fusion(lake_1_dict["wse_r_u"], lake_2_dict["wse_r_u"], 'ratio', ratio)
        lake_merged_dict["wse_std"] = attributes_fusion(lake_1_dict["wse_std"], lake_2_dict["wse_std"], 'ratio', ratio)
        lake_merged_dict["n_wse"] = attributes_fusion(lake_1_dict["n_wse"], lake_2_dict["n_wse"], 'sum', ratio)

        lake_merged_dict["wse2"] = attributes_fusion(lake_1_dict["wse2"], lake_2_dict["wse2"], 'ratio', ratio)
        lake_merged_dict["wse2_u"] = attributes_fusion(lake_1_dict["wse2_u"], lake_2_dict["wse2_u"], 'ratio', ratio)
        lake_merged_dict["wse2_r_u"] = attributes_fusion(lake_1_dict["wse2_r_u"], lake_2_dict["wse2_r_u"], 'ratio', ratio)
        lake_merged_dict["wse2_std"] = attributes_fusion(lake_1_dict["wse2_std"], lake_2_dict["wse2_std"], 'ratio', ratio)
        lake_merged_dict["n_wse2"] = attributes_fusion(lake_1_dict["n_wse2"], lake_2_dict["n_wse2"], 'sum', ratio)

        lake_merged_dict["area_total"] = attributes_fusion(lake_1_dict["area_total"], lake_2_dict["area_total"], 'sum', ratio)
        lake_merged_dict["area_tot_u"] = attributes_fusion(lake_1_dict["area_tot_u"], lake_2_dict["area_tot_u"], 'ratio', ratio)
        lake_merged_dict["area_detct"] = attributes_fusion(lake_1_dict["area_detct"], lake_2_dict["area_detct"], 'sum', ratio)
        lake_merged_dict["area_det_u"] = attributes_fusion(lake_1_dict["area_det_u"], lake_2_dict["area_det_u"], 'ratio', ratio)
        lake_merged_dict["n_area"] = attributes_fusion(lake_1_dict["n_area"], lake_2_dict["n_area"], 'sum', ratio)

        lake_merged_dict["layovr_val"] = attributes_fusion(lake_1_dict["layovr_val"], lake_2_dict["layovr_val"], 'ratio', ratio)
        lake_merged_dict["xtrk_dist"] = attributes_fusion(abs(lake_1_dict["xtrk_dist"]), abs(lake_2_dict["xtrk_dist"]), 'ratio', ratio)
        try:
            if ((ratio > 0.5) and (lake_1_dict["xtrk_dist"] < 0)) or ((ratio < 0.5) and (lake_2_dict["xtrk_dist"] < 0)):
                lake_merged_dict["xtrk_dist"] = -lake_merged_dict["xtrk_dist"]
        except:
            pass

        # storage change
        lake_merged_dict["ds1_l"] = attributes_fusion(lake_1_dict["ds1_l"], lake_2_dict["ds1_l"], 'sum', ratio)
        lake_merged_dict["ds1_l_u"] = attributes_fusion(lake_1_dict["ds1_l_u"], lake_2_dict["ds1_l_u"], 'ratio', ratio)

        lake_merged_dict["ds1_q"] = attributes_fusion(lake_1_dict["ds1_q"], lake_2_dict["ds1_q"], 'sum', ratio)
        lake_merged_dict["ds1_q_u"] = attributes_fusion(lake_1_dict["ds1_q_u"], lake_2_dict["ds1_q_u"], 'ratio', ratio)

        lake_merged_dict["ds2_l"] = attributes_fusion(lake_1_dict["ds2_l"], lake_2_dict["ds2_l"], 'sum', ratio)
        lake_merged_dict["ds2_l_u"] = attributes_fusion(lake_1_dict["ds2_l_u"], lake_2_dict["ds2_l_u"], 'ratio', ratio)

        lake_merged_dict["ds2_q"] = attributes_fusion(lake_1_dict["ds2_q"], lake_2_dict["ds2_q"], 'sum', ratio)
        lake_merged_dict["ds2_q_u"] = attributes_fusion(lake_1_dict["ds2_q_u"], lake_2_dict["ds2_q_u"], 'ratio', ratio)

        lake_merged_dict["quality_f"] = attributes_fusion(lake_1_dict["quality_f"], lake_2_dict["quality_f"], 'max',
                                                          ratio)
        lake_merged_dict["qual_f_b"] = attributes_fusion(lake_1_dict["qual_f_b"], lake_2_dict["qual_f_b"], 'or',
                                                          ratio)
        lake_merged_dict["wse_frac"] = attributes_fusion(lake_1_dict["wse_frac"], lake_2_dict["wse_frac"], 'ratio',
                                                          ratio)
        lake_merged_dict["edge_frac"] = attributes_fusion(lake_1_dict["edge_frac"], lake_2_dict["edge_frac"], 'ratio',
                                                          ratio)
        lake_merged_dict["dark_frac"] = attributes_fusion(lake_1_dict["dark_frac"], lake_2_dict["dark_frac"], 'ratio',
                                                          ratio)
        lake_merged_dict["low_frac"] = attributes_fusion(lake_1_dict["low_frac"], lake_2_dict["low_frac"], 'ratio',
                                                          ratio)
        lake_merged_dict["spec_frac"] = attributes_fusion(lake_1_dict["spec_frac"], lake_2_dict["spec_frac"], 'ratio',
                                                          ratio)

        lake_merged_dict["ice_clim_f"] = attributes_fusion(lake_1_dict["ice_clim_f"], lake_2_dict["ice_clim_f"], 'max',
                                                           ratio)
        lake_merged_dict["ice_dyn_f"] = attributes_fusion(lake_1_dict["ice_dyn_f"], lake_2_dict["ice_dyn_f"], 'max',
                                                          ratio)
        lake_merged_dict["partial_f"] = attributes_fusion(lake_1_dict["partial_f"], lake_2_dict["partial_f"], 'max',
                                                          ratio)
        lake_merged_dict["xovr_cal_q"] = attributes_fusion(lake_1_dict["xovr_cal_q"], lake_2_dict["xovr_cal_q"], 'max',
                                                           ratio)
        
        lake_merged_dict["sig0"] = attributes_fusion(lake_1_dict["sig0"], lake_2_dict["sig0"], 'ratio', ratio)
        lake_merged_dict["sig0_u"] = attributes_fusion(lake_1_dict["sig0_u"], lake_2_dict["sig0_u"], 'ratio', ratio)
        lake_merged_dict["sig0_std"] = attributes_fusion(lake_1_dict["sig0_std"], lake_2_dict["sig0_std"], 'ratio', ratio)
        try:
            if lake_1_dict["pol"] == lake_2_dict["pol"]:
                lake_merged_dict["pol"] = lake_1_dict["pol"]
            else:
                lake_merged_dict["pol"] = "H+V"
        except:
            lake_merged_dict["pol"] = ""
        
        lake_merged_dict["geoid_hght"] = attributes_fusion(lake_1_dict["geoid_hght"], lake_2_dict["geoid_hght"],
                                                           'ratio', ratio)
        lake_merged_dict["delta_mls"] = attributes_fusion(lake_1_dict["delta_mls"], lake_2_dict["delta_mls"],
                                                           'ratio', ratio)
        lake_merged_dict["solid_tide"] = attributes_fusion(lake_1_dict["solid_tide"], lake_2_dict["solid_tide"],
                                                           'ratio', ratio)
        lake_merged_dict["load_tidef"] = attributes_fusion(lake_1_dict["load_tidef"], lake_2_dict["load_tidef"],
                                                           'ratio', ratio)
        lake_merged_dict["load_tideg"] = attributes_fusion(lake_1_dict["load_tideg"], lake_2_dict["load_tideg"],
                                                           'ratio', ratio)
        lake_merged_dict["pole_tide"] = attributes_fusion(lake_1_dict["pole_tide"], lake_2_dict["pole_tide"], 'ratio',
                                                          ratio)

        lake_merged_dict["dry_trop_c"] = attributes_fusion(lake_1_dict["dry_trop_c"], lake_2_dict["dry_trop_c"],
                                                           'ratio', ratio)
        lake_merged_dict["wet_trop_c"] = attributes_fusion(lake_1_dict["wet_trop_c"], lake_2_dict["wet_trop_c"],
                                                           'ratio', ratio)
        lake_merged_dict["iono_c"] = attributes_fusion(lake_1_dict["iono_c"], lake_2_dict["iono_c"], 'ratio', ratio)
        lake_merged_dict["xovr_cal_c"] = attributes_fusion(lake_1_dict["xovr_cal_c"], lake_2_dict["xovr_cal_c"],
                                                           'ratio', ratio)

        # PLD attributes should be identical
        lake_merged_dict["lake_name"] = lake_1_dict["lake_name"]
        lake_merged_dict["p_res_id"] = lake_1_dict["p_res_id"]
        lake_merged_dict["p_lon"] = lake_1_dict["p_lon"]
        lake_merged_dict["p_lat"] = lake_1_dict["p_lat"]
        lake_merged_dict["p_ref_wse"] = lake_1_dict["p_ref_wse"]
        lake_merged_dict["p_ref_area"] = lake_1_dict["p_ref_area"]
        lake_merged_dict["p_date_t0"] = lake_1_dict["p_date_t0"]
        lake_merged_dict["p_ds_t0"] = lake_1_dict["p_ds_t0"]
        lake_merged_dict["p_storage"] = lake_1_dict["p_storage"]

        # 3 - Build feature from merged attributes
        feat_lake_merged = feat_lake_1.Clone()
        feat_lake_merged.SetGeometry(geom_merged)
        feat_lake_merged = self.fill_feature_from_dict(feat_lake_merged, lake_merged_dict)

        return feat_lake_merged


class LakeSPUnassignedProduct(ShapefileProduct):
    """
    Deal with LakeTile_Unassigned file, the shapefile dedicated to unassigned water features
    """
    
    def __init__(self, in_layer_name, in_filename=None, in_laketile_metadata=None, in_proc_metadata=None, flag_create=True):
        """
        Constructor; specific to LakeSP_Unassigned shapefile
        
        :param in_layer_name: name for shapefile layer        
        :type in_layer_name: string
        :param in_filename: full path of the shapefile for the layer; if None, create a memory layer (default)
        :param in_filename: string     
        :param in_laketile_metadata: metadata specific to L2_HR_LakeTile product
        :type in_laketile_metadata: dict
        :param in_proc_metadata: metadata specific to LakeTile processing
        :type in_proc_metadata: dict
        :param flag_create: =True (default) to create a shapefile or a memory layer, =False otherwise
        :type flag_create: boolean
        """
        
        # 1 - Inheritance of ShapefileProduct
        super().__init__(os.path.join(os.path.dirname( __file__ ), "xml/lakesp_unassigned.xml"),
                          in_layer_name, 
                          in_filename=in_filename, 
                          in_inprod_metadata=in_laketile_metadata,
                          in_proc_metadata=in_proc_metadata,
                          flag_create=flag_create)
            

#######################################
        
        
class LakeAvgProduct(ShapefileProduct):
    """
    Deal with LakeAvg shapefile
    """
    
    def __init__(self, in_layer_name, in_filename=None, in_proc_metadata=None, flag_create=True):
        """
        Constructor; specific to LakeAvg shapefile
        
        :param in_layer_name: name for shapefile layer        
        :type in_layer_name: string
        :param in_filename: full path of the shapefile for the layer; if None, create a memory layer (default)
        :param in_filename: string   
        :param in_proc_metadata: metadata specific to LakeAvg processing
        :type in_proc_metadata: dict
        :param flag_create: =True (default) to create a shapefile or a memory layer, =False otherwise
        :type flag_create: boolean
        """
        
        # 1 - Inheritance of ShapefileProduct
        super().__init__(os.path.join(os.path.dirname( __file__ ), "xml/lakeavg.xml"),
                          in_layer_name, 
                          in_filename=in_filename,
                          in_proc_metadata=in_proc_metadata,
                          flag_create=flag_create)


#######################################
        
   
def convert_str_to_type(in_val, in_type):
    """
    Convert string in the format specified in input
    
    :param in_val: string value to convert
    :type in_val: string
    :param in_type: type of the wanted output, among "integer", "float", "string"
    :type in_type: string
    
    :return: retour = converted string
    :rtype: retour = depends on the wanted type
    """
    
    if in_type == "integer":
        retour = int(in_val)
    elif in_type == "float":
        retour = float(in_val)
    else:
        retour = str(in_val)
    return retour
        
   
def convert_wrt_type(in_val, in_type):
    """
    Convert in_val in the format specified in input
    
    :param in_val: value to convert
    :type in_val: depends on the value
    :param in_type: type of the wanted output (default=str)
    :type in_type: type OGR
    
    :return: retour = converted value
    :rtype: retour = depends on the wanted type
    """
    
    if in_type == ogr.OFTInteger:
        retour = int(in_val)
    elif in_type == ogr.OFTReal:
        retour = float(in_val)
    else:
        retour = str(in_val)
    return retour
            
    
#######################################
    

def check_lake_db(in_xml_tree, in_lakedb_file):
    """
    Check lake_db filename: value in LakeTile XML files and in LakeSP command file should be the same
    
    :param IN_xml_tree: XML tree from LakeTile shapefile
    :type IN_xml_tree: etree.parse
    :param in_lakedb_file: lake database filename
    :type in_lakedb_file: string
    """ 
    # Get instance of service config file
    cfg = service_config_file.get_instance()
    flag_write_full_path = cfg.getboolean("OPTIONS", "Write full path")
    # Get logger
    logger = logging.getLogger("locnes_products_shapefile")
    logger.debug("Compare LAKE_DB value in LakeTile XML file to value in LakeSP command file")
    
    # Read PLD filename in the XML file
    try:
        lake_db_xml = in_xml_tree.xpath("//swot_product/global_metadata/xref_prior_lake_db_file")[0].text
    except(IndexError):  # Default is None
        lake_db_xml = None
        logger.error("NO lake_db information in LakeTile XML file")
    
    # Both should be the same
    if (in_lakedb_file is None) and (lake_db_xml == None):
        logger.warning('LAKE_DB not filled => LakeSP product not linked to prior lake database')
    else:
        # Compare full path or basename of PLD
        if flag_write_full_path:
            tmp_compare_cfg = in_lakedb_file
            tmp_compare_xml = lake_db_xml
        else:
            if in_lakedb_file is None:
                tmp_compare_cfg = None
            else:
                tmp_compare_cfg = os.path.basename(in_lakedb_file)
            if lake_db_xml is None:
                tmp_compare_xml = None
            else:
                tmp_compare_xml = os.path.basename(lake_db_xml)
            
        # Test
        if (in_lakedb_file is None) and (tmp_compare_xml is not None):
            message = "ERROR bad lake_db. In XML file: " + lake_db_xml + " WHEREAS in command file: None"
            raise service_error.ProcessingError(message, logger)
        elif (in_lakedb_file is not None) and (tmp_compare_xml is None):
            message = "ERROR bad lake_db. In command file: " + in_lakedb_file + " WHEREAS in XML file: None"
            raise service_error.ProcessingError(message, logger)
        elif (tmp_compare_cfg != tmp_compare_xml):
            pass
            # Suppression erreur suite au changement de PLD à LakeDatabase --> a modifier apres
            #message = "ERROR bad lake_db. In XML file: " + tmp_compare_xml + " WHEREAS in command file: " + tmp_compare_cfg
            #raise service_error.ProcessingError(message, logger)
            
    # Print PLD filename in log
    #logger.debug('LAKE_DB = ' + str(in_lakedb_file))


def check_lake_db_id(in_xml_tree):
    """
    Check lake_db_id value: value in LakeTile XML files and in LakeSP command file should be the same
    
    :param IN_xml_tree: XML tree
    :type IN_xml_tree: etree.parse
    """ 
    # Get instance of service config file
    cfg = service_config_file.get_instance()
    # Get logger
    logger = logging.getLogger("locnes_products_shapefile")
    #logger.debug("Compare LAKE_DB_ID value in LakeTile XML file to value in LakeSP command file")

    # Read lake identifier in LakeTile XML file
    try:
        lake_db_id_xml = in_xml_tree.xpath("//swot_product/processing_parameters/lake_db_id")[0].text
    except(IndexError):
        lake_db_id_xml = None
        #logger.error("NO lake_db_id information in LakeTile XML file")
    
    # Retrieve lake identifier from command file
    lake_db_id = cfg.get('DATABASES', 'LAKE_DB_ID')
    
    # Both should be the same
    if (lake_db_id is None) or (lake_db_id_xml == "None"):
        #logger.debug('LAKE_DB_ID not filled')
        pass
    else:
        if (lake_db_id != lake_db_id_xml):
            message = "ERROR bad lake_db_id. In XML file: " + lake_db_id + " WHEREAS in command file: " + lake_db_id_xml
            raise service_error.ProcessingError(message, logger)
            
    # Print lake identifier in log
    #logger.debug('LAKE_DB_ID = ' + str(cfg.get('DATABASES', 'LAKE_DB_ID')))
    

def check_mls(in_xml_tree, in_mls_file):
    """
    Check mls filename: value in LakeTile XML files and in LakeSP command file should be the same
    
    :param IN_xml_tree: XML tree from LakeTile shapefile
    :type IN_xml_tree: etree.parse
    :param in_mls_file: mean lake surfaces database filename
    :type in_mls_file: string
    """ 
    # Get instance of service config file
    cfg = service_config_file.get_instance()
    flag_write_full_path = cfg.getboolean("OPTIONS", "Write full path")
    # Get logger
    logger = logging.getLogger("locnes_products_shapefile")
    logger.debug("Compare MLS value in LakeTile XML file to value in LakeSP command file")
    
    # Read MLS filename in the XML file
    try:
        mls_xml = in_xml_tree.xpath("//swot_product/global_metadata/xref_mls_file")[0].text
    except(IndexError):  # Default is None
        mls_xml = None
        logger.error("NO mls information in LakeTile XML file")
    
    # Both should be the same
    if (in_mls_file is None) and (mls_xml == None):
        logger.warning('MLS not filled => NO correction of geoid for large lakes applied')
    else:
        # Compare full path or basename of PLD
        if flag_write_full_path:
            tmp_compare_cfg = in_mls_file
            tmp_compare_xml = mls_xml
        else:
            if in_mls_file is None:
                tmp_compare_cfg = None
            else:
                tmp_compare_cfg = os.path.basename(in_mls_file)
            if mls_xml is None:
                tmp_compare_xml = None
            else:
                tmp_compare_xml = os.path.basename(mls_xml)
            
        # Test
        if (in_mls_file is None) and (tmp_compare_xml is not None):
            message = "ERROR bad mls. In XML file: " + mls_xml + " WHEREAS in command file: None"
            raise service_error.ProcessingError(message, logger)
        elif (in_mls_file is not None) and (tmp_compare_xml is None):
            message = "ERROR bad mls. In command file: " + in_mls_file + " WHEREAS in XML file: None"
            raise service_error.ProcessingError(message, logger)
        elif (tmp_compare_cfg != tmp_compare_xml):
            message = "ERROR bad mls. In XML file: " + tmp_compare_xml + " WHEREAS in command file: " + tmp_compare_cfg
            raise service_error.ProcessingError(message, logger)


def set_config_from_xml(in_xml_tree):
    """
    Set parameter variables from an XML tree
    
    :param IN_xml_tree: XML tree (typically from .shp.xml file)
    :type IN_xml_tree: etree.parse
    """
    # Get instance of service config file
    cfg = service_config_file.get_instance()
    # Get logger
    logger = logging.getLogger("locnes_products_shapefile")
    logger.debug("Read LakeTile XML")

    section = "CONFIG_PARAMS"
    # If section doesn't exist: add and init values
    # If exists: check values
    if section not in cfg.sections():
        
        # Add value if not in service_config_file
        cfg.add_section(section)
        
        # List of classif flags to keep for processing: 2=land_near_water  3=water_near_land   4=open_water  5=dark_water 
        classif_list = in_xml_tree.xpath("//swot_product/processing_parameters/classif_list")[0].text
        cfg.set(section, "CLASSIF_LIST", classif_list)
        logger.debug('CLASSIF_LIST = ' + str(cfg.get('CONFIG_PARAMS', 'CLASSIF_LIST')))
        
        # List of classif flags a lake feature must contain at least
        classif_mandatory = in_xml_tree.xpath("//swot_product/processing_parameters/classif_mandatory")[0].text
        cfg.set(section, "CLASSIF_MANDATORY", classif_mandatory)
        logger.debug('CLASSIF_MANDATORY = ' + str(cfg.get('CONFIG_PARAMS', 'CLASSIF_MANDATORY')))
        
        # True/False list of classif flags which contain observed water (1-to-1 correspondance with CLASSIF_LIST)
        classif_water = in_xml_tree.xpath("//swot_product/processing_parameters/classif_water")[0].text
        cfg.set(section, "CLASSIF_WATER", classif_water)
        logger.debug('CLASSIF_WATER = ' + str(cfg.get('CONFIG_PARAMS', 'CLASSIF_WATER')))
        
        # True/False list of classif flags to use for lake boundary construction (1-to-1 correspondance with CLASSIF_LIST)
        classif_4hull = in_xml_tree.xpath("//swot_product/processing_parameters/classif_4hull")[0].text
        cfg.set(section, "CLASSIF_4HULL", classif_4hull)
        logger.debug('CLASSIF_4HULL = ' + str(cfg.get('CONFIG_PARAMS', 'CLASSIF_4HULL')))
        
        # If CLASSIF_LIST contains 3: =True to remove 1 pixel in top/bottom azimuth to build the lake boundary; =False to keep all pixels (default)
        del_1pix_in_az_4hull = in_xml_tree.xpath("//swot_product/processing_parameters/del_1pix_in_az_4hull")[0].text
        cfg.set(section, "DEL_1PIX_IN_AZ_4HULL", del_1pix_in_az_4hull)
        logger.debug('DEL_1PIX_IN_AZ_4HULL = ' + str(cfg.get('CONFIG_PARAMS', 'DEL_1PIX_IN_AZ_4HULL')))
        
        # True/False list of classif flags to use for WSE computation (1-to-1 correspondance with CLASSIF_LIST)
        classif_4wse = in_xml_tree.xpath("//swot_product/processing_parameters/classif_4wse")[0].text
        cfg.set(section, "CLASSIF_4WSE", classif_4wse)
        logger.debug('CLASSIF_4WSE = ' + str(cfg.get('CONFIG_PARAMS', 'CLASSIF_4WSE')))
        
        # True/False list of classif flags to use for lake area computation (1-to-1 correspondance with CLASSIF_LIST)
        classif_4area = in_xml_tree.xpath("//swot_product/processing_parameters/classif_4area")[0].text
        cfg.set(section, "CLASSIF_4AREA", classif_4area)
        logger.debug('CLASSIF_4AREA = ' + str(cfg.get('CONFIG_PARAMS', 'CLASSIF_4AREA')))
        
        # If CLASSIF_LIST contains 3: =True to remove 1 pixel in top/bottom azimuth to compute the lake area; =False to keep all pixels (default)
        del_1pix_in_az_4area = in_xml_tree.xpath("//swot_product/processing_parameters/del_1pix_in_az_4area")[0].text
        cfg.set(section, "DEL_1PIX_IN_AZ_4AREA", del_1pix_in_az_4area)
        logger.debug('DEL_1PIX_IN_AZ_4AREA = ' + str(cfg.get('CONFIG_PARAMS', 'DEL_1PIX_IN_AZ_4AREA')))
        
        # Min cross-track distance to select PixC
        min_xtrack = float(in_xml_tree.xpath("//swot_product/processing_parameters/min_xtrack")[0].text)
        cfg.set(section, "MIN_XTRACK", min_xtrack)
        logger.debug('MIN_XTRACK = ' + str(cfg.get('CONFIG_PARAMS', 'MIN_XTRACK')))
        
        # Max cross-track distance to select PixC
        max_xtrack = float(in_xml_tree.xpath("//swot_product/processing_parameters/max_xtrack")[0].text)
        cfg.set(section, "MAX_XTRACK", max_xtrack)
        logger.debug('MAX_XTRACK = ' + str(cfg.get('CONFIG_PARAMS', 'MAX_XTRACK')))
        
        # Minimum number of pixels for a water feature to start being processed
        min_nb_pixc = float(in_xml_tree.xpath("//swot_product/processing_parameters/min_nb_pixc")[0].text)
        cfg.set(section, "MIN_NB_PIXC", min_nb_pixc)
        logger.debug('MIN_NB_PIXC = ' + str(cfg.get('CONFIG_PARAMS', 'MIN_NB_PIXC')))
        
        # Min size for a lake to generate a lake product (=polygon + attributes) for it
        min_size = float(in_xml_tree.xpath("//swot_product/processing_parameters/min_size")[0].text)
        cfg.set(section, "MIN_SIZE", min_size)
        logger.debug('MIN_SIZE = ' + str(cfg.get('CONFIG_PARAMS', 'MIN_SIZE')))

        # Maximum number of azimuth tiles covered by a lake processed in LakeSP
        max_nb_tiles_full_az = float(in_xml_tree.xpath("//swot_product/processing_parameters/max_nb_tiles_full_az")[0].text)
        cfg.set(section, "MAX_NB_TILES_FULL_AZ", max_nb_tiles_full_az)
        logger.debug('MAX_NB_TILES_FULL_AZ = ' + str(cfg.get('CONFIG_PARAMS', 'MAX_NB_TILES_FULL_AZ')))

        # PLD table to use to assign pixels to a PLD lake
        assignment_layer = in_xml_tree.xpath("//swot_product/processing_parameters/assignment_layer")[0].text
        cfg.set(section, "ASSIGNMENT_LAYER", assignment_layer)
        logger.debug('ASSIGNMENT_LAYER = ' + str(cfg.get('CONFIG_PARAMS', 'ASSIGNMENT_LAYER')))

        # Min percentage of overlapping area to consider a PLD lake linked to an observed feature
        min_overlap = float(in_xml_tree.xpath("//swot_product/processing_parameters/min_overlap")[0].text)
        cfg.set(section, "MIN_OVERLAP", min_overlap)
        logger.debug('MIN_OVERLAP = ' + str(cfg.get('CONFIG_PARAMS', 'MIN_OVERLAP')))
        
        # Method to compute area_total attribute
        # polygon = area_total is the area of the polygon of the water region
        # pixc = area_total is computed using JPL aggregate function (default)
        area_method = in_xml_tree.xpath("//swot_product/processing_parameters/area_method")[0].text
        cfg.set(section, "AREA_METHOD", area_method)
        logger.debug('AREA_METHOD = ' + str(cfg.get('CONFIG_PARAMS', 'AREA_METHOD')))
        
        # Method to compute height-segmentation
        # 1 = Felzenszwalb
        # 2 = SLIC
        # 3 = Unsupervised thresholding
        # 4 = Watershed segmentation method
        # 5 = k-means clustering method
        # 6 = hierarchical clustering
        # 7 (default) = Otsu thresholding
        # 8 = MeanShift
        segmentation_method = float(in_xml_tree.xpath("//swot_product/processing_parameters/segmentation_method")[0].text)
        cfg.set(section, "SEGMENTATION_METHOD", segmentation_method)
        logger.debug('SEGMENTATION_METHOD = ' + str(cfg.get('CONFIG_PARAMS', 'SEGMENTATION_METHOD')))
            
        # To improve PixC geolocation (=True) or not (=False)
        imp_geoloc = in_xml_tree.xpath("//swot_product/processing_parameters/imp_geoloc")[0].text
        cfg.set(section, "IMP_GEOLOC", imp_geoloc)
        logger.debug('IMP_GEOLOC = ' + str(cfg.get('CONFIG_PARAMS', 'IMP_GEOLOC')))
        # Maximum tolerated difference (m) between improved height and initial height
        max_diff_height_geoloc = in_xml_tree.xpath("//swot_product/processing_parameters/max_diff_height_geoloc")[0].text
        cfg.set(section, "MAX_DIFF_HEIGHT_GEOLOC", max_diff_height_geoloc)
        logger.debug('MAX_DIFF_HEIGHT_GEOLOC = ' + str(cfg.get('CONFIG_PARAMS', 'MAX_DIFF_HEIGHT_GEOLOC')))
        # Method to compute lake boundary or polygon hull
        # 0 = convex hull 
        # 1.0 = concave hull computed in ground geometry, based on Delaunay triangulation - using CGAL library
        # 1.1 = concave hull computed in ground geometry, based on Delaunay triangulation (time consuming)
        # 1.2 = concave hull computed in ground geometry, based on Delaunay triangulation with alpha parameter varying along range
        # 2.0 = edge computed in radar geometry, then converted in ground geometry by reordering lon lat (default)
        # 2.1 = edge computed in radar geometry, then converted in ground geometry by checking each segment
        hull_method = float(in_xml_tree.xpath("//swot_product/processing_parameters/hull_method")[0].text)
        cfg.set(section, "HULL_METHOD", hull_method)
        logger.debug('HULL_METHOD = ' + str(cfg.get('CONFIG_PARAMS', 'HULL_METHOD')))
        
        # Max number of pixels for hull computation method 1
        nb_pix_max_delauney = int(in_xml_tree.xpath("//swot_product/processing_parameters/nb_pix_max_delauney")[0].text)
        cfg.set(section, "NB_PIX_MAX_DELAUNEY", nb_pix_max_delauney)
        logger.debug('NB_PIX_MAX_DELAUNEY = ' + str(cfg.get('CONFIG_PARAMS', 'NB_PIX_MAX_DELAUNEY')))

        # Big lakes parameters for improved geoloc
        # Model to deal with big lake processing
        biglake_model = in_xml_tree.xpath("//swot_product/processing_parameters/biglake_model")[0].text
        cfg.set(section, "BIGLAKE_MODEL", biglake_model)
        logger.debug('BIGLAKE_MODEL = ' + str(cfg.get('CONFIG_PARAMS', 'BIGLAKE_MODEL')))
        # Min size for lake to be considered as big
        biglake_min_size = float(in_xml_tree.xpath("//swot_product/processing_parameters/biglake_min_size")[0].text)
        cfg.set(section, "BIGLAKE_MIN_SIZE", biglake_min_size)        
        logger.debug('BIGLAKE_MIN_SIZE = ' + str(cfg.get('CONFIG_PARAMS', 'BIGLAKE_MIN_SIZE')))
        # Grid spacing for lake height smoothing
        biglake_grid_spacing = int(in_xml_tree.xpath("//swot_product/processing_parameters/biglake_grid_spacing")[0].text)
        cfg.set(section, "BIGLAKE_GRID_SPACING", biglake_grid_spacing)
        logger.debug('BIGLAKE_GRID_SPACING = ' + str(cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_SPACING')))
        # Grid resolution for lake height smoothing
        biglake_grid_res = int(in_xml_tree.xpath("//swot_product/processing_parameters/biglake_grid_res")[0].text)
        cfg.set(section, "BIGLAKE_GRID_RES", biglake_grid_res)
        logger.debug('BIGLAKE_GRID_RES = ' + str(cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_RES')))
        
        # Input data for storage change computation
        stocc_input = in_xml_tree.xpath("//swot_product/processing_parameters/stocc_input")[0].text
        cfg.set(section, "STOCC_INPUT", stocc_input)
        logger.debug('STOCC_INPUT = ' + str(cfg.get('CONFIG_PARAMS', 'STOCC_INPUT')))
        
        # Input data for ice_dyn_f computation
        window_4_ice_clim = in_xml_tree.xpath("//swot_product/processing_parameters/window_4_ice_clim")[0].text
        cfg.set(section, "WINDOW_4_ICE_CLIM", window_4_ice_clim)
        logger.debug('WINDOW_4_ICE_CLIM = ' + str(cfg.get('CONFIG_PARAMS', 'WINDOW_4_ICE_CLIM')))
        
        # To add not-observed PLD lakes as empty features to the LakeSP_Prior shapefile (=True) or not (=False)
        add_all = in_xml_tree.xpath("//swot_product/processing_parameters/add_all")[0].text
        cfg.set(section, "ADD_ALL", add_all)
        logger.debug('ADD_ALL = ' + str(cfg.get('CONFIG_PARAMS', 'ADD_ALL')))
        # Nb digits for counter of lakes, used in the obs_id identifier of each observed lake
        nb_digits = int(in_xml_tree.xpath("//swot_product/processing_parameters/nb_digits")[0].text)
        cfg.set(section, "NB_DIGITS", nb_digits)
        logger.debug('NB_DIGITS = ' + str(cfg.get('CONFIG_PARAMS', 'NB_DIGITS')))

    else:
        
        # Check values
        try:
            
            # List of classif flags to keep for processing: 2=land_near_water  3=water_near_land   4=open_water  5=dark_water 
            classif_list_xml = in_xml_tree.xpath("//swot_product/processing_parameters/classif_list")[0].text
            classif_list = cfg.get(section, "CLASSIF_LIST")
            if (classif_list_xml != classif_list):
                message = "ERROR bad classif_list. In XML file: " + str(classif_list_xml) + " WHEREAS previously set at: " + str(classif_list)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('CLASSIF_LIST = ' + str(cfg.get('CONFIG_PARAMS', 'CLASSIF_LIST')))
            
            # List of classif flags a lake feature must contain at least
            classif_mandatory_xml = in_xml_tree.xpath("//swot_product/processing_parameters/classif_mandatory")[0].text
            classif_mandatory = cfg.get(section, "CLASSIF_MANDATORY")
            if (classif_mandatory_xml != classif_mandatory):
                message = "ERROR bad classif_mandatory. In XML file: " + str(classif_mandatory_xml) + " WHEREAS previously set at: " + str(classif_mandatory)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('CLASSIF_MANDATORY = ' + str(cfg.get('CONFIG_PARAMS', 'CLASSIF_MANDATORY')))
            
            # True/False list of classif flags which contain observed water (1-to-1 correspondance with CLASSIF_LIST)
            classif_water_xml = in_xml_tree.xpath("//swot_product/processing_parameters/classif_water")[0].text
            classif_water = cfg.get(section, "CLASSIF_WATER")
            if (classif_water_xml != classif_water):
                message = "ERROR bad classif_water. In XML file: " + str(classif_water_xml) + " WHEREAS previously set at: " + str(classif_water)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('CLASSIF_WATER = ' + str(cfg.get('CONFIG_PARAMS', 'CLASSIF_WATER')))
            
            # True/False list of classif flags to use for lake boundary construction (1-to-1 correspondance with CLASSIF_LIST)
            classif_4hull_xml = in_xml_tree.xpath("//swot_product/processing_parameters/classif_4hull")[0].text
            classif_4hull = cfg.get(section, "CLASSIF_4HULL")
            if (classif_4hull_xml != classif_4hull):
                message = "ERROR bad classif_4hull. In XML file: " + str(classif_4hull_xml) + " WHEREAS previously set at: " + str(classif_4hull)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('CLASSIF_4HULL = ' + str(cfg.get('CONFIG_PARAMS', 'CLASSIF_4HULL')))
            
            # If CLASSIF_LIST contains 3: =True to remove 1 pixel in top/bottom azimuth to build the lake boundary; =False to keep all pixels (default)
            del_1pix_in_az_4hull_xml = in_xml_tree.xpath("//swot_product/processing_parameters/del_1pix_in_az_4hull")[0].text
            del_1pix_in_az_4hull = cfg.get(section, "DEL_1PIX_IN_AZ_4HULL")
            if (del_1pix_in_az_4hull_xml != del_1pix_in_az_4hull):
                message = "ERROR bad classif_4hull. In XML file: " + str(del_1pix_in_az_4hull_xml) + " WHEREAS previously set at: " + str(del_1pix_in_az_4hull)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('DEL_1PIX_IN_AZ_4HULL = ' + str(cfg.get('CONFIG_PARAMS', 'DEL_1PIX_IN_AZ_4HULL')))
            
            # True/False list of classif flags to use for WSE computation (1-to-1 correspondance with CLASSIF_LIST)
            classif_4wse_xml = in_xml_tree.xpath("//swot_product/processing_parameters/classif_4wse")[0].text
            classif_4wse = cfg.get(section, "CLASSIF_4WSE")
            if (classif_4wse_xml != classif_4wse):
                message = "ERROR bad classif_4wse. In XML file: " + str(classif_4wse_xml) + " WHEREAS previously set at: " + str(classif_4wse)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('CLASSIF_4WSE = ' + str(cfg.get('CONFIG_PARAMS', 'CLASSIF_4WSE')))
            
            # True/False list of classif flags to use for lake area computation (1-to-1 correspondance with CLASSIF_LIST)
            classif_4area_xml = in_xml_tree.xpath("//swot_product/processing_parameters/classif_4area")[0].text
            classif_4area = cfg.get(section, "CLASSIF_4AREA")
            if (classif_4area_xml != classif_4area):
                message = "ERROR bad classif_4area. In XML file: " + str(classif_4area_xml) + " WHEREAS previously set at: " + str(classif_4area)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('CLASSIF_4AREA = ' + str(cfg.get('CONFIG_PARAMS', 'CLASSIF_4AREA')))
            
            # If CLASSIF_LIST contains 3: =True to remove 1 pixel in top/bottom azimuth to compute the lake area; =False to keep all pixels (default)
            del_1pix_in_az_4area_xml = in_xml_tree.xpath("//swot_product/processing_parameters/del_1pix_in_az_4area")[0].text
            del_1pix_in_az_4area = cfg.get(section, "DEL_1PIX_IN_AZ_4AREA")
            if (del_1pix_in_az_4area_xml != del_1pix_in_az_4area):
                message = "ERROR bad classif_4hull. In XML file: " + str(del_1pix_in_az_4area_xml) + " WHEREAS previously set at: " + str(del_1pix_in_az_4area)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('DEL_1PIX_IN_AZ_4AREA = ' + str(cfg.get('CONFIG_PARAMS', 'DEL_1PIX_IN_AZ_4AREA')))
            
            # Min cross-track distance to select PixC
            min_xtrack_xml = float(in_xml_tree.xpath("//swot_product/processing_parameters/min_xtrack")[0].text)
            min_xtrack = cfg.getfloat(section, "MIN_XTRACK")
            if (min_xtrack_xml != min_xtrack):
                message = "ERROR bad min_xtrack. In XML file: " + str(min_xtrack_xml) + " WHEREAS previously set at: " + str(min_xtrack)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('MIN_XTRACK = ' + str(cfg.get('CONFIG_PARAMS', 'MIN_XTRACK')))
            
            # Max cross-track distance to select PixC
            max_xtrack_xml = float(in_xml_tree.xpath("//swot_product/processing_parameters/max_xtrack")[0].text)
            max_xtrack = cfg.getfloat(section, "MAX_XTRACK")
            if (max_xtrack_xml != max_xtrack):
                message = "ERROR bad max_xtrack. In XML file: " + str(max_xtrack_xml) + " WHEREAS previously set at: " + str(max_xtrack)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('MAX_XTRACK = ' + str(cfg.get('CONFIG_PARAMS', 'MAX_XTRACK')))
                        
            # Minimum number of pixels for a water feature to start being processed
            min_nb_pixc_xml = float(in_xml_tree.xpath("//swot_product/processing_parameters/min_nb_pixc")[0].text)
            min_nb_pixc = cfg.getint(section, "MIN_NB_PIXC")
            if (min_nb_pixc_xml != min_nb_pixc):
                message = "ERROR bad min_nb_pixc. In XML file: " + str(min_nb_pixc_xml) + " WHEREAS previously set at: " \
                          + str(min_nb_pixc)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('MIN_NB_PIXC = ' + str(cfg.get('CONFIG_PARAMS', 'MIN_NB_PIXC')))
            
            # Min size for a lake to generate a lake product (=polygon + attributes) for it
            min_size_xml = float(in_xml_tree.xpath("//swot_product/processing_parameters/min_size")[0].text)
            min_size = cfg.getfloat(section, "MIN_SIZE")
            if (min_size_xml != min_size):
                message = "ERROR bad min_size. In XML file: " + str(min_size_xml) + " WHEREAS previously set at: " + str(min_size)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('MIN_SIZE = ' + str(cfg.get('CONFIG_PARAMS', 'MIN_SIZE')))

            # Maximum number of azimuth tiles covered by a lake processed in LakeSP
            max_nb_tiles_full_az_xml = float(in_xml_tree.xpath("//swot_product/processing_parameters/max_nb_tiles_full_az")[0].text)
            max_nb_tiles_full_az = cfg.getfloat(section, "MAX_NB_TILES_FULL_AZ")
            if (max_nb_tiles_full_az_xml != max_nb_tiles_full_az):
                message = "ERROR bad max_nb_tiles_full_az. In XML file: " + str(max_nb_tiles_full_az_xml) + " WHEREAS previously set at: " \
                                                                        + str(max_nb_tiles_full_az)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('MAX_NB_TILES_FULL_AZ = ' + str(cfg.get('CONFIG_PARAMS', 'MAX_NB_TILES_FULL_AZ')))

            # Min percentage of overlapping area to consider a PLD lake linked to an observed feature
            assignment_layer_xml = in_xml_tree.xpath("//swot_product/processing_parameters/assignment_layer")[0].text
            assignment_layer = cfg.get(section, "ASSIGNMENT_LAYER")
            if (assignment_layer_xml != assignment_layer):
                message = "ERROR bad assignment_layer. In XML file: " + str(assignment_layer_xml) + " WHEREAS previously set at: " + str(assignment_layer)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('ASSIGNMENT_LAYER = ' + str(cfg.get('CONFIG_PARAMS', 'ASSIGNMENT_LAYER')))

            # Min percentage of overlapping area to consider a PLD lake linked to an observed feature
            min_overlap_xml = float(in_xml_tree.xpath("//swot_product/processing_parameters/min_overlap")[0].text)
            min_overlap = cfg.getfloat(section, "MIN_OVERLAP")
            if (min_overlap_xml != min_overlap):
                message = "ERROR bad min_overlap. In XML file: " + str(min_overlap_xml) + " WHEREAS previously set at: " + str(min_overlap)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('MIN_OVERLAP = ' + str(cfg.get('CONFIG_PARAMS', 'MIN_OVERLAP')))
            
            # Method to compute area_total attribute
            area_method_xml = in_xml_tree.xpath("//swot_product/processing_parameters/area_method")[0].text
            area_method = cfg.get(section, "AREA_METHOD")
            if (area_method_xml != area_method):
                message = "ERROR bad area_method. In XML file: " + str(area_method_xml) + " WHEREAS previously set at: " + str(area_method)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('AREA_METHOD = ' + str(cfg.get('CONFIG_PARAMS', 'AREA_METHOD')))
            
            # Method to compute height-segmentation
            # 1 = Felzenszwalb
            # 2 = SLIC
            # 3 = Unsupervised thresholding
            # 4 = Watershed segmentation method
            # 5 = k-means clustering method
            # 6 = hierarchical clustering
            # 7 (default) = Otsu thresholding
            # 8 = MeanShift
            segmentation_method_xml = float(in_xml_tree.xpath("//swot_product/processing_parameters/segmentation_method")[0].text)
            segmentation_method = cfg.getint(section, "SEGMENTATION_METHOD")
            if (segmentation_method_xml != segmentation_method):
                message = "ERROR bad segmentation_method. In XML file: " + str(segmentation_method_xml) + " WHEREAS previously set at: " \
                          + str(segmentation_method)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('SEGMENTATION_METHOD = ' + str(cfg.get('CONFIG_PARAMS', 'SEGMENTATION_METHOD')))
    
            # To improve PixC geolocation (=True) or not (=False)
            imp_geoloc_xml = in_xml_tree.xpath("//swot_product/processing_parameters/imp_geoloc")[0].text
            imp_geoloc = cfg.getboolean(section, "IMP_GEOLOC")
            if (imp_geoloc_xml != str(imp_geoloc)):
                message = "ERROR bad imp_geoloc. In XML file: " + str(imp_geoloc_xml) + " WHEREAS previously set at: " + str(imp_geoloc)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('IMP_GEOLOC = ' + str(cfg.get('CONFIG_PARAMS', 'IMP_GEOLOC')))
            # Maximum tolerated difference (m) between improved height and initial height
            max_diff_height_geoloc_xml = float(in_xml_tree.xpath("//swot_product/processing_parameters/max_diff_height_geoloc")[0].text)
            max_diff_height_geoloc = cfg.getfloat(section, "MAX_DIFF_HEIGHT_GEOLOC")
            if (max_diff_height_geoloc_xml != max_diff_height_geoloc):
                message = "ERROR bad max_diff_height_geoloc In XML file: " + str(max_diff_height_geoloc_xml) + " WHEREAS previously set at: " + str(max_diff_height_geoloc)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('MAX_DIFF_HEIGHT_GEOLOC = ' + str(cfg.get('CONFIG_PARAMS', 'MAX_DIFF_HEIGHT_GEOLOC')))
            # Method to compute lake boundary or polygon hull
            # 0 = convex hull 
            # 1.0 = concave hull computed in ground geometry, based on Delaunay triangulation - using CGAL library
            # 1.1 = concave hull computed in ground geometry, based on Delaunay triangulation (time consuming)
            # 1.2 = concave hull computed in ground geometry, based on Delaunay triangulation with alpha parameter varying along range
            # 2.0 = edge computed in radar geometry, then converted in ground geometry by reordering lon lat (default)
            # 2.1 = edge computed in radar geometry, then converted in ground geometry by checking each segment
            hull_method_xml = float(in_xml_tree.xpath("//swot_product/processing_parameters/hull_method")[0].text)
            hull_method = cfg.getfloat(section, "HULL_METHOD")
            if (hull_method_xml != hull_method):
                message = "ERROR bad hull_method. In XML file: " + str(hull_method_xml) + " WHEREAS previously set at: " + str(hull_method)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('HULL_METHOD = ' + str(cfg.get('CONFIG_PARAMS', 'HULL_METHOD')))
            
            # Max number of pixels for hull computation method 1
            nb_pix_max_delauney_xml = int(in_xml_tree.xpath("//swot_product/processing_parameters/nb_pix_max_delauney")[0].text)
            nb_pix_max_delauney = cfg.getint(section, "NB_PIX_MAX_DELAUNEY")
            if (nb_pix_max_delauney_xml != nb_pix_max_delauney):
                message = "ERROR bad nb_pix_max_delauney. In XML file: " + str(nb_pix_max_delauney_xml) + " WHEREAS previously set at: " \
                                                                         + str(nb_pix_max_delauney)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('NB_PIX_MAX_DELAUNEY = ' + str(cfg.get('CONFIG_PARAMS', 'NB_PIX_MAX_DELAUNEY')))
    
            # Big lakes parameters for improved geoloc
            # Model to deal with big lake processing
            biglake_model_xml = in_xml_tree.xpath("//swot_product/processing_parameters/biglake_model")[0].text
            biglake_model = cfg.get(section, "BIGLAKE_MODEL")
            if (biglake_model_xml != biglake_model):
                message = "ERROR bad biglake_model. In XML file: " + str(biglake_model_xml) + " WHEREAS previously set at: " + str(biglake_model)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('BIGLAKE_MODEL = ' + str(cfg.get('CONFIG_PARAMS', 'BIGLAKE_MODEL')))
            # Min size for lake to be considered as big
            biglake_min_size_xml = float(in_xml_tree.xpath("//swot_product/processing_parameters/biglake_min_size")[0].text)
            biglake_min_size = cfg.getfloat(section, "BIGLAKE_MIN_SIZE")
            if (biglake_min_size_xml != biglake_min_size):
                message = "ERROR bad biglake_min_size. In XML file: " + str(biglake_min_size_xml) + " WHEREAS previously set at: " \
                          + str(biglake_min_size)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('BIGLAKE_MIN_SIZE = ' + str(cfg.get('CONFIG_PARAMS', 'BIGLAKE_MIN_SIZE')))
            # Grid spacing for lake height smoothing
            biglake_grid_spacing_xml = int(in_xml_tree.xpath("//swot_product/processing_parameters/biglake_grid_spacing")[0].text)
            biglake_grid_spacing = cfg.getint(section, "BIGLAKE_GRID_SPACING")
            if (biglake_grid_spacing_xml != biglake_grid_spacing):
                message = "ERROR bad biglake_grid_spacing. In XML file: " + str(biglake_grid_spacing_xml) + " WHEREAS previously set at: " \
                          + str(biglake_grid_spacing)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('BIGLAKE_GRID_SPACING = ' + str(cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_SPACING')))
            # Grid resolution for lake height smoothing
            biglake_grid_res_xml = int(in_xml_tree.xpath("//swot_product/processing_parameters/biglake_grid_res")[0].text)
            biglake_grid_res = cfg.getint(section, "BIGLAKE_GRID_RES")
            if (biglake_grid_res_xml != biglake_grid_res):
                message = "ERROR bad biglake_grid_res. In XML file: " + str(biglake_grid_res_xml) + " WHEREAS previously set at: " \
                          + str(biglake_grid_res)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('BIGLAKE_GRID_RES = ' + str(cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_RES')))
            
            # Input data for storage change computation
            stocc_input_xml = in_xml_tree.xpath("//swot_product/processing_parameters/stocc_input")[0].text
            stocc_input = cfg.get(section, "STOCC_INPUT")
            if (stocc_input_xml != stocc_input):
                message = "ERROR bad stocc_input. In XML file: " + str(stocc_input_xml) + " WHEREAS previously set at: " \
                          + str(stocc_input)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('STOCC_INPUT = ' + str(cfg.get('CONFIG_PARAMS', 'STOCC_INPUT')))
            
            # Input data for ice_dyn_f computation
            window_4_ice_clim_xml = in_xml_tree.xpath("//swot_product/processing_parameters/window_4_ice_clim")[0].text
            window_4_ice_clim = cfg.get(section, "WINDOW_4_ICE_CLIM")
            if (window_4_ice_clim_xml != window_4_ice_clim):
                message = "ERROR bad window_4_ice_clim. In XML file: " + str(stocc_input_xml) + " WHEREAS previously set at: " \
                          + str(window_4_ice_clim)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('WINDOW_4_ICE_CLIM = ' + str(cfg.get('CONFIG_PARAMS', 'WINDOW_4_ICE_CLIM')))
            
            # To add not-observed PLD lakes as empty features to the LakeSP_Prior shapefile (=True) or not (=False)
            add_all_xml = in_xml_tree.xpath("//swot_product/processing_parameters/add_all")[0].text
            add_all = cfg.getboolean(section, "ADD_ALL")
            if (add_all_xml != str(add_all)):
                message = "ERROR bad add_all. In XML file: " + str(add_all_xml) + " WHEREAS previously set at: " + str(add_all)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('ADD_ALL = ' + str(cfg.get('CONFIG_PARAMS', 'ADD_ALL')))
            # Nb digits for counter of lakes, used in the obs_id identifier of each observed lake
            nb_digits_xml = int(in_xml_tree.xpath("//swot_product/processing_parameters/nb_digits")[0].text)
            nb_digits = cfg.getint(section, "NB_DIGITS")
            if (nb_digits_xml != nb_digits):
                message = "ERROR bad nb_digits. In XML file: " + str(nb_digits_xml) + " WHEREAS previously set at: " \
                          + str(nb_digits)
                raise service_error.ProcessingError(message, logger)
            #logger.debug('NB_DIGITS = ' + str(cfg.get('CONFIG_PARAMS', 'NB_DIGITS')))

        except(IndexError):
            raise


#######################################


def attributes_fusion(attribute1, attribute2, fusion_mode, ratio):
    """
    Merge attributes depending on fusion mode :
    - ratio: weight average
    - max: attribute1 if ratio > 1, attribute2 otherwise
    - sum: attribute1 + attribute2
    - or: attribute1 | attribute2

    :param attribute1: value to merge
    :type attribute1: depends on the value
    :param attribute2: value to merge
    :type attribute2: depends on the value
    :param fusion_mode: ratio, max, sum, or
    :type fusion_mode: string
    :param ratio: weight to apply
    :type ratio: float

    :return: out_value = merged value
    :rtype: out_value = depends on the wanted type
    """
    logger = logging.getLogger("locnes_products_shapefile")
    
    out_value = None
    
    if ratio is None:
        out_value = None
    elif ratio == 1:
        out_value = attribute1
    elif ratio == 0:
        out_value = attribute2
    elif (attribute1 is None) or (attribute2 is None):
        out_value = None
    elif fusion_mode == 'ratio':
        out_value = ratio * attribute1 + (1 - ratio) * attribute2
    elif fusion_mode == 'max':
        if ratio >= 1:
            out_value = attribute1
        else:
            out_value = attribute2
    elif fusion_mode == 'sum':
        out_value = attribute1 + attribute2
    elif fusion_mode == 'or':
        out_value = attribute1 | attribute2
    else:
        message = "ERROR : unknown fusion mode %s with attributes %s %s" % (fusion_mode, str(attribute1), str(attribute2))
        raise service_error.ProcessingError(message, logger)
        
    return out_value
