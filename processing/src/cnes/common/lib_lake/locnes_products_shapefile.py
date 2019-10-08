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
.. module:: locnes_products_shapefile.py
    :synopsis: Deals with SWOT shapefile products
    Created on 02/26/2019

.. moduleauthor: Claire POTTIER - CNES DSO/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""

from collections import OrderedDict
import datetime
import logging
import os
from lxml import etree as ET
from osgeo import ogr, osr

import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_variables as my_var

import cnes.common.service_error as service_error
import cnes.common.service_config_file as service_config_file


def check_lake_db(in_xml_tree):
    """
    Check lake_db file
    
    :param IN_xml_tree: XML tree
    :type IN_xml_tree: etree.parse
    """ 
    # Get instance of service config file
    cfg = service_config_file.get_instance()
    # Get logger
    logger = logging.getLogger("locnes_products_shapefile")
    logger.info("Read xml lake_tile")
    
    # check Lake database
    try:
        lake_db_xml = in_xml_tree.xpath("//swot_product/processing_parameters/lake_db")[0].text
    except(IndexError):
        lake_db_xml = "None"
        logger.error("no lake_db information in lake_tile xml file")
    
    lake_db = cfg.get('DATABASES', 'LAKE_DB')
    if (lake_db is None) and (lake_db_xml == "None"):
        logger.warning('LAKE_DB not filled => LakeTile product not linked to a lake database')
    else:
        if (lake_db is None):
            message = "ERROR bad lake_db. In xml file: " + lake_db_xml + " in config file: None"
            raise service_error.ProcessingError(message, logger)
        elif (os.path.basename(lake_db) != lake_db_xml):
            message = "ERROR bad lake_db. In xml file: " + lake_db_xml + " in config file: " + os.path.basename(lake_db)
            raise service_error.ProcessingError(message, logger)
    logger.debug('LAKE_DB = ' + str(cfg.get('DATABASES', 'LAKE_DB')))


def check_lake_db_id(in_xml_tree):
    """
    Check lake_db file
    
    :param IN_xml_tree: XML tree
    :type IN_xml_tree: etree.parse
    """ 
    # Get instance of service config file
    cfg = service_config_file.get_instance()
    # Get logger
    logger = logging.getLogger("locnes_products_shapefile")
    logger.info("Read xml lake_tile")

    # check lake_db_id
    # Lake identifier attribute name in the database
    try:
        lake_db_id_xml = in_xml_tree.xpath("//swot_product/processing_parameters/lake_db_id")[0].text
    except(IndexError):
        lake_db_id_xml = "None"
        logger.error("no lake_db_id information in lake_tile xml file")
    
    lake_db_id = cfg.get('DATABASES', 'LAKE_DB_ID')
    if (lake_db_id is None) or (lake_db_id_xml == "None"):
        logger.warning('LAKE_DB_ID not filled')
    else:
        if (lake_db_id != lake_db_id_xml):
            message = "ERROR bad lake_db_id. In xml file: " + lake_db_id + " in config file: " + lake_db_id_xml
            raise service_error.ProcessingError(message, logger)
    logger.debug('LAKE_DB_ID = ' + str(cfg.get('DATABASES', 'LAKE_DB_ID')))


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
    logger.info("Read xml lake_tile")
    
    # Check if lake db specify in the command file is the same
    # check Lake database
    check_lake_db(in_xml_tree)
    # check lake_db_id
    check_lake_db_id(in_xml_tree)

    # add CONFIG_PARAMS section
    section = "CONFIG_PARAMS"
    if section not in cfg.sections():
        # add value if not in service_config_file
        cfg.add_section(section)
        # Water flag = 3=water near land edge  4=interior water
        flag_water = in_xml_tree.xpath("//swot_product/processing_parameters/flag_water")[0].text
        cfg.set(section, "FLAG_WATER", flag_water)
        logger.debug('FLAG_WATER = ' + str(cfg.get('CONFIG_PARAMS', 'FLAG_WATER')))
        # Dark water flag = 23=darkwater near land  24=interior dark water
        flag_dark = in_xml_tree.xpath("//swot_product/processing_parameters/flag_dark")[0].text
        cfg.set(section, "FLAG_DARK", flag_dark)
        logger.debug('FLAG_DARK = ' + str(cfg.get('CONFIG_PARAMS', 'FLAG_DARK')))
        # Min size for a lake to generate a lake product (=polygon + attributes) for it
        min_size = float(in_xml_tree.xpath("//swot_product/processing_parameters/min_size")[0].text)
        cfg.set(section, "MIN_SIZE", min_size)
        logger.debug('MIN_SIZE = ' + str(cfg.get('CONFIG_PARAMS', 'MIN_SIZE')))
        # Maximal standard deviation of height inside a lake (-1 = do not compute lake height segmentation)
        std_height_max = float(in_xml_tree.xpath("//swot_product/processing_parameters/std_height_max")[0].text)
        cfg.set(section, "STD_HEIGHT_MAX", std_height_max)
        logger.debug('STD_HEIGHT_MAX = ' + str(cfg.get('CONFIG_PARAMS', 'STD_HEIGHT_MAX')))
            
        # To improve PixC golocation (=True) or not (=False)
        imp_geoloc = bool(in_xml_tree.xpath("//swot_product/processing_parameters/imp_geoloc")[0].text)
        cfg.set(section, "IMP_GEOLOC", imp_geoloc)
        logger.debug('IMP_GEOLOC = ' + str(cfg.get('CONFIG_PARAMS', 'IMP_GEOLOC')))
        # Method to compute lake boundary or polygon hull
        # 0 = convex hull
        # 1.0 = concave hull computed in ground geometry, based on Delaunay triangulation - using CGAL library
        # 1.1 = concave hull computed in ground geometry, based on Delaunay triangulation - with alpha parameter varying across-track
        # 2 = edge computed in radar geometry, then converted in ground geometry (default)
        hull_method = float(in_xml_tree.xpath("//swot_product/processing_parameters/hull_method")[0].text)
        cfg.set(section, "HULL_METHOD", hull_method)
        logger.debug('HULL_METHOD = ' + str(cfg.get('CONFIG_PARAMS', 'HULL_METHOD')))
        # max number of pixel for hull computation 1
        nb_pix_max_delauney = int(in_xml_tree.xpath("//swot_product/processing_parameters/nb_pix_max_delauney")[0].text)
        cfg.set(section, "NB_PIX_MAX_DELAUNEY", nb_pix_max_delauney)
        logger.debug('NB_PIX_MAX_DELAUNEY = ' + str(cfg.get('CONFIG_PARAMS', 'NB_PIX_MAX_DELAUNEY')))
        # max number of contour points for hull computation 2
        nb_pix_max_contour = int(in_xml_tree.xpath("//swot_product/processing_parameters/nb_pix_max_contour")[0].text)
        cfg.set(section, "NB_PIX_MAX_CONTOUR", nb_pix_max_contour)
        logger.debug('NB_PIX_MAX_CONTOUR = ' + str(cfg.get('CONFIG_PARAMS', 'NB_PIX_MAX_CONTOUR')))
        
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

    else:
        # check values
        try:
            # Water flag = 3=water near land edge  4=interior water
            flag_water_xml = in_xml_tree.xpath("//swot_product/processing_parameters/flag_water")[0].text
            flag_water = cfg.get(section, "FLAG_WATER")
            if (flag_water_xml != flag_water):
                message = "ERROR bad flag_water. In xml file: " + str(flag_water_xml) + " in config file: " + str(flag_water)
                raise service_error.ProcessingError(message, logger)
            logger.debug('FLAG_WATER = ' + str(cfg.get('CONFIG_PARAMS', 'FLAG_WATER')))
            # Dark water flag = 23=darkwater near land  24=interior dark water
            flag_dark_xml = in_xml_tree.xpath("//swot_product/processing_parameters/flag_dark")[0].text
            flag_dark = cfg.get(section, "FLAG_DARK")
            if (flag_dark_xml != flag_dark):
                message = "ERROR bad flag_dark. In xml file: " + str(flag_dark_xml) + " in config file: " + str(flag_dark)
                raise service_error.ProcessingError(message, logger)
            logger.debug('FLAG_DARK = ' + str(cfg.get('CONFIG_PARAMS', 'FLAG_DARK')))
            # Min size for a lake to generate a lake product (=polygon + attributes) for it
            min_size_xml = float(in_xml_tree.xpath("//swot_product/processing_parameters/min_size")[0].text)
            min_size = cfg.getfloat(section, "MIN_SIZE")
            if (min_size_xml != min_size):
                message = "ERROR bad min_size. In xml file: " + str(min_size_xml) + " in config file: " + str(min_size)
                raise service_error.ProcessingError(message, logger)
            logger.debug('MIN_SIZE = ' + str(cfg.get('CONFIG_PARAMS', 'MIN_SIZE')))
            # Maximal standard deviation of height inside a lake (-1 = do not compute lake height segmentation)
            std_height_max_xml = float(in_xml_tree.xpath("//swot_product/processing_parameters/std_height_max")[0].text)
            std_height_max = cfg.getfloat(section, "STD_HEIGHT_MAX")
            if (std_height_max_xml != std_height_max):
                message = "ERROR bad std_height_max. In xml file: " + str(std_height_max_xml) + " in config file: " + str(std_height_max)
                raise service_error.ProcessingError(message, logger)
            logger.debug('STD_HEIGHT_MAX = ' + str(cfg.get('CONFIG_PARAMS', 'STD_HEIGHT_MAX')))
    
            # To improve PixC golocation (=True) or not (=False)
            imp_geoloc_xml = bool(in_xml_tree.xpath("//swot_product/processing_parameters/imp_geoloc")[0].text)
            imp_geoloc = cfg.getboolean(section, "IMP_GEOLOC")
            if (imp_geoloc_xml != imp_geoloc):
                message = "ERROR bad imp_geoloc. In xml file: " + str(imp_geoloc_xml) + " in config file: " + str(imp_geoloc)
                raise service_error.ProcessingError(message, logger)
            logger.debug('IMP_GEOLOC = ' + str(cfg.get('CONFIG_PARAMS', 'IMP_GEOLOC')))
            # Method to compute lake boundary or polygon hull
            # 0 = convex hull
            # 1.0 = concave hull computed in ground geometry, based on Delaunay triangulation - using CGAL library
            # 1.1 = concave hull computed in ground geometry, based on Delaunay triangulation - with alpha parameter varying across-track
            # 2 = edge computed in radar geometry, then converted in ground geometry (default)
            hull_method_xml = float(in_xml_tree.xpath("//swot_product/processing_parameters/hull_method")[0].text)
            hull_method = cfg.getfloat(section, "HULL_METHOD")
            if (hull_method_xml != hull_method):
                message = "ERROR bad hull_method. In xml file: " + str(hull_method_xml) + " in config file: " + str(hull_method)
                raise service_error.ProcessingError(message, logger)
            logger.debug('HULL_METHOD = ' + str(cfg.get('CONFIG_PARAMS', 'HULL_METHOD')))
            # max number of pixel for hull computation 1
            nb_pix_max_delauney_xml = int(in_xml_tree.xpath("//swot_product/processing_parameters/nb_pix_max_delauney")[0].text)
            nb_pix_max_delauney = cfg.getint(section, "NB_PIX_MAX_DELAUNEY")
            if (nb_pix_max_delauney_xml != nb_pix_max_delauney):
                message = "ERROR bad nb_pix_max_delauney. In xml file: " + str(nb_pix_max_delauney_xml) + " in config file: " \
                                                                         + str(nb_pix_max_delauney)
                raise service_error.ProcessingError(message, logger)
            logger.debug('NB_PIX_MAX_DELAUNEY = ' + str(cfg.get('CONFIG_PARAMS', 'NB_PIX_MAX_DELAUNEY')))
            # max number of contour points for hull computation 2
            nb_pix_max_contour_xml = int(in_xml_tree.xpath("//swot_product/processing_parameters/nb_pix_max_contour")[0].text)
            nb_pix_max_contour = cfg.getint(section, "NB_PIX_MAX_CONTOUR")
            if (nb_pix_max_contour_xml != nb_pix_max_contour):
                message = "ERROR bad nb_pix_max_contour. In xml file: " + str(nb_pix_max_contour_xml) + " in config file: " \
                                                                        + str(nb_pix_max_contour)
                raise service_error.ProcessingError(message, logger)
            logger.debug('NB_PIX_MAX_CONTOUR = ' + str(cfg.get('CONFIG_PARAMS', 'NB_PIX_MAX_CONTOUR')))
    
            # Big lakes parameters for improved geoloc
            # Model to deal with big lake processing
            biglake_model_xml = in_xml_tree.xpath("//swot_product/processing_parameters/biglake_model")[0].text
            biglake_model = cfg.get(section, "BIGLAKE_MODEL")
            if (biglake_model_xml != biglake_model):
                message = "ERROR bad biglake_model. In xml file: " + str(biglake_model_xml) + " in config file: " + str(biglake_model)
                raise service_error.ProcessingError(message, logger)
            logger.debug('BIGLAKE_MODEL = ' + str(cfg.get('CONFIG_PARAMS', 'BIGLAKE_MODEL')))
            # Min size for lake to be considered as big
            biglake_min_size_xml = float(in_xml_tree.xpath("//swot_product/processing_parameters/biglake_min_size")[0].text)
            biglake_min_size = cfg.getfloat(section, "BIGLAKE_MIN_SIZE")
            if (biglake_min_size_xml != biglake_min_size):
                message = "ERROR bad biglake_min_size. In xml file: " + str(biglake_min_size_xml) + " in config file: " + str(biglake_min_size)
                raise service_error.ProcessingError(message, logger)
            logger.debug('BIGLAKE_MIN_SIZE = ' + str(cfg.get('CONFIG_PARAMS', 'BIGLAKE_MIN_SIZE')))
            # Grid spacing for lake height smoothing
            biglake_grid_spacing_xml = int(in_xml_tree.xpath("//swot_product/processing_parameters/biglake_grid_spacing")[0].text)
            biglake_grid_spacing = cfg.getint(section, "BIGLAKE_GRID_SPACING")
            if (biglake_grid_spacing_xml != biglake_grid_spacing):
                message = "ERROR bad biglake_grid_spacing. In xml file: " + str(biglake_grid_spacing_xml) + " in config file: " \
                          + str(biglake_grid_spacing)
                raise service_error.ProcessingError(message, logger)
            logger.debug('BIGLAKE_GRID_SPACING = ' + str(cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_SPACING')))
            # Grid resolution for lake height smoothing
            biglake_grid_res_xml = int(in_xml_tree.xpath("//swot_product/processing_parameters/biglake_grid_res")[0].text)
            biglake_grid_res = cfg.getint(section, "BIGLAKE_GRID_RES")
            if (biglake_grid_res_xml != biglake_grid_res):
                message = "ERROR bad biglake_grid_res. In xml file: " + str(biglake_grid_res_xml) + " in config file: " + str(biglake_grid_res)
                raise service_error.ProcessingError(message, logger)
            logger.debug('BIGLAKE_GRID_RES = ' + str(cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_RES')))

        except(IndexError):
            raise


#######################################

        
def set_dico_val(in_out_dico, in_values):
    """
    Setter of dictionary values
        
    :param in_out_dico: original dictionary in which the values will be modified
    :type in_out_dico: dict
    :param in_values: values with which we fill the original dictionary
    :type in_values: dict
    """
    
    # Recopy the dictionary in input
    dico_keys = in_out_dico.keys()
    
    # Update values in input dictionary
    for key, value in in_values.items():
        if key in dico_keys:
            in_out_dico[key] = value
        
   
def convert_wrt_type(in_val, in_type):
    """
    Convert in_val in the format specified in input
    
    :param in_val: value to convert
    :type in_val: depends on the value
    :param in_type: type of the wanted output (default=str)
    :type in_type: type OGR
    """
    
    if in_type == ogr.OFTInteger:
        retour = int(in_val)
    elif in_type == ogr.OFTReal:
        retour = float(in_val)
    else:
        retour = str(in_val)
    return retour
            
    
#######################################


class ShapefileProduct(object):
    """
    Deal with SWOT Shapefile products: LakeTile_shp
    """
    
    def __init__(self):
        """
        Constructor: set the general values
        
        Variables of the object:
        - attributes / OrderedDict: dictionary having key=attributes names and value=OrderedDict of shapefile attributes
        - metadata / OrderedDict: dictionary having key=XML attributes and value=their values
        - data_source / OGRDataSource: data source of the shape
        - layer / OGRLayer: data layer
        - layer_defn / OGRLayerDefn: layer definition
        """
        
        # 1 - Init attributes
        self.attributes = OrderedDict()
        
        # 2 - Init metadata
        self.metadata = OrderedDict()
        # 2.1 - General metadata
        general_metadata = OrderedDict()
        general_metadata["conventions"] = "Esri conventions as given in \"ESRI Shapefile Technical Description, an ESRI White Paper, July 1998\" http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf"
        general_metadata["title"] = ""
        general_metadata["institution"] = "CNES"
        general_metadata["source"] = "Ka-band radar interferometer"
        general_metadata["history"] = "%sZ: Creation" % my_tools.swot_timeformat(datetime.datetime.utcnow(), in_format=1)
        general_metadata["mission_name"] = "SWOT"
        general_metadata["references"] = ""
        general_metadata["reference_document"] = ""
        general_metadata["contact"] = "claire.pottier@cnes.fr"
        self.metadata['global_attributes'] = general_metadata
        # 2.2 - Metadata specific to granule specification
        self.metadata['content'] = OrderedDict()
        self.metadata['content'] = self.attributes
        
        # 3 - Shapefile attributes
        self.data_source = None  # Data source
        self.layer = None  # Data layer
        self.layer_defn = None  # Layer definition
        
    def free(self):
        """
        Destroy memory layer
        """
        if self.data_source is not None:
            self.data_source.Destroy()

    #----------------------------------------
    
    def create_mem_layer(self, in_layer_name, in_geom_type, in_spatial_ref=4326):
        """
        Create product memory layer and attributes creation
        
        :param in_layer_name: name for product layer        
        :type in_layer_name: string
        :param in_geom_type: type of geometry
        :type in_geom_type: OGRGeometry
        :param in_spatial_ref: name of the spatial reference (default=4326=WGS84)
        :type in_spatial_ref: int        
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        mem_driver = ogr.GetDriverByName(str('MEMORY'))  # Driver for memory layers

        # 1 - Create memory layer
        logger.info('> Creating memory layer')
        self.data_source = mem_driver.CreateDataSource('memData')
        
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
        
        for key, value in self.attributes.items():
            ogr_type = my_var.FORMAT_OGR[value['type']]
            if ogr_type == 0:
                logger.debug("> Create attribute %s of type ogr.OFTInteger" % key)
            elif ogr_type == 2:
                logger.debug("> Create attribute %s of type ogr.OFTReal" % key)
            elif ogr_type == 4:
                logger.debug("> Create attribute %s of type ogr.OFTString" % key)
            else:
                logger.debug("> Create attribute {} of type {}".format(key, ogr_type))
            self.layer.CreateField(ogr.FieldDefn(str(key), my_var.FORMAT_OGR[value['type']]))

    #----------------------------------------
    
    def add_feature(self, in_geom, in_attributes):
        """
        Add feature with geometry given by in_geom and attribute values given by in_attributes
        
        :param in_geom: geometry to add to feature
        :type in_geom: OGRGeometry
        :param in_attributes: list of attributes and their value to
        :type in_attributes: dict
        """
        
        # 1 - Create the object
        feature = ogr.Feature(self.layer_defn)
        
        # 2 - Add geometry
        feature.SetGeometry(in_geom)
        
        # 3 - Add attribute values
        att_keys = in_attributes.keys()  # List of attributes to value
        for att_name, attr_carac in self.attributes.items():  # Loop over the whole list of the existing attributes
            value_to_write = None  # Init value to write
            if att_name in att_keys:  # If in the list of attributes to modify
                if in_attributes[att_name] is None:  
                    value_to_write = my_var.FV_SHP[attr_carac["type"]]  # Fill to _FillValue if None
                else:
                    value_to_write = convert_wrt_type(in_attributes[att_name], attr_carac["type"])  # Fill with appropriate type
            else:
                value_to_write = my_var.FV_SHP[attr_carac["type"]]  # Fill to _FillValue if not in list
            feature.SetField(str(att_name), value_to_write)  # Fill the field with the wanted value
            
        # 4 - Add feature
        self.layer.CreateFeature(feature)
        
        # 5 - Destroy the feature to free resources
        feature.Destroy()

    #----------------------------------------
    
    def write_metadata_file(self, in_filename):
        """
        Write the metadata file associated to the lake shapefile
        
        :param in_filename: full path for the metadata file
        :type in_filename: string
        """
        
        # 1 - Define root element
        root = ET.Element("swot_product")
        
        # 2 - Creating the tree
        for key1, value1 in self.metadata.items():  # 1st sub-level
            sub1 = ET.SubElement(root, str(key1))
            for key2, value2 in value1.items():  # 2nd sub-level
                if str(value2).startswith("{"):  # Go to 3rd level
                    sub2 = ET.SubElement(sub1, str(key2))
                    for key3, value3 in value2.items():
                        ET.SubElement(sub2, str(key3)).text = str(value3)
                else:
                    ET.SubElement(sub1, str(key2)).text = str(value2)

        # 3 - Write tree element in XML file
        tree = ET.ElementTree(root)
        tree.write(in_filename, pretty_print=True, xml_declaration=True, encoding='utf-8')


#######################################


class LakeSPShpProduct(ShapefileProduct):
    """
    Deal with LakeTile_shp shapefile and LakeSP product
    """
    
    def __init__(self, in_product_type, in_layer_name, in_pixc_metadata=None, in_laketile_metadata=None, in_proc_metadata=None):
        """
        Constructor
        
        :param in_product_type: type of product among "SP"=LakeSP and "TILE"=LakeTile
        :type in_product_type: string
        :param in_layer_name: name for product layer        
        :type in_layer_name: string
        :param in_pixc_metadata: metadata specific to L2_HR_PIXC product
        :type in_pixc_metadata: dict
        :param in_laketile_metadata: metadata specific to L2_HR_LakeTile_shp shapefile
        :type in_laketile_metadata: dict
        :param in_proc_metadata: metadata specific to LakeTile or LakeSP processing
        :type in_proc_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Init NetCDF_product
        super().__init__()
        
        # 2 - Init attributes info 
        
        self.attributes['obs_id'] = {'type': "text", 
                       'long_name': "identifier of the observed lake", 
                       'tag_basic_expert': "Basic", 
                       'comment': "Unique lake identifier within the product. The format of the identifier is CBBTTTSNNNNNN, where C=continent, B=basin, TTT=tile number within the pass, S=swath side, N=lake counter within the tile."}

        self.attributes['lake_id'] = {'type': "text", 
                       'fill_value': "no_data", 
                       'long_name': "lake ID(s) from prior database", 
                       'tag_basic_expert': "Basic", 
                       'comment': "List of identifiers of prior lakes that intersect the observed lake. The format of the identifier is CBBNNNNNNT, where C=continent, B=basin, N=lake counter within the basin, T=type. The different lake identifiers are separated by semicolons."}

        self.attributes['time'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "time (UTC)", 
                       'standard_name': "time", 
                       'calendar': "gregorian", 
                       'tai_utc_difference': "[value of TAI-UTC at time of first record]", 
                       'leap_second': "YYYY-MM-DD hh:mm:ss", 
                       'units': "seconds since 2000-01-01 00:00:00.000", 
                       'tag_basic_expert': "Basic", 
                       'comment': "Time of measurement in seconds in the UTC time scale since 1 Jan 2000 00:00:00 UTC. [tai_utc_difference] is the difference between TAI and UTC reference time (seconds) for the first measurement of the data set. If a leap second occurs within the data set, the attribute leap_second is set to the UTC time at which the leap second occurs."}

        self.attributes['time_tai'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "time (TAI)", 
                       'standard_name': "time", 
                       'calendar': "gregorian", 
                       'units': "seconds since 2000-01-01 00:00:00.000", 
                       'tag_basic_expert': "Basic", 
                       'comment': "Time of measurement in seconds in the TAI time scale since 1 Jan 2000 00:00:00 TAI. This time scale contains no leap seconds. The difference (in seconds) with time in UTC is given by the attribute [time:tai_utc_difference]."}

        self.attributes['time_str'] = {'type': "text", 
                       'fill_value': "no_data", 
                       'long_name': "UTC time", 
                       'standard_name': "time", 
                       'calendar': "gregorian", 
                       'tai_utc_difference': "[value of TAI-UTC at time of first record]", 
                       'leap_second': "YYYY-MM-DD hh:mm:ss", 
                       'tag_basic_expert': "Basic", 
                       'comment': "Time string giving UTC time. The format is YYYY-MM-DDThh:mm:ss.ssssssZ, where the Z suffix indicates UTC time."}

        self.attributes['wse'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "lake-averaged water surface elevation with respect to the geoid", 
                       'units': "m", 
                       'valid_min': -1000, 
                       'valid_max': 100000, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Lake-averaged water surface elevation, relative to the provided model of the geoid (geoid_hght), with corrections for media delays (wet and dry troposphere, and ionosphere), crossover correction, and tidal effects (solid_tide, load_tide1, and pole_tide) applied."}

        self.attributes['wse_u'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "total uncertainty in lake water surface elevation", 
                       'units': "m", 
                       'valid_min': 0, 
                       'valid_max': 100, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Total one-sigma uncertainty (random and systematic) in the lake WSE, including uncertainties of corrections and references."}

        self.attributes['wse_r_u'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "random-only uncertainty in the height water surface elevation", 
                       'units': "m", 
                       'valid_min': 0, 
                       'valid_max': 100, 
                       'tag_basic_expert': "Expert", 
                       'comment': "Random-only component in the lake water surface elevation, including uncertainties of corrections and references, and variation about the fit."}

        self.attributes['wse_std'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "standard deviation of pixels wse", 
                       'units': "m", 
                       'valid_min': -1000, 
                       'valid_max': 100000, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Standard deviation of the water surface elevation of all the pixels composing the lake. Note that this value is therefore with respect to the provided model of the geoid (geoid_hght attribute) whereas the height of pixels is given with respect to the ellipsoid. This parameter is computed only for large lakes (> 5000ha)."}

        self.attributes['area_total'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "total water area with estimate of dark water", 
                       'units': "m^2", 
                       'valid_min': 0, 
                       'valid_max': 2000000, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Total estimated area, including dark water that was not detected as water in the SWOT observation but identified through the use of a prior water likelihood map."}

        self.attributes['area_tot_u'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "uncertainty in total water area", 
                       'units': "m^2", 
                       'valid_min': 0, 
                       'valid_max': 2000000, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Total uncertainty (random and systematic) in the total water area."}

        self.attributes['area_detct'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "area of detected water pixels", 
                       'units': "m^2", 
                       'valid_min': 0, 
                       'valid_max': 1000000000, 
                       'tag_basic_expert': "Expert", 
                       'comment': "Aggregation of used detected pixels area."}

        self.attributes['area_det_u'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "uncertainty in area of detected water", 
                       'units': "m^2", 
                       'valid_min': 0, 
                       'valid_max': 1000000000, 
                       'tag_basic_expert': "Expert", 
                       'comment': "Total uncertainty (random and systematic) in the area of detected water pixels."}

        self.attributes['layovr_val'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "metric of layover effect", 
                       'units': "TBD", 
                       'valid_min': "TBD", 
                       'valid_max': "TBD", 
                       'tag_basic_expert': "Expert", 
                       'comment': "Value indicating an estimate of the height error due to layover."}

        self.attributes['xtrk_dist'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "distance of lake polygon centroid to the satellite ground track", 
                       'units': "m", 
                       'valid_min': -75000, 
                       'valid_max': 75000, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Distance of centroid of polygon delineating lake boundary to the satellite ground track. A negative value indicates the left side of the swath, relative to the spacecraft velocity vector. A positive value indicates the right side of the swath."}

        self.attributes['delta_s_l'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "storage change computed by linear method", 
                       'units': "m^3", 
                       'valid_min': -10000000, 
                       'valid_max': 10000000, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Storage change with regards to the reference area and height from PLD"}

        self.attributes['ds_l_u'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "uncertainty in storage change computed by linear method", 
                       'units': "m^3", 
                       'valid_min': -10000000, 
                       'valid_max': 10000000, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Uncertainty in storage change computed by linear method."}

        self.attributes['delta_s_q'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "storage change computed by quadratic method", 
                       'units': "m^3", 
                       'valid_min': -10000000, 
                       'valid_max': 10000000, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Storage change with regards to the reference area and height from PLD"}

        self.attributes['ds_q_u'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "uncertainty in storage change computed by quadratic method", 
                       'units': "m^3", 
                       'valid_min': -10000000, 
                       'valid_max': 10000000, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Uncertainty in storage change computed by quadratic method."}

        self.attributes['quality_f'] = {'type': "int4", 
                       'fill_value': -999, 
                       'long_name': "summary quality indicator for lake measurement", 
                       'flag_meanings': "good bad", 
                       'flag_values': "0 1", 
                       'valid_min': 0, 
                       'valid_max': 1, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Summary quality flag for the lake measurement. Values of 0 and 1 indicate nominal and off-nominal measurements."}

        self.attributes['dark_frac'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "fractional area of dark water", 
                       'units': 1, 
                       'valid_min': 0, 
                       'valid_max': 1, 
                       'tag_basic_expert': "Expert", 
                       'comment': "Fraction of lake area_total covered by dark water. The value is between 0 and 1."}

        self.attributes['ice_clim_f'] = {'type': "int4", 
                       'fill_value': -999, 
                       'long_name': "climatological ice cover flag", 
                       'source': "UNC", 
                       'flag_meanings': "no_ice_cover partial_ice_cover full_ice_cover not_available", 
                       'flag_values': "0 1 2 255", 
                       'valid_min': 0, 
                       'valid_max': 255, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Climatological ice cover flag indicating whether the lake is ice-covered on the day of the observation based on external climatological information (not the SWOT measurement). Values of 0, 1, and 2 indicate that the lake is not ice covered, partially ice covered, and fully ice covered, respectively. A value of 255 indicates that this flag is not available."}

        self.attributes['ice_dyn_f'] = {'type': "int4", 
                       'fill_value': -999, 
                       'long_name': "dynamical ice cover flag", 
                       'source': "UNC", 
                       'flag_meanings': "no_ice_cover partial_ice_cover full_ice_cover not_available", 
                       'flag_values': "0 1 2 255", 
                       'valid_min': 0, 
                       'valid_max': 255, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Dynamic ice cover flag indicating whether the lake is ice-covered on the day of the observation based on analysis of external satellite optical data. Values of 0, 1, and 2 indicate that the lake is not ice covered, partially ice covered, and fully ice covered, respectively. A value of 255 indicates that this flag is not available."}

        self.attributes['partial_f'] = {'type': "int4", 
                       'fill_value': -999, 
                       'long_name': "partially covered lake flag", 
                       'flag_meanings': "covered partially_covered", 
                       'flag_values': "0 1", 
                       'valid_min': 0, 
                       'valid_max': 1, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Flag that indicates only partial lake coverage.  0= Indicates that the observed lake is entirely covered by the swath. 1= Indicates that the observed lake is partially covered by the swath."}

        self.attributes['xovr_cal_q'] = {'type': "int4", 
                       'fill_value': -999, 
                       'long_name': "quality of the cross-over calibrations", 
                       'flag_meanings': "TBD", 
                       'flag_masks': "TBD", 
                       'flag_values': "TBD", 
                       'valid_min': 0, 
                       'valid_max': "TBD", 
                       'tag_basic_expert': "Basic", 
                       'comment': "Quality of the cross-over calibration."}

        self.attributes['geoid_hght'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "geoid height", 
                       'standard_name': "geoid_height_above_reference_ellipsoid", 
                       'source': "EGM2008", 
                       'institution': "GSFC", 
                       'units': "m", 
                       'valid_min': -150, 
                       'valid_max': 150, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Lake-averaged geoid model height above the reference ellipsoid. The value is computed from the EGM2008 geoid model with a correction to refer the value to the mean tide system (i.e., includes the zero-frequency permanent tide)."}

        self.attributes['solid_tide'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "solid Earth tide height", 
                       'source': "Cartwright and Edden [1973] Corrected tables of tidal harmonics - J. Geophys. J. R. Astr. Soc., 33, 253-264", 
                       'units': "m", 
                       'valid_min': -1, 
                       'valid_max': 1, 
                       'tag_basic_expert': "Expert", 
                       'comment': "Solid-Earth (Body) tide height, averaged over the lake. The zero-frequency permanent tide component is not included. The value is computed from the Cartwright/Taylor model."}

        self.attributes['pole_tide'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "height of pole tide", 
                       'units': "m", 
                       'valid_min': 0, 
                       'valid_max': 0, 
                       'tag_basic_expert': "Expert", 
                       'comment': "Geocentric pole tide height. The sum total of the contribution from the solid-Earth (body) pole tide height and the load pole tide height (i.e., the effect of the ocean pole tide loading of the Earth's crust)."}

        self.attributes['load_tide1'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "geocentric load tide height", 
                       'source': "FES2014", 
                       'institution': "LEGOS/CNES", 
                       'units': "m", 
                       'valid_min': 0, 
                       'valid_max': 0, 
                       'tag_basic_expert': "Expert", 
                       'comment': "Geocentric load tide height. The effect of the ocean tide loading of the Earth's crust. This value is used to compute wse. The value is computed from the FES2014 ocean tide model."}

        self.attributes['load_tide2'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "geocentric load tide height", 
                       'source': "GOT4.10c", 
                       'institution': "GSFC", 
                       'units': "m", 
                       'valid_min': 0, 
                       'valid_max': 0, 
                       'tag_basic_expert': "Expert", 
                       'comment': "Geocentric load tide height. The effect of the ocean tide loading of the Earth's crust. The value is computed from the GOT4.10c ocean tide model."}

        self.attributes['dry_trop_c'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "dry tropospheric vertical correction to WSE", 
                       'source': "European Centre for Medium-Range Weather Forecasting", 
                       'units': "m", 
                       'valid_min': -3, 
                       'valid_max': -1, 
                       'tag_basic_expert': "Expert", 
                       'comment': "Equivalent vertical correction due to dry troposphere delay. Adding the reported correction to the reported lake WSE results in the uncorrected lake WSE. The value is computed from the European Centre for Medium-Range Weather Forecasts (ECMWF) model."}

        self.attributes['wet_trop_c'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "wet tropospheric vertical correction to WSE", 
                       'source': "European Centre for Medium-Range Weather Forecasting", 
                       'units': "m", 
                       'valid_min': -1, 
                       'valid_max': 0, 
                       'tag_basic_expert': "Expert", 
                       'comment': "Equivalent vertical correction due to wet troposphere delay. Adding the reported correction to the reported lake WSE results in the uncorrected lake WSE. The value is computed from the ECMWF model."}

        self.attributes['iono_c'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "ionospheric vertical correction to WSE", 
                       'source': "Global Ionosphere Maps", 
                       'institution': "JPL", 
                       'units': "m", 
                       'valid_min': 0, 
                       'valid_max': 0, 
                       'tag_basic_expert': "Expert", 
                       'comment': "Equivalent vertical correction due to ionosphere delay. Adding the reported correction to the reported lake WSE results in the uncorrected lake WSE. The value is computed from the JPL Global Ionosphere Maps (GIM)."}

        self.attributes['xovr_cal_c'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "crossover calibration height correction", 
                       'units': "m", 
                       'valid_min': -10, 
                       'valid_max': 10, 
                       'tag_basic_expert': "Expert", 
                       'comment': "Equivalent height correction estimated from KaRIn crossover calibration. The correction is applied during processing before geolocation in terms of roll, baseline dilation, etc., but reported as an equivalent height correction. The correction term should be subtracted from the reported WSE to obtain the uncorrected WSE."}

        self.attributes['p_name'] = {'type': "text", 
                       'fill_value': "no_data", 
                       'long_name': "name(s) of the lake", 
                       'comment': "Name(s) of the lake, retrieved from Open Street Map, IGN Carthage, GLWD and vMap0 databases. The different names are separated by semicolons."}

        self.attributes['grand_id'] = {'type': "int9", 
                       'fill_value': -99999999, 
                       'long_name': "reservoir Id from GRanD database", 
                       'source': "https://doi.org/10.1890/100125", 
                       'valid_min': 0, 
                       'valid_max': 9999, 
                       'tag_basic_expert': "Expert", 
                       'comment': "Reservoir ID from the Global Reservoir and Dam (GRanD) database. 0=The lake is not a registered reservoir."}

        self.attributes['p_height'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "reference height", 
                       'units': "m", 
                       'valid_min': -1000, 
                       'valid_max': 100000, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Reference height from the PLD, used to compute the storage change."}

        self.attributes['p_area'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "reference area", 
                       'units': "m", 
                       'valid_min': 0, 
                       'valid_max': 500000000000, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Reference area from the PLD, used to compute the storage change."}

        self.attributes['p_storage'] = {'type': "float", 
                       'fill_value': -999999999999, 
                       'long_name': "maximum water storage", 
                       'units': "m^3", 
                       'valid_min': 0, 
                       'valid_max': 10000000, 
                       'tag_basic_expert': "Basic", 
                       'comment': "Maximum water storage value from the PLD, computed between the minimum (or ground when a bathymetry is available) and maximum observed levels of the lake."}
        
        # 3 - Init metadata depending on product type
        if in_product_type == "TILE":
            self.init_metadata_laketile(in_pixc_metadata=in_pixc_metadata, in_proc_metadata=in_proc_metadata)
        elif in_product_type == "SP":
            self.init_metadata_lakesp(in_laketile_metadata=in_laketile_metadata, in_proc_metadata=in_proc_metadata)
        else:
            message = "ERROR = product type is %s ; should be SP or TILE" % in_product_type
            raise service_error.ProcessingError(message, logger)
            
        # 4 - Create memory layer
        self.create_mem_layer(in_layer_name, ogr.wkbMultiPolygon)
        
    def init_metadata_laketile(self, in_pixc_metadata=None, in_proc_metadata=None):
        """
        Init metadata specific to LakeTile_shp shapefile
        
        :param in_pixc_metadata: metadata specific to L2_HR_PIXC product
        :type in_pixc_metadata: dict
        :param in_proc_metadata: metadata specific to LakeTile processing
        :type in_proc_metadata: dict
        """
        # Get instance of service config file
        cfg = service_config_file.get_instance()
        
        # 1 - Update general metadata
        self.metadata["global_attributes"]["title"] = "Level 2 KaRIn high rate lake tile vector product – LakeTile_shp"
        self.metadata["global_attributes"]["references"] = ""
        self.metadata["global_attributes"]["reference_document"] = "SWOT-TN-CDM-0673-CNES"
        
        # 2 - Metadata retrieved from L2_HR_PIXC product
        self.metadata["global_attributes"]["cycle_number"] = -9999
        self.metadata["global_attributes"]["pass_number"] = -9999
        self.metadata["global_attributes"]["tile_number"] = -9999
        self.metadata["global_attributes"]["swath_side"] = ""
        self.metadata["global_attributes"]["tile_name"] = ""
        self.metadata["global_attributes"]["start_time"] = ""
        self.metadata["global_attributes"]["stop_time"] = ""
        self.metadata["global_attributes"]["inner_first_latitude"] = -9999.0
        self.metadata["global_attributes"]["inner_first_longitude"] = -9999.0
        self.metadata["global_attributes"]["inner_last_latitude"] = -9999.0
        self.metadata["global_attributes"]["inner_last_longitude"] = -9999.0
        self.metadata["global_attributes"]["outer_first_latitude"] = -9999.0
        self.metadata["global_attributes"]["outer_first_longitude"] = -9999.0
        self.metadata["global_attributes"]["outer_last_latitude"] = -9999.0
        self.metadata["global_attributes"]["outer_last_longitude"] = -9999.0
        self.metadata["global_attributes"]["continent"] = "None"
        if in_pixc_metadata is not None:
            set_dico_val(self.metadata["global_attributes"], in_pixc_metadata)
            
        # 3 - Metadata specific to processing
        self.metadata["global_attributes"]["xref_input_l2_hr_pixc_file"] = ""
        self.metadata["global_attributes"]["xref_input_l2_hr_pixc_vec_river_file"] = ""
        self.metadata["global_attributes"]["xref_static_lake_db_file"] = ""
        self.metadata["global_attributes"]["xref_l2_hr_lake_tile_param_file"] = ""
        if in_proc_metadata is not None:
            set_dico_val(self.metadata["global_attributes"], in_proc_metadata)
        
        # 4 - Configuration parameters metadata
        self.metadata['processing_parameters'] = OrderedDict()
        #self.metadata["processing_parameters"]["lake_db"] = str(cfg.get('DATABASES', 'LAKE_DB'))
        self.metadata["processing_parameters"]["lake_db"] = os.path.basename(str(cfg.get('DATABASES', 'LAKE_DB')))
        self.metadata["processing_parameters"]["lake_db_id"] = str(cfg.get('DATABASES', 'LAKE_DB_ID'))
        self.metadata["processing_parameters"]["flag_water"] = str(cfg.get('CONFIG_PARAMS', 'FLAG_WATER'))
        self.metadata["processing_parameters"]["flag_dark"] = str(cfg.get('CONFIG_PARAMS', 'FLAG_DARK'))
        self.metadata["processing_parameters"]["min_size"] = cfg.get('CONFIG_PARAMS', 'MIN_SIZE')
        self.metadata["processing_parameters"]["std_height_max"] = cfg.get('CONFIG_PARAMS', 'STD_HEIGHT_MAX')
        self.metadata["processing_parameters"]["imp_geoloc"] = cfg.get('CONFIG_PARAMS', 'IMP_GEOLOC')
        self.metadata["processing_parameters"]["hull_method"] = cfg.get('CONFIG_PARAMS', 'HULL_METHOD')
        self.metadata["processing_parameters"]["nb_pix_max_delauney"] = cfg.get('CONFIG_PARAMS', 'NB_PIX_MAX_DELAUNEY')
        self.metadata["processing_parameters"]["nb_pix_max_contour"] = cfg.get('CONFIG_PARAMS', 'NB_PIX_MAX_CONTOUR')
        self.metadata["processing_parameters"]["biglake_model"] = cfg.get('CONFIG_PARAMS', 'BIGLAKE_MODEL')
        self.metadata["processing_parameters"]["biglake_min_size"] = cfg.get('CONFIG_PARAMS', 'BIGLAKE_MIN_SIZE')
        self.metadata["processing_parameters"]["biglake_grid_spacing"] = cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_SPACING')
        self.metadata["processing_parameters"]["biglake_grid_res"] = cfg.get('CONFIG_PARAMS', 'BIGLAKE_GRID_RES')
        
    def init_metadata_lakesp(self, in_laketile_metadata=None, in_proc_metadata=None):
        """
        Init metadata specific to LakeSP shapefile
        
        :param in_laketile_metadata: metadata specific to L2_HR_LakeTile_shp shapefile
        :type in_laketile_metadata: dict
        :param in_proc_metadata: metadata specific to LakeSP processing
        :type in_proc_metadata: dict
        """
        
        # 1 - Update general metadata
        self.metadata["global_attributes"]["title"] = "Level 2 KaRIn high rate lake single pass vector product"
        self.metadata["global_attributes"]["references"] = ""
        self.metadata["global_attributes"]["reference_document"] = "SWOT-TN-CDM-0673-CNES"
        
        # 2 - Metadata retrieved from L2_HR_PIXC product
        self.metadata["global_attributes"]["cycle_number"] = -9999
        self.metadata["global_attributes"]["pass_number"] = -9999
        self.metadata["global_attributes"]["start_time"] = ""
        self.metadata["global_attributes"]["stop_time"] = ""
        self.metadata["global_attributes"]["polygon"] = ""
        self.metadata["global_attributes"]["continent"] = "None"
        if in_laketile_metadata is not None:
            set_dico_val(self.metadata["global_attributes"], in_laketile_metadata)
            
        # 3 - Metadata specific to processing
        self.metadata["global_attributes"]["xref_input_l2_hr_pixc_file"] = ""
        self.metadata["global_attributes"]["xref_input_l2_hr_lake_tile_files"] = ""
        self.metadata["global_attributes"]["xref_static_lake_db_file"] = ""
        self.metadata["global_attributes"]["xref_l2_hr_lake_sp_param_file"] = ""
        if in_proc_metadata is not None:
            set_dico_val(self.metadata["global_attributes"], in_proc_metadata)
        
    #----------------------------------------
        
    def update_and_write_metadata(self, in_filename, in_pixc_metadata=None, in_proc_metadata=None):
        """
        Write .shp.xml file
        
        :param in_filename: full path for the metadata file
        :type in_filename: string
        :param in_pixc_metadata: metadata specific to L2_HR_PIXC product
        :type in_pixc_metadata: dict
        :param in_proc_metadata: metadata specific to processing
        :type in_proc_metadata: dict
        """
        
        # Update metadata
        if in_pixc_metadata is not None:
            set_dico_val(self.metadata["global_attributes"], in_pixc_metadata)
        if in_proc_metadata is not None:
            set_dico_val(self.metadata["global_attributes"], in_proc_metadata)

        self.write_metadata_file(in_filename)
        