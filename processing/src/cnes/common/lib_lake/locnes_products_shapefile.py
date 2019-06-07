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
from lxml import etree as ET
from osgeo import ogr, osr

import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_variables as my_var

import cnes.common.service_error as service_error

        
def set_dico_val(in_out_dico, in_values):
    """
    Setter of dictionary values
        
    :param in_out_dico: original dictionary in which the values will be modified
    :type in_out_dico: dict
    :param in_values: values with which filling the original dictionary
    :type in_values: dict
    """
    logger = logging.getLogger("locnes_products_shapefile")
        
    dico_keys = in_out_dico.keys()
    
    for key, value in in_values.items():
        if key in dico_keys:
            in_out_dico[key] = value
        else:
            logger.debug("%s key is not known" % key)
        
   
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


class Shapefile_product(object):
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
        general_metadata["conventions"] = "ESRI shapefile: http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf"
        general_metadata["title"] = ""
        general_metadata["institution"] = "CNES"
        general_metadata["source"] = ""
        general_metadata["history"] = "%sZ: Creation" % my_tools.swot_timeformat(datetime.datetime.utcnow(), in_format=1)
        general_metadata["mission_name"] = "SWOT"
        general_metadata["references"] = ""
        general_metadata["reference_document"] = ""
        general_metadata["contact"] = "claire.pottier@cnes.fr"
        self.metadata['general'] = general_metadata
        # 2.2 - Metadata specific to granule specification
        self.metadata['granule'] = OrderedDict()
        # 2.3 - Metadata specific to processing
        self.metadata['processing'] = OrderedDict()
        
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
            if value['type'] == ogr.OFTReal:
                logger.debug("> Create attribute %s of type ogr.OFTReal(width=%d, precision=%d)" % (key, value['width'], value['precision']))
                self.add_field_real(key, value['width'], value['precision'])
            else:
                if value['type'] == ogr.OFTInteger:
                    logger.debug("> Create attribute %s of type ogr.OFTInteger" % key)
                elif value['type'] == ogr.OFTString:
                    logger.debug("> Create attribute %s of type ogr.OFTString" % key)
                self.layer.CreateField(ogr.FieldDefn(str(key), value['type']))
        
    def add_field_real(self, in_name, in_width, in_precision):
        """
        Add a real field to layer
        
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
        self.layer.CreateField(tmp_field)  # Add field to the layer

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
                    value_to_write = my_var.FV_OGR[attr_carac["type"]]  # Fill to _FillValue if None
                else:
                    value_to_write = convert_wrt_type(in_attributes[att_name], attr_carac["type"])  # Fill with appropriate type
            else:
                value_to_write = my_var.FV_OGR[attr_carac["type"]]  # Fill to _FillValue if not in list
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
        for key1, value1 in self.metadata.items():
            sub1 = ET.SubElement(root, str(key1))  # 1st sub-level
            for key2, value2 in value1.items():
                ET.SubElement(sub1, str(key2)).text = str(value2)

        # 3 - Write tree element in XML file
        tree = ET.ElementTree(root)
        tree.write(in_filename, pretty_print=True, xml_declaration=True, encoding='utf-8')


#######################################


class LakeSPShp_product(Shapefile_product):
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
        
        # 2 - Init attributes info specific to LakeTile_shp file
        
        # 2.1 - Identifiers
        # Object identifier
        self.attributes["lakeobs_id"] = {'type': ogr.OFTString}
        # List of lakes in the a priori DB related to the object
        self.attributes["lakedb_id"] = {'type': ogr.OFTString}
        
        # 2.2 - Mean datetime of observation
        self.attributes["time"] = {'type': ogr.OFTInteger}  # UTC time 
        self.attributes["time_tai"] = {'type': ogr.OFTInteger}  # TAI time
        self.attributes["time_str"] = {'type': ogr.OFTString}  # Time in UTC as a string
        
        # 2.3 - Measured hydrology parameters
        # Mean water surface elevation of the lake and uncertainty
        self.attributes["wse"] = {'type': ogr.OFTReal, 'width': 14, 'precision': 2}
        self.attributes["wse_u"] = {'type': ogr.OFTReal, 'width': 14, 'precision': 2}
        # Water surface elevation standard deviation (only for big lakes)
        self.attributes["wse_std"] = {'type': ogr.OFTReal, 'width': 14, 'precision': 2}
        # Area of detected water pixels and uncertainty
        self.attributes["area_detct"] = {'type': ogr.OFTReal, 'width': 14, 'precision': 2}
        self.attributes["area_det_u"] = {'type': ogr.OFTReal, 'width': 14, 'precision': 2}
        # Total water area and uncertainty
        self.attributes["area_total"] = {'type': ogr.OFTReal, 'width': 14, 'precision': 2}
        self.attributes["area_tot_u"] = {'type': ogr.OFTReal, 'width': 14, 'precision': 2}
        # Area of pixels used to compute water surface elevation
        self.attributes["area_wse"] = {'type': ogr.OFTReal, 'width': 14, 'precision': 2}
        # Metric of layover effect
        self.attributes["layovr_val"] = {'type': ogr.OFTReal, 'width': 14, 'precision': 2}
        # Average distance from polygon centroid to the satellite ground track
        self.attributes["xtrk_dist"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 3}
        
        # 2.4 - Storage change and uncertainty
        # Linear model
        self.attributes["delta_s_L"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 3}
        self.attributes["ds_L_u"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 3}
        # Quadratic model
        self.attributes["delta_s_Q"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 3}
        self.attributes["ds_Q_u"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 3}
        
        # 2.5 - Quality indicators
        # Summary quality indicator
        self.attributes["quality_f"] = {'type': ogr.OFTInteger}
        # Fractional area of dark water
        self.attributes["dark_frac"] = {'type': ogr.OFTReal, 'width': 8, 'precision': 6}
        # Ice cover flags
        self.attributes["ice_clim_f"] = {'type': ogr.OFTInteger}
        self.attributes["ice_dyn_f"] = {'type': ogr.OFTInteger}
        # Partial flag: =1 if the lake is partially covered by the swath, 0 otherwise
        self.attributes["partial_f"] = {'type': ogr.OFTInteger}
        # Quality of cross-over calibrations
        self.attributes["xovr_cal_f"] = {'type': ogr.OFTInteger}
        
        # 2.6 - KaRIn sigma0 information
        self.attributes["sig0"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 6}
        self.attributes["sig0_u"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 6}
        
        # 2.7 - Geophysical references
        # Geoid model height
        self.attributes["geoid"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 3}
        # Solid earth tide
        self.attributes["earth_tide"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 3}
        # Pole tide
        self.attributes["pole_tide"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 3}
        # Load tide
        self.attributes["load_tide1"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 3}  # GOT4.10
        self.attributes["load_tide2"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 3}  # FES2014
        
        # 2.8 - Geophysical range corrections
        # Dry tropo corr
        self.attributes["dry_trop_c"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 3}
        # Wet tropo corr
        self.attributes["wet_trop_c"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 3}
        # Iono corr
        self.attributes["iono_c"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 3}
        
        # 2.9 - Instrument corrections
        # KaRIn correction from crossover cal processing evaluated for lake 
        self.attributes["xovr_cal_c"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 3}
        
        # 2.10 - Prior lake database information
        # Lake name
        self.attributes["p_name"] = {'type': ogr.OFTString}
        # Dam identifier from GRanD database
        self.attributes["grand_id"] = {'type': ogr.OFTString}
        # Reference water surface elevation, used to compute storage change
        self.attributes["p_wse"] = {'type': ogr.OFTReal, 'width': 13, 'precision': 3}
        # Reference area, used to compute storag change
        self.attributes["p_area"] = {'type': ogr.OFTReal, 'width': 14, 'precision': 2}
        
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
        
        # 1 - Update general metadata
        self.metadata["general"]["title"] = "Level 2 KaRIn high rate lake tile vector product – LakeTile_shp"
        self.metadata["general"]["references"] = ""
        self.metadata["general"]["reference_document"] = "SWOT-TN-CDM-0673-CNES"
        
        # 2 - Metadata retrieved from L2_HR_PIXC product
        self.metadata["granule"]["cycle_number"] = -9999
        self.metadata["granule"]["pass_number"] = -9999
        self.metadata["granule"]["tile_number"] = -9999
        self.metadata["granule"]["swath_side"] = ""
        self.metadata["granule"]["tile_name"] = ""
        self.metadata["granule"]["start_time"] = ""
        self.metadata["granule"]["stop_time"] = ""
        self.metadata["granule"]["inner_first_latitude"] = -9999.0
        self.metadata["granule"]["inner_first_longitude"] = -9999.0
        self.metadata["granule"]["inner_last_latitude"] = -9999.0
        self.metadata["granule"]["inner_last_longitude"] = -9999.0
        self.metadata["granule"]["outer_first_latitude"] = -9999.0
        self.metadata["granule"]["outer_first_longitude"] = -9999.0
        self.metadata["granule"]["outer_last_latitude"] = -9999.0
        self.metadata["granule"]["outer_last_longitude"] = -9999.0
        self.metadata["granule"]["continent"] = ""
        self.metadata["granule"]["ellipsoid_semi_major_axis"] = ""
        self.metadata["granule"]["ellipsoid_flattening"] = ""
        if in_pixc_metadata is not None:
            set_dico_val(self.metadata["granule"], in_pixc_metadata)
            
        # 3 - Processing metadata
        self.metadata["processing"]["xref_static_lake_db_file"] = ""
        self.metadata["processing"]["xref_input_l2_hr_pixc_file"] = ""
        self.metadata["processing"]["xref_input_l2_hr_pixc_vec_river_file"] = ""
        self.metadata["processing"]["xref_l2_hr_lake_tile_param_file"] = ""
        if in_proc_metadata is not None:
            set_dico_val(self.metadata["processing"], in_proc_metadata)
        
    def init_metadata_lakesp(self, in_laketile_metadata=None, in_proc_metadata=None):
        """
        Init metadata specific to LakeSP shapefile
        
        :param in_laketile_metadata: metadata specific to L2_HR_LakeTile_shp shapefile
        :type in_laketile_metadata: dict
        :param in_proc_metadata: metadata specific to LakeSP processing
        :type in_proc_metadata: dict
        """
        
        # 1 - Update general metadata
        self.metadata["general"]["title"] = "Level 2 KaRIn high rate lake single pass vector product"
        self.metadata["general"]["references"] = ""
        self.metadata["general"]["reference_document"] = "SWOT-TN-CDM-0673-CNES"
        
        # 2 - Metadata retrieved from L2_HR_PIXC product
        self.metadata["granule"]["cycle_number"] = -9999
        self.metadata["granule"]["pass_number"] = -9999
        self.metadata["granule"]["start_time"] = ""
        self.metadata["granule"]["stop_time"] = ""
        self.metadata["granule"]["polygon"] = ""
        self.metadata["granule"]["continent"] = ""
        self.metadata["granule"]["ellipsoid_semi_major_axis"] = ""
        self.metadata["granule"]["ellipsoid_flattening"] = ""
        if in_laketile_metadata is not None:
            set_dico_val(self.metadata["granule"], in_laketile_metadata)
            
        # 3 - Processing metadata
        self.metadata["processing"]["xref_static_lake_db_file"] = ""
        self.metadata["processing"]["xref_input_l2_hr_pixc_file"] = ""
        self.metadata["processing"]["xref_input_l2_hr_pixc_vec_river_file"] = ""
        self.metadata["processing"]["xref_l2_hr_lake_sp_param_file"] = ""
        if in_proc_metadata is not None:
            set_dico_val(self.metadata["processing"], in_proc_metadata)

    #----------------------------------------
        
    def update_and_write_metadata(self, in_filename, in_pixc_metadata=None, in_proc_metadata=None):
        """
        Write .shp.xml file
        
        :param in_filename: full path for the metadata file
        :type in_filename: string
        :param in_pixc_metadata: metadata specific to L2_HR_PIXC product
        :type in_pixc_metadata: dict
        :param in_proc_metadata: metadata specific to LakeTile processing
        :type in_proc_metadata: dict
        """

        # Update metadata
        if in_pixc_metadata is not None:
            set_dico_val(self.metadata["granule"], in_pixc_metadata)
        if in_proc_metadata is not None:
            set_dico_val(self.metadata["processing"], in_proc_metadata)

        self.write_metadata_file(in_filename)
        
