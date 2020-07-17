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
.. module:: locnes_products_netcdf.py
    :synopsis: Deals with SWOT NetCDF products
     Created on 02/18/2019

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2019 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
"""
from collections import OrderedDict
import logging
from lxml import etree as ET
import numpy as np
import os

import cnes.common.lib.my_netcdf_file as my_nc
import cnes.common.lib_lake.locnes_product_general as my_prod


class NetcdfProduct(my_prod.LocnesProduct):
    """
    Deal with some of SWOT NetCDF products: LakeTile_pixcvec, LakeTile_edge and PIXCVec
    """
    
    def __init__(self, in_xml_file, in_inprod_metadata=None, in_proc_metadata=None):
        """
        Constructor: set the general values
        
        :param in_xml_file: XML file to populate the variables of the object
        :type in_xml_file: string
        :param in_inprod_metadata: metadata specific to input data product
        :type in_inprod_metadata: dict
        :param in_proc_metadata: metadata specific to current processing
        :type in_proc_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # 1 - Inheritance of LocnesProduct
        super().__init__(in_xml_file,
                         in_inprod_metadata=in_inprod_metadata,
                         in_proc_metadata=in_proc_metadata)
        
        # 2 - Update global metadata attributes
        self.global_metadata["Conventions"]["value"] = "CF-1.7"
        
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
        
        # 2 - Scan the available dimensions if there are some
        for element_niv1 in root:
            if element_niv1.tag == "shape":
                cur_shape_name = element_niv1.get("name")
                if cur_shape_name.endswith("_shape"):
                    self.list_shapes_names.append(cur_shape_name)
                    for element_niv2 in element_niv1:
                        dim_name = element_niv2.get("name")
                        if dim_name not in self.dims.keys():
                            self.dims[dim_name] = int(element_niv2.get("extent"))
            elif element_niv1.tag == "science":
                self.file_type = element_niv1.get("uid")
                        
        # 3 - Scan variables and metadata
        level_nodes = root[0][0]
        for element in level_nodes:
            cur_item = element.get("name")
            
            # 3.1 - The item is a global metadata element
            if cur_item.startswith("/@"):
                
                # Retrieve the exact name
                tmp_metadata = cur_item[2:]
                
                # Test existence (should not already be in the list)
                if tmp_metadata in self.global_metadata:
                    logger.warning("Metadata %s is repeated in the XML %s" % (tmp_metadata, in_xml_file))
                    
                else:  # Store information in a dictionary
                    
                    # Init
                    self.global_metadata[tmp_metadata] = OrderedDict()
                    # Type attributes
                    self.global_metadata[tmp_metadata]["type"] = str(element.tag)
                    self.global_metadata[tmp_metadata]["shape"] = element.get("shape")
                    self.global_metadata[tmp_metadata]["width"] = int(element.get("width"))
                    if "signed" in element.keys():
                        self.global_metadata[tmp_metadata]["signed"] = convert_2_boolean(element.get("signed"))
                    # Description
                    annotation = element[0]
                    self.global_metadata[tmp_metadata]["description"] = annotation.get("description")
                    # Init value
                    if tmp_metadata in ["title", "platform", "reference_document"]:
                        self.global_metadata[tmp_metadata]["value"] = self.global_metadata[tmp_metadata]["description"]
                    else:
                        if self.global_metadata[tmp_metadata]["type"] == "string":
                            self.global_metadata[tmp_metadata]["value"] = ""
                        else:
                            self.global_metadata[tmp_metadata]["value"] = -9999
            
            # 3.2 - The item is a variable
            else:
                
                # On détermine le nom de la variable
                cur_variable = cur_item[1:]
                
                # On initialise le dictionnaire associé s'il n'existe pas
                if cur_variable in self.attribute_metadata.keys():
                    logger.warning("Variable %s is repeated in the XML %s" % (cur_variable, in_xml_file))
                else:
                    self.attribute_metadata[cur_variable] = OrderedDict()
                    
                # Type
                self.attribute_metadata[cur_variable]["type"] = str(element.tag)
                self.attribute_metadata[cur_variable]["shape"] = element.get("shape")
                self.attribute_metadata[cur_variable]["width"] = int(element.get("width"))
                if "signed" in element.keys():
                    self.attribute_metadata[cur_variable]["signed"] = convert_2_boolean(element.get("signed"))
                    
                # Compute associated dtype
                tmp_type = self.attribute_metadata[cur_variable]["type"]
                tmp_width = self.attribute_metadata[cur_variable]["width"]
                if "signed" in element.keys():
                    tmp_signed = self.attribute_metadata[cur_variable]["signed"]
                else:
                    tmp_signed = False
                tmp_dtype = convert_type_xml_to_dtype(tmp_type, tmp_width, signed=tmp_signed)
                # Save the value
                if tmp_dtype is None:
                    logger.error("Triplet (type, width, signed)=(%s, %s, %s) is UNKNOWN" % (tmp_type, tmp_width, tmp_signed), exc_info=True)
                    raise
                self.attribute_metadata[cur_variable]["dtype"] = tmp_dtype
                
                # Other
                annotation = element[0]
                for key, value in annotation.items():
                    if key == "units":
                        try:
                            self.attribute_metadata[cur_variable][key] = np.int(value)
                        except:
                            self.attribute_metadata[cur_variable][key] = value
                    elif key == "flag_values":
                        try:
                            self.attribute_metadata[cur_variable][key] = np.array(value.split(" ")).astype(self.attribute_metadata[cur_variable]["dtype"])
                        except:
                            self.attribute_metadata[cur_variable][key] = value
                    else:
                        try:
                            if "signed" in self.attribute_metadata[cur_variable].keys():
                                self.attribute_metadata[cur_variable][key] = convert_str_to_dtype(value, 
                                                                                                   self.attribute_metadata[cur_variable]["type"],
                                                                                                   self.attribute_metadata[cur_variable]["width"],
                                                                                                   signed=self.attribute_metadata[cur_variable]["signed"])
                            else:
                                self.attribute_metadata[cur_variable][key] = convert_str_to_dtype(value, 
                                                                                                   self.attribute_metadata[cur_variable]["type"],
                                                                                                   self.attribute_metadata[cur_variable]["width"])
                        except:
                            self.attribute_metadata[cur_variable][key] = value
                            
    #----------------------------------------
        
    def write_product(self, in_out_file, in_size, in_variables):
        """
        Write NetCDF product
        
        :param in_out_file: output file full path
        :type in_out_file: string
        :param in_size: number of values in 1-D arrays stored in in_variables
        :type in_size: int
        :param in_variables: dictionary with key=variable name and value=variable value
        :type in_param: OrderedDict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # Open file in writing mode
        nc_writer = my_nc.MyNcWriter(in_out_file)
        
        # 1 - Write dimension(s)
        self.dims["points"] = in_size
        for key, value in self.dims.items():
            nc_writer.add_dimension(key, value)
        
        # 2 - Write variables
        if in_size == 0:
            logger.info("Empty NetCDF file generated")
        else:
            self.write_variables(nc_writer, in_variables)
        
        # 3 - Write global attributes
        self.write_metadata(nc_writer)
        
        # Close file
        nc_writer.close()
        
    def write_variables(self, in_nc_writer, in_attributes_value):
        """
        Write NetCDF variables listed in input
        
        :param in_nc_writer: NetCDF file writer
        :type in_nc_writer: my_netcdf_file.MyNcWriter
        :param in_attributes_value: dictionary with key=attribute name and value=attribute value
        :type in_attributes_value: OrderedDict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        var_keys = self.attribute_metadata.keys()
        
        for key, value in in_attributes_value.items():
            if key in var_keys:
                # Compute tuple about dimension(s)
                tmp_tuple_dims = ()
                for cur_dim in self.dims.keys():
                    if cur_dim in self.attribute_metadata[key]["shape"]:
                        tmp_tuple_dims += (cur_dim, )
                # Add variable
                tmp_var_metadata = dict()
                for tmp_key, tmp_value in self.attribute_metadata[key].items():
                    if tmp_key not in ["dtype", "shape", "type", "width", "signed", "value"]:
                        tmp_var_metadata[tmp_key] = tmp_value                
                in_nc_writer.add_variable(key, self.attribute_metadata[key]["dtype"], tmp_tuple_dims, in_attributes=tmp_var_metadata)
                # Fill variable
                in_nc_writer.fill_variable(key, value)
            else:
                logger.debug("Variable %s key is not known in the product" % key)
    
    def write_metadata(self, in_nc_writer):
        """
        Write global attributes
        
        :param in_nc_writer: NetCDF file writer
        :type in_nc_writer: my_netcdf_file.MyNcWriter
        """          
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        for name, att_dict in self.global_metadata.items():
            if "signed" in att_dict.keys():
                converted_value = convert_str_to_dtype(att_dict["value"], 
                                                       att_dict["type"],
                                                       att_dict["width"], 
                                                       signed=att_dict["signed"])
            else:
                converted_value = convert_str_to_dtype(att_dict["value"],
                                                       att_dict["type"],
                                                       att_dict["width"])
            in_nc_writer.add_global_attribute(name, converted_value)


#######################################


class LakeTilePixcvecProduct(NetcdfProduct):
    """
    Deal with LakeTile_PIXCVec NetCDF file
    """
    
    def __init__(self, in_pixcvecriver_metadata=None, in_proc_metadata=None):
        """
        Constructor; specific to LakeTile_PIXCVec NetCDF file
        
        :param in_pixcvecriver_metadata: metadata specific to L2_HR_PIXCVecRiver product
        :type in_pixcvecriver_metadata: dict
        :param in_proc_metadata: metadata specific to LakeTile processing
        :type in_proc_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Inheritance of NetcdfProduct
        super().__init__(os.path.join(os.path.dirname( __file__ ), "xml/laketile_pixcvec.xml"),
                          in_inprod_metadata=in_pixcvecriver_metadata,
                          in_proc_metadata=in_proc_metadata)


class LakeTileEdgeProduct(NetcdfProduct):
    """
    Deal with LakeTile_Edge NetCDF file
    """
    
    def __init__(self, in_pixc_metadata=None, in_proc_metadata=None):
        """
        Constructor; specific to LakeTile_Edge NetCDF file
        
        :param in_pixc_metadata: metadata specific to L2_HR_PIXC product
        :type in_pixc_metadata: dict
        :param in_proc_metadata: metadata specific to LakeTile processing
        :type in_proc_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Inheritance of NetcdfProduct
        super().__init__(os.path.join(os.path.dirname( __file__ ), "xml/laketile_edge.xml"),
                          in_inprod_metadata=in_pixc_metadata,
                          in_proc_metadata=in_proc_metadata)


class PixcvecProduct(NetcdfProduct):
    """
    Deal with PIXCVec NetCDF file
    """
    
    def __init__(self, in_laketile_pixcvec_metadata=None, in_proc_metadata=None):
        """
        Constructor; specific to PIXCVec NetCDF file
        
        :param in_laketile_pixcvec_metadata: metadata specific to L2_HR_LakeTile_Pixcvec product
        :type in_laketile_pixcvec_metadata: dict
        :param in_proc_metadata: metadata specific to LakeTile processing
        :type in_proc_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Inheritance of NetcdfProduct
        super().__init__(os.path.join(os.path.dirname( __file__ ), "xml/pixcvec.xml"),
                          in_inprod_metadata=in_laketile_pixcvec_metadata,
                          in_proc_metadata=in_proc_metadata)
            

#######################################
                

def convert_str_to_dtype(in_str_number, in_type_name, in_width, signed=True):
    """
    Convert in_str_number to approriate Numpy.dtype given type name, width and signed/unsigned information
    
    :param in_str_number: number in string to convert in the appropriate format
    :type in_str_number: str
    :param in_type_name: name of the type, among "integer", "real", "string"
    :type in_type_name: string
    :param in_width: width of type
    :type in_width: string
    :param signed: =True if the type is signed, =False otherwise
    :type signed: boolean
    
    :return: out_number = converted number
    :rtype: out_number = Numpy.dtype
    """
    logger = logging.getLogger("convert_str_to_dtype")
    
    # 0 - Init the output value
    if (in_str_number == "") and (in_type_name != "string"):
        tmp_number = -9999
    else:
        tmp_number = in_str_number
        
    # 1 - Compute the conversion to the wanted type, given its name, width and signed/unsigned flag
    if in_type_name == "integer":
        
        if signed:
            if in_width == 8:
                #out_number = np.int8(tmp_number)
                out_number = np.byte(tmp_number)
            elif in_width == 16:
                #out_number = np.int16(tmp_number)
                out_number = np.short(tmp_number)
            elif in_width == 64:
                out_number = np.int64(tmp_number)
            else:
                out_number = np.int32(tmp_number)
                
        else:
            if in_width == 8:
                #out_number = np.uint8(tmp_number)
                out_number = np.ubyte(tmp_number)
            elif in_width == 16:
                #out_number = np.uint16(tmp_number)
                out_number = np.ushort(tmp_number)
            elif in_width == 64:
                out_number = np.uint64(tmp_number)
            else:
                out_number = np.uint32(tmp_number)
                
    elif in_type_name == "real":
                
        if in_width == 64:
            out_number = np.float64(tmp_number)
        else:
            out_number = np.float32(tmp_number)
        
    else:
        
        if in_type_name != "string":
            logger.warning("Type %s UNKNOWN - Keep it as a string" % in_type_name)
        out_number = tmp_number
                
    return out_number
           

def convert_type_xml_to_dtype(in_type_name, in_width, signed=True):
    """
    Compute Numpy.dtype given type name, width and signed/unsigned information
    
    :param in_type_name: name of the type, among "integer", "real", "string"
    :type in_type_name: string
    :param in_width: width of type
    :type in_width: string
    :param signed: =True if the type is signed, =False otherwise
    :type signed: boolean
    
    :return: out_dtype = type in Numpy format given the input
    :rtype: out_dtype = Numpy.dtype
    """
    logger = logging.getLogger("convert_type_xml_to_type_netcdf")
    
    # 0 - Init the output value
    out_dtype = None
    
    # 1 - Compute the wanted Numpy.dtype, given the type name, width and signed/unsigned flag
    if in_type_name == "integer":
        
        if signed:
            if in_width == 8:
                #out_dtype = np.int8
                out_dtype = np.byte
            elif in_width == 16:
                #out_dtype = np.int16
                out_dtype = np.short
            elif in_width == 64:
                out_dtype = np.int64
            else:
                out_dtype = np.int32
                
        else:
            if in_width == 8:
                #out_dtype = np.uint8
                out_dtype = np.ubyte
            elif in_width == 16:
                #out_dtype = np.uint16
                out_dtype = np.ushort
            elif in_width == 64:
                out_dtype = np.uint64
            else:
                out_dtype = np.uint32
                
    elif in_type_name == "real":
        if in_width == 64:
            out_dtype = np.float64
        else:
            out_dtype = np.float32
            
    elif in_type_name == "string":
        out_dtype = "c"
        
    else:
        logger.error("Type %s UNKNOWN" % in_type_name)      
                
    return out_dtype
    


#######################################
    

def convert_2_boolean(in_str):
    """
    Convert string "true" or "false" into boolean True or False
    
    :param in_str: string ton convert into boolean
    :type in_str: string
    
    :return: out_boolean = converted string
    :rtype: out_boolean = boolean
    """
    
    out_boolean = False
    if in_str == "true":
        out_boolean = True
        
    return out_boolean
