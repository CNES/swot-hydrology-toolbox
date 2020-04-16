# -*- coding: utf8 -*-
"""
.. module product_netcdf.py
    :synopsis: Deal with NetCDF products
    Created on 04/09/2020

.. module author: Claire POTTIER - CNES DSO/SI/TR

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""

from collections import OrderedDict
from datetime import datetime
from lxml import etree as ET
import numpy as np
import os
import lib.my_api as my_api

import cnes.common.lib.my_netcdf_file as my_nc


class GroupNetCDF(object):
    """
    Class dealing with a NetCDF group
    """
    
    def __init__(self):
        """
        Constructor: set the general values
        
        Variables of the object:
        - list_shapes / list: list of shapes used in the group
        - variables / OrderedDict: dictionary having key=variables names and value=OrderedDict of NetCDF variable attributes
        - metadata / OrderedDict: dictionary having key=global attributes and value=their values
        """
        my_api.printInfo("[GroupNetCDF] == INIT ==")
        
        # 1 - Init list of shapes
        self.list_shapes = []
        
        # 2 - Init variables
        self.variables = OrderedDict()
        
        # 3 - Init metadata
        self.metadata = OrderedDict()


class NetCDFProduct(object):
    """
    Main class dealing with NetCDF products
    """
    
    def __init__(self, in_xml_file=None):
        """
        Constructor: set the general values
        
        :param in_xml_file: XML file to populate the variables of the object
        :type in_xml_file: string
        
        Variables of the object:
        - list_shapes / dict of dict: dictionary having key=shape name and value=dictionary of dimensions related to the shape
        - group_names / list of string: list of name of each NetCDF group; it is empty of no group in the NetCDF
        - list_groups / dict of GroupNetCDF: dictionary having key=group name and value=GroupNetCDF object
        - metadata / OrderedDict: dictionary having key=global attributes and value=their values
        """
        my_api.printInfo("[NetCDFProduct] == INIT ==")
        
        # 1 - Init dict of shapes
        self.list_shapes = {}
        
        # 2 - Init list of group names
        self.group_names = []
        
        # 3 - Init dict of groups attributes
        self.list_groups = {}
        
        # 4 - Init metadata
        self.metadata = OrderedDict()
        
        # 5 - Populate the variables from the XML file
        if in_xml_file is not None:
            self.read_content_from_xml(in_xml_file)
        
    #----------------------------------------
    
    def read_content_from_xml(self, in_xml_file):
        """
        Populate the variables of the object from an XML file (generated from the Product Description Document)
        
        :param in_xml_file: XML file to populate the variables of the object
        :type in_xml_file: string
        """
        my_api.printInfo("[NetCDFProduct] == read_content_from_xml ==")
        
        # 1 - Load the XML file
        xml_reader = ET.parse(in_xml_file)
        root = xml_reader.getroot()
        
        # 2 - Scan the available dimensions
        self.list_shapes = {}
        for element_niv1 in root:
            if element_niv1.tag == "shape":
                cur_shape_name = element_niv1.get("name")
                if cur_shape_name.endswith("_shape"):
                    self.list_shapes[cur_shape_name] = {}
                    for element_niv2 in element_niv1:
                        self.list_shapes[cur_shape_name][element_niv2.get("name")] = int(element_niv2.get("extent"))
                        
        # 3 - Scan variables and metadata
        level_nodes = root[0][0]
        for element in level_nodes:
            cur_item = element.get("name")
            
            # 3.1 - The item is a global metadata element
            if cur_item.startswith("/@"):
                # Retrieve the exact name
                tmp_metadata = cur_item[2:]
                # Test existence (should not already be in the list)
                if tmp_metadata in self.metadata:
                    my_api.exitWithError("[ERROR] Metadata %s is repeated in the XML %s" % (tmp_metadata, in_xml_file))
                else:
                    # Add global metadata element to the list
                    self.metadata[tmp_metadata] = ""  
            
            # 3.2 - The item is a variable
            elif cur_item.startswith("/"):
                
                # Test if there are groups in the NetCDF
                if "/" in cur_item[1:]:
                    tmp_split_item = cur_item[1:].split("/")
                    # 3.2.1 - Group information
                    cur_group = tmp_split_item[0]  # Retrieve the exact name of the group
                    # 3.2.2 - Item information
                    cur_variable = tmp_split_item[1]  # Retrieve the exact name
                else:
                    # 3.2.1 - Group information
                    cur_group = "no_group"  # Use a fake group
                    # 3.2.2 - Item information
                    cur_variable = cur_item[1:]  # Retrieve the exact name
                    
                # If this is a new group, declare it
                if cur_group not in self.group_names:
                    self.group_names.append(cur_group)  # Add the name to the list
                    self.list_groups[cur_group] = GroupNetCDF()
                
                # The name corresponds to a metadata
                if cur_variable.startswith("@"):
                    cur_metadata = cur_variable[1:]
                    if cur_metadata in self.list_groups[cur_group].metadata:
                        my_api.exitWithError("[ERROR] Group metadata %s/%s is repeated in the XML %s" % (cur_group, cur_metadata, in_xml_file))
                    else:
                        self.list_groups[cur_group].metadata[cur_metadata] = ""
                        
                # The name corresponds to a variable in the group
                else:            
                    if cur_variable in self.list_groups[cur_group].variables.keys():
                        my_api.exitWithError("[ERROR] Group variable %s/%s is repeated in the XML %s" % (cur_group, cur_variable, in_xml_file))
                    else:
                        self.list_groups[cur_group].variables[cur_variable] = OrderedDict()
                    # Dimensions
                    var_shape = element.get("shape")
                    self.list_groups[cur_group].variables[cur_variable]["shape"] = var_shape
                    if var_shape not in self.list_groups[cur_group].list_shapes:
                        self.list_groups[cur_group].list_shapes.append(var_shape)
                    # Type
                    flag_signed = False
                    if element.get("signed"):
                        flag_signed = True
                    var_type = find_dtype(element.tag, int(element.get("width")), flag_signed)
                    self.list_groups[cur_group].variables[cur_variable]["dtype"] = var_type
                    # Other
                    annotation = element[0]
                    for key, value in annotation.items():
                        if key not in ["app", "_FillValue"]:
                            self.list_groups[cur_group].variables[cur_variable][key] = value
        
    #----------------------------------------
        
    def set_dim_val(self, in_dim_name, in_dim_value):
        """
        Setter of dimension value
        
        :param in_dim_name: name of the dimension to set
        :type in_dim_name: str
        :param in_dim_value: value to attribute to the dimension in_name
        :type in_dim_value: int
        """
        
        for shape_name, shape_value in self.list_shapes.items():
            if in_dim_name in shape_value.keys():
                self.list_shapes[shape_name][in_dim_name] = in_dim_value
        
    def set_metadata_val(self, in_metadata, group=None):
        """
        Setter of metadata value
        
        :param in_metadata: metadata stored as a dictionary with key=metadata name and value=metadata value
        :type in_metadata: dict
        :param group: name of the group (optionnal)
        :type group: string
        """
        
        # 1 - Set pointer on metadata
        if group is None:
            metadata_pointer = self.metadata
        else:
            metadata_pointer = self.list_groups[group].metadata
        
        # 2 - Update metadata
        for key, value in in_metadata.items():
            if key in metadata_pointer.keys():
                metadata_pointer[key] = value

    #----------------------------------------
        
    def write_product(self, in_out_file, in_vars_value=None):
        """
        Write NetCDF product
        
        :param in_out_file: output file full path
        :type in_out_file: string
        :param in_vars_value: variables value stored as a dictionary with key=variable name and value=variable value
        :type in_vars_value: dict
        """
        my_api.printInfo("[NetCDFProduct] == write_product : %s ==" % in_out_file)
        
        # Open file in writing mode
        nc_writer = my_nc.MyNcWriter(in_out_file)
        
        # 1 - Write global attributes
        self.write_metadata(nc_writer)
        
        # 2 - Write groups, variables and metadata
        # If there is no group in the NetCDF, there is one single group object named "no_group"
        for cur_group_name in self.group_names:
            
            # 2.1 - Add group
            if cur_group_name == "no_group":
                cur_group = None
                vars_value = in_vars_value
            else:
                cur_group = nc_writer.add_group(cur_group_name)
                if in_vars_value is None:
                    vars_value = None
                else:
                    vars_value = in_vars_value[cur_group_name]
            
            # 2.2 - Add group attributes
            self.write_metadata(nc_writer, group=cur_group) 
            
            # 2.3 - Add group dimensions
            self.write_dimensions(nc_writer, group=cur_group)
            
            # 2.4 - Add group variables
            self.write_variables(nc_writer, vars_value, group=cur_group)
        
        # 3 - Close file
        nc_writer.close()
        
    def write_dimensions(self, in_nc_writer, group=None):
        """
        Write dimensions, dedicated to a group of not
        
        :param in_nc_writer: NetCDF file writer
        :type in_nc_writer: my_netcdf_file.MyNcWriter
        :param group: NetCDF group in which the metadata will be added; if None: metadata are considered as global
        :type group: netCDF4._netCDF4.Group
        """
        
        # 0 - Init variables
        dims_to_write = {}
        tmp_list_shapes = None
        
        # 1 - Set the list of shapes to consider, given a group or not
        if group is None:
            tmp_list_shapes = self.list_shapes
        else:
            tmp_list_shapes = self.list_groups[group.name].list_shapes
            
        # 2 - Compute the list of dimensions to write
        for cur_shape in tmp_list_shapes:
            for key, value in self.list_shapes[cur_shape].items():
                if key not in dims_to_write.keys():
                    dims_to_write[key] = value
                    
        # 3 - Write dimensions
        for key, value in dims_to_write.items():
            in_nc_writer.add_dimension(key, value, in_group=group)
        
    def write_variables(self, in_nc_writer, in_vars_value=None, group=None):
        """
        Write NetCDF variables listed in input
        
        :param in_nc_writer: NetCDF file writer
        :type in_nc_writer: my_netcdf_file.MyNcWriter
        :param in_vars_value: variables value stored as a dictionary with key=variable name and value=variable value
        :type in_vars_value: dict
        :param group: NetCDF group in which the metadata will be added; if None: metadata are considered as global
        :type group: netCDF4._netCDF4.Group
        """
        
        # 1 - Set the list of variables to consider, given a group or not
        if group is None:
            tmp_list_vars = self.list_groups["no_group"].variables
        else:
            tmp_list_vars = self.list_groups[group.name].variables
        
        for var_name, var_att in tmp_list_vars.items():
            
            # 2 - Determine the dimension to use
            var_shape = var_att["shape"]
            var_dims = tuple(self.list_shapes[var_shape].keys())
            
            # 3 - Create the variable
            in_nc_writer.add_variable(var_name, tmp_list_vars[var_name]['dtype'], var_dims, in_group=group, in_attributes=tmp_list_vars[var_name])
            
            # 4 - Fill variable
            if in_vars_value is not None:
                if var_name in in_vars_value.keys():
                    in_nc_writer.fill_variable(var_name, in_vars_value[var_name], in_group=group)
                else:
                    my_api.printInfo("[product_netcdf/NetCDFProduct/write_variables] Variable %s not populated" % var_name)
    
    def write_metadata(self, in_nc_writer, group=None):
        """
        Write global attributes, dedicated to a group of not
        
        :param in_nc_writer: NetCDF file writer
        :type in_nc_writer: my_netcdf_file.MyNcWriter
        :param group: NetCDF group in which the metadata will be added; if None: metadata are considered as global
        :type group: netCDF4._netCDF4.Group
        """
        
        # 1 - Set the metadata to consider, given a group or not
        if group is None:
            metadata_pointer = self.metadata
        else:
            metadata_pointer = self.list_groups[group.name].metadata
        
        # 2 - Add the metadata
        for key, value in metadata_pointer.items():
            in_nc_writer.add_global_attribute(key, value, in_group=group)


#######################################


class PixcProduct(NetCDFProduct):
    """
    Deal with PIXC NetCDF file
    """
    
    def __init__(self):
        """
        Constructor; specific to PIXC NetCDF file
        """
        
        # 1 - Init CommonPixc
        super().__init__(os.path.expandvars("$SWOT_HYDROLOGY_TOOLBOX/sisimp/data/pixc.xml"))
                
        # 2 - Update general metadata specific to PIXC file
        tmp_metadata = {}
        tmp_metadata['conventions'] = 'CF-1.7'
        tmp_metadata['title'] = 'Level 2 KaRIn High Rate Water Mask Pixel Clould Data Product'
        tmp_metadata['institution'] = 'CNES - Large scale simulator'
        tmp_metadata['source'] = 'Ka-band radar interferometer'
        tmp_metadata['history'] = "%sZ: Creation" % datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
        tmp_metadata['platform'] = "SWOT"
        tmp_metadata['references'] = 'Large scale simulator'
        tmp_metadata['reference_document'] = 'JPL D-56411 - Initial release - February 11, 2019'
        tmp_metadata['wavelength'] = 0.008385803020979
        self.set_metadata_val(tmp_metadata)
        
        # 3 - Update group metadata
        # 3.1 - pixel_cloud group
        pixel_cloud_metadata = {}
        pixel_cloud_metadata['description'] = 'cloud of geolocated interferogram pixels'      
        pixel_cloud_metadata['looks_to_efflooks'] = 1.75
        self.set_metadata_val(pixel_cloud_metadata, group="pixel_cloud")
        # 3.2 - tvp group
        tvp_metadata = {}
        tvp_metadata['description'] = 'Time varying parameters group including spacecraft attitude, position, velocity, and antenna position information'
        self.set_metadata_val(tvp_metadata, group="tvp")
        # 3.3 - noise group
        noise_metadata = {}
        noise_metadata['description'] = 'Measured noise power for each recieve echo of the plus_y and minus_y SLC channels'
        self.set_metadata_val(noise_metadata, group="noise")


#######################################


class PixcVecRiverProduct(NetCDFProduct):
    """
    Deal with PIXCVecRiver NetCDF file
    """
    
    def __init__(self):
        """
        Constructor; specific to PIXCVecRiver NetCDF file
        """
        
        # 1 - Init NetCDF_product
        super().__init__(os.path.expandvars("$SWOT_HYDROLOGY_TOOLBOX/sisimp/data/pixcvecriver.xml"))
                
        # 2 - Update general metadata specific to PIXC file
        tmp_metadata = {}
        tmp_metadata['Conventions'] = 'CF-1.7'
        tmp_metadata['title'] = 'Level 2 KaRIn high rate pixel cloud vector river product'
        tmp_metadata['institution'] = 'CNES - Large scale simulator'
        tmp_metadata['source'] = 'Ka-band radar interferometer'
        tmp_metadata['history'] = "%sZ: Creation" % datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
        tmp_metadata['mission_name'] = "SWOT"
        tmp_metadata['references'] = 'Large scale simulator'
        tmp_metadata['reference_document'] = 'JPL D-56415 - Initial release - February 02, 2020'
        self.set_metadata_val(tmp_metadata)


#######################################
            
            
def find_dtype(in_name, in_width, in_flag_unsigned=False):
    """
    Find the Numpy type given the type name, width and signed/unsigned information
    
    :param in_name: name of the variable
    :type in_name: string
    :param in_width: width of the type
    :type in_width: integer
    :param in_flag_unsigned: flag indicating if the type is unsigned or signed(=default)
    :type in_flag_unsigned: boolean
    
    :return: out_numpy_type = Numpy type
    :rtype: out_numpy_type = Numpy.dtype
    """
    
    # 0 - Init the output value
    out_numpy_type = None
    
    # 1 - Set the letter corresponding to the type, given its name and signed/unsigned flag
    if in_name == "integer":
        if in_flag_unsigned:
            tmp_type = "u"
        else:
            tmp_type = "i"
    elif in_name == "real":
        tmp_type = "f"
    else:
        my_api.exitWithError("[ERROR] Type %s unknown" % in_name)
    
    # 2 - Set the length
    tmp_length = str(int(in_width/8.0))
    
    # 3 - Compute the output value
    out_numpy_type = np.dtype(tmp_type+tmp_length)
    # TODO: to improve; used to remedy to the pb 
    # ProcessingError: datatype not recognized : uint64
    # that we have for reach_id and node_id
    if out_numpy_type == "u8":
        out_numpy_type = "u4"

    return out_numpy_type
