"""
.. module:: locnes_products_general.py
    :synopsis: Deals with SWOT products
     Created on 06/10/2020

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2019 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
"""
from collections import OrderedDict
import datetime
import logging
import os

import cnes.common.service_config_file as service_config_file

import cnes.common.lib.my_tools as my_tools


class LocnesProduct(object):
    """
    Main class dealing with the SWOT products in output of LOCNES 
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
        
        Variables of the object:
            - file_type / string: file type
            - global_metadata / OrderedDict: global metadata attributes; key=attribute name; 
                value=OrderedDict with (key, value)=(name, value) of each metadata
            - list_shapes_names / set: list of shape names (optionnal)
            - dims / dict: size of each dimension (optionnal)
            - attribute_metadata / OrderedDict: dictionary having key=attribute name and 
                value=OrderedDict with (key, value)=(name, value) of each metadata
        """
        # Get instance of service config file
        cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # 0 - Init variables of the object
        self.file_type = None
        self.global_metadata = OrderedDict()
        self.list_shapes_names = list()
        self.dims = dict()
        self.attribute_metadata = OrderedDict()
        
        # 1 - Set value from XML file
        self.set_from_xml(in_xml_file)
        
        # 2 - Update global metadata attributes
        self.global_metadata["institution"]["value"] = cfg.get('FILE_INFORMATION', 'PRODUCER')
        self.global_metadata["source"]["value"] = cfg.get('FILE_INFORMATION', 'SOURCE')
        self.global_metadata["history"]["value"] = "%sZ: Creation" % my_tools.swot_timeformat(datetime.datetime.utcnow(), in_format=1)
        self.global_metadata["references"]["value"] = cfg.get('FILE_INFORMATION', 'SOFTWARE_VERSION')
        self.global_metadata["contact"]["value"] = cfg.get('FILE_INFORMATION', 'CONTACT')
        
        # 3.1 - Update metadata retrieved from input data product
        if in_inprod_metadata is not None:
            self.set_metadata_val(in_inprod_metadata)
        # 3.2 - Update processing metadata
        if in_proc_metadata is not None:
            self.set_metadata_val(in_proc_metadata)
            
        # 4 - Update file path if necessary
        self.fullpath_or_basename()
        
    #----------------------------------------
        
    def set_metadata_val(self, in_metadata):
        """
        Setter of metadata value
        
        :param in_metadata: metadata stored as a dictionary with key=metadata key and value=metadata value
        :type in_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        metadata_keys = self.global_metadata.keys()
        
        for key, value in in_metadata.items():
            if key in metadata_keys:
                if value not in ["", "None", None]:
                    self.global_metadata[key]["value"] = value
            else:
                logger.debug("%s key is not used as metadata" % key)
                            
    #----------------------------------------
    
    def fullpath_or_basename(self):
        """
        If LOCNES is not run in simulation environment, convert xref_[...] global metadata full path into basename
        """
        # Get instance of service config file
        cfg = service_config_file.get_instance()
        
        if cfg.get("FILE_INFORMATION", "SOURCE") != "Simulation":
            for name, att_dict in self.global_metadata.items():
                if name.startswith("xref_"):
                    tmp_desc = os.path.basename(self.global_metadata[name]["value"])
                    self.global_metadata[name]["value"] = tmp_desc
                    