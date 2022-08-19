# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:3.0.0:DM:#91:2021/03/12:Poursuite industrialisation
# VERSION:3.2.0:DM:#91:2021/10/27:Poursuite industrialisation
# VERSION:4.0.0:DM:#91:2022/05/05:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
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
import logging
import os

import cnes.common.service_config_file as service_config_file

import cnes.common.lib.my_variables as my_var


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
        list_fileinfo = cfg.options("FILE_INFORMATION")
        if "institution" in list_fileinfo:
            self.global_metadata["institution"]["value"] = cfg.get('FILE_INFORMATION', 'INSTITUTION')
        if "product_version" in list_fileinfo:
            self.global_metadata["product_version"]["value"] = cfg.get('FILE_INFORMATION', 'PRODUCT_VERSION')
        if "crid" in list_fileinfo:
            self.global_metadata["crid"]["value"] = cfg.get("FILE_INFORMATION", "CRID")
        if "crid_laketile" in list_fileinfo:
            self.global_metadata["crid"]["value"] = cfg.get("FILE_INFORMATION", "CRID_LAKETILE")
        if "crid_lakesp" in list_fileinfo:
            self.global_metadata["crid"]["value"] = cfg.get("FILE_INFORMATION", "CRID_LAKESP")
        if "crid_lakeavg" in list_fileinfo:
            self.global_metadata["crid"]["value"] = cfg.get("FILE_INFORMATION", "CRID_LAKEAVG")
        if "pge_version" in list_fileinfo:
            self.global_metadata["pge_version"]["value"] = cfg.get('FILE_INFORMATION', 'PGE_VERSION')
        if "contact" in list_fileinfo:
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
                if value in ["", "None", None]:
                    #self.global_metadata[key]["value"] = "None"
                    self.global_metadata[key]["value"] = ""
                else:
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
        flag_write_full_path = cfg.getboolean("OPTIONS", "Write full path")
        
        if not flag_write_full_path:
            for name, att_dict in self.global_metadata.items():
                if name.startswith("xref_"):
                    if ";" in self.global_metadata[name]["value"]:
                        tmp_split = self.global_metadata[name]["value"].split(";")
                        tmp_set = set()
                        for tmp_fullpath in tmp_split:
                            tmp_set.add(os.path.basename(tmp_fullpath))
                        self.global_metadata[name]["value"] = my_var.SEP_FILES.join(tmp_set)
                    else:
                        self.global_metadata[name]["value"] = os.path.basename(self.global_metadata[name]["value"])
                    
