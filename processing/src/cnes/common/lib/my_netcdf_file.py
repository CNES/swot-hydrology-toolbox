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
.. module:: my_netcdf_file.py
    :synopsis: Deals with NetCDF files (reader and writer)
     Created on 2016/08/04

.. moduleauthor:: Claire POTTIER - CNES DCT/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from __future__ import absolute_import, division, print_function, unicode_literals 

from netCDF4 import Dataset
import numpy
import logging

import cnes.common.lib.my_variables as my_var
import cnes.common.service_error as service_error


class myNcReader(object):
    """
        class myNcReader
    """
    def __init__(self, in_filename, mode='r'):
        """
        Constructor
        
        :param in_filename: filename of the NetCDF file
        :type in_filename: string
        :param mode: open mode of filename
        :type mode: string


        Variables of the object:   
          - filename / string: filename of the NetCDF file
          - content / netCDF4.Dataset: content of the NetCDF file
        """
        
        # Filename
        self.filename = in_filename
        
        # Content
        self.content = Dataset(in_filename, mode)
        
    def close(self):
        """
        Close the NetCDF file
        """
        self.content.close()

    #----------------------------------------
    
    def getListDim(self, in_group=None):
        """
        Get the dictionnary of the dimensions

        :param in_group: group of variables to study
        :type in_group: netCDF4.Group

        :return: dictionnary of dimensions
        :rtype: dict
        """
        if in_group is None:
            retour = self.content.dimensions
        else:
            retour = in_group.dimensions
        return retour
                    
    def printListDim(self, in_group=None):
        """
        Print list of dimensions

        :param in_group: group of variables to study
        :type in_group: netCDF4.Group
        """
        for value in self.getListDim(in_group=in_group):
            print("%s - size = %d" % (value, self.getDimValue(value, in_group=in_group)))
    
    def getDimValue(self, in_name, in_group=None):
        """
        Get the size of the dimension named in_name
        
        :param in_name: name of the dimension
        :type in_name: string
        :param in_group: group containing the dimension in_name
        :type in_group: netCDF4.Group
        
        :return: size of the dimension
        :rtype: int
        """
        if in_group is None:
            retour = len(self.content.dimensions[in_name])
        else:
            retour = len(in_group.dimensions[in_name])
        return retour
    
    #----------------------------------------
    
    def getListAtt(self, in_group=None):
        """
        Get the list of all the global attributes included in the NetCDF file

        :param in_group: group to study
        :type in_group: netCDF4.Group
        
        :return: global attribute names
        :rtype: list of string
        """
        if in_group is None:
            retour = self.content.ncattrs()
        else:
            retour = in_group.ncattrs()
        return retour
        
    def printListAtt(self, in_group=None):
        """
        Print list of global attributes

        :param in_group: group to study
        :type in_group: netCDF4.Group
        """
        for value in self.getListAtt(in_group=in_group):
            print("%s = %s" % (value, str(self.getAttValue(value, in_group=in_group))))
    
    def getAttValue(self, in_name, in_group=None):
        """
        Get the value associated to the global attribute named in_name
        
        :param in_name: name of the attribute
        :type in_name: string
        :param in_group: group containing the global attribute in_name
        :type in_group: netCDF4.Group
        
        :return: value associated to the attribute
        :rtype: string
        """
        if in_group is None:
            retour = self.content.getncattr(in_name)
        else:
            retour = in_group.getncattr(in_name)
        return retour
    
    #----------------------------------------
    
    def getListVar(self, in_group=None):
        """
        Get the list of variables

        :param in_group: group to study
        :type in_group: netCDF4.Group
        
        :return:out_list_var = global or group variables
        :rtype: list of string
        """
        if in_group is None:
            out_list_var = sorted(self.content.variables.keys())
        else:
            out_list_var = sorted(in_group.variables.keys())
            
        return out_list_var
        
    def printListVar(self, in_group=None):
        """
        Print list of variables

        :param in_group: group to study
        :type in_group: netCDF4.Group
        """

        # 1 - Get list of variables
        list_var = self.getListVar(in_group=in_group)

        # 2 - Print list
        for value in list_var:
            print(value + " - units = " + self.getVarUnit(value, in_group=in_group))
    
    def getVarValue(self, in_name, in_group=None):
        """
        Get the data associated to the variable in_name
        The multiplication by the scale_factor is done if there is a scale_factor attribute 
        _FillValue values are let as is (not filled by NaN)
        
        :param in_name: name of the variable
        :type in_name: string
        :param in_group: group containing the variable in_name
        :type in_group: netCDF4.Group
        
        :return: out_data = formatted data
        :rtype: [depend on the inumpy.t variable]
        """
        logger = logging.getLogger(self.__class__.__name__)

        # 0 - Select the group to study
        if in_group is None:
            cur_content = self.content
        else:
            cur_content = in_group
        
        # 1 - Get data
        try:
            out_data = numpy.copy(cur_content.variables[in_name][:])
        except KeyError:
            message = "Variable %s does not exist in netCDF file" % in_name
            raise service_error.ProcessingError(message, logger)
        # TODO: remove when inumpy.t variables corrected (no NaN values anymore)
        out_type = str(out_data.dtype)
        if out_type.startswith("float") or out_type.startswith("double"):
            nan_idx = numpy.argwhere(numpy.isnan(out_data))
            nb_nan = len(nan_idx)
            if nb_nan > 0:
                try:
                    logger.warning("{} NaN values remaining in {} variable => replaced by {}".format(nb_nan, in_name, cur_content.variables[in_name]._FillValue))
                    out_data[nan_idx] = cur_content.variables[in_name]._FillValue
                except AttributeError:
                    logger.warning("{} NaN values remaining in {} variable => replaced by {} (_FillValue unknown)".format(nb_nan, in_name, my_var.FV_NETCDF[out_type]))
                    out_data[nan_idx] = my_var.FV_NETCDF[out_type]
        # == END-TODO ==
        
        # 2 - Get not _FillValue indices
        fv_flag = False
        # TODO Change that try/except/pass without error class is not allowed
        try:
            not_nan_idx = numpy.where(out_data < cur_content.variables[in_name]._FillValue)
            if len(not_nan_idx) < len(out_data):
                fv_flag = True
        except:
            pass
                
        # 3 - Multiplication by the scale factor where not NaN
        try:
            if fv_flag:
                out_data[not_nan_idx] *= cur_content.variables[in_name].scale_factor
            else:
                out_data *= cur_content.variables[in_name].scale_factor
        except:
            pass
            
        return out_data
    
    def getVarValue_orEmpty(self, in_name, in_group=None):
        """
        Get the data associated to the variable in_name if it exists
        If not, return None
        
        :param in_name: name of the variable
        :type in_name: string
        :param in_group: group containing the variable in_name
        :type in_group: netCDF4.Group
        
        :return: out_data = formatted data or None
        :rtype: [depend on the inumpy.t variable]
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        list_vars = self.getListVar(in_group=in_group)
        
        if in_name in list_vars:
            out_data = self.getVarValue(in_name, in_group)
        else:
            logger.info("Variable %s doesn't exit => return empty array" % in_name)
            # Get all dimensions
            list_dims = self.getListDim(in_group=in_group)
            # Compute greatest one
            nb_records = 0
            for key, value in list_dims.items():
                if value.size > nb_records:
                    nb_records = value.size
            # Make empty array of int with _FillValue
            out_data = numpy.zeros(nb_records, dtype=numpy.int32) + 2147483647
            
        return out_data
    
    def getVarUnit(self, in_name, in_group=None):
        """
        Get the unit of variable named in_name
        
        :param in_name: name of the variable
        :type in_name: string
        :param in_group: group containing the variable in_name
        :type in_group: netCDF4.Group
        
        :return: unit
        :rtype: string
        """

        # 0 - Select the group to study
        if in_group is None:
            cur_content = self.content
        else:
            cur_content = in_group

        try:
            retour = cur_content.variables[in_name].units
        except AttributeError:
            retour = "-"

        return retour
    
            
#######################################


class myNcWriter(object):
    """
        class myNcWriter
    """
    def __init__(self, in_filename):
        """
        Constructor
        
        :param in_filename: full path of the NetCDF file
        :type in_filename: string
        

        Variables of the object:
         - filename / string: full path of the NetCDF file
         - content / netCDF4.Dataset: content of the NetCDF file
        """
        
        # Filename
        self.filename = in_filename
        
        # Content
        self.content = Dataset(in_filename, 'w')
        
    def close(self):
        """
        Close the NetCDF file
        """
        self.content.close()
    
    #----------------------------------------

    def add_group(self, in_name):
        """
        Add a group
        
        :param in_name: the name of the group
        :type in_name: string
        """
        return self.content.createGroup(in_name)
        
    def add_dimension(self, in_name, in_size, in_group=None):
        """
        Set the size of the dimension with the given name

        :param in_name: the name of the dimension
        :type in_name: string
        :param in_size: the size of the dimension
        :type in_size: int
        :param in_group: group which will contain the dimension in_name
        :type in_group: netCDF4.Group
        """
        if in_group is None:
            self.content.createDimension(in_name, in_size)
        else:
            in_group.createDimension(in_name, in_size)

    #----------------------------------------
        
    def add_global_attribute(self, in_name, in_value, in_group=None):
        """
        Set the value of a global attribute with the given name

        :param in_name: the name of the attribute
        :type in_name: string
        :param in_value: the value to store
        :type in_value: unknown
        :param in_group: group which will contain the global attribute in_name
        :type in_group: netCDF4.Group
        """
        
        # Fix python_netcdf4 bug with in_value.encode('ascii')
        if in_group is None:
            if type(in_value) is str:
                self.content.setncattr(in_name, in_value.encode('utf-8'))
            else:
                self.content.setncattr(in_name, in_value)
        else:
            if type(in_value) is str:
                in_group.setncattr(in_name, in_value.encode('utf-8'))
            else:
                in_group.setncattr(in_name, in_value)

    #----------------------------------------
        
    def add_variable(self, in_name, in_datatype, in_dimensions, in_group=None, in_attributes=None, in_compress=True):
        """
        Add the data content of the variable
        
        :param in_name: the name of the variable
        :type in_name: string
        :param in_datatype: the type of the variable
        :type in_datatype: ex numpy.int32, ...
        :param in_dimensions: the name of the dimensions of the variable
        :type in_dimensions: string
        :param in_group: group which will contain the variable in_name
        :type in_group: netCDF4.Group
        :param in_attributes: variable attributes stored as a dictionary key=attribute name value=attribute value
        :type in_attributes: dict
        :param in_compress: true to compress the content of the variable (save disk space), else false
        :type in_compress: boolean
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        # Select the group to study
        if in_group is None:
            cur_content = self.content
        else:
            cur_content = in_group
        """        
<<<<<<< HEAD
        if (in_datatype is numpy.double) or (in_datatype is numpy.float64):
            # double netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, in_compress, 2, fill_value=my_variables.FV_DOUBLE)
        elif (in_datatype is numpy.float):
            # float netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, in_compress, 2, fill_value=my_variables.FV_FLOAT)
        elif (in_datatype is numpy.int32) or (in_datatype is int):
            # int netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, in_compress, 2, fill_value=my_variables.FV_INT)
        elif (in_datatype is numpy.uint32):
            # unsigned int netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, in_compress, 2, fill_value=my_variables.FV_UINT)
        elif (in_datatype is numpy.int16):
            # short netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, in_compress, 2, fill_value=my_variables.FV_SHORT)
        elif (in_datatype is numpy.uint16):
            # unsigned short netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, in_compress, 2, fill_value=my_variables.FV_USHORT)
        elif (in_datatype is numpy.int8):
            # byte netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, in_compress, 2, fill_value=my_variables.FV_BYTE)
        elif (in_datatype is numpy.uint8):
            # unsigned byte netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, in_compress, 2, fill_value=my_variables.FV_UBYTE)
        elif (in_datatype is str):
            # char netCDF variable command below is not correctly support in python netCDF...
            #cur_content.createVariable(in_name, in_datatype, in_dimensions, in_compress, 2, fill_value=my_variables.FV_STRING)
            cur_content.createVariable(in_name, in_datatype, in_dimensions, in_compress, 2)
=======
        """
        # Create variable depending on its type
        if in_datatype is str:
            # string type is not allowed in netCDF, force NC_CHAR type
            #cur_content.createVariable(in_name, 'c', in_dimensions, in_compress, 2)
            cur_content.createVariable(in_name, in_datatype, in_dimensions, in_compress, 2)
        elif numpy.dtype(in_datatype).name in my_var.FV_NETCDF:
            cur_content.createVariable(in_name, in_datatype, in_dimensions, in_compress, 2, fill_value=my_var.FV_NETCDF[numpy.dtype(in_datatype).name])
        else:
            # datatype not recognized !
            message = "datatype not recognized : %s %s" % str(numpy.dtype(in_datatype).name)
            print(str(my_var.FV_NETCDF))
            raise service_error.ProcessingError(message, logger)
            
        # Add variable attributes
        if in_attributes is not None:
            for key, value in in_attributes.items():
                if key != "dtype":
                    self.add_variable_attribute(in_name, key, value, in_group=in_group)
        
    def add_variable_attribute(self, in_varname, in_attname, in_value, in_group=None):
        """
        Set the value of a variable attribute with the given name

        :param in_varname: the name of the variable
        :type in_varname: string
        :param in_attname: the name of the attribute of the variable
        :type in_attname: string
        :param in_value: the value to store
        :type in_value: unknown
        :param in_group: group which will contain the variable in_name
        :type in_group: netCDF4.Group
        """
        logger = logging.getLogger(self.__class__.__name__)

        # 0 - Select the group to study
        if in_group is None:
            cur_content = self.content
        else:
            cur_content = in_group

        if in_varname in cur_content.variables:
            if type(in_value) is str:
                cur_content.variables[in_varname].setncattr(in_attname, in_value.encode('ascii'))
            else:
                cur_content.variables[in_varname].setncattr(in_attname, in_value)
        else:
            # Variable not recognized !
            logger.debug("Could not add variable attribute %s because %s variable does not exist" % (in_attname, in_varname))

    #----------------------------------------
        
    def fill_variable(self, in_name, in_data, in_group=None):
        """
        Write the given data into the named variable 

        :param in_name: the name of the variable
        :type in_name: string
        :param in_data: the data to store
        :type in_data: unknown
        :param in_group: group which will contain the variable in_name
        :type in_group: netCDF4.Group
        """
        logger = logging.getLogger(self.__class__.__name__)

        # 0 - Select the group to study
        if in_group is None:
            cur_content = self.content
        else:
            cur_content = in_group

        if in_name in cur_content.variables:
            # Test existence of NaNs and replace them by _FillValue
            data_type = str(in_data.dtype)
            if data_type.startswith("float") or data_type.startswith("double"):
                nan_idx = numpy.argwhere(numpy.isnan(in_data))
                nb_nan = len(nan_idx)
                if nb_nan > 0:
                    try:
                        logger.warning("{} NaN values remaining in {} variable => replaced by {}".format(nb_nan, in_name, cur_content.variables[in_name]._FillValue))
                        in_data[nan_idx] = cur_content.variables[in_name]._FillValue
                    except AttributeError:
                        logger.warning("{} NaN values remaining in {} variable => replaced by {} (_FillValue unknown)".format(nb_nan, in_name, my_var.FV_NETCDF[data_type]))
                        in_data[nan_idx] = my_var.FV_NETCDF[data_type]
            # Write the whole array
            cur_content.variables[in_name][:] = in_data
        else:
            # Variable not recognized !
            logger.debug("Could not fill variable %s because it does not exist" % in_name)
            
