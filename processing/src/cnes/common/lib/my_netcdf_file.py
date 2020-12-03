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

import logging
import netCDF4
from netCDF4 import Dataset
import numpy
import re
import cnes.common.service_error as service_error

import cnes.common.lib.my_variables as my_var


class MyNcReader(object):
    """
        class MyNcReader
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
    
    def get_list_dim(self, in_group=None):
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
                    
    def print_list_dim(self, in_group=None):
        """
        Print list of dimensions

        :param in_group: group of variables to study
        :type in_group: netCDF4.Group
        """
        for value in self.get_list_dim(in_group=in_group):
            print("%s - size = %d" % (value, self.get_dim_value(value, in_group=in_group)))
    
    def get_dim_value(self, in_name, in_group=None):
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
    
    def get_list_att(self, in_group=None):
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
        
    def print_list_att(self, in_group=None):
        """
        Print list of global attributes

        :param in_group: group to study
        :type in_group: netCDF4.Group
        """
        for value in self.get_list_att(in_group=in_group):
            print("%s = %s" % (value, str(self.get_att_value(value, in_group=in_group))))
    
    def get_att_value(self, in_name, in_group=None):
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
    
    def get_list_var(self, in_group=None):
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
        
    def print_list_var(self, in_group=None):
        """
        Print list of variables

        :param in_group: group to study
        :type in_group: netCDF4.Group
        """

        # 1 - Get list of variables
        list_var = self.get_list_var(in_group=in_group)

        # 2 - Print list
        for value in list_var:
            print(value + " - units = " + self.get_var_unit(value, in_group=in_group))
    
    def get_var_value(self, in_name, in_group=None):
        """
        Get the data associated to the variable in_name
        _FillValue values are converted to numpy.nan
        The multiplication by the scale_factor is done if there is a scale_factor attribute 
        
        :param in_name: name of the variable
        :type in_name: string
        :param in_group: group containing the variable in_name (optionnal)
        :type in_group: netCDF4.Group
        
        :return: out_data = formatted data
        :rtype: numpy.array
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
            message = "Variable %s does not exist in NetCDF file" % in_name
            raise service_error.ProcessingError(message, logger)
        
        # 2 - !!! SHOULD NOT OCCUR => IN CASE _FillValue != my_var.FV_NETCDF
        # Replace FillValue from NetCDF file with FillValue from my_variables
        type_out_data = str(out_data.dtype)
        try:
            # Test if FillValue are different
            flag_nok = False
            if type_out_data.startswith("float"):  
                # Particular case of float type because it is loaded in double variable. 
                # The difference test must be truncated.
                if abs(cur_content.variables[in_name]._FillValue - my_var.FV_NETCDF[type_out_data])/1e30 > .1:
                    flag_nok = True
            else:
                if cur_content.variables[in_name]._FillValue != my_var.FV_NETCDF[type_out_data]:
                    flag_nok = True
            # Replace if needed
            if flag_nok:
                logger.warning("FillValue for type %s from NetCDF file is different from my_variable.FV_NETCDF: %s w.r.t. %s " \
                                   %(type_out_data, str(cur_content.variables[in_name]._FillValue), str(my_var.FV_NETCDF[type_out_data])))
                out_data[numpy.where(out_data == cur_content.variables[in_name]._FillValue)] = my_var.FV_NETCDF[type_out_data]
        except:
            pass
        
        # 3 - Replace FillValue with numpy.nan
        # NB: use of NaN numpy library after in the code
        if type_out_data.startswith("float") or type_out_data.startswith("double"):
            fillvalue_idx = numpy.where(out_data > 1e30)[0]
            nb_fillvalue = len(fillvalue_idx)
            if nb_fillvalue > 0:
                logger.warning("{} _FillValue in variable {} => replaced by numpy.nan".format(nb_fillvalue, in_name))
                out_data[fillvalue_idx] = numpy.nan
                
        # 4 - Multiplication by the scale factor if it exists
        try:
            scale_factor = cur_content.variables[in_name].scale_factor
            logger.info("Apply scale_factor={} for variable {}" % (scale_factor, in_name))
            out_data *= scale_factor
        except:
            logger.info("No scale_factor for variable %s" % in_name)

        return out_data
    
    def get_var_unit(self, in_name, in_group=None):
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


class MyNcWriter(object):
    """
        class MyNcWriter
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
            
        # Create variable depending on its type
        if numpy.dtype(in_datatype).char == 'S' or numpy.dtype(in_datatype).char == 'U' or numpy.dtype(in_datatype).char == 'c':
            cur_content.createVariable(in_name, 'c', in_dimensions, in_compress, 2)
        elif numpy.dtype(in_datatype).name in my_var.FV_NETCDF:
            cur_content.createVariable(in_name, in_datatype, in_dimensions, in_compress, 2, \
                                       fill_value=my_var.FV_NETCDF[numpy.dtype(in_datatype).name])
        else:
            # datatype not recognized !
            message = "datatype not recognized : %s" % str(numpy.dtype(in_datatype).name)
            raise service_error.ProcessingError(message, logger)
            
        # Add variable attributes
        if in_attributes is not None:
            for key, value in in_attributes.items():
                if key not in ["dtype", "shape", "type", "width", "signed", "value", "_FillValue"]:
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
                tmp = re.sub(' +', ' ', in_value)
                cur_content.variables[in_varname].setncattr(in_attname, tmp.encode('ascii'))
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
            data_type = str(in_data.dtype)
            
            # == Case of float and double
            # Test existence of NaNs and replace them by _FillValue
            if data_type.startswith("float") or data_type.startswith("double"):
                nan_idx = numpy.argwhere(numpy.isnan(in_data))
                nb_nan = len(nan_idx)
                if nb_nan > 0:
                    try:
                        logger.warning("{} NaN values remaining in {} variable => replaced by {}".format(nb_nan, in_name, \
                                                                                                         cur_content.variables[in_name]._FillValue))
                        in_data[nan_idx] = cur_content.variables[in_name]._FillValue
                    except AttributeError:
                        logger.warning("{} NaN values remaining in {} variable => replaced by {} (_FillValue unknown)".format(nb_nan, in_name, \
                                                                                                                         my_var.FV_NETCDF[data_type]))
                        in_data[nan_idx] = my_var.FV_NETCDF[data_type]
            
            # == In case of table of char
            if data_type.startswith("|S") or data_type.startswith("<U") or data_type.startswith(">U") or data_type.startswith("|U"):
                # Data needs to be converted before being written in NetCDF file
                cur_content.variables[in_name][:] = netCDF4.stringtochar(in_data.astype('S'))
            elif len(in_data.shape) == 2:
                cur_content.variables[in_name][:,:] = in_data
            else:
                cur_content.variables[in_name][:] = in_data
                
        else:
            # Variable not recognized !
            logger.warning("Could not fill variable %s because it does not exist" % in_name)
            