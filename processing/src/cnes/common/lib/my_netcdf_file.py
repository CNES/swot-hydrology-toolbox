# -*- coding: utf8 -*-
"""
.. module:: my_netcdf_file.py
    :synopsis: Deals with NetCDF files (reader and writer)
    Created on 08/04/2016

.. moduleauthor:: Claire POTTIER - CNES DCT/SI/TR

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from __future__ import absolute_import, division, print_function, unicode_literals 

from netCDF4 import Dataset
import numpy
import logging

import cnes.common.lib.my_variables as my_variables
import cnes.common.service_error as service_error

class myNcReader(object):

    def __init__(self, in_filename):
        """
        Constructor
        
        :param in_filename: filename of the NetCDF file
        :type in_filename: string
        
        Variables:
            filename / string: filename of the NetCDF file
            content / netCDF4.Dataset: content of the NetCDF file
        """
        
        # Filename
        self.filename = in_filename
        
        # Content
        self.content = Dataset(in_filename, 'r')
        
    def close(self):
        """
        Close the NetCDF file
        """
        self.content.close()

    #----------------------------------------
    
    def getListDim(self, IN_group=None):
        """
        Get the dictionnary of the dimensions

        :param IN_group: group of variables to study
        :type IN_group: netCDF4.Group

        :return: dictionnary of dimensions
        :rtype: dict
        """
        if IN_group is None:
            retour = self.content.dimensions
        else:
            retour = IN_group.dimensions
        return retour
                    
    def printListDim(self, IN_group=None):
        """
        Print list of dimensions

        :param IN_group: group of variables to study
        :type IN_group: netCDF4.Group
        """
        for value in self.getListDim(IN_group=IN_group):
            print("%s - size = %d" % (value, self.getDimValue(value, IN_group=IN_group)))
    
    def getDimValue(self, in_name, IN_group=None):
        """
        Get the size of the dimension named in_name
        
        :param in_name: name of the dimension
        :type in_name: string
        :param IN_group: group containing the dimension in_name
        :type IN_group: netCDF4.Group
        
        :return: size of the dimension
        :rtype: int
        """
        if IN_group is None:
            retour = len(self.content.dimensions[in_name])
        else:
            retour = len(IN_group.dimensions[in_name])
        return retour
    
    #----------------------------------------
    
    def getListAtt(self, IN_group=None):
        """
        Get the list of all the global attributes included in the NetCDF file

        :param IN_group: group to study
        :type IN_group: netCDF4.Group
        
        :return: global attribute names
        :rtype: list of string
        """
        if IN_group is None:
            retour = self.content.ncattrs()
        else:
            retour = IN_group.ncattrs()
        return retour
        
    def printListAtt(self, IN_group=None):
        """
        Print list of global attributes

        :param IN_group: group to study
        :type IN_group: netCDF4.Group
        """
        for value in self.getListAtt(IN_group=IN_group):
            print("%s = %s" % (value, str(self.getAttValue(value, IN_group=IN_group))))
    
    def getAttValue(self, in_name, IN_group=None):
        """
        Get the value associated to the global attribute named in_name
        
        :param in_name: name of the attribute
        :type in_name: string
        :param IN_group: group containing the global attribute in_name
        :type IN_group: netCDF4.Group
        
        :return: value associated to the attribute
        :rtype: string
        """
        if IN_group is None:
            retour = self.content.getncattr(in_name)
        else:
            retour = IN_group.getncattr(in_name)
        return retour
    
    #----------------------------------------
        
    def printListVar(self, IN_group=None):
        """
        Print list of variables

        :param IN_group: group to study
        :type IN_group: netCDF4.Group
        """

        # 1 - Get list of variables
        if IN_group is None:
            list_var = sorted(self.content.variables.keys())
        else:
            list_var = sorted(IN_group.variables.keys())

        # 2 - Print list
        for value in list_var:
            print(value + " - units = " + self.getVarUnit(value, IN_group=IN_group))
    
    def getVarValue(self, in_name, IN_group=None):
        """
        Get the data associated to the variable in_name
        _FillValue value are filled by NaN and the multiplication by the scale_factor is done if needed
        
        :param in_name: name of the variable
        :type in_name: string
        :param IN_group: group containing the variable in_name
        :type IN_group: netCDF4.Group
        
        :return: formatted data
        :rtype: [depend on the input variable]
        """

        # 0 - Select the group to study
        if IN_group is None:
            cur_content = self.content
        else:
            cur_content = IN_group
        
        # 1 - Get data
        out_data = numpy.copy(cur_content.variables[in_name][:])
            
        # 2 - If data values are not "int", change _FillValue by NaN
        out_type = str(out_data.dtype)
        if not out_type.startswith("int"):  # NaN ne s'applique pas aux tableaux d'entiers
            try:
                out_data[out_data == cur_content.variables[in_name]._FillValue] = numpy.nan
            except AttributeError:
                pass
                # print "[my_netcdf_file] %s: no _FillValue" % ( in_name )
                
        # 3 - Multiplication by the scale factor
        try:
            out_data = out_data * cur_content.variables[in_name].scale_factor
        except AttributeError:
            pass
            # print('[my_netcdf_file] %s: no scale_factor') % ( in_name )
            
        return out_data
    
    def getVarUnit(self, in_name, IN_group=None):
        """
        Get the unit of variable named in_name
        
        :param in_name: name of the variable
        :type in_name: string
        :param IN_group: group containing the variable in_name
        :type IN_group: netCDF4.Group
        
        :return: unit
        :rtype: string
        """

        # 0 - Select the group to study
        if IN_group is None:
            cur_content = self.content
        else:
            cur_content = IN_group

        try:
            retour = cur_content.variables[in_name].units
        except AttributeError:
            retour = "-"

        return retour
            
#######################################


class myNcWriter(object):

    def __init__(self, in_filename):
        """
        Constructor
        
        :param in_filename: full path of the NetCDF file
        :type in_filename: string
        
        Variables:
            filename / string: full path of the NetCDF file
            content / netCDF4.Dataset: content of the NetCDF file
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
        
    def add_dimension(self, in_name, in_size, IN_group=None):
        """
        Set the size of the dimension with the given name

        :param in_name: the name of the dimension
        :type in_name: string
        :param in_size: the size of the dimension
        :type in_size: int
        :param IN_group: group which will contain the dimension in_name
        :type IN_group: netCDF4.Group
        """
        if IN_group is None:
            self.content.createDimension(in_name, in_size)
        else:
            IN_group.createDimension(in_name, in_size)

    #----------------------------------------
        
    def add_global_attribute(self, in_name, in_value, IN_group=None):
        """
        Set the value of a global attribute with the given name

        :param in_name: the name of the attribute
        :type in_name: string
        :param in_value: the value to store
        :type in_value: unknown
        :param IN_group: group which will contain the global attribute in_name
        :type IN_group: netCDF4.Group
        """
        if IN_group is None:
            self.content.setncattr(in_name, in_value)
        else:
            IN_group.setncattr(in_name, in_value)

    #----------------------------------------
        
    def add_variable(self, in_name, in_datatype, in_dimensions, IN_fill_value=None, IN_group=None, IN_compress=True):
        """
        Add the data content of the variable
        
        :param in_name: the name of the variable
        :type in_name: string
        :param in_datatype: the type of the variable
        :type in_datatype: ex np.int32, ...
        :param in_dimensions: the name of the dimensions of the variable
        :type in_dimensions: string
        :param IN_fill_value: the value used when no data is saved
        :type IN_fill_value: int or float
        :param IN_group: group which will contain the variable in_name
        :type IN_group: netCDF4.Group
        :param IN_compress: true to compress the content of the variable (save disk space), else false
        :type IN_compress: boolean
        """
        
        logger = logging.getLogger(self.__class__.__name__)
        # 0 - Select the group to study
        if IN_group is None:
            cur_content = self.content
        else:
            cur_content = IN_group
        
        if (in_datatype is numpy.double) or (in_datatype is numpy.float64):
            # double netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, IN_compress, 2, fill_value=my_variables.FV_DOUBLE)
        elif (in_datatype is numpy.float):
            # float netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, IN_compress, 2, fill_value=my_variables.FV_FLOAT)
        elif (in_datatype is numpy.int32) or (in_datatype is int):
            # int netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, IN_compress, 2, fill_value=my_variables.FV_INT)
        elif (in_datatype is numpy.uint32):
            # unsigned int netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, IN_compress, 2, fill_value=my_variables.FV_UINT)
        elif (in_datatype is numpy.int16):
            # short netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, IN_compress, 2, fill_value=my_variables.FV_SHORT)
        elif (in_datatype is numpy.uint16):
            # unsigned short netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, IN_compress, 2, fill_value=my_variables.FV_USHORT)
        elif (in_datatype is numpy.int8):
            # byte netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, IN_compress, 2, fill_value=my_variables.FV_BYTE)
        elif (in_datatype is numpy.uint8):
            # unsigned byte netCDF variable
            cur_content.createVariable(in_name, in_datatype, in_dimensions, IN_compress, 2, fill_value=my_variables.FV_UBYTE)
        elif (in_datatype is str):
            # char netCDF variable command below is not correctly support in python netCDF...
            #cur_content.createVariable(in_name, in_datatype, in_dimensions, IN_compress, 2, fill_value=my_variables.FV_STRING)
            cur_content.createVariable(in_name, in_datatype, in_dimensions, IN_compress, 2)
        else:
            # datatype not recognized !
            message = "datatype not recognized : %s" % str(in_datatype)
            raise service_error.ProcessingError(message, logger)
        
        
        
    def add_variable_attribute(self, in_varname, in_attname, in_value, IN_group=None):
        """
        Set the value of a variable attribute with the given name

        :param in_varname: the name of the variable
        :type in_varname: string
        :param in_attname: the name of the attribute of the variable
        :type in_attname: string
        :param in_value: the value to store
        :type in_value: unknown
        :param IN_group: group which will contain the variable in_name
        :type IN_group: netCDF4.Group
        """

        # 0 - Select the group to study
        if IN_group is None:
            cur_content = self.content
        else:
            cur_content = IN_group

        if in_varname in cur_content.variables:
            cur_content.variables[in_varname].setncattr(in_attname, in_value)
        else:
            exit("[my_netcdf_file/add_variable_attribute] Could not add variable attribute ; %s variable does not exist" % in_varname)

    #----------------------------------------
        
    def fill_variable(self, in_name, in_data, IN_group=None):
        """
        Write the given data into the named variable 

        :param in_name: the name of the variable
        :type in_name: string
        :param in_data: the data to store
        :type in_data: unknown
        :param IN_group: group which will contain the variable in_name
        :type IN_group: netCDF4.Group
        """

        # 0 - Select the group to study
        if IN_group is None:
            cur_content = self.content
        else:
            cur_content = IN_group

        if in_name in cur_content.variables:
            # Write the whole array
            cur_content.variables[in_name][:] = in_data
        else:
            exit("[my_netcdf_file/fill_variable] Could not fill variable ; %s variable does not exist" % in_name)


#######################################


if __name__ == "__main__":
    
    pixc_main = myNcReader("C:\\Users\\pottierc\\Documents\\workspace_eclipse\\swot_py\\res\\sisimp_real_pixc\\outputs\\SWOT_L2_HR_PIXC_001_060_43N-L_main.nc")
    #pixc_main.printListDim()
    #pixc_main.printListAtt()
    #pixc_main.printListVar()
    #"""
    print("== Retrieve some data ==")
    az_index = pixc_main.getVarValue('azimuth_index')
    print("azimuth_index[0] = %d" % az_index[0])
    r_index = pixc_main.getVarValue("range_index")
    print("range_index[0] = %d" % r_index[0])
    pix_area = pixc_main.getVarValue('pixel_area')
    print("pixel_area[0] = %.6f" % pix_area[0])
    #"""

