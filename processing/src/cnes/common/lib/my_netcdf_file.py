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


class myNcReader(object):

    def __init__(self, IN_filename):
        """
        Constructor
        
        :param IN_filename: filename of the NetCDF file
        :type IN_filename: string
        
        Variables:
            filename / string: filename of the NetCDF file
            content / netCDF4.Dataset: content of the NetCDF file
        """
        
        # Filename
        self.filename = IN_filename
        
        # Content
        self.content = Dataset(IN_filename, 'r')
        
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
            return self.content.dimensions
        return IN_group.dimensions
                    
    def printListDim(self, IN_group=None):
        """
        Print list of dimensions

        :param IN_group: group of variables to study
        :type IN_group: netCDF4.Group
        """
        for value in self.getListDim(IN_group=IN_group):
            print("%s - size = %d" % (value, self.getDimValue(value, IN_group=IN_group)))
    
    def getDimValue(self, IN_name, IN_group=None):
        """
        Get the size of the dimension named IN_name
        
        :param IN_name: name of the dimension
        :type IN_name: string
        :param IN_group: group containing the dimension IN_name
        :type IN_group: netCDF4.Group
        
        :return: size of the dimension
        :rtype: int
        """
        if IN_group is None:
            return len(self.content.dimensions[IN_name])
        return len(IN_group.dimensions[IN_name])
    
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
            return self.content.ncattrs()
        return IN_group.ncattrs()
        
    def printListAtt(self, IN_group=None):
        """
        Print list of global attributes

        :param IN_group: group to study
        :type IN_group: netCDF4.Group
        """
        for value in self.getListAtt(IN_group=IN_group):
            print("%s = %s" % (value, str(self.getAttValue(value, IN_group=IN_group))))
    
    def getAttValue(self, IN_name, IN_group=None):
        """
        Get the value associated to the global attribute named IN_name
        
        :param IN_name: name of the attribute
        :type IN_name: string
        :param IN_group: group containing the global attribute IN_name
        :type IN_group: netCDF4.Group
        
        :return: value associated to the attribute
        :rtype: string
        """
        if IN_group is None:
            return self.content.getncattr(IN_name)
        return IN_group.getncattr(IN_name)
    
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
    
    def getVarValue(self, IN_name, IN_group=None):
        """
        Get the data associated to the variable IN_name
        _FillValue value are filled by NaN and the multiplication by the scale_factor is done if needed
        
        :param IN_name: name of the variable
        :type IN_name: string
        :param IN_group: group containing the variable IN_name
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
        OUT_data = numpy.copy(cur_content.variables[IN_name][:])
            
        # 2 - If data values are not "int", change _FillValue by NaN
        OUT_type = str(OUT_data.dtype)
        if not OUT_type.startswith("int"):  # NaN ne s'applique pas aux tableaux d'entiers
            try:
                OUT_data[OUT_data == cur_content.variables[IN_name]._FillValue] = numpy.nan
            except AttributeError:
                pass
                # print "[my_netcdf_file] %s: no _FillValue" % ( IN_name )
                
        # 3 - Multiplication by the scale factor
        try:
            OUT_data = OUT_data * cur_content.variables[IN_name].scale_factor
        except AttributeError:
            pass
            # print('[my_netcdf_file] %s: no scale_factor') % ( IN_name )
            
        return OUT_data
    
    def getVarUnit(self, IN_name, IN_group=None):
        """
        Get the unit of variable named IN_name
        
        :param IN_name: name of the variable
        :type IN_name: string
        :param IN_group: group containing the variable IN_name
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
            return cur_content.variables[IN_name].units
        except AttributeError:
            return "-"

            
#######################################


class myNcWriter(object):

    def __init__(self, IN_filename):
        """
        Constructor
        
        :param IN_filename: full path of the NetCDF file
        :type IN_filename: string
        
        Variables:
            filename / string: full path of the NetCDF file
            content / netCDF4.Dataset: content of the NetCDF file
        """
        
        # Filename
        self.filename = IN_filename
        
        # Content
        self.content = Dataset(IN_filename, 'w')
        
    def close(self):
        """
        Close the NetCDF file
        """
        self.content.close()
    
    #----------------------------------------

    def add_group(self, IN_name):
        """
        Add a group
        
        :param IN_name: the name of the group
        :type IN_name: string
        """
        return self.content.createGroup(IN_name)
        
    def add_dimension(self, IN_name, IN_size, IN_group=None):
        """
        Set the size of the dimension with the given name

        :param IN_name: the name of the dimension
        :type IN_name: string
        :param IN_size: the size of the dimension
        :type IN_size: int
        :param IN_group: group which will contain the dimension IN_name
        :type IN_group: netCDF4.Group
        """
        if IN_group is None:
            self.content.createDimension(IN_name, IN_size)
        else:
            IN_group.createDimension(IN_name, IN_size)

    #----------------------------------------
        
    def add_global_attribute(self, IN_name, IN_value, IN_group=None):
        """
        Set the value of a global attribute with the given name

        :param IN_name: the name of the attribute
        :type IN_name: string
        :param IN_value: the value to store
        :type IN_value: unknown
        :param IN_group: group which will contain the global attribute IN_name
        :type IN_group: netCDF4.Group
        """
        if IN_group is None:
            self.content.setncattr(IN_name, IN_value)
        else:
            IN_group.setncattr(IN_name, IN_value)

    #----------------------------------------
        
    def add_variable(self, IN_name, IN_datatype, IN_dimensions, IN_fill_value=None, IN_group=None, IN_compress=True):
        """
        Add the data content of the variable
        
        :param IN_name: the name of the variable
        :type IN_name: string
        :param IN_datatype: the type of the variable
        :type IN_datatype: ex np.int32, ...
        :param IN_dimensions: the name of the dimensions of the variable
        :type IN_dimensions: string
        :param IN_fill_value: the value used when no data is saved
        :type IN_fill_value: int or float
        :param IN_group: group which will contain the variable IN_name
        :type IN_group: netCDF4.Group
        :param IN_compress: true to compress the content of the variable (save disk space), else false
        :type IN_compress: boolean
        """

        # 0 - Select the group to study
        if IN_group is None:
            cur_content = self.content
        else:
            cur_content = IN_group

        cur_content.createVariable(IN_name, IN_datatype, IN_dimensions, IN_compress, 2, fill_value=IN_fill_value)
        
    def add_variable_attribute(self, IN_varname, IN_attname, IN_value, IN_group=None):
        """
        Set the value of a variable attribute with the given name

        :param IN_varname: the name of the variable
        :type IN_varname: string
        :param IN_attname: the name of the attribute of the variable
        :type IN_attname: string
        :param IN_value: the value to store
        :type IN_value: unknown
        :param IN_group: group which will contain the variable IN_name
        :type IN_group: netCDF4.Group
        """

        # 0 - Select the group to study
        if IN_group is None:
            cur_content = self.content
        else:
            cur_content = IN_group

        if IN_varname in cur_content.variables:
            cur_content.variables[IN_varname].setncattr(IN_attname, IN_value)
        else:
            exit("[my_netcdf_file/add_variable_attribute] Could not add variable attribute ; %s variable does not exist" % IN_varname)

    #----------------------------------------
        
    def fill_variable(self, IN_name, IN_data, IN_group=None):
        """
        Write the given data into the named variable 

        :param IN_name: the name of the variable
        :type IN_name: string
        :param IN_data: the data to store
        :type IN_data: unknown
        :param IN_group: group which will contain the variable IN_name
        :type IN_group: netCDF4.Group
        """

        # 0 - Select the group to study
        if IN_group is None:
            cur_content = self.content
        else:
            cur_content = IN_group

        if IN_name in cur_content.variables:
            # Write the whole array
            cur_content.variables[IN_name][:] = IN_data
        else:
            exit("[my_netcdf_file/fill_variable] Could not fill variable ; %s variable does not exist" % IN_name)


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

