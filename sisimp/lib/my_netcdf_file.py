# -*- coding: utf8 -*-
'''
.. module my_netcdf_file.py
    :synopsis: Deals with NetCDF files (reader and writer)
    Created on 08/04/2016

.. module author: Claire POTTIER - CNES DCT/SI/TR

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


'''
from __future__ import absolute_import, division, print_function, unicode_literals 

from netCDF4 import Dataset
import numpy

class myNcReader():

    def __init__(self, IN_filename, dim=1):
        '''
        Constructor
        
        :param IN_filename: filename of the NetCDF file
        :type IN_filename: string
        
        Variables:
            filename / string: filename of the NetCDF file
            content / netCDF4.Dataset: content of the NetCDF file
        '''
        
        # Filename
        self.filename = IN_filename
        
        # Content
        if dim == 1 :
            self.content = Dataset(IN_filename, 'r')
        elif dim == 2 :
            print("Ouverture fichier 2D geopanda")


            self.content = xr.open_dataset(IN_filename)
            self.content.attrs
        else :
            raise("Dimension trop grande")
        
    def close(self):
        '''
        Close the NetCDF file
        '''
        self.content.close()
    
    #----------------------------------------
    
    def getListDim(self, group = None):
        '''
        Get the dictionnary of the dimensions

        :return: dictionnary of dimensions
        :rtype: dict
        '''
        if group == None:
            return self.content.dimensions
        else:
            return group.dimensions
                    
    def printListDim(self, group = None):
        '''
        Print list of dimensions
        '''
        
        if group == None:
            for value in self.getListDim():
                print(value + " - size = " + str(self.getDimValue(value)))
        else:
            for value in self.getListDim(group = group):
                print(value + " - size = " + str(self.getDimValue(value, group = group)))        
    
    def getDimValue(self, IN_name, group = None):
        '''
        Get the size of the dimension named IN_name
        
        :param IN_name: name of the dimension
        :type IN_name: string
        
        :return: size of the dimension
        :rtype: int
        '''
        if group == None:
            return len(self.content.dimensions[IN_name])
        else:
            return len(group.dimensions[IN_name])
               
    
    #----------------------------------------
    
    def getListAtt(self, group = None):
        '''
        Get the list of all the global attributes included in the NetCDF file
        
        :return: global attribute names
        :rtype: list of string
        '''
        if group == None:
            return self.content.ncattrs()
        else:
            return group.ncattrs()
            
        
    def printListAtt(self, group = None):
        '''
        Print list of global attributes
        '''
        
        if group == None:
            for value in self.getListAtt():
                print(value + " = " + self.getAttValue(value))
        else:
            for value in self.getListAtt(group = group):
                print(value + " = " + self.getAttValue(value, group = group))            
    
    def getAttValue(self, IN_name, group = None):
        '''
        Get the value associated to the global attribute named IN_name
        
        :param IN_name: name of the attribute
        :type IN_name: string
        
        :return: value associated to the attribute
        :rtype: string
        '''
        if group == None:
            return self.content.getncattr(IN_name)
        else:
            return group.getncattr(IN_name)
            
    
    #----------------------------------------
        
    def printListVar(self, group = None):
        '''
        Print list of variables
        '''
        
        if group == None:
            list_var = sorted(self.content.variables.keys())
            for value in list_var:
                print(value + " - units = " + self.getVarUnit(value))
                
        else:
            list_var = sorted(group.variables.keys())
            for value in list_var:
                print(value + " - units = " + self.getVarUnit(value, group = group))                
    
    def getVarValue(self, IN_name, group = None):
        '''
        Get the data associated to the variable IN_name
        _FillValue value are filled by NaN and the multiplication by the scale_factor is done if needed
        
        :param IN_name: name of the variable
        :type IN_name: string
        
        :return: formatted data
        :rtype: [depend on the input variable]
        '''
        
        if group == None:
            
            # 1 - Get data
            OUT_data = numpy.copy(self.content.variables[IN_name][:])
            
            # 2 - If data values are not "int", change _FillValue by NaN
            OUT_type = str(OUT_data.dtype)
            if not OUT_type.startswith("int"): # NaN ne s'applique pas aux tableaux d'entiers
                try:
                    OUT_data[OUT_data == self.content.variables[IN_name]._FillValue] = numpy.nan
                except AttributeError:
                    pass
                    #print "[my_netcdf_file] %s: no _FillValue" % ( IN_name )
                
            # 3 - Multiplication by the scale factor
            try:
                OUT_data = OUT_data * self.content.variables[IN_name].scale_factor
            except AttributeError:
                pass
                #print('[my_netcdf_file] %s: no scale_factor') % ( IN_name )
            
            return OUT_data
            
        else:
            
            # 1 - Get data
            OUT_data = numpy.copy(group.variables[IN_name][:])
            
            # 2 - If data values are not "int", change _FillValue by NaN
            OUT_type = str(OUT_data.dtype)
            if not OUT_type.startswith("int"): # NaN ne s'applique pas aux tableaux d'entiers
                try:
                    OUT_data[OUT_data == group.variables[IN_name]._FillValue] = numpy.nan
                except AttributeError:
                    pass
                    #print "[my_netcdf_file] %s: no _FillValue" % ( IN_name )
                
            # 3 - Multiplication by the scale factor
            try:
                OUT_data = OUT_data * group.variables[IN_name].scale_factor
            except AttributeError:
                pass
                #print('[my_netcdf_file] %s: no scale_factor') % ( IN_name )
            
            return OUT_data            
    
    def getVarUnit(self, IN_name, group = None):
        '''
        Get the unit of variable named IN_name
        
        :param IN_name: name of the variable
        :type IN_name: string
        
        :return: unit
        :rtype: string
        '''
        if group == None:
            try:
                return self.content.variables[IN_name].units
            except AttributeError:
                return "-"
        else:
            try:
                return group.variables[IN_name].units
            except AttributeError:
                return "-"            
            
            
#######################################


class myNcWriter():

    def __init__(self, IN_filename):
        '''
        Constructor
        
        :param IN_filename: full path of the NetCDF file
        :type IN_filename: string
        
        Variables:
            filename / string: full path of the NetCDF file
            content / netCDF4.Dataset: content of the NetCDF file
        '''
        
        # Filename
        self.filename = IN_filename
        
        # Content
        self.content = Dataset(IN_filename, 'w')
        
    def close(self):
        '''
        Close the NetCDF file
        '''
        self.content.close()
    
    #----------------------------------------

    def add_group(self, IN_name):
        '''
        Add the data content of the variable
        
        :param IN_name: the name of the group
        :type IN_name: string
        '''
        return self.content.createGroup(IN_name)
        
        
    def add_dimension(self, IN_name, IN_size, group = None):
        '''
        Set the size of the dimension with the given name

        :param IN_name: the name of the dimension
        :type IN_name: string
        :param IN_size: the size of the dimension
        :type IN_size: int
        '''
        if group == None:
            self.content.createDimension(IN_name, IN_size)
        else:
            group.createDimension(IN_name, IN_size)  
            
    #----------------------------------------
        
    def add_global_attribute(self, IN_name, IN_value, group = None):
        '''
        Set the value of a global attribute with the given name

        :param IN_name: the name of the attribute
        :type IN_name: string
        :param IN_value: the value to store
        :type IN_value: unknown
        '''
        if group == None:
            self.content.setncattr(IN_name, IN_value)
        else:
            group.setncattr(IN_name, IN_value)   
            
    #----------------------------------------          
        
    def add_variable(self, IN_name, IN_datatype, IN_dimensions, IN_fill_value = None, IN_compress = True, group = None):
        '''
        Add the data content of the variable
        
        :param IN_name: the name of the variable
        :type IN_name: string
        :param IN_datatype: the type of the variable
        :type IN_datatype: ex np.int32, ...
        :param IN_dimensions: the name of the dimensions of the variable
        :type IN_dimensions: string
        :param IN_fill_value: the value used when no data is saved
        :type IN_fill_value: int or float
        :param IN_compress: true to compress the content of the variable (save disk space), else false
        :type IN_compress: boolean
        '''
        if group == None:
            self.content.createVariable(IN_name, IN_datatype, IN_dimensions, IN_compress, 4, fill_value = IN_fill_value)
        else:
            group.createVariable(IN_name, IN_datatype, IN_dimensions, IN_compress, 4, fill_value = IN_fill_value)
        
    def add_variable_attribute(self, IN_varname, IN_attname, IN_value, group = None):
        '''
        Set the value of a variable attribute with the given name

        :param IN_varname: the name of the variable
        :type IN_varname: string
        :param IN_attname: the name of the attribute of the variable
        :type IN_attname: string
        :param IN_value: the value to store
        :type IN_value: unknown
        '''
        if group == None:
            if IN_varname in self.content.variables:
                self.content.variables[IN_varname].setncattr(IN_attname, IN_value)
            else:
                exit("[my_netcdf_file/add_variable_attribute] Could not add variable attribute ; %s variable does not exist" % ( IN_varname ))
        
        else:
            if IN_varname in group.variables:
                group.variables[IN_varname].setncattr(IN_attname, IN_value)
            else:
                exit("[my_netcdf_file/add_variable_attribute] Could not add variable attribute ; %s variable does not exist" % ( IN_varname ))                
    
    #----------------------------------------
        
    def fill_variable(self, IN_name, IN_data, group = None):
        '''
        Write the given data into the named variable 

        :param IN_name: the name of the variable
        :type IN_name: string
        :param IN_data: the data to store
        :type IN_data: unknown
        '''
        if group == None: 
            if IN_name in self.content.variables:
                # Write the whole array
                self.content.variables[IN_name][:] = IN_data
            else:
                exit("[my_netcdf_file/fill_variable] Could not fill variable ; %s variable does not exist" % ( IN_name ))

        else: 
            if IN_name in group.variables:
                # Write the whole array
                group.variables[IN_name][:] = IN_data
            else:
                exit("[my_netcdf_file/fill_variable] Could not fill variable ; %s variable does not exist" % ( IN_name ))                


#######################################


if __name__ == "__main__":
    
    '''
    coord = myNcReader("test.gdem_orbit_cycle_0001_pass_0059.nc")
    print(coord.printListDim())
    #print coord.printListAtt()
    #print coord.printListVar()
    print("== Retrieve some data ==")
    lon = coord.getVarValue('longitude')
    print(lon[0])
    lat = coord.getVarValue('latitude')
    print(lat[0])
    alt = coord.getVarValue("altitude")
    print(alt[0])
    x = coord.getVarValue("x")
    print(x[0])
    y = coord.getVarValue("y")
    print(y[0])
    z = coord.getVarValue("z")
    print(z[0])
    '''
    
    data = myNcWriter("SWOT_L2_HR_PIXC_I_000_017_043N_L_20140101T133024_T_000.nc")
    