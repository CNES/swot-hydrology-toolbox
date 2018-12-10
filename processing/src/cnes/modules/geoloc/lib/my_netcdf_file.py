# -*- coding: utf8 -*-
'''
.. module my_netcdf_file.py
    :synopsis: Deals with NetCDF files (reader and writer)
    Created on 08/04/2016

.. module author: Claire POTTIER - CNES DCT/SI/TR

Copyright (c) 2016 CNES. All rights reserved.
'''

from netCDF4 import Dataset
import numpy
import xarray as xr



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
    
    def getListDim(self):
        '''
        Get the dictionnary of the dimensions

        :return: dictionnary of dimensions
        :rtype: dict
        '''
        return self.content.dimensions
        
    def printListDim(self):
        '''
        Print list of dimensions
        '''
        
        print("== DIMENSIONS ==")
        
        for value in self.getListDim():
            print(value + " - size = " + str(self.getDimValue(value)))
    
    def getDimValue(self, IN_name):
        '''
        Get the size of the dimension named IN_name
        
        :param IN_name: name of the dimension
        :type IN_name: string
        
        :return: size of the dimension
        :rtype: int
        '''
        return len(self.content.dimensions[IN_name])
    
    #----------------------------------------
    
    def getListAtt(self):
        '''
        Get the list of all the global attributes included in the NetCDF file
        
        :return: global attribute names
        :rtype: list of string
        '''
        return self.content.ncattrs()
        
    def printListAtt(self):
        '''
        Print list of global attributes
        '''
        
        print("== GLOBAL ATTRIBUTES ==")
        for value in self.getListAtt():
            print(value + " = " + self.getAttValue(value))
    
    def getAttValue(self, IN_name):
        '''
        Get the value associated to the global attribute named IN_name
        
        :param IN_name: name of the attribute
        :type IN_name: string
        
        :return: value associated to the attribute
        :rtype: string
        '''
        return self.content.getncattr(IN_name)
    
    #----------------------------------------
        
    def printListVar(self):
        '''
        Print list of variables
        '''
        
        list_var = sorted(self.content.variables.keys())
        
        print("== VARIABLES ==")
        for value in list_var:
            print(value + " - units = " + self.getVarUnit(value))
    
    def getVarValue(self, IN_name):
        '''
        Get the data associated to the variable IN_name
        _FillValue value are filled by NaN and the multiplication by the scale_factor is done if needed
        
        :param IN_name: name of the variable
        :type IN_name: string
        
        :return: formatted data
        :rtype: [depend on the input variable]
        '''
        
        # 1 - Get data
        OUT_data = numpy.copy(self.content.variables[IN_name][:])
        
        # 2 - If data values are not "int", change _FillValue by NaN
        OUT_type = str(OUT_data.dtype)
        if not OUT_type.startswith("int"): # NaN ne s'applique pas aux tableaux d'entiers
            try:
                OUT_data[OUT_data == self.content.variables[IN_name]._FillValue] = numpy.nan
            except AttributeError:
                print("[my_netcdf_file] %s: no _FillValue" % ( IN_name ))
            
        # 3 - Multiplication by the scale factor
        try:
            OUT_data = OUT_data * self.content.variables[IN_name].scale_factor
        except AttributeError:
            print(('[my_netcdf_file] %s: no scale_factor') % ( IN_name ))
        
        return OUT_data
    
    def getVarUnit(self, IN_name):
        '''
        Get the unit of variable named IN_name
        
        :param IN_name: name of the variable
        :type IN_name: string
        
        :return: unit
        :rtype: string
        '''
        
        try:
            return self.content.variables[IN_name].units
        except AttributeError:
            return "-"

    def getVarValue2d(self, keyword_1, keyword_2):
    
        df = self.content.to_dataframe()
        df = df.reset_index()

        nb_lon = self.content.dims[keyword_1]
        nb_lat = self.content.dims[keyword_2]
        
        return df, nb_lon, nb_lat  
              
#######################################

class myNcWriter():

    def __init__(self, IN_filename):
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
        self.content = Dataset(IN_filename, 'w')
        
    def close(self):
        '''
        Close the NetCDF file
        '''
        self.content.close()
    
    #----------------------------------------

    def add_dimension(self, IN_name, IN_size):
        '''
        Set the size of the dimension with the given name

        :param IN_name: the name of the dimension
        :type IN_name: string
        :param IN_size: the size of the dimension
        :type IN_size: int
        '''
        self.content.createDimension(IN_name, IN_size)
    
    #----------------------------------------
        
    def add_global_attribute(self, IN_name, IN_value):
        '''
        Set the value of a global attribute with the given name

        :param IN_name: the name of the attribute
        :type IN_name: string
        :param IN_value: the value to store
        :type IN_value: unknown
        '''
        self.content.setncattr(IN_name, IN_value)
    
    #----------------------------------------

    def add_variable(self, IN_name, IN_datatype, IN_dimensions, IN_fill_value = None, IN_compress = True):
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
        self.content.createVariable(IN_name, IN_datatype, IN_dimensions, IN_compress, 4, fill_value = IN_fill_value)
        
    def add_variable_attribute(self, IN_varname, IN_attname, IN_value):
        '''
        Set the value of a variable attribute with the given name

        :param IN_varname: the name of the variable
        :type IN_varname: string
        :param IN_attname: the name of the attribute of the variable
        :type IN_attname: string
        :param IN_value: the value to store
        :type IN_value: unknown
        '''
        
        if IN_varname in self.content.variables:
            self.content.variables[IN_varname].setncattr(IN_attname, IN_value)
        else:
            exit("[my_netcdf_file/add_variable_attribute] Could not add variable attribute ; %s variable does not exist" % ( IN_varname ))
    
    #----------------------------------------
        
    def fill_variable(self, IN_name, IN_data):
        '''
        Write the given data into the named variable 

        :param IN_name: the name of the variable
        :type IN_name: string
        :param IN_data: the data to store
        :type IN_data: unknown
        '''
        
        if IN_name in self.content.variables:
            # Write the whole array
            self.content.variables[IN_name][:] = IN_data
        else:
            exit("[my_netcdf_file/fill_variable] Could not fill variable ; %s variable does not exist" % ( IN_name ))
            
#######################################

if __name__ == "__main__":
    
    coord = myNcReader("D:\\Utilisateurs\\pottierc\\Documents\\workspace_eclipse\\swot_py\\res\\test.gdem_orbit_cycle_0001_pass_0059.nc")
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
