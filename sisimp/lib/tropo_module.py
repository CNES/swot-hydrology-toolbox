'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


import netCDF4 as nc
import numpy
import scipy.interpolate

import lib.height_model as height_model
import lib.my_api as my_api


class Tropo_module(object):

    def __init__(self, model):
        self.model = model
        
        
    def calculate_tropo_error_gaussian(self, az, r, stdv_error, mean_error, correlation_length):
        rmin = r.min()
        rmax = r.max()
        azmin = az.min()
        azmax = az.max()
        my_api.printInfo("Computing random tropo_error field on %d x %d image" % (azmax-azmin, rmax-rmin))

        tropo_map_rg_az = height_model.generate_2d_profile_gaussian(1, azmin, azmax+1, 1, rmin, rmax+1, stdv_error, lcorr = correlation_length)+mean_error
        tropo_error = tropo_map_rg_az[az-azmin,r-rmin]
        
        return tropo_error
        
        
    def calculate_tropo_error_map(self, latmin, az, r, map_file, correlation_length):
        rmin = r.min()
        rmax = r.max()
        azmin = az.min()
        azmax = az.max()
        my_api.printInfo("Computing map tropo_error field on %d x %d image" % (azmax-azmin, rmax-rmin))

        ds = nc.Dataset(map_file, 'r')
        delta_wtc_MEAN = ds.variables['delta_wtc_MEAN'][:]
        delta_wtc_STD = ds.variables['delta_wtc_STD'][:]
        latitude = ds.variables['latitude'][:]

        f=scipy.interpolate.interp1d(latitude, delta_wtc_MEAN,kind='linear')
        delta_wtc_MEAN_local=f(latmin)        
        f=scipy.interpolate.interp1d(latitude, delta_wtc_STD,kind='linear')
        delta_wtc_STD_local=f(latmin) 
               
        my_api.printInfo("%f cm mean biais and %f cm stv  estimated at latitude %f" % (delta_wtc_MEAN_local, delta_wtc_STD_local, latmin))

        tropo_map_rg_az = height_model.generate_2d_profile_gaussian(1, azmin, azmax+1, 1, rmin, rmax+1, delta_wtc_STD_local*0.01, lcorr = correlation_length)+delta_wtc_MEAN_local*0.01
        tropo_error = tropo_map_rg_az[az-azmin,r-rmin]
        
        return tropo_error
