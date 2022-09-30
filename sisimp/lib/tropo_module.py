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

    def __init__(self, tropo_model, rmin, rmax, azmin, azmax, tropo_error_stdv,
                 tropo_error_mean, tropo_error_correlation,
                 tropo_error_map_file, seed=None):
        self.model = tropo_model
        self.rmin = rmin
        self.rmax = rmax
        self.azmin = azmin
        self.azmax = azmax
        self.tropo_error_stdv = tropo_error_stdv
        self.tropo_error_mean = tropo_error_mean
        self.tropo_error_correlation = tropo_error_correlation
        self.tropo_error_map_file = tropo_error_map_file
        self.seed = seed

    def calculate_tropo_error_gaussian(self):
        my_api.printInfo("Computing random tropo_error field on %d x %d image"
                         % (self.azmax-self.azmin, self.rmax-self.rmin))
        self.tropo_map_rg_az = height_model.generate_2d_profile_gaussian(
            1, self.azmin, self.azmax+1, 1, self.rmin, self.rmax+1,
            self.tropo_error_stdv, lcorr=self.tropo_error_correlation,
            seed=self.seed) + self.tropo_error_mean

    def calculate_tropo_error_map(self, latmin):

        my_api.printInfo("Computing map tropo_error field on %d x %d image" % (self.azmax-self.azmin, self.rmax-self.rmin))

        ds = nc.Dataset(self.tropo_error_map_file, 'r')
        delta_wtc_MEAN = ds.variables['delta_wtc_MEAN'][:]
        delta_wtc_STD = ds.variables['delta_wtc_STD'][:]
        latitude = ds.variables['latitude'][:]

        f=scipy.interpolate.interp1d(latitude, delta_wtc_MEAN,kind='linear')
        delta_wtc_MEAN_local=f(latmin)        
        f=scipy.interpolate.interp1d(latitude, delta_wtc_STD,kind='linear')
        delta_wtc_STD_local=f(latmin) 
               
        my_api.printInfo("%f cm mean biais and %f cm stv  estimated at latitude %f" % (delta_wtc_MEAN_local, delta_wtc_STD_local, latmin))

        self.tropo_map_rg_az = height_model.generate_2d_profile_gaussian(
            1, self.azmin, self.azmax+1, 1, self.rmin, self.rmax+1,
            delta_wtc_STD_local*0.01, lcorr=self.tropo_error_correlation,
            seed=self.seed) + delta_wtc_MEAN_local*0.01

    def apply_tropo_error_on_pixels(self, az, r):
        self.tropo_2d_field = self.tropo_map_rg_az[az-self.azmin,r-self.rmin]
        

    def generate_tropo_field_over_pass(self, latmin):
        if self.model == 'gaussian':
            self.calculate_tropo_error_gaussian()
        elif self.model == 'map':
            self.calculate_tropo_error_map(latmin)
        else:
            self.tropo_map_rg_az = None
        
        
