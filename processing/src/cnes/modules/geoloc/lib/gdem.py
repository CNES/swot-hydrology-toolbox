'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


import cnes.modules.geoloc.lib.my_netcdf_file as myCdf
import cnes.modules.geoloc.lib.tools as lib

import numpy as np
from scipy.interpolate import griddata



class Gdem(object):
    
    def __init__(self, IN_gdem_file, list_landtype_flags):
      
        gdem_main = myCdf.myNcReader(IN_gdem_file, dim=2)
        df, nb_lon, nb_lat = gdem_main.getVarValue2d('longitude', 'latitude')
        self.dataframe = df
        self.nb_lon = nb_lon
        self.nb_lat = nb_lat
        
        print("Calcul des longitudes et latitudes du gdem")
                
        self.elevation = self.dataframe.elevation.values.reshape(self.nb_lat,  self.nb_lon)
        self.landtype = self.dataframe.landtype.values.reshape(self.nb_lat,  self.nb_lon)

        self.longitude = self.dataframe.longitude.values
        self.latitude = self.dataframe.latitude.values
                
    def compute_input_gdem_dem(self, pixc, dem_res_x, dem_res_y):
        '''
        Compute an extract of the gdem with dem_res_x, dem_res_y resolution on the same area of a pixel cloud pixc
        '''

        
        print("Determination de la zone d interet sur le gdem")
   
        self.nb_lon_crop = int(((self.dataframe.loc[(self.dataframe.longitude > pixc.lonmin) & (self.dataframe.longitude < pixc.lonmax)]).elevation.size)/self.nb_lat)
        self.nb_lat_crop = int(((self.dataframe.loc[(self.dataframe.latitude > pixc.latmin) & (self.dataframe.latitude < pixc.latmax)]).elevation.size)/self.nb_lon)
        
        self.crop_dataframe = self.dataframe.loc[(self.dataframe.longitude > pixc.lonmin) & (self.dataframe.longitude < pixc.lonmax) & (self.dataframe.latitude > pixc.latmin) & (self.dataframe.latitude < pixc.latmax) ]
    
        self.crop_elevation = self.crop_dataframe.elevation.values.reshape(self.nb_lat_crop , self.nb_lon_crop)
        self.crop_landtype = self.crop_dataframe.landtype.values.reshape(self.nb_lat_crop , self.nb_lon_crop)

        self.crop_longitude = self.crop_dataframe.longitude.values
        self.crop_latitude = self.crop_dataframe.latitude.values

        print("Conversion et reechantillonnage du GDEM en xyz")

        
        xyz = lib.llh2xyz(self.crop_longitude.reshape(self.nb_lat_crop , self.nb_lon_crop), self.crop_latitude.reshape(self.nb_lat_crop , self.nb_lon_crop), self.crop_elevation, IN_flag_rad=False)
        
        size = xyz[0].size
        values = np.zeros([size,4])
        values[:,0] = np.reshape(xyz[0],size)
        values[:,1] = np.reshape(xyz[1],size)
        values[:,2] = np.reshape(self.crop_elevation,size)
        values[:,3] = np.reshape(self.crop_landtype,size)

        values[:,0] = np.where(np.isnan(values[:,0]), 0, values[:,0])
        values[:,1] = np.where(np.isnan(values[:,1]), 0, values[:,1])

        grid_x, grid_y = np.mgrid[pixc.xmin:pixc.xmax:dem_res_x,pixc.ymin:pixc.ymax:dem_res_y]
        self.crop_dem_xyz = griddata(values[::100,0:2], values[::100,2],(grid_x, grid_y), method='linear')
        self.crop_landtype_xyz = griddata(values[::100,0:2], values[::100,3],(grid_x, grid_y), method='linear')
