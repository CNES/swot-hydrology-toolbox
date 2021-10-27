# -*- coding: utf8 -*-
'''
Reader for netcdf file
Convert netcdf files to geopandas dataframe
Pixel cloud reader using geopandas

Copyright (c) 2018 CNES. All rights reserved.
'''

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import xarray as xr
import math
import numpy as np
from netCDF4 import Dataset
import osr
import textwrap

EARTH_RADIUS = 6371000.
class PixcReader():

    def __init__(self, pixc: str, vec: str = None):
        '''
        Read pixel cloud file and vec file.
        
        :param pixc: filename of the pixel cloud file
        :param vec: filename of the vec file
        '''
        
        # Read pixel cloud and convert to dataframe
        dnc = xr.open_dataset(pixc,group="pixel_cloud")  
        df = dnc.to_dataframe()
        df = df.loc[df.index.get_level_values('complex_depth') == 0]

        
        try:        
            df = df[['classification','pixel_area',
                     'longitude','latitude', 'height', 
                     'range_index', 'azimuth_index',
                     'inc', 'cross_track', 'pole_tide', 'load_tide_got', 'load_tide_fes', 'solid_earth_tide', 'geoid']]
                     
            tide_correctly_loaded = True
        except:
            # should be removed, kept to be compatible with old simulations...
            df = df[['classification','pixel_area',
                'longitude','latitude', 'height', 
                'range_index', 'azimuth_index',
                'inc', 'cross_track', 'geoid']]
            tide_correctly_loaded = False
         
        df = df.loc[(df.classification > 1.0)] # Keep only water pixel point
        
        
        # Read vec file and convert to dataframe
        if vec:
            vec_dnc = xr.open_dataset(vec,decode_times=False)  
            vec_df = vec_dnc.to_dataframe()
            vec_df[['index_i']] = vec_df[['range_index']].astype(int)
            vec_df[['index_j']] = vec_df[['azimuth_index']].astype(int)
            vec_df = vec_df[['index_i', 'index_j', 'longitude_vectorproc',
                             'latitude_vectorproc', 'height_vectorproc']]
            vec_df.rename(columns={'longitude_vectorproc': 'longitude', 
                    'latitude_vectorproc': 'latitude',
                    'height_vectorproc': 'height'}, inplace=True)
            vec_df = vec_df.set_index(['index_i', 'index_j'])
            
            df = df.drop(['longitude','latitude', 'height'],axis=1)
            df[['index_i']] = df[['range_index']].astype(int)
            df[['index_j']] = df[['azimuth_index']].astype(int)
            df = df.set_index(['index_i', 'index_j'])
            
            df = pd.concat([df, vec_df], axis=1,join='inner')
            df = df.reset_index()
            df = df.drop(['index_i','index_j'],axis=1)
        
        geom = [Point(x,y) for x, y in zip(df['longitude'], df['latitude'])]

        # Create geodataframe
        self.data = gpd.GeoDataFrame(df, geometry=geom)
        # Drop Nan values
        self.data.dropna(inplace=True)
        # Longitude must be set between -180 and 180
        self.data['longitude'] = self.data.apply(lambda row: row.longitude if row.longitude < 180.0 else row.longitude - 360.0,
                                                 axis=1)
        self.data[['range_index']] = self.data[['range_index']].astype(int)
        self.data[['azimuth_index']] = self.data[['azimuth_index']].astype(int)
         
        # Add information (range spacing)
        attrs_globals = xr.open_dataset(pixc).attrs
        nominal_range_spacing = attrs_globals["nominal_slant_range_spacing"]
        near_range = attrs_globals["near_range"]
        try:
            height = xr.open_dataset(pixc,group="tvp")["altitude"]
        except:
            # should be removed, kept to be compatible with old simulations...
            height = xr.open_dataset(pixc,group="tvp")["height"]

        #to be cleaned
      
        self.data['cos_alpha'] = self.data.apply(lambda row:
                                           float(height[row['azimuth_index']])/(near_range + row['range_index'] * nominal_range_spacing), 
                                            axis=1)
        self.data['cos_alpha']  = self.data.apply(lambda row: min(row['cos_alpha'], 0.99999),axis=1)
        self.data['range_spacing'] = self.data.apply(lambda row: nominal_range_spacing/math.sin(math.acos(row['cos_alpha'])),axis=1)
 
                                                   
        # ~ self.data['range_spacing'] = self.data.apply(lambda row:
                                           # ~ nominal_range_spacing / 
                                           # ~ math.sin(math.acos(
                                           # ~ float(height[row['azimuth_index']])/(near_range + row['range_index'] * nominal_range_spacing))), 
                                           # ~ axis=1)
        
        
        # Compute wse using heights and corrections 
        ## TODO : Deal with nan correction and consider them as 0.
        if tide_correctly_loaded:
            self.data['elevation'] = self.data.apply(lambda row: row['height']-row['pole_tide']-
                                                row['load_tide_got']-row['load_tide_fes']-
                                                row['solid_earth_tide']-row['geoid'], axis=1)
        else:
            # should be removed, kept to be compatible with old simulations...
            self.data['elevation'] = self.data.apply(lambda row: row['height']-row['geoid'], axis=1)            
  
        # Get attributes
        self.range_max = dnc.attrs['interferogram_size_range']
        self.azimuth_max = dnc.attrs['interferogram_size_azimuth']
        
        # Read orbit and convert to dataframe
        dnc_trj = xr.open_dataset(pixc,group="tvp")  
        df_trj = dnc_trj.to_dataframe()    
            
        d_az_trj =  np.sqrt((df_trj[['x']].values[0]-df_trj[['x']].values[1])**2+(df_trj[['y']].values[0]-df_trj[['y']].values[1])**2+(df_trj[['z']].values[0]-df_trj[['z']].values[1])**2)
        alt = np.sqrt(df_trj[['x']].values[0]**2+df_trj[['y']].values[0]**2+df_trj[['z']].values[0]**2)
        d_az_ground = d_az_trj*EARTH_RADIUS/alt
        
        self.along_track_sampling = d_az_ground
        
    def get_data(self) -> gpd.GeoDataFrame:
        return self.data.copy()


def textjoin(text):
    """ Dedent join and strip text """
    text = textwrap.dedent(text)
    text = text.replace('\n', ' ')
    text = text.strip()
    return text
    
def from_file(filename: str) -> gpd.GeoDataFrame:
        '''
        Read netcdf file and convert to geodataframe.
        Require latitie and longitude columns
        
        :param filename: netcdf filename
        :return geodataframe
        '''

        # Read pixel cloud and convert to dataframe
        dnc = xr.open_dataset(filename)
        df = dnc.to_dataframe()
        # Geometry
        geom = [Point(x,y) for x, y in zip(df['longitude'], df['latitude'])]
        # Create geodataframe
        return gpd.GeoDataFrame(df, geometry=geom)

def write_raster_gridded(filename: str, mode: str, x: np.array, y: np.array, out_image: np.array, out_dist_min_2d: np.array, \
                out_dist_mean_2d: np.array, qual_flag: np.array, resolution: float, espg: str):
        '''
        write output raster netcdf file 
        '''

        ds = Dataset(filename, 'w')
        ds.Conventions = "CF-1.7" 
        ds.title = "Level 2 KaRIn High Rate FPDEM Gridded Data Product"
        ds.institution = "CNES"
        ds.source = "Large scale Simulator"
        ds.history = "None"
        ds.mission_name ="SWOT"
        ds.references = "None"
        ds.reference_document = "None"
        ds.contact = "damien.desroches@cnes.fr"
        ds.coordinate_reference_system = espg
        ds.sampling = resolution
        ds.short_name = "L2_HR_FPDEM_Gridded"
        ds.descriptor_string = filename.split('_')[5]
        ds.crid = "Dx0000"
        ds.product_version = "1"
        ds.pge_name = "Pge_v0"
        ds.pge_version = "1"
        ds.time_coverage_start = filename.split('_')[6]
        ds.time_coverage_end = filename.split('_')[7]
        ds.geospatial_lon_min = np.min(y)
        ds.geospatial_lon_max = np.max(y)
        ds.geospatial_lat_min = np.min(x)
        ds.geospatial_lat_max = np.max(x)
        
        

        
        # ~ if mode == 'utm':
            # ~ x_dim    = ds.createDimension('x',len(x))
            # ~ y_dim    = ds.createDimension('y',len(y))
            # ~ x_var    = ds.createVariable("x", "float64", ("x"), fill_value=-9999.)
            # ~ y_var    = ds.createVariable("y", "float64", ("y"), fill_value=-9999.)
            # ~ z_var    = ds.createVariable("elevation", "float64", ("y", "x"), fill_value=-9999.)   
            # ~ z_var_u  = ds.createVariable("elevation_uncert", "float64", ("y", "x"), fill_value=-9999.)   
            # ~ min_dist_var = ds.createVariable("distance_to_closest", "float64", ("y", "x"), fill_value=-9999.)   
            # ~ mean_dist_var = ds.createVariable("mean_distance", "float64", ("y", "x"), fill_value=-9999.)   
            # ~ qual_var = ds.createVariable("fpdem_qual", "float64", ("y", "x"), fill_value=-9999.)   

            # ~ ds.projection = "UTM Coordinates"
                        
        if mode == 'latlon':
            x_dim = ds.createDimension('latitude',len(x))
            y_dim = ds.createDimension('longitude',len(y))

            
            coordinate_system = osr.SpatialReference()
            coordinate_system.ImportFromEPSG(4326)
    
            crs = ds.createVariable("crs", 'S1')
            crs.long_name = 'CRS Definition'
            crs.grid_mapping_name = 'latitude_longitude'
            crs.geographic_crs_name = 'WGS 84'
            crs.reference_ellipsoid_name = 'WGS 84'
            crs.horizontal_datum_name = 'WGS_1984'
            crs.prime_meridian_name = 'Greenwich'
            crs.longitude_of_prime_meridian = 0.
            crs.semi_major_axis = 6378137.
            crs.inverse_flattening = 298.257223563
            crs.crs_wkt = coordinate_system.ExportToWkt()
            crs.spatial_ref = coordinate_system.ExportToWkt()
            crs.comment = 'Geodetic lat/lon coordinate reference system.'
            
                    
            x_var = ds.createVariable("latitude", "float64", ("latitude"), fill_value=-9999.)
            x_var.long_name = 'latitude (positive N, negative S)'
            x_var.standard_name = 'latitude'
            x_var.units = 'degrees_north'
            x_var.valid_min = -80
            x_var.valid_max = 80
            x_var.comment = textjoin("""
                    Geodetic latitude [-80,80] (degrees north of equator) of
                    the pixel.""")
                    
            y_var = ds.createVariable("longitude", "float64", ("longitude"), fill_value=-9999.)
            y_var.long_name = 'longitude (degrees East)'
            y_var.standard_name = 'longitude'
            y_var.units = 'degrees_east'
            y_var.valid_min = -180
            y_var.valid_max = 180
            y_var.comment = textjoin("""
                    Geodetic longitude [-180,180] (east of the Greenwich meridian) of
                    the pixel.""")            
            
            z_var    = ds.createVariable("elevation", "float64", ("longitude", "latitude"), fill_value=-9999.)  
            z_var.long_name = 'elevation (meters)'
            z_var.standard_name = 'elevation'
            z_var.units = 'meters'
            z_var.valid_min = -9999
            z_var.valid_max = 9999
            z_var.comment = textjoin("""
                    Elevation [-9999,9999] (relative  to the geoid) of
                    the pixel.""")             
            
            z_var_u  = ds.createVariable("elevation_uncert", "float64", ("longitude", "latitude"), fill_value=-9999.)   
            z_var_u.long_name = 'elevation uncertainty (meters)'
            z_var_u.standard_name = 'elevation uncertainty'
            z_var_u.units = 'meters'
            z_var_u.valid_min = -9999
            z_var_u.valid_max = 9999
            z_var_u.comment = textjoin("""
                    Elevation uncertainty [-9999,9999] (relative  to the geoid) of
                    the pixel.""") 
                    
            min_dist_var = ds.createVariable("distance_to_closest", "float64", ("longitude", "latitude"), fill_value=-9999.)   
            min_dist_var.long_name = 'minimum distance to closest pixel (meters)'
            min_dist_var.standard_name = 'distance to closest'
            min_dist_var.units = 'meters'
            min_dist_var.valid_min = 0
            min_dist_var.valid_max = 9999
            min_dist_var.comment = textjoin("""
                    Distance to closest pixel [0,9999] from each
                    the pixel.""") 
                                
            mean_dist_var = ds.createVariable("mean_distance", "float64", ("longitude", "latitude"), fill_value=-9999.)   
            mean_dist_var.long_name = 'mean distance from pixel to aggregated pixels (meters)'
            mean_dist_var.standard_name = 'mean distance'
            mean_dist_var.units = 'meters'
            mean_dist_var.valid_min = 0
            mean_dist_var.valid_max = 9999
            mean_dist_var.comment = textjoin("""
                    Mean distance from pixels to aggregated pixels [0,9999].""") 
                    
            qual_var = ds.createVariable("fpdem_qual", "float64", ("longitude", "latitude"), fill_value=9999.)   
            qual_var.long_name = 'Quality flag for each pixel'
            qual_var.standard_name = 'Quality flag'
            qual_var.units = 'None'
            qual_var.valid_min = 0
            qual_var.valid_max = 10
            qual_var.comment = textjoin("""
                    Quality flag [0,10] of
                    the pixel.""") 
                      
             
        x_var[:] = x
        y_var[:] = y
        z_var[:,:] = out_image  
        min_dist_var[:,:] = out_dist_min_2d  
        mean_dist_var[:,:] = out_dist_mean_2d  
        qual_var[:,:] = qual_flag  
         

        ds.close

def write_raster_ungridded(data, filename):
    
    ds = Dataset(filename, 'w')
    ds.Conventions = "CF-1.7" 
    ds.title = "Level 2 KaRIn High Rate FPDEM Ungridded Data Product"
    ds.institution = "CNES"
    ds.source = "Large scale Simulator"
    ds.history = "None"
    ds.mission_name ="SWOT"
    ds.references = "None"
    ds.reference_document = "None"
    ds.contact = "damien.desroches@cnes.fr"
    ds.coordinate_reference_system = data.attrs['espg']
    ds.short_name = "L2_HR_FPDEM_Ungridded"
    ds.descriptor_string = filename.split('_')[5]
    ds.crid = "Dx0000"
    ds.product_version = "1"
    ds.pge_name = "Pge_v0"
    ds.pge_version = "1"
    ds.espg = data.attrs["espg"]
    ds.time_coverage_start = filename.split('_')[6]
    ds.time_coverage_end = filename.split('_')[7]
    ds.geospatial_lon_min = np.min(data.variables["longitude"])
    ds.geospatial_lon_max = np.max(data.variables["longitude"])
    ds.geospatial_lat_min = np.min(data.variables["latitude"])
    ds.geospatial_lat_max = np.max(data.variables["latitude"])

    index_dim = ds.createDimension('index',len(data.variables['longitude']))
    
    coordinate_system = osr.SpatialReference()
    coordinate_system.ImportFromEPSG(4326)

    
    x_var = ds.createVariable("latitude", "float64", ("index"), fill_value=-9999.)
    x_var.long_name = 'latitude (positive N, negative S)'
    x_var.standard_name = 'latitude'
    x_var.units = 'degrees_north'
    x_var.valid_min = -80
    x_var.valid_max = 80
    x_var.comment = textjoin("""
            Geodetic latitude [-80,80] (degrees north of equator) of
            the pixel.""")
            
    y_var = ds.createVariable("longitude", "float64", ("index"), fill_value=-9999.)
    y_var.long_name = 'longitude (degrees East)'
    y_var.standard_name = 'longitude'
    y_var.units = 'degrees_east'
    y_var.valid_min = -180
    y_var.valid_max = 180
    y_var.comment = textjoin("""
            Geodetic longitude [-180,180] (east of the Greenwich meridian) of
            the pixel.""")            
                    
    z_var    = ds.createVariable("elevation", "float64", ("index"), fill_value=-9999.)  
    z_var.long_name = 'elevation (meters)'
    z_var.standard_name = 'elevation'
    z_var.units = 'meters'
    z_var.valid_min = -9999
    z_var.valid_max = 9999
    z_var.comment = textjoin("""
            Elevation [-9999,9999] (relative  to the geoid) of
            the pixel.""")      

    data_valid_var  = ds.createVariable("data_valid", "int64", ("index"), fill_value=-9999)  
    data_valid_var.long_name = 'data validation flag'
    data_valid_var.standard_name = 'data valid'
    data_valid_var.units = ''
    data_valid_var.valid_min = -9999
    data_valid_var.valid_max = 9999
    data_valid_var.comment = textjoin("""
            Data validation flag [0, 1] of
            the pixel.""")  
                        
    x_var[:] = data.variables['latitude']
    y_var[:] = data.variables['longitude']
    z_var[:] = data.variables['elevation'] 
    data_valid_var[:] = data.variables["data_valid"]
