# -*- coding: utf8 -*-
'''
Create raster from floodplain dem pixel cloud

Copyright (c) 2018, CNES
'''
import os
import numpy as np
import xarray as xr
import lib.my_rdf_file as my_rdf
import lib.idw as idw
import lib.rasterization as rasterization
import floodplain.io.shp as shp

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
    
import pyproj
import shapely.geometry as shpgeo
import shapefile
  
import rasterio.mask
from rasterio.io import MemoryFile
from rasterio.transform import from_origin, from_bounds
from floodplain.io.nc import write_raster_gridded
from lib.constants import EARTH_RADIUS

from scipy.spatial import cKDTree
from scipy import interpolate


class FPDEM_Raster(object):
    """
    Class FPDEM_Raster
    Main class to compute HR_FPDEM prodcut
    """
    def __init__(self, param, input_file = None, output_file = None, mask = None):
        """
        Constructor: initialize variables
        
        :param in_params: input parameters to run the processor
        :type in_params: dictionary
        """
        
        if input_file == None:
            self.input_file = param.getValue("input_file")
        else:
            self.input_file = input_file 

        if input_file == None:
            self.output_file = param.getValue("output_file")
        else:
            self.output_file = output_file 

        if mask == None:
            self.mask = param.getValue("mask_file")
        else:
            self.mask = mask 
            

        self.method_raster = param.getValue("method_raster")
        
        if self.method_raster == 'cars':
            self.resolution = float(param.getValue("resolution"))
            self.sigma = float(param.getValue("sigma"))
            self.radius = float(param.getValue("radius"))
            print('This method is currently not operational')
            print("Please select idw raster_method")

        if self.method_raster == 'idw':
            self.resolution = float(param.getValue("resolution"))
            self.number_of_neighbors_considered = int(param.getValue("number_of_neighbors_considered"))
            self.mode = param.getValue("mode")

        if self.method_raster == 'exotic':
            self.resolution = float(param.getValue("resolution"))
            self.mode = param.getValue("mode")
                    
        self.plot = param.getValue("plot")
        
    def compute_fpdem_raster(self):
        if self.method_raster == 'idw':
            self.compute_raster_from_pixc_idw()
        if self.method_raster == 'exotic':
            self.compute_raster_exotic()
            
        else:
            print('This method is currently not operational')
            print("Please select idw raster_method")
        

    def compute_raster_from_pixc_cars(self, input_file=None, output_file=None, mask=None):

        cloud_xr = xr.open_dataset(self.input_file)
        cloud_df = cloud_xr.to_dataframe()

        xstart, ystart, xsize, ysize = \
            rasterization.compute_xy_starts_and_sizes(self.resolution, cloud_df)
        raster = rasterization.\
            rasterize(
                cloud_df, self.resolution, cloud_xr.attrs['espg'], xstart, ystart, xsize, ysize, self.sigma, self.radius,
                hgt_no_data=np.nan, color_no_data=np.nan)
                   
        raster.to_netcdf(self.output_file)
    
    def load_input_fpdem_raster(self):
        cloud_xr = xr.open_dataset(self.input_file)
        self.cloud_df_raster = cloud_xr.to_dataframe()
        self.espg = cloud_xr.attrs['espg']
    
    def compute_raster_from_pixc_idw(self):
                        
        # Extract xyz
        if self.mode == 'utm':
            nb = len(self.cloud_df_raster['x'])
            xyz_data = np.ndarray((nb,3))
            xyz_data[:,0] = self.cloud_df_raster['x']
            xyz_data[:,1] = self.cloud_df_raster['y']

        # Extract latlon
        if self.mode == 'latlon':
            nb = len(self.cloud_df_raster['longitude'])
            lonlat_data = np.ndarray((nb,3))
            lonlat_data[:,0] = self.cloud_df_raster['longitude']
            lonlat_data[:,1] = self.cloud_df_raster['latitude']
        
        # Train
        if self.mode == 'utm':
            idw_tree = idw.tree(xyz_data[:,0:2], self.cloud_df_raster['elevation'])
        if self.mode == 'latlon':
            idw_tree = idw.tree(lonlat_data[:,0:2], self.cloud_df_raster['elevation'])
            
        # Compute grid
        if self.mode == 'utm':
            # Interpolation
            x0 = min(xyz_data[:,0])
            y0 = min(xyz_data[:,1])
            nx = abs(int((x0-max(xyz_data[:,0]))/self.resolution))+1
            ny = abs(int((y0-max(xyz_data[:,1]))/self.resolution))+1
        if self.mode =='latlon':
            x0 = min(lonlat_data[:,0])
            y0 = min(lonlat_data[:,1])
            nx = abs(int((x0-max(lonlat_data[:,0]))/self.resolution))+1
            ny = abs(int((y0-max(lonlat_data[:,1]))/self.resolution))+1 
                
        x1 = x0 + self.resolution*(nx-1)
        y1 = y0 + self.resolution*(ny-1)
        
        x = np.linspace(x0, x1, nx)
        y = np.linspace(y0, y1, ny)
        
        GRID = np.meshgrid(x, y)
        grid_shape = GRID[0].shape
        GRID = np.reshape(GRID, (2, -1)).T    
        
        # Compute idw interpolation
        Z, dist_min, dist_mean, z_rel = idw_tree(GRID, k=self.number_of_neighbors_considered)

        if self.mode == 'latlon':
            dist_min = 2*np.pi*EARTH_RADIUS*dist_min/360.
            dist_mean = 2*np.pi*EARTH_RADIUS*dist_mean/360.
            
        mask = shapefile.Reader(self.mask)
            
        proj = pyproj.Proj(init='EPSG:'+str(self.espg))    

        # Reproject into 2d grid
        Z2d = Z.reshape(grid_shape)
        dist_min_2d = dist_min.reshape(grid_shape)
        dist_mean_2d = dist_mean.reshape(grid_shape)
        z_rel_2d = z_rel.reshape(grid_shape)
        
        # Compute quality flag
        qual_flag_2d = compute_qual_flag(dist_mean_2d, z_rel_2d)
        
        # Compute mask in latlon and utm coordinates
        extracted_area_utm = []     
        extracted_area_latlon = []
        
        for extracted_area in mask.shapeRecords():
            extracted_area_utm.append(toFromUTM(extracted_area.shape, proj))
            extracted_area_latlon.append(extracted_area.shape)
            
        transform = rasterio.transform.from_origin(min(x), min(y), self.resolution, -self.resolution)
        
        #Filter data out of mask
        if self.mode == 'utm':
            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=Z2d.shape[0], width=Z2d.shape[1], count=1, dtype=Z2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(Z2d, 1)
                    out_image, out_transform = rasterio.mask.mask(src, extracted_area_utm, crop=False, indexes = 1, nodata=np.nan)

            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=dist_min_2d.shape[0], width=dist_min_2d.shape[1], count=1, dtype=dist_min_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(dist_min_2d, 1)
                    out_dist_min_2d, out_transform = rasterio.mask.mask(src, extracted_area_utm, crop=False, indexes = 1, nodata=np.nan)

            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=dist_mean_2d.shape[0], width=dist_mean_2d.shape[1], count=1, dtype=dist_mean_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(dist_mean_2d, 1)
                    out_dist_mean_2d, out_transform = rasterio.mask.mask(src, extracted_area_utm, crop=False, indexes = 1, nodata=np.nan)

            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=z_rel_2d.shape[0], width=z_rel_2d.shape[1], count=1, dtype=z_rel_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(z_rel_2d, 1)
                    out_z_rel_2d, out_transform = rasterio.mask.mask(src, extracted_area_utm, crop=False, indexes = 1, nodata=np.nan)
                    
            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=qual_flag_2d.shape[0], width=qual_flag_2d.shape[1], count=1, dtype=qual_flag_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(qual_flag_2d, 1)
                    out_qual_flag_2d, out_transform = rasterio.mask.mask(src, extracted_area_utm, crop=False, indexes = 1, nodata=np.nan)
                    
        if self.mode == 'latlon':
            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=Z2d.shape[0], width=Z2d.shape[1], count=1, dtype=Z2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(Z2d, 1)
                    out_image, out_transform = rasterio.mask.mask(src, extracted_area_latlon, crop=False, indexes = 1, nodata=np.nan)        

            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=dist_min_2d.shape[0], width=dist_min_2d.shape[1], count=1, dtype=dist_min_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(dist_min_2d, 1)
                    out_dist_min_2d, out_transform = rasterio.mask.mask(src, extracted_area_latlon, crop=False, indexes = 1, nodata=np.nan)

            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=dist_mean_2d.shape[0], width=dist_mean_2d.shape[1], count=1, dtype=dist_mean_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(dist_mean_2d, 1)
                    out_dist_mean_2d, out_transform = rasterio.mask.mask(src, extracted_area_latlon, crop=False, indexes = 1, nodata=np.nan)

            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=z_rel_2d.shape[0], width=z_rel_2d.shape[1], count=1, dtype=z_rel_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(z_rel_2d, 1)
                    out_z_rel_2d, out_transform = rasterio.mask.mask(src, extracted_area_latlon, crop=False, indexes = 1, nodata=np.nan)
                                    
            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=qual_flag_2d.shape[0], width=qual_flag_2d.shape[1], count=1, dtype=qual_flag_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(qual_flag_2d, 1)
                    out_qual_flag_2d, out_transform = rasterio.mask.mask(src, extracted_area_latlon, crop=False, indexes = 1, nodata=np.nan)

            
        self.x = x
        self.y = y
        self.out_image = out_image
        self.out_dist_min_2d = out_dist_min_2d
        self.out_dist_mean_2d = out_dist_mean_2d
        self.out_qual_flag_2d = out_qual_flag_2d
        
        print("Output raster shape = ", Z2d.shape)
        
        if self.plot =='yes':
            
            fig, [[ax1, ax2],[ax3, ax4]] = plt.subplots(2,2, sharex=True, sharey=True, figsize=(18,18))
            if self.mode == 'utm':
                mat1 = ax1.scatter(xyz_data[:,0], xyz_data[:,1], s=2, c=self.cloud_df_raster['elevation'], linewidths=0,  cmap="RdBu", vmin=np.min(Z), vmax=np.max(Z))
            if self.mode =='latlon':
                mat1 = ax1.scatter(lonlat_data[:,0], lonlat_data[:,1], s=2, c=self.cloud_df_raster['elevation'], linewidths=0,  cmap="RdBu", vmin=np.min(Z), vmax=np.max(Z))            
            mat2 = ax2.contourf(x, y, out_image, 100, cmap="RdBu", vmin=np.nanmin(out_image), vmax=np.nanmax(out_image))
            mat3 = ax3.contourf(x, y, out_z_rel_2d, 100, cmap="RdBu", vmin=np.nanmin(out_z_rel_2d), vmax=np.nanmax(out_z_rel_2d))
            mat4 = ax4.contourf(x, y, out_qual_flag_2d, 100, cmap="RdBu", vmin=np.nanmin(out_qual_flag_2d), vmax=np.nanmax(out_qual_flag_2d))
            ax1.set_title('Point cloud')
            ax2.set_title('Raster (iwd interpolation)')
            ax3.set_title('Relative height variance')
            ax4.set_title('Flag')
            fig.colorbar(mat1,label="Height (m)", orientation="vertical", ax=ax1)        
            fig.colorbar(mat2,label="Height (m)", orientation="vertical", ax=ax2)
            fig.colorbar(mat3,label="Relative height variance", orientation="vertical", ax=ax3)
            fig.colorbar(mat4,label="Flag", orientation="vertical", ax=ax4)
            plt.show()


    def compute_raster_exotic(self):
                        
        # Extract xyz
        if self.mode == 'utm':
            nb = len(self.cloud_df_raster['x'])
            xyz_data = np.ndarray((nb,3))
            xyz_data[:,0] = self.cloud_df_raster['x']
            xyz_data[:,1] = self.cloud_df_raster['y']

        # Extract latlon
        if self.mode == 'latlon':
            nb = len(self.cloud_df_raster['longitude'])
            lonlat_data = np.ndarray((nb,3))
            lonlat_data[:,0] = self.cloud_df_raster['longitude']
            lonlat_data[:,1] = self.cloud_df_raster['latitude']
            
        elevation = self.cloud_df_raster['elevation']
        
        # Train
        if self.mode == 'utm':
            idw_tree = cKDTree(xyz_data[:,0:2])
        if self.mode == 'latlon':
            idw_tree = cKDTree(lonlat_data[:,0:2])
            
        # Compute grid
        if self.mode == 'utm':
            # Interpolation
            x0 = min(xyz_data[:,0])
            y0 = min(xyz_data[:,1])
            nx = abs(int((x0-max(xyz_data[:,0]))/self.resolution))+1
            ny = abs(int((y0-max(xyz_data[:,1]))/self.resolution))+1
        if self.mode =='latlon':
            x0 = min(lonlat_data[:,0])
            y0 = min(lonlat_data[:,1])
            nx = abs(int((x0-max(lonlat_data[:,0]))/self.resolution))+1
            ny = abs(int((y0-max(lonlat_data[:,1]))/self.resolution))+1 
                
        x1 = x0 + self.resolution*(nx-1)
        y1 = y0 + self.resolution*(ny-1)
        
        x = np.linspace(x0, x1, nx)
        y = np.linspace(y0, y1, ny)
        
        GRID = np.meshgrid(x, y)
        grid_shape = GRID[0].shape
        GRID = np.reshape(GRID, (2, -1)).T    
        
        # Compute first interpolation in order to reduce to raster grid without interpolation (todo, inverse distance weighting ?)
        results = idw_tree.query_ball_point(GRID, self.resolution)
        Z = np.zeros_like(results, dtype = float)
        
        for i, k in enumerate(results):
            if len(k) != 0:
                Z[i] = np.mean(elevation[k])
                
        Z2d = Z.reshape(grid_shape)
        
        
        # Determine good (no empty) points to use
        good_points = np.where(Z2d != 0)
        points = np.dstack((good_points[0], good_points[1]))
        grid_x, grid_y = np.mgrid[range(Z2d.shape[0]), range(Z2d.shape[1])]
        
        # Compute cubic interpolation
        Z2d = interpolate.griddata(points[0], Z2d[np.where(Z2d != 0)], (grid_x, grid_y), method='cubic')
        

        if self.mode == 'latlon':
            # ~ dist_min = 2*np.pi*EARTH_RADIUS*dist_min/360.
            # ~ dist_mean = 2*np.pi*EARTH_RADIUS*dist_mean/360.
                        
            dist_min = np.zeros_like(results, dtype = float)
            dist_mean = np.zeros_like(results, dtype = float)
            z_rel = np.zeros_like(results, dtype = float)
            
            
        mask = shapefile.Reader(self.mask)
            
        proj = pyproj.Proj(init='EPSG:'+str(self.espg))    

        # Reproject into 2d grid
        # ~ Z2d = Z.reshape(grid_shape)
        dist_min_2d = dist_min.reshape(grid_shape)
        dist_mean_2d = dist_mean.reshape(grid_shape)
        z_rel_2d = z_rel.reshape(grid_shape)
        
        # Compute quality flag
        qual_flag_2d = compute_qual_flag(dist_mean_2d, z_rel_2d)
        
        # Compute mask in latlon and utm coordinates
        extracted_area_utm = []     
        extracted_area_latlon = []
        
        for extracted_area in mask.shapeRecords():
            extracted_area_utm.append(toFromUTM(extracted_area.shape, proj))
            extracted_area_latlon.append(extracted_area.shape)
            
        transform = rasterio.transform.from_origin(min(x), min(y), self.resolution, -self.resolution)
        
        #Filter data out of mask
        if self.mode == 'utm':
            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=Z2d.shape[0], width=Z2d.shape[1], count=1, dtype=Z2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(Z2d, 1)
                    out_image, out_transform = rasterio.mask.mask(src, extracted_area_utm, crop=False, indexes = 1, nodata=np.nan)

            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=dist_min_2d.shape[0], width=dist_min_2d.shape[1], count=1, dtype=dist_min_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(dist_min_2d, 1)
                    out_dist_min_2d, out_transform = rasterio.mask.mask(src, extracted_area_utm, crop=False, indexes = 1, nodata=np.nan)

            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=dist_mean_2d.shape[0], width=dist_mean_2d.shape[1], count=1, dtype=dist_mean_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(dist_mean_2d, 1)
                    out_dist_mean_2d, out_transform = rasterio.mask.mask(src, extracted_area_utm, crop=False, indexes = 1, nodata=np.nan)

            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=z_rel_2d.shape[0], width=z_rel_2d.shape[1], count=1, dtype=z_rel_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(z_rel_2d, 1)
                    out_z_rel_2d, out_transform = rasterio.mask.mask(src, extracted_area_utm, crop=False, indexes = 1, nodata=np.nan)
                    
            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=qual_flag_2d.shape[0], width=qual_flag_2d.shape[1], count=1, dtype=qual_flag_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(qual_flag_2d, 1)
                    out_qual_flag_2d, out_transform = rasterio.mask.mask(src, extracted_area_utm, crop=False, indexes = 1, nodata=np.nan)
                    
        if self.mode == 'latlon':
            with MemoryFile() as memfile:                
                with memfile.open(driver='GTiff', height=Z2d.shape[0], width=Z2d.shape[1], count=1, dtype=Z2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(Z2d, 1)
                    out_image, out_transform = rasterio.mask.mask(src, extracted_area_latlon, crop=False, indexes = 1, nodata=np.nan)        

            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=dist_min_2d.shape[0], width=dist_min_2d.shape[1], count=1, dtype=dist_min_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(dist_min_2d, 1)
                    out_dist_min_2d, out_transform = rasterio.mask.mask(src, extracted_area_latlon, crop=False, indexes = 1, nodata=np.nan)

            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=dist_mean_2d.shape[0], width=dist_mean_2d.shape[1], count=1, dtype=dist_mean_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(dist_mean_2d, 1)
                    out_dist_mean_2d, out_transform = rasterio.mask.mask(src, extracted_area_latlon, crop=False, indexes = 1, nodata=np.nan)

            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=z_rel_2d.shape[0], width=z_rel_2d.shape[1], count=1, dtype=z_rel_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(z_rel_2d, 1)
                    out_z_rel_2d, out_transform = rasterio.mask.mask(src, extracted_area_latlon, crop=False, indexes = 1, nodata=np.nan)
                                    
            with MemoryFile() as memfile:
                with memfile.open(driver='GTiff', height=qual_flag_2d.shape[0], width=qual_flag_2d.shape[1], count=1, dtype=qual_flag_2d.dtype, crs='EPSG:'+str(self.espg), transform=transform) as src:
                    src.write(qual_flag_2d, 1)
                    out_qual_flag_2d, out_transform = rasterio.mask.mask(src, extracted_area_latlon, crop=False, indexes = 1, nodata=np.nan)

            
        self.x = x
        self.y = y
        self.out_image = out_image
        self.out_dist_min_2d = out_dist_min_2d
        self.out_dist_mean_2d = out_dist_mean_2d
        self.out_qual_flag_2d = out_qual_flag_2d
        
        print("Output raster shape = ", Z2d.shape)
        
        if self.plot =='yes':
            
            fig, [[ax1, ax2],[ax3, ax4]] = plt.subplots(2,2, sharex=True, sharey=True, figsize=(18,18))
            if self.mode == 'utm':
                mat1 = ax1.scatter(xyz_data[:,0], xyz_data[:,1], s=2, c=self.cloud_df_raster['elevation'], linewidths=0,  cmap="RdBu", vmin=np.min(Z), vmax=np.max(Z))
            if self.mode =='latlon':
                mat1 = ax1.scatter(lonlat_data[:,0], lonlat_data[:,1], s=2, c=self.cloud_df_raster['elevation'], linewidths=0,  cmap="RdBu", vmin=np.min(Z), vmax=np.max(Z))            
            mat2 = ax2.contourf(x, y, out_image, 100, cmap="RdBu", vmin=np.nanmin(out_image), vmax=np.nanmax(out_image))
            mat3 = ax3.contourf(x, y, out_z_rel_2d, 100, cmap="RdBu", vmin=np.nanmin(out_z_rel_2d), vmax=np.nanmax(out_z_rel_2d))
            mat4 = ax4.contourf(x, y, out_qual_flag_2d, 100, cmap="RdBu", vmin=np.nanmin(out_qual_flag_2d), vmax=np.nanmax(out_qual_flag_2d))
            ax1.set_title('Point cloud')
            ax2.set_title('Raster (iwd interpolation)')
            ax3.set_title('Relative height variance')
            ax4.set_title('Flag')
            fig.colorbar(mat1,label="Height (m)", orientation="vertical", ax=ax1)        
            fig.colorbar(mat2,label="Height (m)", orientation="vertical", ax=ax2)
            fig.colorbar(mat3,label="Relative height variance", orientation="vertical", ax=ax3)
            fig.colorbar(mat4,label="Flag", orientation="vertical", ax=ax4)
            plt.show()
  
    def write_fpdem_raster(self):
        # Write output file
        write_raster_gridded(self.output_file, self.mode, self.x, self.y, self.out_image, self.out_dist_min_2d, self.out_dist_mean_2d, self.out_qual_flag_2d, self.resolution, self.espg)
    


def compute_qual_flag(mean_dist, mean_z_rel):
    flag = np.zeros_like(mean_z_rel)
    flag[:] = 1
    flag[np.where(mean_z_rel > 0.1)] = 2
    flag[np.where(mean_dist > 500.)] = 3
    flag[np.logical_and(mean_z_rel > 0.30,mean_dist > 500.)] = 4
    
    return flag

def toFromUTM(shp, proj, inv=False):
    geoInterface = shp.__geo_interface__

    shpType = geoInterface['type']
    coords = geoInterface['coordinates']
    if shpType == 'Polygon':
        newCoord = [[proj(*point, inverse=inv) for point in linring] for linring in coords]
    elif shpType == 'MultiPolygon':
        newCoord = [[[proj(*point, inverse=inv) for point in linring] for linring in poly] for poly in coords]

    return shpgeo.shape({'type': shpType, 'coordinates': tuple(newCoord)})         
            

# Main program
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute fpdem raster product")
    parser.add_argument("parameter_file", help="parameter_file (*.rdf)")
    args = parser.parse_args()
    parameters = my_rdf.myRdfReader(args.parameter_file)
    
    fpdem_raster = FPDEM_Raster(parameters)
    fpdem_raster.load_input_fpdem_raster()
    fpdem_raster.compute_fpdem_raster()
    fpdem_raster.write_fpdem_raster()
