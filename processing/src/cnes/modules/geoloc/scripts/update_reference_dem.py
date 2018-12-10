#!/usr/bin/env python

import sys
import os
import os.path

import cnes.modules.geoloc.lib.tools as lib
import cnes.modules.geoloc.lib.my_netcdf_file as myCdf
import cnes.modules.geoloc.lib.gdem as gdem
import cnes.modules.geoloc.lib.function_mean_h as function_mean_h
import cnes.modules.geoloc.lib.pixel_cloud as pixel_cloud
import utm

import numpy as np

import matplotlib.pyplot as plt
from scipy import spatial
from timeit import default_timer as timer
import itertools as itertools

import dask.bag as db
import numpy as np
import dask.array as da
import dask.multiprocessing
import dask.threaded
from operator import add
from dask import delayed
from dask.distributed import Client, progress

     
class Update_reference_dem(object):
    
    def __init__(self):
        pass

    def utm_from_latlon(self, tablat, tablon):
        tab = np.zeros([2,len(tablat)], float)
        for i in np.arange(len(tablat)):
            tab[:,i] = utm.from_latlon(tablat[i], tablon[i])[0:2]
        return tab
        
    def create_grid(self, pixc, grid_size_x, grid_size_y):
        '''
        Calculate the raster grid with grid_size_x, grid_size_y resolution corresponding to the pixel cloud area
        '''   
        self.xres = grid_size_x
        self.yres = grid_size_y
        
        pixc.latmin = min(pixc.latitude)
        pixc.latmax = max(pixc.latitude)
        pixc.lonmin = min(pixc.longitude)
        pixc.lonmax = max(pixc.longitude)
        
        c0 = utm.from_latlon(pixc.latmin, pixc.lonmin)
        c1 = utm.from_latlon(pixc.latmin, pixc.lonmax)
        c2 = utm.from_latlon(pixc.latmax, pixc.lonmin)
        c3 = utm.from_latlon(pixc.latmax, pixc.lonmax)
        
        #~ c0 = lib.llh2xyz(pixc.lonmin, pixc.latmin, 0, IN_flag_rad=False)
        #~ c1 = lib.llh2xyz(pixc.lonmax, pixc.latmin, 0, IN_flag_rad=False)
        #~ c2 = lib.llh2xyz(pixc.lonmin, pixc.latmax, 0, IN_flag_rad=False)
        #~ c3 = lib.llh2xyz(pixc.lonmax, pixc.latmax, 0, IN_flag_rad=False)
    
        pixc.xmin = min(c0[0], c1[0], c2[0], c3[0])
        pixc.xmax = max(c0[0], c1[0], c2[0], c3[0])
        pixc.ymin = min(c0[1], c1[1], c2[1], c3[1])
        pixc.ymax = max(c0[1], c1[1], c2[1], c3[1])
        
        print("ymin ymax xmin xmax = ", pixc.ymin, pixc.ymax, pixc.xmin, pixc.xmax)
        
        self.nx = int((pixc.xmax-pixc.xmin)/self.xres)+1
        self.ny = int((pixc.ymax-pixc.ymin)/self.yres)+1
        
        self.dem = np.zeros([self.nx, self.ny, 3], float)
        self.dem[:,:,0] = np.array([pixc.xmin + self.xres*np.arange(self.nx),]*self.ny).transpose()
        self.dem[:,:,1] = np.array([pixc.ymin + self.yres*np.arange(self.ny),]*self.nx)

        print("Size DEM = ", len(self.dem[:,:,0]), len(self.dem[0,:,0]))

        print("nx ny = ", self.nx, self.ny)
        
    def estimate_averaged_mean_height(self, pixc):
        '''
        Calculate the height averaged on a GDEM grid using pixel cloud height coordinates
        '''    
              
        htot = np.zeros([self.nx,self.ny], float)
        count = np.zeros([self.nx,self.ny], float) 
        
        for i in self.kept_labels:
            self.kept_labels_hmean = np.mean(self.height[np.where(pixc.labels == i)])
            
            h = pixc.height[np.where(pixc.labels == i)]
            lat = pixc.latitude[np.where(pixc.labels == i)]
            lon = pixc.longitude[np.where(pixc.labels == i)]
            
            #~ xyz = lib.llh2xyz(lon, lat, h, IN_flag_rad=False)
            xyz = self.utm_from_latlon(lat, lon)
            
            
            case_x = np.floor((abs(self.dem[self.nx-1,0,0] - xyz[0]))/self.xres)
            case_y = np.floor((abs(self.dem[0,self.ny-1,1] - xyz[1]))/self.yres)
               

            for k in np.arange(len(h)):
                            
                htot[int(case_x[k]),int(case_y[k])] += h[k]
                
                count[int(case_x[k]),int(case_y[k])] += 1.
                
        count = np.where(count == 0, 1, count)
        self.dem[:,:,2] = htot/count
        self.dem[:,:,2] = np.where(self.dem[:,:,2] != 0, self.dem[:,:,2], np.NaN)
        
    def estimate_averaged_mean_height_quadtree(self, pixc, nb_points, dist_max):
        htot = np.zeros([self.nx,self.ny], float)
        count = np.zeros([self.nx,self.ny], float) 
        
        #~ xyz = lib.llh2xyz(pixc.longitude, pixc.latitude, pixc.height, IN_flag_rad=False)
        xyz = self.utm_from_latlon(pixc.latitude, pixc.longitude)
        
        for i in self.kept_labels:
            
            
            print("Calcul moyenne de hauteur glissante pour labels: ", i)
            ind_lab = np.where(pixc.labels == i)
            pts = list(zip(xyz[0][ind_lab], xyz[1][ind_lab]))
            h = pixc.height[ind_lab]
            tree = spatial.KDTree(pts)
            
            for m in range(self.nx):
                for n in range(self.ny):
                     xi, yi = self.dem[m, n, 0], self.dem[m, n, 1]
                     if tree.query([xi, yi],k=1)[0] < dist_max:
                        dist, ind = tree.query([xi, yi],k = nb_points)
                        self.dem[m, n, 2] = np.mean(h[ind])


        self.dem[:,:,2] = np.where(self.dem[:,:,2] != 0, self.dem[:,:,2], np.NaN)


    def estimate_averaged_mean_height_quadtree_dask(self, pixc, nb_points, dist_max):
        htot = np.zeros([self.nx,self.ny], float)
        count = np.zeros([self.nx,self.ny], float)
         
        #~ xyz = lib.llh2xyz(pixc.longitude, pixc.latitude, pixc.height, IN_flag_rad=False)
        xyz = self.utm_from_latlon(pixc.latitude, pixc.longitude)

        points_list=[]
        for i in pixc.kept_labels:
            print("Calcul moyenne de hauteur glissante pour labels: ", i)
            ind_lab = np.where(pixc.labels == i)
            pts = list(zip(xyz[0][ind_lab], xyz[1][ind_lab]))
            h = pixc.height[ind_lab]
            points_list.append([pts,h])

        points_bag = db.from_sequence(points_list)
            
        f = function_mean_h.Function_mean_h(self.nx, self.ny, self.dem, nb_points)
        dem = points_bag.map(f.calcul_h, dist_max)
        dems = dem.compute()        
        self.dem[:,:,2] = sum(dems[:])[:,:]
        self.dem[:,:,2] = np.where(self.dem[:,:,2] != 0, self.dem[:,:,2], np.NaN)
        

    def estimate_averaged_mean_height_quadtree_daskarray(self, pixc, nb_points, dist_max):
        htot = np.zeros([self.nx,self.ny], float)
        count = np.zeros([self.nx,self.ny], float) 
        
        #~ xyz = lib.llh2xyz(pixc.longitude, pixc.latitude, pixc.height, IN_flag_rad=False)
        xyz = self.utm_from_latlon(pixc.latitude, pixc.longitude)
        
        points_list=[]
        
                    
        for i in self.kept_labels:
            print("Calcul moyenne de hauteur glissante pour labels: ", i)
            ind_lab = np.where(pixc.labels == i)
            pts = list(zip(xyz[0][ind_lab], xyz[1][ind_lab]))
            h = pixc.height[ind_lab]
            tree = spatial.KDTree(pts)
            daskarray = da.from_array(self.dem[:,:,0:2], chunks=(20,20,2))
            f = function_mean_h.Function_mean_h(self.nx, self.ny, self.dem, nb_points, tree)
            dem = daskarray.map_blocks(f.calcul_h_array, h, dist_max, dtype=float, chunks=(20,20), drop_axis=2).compute()
            #~ dem = daskarray.map_blocks(f.calcul_h_array, h, dist_max, dtype=float, chunks=(20,20), drop_axis=2).compute(get=dask.threaded.get)
            #~ dem = daskarray.map_blocks(f.calcul_h_array, h, dist_max, dtype=float, chunks=(20,20), drop_axis=2).compute(get=dask.multiprocessing.get)
            
            self.dem[:,:,2]+=dem[:,:]

        
    def estimate_averaged_mean_height_quadtree_daskmap(self, pixc, nb_points, dist_max):
        
        htot = np.zeros([self.nx,self.ny], float)
        count = np.zeros([self.nx,self.ny], float) 
        
        #~ xyz = lib.llh2xyz(pixc.longitude, pixc.latitude, pixc.height, IN_flag_rad=False)
        xyz = self.utm_from_latlon(pixc.latitude, pixc.longitude)
        
        points_list=[]
        
        def calcul_h_map(array, h, tree, nb_points, dist_max):
            xi, yi = array[0], array[1]
            if tree.query([xi, yi],k=1)[0] < dist_max:
                dist, ind = tree.query([xi, yi],k = nb_points)
                dem = np.mean(h[ind])

            return dem
            
                    
        for i in self.kept_labels:
            print("Calcul moyenne de hauteur glissante pour labels: ", i)
            ind_lab = np.where(pixc.labels == i)
            pts = list(zip(xyz[0][ind_lab], xyz[1][ind_lab]))
            h = pixc.height[ind_lab]
            tree = spatial.KDTree(pts)
            dem_reshape = np.reshape(self.dem[:,:,0:2], [self.nx * self.ny, 2])
            daskbag = db.from_sequence(dem_reshape)
            f = function_mean_h.Function_mean_h(self.nx, self.ny, self.dem, nb_points)
            dem = daskbag.map(f.calcul_h_map, h, tree, dist_max).compute()
            #~ dem = da.map_blocks(calcul_h_array, daskarray, h, tree, nb_points, dist_max, dtype=float, chunks=(20,20), drop_axis=2).compute(get=dask.threaded.get)
            #~ dem = da.map_blocks(calcul_h_array, daskarray, h, tree, nb_points, dist_max, dtype=float, chunks=(20,20), drop_axis=2).compute(get=dask.multiprocessing.get)

            self.dem[:,:,2]+=reshape(dem, [self.nx, self.ny])
        
                
    def estimate_averaged_mean_height_quadtree_superdask(self, pixc, nb_points, dist_max):
        
        def points_to_multipoints(points, nb_voisins):
            x=points[0]
            y=points[1]
            h=points[2]
            lab=points[3]
            indx = int(abs((x - pixc.xmin)//self.xres))
            indy = int(abs((y - pixc.ymin)//self.yres))

            return [((i, j, lab), (x, y, h)) for i in range(max(0, indx-nb_voisins), min(indx+nb_voisins, self.nx-1))
                for j in range(max(0, indy-nb_voisins), min(indy+nb_voisins, self.ny-1)) ]

        def cle(x):
            return tuple(x[0])

                        
        htot = np.zeros([self.nx,self.ny], float)
        count = np.zeros([self.nx,self.ny], float) 
        
        #~ xyz = lib.llh2xyz(pixc.longitude, pixc.latitude, pixc.height, IN_flag_rad=False)
        xyz = self.utm_from_latlon(pixc.latitude, pixc.longitude)
        
        
        points_list=[]
        f = function_mean_h.Function_mean_h(self.nx, self.ny, self.dem, nb_points)
        
        for i in pixc.kept_labels:
            ind_lab = np.where(pixc.labels == i)
            pts = list(zip(*[xyz[0][ind_lab], xyz[1][ind_lab], pixc.height[ind_lab], np.full(np.size(ind_lab), i)]))
            points_list+=pts
        points_bag = db.from_sequence(points_list[0:], npartitions=1)
        #~ points_bag = db.from_sequence(points_list[0:])

        indexed_points = points_bag.map(points_to_multipoints, nb_voisins= int(5000/min(self.xres, self.yres)/2))
        indexed_points_flat = indexed_points.flatten()
        indexed_points_fb = indexed_points_flat.foldby(cle, lambda l, e : l + [e[1]], [], add)
        mean_height_by_id = indexed_points_fb.map(f.calcul_h_map_tuple, dist_max).compute()

        
        for i in mean_height_by_id:
            self.dem[i[0][0],i[0][1],2] = i[1] 
        
        
    def compute_upgraded_dem(self, dem_res_x, dem_res_y):
        
        size = self.dem[:,:,0].size
        values = np.zeros([size,3])
        values[:,0] = np.reshape(self.dem[:,:,0],size)
        values[:,1] = np.reshape(self.dem[:,:,1],size)
        values[:,2] = np.reshape(self.dem[:,:,2],size)

        values0=np.delete(values[:,0],np.where(np.isnan(values[:,2])))
        values1=np.delete(values[:,1],np.where(np.isnan(values[:,2])))
        values2=np.delete(values[:,2],np.where(np.isnan(values[:,2])))

        values = np.zeros([values0.size,3])
        values[:,0] = values0
        values[:,1] = values1
        values[:,2] = values2
        
        grid_x, grid_y = np.mgrid[self.xmin:self.xmax:dem_res_x,self.ymin:self.ymax:dem_res_y]
        self.upgraded_dem = griddata(values[:,0:2], values[:,2],(grid_x, grid_y), method='linear')

        
        
if __name__ == "__main__":

    start = timer()

    # Input Pixel Cloud and GDEM
    #~ input_pixc_path = "/work/ALT/swot/swotdev/SIMU_JPL/results_proc_simu_US_land_-5dB_no_dark_water/proc_simu9/cycle_0001_pass_0385_LeftSwath_nlcd50_darkWater_25m_CBE/output/cycle_0001_pass_0385_LeftSwath_nlcd50_darkWater_25m_CBE/input/slc.Noise.LeftSwath.Flat_pixc.nc"
    #~ input_gdem_path = "/work/ALT/swot/swotdev/outils_swot/DATA/gdem_dem_3426_NLCD.nc"
    #~ input_pixc_path = "/work/ALT/swot/swotdev/SIMU_JPL/results_proc_simu_US_land_-5dB_no_dark_water/proc_simu14/cycle_0001_pass_0051_RightSwath_nlcd50_darkWater_25m_CBE/output/cycle_0001_pass_0051_RightSwath_nlcd50_darkWater_25m_CBE/input/slc.Noise.RightSwath.Flat_pixc.nc"
    #~ input_gdem_path = "/work/ALT/swot/swotdev//Livraison_JPL_11_10/GDEM_DEM/gdem_dem_NLCD_2/gdem_dem_3420_NLCD.nc"
    input_pixc_path = "/work/ALT/swot/swotdev/SIMU_JPL/results_proc_simu_US_land_-5dB_no_dark_water/proc_simu15/cycle_0001_pass_0079_RightSwath_nlcd50_darkWater_25m_CBE/output/cycle_0001_pass_0079_RightSwath_nlcd50_darkWater_25m_CBE/input/slc.Noise.RightSwath.Flat_pixc.nc"
    input_gdem_path = "/work/ALT/swot/swotdev//Livraison_JPL_11_10/GDEM_DEM/gdem_dem_NLCD_2/gdem_dem_3465_NLCD.nc"

    # Output Directory (NOT USED YET)
    output_dir = "/work/ALT/swot/swotdev/outils_swot/DATA/"
    # Grid size of the computed new GDEM from pixel cloud heights
    grid_size=[100,100]
    # Classification value from pixel cloud used to compute GDEM water height
    classification_values_pixc = [4,14]
    # Classification value from initial gdem (used to validate results)
    classification_values_gdem = [1]
    # Size of kept water entities
    size_entities = 2000.*2000.
    # Number of pixels used during averaging
    nb_points_to_average = 1000
    # Distance max between gdem center of pixel and any pixel from pixel cloud to calculate a water averaged value for the gdem pixel
    dist_max = 50.
    

    objPixc = pixel_cloud.PixelCloud(input_pixc_path, classification_values_pixc)
    objPixc.computeSeparateEntities()
    objPixc.filterShorterEntities(size_entities)
    objUpdateRefDem = Update_reference_dem()
    objUpdateRefDem.create_grid(objPixc, grid_size[0],grid_size[1])
    
    objUpdateRefDem.estimate_averaged_mean_height_quadtree_dask(objPixc, nb_points_to_average, dist_max)
    #~ objUpdateRefDem.estimate_averaged_mean_height_quadtree(objPixc, nb_points_to_average, dist_max)
    #~ objUpdateRefDem.estimate_averaged_mean_height_quadtree_daskarray(objPixc, nb_points_to_average, dist_max)
    #~ objUpdateRefDem.estimate_averaged_mean_height_quadtree_daskmap(objPixc, nb_points_to_average, dist_max)
    #~ objUpdateRefDem.estimate_averaged_mean_height_quadtree_superdask(objPixc, nb_points_to_average, dist_max)    
    

    ### Obsolete
    ##~ objPixc.estimate_averaged_mean_height()
    ##~ objPixc.compute_upgraded_dem(final_gdem_res[0],final_gdem_res[1])
    ###
    
    #~ objGdem = gdem.Gdem(input_gdem_path, classification_values_gdem)
    #~ objGdem.compute_input_gdem_dem(objPixc, grid_size[0],grid_size[1])


    end = timer()
    print("Temps = ", end - start)

    plt.figure()
    plt.imshow(objPixc.sepEntities)
    plt.figure()
    plt.imshow(((objUpdateRefDem.dem[::,::,2]).T)[::-1,::])
    #~ plt.figure()
    #~ plt.imshow(objGdem.crop_elevation)
    #~ plt.figure()
    #~ plt.imshow(objGdem.crop_landtype)
    #~ plt.figure()
    #~ plt.imshow(objGdem.crop_dem_xyz.T)
    #~ plt.figure()
    #~ plt.imshow(objGdem.crop_landtype_xyz.T)
    #~ plt.figure()
    #~ plt.imshow(np.where(objGdem.crop_landtype_xyz == 1, objGdem.crop_dem_xyz, np.NaN).T)
    #~ plt.figure()
    #~ plt.imshow(np.where(objGdem.crop_landtype_xyz == 1, objGdem.crop_dem_xyz, np.NaN).T - objUpdateRefDem.dem[::,::, 2].T)


    plt.show()

