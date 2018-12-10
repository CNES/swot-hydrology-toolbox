import cnes.modules.geoloc.lib.pixel_cloud as pixel_cloud
import cnes.modules.geoloc.lib.pixel_cloud_SGE as pixel_cloud_sge
import cnes.modules.geoloc.lib.sensor_SGE as sensor_sge

import warnings
import numpy as np
import matplotlib.pyplot as plt
import math 
from scipy.interpolate import griddata
from scipy import spatial
#~ from scipy.signal import convolve2d as convolve
from scipy.signal import convolve as convolve

import scipy.signal

def taille_pix(col, nr, h, pd, alt):
    # function to estimate pixel size (TO BE CHANGED USING PIXEL CLOUD PIXEL SIZE)
    # col :  radar column
    # nr : Near range value (m)
    # h : satellite altitude (m)
    # pd : range gate (m)
    # alt : groung mean height (m)
    theta = math.acos((alt-h)/(nr+col))
    return pd/math.sin(theta)

def mean_classif(tab,classif):
    # Estimate the averared value of an array
    classif = np.where(np.isnan(tab) == False, classif, 0)
    classif = np.where(classif != 0, 1, 0)
    tab = np.where(np.isnan(tab) == False, tab, 0)
    tot = np.sum(tab*classif)

    
    if tot != 0.:
        return tot/np.sum(classif)
    else:
        return 0. 

def mean_classif_slidding(tab,classif,classif_center):
    # Estimate the averared value of an array
    classif = np.where(np.isnan(tab) == False, classif, 0)
    classif = np.where(classif != 0, 1, 0)
    
    #~ classif_center= np.where(np.isnan(tab) == False, classif_center, 0)
    classif_center = np.where(classif_center != 0, 1, 0)

    tab = np.where(np.isnan(tab) == False, tab, 0)
    tot = np.sum(tab*classif)
    
    tot_classif_center = np.sum(classif_center)
    tot_classif = np.sum(classif)
    
    #~ if tot_classif > 1500. :    
    if tot_classif_center != 0. :
    #~ if tot_classif != 0. :
        return tot/np.sum(classif)
    else:
        return 0. 
                

def local_mean(image, kernel, nodata_mask, mode="same", no_data_output_value=np.nan):
    """
    Local average of image, with a nodata mask

    Parameters:
        image: input image
        kernel: averaging kernel, typically np.ones(N, N)
        nodata_mask: binary mask (same shape as image), True where there is nodata
        mode: mode for scipy's convolve2d

    Returns:
        Averaged image
    """

    # To handle no-data, use 0 inplace of invalid values,
    # then convolve numerator and denominator separatly
    # essentially dividing by the number of valid pixels in the area (instead of simply kernel.sum())
    # originally from https://stackoverflow.com/a/42029965/5815110
    numerator = convolve(np.where(nodata_mask, 0, image), kernel, mode=mode)
    denominator = convolve(np.logical_not(nodata_mask), kernel, mode=mode)

    # Silence the divide by zero warning
    # No big deal because we replace by NO_DATA_OUTPUT_VALUE right after
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = np.where(denominator != 0, numerator / denominator, no_data_output_value)

    # Reset nodata pixels to nodata value
    result = np.where(nodata_mask, no_data_output_value, result)
    return result
    
    
class Geoloc_big_lake:

    def __init__(self):
        pass
        
    def construc_lakes_radar_geom(self, pixc):
        
        nblabels = len(pixc.kept_labels)
        # number of colums of radar image (max(range indice)+1)
        nc = pixc.nr_pixels
        # number of lines of radar image (max(range line)+1)
        nl = pixc.nr_lines
            
        # Initialisation of classif, height, taille_pix arrays
        self.classif = np.zeros([nl,nc,nblabels],int)
        self.height = np.zeros([nl,nc,nblabels], float)
        self.taille_pix_tab = np.zeros([nl,nc,nblabels], float)
        self.height_averaged = np.zeros([nl,nc,nblabels],float)

        self.x = np.zeros([nl,nc],float)
        self.y = np.zeros([nl,nc],float)
        
        for i in range(nl):
            self.y[i,:]= np.arange(nc)
        for i in range(nc):
            self.x[:,i]= np.arange(nl)
            
        for n in np.arange(len(pixc.kept_labels)):
            
            ind_lab = np.where(pixc.labels == pixc.kept_labels[n])
            print('labels = ', pixc.kept_labels[n])
            for i in np.arange(len(ind_lab)):
                self.classif[pixc.azimuth_idx[ind_lab[i]],pixc.range_idx[ind_lab[i]], n] = 1
                self.height[pixc.azimuth_idx[ind_lab[i]],pixc.range_idx[ind_lab[i]], n] = pixc.height[ind_lab[i]]
                self.taille_pix_tab[pixc.azimuth_idx[ind_lab[i]],pixc.range_idx[ind_lab[i]], n] = pixc.pixel_area[ind_lab[i]]

            
    def compute_averaged_grid(self, pixc, sensor, av_pixel_size):
    
        ### Calculate the averaged grid ###
        
        # number of colums of radar image (max(range indice)+1)
        nc = pixc.nr_pixels
        # number of lines of radar image (max(range line)+1)
        nl = pixc.nr_lines
        
        ## Mean altitude
        altitude = sensor.altitude
        altitude_mean = np.mean(altitude)
        
        ## Pixel azimuth size (TO BE CHANGED USING METADATA, not in actual big scale simulator)
        azimuth_spacing = sensor.azimuth_spacing

        ## Range gate value (TO BE CHANGED USING METADATA, not in actual big scale simulator)
        range_spacing = 0.75
        
        ## Mean Near-range (TO BE CHANGED USING METADATA, not in actual big scale simulator)
        near_range_mean = altitude_mean + 50.

        self.tab_ind_l=[]
        self.tab_ind_c=[]

        for n in np.arange(len(pixc.kept_labels)):
            
            ## hmean (TO BE CHANGED USING MORE PRECISE VALUE IF NEEDED)
            hmean = mean_classif(self.height[:,:,n],self.classif[:,:,n])
            print(hmean)
            # range indice for a size of av_pixel_size[1]
            total = 0.
            k=1
            ind_c=[0]
            for i in range(nc):
                total += np.sum(taille_pix(i, near_range_mean, hmean, range_spacing, altitude_mean))
                if total > k*av_pixel_size[1] :
                    k+=1.
                    ind_c.append(i)
            ind_c.append(nc)

            # azimuth indice for a size of av_pixel_size[0]
            total = 0.
            k=1
            ind_l=[0]
            for i in range(nl):
                total += azimuth_spacing
                if total > k*av_pixel_size[0] :
                    k+=1.
                    ind_l.append(i)
            ind_l.append(nl)    
                
            self.tab_ind_l.append(ind_l)
            self.tab_ind_c.append(ind_c)

        
    def compute_averaged_height(self, pixc, k_window):

        # number of colums of radar image (max(range indice)+1)
        nc = pixc.nr_pixels
        # number of lines of radar image (max(range line)+1)
        nl = pixc.nr_lines
        
        for n in np.arange(len(pixc.kept_labels)):

            '''      
            ind_l = self.tab_ind_l[n]
            ind_c = self.tab_ind_c[n]

            h_km = np.zeros([len(ind_l)-1,len(ind_c)-1])
            x_km = np.zeros([len(ind_l)-1,len(ind_c)-1])
            y_km = np.zeros([len(ind_l)-1,len(ind_c)-1])
            h_km_list = []
            x_km_list = []
            y_km_list = []

           
            # loop on x
            for i in range(len(ind_l)-1):
                # loop on y
                for j in range(len(ind_c)-1):
                    # Averaging of x, y on initial window without overlap
                    x_km[i,j] = mean_classif(self.x[ind_l[i]:ind_l[i+1],ind_c[j]:ind_c[j+1]],self.classif[ind_l[i]:ind_l[i+1],ind_c[j]:ind_c[j+1],n])
                    y_km[i,j] = mean_classif(self.y[ind_l[i]:ind_l[i+1],ind_c[j]:ind_c[j+1]],self.classif[ind_l[i]:ind_l[i+1],ind_c[j]:ind_c[j+1],n])
                    
                    ## Averaging of heights with overlap
                 
                    h_km[i,j] = mean_classif_slidding(self.height[ind_l[max(i-k_window,0)]:ind_l[min(i+1+k_window,len(ind_l)-1)],ind_c[max(j-k_window,0)]:ind_c[min(j+1+k_window,len(ind_c)-1)],n],self.classif[ind_l[max(i-k_window,0)]:ind_l[min(i+1+k_window,len(ind_l)-1)],ind_c[max(j-k_window,0)]:ind_c[min(j+1+k_window,len(ind_c)-1)],n], self.classif[ind_l[i]:ind_l[i+1],ind_c[j]:ind_c[j+1],n])
                   
                    # Adding results in list for interpolation
                    
                    if np.isnan(h_km[i,j]) == False:
                        h_km_list.append(h_km[i,j])       
                        x_km_list.append(x_km[i,j])
                        y_km_list.append(y_km[i,j])
            # x,y grid at radar sampling (azimuth, range)        
            grid_x, grid_y = np.mgrid[0:nl, 0:nc]

            # 2D tab with x,y coordinates associated to averaged heights
            values = np.zeros([len(h_km_list),2])
            values[:,0] = x_km_list
            values[:,1] = y_km_list
            
            # Interpolation 
            height_averaged = griddata(values, h_km_list, (grid_x, grid_y), method='cubic', fill_value = 'nan')
            height_averaged = np.where(self.classif[:,:,n] == 0, np.nan, height_averaged)
            self.height_averaged[:,:,n]= height_averaged
            self.h_km = h_km
            '''
        
            self.height_averaged[:,:,n] = local_mean(self.height[:,:,n], np.ones([501,501]), np.logical_not(self.classif[:,:,n]), mode="same")
        
    def add_height_variation_on_pixel_cloud(self, pixc):
        
        for n in np.arange(len(pixc.kept_labels)):
            self.height[:,:,n] += np.array([np.arange(pixc.nr_pixels),]*pixc.nr_lines)/100.
            self.height[:,:,n] = np.where(self.classif[:,:,n] == 0, np.nan, self.height[:,:,n])


# Open sensor file
##~ input_sensor_path= "/work/ALT/swot/swotdev/outils_swot/DATA/SWOT_L2_HR_PIXC_001_004_46N-R_sensor.nc"
#~ input_sensor_path= "/work/ALT/swot/swothr/users/desrochesd/runs/lacs.sis/sisimp_pos_2.00/outputs/SWOT_L2_HR_PIXC_000_060_45N-R_sensor.nc"
#~ input_sensor_path= "/work/ALT/swot/swothr/users/desrochesd/runs/leman.sis/sisimp_pos_2.00/outputs/SWOT_L2_HR_PIXC_000_004_46N-R_sensor.nc"
input_sensor_path = "/home/poughov/work/workdir/sisimp_test/mysterious_lake/output/SWOT_L2_HR_PIXC_000_060_45N-R_sensor.nc"

# Open pixel cloud file
##~ input_pixc_path= "/work/ALT/swot/swotdev/outils_swot/DATA/SWOT_L2_HR_PIXC_001_004_46N-R_main.nc"
#~ input_pixc_path= "/work/ALT/swot/swothr/users/desrochesd/runs/lacs.sis/sisimp_pos_2.00/outputs/SWOT_L2_HR_PIXC_000_060_45N-R_main.nc"
#~ input_pixc_path= "/work/ALT/swot/swothr/users/desrochesd/runs/leman.sis/sisimp_pos_2.00/outputs/SWOT_L2_HR_PIXC_000_004_46N-R_main.nc"
input_pixc_path = "/home/poughov/work/workdir/sisimp_test/mysterious_lake/output/SWOT_L2_HR_PIXC_000_060_45N-R_main.nc"






classification_values_pixc=[4]
size_entities=1000*1000
av_pixel_size = [1000, 1000]
k_window = 2

objPixc = pixel_cloud_sge.PixelCloudSGE(input_pixc_path, classification_values_pixc)
objSens = sensor_sge.SensorSGE(input_sensor_path)

objPixc.computeSeparateEntities()
objPixc.filterShorterEntities(size_entities)

objBigLake = Geoloc_big_lake()
objBigLake.construc_lakes_radar_geom(objPixc)
objBigLake.compute_averaged_grid(objPixc, objSens, av_pixel_size)

#~ objBigLake.add_height_variation_on_pixel_cloud(objPixc)

objBigLake.compute_averaged_height(objPixc, k_window)


fig = plt.figure(figsize=[12,8])
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
im1 = ax1.imshow(objBigLake.height[4800:5100,2500:,0])
#~ im1 = ax1.imshow(objBigLake.classif[4800:5100,2500:,0])
im2 = ax2.imshow(objBigLake.height_averaged[4800:5100,2500:,0])

from scipy.misc import imsave, toimage
print(objBigLake.height_averaged[4800:5100,2500:,0].shape)
toimage(objBigLake.height_averaged[4800:5100,2500:,0], cmin=-3, cmax=-2).save("h_ave.png")

fig.colorbar(im2)


#~ fig = plt.figure(figsize=[8,8])
#~ im1 = plt.imshow(objBigLake.h_km[100:,40:])
#~ fig.colorbar(im1)

#~ fig = plt.figure()
#~ im1 = plt.plot(np.nanmean(objBigLake.height_averaged[:,:,0],axis=1))
#~ im2 = plt.plot(np.nanmean(objBigLake.height[:,:,0],axis=1))

plt.show()
