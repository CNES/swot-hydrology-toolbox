# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 14:35:20 2018

@author: plousser
"""
import numpy as np
import random
import lib.height_model as height_model

from skimage.measure import label
import random

# to move in other file and import
def dark_water_simulation(taille_2D,fact_echelle_dw, pourcent_dw,seedvalue=None):
        #create a 2D gaussian profile
        if seedvalue is not None:
            np.random.seed(int(seedvalue))
        profile_2d = height_model.generate_2d_profile_gaussian(taille_2D,0.,'Default',1,fact_echelle_dw)
        ### Create Dark water mask
        # dfine the threshold value to keep the desired % of dark water
        threshold_value = np.percentile(profile_2d,100-pourcent_dw)
        #initiate the dw_mask with zeros and set to 1 the dark water values
        mask_dw = np.zeros(profile_2d.shape)
        mask_dw[np.where(profile_2d>threshold_value)] = 1
        return mask_dw
    
    

def dark_water_randomerase(mask_dw,water_pixmatrix,taille_2D,seedvalue=None):
            #label dark_water regions
            mask_regions=label(mask_dw)

            #reshape to water extent
            mask_regions=mask_regions[taille_2D]
            mask_dw=mask_dw[taille_2D]
            ### Delete random regions
            #create a vector with each region label
            region_list=[i for i in range(1,np.max(mask_regions))]
            
            #Randomly select the regions to be deleted. the number of deleted regions is also random
            if seedvalue is not None:
                random.seed(seedvalue)
            region_todelete = random.sample(region_list,np.random.choice(region_list))
            print("Dark water non detected pixels : %d" % len(np.where(np.isin(mask_regions, region_todelete))[0]))

            #Erase non detected dark_water pixels
            water_pixmatrix[np.where(np.isin(mask_regions, region_todelete))]=0
            mask_dw=np.delete(mask_dw,np.where(np.isin(mask_regions, region_todelete)))
            #mask_regions=np.delete(mask_regions,np.where(np.isin(mask_regions, region_todelete)))
            
            
            return mask_dw, water_pixmatrix


def dark_water_non_detected_simulation(mask_dw,fact_echelle_nondetected_dw,percent_detected_dw,seedvalue=None):
    # label dark water regions
    mask_regions=label(mask_dw)
    region_list=[i for i in range(1,np.max(mask_regions+1))]

    #initialize the array of non detected dark water mask
    non_detected_dw_mask=np.zeros(mask_dw.shape)
    #set the value to a number different from the generated 2D profile value
    non_detected_dw_mask.fill(-999)
    # test if at least one region is present
    if len(region_list)>0 :
        for i in region_list :
            #Locate the pixels corresponding to the region i
            region_ind = np.where(mask_regions==i)
            # Get the bounding lines and columns of the region
            minx, maxx, miny, maxy = min(region_ind[0]),max(region_ind[0]),min(region_ind[1]),max(region_ind[1])
            size_x = maxx-minx+1
            size_y=maxy-miny+1

            #Simulate non detected dark water
            profile_2d = height_model.generate_2d_profile_gaussian([size_x,size_y],0.,'Default',1,fact_echelle_nondetected_dw)
            profile_2d[np.where(mask_regions[minx:maxx+1,miny:maxy+1]!=i)]=-999
            non_detected_dw_mask[minx:maxx+1,miny:maxy+1]=profile_2d
        # Define the threshold value to keep the percentage of non_detected dark water inside detected dark_water regions
        threshold_value=np.percentile(non_detected_dw_mask[np.where(non_detected_dw_mask>-999)],percent_detected_dw)
        ind_inf_threshold=np.where(non_detected_dw_mask<=threshold_value)
        non_detected_dw_mask[np.where(non_detected_dw_mask>threshold_value)]=1
        non_detected_dw_mask[ind_inf_threshold]=0
        mask_dw=mask_dw+non_detected_dw_mask
        non_detected_dw_mask = None 
        mask_regions= None

    return mask_dw
