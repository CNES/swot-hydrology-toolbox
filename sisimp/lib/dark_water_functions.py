# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 14:35:20 2018

@author: plousser
"""

import numpy as np
import random
from skimage.measure import label

import lib.height_model as height_model
import lib.my_api as my_api


# To move in other file and import
def dark_water_simulation(taille_2D, fact_echelle_dw, pourcent_dw, seedvalue=None):
    
    # Create a 2D gaussian profile
    if seedvalue is not None:
        np.random.seed(int(seedvalue))
    profile_2d = height_model.generate_2d_profile_gaussian(taille_2D, 0., 'Default', 1, fact_echelle_dw)
    
    ### Create Dark water mask
    # Define the threshold value to keep the desired % of dark water
    threshold_value = np.percentile(profile_2d, 100-pourcent_dw)
    # Initiate the dw_mask with zeros and set to 1 the dark water values
    mask_dw = np.zeros(profile_2d.shape)
    mask_dw[np.where(profile_2d > threshold_value)] = 1
    
    return mask_dw
    

def dark_water_randomerase(mask_dw, water_pixmatrix, taille_2D, seedvalue=None):
    
    # Label dark_water regions
    mask_regions=label(mask_dw)

    # Reshape to water extent
    mask_regions=mask_regions[taille_2D]
    mask_dw=mask_dw[taille_2D]
    
    ### Delete random regions
    #Create a vector with each region label
    region_list=[i for i in range(1,np.max(mask_regions))]
            
    # Randomly select the regions to be deleted. The number of deleted regions is also random
    if seedvalue is not None:
        random.seed(seedvalue)
    region_todelete = random.sample(region_list,np.random.choice(region_list))
    my_api.printDebug("Dark water non detected pixels : %d" % len(np.where(np.isin(mask_regions, region_todelete))[0]))

    # Erase non detected dark_water pixels
    water_pixmatrix[np.where(np.isin(mask_regions, region_todelete))]=0
    mask_dw = np.delete(mask_dw,np.where(np.isin(mask_regions, region_todelete)))
    #mask_regions=np.delete(mask_regions,np.where(np.isin(mask_regions, region_todelete)))
    
    return mask_dw, water_pixmatrix
