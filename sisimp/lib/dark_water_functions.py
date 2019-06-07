# -*- coding: utf-8 -*-
'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''

import numpy as np
import random
from skimage.measure import label
import lib.height_model as height_model
import lib.my_api as my_api


def dark_water_simulation(dlat, latmin, latmax, dlon, lonmin, lonmax, pourcent_dw, seedvalue, lcorr = 50):
    """Simulates dark water areas over water areas using a 2d gaussian profile
        Arguments : 
            dlat (float) : step between two lines of the matrix 
            latmin (float) : minimum line/latitude index of the water area
            latmax (float): maximum line/latitude index of the water area
            dlon (float) : step between two columns of the matrix
            lonmin (float) : minimum column/longitude index of the water area
            lonmax (float): maximum column/longitude index of the water area
            pourcent_dw (float) : percentage of dark_water to add in the water area 
            seedvalue (float): seedvalue for reproducible simulation
            lcorr (float) : correlation length of dark water region
        Returns a binary mask of dark water (1 = dark water pixels)

    """
    # Create a 2D gaussian profile
    profile_2d = height_model.generate_2d_profile_gaussian(dlat, latmin, latmax, dlon, lonmin, lonmax, 1, plot = False, lcorr = lcorr, seed = seedvalue)
    
    ### Create Dark water mask
    # Define the threshold value to keep the desired % of dark water
    threshold_value = np.percentile(profile_2d, 100-pourcent_dw)
    # Initiate the dw_mask with zeros and set to 1 the dark water values
    mask_dw = np.zeros(profile_2d.shape)
    mask_dw[np.where(profile_2d > threshold_value)] = 1
    
    return mask_dw




def dark_water_non_detected_simulation(mask_dw, dlat, latmin, latmax, dlon, lonmin, lonmax,percent_detected_dw,seedvalue, scale_factor = 0.5):
    """Simulates non detected dark water over dark water areas using a 2d gaussian profile
        Arguments : 
            mask_dw (2d array): dark water mask
            dlat (float): step between two lines of the matrix
            latmin (float) : minimum line/latitude index of the dark water area
            latmax (float): maximum line/latitude index of the dark water area
            dlon (float) : step between two columns of the matrix
            lonmin (float) : minimum column/longitude index of the dark water area
            lonmax (float): maximum column/longitude index of the dark water area
            percent_detected_dw (float) : percentage of detected dark_water to keep in the dark water mask
            seedvalue (float): seedvalue for reproducible simulation
            scale factor (float) : scale factor of dark water region

    """
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
            size_y = maxy-miny+1
            size_x = maxx-minx+1

            #Simulate non detected dark water
            
            lcorr_ = np.int(scale_factor*(size_y+size_x)/2.)
            if lcorr_ == 0:
                lcorr_ = 1
                
            profile_2d = height_model.generate_2d_profile_gaussian(1, minx, maxx+1, 1, miny, maxy+1,1, plot = False, lcorr = lcorr_, seed = seedvalue)
            
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

