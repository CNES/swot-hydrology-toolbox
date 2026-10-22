# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2026/03/25:version initiale.
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: pixel_ice.py
   :synopsis: Compute ice_dyn_f 
    Created on 2026/03/25

.. moduleauthor:: Mélanie Trudel and Victoria Litalien - UDES 

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""

from datetime import datetime, timedelta
import joblib
import numpy as np
import pandas as pd

import cnes.common.lib.my_variables as my_var


def load_data(in_fic_model):
    """
    Load pre-trained classification model used to compute ice_dyn_f

    Parameters
    ----------
    in_fic_model : string
        Full path of the file storing the model
        
    Returns
    -------
    model : 
        Pre-trained classification model 
    """
    model = joblib.load(in_fic_model)
    return model


#######################################

    
def compute_by_pixel(model, inc, sig0, phase_noise_std, interferogram_qual, sig0_qual):
    """
    Compute ice_dyn_f for each pixels.

    Parameters
    ----------
    model : 
        Pre-trained classification model used to infer ice_dyn_f at the pixel level
    inc : 1D-array of float
        Incidence angle
    sig0 : 1D array of float
        sigma0
    phase_noise_std : 1D array of float
        Phase noise standard deviation
    interferogram_qual : 1D array of uint
        Flag that indicates the quality of the interferogram quantities 
    sig0_qual : 1D array of uint
        Quality flag for sigma0

    Returns
    -------
    out_ice_dyn_f : numpy.ndarray (int8)
        Array of ice_dyn_f classification for each pixel.
    """
    
    # Init output array
    out_ice_dyn_f = np.full(len(inc), my_var.FV_NETCDF['int8'], dtype=np.int8)
    
    if model is None:
        return out_ice_dyn_f

    # Transform sig0 into dB
    sig0 = np.where(sig0 <= 0, np.nan, sig0)
    sig0_db = 10 * np.log10(sig0)

    # Construct a dataframe with the 5 variables
    df = pd.DataFrame({
        'sig0_dB': sig0_db,
        'inc': inc,
        'phase_std': phase_noise_std,
        'qual_sig0': sig0_qual,
        'qual_class': interferogram_qual
    })
    
    mask_good = (
        (df['qual_sig0'] < 262144) &
        (df['qual_class'] < 262144) &
        (df['inc'] > 0.5) &
        (df['inc'] < 4.3) &
        (np.isfinite(df['sig0_dB']))
    )

    features = ['sig0_dB', 'phase_std', 'inc']
    X = df[features]

    # Making the prediction
    y_pred = model.predict(X)
    ice_dyn_f = y_pred.copy()

    # Assign output values
    out_ice_dyn_f[mask_good] = ice_dyn_f[mask_good]

    return out_ice_dyn_f


#######################################


def threshold_ice_dyn_f(p_0, p_1, p_2, std_sig0):
    """
    Compute ice_dyn_f at the lake level wrt different thresholds of input parameters

    Parameters
    ----------
    p_0 : float
        Proportion of lake pixels with ice_dyn_f=0
    p_1 : float
        Proportion of lake pixels with ice_dyn_f=1
    p_2 : float
        Proportion of lake pixels with ice_dyn_f=2
    std_sig0 : float
        Standard deviation of sigma0 of the lake pixels

    Returns
    -------
    ice_dyn_f : int8
       ice_dyn_f value.
    """
    
    # Init default output value
    ice_dyn_f = my_var.FV_SHP['int8']  # Default output value
    
    if p_1 < 0.8: 
        
        if p_2 <= 0.30:
            
            if p_0 <= 0.20:
                if std_sig0 <= 7.5:
                    ice_dyn_f = 1
                else:
                    ice_dyn_f = 0
                    
            else: 
                ice_dyn_f = 0
                
        else:  # p_2 > 30%
            
            if std_sig0 <= 3.5:
                if p_0 <= 0.55:
                    ice_dyn_f = 2
                else: 
                    ice_dyn_f = 0
                    
            else: # std_sig0 > 3.5
                if std_sig0 <= 4.5:
                    ice_dyn_f = 2
                else:
                    ice_dyn_f = 1
                    
    return ice_dyn_f

    
def compute_ice_dyn_f_4lake(pix_lake, sig0):
    """
    Compute ice_dyn_f for Obs lake

    Parameters
    ----------
    pix_lake : array-like
        Table of pixel classes (0, 1, 2)
    sig0 : 1D array of float
        sigma0

    Returns
    -------
    int
        Value of ice_dyn_f (0, 1, 2, or -999)
    """

    pix_lake = np.asarray(pix_lake)

    # Total number of valid pixels for a lake
    n = len(pix_lake) - np.sum(pix_lake == my_var.FV_NETCDF['int8'])
    if (n < 5) or (not np.any(pix_lake < 5)):
        return my_var.FV_SHP['int8']

    # Transform sig0 into dB
    sig0 = np.where(sig0 <= 0, np.nan, sig0)
    sig0_db = 10 * np.log10(sig0)
    
    mask_good = (pix_lake != my_var.FV_NETCDF['int8'])
    std_sig0 = np.nanstd(sig0_db[mask_good])
    if std_sig0 == np.nan:
        std_sig0 = 0.0
    
    # Calculate proportions
    p_0 = np.sum(pix_lake == 0) / n
    p_1 = np.sum(pix_lake == 1) / n
    p_2 = np.sum(pix_lake == 2) / n
    
    # Applying thresholds
    ice_dyn_f = threshold_ice_dyn_f(p_0, p_1, p_2, std_sig0)
       
    return ice_dyn_f

    
def compute_ice_dyn_f_4prior(pix_lake, sig0, ice_start, ice_end, date_acquisition, time_window=20):
    
    """
    Compute ice_dyn_f for Prior lake

    Parameters
    ----------
    pix_lake : array-like
        Table of pixel classes (0, 1, 2)
    sig0 : 1D array of float
        sigma0
    ice_start : str (format 'MMDD')
    ice_end : str (format 'MMDD')
    date_acquisition : str (format 'YYYY-MM-DD')
    time_window : int
        Time window (in days) around the climatological ice cover period (ie from t-N days until t+N days) to consider ice_dyn_f=2 possible

    Returns
    -------
    int
        Value of ice_dyn_f (0, 1, 2, or -999)
    """

    pix_lake = np.asarray(pix_lake)

    # Total number of pixels for a lake
    n = len(pix_lake) - np.sum(pix_lake == my_var.FV_NETCDF['int8'])
    if (n < 5) or (not np.any(pix_lake < 5)):
        return my_var.FV_SHP['int8']
        
    # Acquisition date
    dt_acq = datetime.strptime(date_acquisition[:10], "%Y-%m-%d")
    year = dt_acq.year
        
    # Case 1 : ice_start and ice_end in the same year
    if ice_end >= ice_start:
        start = datetime.strptime(f"{year}{ice_start}", "%Y%m%d")
        end   = datetime.strptime(f"{year}{ice_end}", "%Y%m%d")

    # Case 2 : ice_start and ice_end not in the same year (ex: Nov 2024 to April 2025)
    else:
        # If acquisition is at the end of the year (nov/dec)
        if dt_acq.strftime("%m%d") >= ice_start:
            start = datetime.strptime(f"{year}{ice_start}", "%Y%m%d")
            end   = datetime.strptime(f"{year+1}{ice_end}", "%Y%m%d")
        else:
            # If acquisition is at the begining of the year (jan/fev)
            start = datetime.strptime(f"{year-1}{ice_start}", "%Y%m%d")
            end   = datetime.strptime(f"{year}{ice_end}", "%Y%m%d")
            
    # Apply window around the climatological ice cover period
    start_window = start - timedelta(days=time_window)
    end_window = end + timedelta(days=time_window)
      
    # Transform sig0 into dB
    sig0 = np.where(sig0 <= 0, np.nan, sig0)
    sig0_db = 10 * np.log10(sig0)
    
    mask_good = (pix_lake != my_var.FV_NETCDF['int8'])
    std_sig0 = np.nanstd(sig0_db[mask_good])
    if std_sig0 == np.nan:
        std_sig0 = 0.0
    
    # Calculate proportions
    p_0 = np.sum(pix_lake == 0) / n
    p_1 = np.sum(pix_lake == 1) / n
    p_2 = np.sum(pix_lake == 2) / n
    
    # Applying thresholds
    ice_dyn_f = threshold_ice_dyn_f(p_0, p_1, p_2, std_sig0)
        
    # Temporal filter
    if ice_dyn_f == 2:
        if (dt_acq < start_window) or (dt_acq > end_window):
            ice_dyn_f = my_var.FV_SHP['int8']
            
    # LOG FILE 
    #with open(log_file, "a") as f:
        #f.write(
            #f"{date_acquisition},"
            #f"{ice_start},{ice_end},"
            #f"{start_window.strftime('%Y-%m-%d')},"
            #f"{end_window.strftime('%Y-%m-%d')},"
            #f"{p_0:.3f},{p_1:.3f},{p_ice:.3f},"
            #f"{ice_dyn_f}\n"
        #)

    return ice_dyn_f