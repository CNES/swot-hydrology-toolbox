#!/usr/bin/env python
"""
module write_polygons.py

module author : Capgemini

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from __future__ import absolute_import, division, print_function, unicode_literals

from osgeo import gdal, ogr
from osgeo.gdalconst import GDT_Float32

import numpy as np
import sys
import os
import utm 
import pyproj
from scipy import ndimage
from scipy.spatial import cKDTree
import time

import lib.dark_water_functions as dark_water
import lib.my_api as my_api
import lib.my_shp as my_shp
import lib.my_tools as my_tools
from lib.my_lacs import Constant_Lac, Reference_height_Lac, Gaussian_Lac, Polynomial_Lac, Height_in_file_Lac
from lib.my_variables import RAD2DEG, DEG2RAD, GEN_APPROX_RAD_EARTH
from lib.roll_module import Roll_module
from lib.tropo_module import Tropo_module

from inversion_algo import inversionCore

import mathematical_function as math_fct
import proc_real_pixc
import proc_real_pixc_vec_river

import pickle

GEN_RAD_EARTH = 6378137.0
GEN_RAD_EARTH_POLE = 6356752.31425

class orbitAttributes:
    
    def __init__(self):
        
        # 1 - Attributes from parameters file or default
    
        # 1.1 - Files and directories
        self.out_dir = None  # Output directory
        self.shapefile_path = None  # Full path of shapefile of input water bodies
    
        # 1.2 - Instrument parameters
        self.swath_width = None  # Swath width
        self.nr_cross_track = None  # NR cross track
        self.sensor_wavelength = None  # Sensor wavelength
        self.nb_pix_range = None  # Number of pixels in range
        self.range_sampling = None  # Range sampling
    
        # 1.3 - Orbit parameters
        self.orbit_jitter = None  # Orbit jitter
        self.mission_start_time = None  # Launch time
        self.cycle_duration = None  # Cycle duration
        self.multi_orbit_option = None  # Multiple orbit option
        self.orbit_number = None  # Orbit number to process
        self.passplan_path = None  # Full path of passplan.txt file
    
        # 1.4 - Height parameters
        self.trueheight_file = None  # True height file
        self.height_model = None  # Height model
        self.height_name = None  # Height field name
        self.height_model_min_area = None  # Optionnal argument to add complex 2D height model
        # Constant height model (always applied, set HEIGHT_MODEL_A to 0 to desactivate)
        self.height_model_a = None  # Height model A
        self.height_model_a_tab = None  # Height model from shapefile
        self.height_model_t0 = None  # Height model t0
        self.height_model_period = None  # Height model period
        # Gaussian parameter for gaussian model
        self.height_model_stdv = None  # Height model standard deviation (for gaussian model)
        
        # Dark water
        self.dark_water = None
        self.fact_echelle_dw = None
        self.dw_pourcent = None
        self.dw_seed = None
        self.darkwater_flag = None
        self.scale_factor_non_detected_dw = None
        self.dw_detected_percent = None
        self.dw_detected_noise_factor = None

        # 1.5 - Water flag
        self.water_flag = None
    
        # 1.6 - Noise parameters
        self.geolocalisation_improvement = None  # No noise applied on geolocation
        self.noise_multiplier_factor = None  # Noise multiplier factor
        self.height_bias_std = None  # Height bias std
    
        # 1.7 - Cross-over residual roll error
        self.roll_file = None  # Full path
    
        # 2 - Attributes for computation configuration
        self.create_shapefile = None
        self.create_pixc_vec_river = None
    
        # 3 - Working variables init
        self.sisimp_filenames = None  # Filenames specific to SISIMP
        self.noise_height = None  # Noise tab
        self.orbit_list = []  # List of triplets (cycle number, pass number, orbit file) to process
        self.compute_pixc_vec_river = None  # Flag for PIXCVecRiver file computation
        self.near_range = None
        self.tile_coords = {}  # Dictionnary for storing tile coordinates for swath R and L

        self.dw_detected_noise_height = None  # dw detected noise tab

        # 3.1 - From orbit file
        self.azimuth_spacing = None  # Azimuth spacing
        self.azimuth_layover = {} # Azimuth layover when using multi-tile
        self.lon = None
        self.lon_orbit = None
        self.lon_init = None
        self.lat = None
        self.lat_orbit = None
        self.lat_init = None
        self.alt = None
        self.heading = None
        self.heading_init = None
        self.orbit_time = None
        self.x = None
        self.y = None
        self.z = None
        self.cosphi_init = None
        self.sinphi_init = None
        self.costheta_init = None
        self.sintheta_init = None
        self.cospsi_init = None
        self.sinpsi_init = None

        # 4 - List of each water body
        self.liste_lacs = None

        # number of pixels overlapping top and bottom tile in tilling process
        self.nb_pix_overlap_begin = 0
        self.nb_pix_overlap_end = 0


#######################################
        

def compute_pixels_in_water(IN_fshp_reproj, IN_pixc_vec_only, IN_attributes):
    """
    Compute the position of the radar pixels that are inside a water body

    Method: use the GDAL rasterize function to mark the pixels that are
    contained in some polygons. This is very fast. Then we extract the pixels
    that are marked (value=1) and write them in an text file for check with qgis

    Important note: only polygons in one swath (left or right) must be in
                    fshp_reproj

    :param IN_fshp_reproj: full path of the shapefile with polygons in radar projection
    :type IN_fshp_reproj: string
    :param IN_pixc_vec_only: if set, deal only with polygons with field RIV_FLAG != 0
    :type IN_pixc_vec_only: boolean
    :param IN_attributes:
    :type IN_attributes:

    :return OUT_burn_data: the radar pixels that are inside a water body
    :rtype OUT_burn_data: 2D-array of int (0=land 1=water)
    :rettun OUT_height_data : the height pixels that are inside a water body
    :rtype OUT_height_data : 2D-array of int (height of each water body pixel)
    """
    if IN_pixc_vec_only:
        my_api.printInfo("[write_polygons] [compute_pixels_in_water] == compute_pixels_in_water / river polygons only ==")
    else:
        my_api.printInfo("[write_polygons] [compute_pixels_in_water] == compute_pixels_in_water / all polygons ==")

    # 1 - Read the reprojected shapefile
    layer, da_shapefile = my_shp.open_shp(IN_fshp_reproj)

    if IN_pixc_vec_only:
        layer.SetAttributeFilter(str("RIV_FLAG != '0'"))
        my_api.printInfo("[write_polygons] [compute_pixels_in_water] compute_pixels_in_water / river pixels only - %d features to deal with" % layer.GetFeatureCount())

    # Create a GDAL raster in memory
    nx = len(IN_attributes.lon)
    ny = IN_attributes.nb_pix_range
    
    cover = np.zeros((ny, nx), dtype='float64')
    cover_height = np.zeros((ny, nx), dtype='float64')
    cover_code = np.zeros((ny, nx), dtype='float64')
    cover_ind_lac = np.zeros((ny, nx), dtype='int')
    
    # Create 2 raster band if heights value come from shapefile
    #~ nb_rasterBand = [1, 2][IN_attributes.height_model == "reference_height" or IN_attributes.height_model == "gaussian" or IN_attributes.height_model == "polynomial"]
    nb_rasterBand = 4

    ds = gdal.GetDriverByName(str('MEM')).Create('', nx, ny, nb_rasterBand, GDT_Float32)
    ds.SetGeoTransform([-0.5, 1, 0, -0.5, 0, 1])
    ds.GetRasterBand(1).WriteArray(cover)
    
    # Burn with a value 1 the pixels that are in a polygon
    gdal.RasterizeLayer(ds, [1], layer, None, None, burn_values=[1])
    #~ #                                        , options=['ALL_TOUCHED=TRUE'])

    # Convert the radar pixels touching water in lon-lat
    OUT_burn_data = ds.GetRasterBand(1).ReadAsArray()
    
    # Analyse value in shapefile
    OUT_height_data = None
    OUT_code_data = None
    OUT_ind_lac_data = None
    
    #not used, commented to improve speed
    
    if IN_attributes.height_model == "reference_height":
        ds.GetRasterBand(2).WriteArray(cover_height)
        # Burn with height value pixels in the associated polygon
        gdal.RasterizeLayer(ds, [2], layer, None, options=["ATTRIBUTE=HEIGHT"])
        # Get height pixels in lon-lat
        OUT_height_data = ds.GetRasterBand(2).ReadAsArray()
    
       
    if IN_attributes.height_model == "gaussian" or IN_attributes.height_model == "polynomial":
        ds.GetRasterBand(3).WriteArray(cover_code)
        # Burn with height value pixels in the associated polygon
        gdal.RasterizeLayer(ds, [3], layer, None, options=["ATTRIBUTE=CODE"])
        # Get height pixels in lon-lat
        OUT_code_data = ds.GetRasterBand(3).ReadAsArray().astype(int)
 
    ds.GetRasterBand(4).WriteArray(cover_ind_lac)
    gdal.RasterizeLayer(ds, [4], layer, None, options=["ATTRIBUTE=IND_LAC"])
    OUT_ind_lac_data = ds.GetRasterBand(4).ReadAsArray().astype(int)

    # Close the raster
    ds = None

    if not IN_pixc_vec_only:
        my_api.printInfo("[write_polygons] [compute_pixels_in_water] Compute lakes labels")

        nb_lakes = len(IN_attributes.liste_lacs)

        if nb_lakes > 100 : # Quick way to compute labels for a large amount of lakes
            labels_coords = my_tools.coords_from_labels(OUT_ind_lac_data)

            for i, lake in enumerate(IN_attributes.liste_lacs):
                if lake.num in labels_coords:
                    lake.set_pixels_coods(np.array(labels_coords[lake.num]).transpose())
                else :
                    # lakes without pixels
                    lake.set_pixels_coods([[], []])
        else :
            for lake in IN_attributes.liste_lacs:
                lake.compute_pixels_in_given_lac(OUT_ind_lac_data)

    return OUT_burn_data, OUT_height_data, OUT_code_data, OUT_ind_lac_data, IN_attributes
                

#######################################
    

def write_water_pixels_realPixC(IN_water_pixels, IN_swath, IN_cycle_number, IN_orbit_number, IN_attributes, swath_polygon):
    """
    Check what pixels are marked as water and write their position.
    Real PixC files (Version 08/2018) are produced.
    Associated shapefiles are also produced if asked.
    
    :param IN_water_pixels: 2D-array of water pixels vs land pixels
    :type IN_water_pixels: 2D-array of int
    :param IN_swath: Left or Right swath
    :type IN_swath: string
    :param IN_cycle_number: cycle number
    :type IN_cycle_number: int
    :param IN_orbit_number: orbit number
    :type IN_orbit_number: int
    :param IN_attributes:
    :type IN_attributes:
    """
    my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] == write_water_pixels_realPixC ==")
    
    ################################
    # Variables for PixC main file #
    ################################
    
    # Print number of water pixels
    size_of_tabs = np.count_nonzero(IN_water_pixels) 
    my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] Nb water pixels: %d" %(size_of_tabs))

    # 1 - Get range and azimuth indices of all water pixels
    ind = np.nonzero(IN_water_pixels)  # Get indices 1=lake and 2=river (remove 0=land)
    r, az = [ind[0], ind[1]]  # Range index
    #az = ind[1]  # Azimuth index
    river_flag = IN_water_pixels[ind]  # 1=lake and 2=river
    # Dark Water
    
    # Check if dark water is simulated or not
    if IN_attributes.dark_water.lower() == "yes":
        # Simulate dark water if water pixel are present
        # TODO: integrer dans le if suivant après quand ça marche pour éviter les répétitions
        if size_of_tabs != 0.:
            rmin, rmax, azmin, azmax = r.min(), r.max(), az.min(), az.max()

            # Simulate dark_water
            dw_mask = dark_water.dark_water_simulation(1, azmin, azmax+1, 1, rmin, rmax+1, IN_attributes.dw_pourcent, IN_attributes.dw_seed, lcorr=IN_attributes.dw_correlation_length)
    
            # Get water extent
            indice_r = np.array(ind[0]-rmin)
            indice_az = np.array(ind[1]-azmin)
            # Randomly classify or erase dark water regions
            # Randomly erase DW regions in DW mask
            dw_mask = dark_water.dark_water_non_detected_simulation(dw_mask, 1, azmin, azmax+1, 1, rmin, rmax+1, IN_attributes.dw_detected_percent, IN_attributes.dw_seed, scale_factor=IN_attributes.scale_factor_non_detected_dw)
            # Reshape dark_water to water extent
            dw_mask = dw_mask[indice_az, indice_r]

            # Update IN_water_pixels with the deleted water pixels
            # Check if dw pixels are deleted = dw_mask pixels set to 2
            if np.where(dw_mask == 2)[0].size > 0:
                my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] Nb detected dark water pixels : %d" % np.where(dw_mask == 1)[0].size)
                my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] Nb non detected dark water pixels : %d" % np.where(dw_mask == 2)[0].size)
                # Update river_flag and IN_water_pixels with the deleted water pixels
                river_flag[np.where(dw_mask == 2)] = 0
                IN_water_pixels[ind] = river_flag
                # Delete the corresponding pixels in the dw mask to update indices values
                dw_mask = np.delete(dw_mask, np.where(dw_mask == 2))
                # Update size of tabs etc... because water pixels were deleted
                size_of_tabs = np.count_nonzero(IN_water_pixels)
                my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] Nb water pixels: %d" % size_of_tabs)
                ind = np.nonzero(IN_water_pixels)
                r, az = [ind[0], ind[1]]
                river_flag = IN_water_pixels[ind]
            # Build the classification array with water flag Dark_water flag
            classification_tab = np.ones(size_of_tabs) * IN_attributes.water_flag  # Classification as water
            # Locate DW pixels
            dark_water_loc = np.where(dw_mask == 1)
            # Update classification value for dark water pixels with DW flag
            classification_tab[dark_water_loc] = IN_attributes.darkwater_flag
        else:
            classification_tab = np.ones(size_of_tabs) * IN_attributes.water_flag
    else:
        classification_tab = np.ones(size_of_tabs) * IN_attributes.water_flag
        my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] No Dark Water will be simulated")

    if IN_attributes.height_model_a_tab is not None:
        height_flag = IN_attributes.height_model_a_tab[ind]
        
    if IN_attributes.code is not None:
        code_flag = IN_attributes.code[ind]
        
    # 2 - Compute radar variables for water pixels

    ri = IN_attributes.near_range + r * IN_attributes.range_sampling  # Radar-to-ground distance
    
    Hi = IN_attributes.alt[az]  # Altitude

    elevation_tab = np.zeros(len(az))
    water_frac = np.zeros(len(az))

    processing = np.round(np.linspace(0, len(IN_attributes.liste_lacs), 11), 0)
    for i, lac in enumerate(IN_attributes.liste_lacs):
        if i in processing:
            my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] Processing %d%%" % (int(100 * (i + 1) / len(IN_attributes.liste_lacs))))


        if lac.nb_pix > 0 :

            #~ indice = np.logical_and(np.isin(r, lac.pixels[0]), np.isin(az, lac.pixels[1]))
            def merge(a,b):
                return a*100000+b
            titi = merge(lac.pixels[0], lac.pixels[1])
            toto = merge(r, az)
            indice = np.isin(toto,titi)
            
            lon, lat = math_fct.lonlat_from_azy(az, ri, IN_attributes, IN_swath, IN_unit="deg", h=lac.hmean)
            elevation_tab[indice] = lac.compute_h(lat[indice], lon[indice])
            elevation_tab[indice] += lac.h_ref
            # TODO : Calcul water_frac for all pts
            #for i, val in enumerate(indice):
            #    if val:
            #        ring = ogr.Geometry(ogr.wkbLinearRing)
            #        ring.AddPoint( lon[i] + (IN_attributes.range_sampling / GEN_RAD_EARTH)*(180./np.pi)/np.cos(lat[i]*np.pi/180.), lat[i] + (IN_attributes.azimuth_spacing / GEN_RAD_EARTH)*(180./np.pi))
            #        ring.AddPoint( lon[i] - (IN_attributes.range_sampling / GEN_RAD_EARTH)*(180./np.pi)/np.cos(lat[i]*np.pi/180.), lat[i] + (IN_attributes.azimuth_spacing / GEN_RAD_EARTH)*(180./np.pi))
            #        ring.AddPoint( lon[i] + (IN_attributes.range_sampling / GEN_RAD_EARTH)*(180./np.pi)/np.cos(lat[i]*np.pi/180.), lat[i] - (IN_attributes.azimuth_spacing / GEN_RAD_EARTH)*(180./np.pi))
            #        ring.AddPoint( lon[i] - (IN_attributes.range_sampling / GEN_RAD_EARTH)*(180./np.pi)/np.cos(lat[i]*np.pi/180.), lat[i] - (IN_attributes.azimuth_spacing / GEN_RAD_EARTH)*(180./np.pi))
            #        ring.CloseRings()
            #    
            #        poly_point = ogr.Geometry(ogr.wkbPolygon)
            #        poly_point.AddGeometry(ring)
            #        poly_intersection = swath_polygon.Intersection(poly_point)
            #        if poly_intersection==None:
            #            water_frac[i]=1.0
            #        else:
            #            try:
            #                water_frac[i] = (poly_point.GetArea() - (poly_point.GetArea() - poly_intersection.GetArea()))/poly_point.GetArea()
            #            except ZeroDivisionError:
            #                water_frac[i] = poly_point.GetArea()

    # 3 - Build cross-track distance array
    # Compute theorical cross-track distance for water pixels
    # ~ sign = [-1, 1][IN_swath.lower() == 'right']
    # ~ y = sign * np.sqrt((ri + Hi - elevation_tab) * (ri - (Hi-elevation_tab)) / (1. + (Hi-elevation_tab) / GEN_APPROX_RAD_EARTH))
    # ~ angles = np.arccos((Hi - elevation_tab)/ri)  # Look angles
    
    angles = np.arccos((ri**2 + (GEN_APPROX_RAD_EARTH+Hi)**2 - (GEN_APPROX_RAD_EARTH+elevation_tab)**2)/(2*ri*(GEN_APPROX_RAD_EARTH+Hi)))
    theta = np.arccos(((GEN_APPROX_RAD_EARTH+Hi)**2 + (GEN_APPROX_RAD_EARTH+elevation_tab)**2 - ri**2)/(2*(GEN_APPROX_RAD_EARTH+elevation_tab)*(GEN_APPROX_RAD_EARTH+Hi)))
    y = theta*GEN_APPROX_RAD_EARTH
    
    if IN_swath.upper()[0] == 'L':
        y = -y
    
    
    pixel_area = IN_attributes.azimuth_spacing * IN_attributes.range_sampling / np.sin(angles)  # Pixel area
    
    ## 4.1 - Compute noise over height
    #if IN_attributes.dark_water.lower() == "yes":
    #    noise_seed = int(str(time.time()).split('.')[1])
    #    delta_h = np.zeros(elevation_tab.shape)
    #    phase_noise_std = np.zeros(elevation_tab.shape)
    #    dh_dphi = np.zeros(elevation_tab.shape)
    #    water_pixels = np.where(classification_tab == IN_attributes.water_flag)
    #    dw_pixels = np.where(classification_tab == IN_attributes.darkwater_flag)
    #    delta_h[water_pixels], phase_noise_std[water_pixels], dh_dphi[water_pixels] = math_fct.calc_delta_h(angles[water_pixels], IN_attributes.noise_height, IN_attributes.height_bias_std, IN_attributes.sensor_wavelength, IN_attributes.baseline, IN_attributes.near_range, seed=noise_seed)
    #    delta_h[dw_pixels], phase_noise_std[dw_pixels], dh_dphi[dw_pixels] = math_fct.calc_delta_h(angles[dw_pixels], IN_attributes.dw_detected_noise_height, IN_attributes.height_bias_std, IN_attributes.sensor_wavelength, IN_attributes.baseline, IN_attributes.near_range, seed=noise_seed)
    #else:
    #    delta_h, phase_noise_std, dh_dphi = math_fct.calc_delta_h(angles, IN_attributes.noise_height, IN_attributes.height_bias_std, IN_attributes.sensor_wavelength, IN_attributes.baseline, IN_attributes.near_range, seed=noise_seed)

    # 4.1 bis - Compute mean noise over points
    noise_seed = int(str(time.time()).split('.')[1])

    angles_pixels = np.zeros((len(IN_water_pixels), len(IN_water_pixels[0])))
    r_pixels = np.arange(0, len(IN_water_pixels)-1)
    ri_pixels = IN_attributes.near_range + r_pixels * IN_attributes.range_sampling
    Hi_pixels = IN_attributes.alt[np.arange(0, len(IN_water_pixels[0])-1)]
    for az_ind in range(len(IN_water_pixels[0])-1):
        lon_pixels, lat_pixels = math_fct.lonlat_from_azy(np.ones(len(ri_pixels), dtype=np.int)*i, ri_pixels, IN_attributes, IN_swath, IN_unit="deg", h=lac.hmean)
        elevation_az = lac.compute_h(lat_pixels, lon_pixels)
        Hi_pixels = IN_attributes.alt[i]
        angles_az = np.arccos((ri_pixels**2 + (GEN_APPROX_RAD_EARTH+Hi_pixels)**2 - (GEN_APPROX_RAD_EARTH+elevation_az)**2)/(2*ri_pixels*(GEN_APPROX_RAD_EARTH+Hi_pixels)))
        angles_az[np.isnan(angles_az)]=0.
        for r_ind in range(len(angles_az)-1):
            angles_pixels[r_ind][az_ind]=angles_az[r_ind]

    delta_h, phase_noise_std, dh_dphi = math_fct.calc_delta_h(IN_water_pixels, angles_pixels, angles, IN_attributes.noise_height, IN_attributes.height_bias_std, IN_attributes.sensor_wavelength, IN_attributes.baseline, IN_attributes.near_range, seed=noise_seed)

    conv_delta_h = ndimage.convolve(delta_h, np.array([[1/9, 1/9, 1/9], [1/9, 1/9, 1/9], [1/9, 1/9, 1/9]]))

    ind=np.where(IN_water_pixels!=0)
    delta_h = conv_delta_h[ind] # Keep only water points


    # 4.2 Add residual roll error
    try:
        if IN_attributes.roll_repo is not None:
            
            my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] Applying roll residual error")

            roll = Roll_module(IN_attributes.roll_repo)
            roll.get_roll_file_associated_to_orbit_and_cycle(IN_orbit_number, IN_cycle_number)
            roll.interpolate_roll_on_sensor_grid(IN_attributes.orbit_time)
            
            # Apply roll for each pixel
            roll.interpolate_roll_on_pixelcloud(IN_attributes.orbit_time, IN_attributes.orbit_time[az], y)
            delta_h_roll = roll.roll1_err_cloud
            
            delta_h += delta_h_roll

        else:
            my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] No roll error applied")

    except:
        my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] No roll error applied")
         
    # 4.3 Add tropospheric delay
    
    tropo = Tropo_module(IN_attributes.tropo_model, min(r), max(r), min(az), max(az), \
    IN_attributes.tropo_error_stdv, IN_attributes.tropo_error_mean, IN_attributes.tropo_error_correlation, \
    IN_attributes.tropo_error_map_file)
        
    if IN_attributes.tropo_map_rg_az is None:

        if tropo.model == 'gaussian':
            my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] Applying wet tropo gaussian model")
            tropo.calculate_tropo_error_gaussian(IN_attributes.tropo_error_stdv, IN_attributes.tropo_error_mean, IN_attributes.tropo_error_correlation) 
            tropo.apply_tropo_error_on_pixels(az, r)
            tropo_2d_field = tropo.tropo_2d_field
            delta_h += tropo_2d_field
            
        if tropo.model == 'map':
            my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] Applying wet tropo map gaussian model")
            tropo.calculate_tropo_error_map(np.mean(IN_attributes.lat), IN_attributes.tropo_error_map_file, IN_attributes.tropo_error_correlation)
            tropo.apply_tropo_error_on_pixels(az, r)
            tropo_2d_field = tropo.tropo_2d_field
            delta_h += tropo_2d_field

    else :
        tropo.tropo_map_rg_az = IN_attributes.tropo_map_rg_az
        tropo.apply_tropo_error_on_pixels(az, r)
        tropo_2d_field = tropo.tropo_2d_field
        delta_h += tropo_2d_field        
   
    
    # 5.3 - Compute final noisy heights (elevation + thermal noise + roll error + height model) 
    elevation_tab_noisy = elevation_tab + delta_h           
    
    # 7 - Build velocity arrays
    nb_pix_nadir = IN_attributes.x.size  # Nb pixels at nadir
    # Init velocity arrays
    vx = np.zeros(nb_pix_nadir)
    vy = np.zeros(nb_pix_nadir)
    vz = np.zeros(nb_pix_nadir)
    # Compute first value
    vx[0], vy[0], vz[0] = [(IN_attributes.x[1] - IN_attributes.x[0]), (IN_attributes.y[1] - IN_attributes.y[0]), (IN_attributes.z[1] - IN_attributes.z[0])] / (IN_attributes.orbit_time[1] - IN_attributes.orbit_time[0])
    # Compute last value
    vx[-1], vy[-1], vz[-1] = [(IN_attributes.x[-1] - IN_attributes.x[-2]), (IN_attributes.y[-1] - IN_attributes.y[-2]), (IN_attributes.z[-1] - IN_attributes.z[-2])] / (IN_attributes.orbit_time[-1] - IN_attributes.orbit_time[-2])
    # Compute middle values
    for indp in range(1, nb_pix_nadir-1):
        vx[indp], vy[indp], vz[indp] = [(IN_attributes.x[indp+1] - IN_attributes.x[indp-1]), (IN_attributes.y[indp+1] - IN_attributes.y[indp-1]), (IN_attributes.z[indp+1] - IN_attributes.z[indp-1])] / (IN_attributes.orbit_time[indp+1] - IN_attributes.orbit_time[indp-1])
 
    # Convert nadir latitudes/longitudes in degrees
    nadir_lat_deg = IN_attributes.lat * RAD2DEG
    nadir_lon_deg = IN_attributes.lon * RAD2DEG
    
    # Remove 1st and last values because just here for extrapolators (cf. read_orbit)
    nadir_alt = IN_attributes.alt
    nadir_heading = IN_attributes.heading

    lon_noisy, lat_noisy = math_fct.lonlat_from_azy(az, ri, IN_attributes, IN_swath, h=elevation_tab+delta_h)
    lon_noisy *= RAD2DEG  # Conversion in degrees
    lat_noisy *= RAD2DEG  # Conversion in degrees
    
    dlon_dphi = (lon_noisy-lon)/delta_h*dh_dphi
    dlat_dphi = (lat_noisy-lat)/delta_h*dh_dphi
    
    
    ######################
    # Write output files #
    ######################
    
    # Cut arrays in order to write real PixC files
    # Tiles correspond to theoretical tiles, with 60km length at nadir
    if IN_attributes.lat.size != 0:  
         
        my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] == Dealing with tile number %03d" % IN_attributes.tile_number)
                
        # Get azimuth indices corresponding to this integer value of latitude
        nadir_az = np.arange(1, len(IN_attributes.alt)-1)
                    
        az_min = np.sort(nadir_az)[0]  # Min azimuth index, to remove from tile azimuth ivector
        my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] = %d pixels in azimuth (index %d put to 0)" % (nadir_az.size, az_min))

        # Get overlapping pixels du to tilling process
        nb_pix_overlap_begin = IN_attributes.nb_pix_overlap_begin
        nb_pix_overlap_end = IN_attributes.nb_pix_overlap_end

        # Get pixel indices of water pixels corresponding to this latitude interval
        az_indices = np.where((az >= min(nadir_az) + nb_pix_overlap_begin -1 ) & (az <= max(nadir_az) - nb_pix_overlap_end +1 ))[0]
        
        nb_pix = az_indices.size  # Number of water pixels for this latitude interval
        my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] = %d water pixels" % nb_pix)
        
        if az_indices.size != 0:  # Write water pixels at this latitude
            
            sub_az, sub_r = [az[az_indices], r[az_indices]]
            
            my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] Min r ind = %d - Max r ind = %d" % (np.min(sub_r), np.max(sub_r)))
            my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] Min az ind = %d - Max az ind = %d" % (np.min(sub_az), np.max(sub_az)))

            # Get output filename

            # Left / right swath flag
            left_or_right = IN_swath.upper()[0]
            
            # General tile reference
            tile_ref = "%03d%s" % (IN_attributes.tile_number, left_or_right)

            sub_az = sub_az - az_min - nb_pix_overlap_begin + 1

            # remove first and last orbit point, added in read_orbit
            if nb_pix_overlap_begin == 0:
                nb_pix_overlap_begin = 1
                nb_pix_overlap_end = nb_pix_overlap_end - nb_pix_overlap_begin
            if nb_pix_overlap_end == 0:
                nb_pix_overlap_end = 1

            #TODO:keep only indices + negative values from azr_from_lonlat
            tile_nadir_lat_deg = nadir_lat_deg[nb_pix_overlap_begin:-nb_pix_overlap_end]
            tile_nadir_lon_deg = nadir_lon_deg[nb_pix_overlap_begin:-nb_pix_overlap_end]
            tile_orbit_time = IN_attributes.orbit_time[nb_pix_overlap_begin:-nb_pix_overlap_end]
            tile_nadir_alt = nadir_alt[nb_pix_overlap_begin:-nb_pix_overlap_end]
            tile_nadir_heading = nadir_heading[nb_pix_overlap_begin:-nb_pix_overlap_end]

            tile_x = IN_attributes.x[nb_pix_overlap_begin:-nb_pix_overlap_end]
            tile_y = IN_attributes.y[nb_pix_overlap_begin:-nb_pix_overlap_end]
            tile_z = IN_attributes.z[nb_pix_overlap_begin:-nb_pix_overlap_end]
            tile_vx = vx[nb_pix_overlap_begin:-nb_pix_overlap_end]
            tile_vy = vy[nb_pix_overlap_begin:-nb_pix_overlap_end]
            tile_vz = vz[nb_pix_overlap_begin:-nb_pix_overlap_end]
            
            # Calcul x, y, z of satellite on earth
            x_earth, y_earth, z_earth = inversionCore.convert_llh2ecef(tile_nadir_lat_deg, tile_nadir_lon_deg, np.zeros(len(tile_nadir_lat_deg)), GEN_RAD_EARTH, GEN_RAD_EARTH_POLE)
            nadir_vx = tile_x - x_earth
            nadir_vy = tile_y - y_earth
            nadir_vz = tile_z - z_earth

            # Calcul norm for each vector
            norm_sat = [np.sqrt(x*x + y*y + z*z) for x, y, z in zip(tile_vx, tile_vy, tile_vz)]
            norm_nadir = [np.sqrt(x*x + y*y + z*z) for x, y, z in zip(nadir_vx, nadir_vy, nadir_vz)]

            # Normalize vectors
            norm_vx = [val / norm_val for val, norm_val in zip(tile_vx, norm_sat)]
            norm_vy = [val / norm_val for val, norm_val in zip(tile_vy, norm_sat)]
            norm_vz = [val / norm_val for val, norm_val in zip(tile_vz, norm_sat)]

            norm_nadir_vx = [val/norm_val for val, norm_val in zip(nadir_vx, norm_nadir)]
            norm_nadir_vy = [val/norm_val for val, norm_val in zip(nadir_vy, norm_nadir)]
            norm_nadir_vz = [val/norm_val for val, norm_val in zip(nadir_vz, norm_nadir)]

            # Calcul normal vector
            normal_vector = [np.cross([vx, vy, vz], [nadir_vx, nadir_vy, nadir_vz]) for vx, vy, vz, nadir_vx, nadir_vy, nadir_vz in zip(norm_vx, norm_vy, norm_vz, norm_nadir_vx, norm_nadir_vy, norm_nadir_vz)]
           
            # Find sensors coordinates
            sensor_plus_y_x = [x + nv[0]*5 for x, nv in zip(tile_x, normal_vector)]
            sensor_plus_y_y = [y + nv[1]*5 for y, nv in zip(tile_y, normal_vector)]
            sensor_plus_y_z = [z + nv[2]*5 for z, nv in zip(tile_z, normal_vector)]

            sensor_minus_y_x = [x - nv[0]*5 for x, nv in zip(tile_x, normal_vector)]
            sensor_minus_y_y = [y - nv[1]*5 for y, nv in zip(tile_y, normal_vector)]
            sensor_minus_y_z = [z - nv[2]*5 for z, nv in zip(tile_z, normal_vector)]

            # Calcul interferogram for all water pixel 
            interf_2d = []
            x_water, y_water, z_water = inversionCore.convert_llh2ecef(lat_noisy[az_indices], lon_noisy[az_indices], elevation_tab_noisy[az_indices], GEN_RAD_EARTH, GEN_RAD_EARTH_POLE)
            for i, (x_w, y_w, z_w) in enumerate(zip(x_water, y_water, z_water), 0):
                interferogram = my_tools.compute_interferogram(sensor_plus_y_x[sub_az[i]], sensor_plus_y_y[sub_az[i]], sensor_plus_y_z[sub_az[i]],\
                    sensor_minus_y_x[sub_az[i]], sensor_minus_y_y[sub_az[i]], sensor_minus_y_z[sub_az[i]], x_w, y_w, z_w) 
                interf_2d.append([interferogram[0], interferogram[1]])

            #  Calcul water frac for each point 

            
            
            nadir_az_size = tile_nadir_lat_deg.size

            tile_coords = compute_tile_coords(IN_attributes, IN_swath, nadir_az_size, nb_pix_overlap_begin)
            IN_attributes.tile_coords[left_or_right] = tile_coords


            # Init L2_HR_PIXC object
            my_pixc = proc_real_pixc.l2_hr_pixc(sub_az, sub_r, classification_tab[az_indices], pixel_area[az_indices],
                                                lat_noisy[az_indices], lon_noisy[az_indices], elevation_tab_noisy[az_indices], phase_noise_std[az_indices],
                                                dh_dphi[az_indices], dlon_dphi[az_indices], dlat_dphi[az_indices], y[az_indices],
                                                tile_orbit_time, tile_nadir_lat_deg, tile_nadir_lon_deg, tile_nadir_alt, tile_nadir_heading,
                                                tile_x, tile_y, tile_z, tile_vx, tile_vy, tile_vz,
                                                IN_attributes.near_range, IN_attributes.mission_start_time, IN_attributes.cycle_duration, IN_cycle_number,
                                                IN_orbit_number, tile_ref, IN_attributes.nb_pix_range, nadir_az_size, IN_attributes.azimuth_spacing,
                                                IN_attributes.range_sampling, IN_attributes.near_range, tile_coords, interf_2d)
            
            # Update filenames with tile ref
            IN_attributes.sisimp_filenames.updateWithTileRef(tile_ref, IN_attributes.orbit_time[nadir_az[0]], IN_attributes.orbit_time[nadir_az[-1]])
            
            # Write main file
            my_pixc.write_pixc_file(IN_attributes.sisimp_filenames.pixc_file+".nc", compress=True)
            
            # Write annotation file
            my_pixc.write_annotation_file(IN_attributes.sisimp_filenames.file_annot_file, IN_attributes.sisimp_filenames.pixc_file+".nc")
            
            # Write shapefiles if asked
            if IN_attributes.create_shapefile:
                my_pixc.write_pixc_asShp(IN_attributes.sisimp_filenames.pixc_file+"_pixc.shp")
                my_pixc.write_tvp_asShp(IN_attributes.sisimp_filenames.pixc_file+"_tvp.shp")
                
            # Write PIXCVec files if asked
            if IN_attributes.create_pixc_vec_river:
                # Init PIXCVec product
                my_pixc_vec = proc_real_pixc_vec_river.l2_hr_pixc_vec_river(sub_az, sub_r, IN_attributes.mission_start_time, IN_attributes.cycle_duration, IN_cycle_number, IN_orbit_number, tile_ref, IN_attributes.nb_pix_range, nadir_az.size)
                # Set improved geoloc
                lon_no_noisy, lat_no_noisy = math_fct.lonlat_from_azy(az, ri, IN_attributes, IN_swath, h=elevation_tab)
                my_pixc_vec.set_vectorproc(lat_no_noisy[az_indices], lon_no_noisy[az_indices], elevation_tab[az_indices])
                # Compute river_flag
                my_pixc_vec.set_river_lake_tag(river_flag[az_indices]-1)  # -1 to have 0=lake and 1=river
                # Write PIXCVec file
                my_pixc_vec.write_file(IN_attributes.sisimp_filenames.pixc_vec_river_file+".nc", compress=True)
                # Write as shapefile if asked
                if IN_attributes.create_shapefile:
                    my_pixc_vec.write_file_asShp(IN_attributes.sisimp_filenames.pixc_vec_river_file+".shp")
            
    else:  
        my_api.printInfo("[write_polygons] [write_water_pixels_realPixC] No output data file to write")
                

#######################################


def reproject_shapefile(IN_filename, IN_swath, IN_driver, IN_attributes, IN_cycle_number):
    """
    Read the water polygon shapefile and compute polygons in radar coordinates. 
    Save the reprojected polygons in a new shapefile and return its name.

    Compute the part of water bodies polygons that is in the swath.
    This is needed to avoid folding of one swath on the other one due to
    left-range ambiguity caused by the transformation in azimuth-range.

    :param IN_filename: full path of the given shapefile
    :type IN_filename: string 
    :param IN_swath: the name of the swath 
    :type IN_swath: string ("Left" or "Right")
    :param IN_driver: OGR driver
    :type IN_driver: -
    :param IN_attributes:
    :type IN_attributes:
    :param IN_cycle_number:
    :type IN_cycle_number:
    
    :return OUT_file/Getname: full path of the shapefile with the water body polygons in radar coordinates
    :rtype OUT_filename: string
    :return OUT_swath_polygons
    :rtype OUT_swath_polygons
    """
    my_api.printInfo("[write_polygons] [reproject_shapefile] == reproject_shapefile ==")

    # 1 - Make swath polygons
    swath_polygon = make_swath_polygon(IN_swath, IN_attributes)

    # 2 - Select water bodies polygons in the swath
    layer, da_shapefile = my_shp.open_shp(IN_filename, swath_polygon)

    # Compute near_range 
    IN_attributes.near_range, hmean = compute_near_range(IN_attributes, layer, cycle_number = IN_cycle_number)
    layer.ResetReading()
    
    ## LOOP if hmean different from 0: 
    if hmean != 0:
        # 1 - Make swath polygons
        swath_polygon = make_swath_polygon(IN_swath, IN_attributes, hmean=hmean)
        # 2 - Select water bodies polygons in the swath
        layer, da_shapefile = my_shp.open_shp(IN_filename, swath_polygon)
        # Compute near_range (assuming height is nul, should be done differrently to handle totpography)
        IN_attributes.near_range, hmean = compute_near_range(IN_attributes, layer, cycle_number = IN_cycle_number)
        layer.ResetReading()    
        ## LOOP
    
    nb_features = layer.GetFeatureCount()

    my_api.printInfo("[write_polygons] [reproject_shapefile] There are %d feature(s) crossing %s swath" % (nb_features, IN_swath))
    sys.stdout.flush()

    # Exit if no feature to deal with
    if nb_features == 0:
        return None, IN_attributes, None

    # 3 - Create the output shapefile
    swath_t = '%s_swath' % IN_swath
    OUT_filename = os.path.join(IN_attributes.out_dir, os.path.splitext(os.path.split(IN_filename)[1])[0] + '_tmp_radarproj_%s.shp' % swath_t)
    if os.path.exists(OUT_filename):
        IN_driver.DeleteDataSource(OUT_filename)
    dataout = IN_driver.CreateDataSource(OUT_filename)
    if dataout is None:
        my_api.printError("[write_polygons] [reproject_shapefile] Could not create file")
    layerout = dataout.CreateLayer(str(os.path.splitext(os.path.split(OUT_filename)[1])[0]), None, geom_type=ogr.wkbPolygon)
    # Create necessary output fields 
    layerout.CreateField(ogr.FieldDefn(str('RIV_FLAG'), ogr.OFTInteger))
    layerout.CreateField(ogr.FieldDefn(str('HEIGHT'), ogr.OFTReal))
    layerout.CreateField(ogr.FieldDefn(str('CODE'), ogr.OFTInteger64))
    layerout.CreateField(ogr.FieldDefn(str('IND_LAC'), ogr.OFTInteger64))

    floutDefn = layerout.GetLayerDefn()
    feature_out = ogr.Feature(floutDefn)

    # 4 - Convert coordinates for each water body polygon
    range_tab = []

    liste_lac = []
    indmax = 0

    nb_lakes = layer.GetFeatureCount()
    processing = np.round(np.linspace(0, nb_lakes, 11), 0)

    
    for ind, polygon_index in enumerate(layer):
        if ind in processing :
            my_api.printInfo("[write_polygons] [reproject_shapefile] Processing %d %%" %( int(100 * (ind+1)/nb_lakes) ) )

        geom = polygon_index.GetGeometryRef()
        
        if geom is not None:  # Test geom.IsValid() not necessary

            # 4.1 - Fill RIV_FLAG flag
            if IN_attributes.compute_pixc_vec_river:
                riv_flag = polygon_index.GetField(str("RIV_FLAG"))
            else:
                riv_flag = 0

            # 4.2 - Create the output polygon in radar projection
            # 4.2.1 - Init output geometry
            geom_out = ogr.Geometry(ogr.wkbPolygon)

            # 4.2.2 - Compute the zone resulting of the intersection between polygon and swath
            intersection = geom.Intersection(swath_polygon)
            
            save_field = False
            try:
                id_lake = polygon_index.GetField("code")
            except ValueError:
                id_lake = polygon_index.GetField("id")
            # Test area of intersect zones

            if intersection is None or id_lake!=25 or intersection.GetArea() != geom.GetArea():
                save_field = True
                intersection=geom
            
# ~ =======
            # ~ if not intersection:
                # ~ geom_buff = geom.Clone().Buffer(0)
                # ~ intersection =  geom_buff.Intersection(swath_polygon)
# ~ >>>>>>> develop
            # 4.2.3 - Convert polygons coordinates
            add_ring = False

            # Oversample the polygon to correctly reconstruct it in SAR geometry
            intersection.Segmentize(0.01)

            for ring in all_linear_rings(intersection):
                npoints = ring.GetPointCount()

                if npoints == 0:
                    continue  # ignore polygons completely outside the swath

                points = np.transpose(np.array(ring.GetPoints()))

                lon, lat = points[0], points[1]

                x_c, y_c, zone_number, zone_letter = utm.from_latlon(lat[0], lon[0])
                latlon = pyproj.Proj(init="epsg:4326")
                utm_proj = pyproj.Proj("+proj=utm +zone={}{} +ellps=WGS84 +datum=WGS84 +units=m +no_defs".format(zone_number, zone_letter))
                X, Y = pyproj.transform(latlon, utm_proj, lon, lat)

                ring_xy = ogr.Geometry(ogr.wkbLinearRing)
                for i in range(len(X)):
                    ring_xy.AddPoint(X[i], Y[i])
                poly_xy = ogr.Geometry(ogr.wkbPolygon)
                poly_xy.AddGeometry(ring_xy)

                # area in ha
                area = poly_xy.GetArea()/10000.

                layerDefn = layer.GetLayerDefn()

                lon = lon * DEG2RAD
                lat = lat * DEG2RAD


                if IN_attributes.height_model is None:
                    lac = Constant_Lac(ind+1, IN_attributes, lat, IN_cycle_number, id_lake)

                if IN_attributes.height_model == 'reference_height':
                    lac = Reference_height_Lac(ind+1, polygon_index, layerDefn, IN_attributes, id_lake)

                if IN_attributes.height_model == 'gaussian':
                    if area > IN_attributes.height_model_min_area:
                        my_api.printInfo(str("Gaussian model applied for big water body of size %f ha" % area))
                        lac = Gaussian_Lac(ind+1, IN_attributes, lat * RAD2DEG, lon * RAD2DEG, IN_cycle_number, id_lake)
                    else:
                        lac = Constant_Lac(ind+1, IN_attributes, lat* RAD2DEG, IN_cycle_number, id_lake)

                if IN_attributes.height_model == 'polynomial':
                    if area > IN_attributes.height_model_min_area:
                        my_api.printInfo(str("Polynomial model applied for big water body of size %f ha" % area))
                        lac = Polynomial_Lac(ind+1, IN_attributes, lat* RAD2DEG, lon* RAD2DEG, IN_cycle_number, id_lake)
                    else:
                        lac = Constant_Lac(ind+1, IN_attributes, lat* RAD2DEG, IN_cycle_number, id_lake)

                if IN_attributes.height_model == "reference_file":
                    if IN_attributes.trueheight_file is not None:
                        lac = Height_in_file_Lac(ind+1, IN_attributes, id_lake)
                    else:
                        lac = Constant_Lac(ind+1, IN_attributes, lat, IN_cycle_number, id_lake)

                lac.set_hmean(np.mean(lac.compute_h(lat* RAD2DEG, lon* RAD2DEG)))
                try:
                    lac.h_ref=float(polygon_index.GetField("ref_height"))
                    if lac.h_ref==None:
                        lac.h_ref=0
                except:
                    lac.h_ref=0

                if IN_attributes.height_model == 'polynomial' and area > IN_attributes.height_model_min_area:
                    # Create database to save all parameters 
                    if save_field:
                        # Read DB
                        try:
                            lake_polynomial_file=open('lake_polynomial_param.pkl', 'rb')
                            polynomial_param_db=pickle.load(lake_polynomial_file)
                            try:
                                # Get saved parameters for this lake
                                lac.X0 = polynomial_param_db[id_lake][0].copy()
                                lac.Y0 = polynomial_param_db[id_lake][1].copy()
                                lac.COEFF_X2 = polynomial_param_db[id_lake][2].copy()
                                lac.COEFF_Y2 = polynomial_param_db[id_lake][3].copy()
                                lac.COEFF_X = polynomial_param_db[id_lake][4].copy()
                                lac.COEFF_Y = polynomial_param_db[id_lake][5].copy()
                                lac.COEFF_XY = polynomial_param_db[id_lake][6].copy()
                                lac.COEFF_CST = polynomial_param_db[id_lake][7].copy()
                            except KeyError:
                                # If lake seen for first time - no layover
                                polynomial_param_db[id_lake] = [lac.X0, lac.Y0, lac.COEFF_X2, lac.COEFF_Y2, lac.COEFF_X, lac.COEFF_Y, lac.COEFF_XY, lac.COEFF_CST] 
                            lake_polynomial_file.close() 
                            lake_polynomial_file=open(os.path.join(IN_attributes.out_dir,'lake_polynomial_param.pkl'), 'wb')
                            # Add azimuth layover of lake's current tile layover
                            pickle.dump(polynomial_param_db, lake_polynomial_file, pickle.HIGHEST_PROTOCOL)
                            lake_polynomial_file.close()
                        # Create DB if not existent    
                        except:
                            lake_polynomial_file=open(os.path.join(IN_attributes.out_dir,'lake_polynomial_param.pkl'), 'wb')
                            pickle.dump({id_lake:[lac.X0, lac.Y0, lac.COEFF_X2, lac.COEFF_Y2, lac.COEFF_X, lac.COEFF_Y, lac.COEFF_XY, lac.COEFF_CST]}, lake_polynomial_file, pickle.HIGHEST_PROTOCOL)
                            lake_polynomial_file.close()
                    # Reset h_mean of the lake
                    lac.set_hmean(np.mean(lac.compute_h(lat*RAD2DEG, lon*RAD2DEG)))
                    az, r = azr_from_lonlat(lon, lat, IN_attributes, heau=lac.compute_h(lat* RAD2DEG, lon* RAD2DEG))

                elif IN_attributes.height_model=='gaussian' and area > IN_attributes.height_model_min_area:
                    # Create DB to save gaussian parameters
                    try:
                        lake_gaussian_file=open(os.path.join(IN_attributes.out_dir,'lake_gaussian_param.pkl'), 'rb')
                        gaussian_param_db=pickle.load(lake_gaussian_file)
                        try:
                            lac.h_interp=gaussian_param_db[id_lake]
                        except KeyError:
                            gaussian_param_db[id_lake]=lac.h_interp
                        lake_gaussian_file.close()
                        lake_gaussian_file=open(os.path.join(IN_attributes.out_dir,'lake_gaussian_param.pkl'),'wb')
                        pickle.dump(gaussian_param_db, lake_gaussian_file, pickle.HIGHEST_PROTOCOL)
                        lake_gaussian_file.close()
                    except:
                        lake_gaussian_file=open(os.path.join(IN_attributes.out_dir,'lake_gaussian_param.pkl'), 'wb')
                        pickle.dump({id_lake:lac.h_interp}, lake_gaussian_file, pickle.HIGHEST_PROTOCOL)
                        lake_gaussian_file.close()
                    lac.set_hmean(np.mean(lac.compute_h(lat*RAD2DEG, lon*RAD2DEG)))
                    az, r = azr_from_lonlat(lon, lat, IN_attributes, heau=lac.hmean)

                else:
                    az, r = azr_from_lonlat(lon, lat, IN_attributes, heau=lac.hmean)

                range_tab = np.concatenate((range_tab, r), -1)
                npoints = len(az)
                if len(az) != len(lon):
                    my_api.printDebug("[write_polygons] [reproject_shapefile] Ignore polygons crossing the swath")
                    continue  # Ignore polygons crossing the swath

                for p in range(npoints):  # no fonction ring.SetPoints()
                    ring.SetPoint(p, az[p], r[p])
                #ring.Intersection(swath_polygon)
                ring.CloseRings()

                add_ring = True
                # Add the reprojected ring to the output geometry
                geom_out.AddGeometry(ring)
            
            #geom_out=geom_out.Intersection(swath_polygon)

            # 4.2.4 - Add Output geometry
            if add_ring:

                # Add the output reprojected polygon to the output feature
                feature_out.SetGeometry(geom_out)
                # Set the RIV_FLAG field
                feature_out.SetField(str("RIV_FLAG"), riv_flag)
                # Set the HEIGHT field
                height_from_shp = False
                layerDefn = layer.GetLayerDefn()
                for i in range(layerDefn.GetFieldCount()):
                    if IN_attributes.height_model == 'polynomial' or IN_attributes.height_model == 'gaussian':
                        try:
                            feature_out.SetField(str("CODE"), polygon_index.GetField("code"))
                        except ValueError:
                            feature_out.SetField(str("CODE"), polygon_index.GetField("id"))
                if not height_from_shp:
                    IN_attributes.height_model_a_tab = None

                feature_out.SetField(str("IND_LAC"), ind+1)

                #~ # Not used, commented to improve speed
                if IN_attributes.height_model == 'reference_height':
                    feature_out.SetField(str("HEIGHT"), lac.height)

                # Add the output feature to the output layer
                layerout.CreateFeature(feature_out)
                liste_lac.append(lac)
                indmax += 1

    dataout.Destroy()

    safe_flag_layover = True
    if (safe_flag_layover) and (IN_attributes.height_model == 'reference_height'):
        liste_lac = intersect(OUT_filename, OUT_filename, indmax, IN_attributes) + liste_lac

    IN_attributes.liste_lacs = liste_lac

    return OUT_filename, IN_attributes, swath_polygon
                

#######################################
        
        
def make_swath_polygon(IN_swath, IN_attributes, hmean=0):
    """
    Make left of right swath polygon
    
    :param IN_swath 
    :type IN_swath
    :param IN_attributes:
    :type IN_attributes:
    """
    n = len(IN_attributes.lon_init) - 4
    az = np.arange(2, n-2 , 1)

    # # Old computing
    # sign = [-1, 1][IN_swath.lower() == 'right']
    # ymin = sign * IN_attributes.nr_cross_track
    # ymax = sign * IN_attributes.swath_width/2
    # y = ymin * np.ones(len(az))
    # lon1_, lat1_ = math_fct.lonlat_from_azy_old(az, y, IN_attributes.lat_init, IN_attributes.lon_init, IN_attributes.heading_init)
    # y = ymax * np.ones(len(az))
    # lon2_, lat2_ = math_fct.lonlat_from_azy_old(az, y, IN_attributes.lat_init, IN_attributes.lon_init, IN_attributes.heading_init)

    # ~ IN_ri = np.sqrt((IN_attributes.alt[az] + (IN_attributes.nr_cross_track ** 2) / (2 * GEN_APPROX_RAD_EARTH)) ** 2 + IN_attributes.nr_cross_track ** 2)
    
    theta = IN_attributes.nr_cross_track/(GEN_APPROX_RAD_EARTH+hmean)
    IN_ri = np.sqrt((GEN_APPROX_RAD_EARTH+IN_attributes.alt[az])**2+(GEN_APPROX_RAD_EARTH+hmean)**2-2*np.cos(theta)*(GEN_APPROX_RAD_EARTH+hmean)*(GEN_APPROX_RAD_EARTH+IN_attributes.alt[az]))
    
    lon1, lat1 = math_fct.lonlat_from_azy(az, IN_ri, IN_attributes, IN_swath, h=hmean, IN_unit="deg")
    lon2, lat2 = math_fct.lonlat_from_azy(az, IN_ri + IN_attributes.nb_pix_range*IN_attributes.range_sampling, IN_attributes, IN_swath, h=hmean, IN_unit="deg")
    
    lonswath = np.concatenate((lon1, lon2[::-1]))
    latswath = np.concatenate((lat1, lat2[::-1]))

    swath_polygon = ogr.Geometry(ogr.wkbPolygon)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    n = len(lonswath)
    for i in range(n):
        ring.AddPoint(lonswath[i], latswath[i])
    ring.CloseRings()
    swath_polygon.AddGeometry(ring)

    return swath_polygon
                

#######################################
    

def azr_from_lonlat(IN_lon, IN_lat, IN_attributes, heau=0.):
    """
    Convert coordinates from lon-lat to azimuth-range for a given track
    
    :param IN_lon: longitude of points
    :type IN_lon: 1D-array of float
    :param IN_lat: latitude of points
    :type IN_lat: 1D-array of float
    :param IN_attributes
    :type IN_attributes
    :param heau:
    :type heau:
    
    :return OUT_az: azimuth coordinate of given points
    :rtype OUT_az: 1D-array of float
    :return OUT_r: range coordinate of given points
    :rtype OUT_r: 1D-array of float
    :return OUT_near_range
    :rtype OUT_near_range
    """
           
    # -----------
    # Preparation
    # -----------
    nr = IN_attributes.nr_cross_track
    dr = IN_attributes.range_sampling
    alt = IN_attributes.alt
    r0 = IN_attributes.near_range
    # ---------------------------------------------------------
    # Compute along-track (az) and across-track (y) coordinates
    # ---------------------------------------------------------
    # 1st iteration
    du = GEN_APPROX_RAD_EARTH * np.cos(IN_lat) * (IN_lon - math_fct.linear_extrap(IN_lat, IN_attributes.lat_init, IN_attributes.lon_init))
    lat0 = IN_lat + (du * np.sin(math_fct.linear_extrap(IN_lat, IN_attributes.lat_init, IN_attributes.heading_init)) * np.cos(math_fct.linear_extrap(IN_lat, IN_attributes.lat_init, IN_attributes.heading_init))) / GEN_APPROX_RAD_EARTH
    psi = math_fct.linear_extrap(lat0, IN_attributes.lat_init, IN_attributes.heading_init)
    y = du * np.cos(psi)  # eq (3)
    OUT_azcoord = (math_fct.linear_extrap(lat0, IN_attributes.lat_init, np.arange(len(IN_attributes.lat_init))))
    nb_points = 50
    gamma = np.zeros([len(IN_lat), nb_points])
    beta = np.zeros([len(IN_lat), nb_points])
  
    for i in range(nb_points):
        k = OUT_azcoord.astype('i4')+i-int(nb_points/2)
        bad_ind = np.logical_or((k < 0), (k > len(IN_attributes.costheta_init)-1))
        k[bad_ind] = 0.

        theta = np.pi/2. - IN_lat
        phi = IN_lon
        
        costheta_0 = IN_attributes.costheta_init[k, ]
        sintheta_0 = IN_attributes.sintheta_init[k, ]
        cosphi_0 = IN_attributes.cosphi_init[k, ]
        sinphi_0 = IN_attributes.sinphi_init[k, ]
        cospsi_0 = IN_attributes.cospsi_init[k, ]
        sinpsi_0 = IN_attributes.sinpsi_init[k, ]
        
        gamma[:, i] = (GEN_APPROX_RAD_EARTH+heau)*(np.sin(theta)*np.cos(phi)*(-cospsi_0*costheta_0*cosphi_0-sinpsi_0*sinphi_0) \
                                                   + np.sin(theta)*np.sin(phi)*(-cospsi_0*costheta_0*sinphi_0+sinpsi_0*cosphi_0) \
                                                   + np.cos(theta)*(cospsi_0*sintheta_0))
                
        beta[:, i] = (GEN_APPROX_RAD_EARTH+heau)*(np.sin(theta)*np.cos(phi)*(sinpsi_0*costheta_0*cosphi_0-cospsi_0*sinphi_0) \
                                                 + np.sin(theta)*np.sin(phi)*(sinpsi_0*costheta_0*sinphi_0+cospsi_0*cosphi_0) \
                                                 + np.cos(theta)*(-sinpsi_0*sintheta_0))

        gamma[bad_ind, i] = 9.99e20
        beta[bad_ind, i] = 9.99e20
                
    ind = np.zeros(len(IN_lat), int)   
    y = np.zeros(len(IN_lat), float)   
             
    for i in range(len(IN_lat)):
        indice = np.argmin(np.abs(gamma[i, :]))
        y[i] = beta[i, indice]
        ind[i] = indice - int(nb_points/2)
    OUT_azcoord2 = OUT_azcoord + ind #-0.5
    #~ # Compute range coordinate (across track)
    OUT_azcoord2[np.where(OUT_azcoord2 < 0)] = 0
    OUT_azcoord2[np.where(OUT_azcoord2 > len(alt)-1)] = len(alt)-1
    
    H = alt[OUT_azcoord2.astype('i4')]
    
    rr = np.sqrt((H-heau + (y ** 2) / (2 * GEN_APPROX_RAD_EARTH)) ** 2 + y ** 2)  # eq (5b)
    OUT_rcoord = (rr - r0) / dr  # eq (4)
    
    return OUT_azcoord2, OUT_rcoord
                

#######################################


def all_linear_rings(geom):
    """ Generator for all linear rings in a geometry """

    if geom.GetGeometryName() == 'MULTIPOLYGON':
        for polygon_index in range(geom.GetGeometryCount()):
            polygon = geom.GetGeometryRef(polygon_index)
            for line_index in range(polygon.GetGeometryCount()):
                yield polygon.GetGeometryRef(line_index)

    if geom.GetGeometryName() == 'POLYGON':
        for line_index in range(geom.GetGeometryCount()):
            yield geom.GetGeometryRef(line_index)


def intersect(input, output, indmax, IN_attributes, overwrite=False):

    # Function to modify polygon shapefile in radar geometry in order to detect overlapped polygons
    # Only intersection between 2 polygons are considered, not multi

    lyr, ds = my_shp.open_shp(input)
    lyr_bis, ds_bis = my_shp.open_shp(input)

    liste_lac = []

    out_ds, out_lyr = my_shp.overwrite_shapefile(output, ds.GetDriver().GetName(), lyr.GetGeomType(),
                                                 lyr.GetSpatialRef(), overwrite)

    defn = out_lyr.GetLayerDefn()
    multi = ogr.Geometry(ogr.wkbMultiPolygon)
    out_lyr.CreateField(ogr.FieldDefn(str('IND_LAC'), ogr.OFTInteger64))

    polygons = []

    for i, polygon_index_1 in enumerate(lyr):

        geom_1 = polygon_index_1.GetGeometryRef()
        geom_1 = geom_1.Buffer(0)

        flag_intersect = 0

        if not geom_1.IsValid():

            geom_1 = geom_1.Buffer(0)


        lyr_bis.SetSpatialFilter(geom_1)

        for polygon_index_2 in lyr_bis:

            # If same feature
            if polygon_index_2.GetFID() == polygon_index_1.GetFID():
                continue

            geom_2 = polygon_index_2.GetGeometryRef()
            geom_2 = geom_2.Buffer(0)

            if geom_1.Intersects(geom_2):
                intersection = geom_1.Intersection(geom_2)
                diff1 = geom_1.Difference(intersection)
                diff2 = geom_2.Difference(intersection)

                h1 = polygon_index_1.GetField(str("HEIGHT"))
                h2 = polygon_index_2.GetField(str("HEIGHT"))
                h = (h1 + h2) / 2

                out_feat = ogr.Feature(defn)

                out_feat.SetGeometry(intersection)
                out_feat.SetField(str("IND_LAC"), indmax + 1)
                out_lyr.CreateFeature(out_feat)

                out_feat1 = ogr.Feature(defn)
                out_feat1.SetGeometry(diff1)
                out_feat1.SetField(str("IND_LAC"), polygon_index_1.GetField(str("IND_LAC")))
                out_lyr.CreateFeature(out_feat1)
                out_feat2 = ogr.Feature(defn)
                out_feat2.SetGeometry(diff2)
                out_feat2.SetField(str("IND_LAC"), polygon_index_2.GetField(str("IND_LAC")))
                out_lyr.CreateFeature(out_feat2)

                lac = Reference_height_Lac(indmax + 1, intersection, defn, IN_attributes)
                lac.height = h
                lac.set_hmean(np.mean(lac.compute_h()))
                liste_lac.append(lac)
                indmax += 1

            else:

                out_feat = ogr.Feature(defn)
                out_feat.SetGeometry(polygon_index_1.GetGeometryRef())
                out_feat.SetField(str("IND_LAC"), polygon_index_1.GetField(str("IND_LAC")))
                out_lyr.CreateFeature(out_feat)

    out_ds.Destroy()
    ds.Destroy()
    ds_bis.Destroy()

    return liste_lac

def compute_near_range(IN_attributes, layer, cycle_number=0):
    hmean=0.
    if IN_attributes.height_model == 'reference_height':
        count, height_tot = 0, 0
        fields_count = layer.GetLayerDefn().GetFieldCount()
        for ind, polygon_index in enumerate(layer):
            for i in range(fields_count):
                if layer.GetLayerDefn().GetFieldDefn(i).GetName() == IN_attributes.height_name:
                    if polygon_index.GetField(str(IN_attributes.height_name)) is not None :
                        height_tot += np.float(polygon_index.GetField(str(IN_attributes.height_name)))
                        count +=1
        if count != 0:
            hmean = height_tot/count
            near_range = np.mean(np.sqrt((IN_attributes.alt - hmean + (IN_attributes.nr_cross_track ** 2) / (2 * GEN_APPROX_RAD_EARTH)) ** 2 + IN_attributes.nr_cross_track ** 2))
        else:
            near_range = 0
    if IN_attributes.height_model == 'gaussian':
        near_range = np.mean(np.sqrt((IN_attributes.alt + (IN_attributes.nr_cross_track ** 2) / (2 * GEN_APPROX_RAD_EARTH)) ** 2 + IN_attributes.nr_cross_track ** 2))
    if IN_attributes.height_model == 'polynomial':
        near_range = np.mean(np.sqrt((IN_attributes.alt + (IN_attributes.nr_cross_track ** 2) / (2 * GEN_APPROX_RAD_EARTH)) ** 2 + IN_attributes.nr_cross_track ** 2))
    if IN_attributes.height_model == 'reference_file':
        near_range = np.mean(np.sqrt((IN_attributes.alt + (IN_attributes.nr_cross_track ** 2) / (2 * GEN_APPROX_RAD_EARTH)) ** 2 + IN_attributes.nr_cross_track ** 2))
    if IN_attributes.height_model == None:
        hmean = IN_attributes.height_model_a + IN_attributes.height_model_a * \
        np.sin(2*np.pi * (np.mean(IN_attributes.orbit_time) + cycle_number * IN_attributes.cycle_duration) - IN_attributes.height_model_t0) / IN_attributes.height_model_period
        near_range = np.mean(np.sqrt((IN_attributes.alt - hmean + (IN_attributes.nr_cross_track ** 2) / (2 * GEN_APPROX_RAD_EARTH)) ** 2 + IN_attributes.nr_cross_track ** 2))
    my_api.printInfo("[write_polygons] [reproject_shapefile] Near_range =  %d" % (near_range))
    
    return near_range, hmean

def compute_tile_coords(IN_attributes, IN_swath, nadir_az_size, nb_pix_overlap_begin):

    inner_first_az = [nb_pix_overlap_begin]
    inner_first_rg = np.sqrt((IN_attributes.alt[inner_first_az] + (IN_attributes.nr_cross_track ** 2) / (2 * GEN_APPROX_RAD_EARTH)) ** 2 + IN_attributes.nr_cross_track ** 2)
    inner_first_lon, inner_first_lat = math_fct.lonlat_from_azy(inner_first_az, inner_first_rg, IN_attributes, IN_swath, h=0, IN_unit="deg")
    inner_first = (inner_first_lon[0], inner_first_lat[0])

    outer_first_az = [nb_pix_overlap_begin]
    outer_first_rg = np.sqrt((IN_attributes.alt[outer_first_az] + (IN_attributes.nr_cross_track ** 2) / (2 * GEN_APPROX_RAD_EARTH)) ** 2 + IN_attributes.nr_cross_track ** 2) + IN_attributes.nb_pix_range * IN_attributes.range_sampling
    outer_first_lon, outer_first_lat = math_fct.lonlat_from_azy(outer_first_az, outer_first_rg, IN_attributes, IN_swath, h=0, IN_unit="deg")
    outer_first = (outer_first_lon[0], outer_first_lat[0])

    inner_last_az = [nadir_az_size + nb_pix_overlap_begin-1]
    inner_last_rg = np.sqrt((IN_attributes.alt[inner_last_az] + (IN_attributes.nr_cross_track ** 2) / (2 * GEN_APPROX_RAD_EARTH)) ** 2 + IN_attributes.nr_cross_track ** 2)
    inner_last_lon, inner_last_lat = math_fct.lonlat_from_azy(inner_last_az, inner_last_rg, IN_attributes, IN_swath, h=0, IN_unit="deg")
    inner_last = (inner_last_lon[0], inner_last_lat[0])

    outer_last_az = [nadir_az_size + nb_pix_overlap_begin-1]
    outer_last_rg = np.sqrt((IN_attributes.alt[outer_last_az] + (IN_attributes.nr_cross_track ** 2) / (2 * GEN_APPROX_RAD_EARTH)) ** 2 + IN_attributes.nr_cross_track ** 2) + IN_attributes.nb_pix_range * IN_attributes.range_sampling
    outer_last_lon, outer_last_lat = math_fct.lonlat_from_azy(outer_last_az, outer_last_rg, IN_attributes, IN_swath, h=0, IN_unit="deg")
    outer_last = (outer_last_lon[0], outer_last_lat[0])

    inner_first = inner_first
    inner_last = inner_last
    outer_first = outer_first
    outer_last = outer_last
    tile_coords = (inner_first, inner_last, outer_first, outer_last)
    return tile_coords
