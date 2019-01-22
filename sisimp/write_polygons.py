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

import lib.my_api as my_api

from lib.my_variables import COEFF_X2, COEFF_Y2, COEFF_X, COEFF_Y, COEFF_XY, COEFF_CST, FACT_ECHELLE, RAD2DEG, DEG2RAD, GEN_APPROX_RAD_EARTH

import lib.height_model as height_model
import lib.true_height_model as true_height_model
from lib.roll_module import Roll_module
from lib.my_tools import llh2xyz
from lib.my_tools import xyz2llh
from cnes.modules.geoloc.lib.geoloc import pointcloud_height_geoloc_vect
import lib.dark_water_functions as dark_water
import proc_real_pixc
import proc_real_pixc_vec_river

import mathematical_function as math_fct


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
        self.height_name = None # Height field name
        self.height_model_min_area = None  # Optionnal argument to add complex 2D height model
        # Constant height model (always applied, set HEIGHT_MODEL_A to 0 to desactivate)
        self.height_model_a = None  # Height model A
        self.height_model_a_tab = None # Height model from shapefile
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
        self.swath_polygons = {}  # Dictionnary for storing swath polygons
        
        self.dw_detected_noise_height = None # dw detected noise tab

        # 3.1 - From orbit file
        self.azimuth_spacing = None  # Azimuth spacing
        self.lon = None
        self.lon_init = None
        self.lat = None
        self.lat_init = None
        self.alt = None
        self.heading = None
        self.heading_init = None
        self.orbit_time = None
        self.x = None
        self.y = None
        self.z = None


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

    :return OUT_burn_data: the radar pixels that are inside a water body
    :rtype OUT_burn_data: 2D-array of int (0=land 1=water)
    :rettun OUT_height_data : the height pixels that are inside a water body
    :rtype OUT_height_data : 2D-array of int (height of each water body pixel)
    """
    if IN_pixc_vec_only:
        my_api.printInfo("[write_polygons] == compute_pixels_in_water / river polygons only ==")
    else:
        my_api.printInfo("[write_polygons] == compute_pixels_in_water / all polygons ==")

    # 1 - Read the reprojected shapefile
    driver = ogr.GetDriverByName(str("ESRI Shapefile"))
    da_shapefile = driver.Open(IN_fshp_reproj, 0)  # 0 = read only
    layer = da_shapefile.GetLayer()
    if IN_pixc_vec_only:
        layer.SetAttributeFilter(str("RIV_FLAG != '0'"))
        my_api.printInfo("compute_pixels_in_water / river pixels only - %d features to deal with" % layer.GetFeatureCount())

    # Create a GDAL raster in memory
    nx = len(IN_attributes.lon)
    ny = IN_attributes.nb_pix_range
    
    cover = np.zeros((ny, nx), dtype='float64')
    cover_height = np.zeros((ny, nx), dtype='float64')
    cover_code = np.zeros((ny, nx), dtype='float64')
    # Create 2 raster band if heights value come from shapefile
    #~ nb_rasterBand = [1, 2][IN_attributes.height_model == "reference_height" or IN_attributes.height_model == "gaussian" or IN_attributes.height_model == "polynomial"]
    nb_rasterBand = 3

    ds = gdal.GetDriverByName(str('MEM')).Create('', nx, ny, nb_rasterBand, GDT_Float32)
    ds.SetGeoTransform([-0.5, 1, 0, -0.5, 0, 1])
    ds.GetRasterBand(1).WriteArray(cover)

    # Burn with a value 1 the pixels that are in a polygon
    gdal.RasterizeLayer(ds, [1], layer, None, None, burn_values=[1])
    #                                        , options=['ALL_TOUCHED=TRUE'])

    # Convert the radar pixels touching water in lon-lat
    OUT_burn_data = ds.GetRasterBand(1).ReadAsArray()
    
    # Analyse value in shapefile
    OUT_height_data = None
    OUT_code_data = None
    
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
 
    # Close the raster
    ds = None
 
    return OUT_burn_data, OUT_height_data, OUT_code_data

def write_water_pixels_realPixC(IN_water_pixels, IN_swath, IN_cycle_number, IN_orbit_number, IN_attributes):
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
    """
    my_api.printInfo("[write_polygons] == write_water_pixels_realPixC ==")  
    
    ################################
    # Variables for PixC main file #
    ################################
    
    # Print number of water pixels
    size_of_tabs = np.count_nonzero(IN_water_pixels) 
    my_api.printInfo(str("Nb water pixels: %d" % size_of_tabs))

    # 1 - Get range and azimuth indices of all water pixels
    ind = np.nonzero(IN_water_pixels)  # Get indices 1=lake and 2=river (remove 0=land)
    r, az = [ind[0], ind[1]]  # Range index
    #az = ind[1]  # Azimuth index
    river_flag = IN_water_pixels[ind]  # 1=lake and 2=river
    
    # Dark Water
    
    ## check if dark water is simulated or not
    if IN_attributes.dark_water.lower() == "yes" :
        ### Simulate dark water if water pixel are present
        ### to do integrer dans le if suivant après quand ça marche pour éviter les répétitions
        if size_of_tabs != 0.:
            rmin, rmax, azmin, azmax = r.min(), r.max(), az.min(), az.max()
            taille_r, taille_az = rmax-rmin+1, azmax-azmin+1

            # Simulate dark_water
            dw_mask=dark_water.dark_water_simulation([taille_az,taille_r],IN_attributes.fact_echelle_dw, IN_attributes.dw_pourcent,IN_attributes.dw_seed)

            ## Get water extent
            indice_r = np.array(ind[0]-rmin)
            indice_az = np.array(ind[1]-azmin)

            ## Randomly classify or erase dark water regions
            ### Randomly erase DW regions in DW mask
            dw_mask = dark_water.dark_water_non_detected_simulation(dw_mask,IN_attributes.scale_factor_non_detected_dw,IN_attributes.dw_detected_percent,IN_attributes.dw_seed)
            #reshape dark_water to water extent
            dw_mask=dw_mask[indice_az,indice_r]
            ### Update IN_water_pixels with the deleted water pixels
            # check if dw pixels are deleted = dw_mask pixels set to 2
            if np.where(dw_mask==2)[0].size > 0:
                my_api.printInfo(str("Nb detected dark water pixels : %d" % np.where(dw_mask==1)[0].size ))
                my_api.printInfo(str("Nb non detected dark water pixels : %d" % np.where(dw_mask==2)[0].size ))
                ### Update river_flag and IN_water_pixels with the deleted water pixels
                river_flag[np.where(dw_mask==2)]=0
                IN_water_pixels[ind]=river_flag
                # delete the corresponding pixels in the dw mask to update indices values
                dw_mask=np.delete(dw_mask,np.where(dw_mask==2))
                ## Update size of tabs etc... because water pixels were deleted
                size_of_tabs = np.count_nonzero(IN_water_pixels)
                my_api.printInfo(str("Nb water pixels: %d" % size_of_tabs))
                ind = np.nonzero(IN_water_pixels)
                r, az = [ind[0], ind[1]]
                river_flag = IN_water_pixels[ind]
            ###  Build the classification array with water flag Dark_water flag
            classification_tab = np.ones(size_of_tabs) * IN_attributes.water_flag  # Classification as water
            #locate DW pixels
            dark_water_loc = np.where(dw_mask==1)
            #update classification value for dark water pixels with DW flag
            classification_tab[dark_water_loc]= IN_attributes.darkwater_flag
        else :
            classification_tab = np.ones(size_of_tabs) * IN_attributes.water_flag
    else :
        classification_tab = np.ones(size_of_tabs) * IN_attributes.water_flag
        my_api.printInfo(str("No Dark Water will be simulated"))

    if IN_attributes.height_model_a_tab is not None:
        height_flag = IN_attributes.height_model_a_tab[ind]
        
    if IN_attributes.code is not None:
        code_flag = IN_attributes.code[ind]
        
    # 2 - Compute radar variables for water pixels
    r0 = np.sqrt((IN_attributes.alt + (IN_attributes.nr_cross_track ** 2) / (2 * GEN_APPROX_RAD_EARTH)) ** 2 + IN_attributes.nr_cross_track ** 2)  # Radar to ground distance for near range pixels
    ri = r0[az] + r * IN_attributes.range_sampling  # Radar-to-ground distance
    Hi = IN_attributes.alt[az]  # Altitude
    angles = np.arccos(Hi/ri)  # Look angles
    pixel_area = IN_attributes.azimuth_spacing * IN_attributes.range_sampling / np.sin(angles)  # Pixel area
    
    # 3 - Build cross-track distance array
    # Compute theorical cross-track distance for water pixels
    sign = [-1, 1][IN_swath.lower() == 'right']
    y = sign * np.sqrt((ri + Hi) * (ri - Hi) / (1. + Hi / GEN_APPROX_RAD_EARTH))
    lon, lat = math_fct.lonlat_from_azy(az, y, IN_attributes.lat_init, IN_attributes.lon_init, IN_attributes.heading_init, IN_unit="deg")
    
    # 4 - Height model    
    ## TBD : Separate height model for each water body !!!
    if size_of_tabs != 0.:
        
        # 4.1 - Constant elevation model
        if IN_attributes.height_model is None:
            # Compute theorical constant elevation for water pixels
            elevation_tab = make_elevation_tab(az, IN_cycle_number, IN_attributes)
      
        # 4.2 - Gaussian model
        elif IN_attributes.height_model == 'gaussian':
            
            # Compute theorical constant elevation for water pixels
            elevation_tab = make_elevation_tab(az, IN_cycle_number, IN_attributes)
            
            # Add gaussian model over big lakes
            for i in np.unique((IN_attributes.code[ind])):
                indice=np.where(IN_attributes.code[ind]==i)
                size_water_body = pixel_area[indice].sum()
                if size_water_body > IN_attributes.height_model_min_area*10000.:
                    ### First height model : random field convoluted with gaussian
                    my_api.printInfo(str("Gaussian model applied for big water body of size %d ha" % int(size_water_body/10000)))
                    rmin, rmax, azmin, azmax = r[indice].min(), r[indice].max(), az[indice].min(), az[indice].max()
                    taille_r, taille_az = rmax-rmin+1, azmax-azmin+1
                    indice_az = np.array(az[indice]-azmin)
                    indice_r = np.array(r[indice]-rmin)
                    height = height_model.generate_2d_profile_gaussian([taille_az, taille_r], 0., "Default", IN_attributes.height_model_stdv, FACT_ECHELLE)
                    height_water = height[indice_az, indice_r]
                    elevation_tab[indice] += height_water

        # 4.3 - Polynomial model
        elif IN_attributes.height_model == 'polynomial':
            
            # Compute theorical constant elevation for water pixels
            elevation_tab = make_elevation_tab(az, IN_cycle_number, IN_attributes)
            
            # Add polynomial model over big lakes
            for i in np.unique((IN_attributes.code[ind])):
                indice=np.where(IN_attributes.code[ind]==i)
                size_water_body = pixel_area[indice].sum()
                if size_water_body > IN_attributes.height_model_min_area*10000.:
        
                    ### Second height model : 2D polynomial model (center in x0, y0)
                    my_api.printInfo(str("Polynomial model applied for big water body of size %d ha" % int(size_water_body/10000)))
                    x_c, y_c, zone_number, zone_letter = utm.from_latlon(lat[0], lon[0])
                    # Convert pixel cloud to UTM (zone of the centroid)
                    latlon = pyproj.Proj(init="epsg:4326")
                    utm_proj = pyproj.Proj("+proj=utm +zone={}{} +ellps=WGS84 +datum=WGS84 +units=m +no_defs".format(zone_number, zone_letter))
                    X, Y = pyproj.transform(latlon, utm_proj, lon[indice], lat[indice])                    
                    k0 = np.random.randint(len(indice[0]))
                    X0, Y0 = X[k0], Y[k0]
                    height_water = height_model.generate_2d_profile_2nd_order_list(X0, Y0, X, Y, COEFF_X2, COEFF_Y2, COEFF_X, COEFF_Y, COEFF_XY, COEFF_CST)
                    elevation_tab[indice] += height_water

        # 4.4 - Height given by an attribute in input shapefile
        elif IN_attributes.height_model == "reference_height" and IN_attributes.height_model_a_tab is not None:
            elevation_tab = height_flag
    
        # 4.5 - Height given in a dedicated file
        elif IN_attributes.height_model == "reference_file" and IN_attributes.trueheight_file is not None:
            # Process true height model from Kevin Larnier
            # TDB : Add security for lat lon boundaries
            # TBD : Add specific model for 1D model (river)
            true_height_model_inst = true_height_model.TrueHeightModel(IN_attributes.trueheight_file, lat, lon, verbose=True)
            true_height_model_inst.apply_model()
            elevation_tab = true_height_model_inst.final_height
    


    # 5 - Error model


    # 4.1 - Compute noise over height
    if IN_attributes.dark_water.lower() == "yes" :
        delta_h=np.zeros(elevation_tab.shape)
        water_pixels=np.where(classification_tab==IN_attributes.water_flag)
        dw_pixels=np.where(classification_tab==IN_attributes.darkwater_flag)
        delta_h[water_pixels]=math_fct.calc_delta_h(angles[water_pixels],IN_attributes.noise_height,IN_attributes.height_bias_std)
        delta_h[dw_pixels]=math_fct.calc_delta_h(angles[dw_pixels],IN_attributes.dw_detected_noise_height,IN_attributes.height_bias_std)
        #delta_h = math_fct.calc_delta_h(angles, IN_attributes.noise_height, IN_attributes.height_bias_std)
        #~ print('delta_h[water_pixels]',delta_h[water_pixels])
        #~ print('delta_h[dw_pixels]',delta_h[dw_pixels])
    else : 
        delta_h = math_fct.calc_delta_h(angles, IN_attributes.noise_height, IN_attributes.height_bias_std)

    # 4.2 Add residual roll error
    try:
        
        roll = Roll_module(IN_attributes.roll_file)
        roll.interpolate_roll_on_sensor_grid(IN_attributes.orbit_time)
        
        ## Change roll values to simulate random acquisitions
        ## TBD ##
        
        # Apply roll for each pixel
        pixel_cloud_time = IN_attributes.orbit_time[az]
        roll.interpolate_roll_on_pixelcloud(IN_attributes.orbit_time, pixel_cloud_time)

        if IN_swath.lower() == 'right' :  # Sign depends on left / right swath
            delta_h_roll = (roll.roll2_err_cloud-roll.roll2_cor_cloud)*y*1e-6
            delta_h += delta_h_roll
            
        if IN_swath.lower() == 'left' :  # Sign depends on left / right swath
            delta_h_roll = (roll.roll1_err_cloud-roll.roll1_cor_cloud)*y*1e-6
            delta_h += delta_h_roll
    
        ## Check what is the better value from roll_module to use as error

    except:
        my_api.printInfo("No roll error applied")
       
    # 5.3 - Compute final noisy heights (elevation + thermal noise + roll error + height model) 
    elevation_tab_noisy = elevation_tab + delta_h           
       
    # 5.4 - Compute noise over geolocation
    
    #~ delta_y = math_fct.calc_delta_sensor(delta_h, Hi, y, ri)
    delta_y = math_fct.calc_delta_sensor(delta_h, Hi, y)
    
    # 5.5 - Add noise over geolocation if asked
    if IN_attributes.geolocalisation_improvement in ['no', 'non', 'not', 'false', 'nope']:
        my_api.printInfo("Add noise to cross track")
        y_noisy = y + delta_y
    else:
        my_api.printInfo("Geolocalisation improvement")
        y_noisy = y   
         
    # 6 - Build geolocation arrays
    # 6.1 - With no noise
    lon, lat = math_fct.lonlat_from_azy(az, y, IN_attributes.lat_init, IN_attributes.lon_init, IN_attributes.heading_init) 
    lon *= RAD2DEG
    lat *= RAD2DEG
    
    # 7 - Build velocity arrays
    nb_pix_nadir = IN_attributes.x.size  # Nb pixels at nadir
    # Init velocity arrays
    vx = np.zeros(nb_pix_nadir)
    vy = np.zeros(nb_pix_nadir)
    vz = np.zeros(nb_pix_nadir)
    # Compute first value
    vx[0], vy[0], vz[0] = [(IN_attributes.x[1] - IN_attributes.x[0]), (IN_attributes.y[1] - IN_attributes.y[0]), (IN_attributes.z[1] - IN_attributes.z[0])]  / (IN_attributes.orbit_time[1] - IN_attributes.orbit_time[0])
    # Compute last value
    vx[-1], vy[-1], vz[-1] = [(IN_attributes.x[-1] - IN_attributes.x[-2]), (IN_attributes.y[-1] - IN_attributes.y[-2]), (IN_attributes.z[-1] - IN_attributes.z[-2])] / (IN_attributes.orbit_time[-1] - IN_attributes.orbit_time[-2])
    # Compute middle values
    for indp in range(1, nb_pix_nadir-1):
        vx[indp], vy[indp], vz[indp] = [(IN_attributes.x[indp+1] - IN_attributes.x[indp-1]), (IN_attributes.y[indp+1] - IN_attributes.y[indp-1]), (IN_attributes.z[indp+1] - IN_attributes.z[indp-1])] / (IN_attributes.orbit_time[indp+1] - IN_attributes.orbit_time[indp-1])
 
    # 8 - Compute noisy geolocation and height
    complex_latlon_noise = False
   
    # Convert nadir latitudes/longitudes in degrees
    nadir_lat_deg = IN_attributes.lat[1:-1] * RAD2DEG
    nadir_lon_deg = IN_attributes.lon[1:-1] * RAD2DEG
    nadir_x = IN_attributes.x[1:-1]
    nadir_y = IN_attributes.y[1:-1]
    nadir_z = IN_attributes.z[1:-1]
        
    # Remove 1st and last values because just here for extrapolators (cf. read_orbit)
    nadir_alt = IN_attributes.alt[1:-1]
    nadir_heading = IN_attributes.heading[1:-1]
     
    
    if (size_of_tabs != 0.) & (complex_latlon_noise == True):
        my_api.printInfo("Complex latlon noise applied from improved geoloc module")
        
        # Update positions
        
        p = project_array((np.vstack([lat, lon, elevation_tab])).T, srcp='latlon', dstp='geocent')
        s1 = project_array((np.vstack([nadir_lon_deg[az], nadir_lat_deg[az], nadir_alt[az]])).T, srcp='latlon', dstp='geocent')
        s = np.transpose([nadir_x[az], nadir_y[az], nadir_z[az]])


        ri_new = np.sqrt((p[:,0]-s[:,0])**2+(p[:,1]-s[:,1])**2+(p[:,2]-s[:,2])**2)
                                                   
        p_final, p_final_llh, h_mu, (iter_grad,nfev_minimize_scalar) = pointcloud_height_geoloc_vect(p, elevation_tab,
                                                   s, np.transpose(np.array([vx[az], vy[az], vz[az]])), 
                                            ri_new, elevation_tab_noisy, 
                                            recompute_Doppler=True, 
                                            #if False uses 0 Doppler, else the value computed from p, s, vs
                                            recompute_R=True,
                                            verbose = False,
                                            max_iter_grad=1, height_goal = 1.e-3, safe_flag=True)
     

        lat_noisy, lon_noisy, elevation_tab_noisy = p_final_llh[:,0], p_final_llh[:,1], p_final_llh[:,2]

    else:
        lon_noisy, lat_noisy = math_fct.lonlat_from_azy(az, y_noisy, IN_attributes.lat_init, IN_attributes.lon_init, IN_attributes.heading_init)
        lon_noisy *= RAD2DEG  # Conversion in degrees
        lat_noisy *= RAD2DEG  # Conversion in degrees 
        
         
    ######################
    # Write output files #
    ######################


    
    # Cut arrays in order to write real PixC files
    # Tiles correspond to 1deg of latitude at nadir
    if lat.size != 0:        
    
        # Get min/max nadir latitudes values
        nadir_lat_max = int(np.max(nadir_lat_deg))
        nadir_lat_min = int(np.min(nadir_lat_deg))

        tile_db = IN_attributes.tile_database
        
        print(tile_db[:,2:3])
        exit()

        for cur_nadir_lat in range(nadir_lat_min, nadir_lat_max + 1):  # Loop on nadir latitude integer intervals

            my_api.printInfo("Dealing with latitudes >= %d AND < %d" % (cur_nadir_lat, cur_nadir_lat+1))
            
            # Get azimuth indices corresponding to this integer value of latitude
            nadir_az = np.where(nadir_lat_deg.astype(int) == cur_nadir_lat)[0]
            az_min = np.sort(nadir_az)[0]  # Min azimuth index, to remove from tile azimuth indices vector
            my_api.printInfo("= %d pixels in azimuth (index %d put to 0)" % (nadir_az.size, az_min))
            
            # Get pixel indices of water pixels corresponding to this latitude interval
            az_indices = np.where((az >= min(nadir_az)) & (az <= max(nadir_az)))[0]
            nb_pix = az_indices.size  # Number of water pixels for this latitude interval
            my_api.printInfo("= %d water pixels" % nb_pix)
            
            if az_indices.size != 0:  # Write water pixels at this latitude
                
                sub_az, sub_r = [az[az_indices], r[az_indices]]
                
                my_api.printInfo("Min r ind = %d - Max r ind = %d" % (np.min(sub_r), np.max(sub_r)))
                my_api.printInfo("Min az ind = %d - Max az ind = %d" % (np.min(sub_az), np.max(sub_az)))
                
                # Get output filename
                # North / south lat flag
                nord_or_south = ["S", "N"][cur_nadir_lat > 0]

                # Left / right swath flag
                left_or_right = IN_swath.upper()[0]
                
                # General tile reference
                tile_ref = str(abs(cur_nadir_lat)).zfill(2) + nord_or_south + "-" + left_or_right
                
                # Init L2_HR_PIXC object
                my_pixc = proc_real_pixc.l2_hr_pixc(sub_az-az_min, sub_r, classification_tab[az_indices], pixel_area[az_indices],
                                                    lat_noisy[az_indices], lon_noisy[az_indices], elevation_tab_noisy[az_indices], y[az_indices],
                                                    IN_attributes.orbit_time[nadir_az], nadir_lat_deg[nadir_az], nadir_lon_deg[nadir_az], nadir_alt[nadir_az], nadir_heading[nadir_az],
                                                    IN_attributes.x[nadir_az], IN_attributes.y[nadir_az], IN_attributes.z[nadir_az], vx[nadir_az], vy[nadir_az], vz[nadir_az], r0[nadir_az],
                                                    IN_attributes.mission_start_time, IN_attributes.cycle_duration, IN_cycle_number, IN_orbit_number, tile_ref, IN_attributes.nb_pix_range, nadir_az.size, IN_attributes.azimuth_spacing, IN_attributes.range_sampling, IN_attributes.near_range)
                
                # Update filenames with tile ref
                IN_attributes.sisimp_filenames.updateWithTileRef(tile_ref, IN_attributes.orbit_time[nadir_az[0]], IN_attributes.orbit_time[nadir_az[-1]])
                
                # Write main file
                my_pixc.write_pixc_file(IN_attributes.sisimp_filenames.pixc_file+".nc", None, True)
                
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
                    my_pixc_vec.set_vectorproc(lat[az_indices], lon[az_indices], elevation_tab[az_indices])
                    # Compute river_flag
                    my_pixc_vec.set_river_lake_tag(river_flag[az_indices]-1)  # -1 to have 0=lake and 1=river
                    # Write PIXCVec file
                    my_pixc_vec.write_file(IN_attributes.sisimp_filenames.pixc_vec_river_file+".nc", None, True)
                    # Write as shapefile if asked
                    if IN_attributes.create_shapefile:
                        my_pixc_vec.write_file_asShp(IN_attributes.sisimp_filenames.pixc_vec_river_file+".shp")
            
    else:  
        my_api.printInfo("No output data file to write")   


def reproject_shapefile(IN_filename, IN_swath, IN_driver, IN_attributes):
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
    
    :return OUT_file/Getname: full path of the shapefile with the water body polygons in radar coordinates
    :rtype OUT_filename: string
    :return OUT_swath_polygons
    :rtype OUT_swath_polygons
    """
    my_api.printInfo("[write_polygons] == reproject_shapefile ==")

    # 1 - Make swath polygons
    swath_polygon = make_swath_polygon(IN_swath, IN_attributes)
    OUT_swath_polygons = IN_attributes.swath_polygons
    OUT_swath_polygons[IN_swath] = swath_polygon

    # 2 - Select water bodies polygons in the swath
    da_shapefile = IN_driver.Open(IN_filename, 0)  # 0 = read only
    layer = da_shapefile.GetLayer()
    layer.SetSpatialFilter(swath_polygon)
    nb_features = layer.GetFeatureCount()
    
    my_api.printInfo("There are %d feature(s) crossing %s swath" % (nb_features, IN_swath))
    sys.stdout.flush()
    
    # Exit if no feature to deal with
    if nb_features == 0:
        return None, IN_attributes
    
    # 3 - Create the output shapefile
    swath_t = '%s_swath' % IN_swath
    OUT_filename = os.path.join(IN_attributes.out_dir, os.path.splitext(os.path.split(IN_filename)[1])[0] + '_tmp_radarproj_%s.shp' % swath_t)
    if os.path.exists(OUT_filename):
        IN_driver.DeleteDataSource(OUT_filename)
    dataout = IN_driver.CreateDataSource(OUT_filename)
    if dataout is None:
        my_api.printError("Could not create file")
    layerout = dataout.CreateLayer(str(os.path.splitext(os.path.split(OUT_filename)[1])[0]), None, geom_type=ogr.wkbPolygon)
    # Create necessary output fields 
    layerout.CreateField(ogr.FieldDefn(str('RIV_FLAG'), ogr.OFTInteger))
    layerout.CreateField(ogr.FieldDefn(str('HEIGHT'), ogr.OFTReal))
    layerout.CreateField(ogr.FieldDefn(str('CODE'), ogr.OFTInteger64))

    floutDefn = layerout.GetLayerDefn()
    feature_out = ogr.Feature(floutDefn)

    # 4 - Convert coordinates for each water body polygon
    range_tab = []
    OUT_near_range = None
    for polygon_index in layer:
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
            # 4.2.3 - Convert polygons coordinates
            add_ring = False
            for ring in all_linear_rings(intersection):
                npoints = ring.GetPointCount()

                if npoints == 0:
                    continue  # ignore polygons completely outside the swath

                points = np.transpose(np.array(ring.GetPoints()))
                lon = points[0] * DEG2RAD
                lat = points[1] * DEG2RAD

                az, r, IN_attributes.near_range = azr_from_lonlat(lon, lat, IN_attributes)
                range_tab = np.concatenate((range_tab, r), -1)
                npoints = len(az)
                if len(az) != len(lon):
                    my_api.printDebug("Ignore polygons crossing the swath")
                    exit()
                    continue  # Ignore polygons crossing the swath
                for p in range(npoints):  # no fonction ring.SetPoints()
                    ring.SetPoint(p, az[p], r[p])
                ring.CloseRings()
                add_ring = True
                # Add the reprojected ring to the output geometry
                geom_out.AddGeometry(ring)
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
                    # Test 'HEIGHT' parameter in input shapefile fields
                    if layerDefn.GetFieldDefn(i).GetName() == IN_attributes.height_name:
                        height_from_shp = True
                        feature_out.SetField(str("HEIGHT"), polygon_index.GetField(str(IN_attributes.height_name)))
                        
                    if IN_attributes.height_model == 'polynomial' or IN_attributes.height_model == 'gaussian':
                        feature_out.SetField(str("CODE"),polygon_index.GetField("code"))
                if not height_from_shp:
                    IN_attributes.height_model_a_tab = None
                # Add the output feature to the output layer
                layerout.CreateFeature(feature_out)
    IN_attributes.swath_polygons = OUT_swath_polygons
    
    return OUT_filename, IN_attributes

def project_array(coordinates, srcp='latlon', dstp='geocent'):
    """
    Project a numpy (n,2) array in projection srcp to projection dstp
    Returns a numpy (n,2) array.
    """
    p1 = pyproj.Proj(proj=srcp, datum='WGS84')
    p2 = pyproj.Proj(proj=dstp, datum='WGS84')
    fx, fy, fz = pyproj.transform(p1, p2, coordinates[:,1], coordinates[:,0], coordinates[:,2])
    # Re-create (n,2) coordinates
    # Inversion of lat and lon !
    return np.dstack([fx, fy, fz])[0]
        
def make_elevation_tab(IN_az, IN_cycle_number, IN_attributes):
    """
    Make the elevation array without noise 
    
    :param IN_az: the azimuth indices
    :type IN_az: 1D-array of int 
    :param IN_cycle_number: cycle number of the currently processed orbit
    :type IN_cycle_number: int
    :param IN_attributes
    :type IN_attributes

    :return OUT_theorical_height: associated elevations
    :rtype OUT_theorical_height: 1D-array of float
    """
    return IN_attributes.height_model_a * np.sin(2*np.pi * ((IN_attributes.orbit_time[IN_az] + IN_cycle_number * IN_attributes.cycle_duration) - IN_attributes.height_model_t0) / IN_attributes.height_model_period) 


def make_swath_polygon(IN_swath, IN_attributes):
    """Make left of right swath polygon
    
    :param IN_swath 
    :type IN_swath
    """
    sign = [-1, 1][IN_swath.lower() == 'right']
    ymin = sign * IN_attributes.nr_cross_track
    ymax = sign * IN_attributes.swath_width/2

    n = len(IN_attributes.lon_init) - 4
    az = np.arange(2, n + 2, 10)
    y = ymin * np.ones(len(az))
    lon1, lat1 = math_fct.lonlat_from_azy(az, y, IN_attributes.lat_init, IN_attributes.lon_init, IN_attributes.heading_init)
    y = ymax * np.ones(len(az))
    lon2, lat2 = math_fct.lonlat_from_azy(az, y, IN_attributes.lat_init, IN_attributes.lon_init, IN_attributes.heading_init)
    lonswath = np.concatenate((lon1, lon2[::-1]))
    latswath = np.concatenate((lat1, lat2[::-1]))
    lonswath *= RAD2DEG
    latswath *= RAD2DEG

    swath_polygon = ogr.Geometry(ogr.wkbPolygon)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    n = len(lonswath)
    for i in range(n):
        ring.AddPoint(lonswath[i], latswath[i])
    ring.CloseRings()
    swath_polygon.AddGeometry(ring)

    return swath_polygon


def azr_from_lonlat(IN_lon, IN_lat, IN_attributes):
    """
    Convert coordinates from lon-lat to azimuth-range for a given track
    
    :param IN_lon: longitude of points
    :type IN_lon: 1D-array of float
    :param IN_lat: latitude of points
    :type IN_lat: 1D-array of float
    :param IN_attributes
    :type IN_attributes
    
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
    daz = IN_attributes.azimuth_spacing
    alt = IN_attributes.alt

    # ---------------------------------------------------------
    # Compute along-track (az) and across-track (y) coordinates
    # ---------------------------------------------------------
    # 1st iteration
    du = GEN_APPROX_RAD_EARTH * np.cos(IN_lat) * (IN_lon - math_fct.linear_extrap(IN_lat, IN_attributes.lat_init, IN_attributes.lon_init))
    lat0 = IN_lat + (du * np.sin(math_fct.linear_extrap(IN_lat, IN_attributes.lat_init, IN_attributes.heading_init)) * np.cos(math_fct.linear_extrap(IN_lat, IN_attributes.lat_init, IN_attributes.heading_init))) / GEN_APPROX_RAD_EARTH
    psi = math_fct.linear_extrap(lat0, IN_attributes.lat_init, IN_attributes.heading_init)
    y = du * np.cos(psi)  # eq (3)
    OUT_azcoord = math_fct.linear_extrap(lat0, IN_attributes.lat_init, np.arange(len(IN_attributes.lat_init)))
    lon_prec, lat_prec = math_fct.lonlat_from_azy(OUT_azcoord, y, IN_attributes.lat_init, IN_attributes.lon_init, IN_attributes.heading_init)
    
    # Next iterations
    precision = 10
    while precision >= 0.5:
        du = GEN_APPROX_RAD_EARTH * (lon_prec - IN_lon) * np.cos(IN_lat)
        dv = GEN_APPROX_RAD_EARTH * (lat_prec - IN_lat)
        az = OUT_azcoord - (du * np.sin(psi) + dv * np.cos(psi)) / daz  
        psi = math_fct.linear_extrap(az, np.arange(len(IN_attributes.lat_init)), IN_attributes.heading_init)
        last_az = OUT_azcoord
        OUT_azcoord = OUT_azcoord - (du * np.sin(psi) + dv * np.cos(psi)) / daz  
        y = y - du * np.cos(psi) + dv * np.sin(psi)
        lon_prec, lat_prec = math_fct.lonlat_from_azy(OUT_azcoord, y, IN_attributes.lat_init, IN_attributes.lon_init, IN_attributes.heading_init)
        precision = np.max(np.abs(OUT_azcoord - last_az))  # iterate until convergence
        
    # Compute range coordinate (across track)
    H = alt[OUT_azcoord.astype('i4')]
    r0 = np.sqrt((H + (nr ** 2) / (2 * GEN_APPROX_RAD_EARTH)) ** 2 + nr ** 2)
    OUT_rcoord = np.sqrt((H + (y ** 2) / (2 * GEN_APPROX_RAD_EARTH)) ** 2 + y ** 2)  # eq (5b)
    OUT_rcoord = (OUT_rcoord - r0) / dr  # eq (4)
    
    OUT_near_range = r0
    
    return OUT_azcoord, OUT_rcoord, OUT_near_range


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
