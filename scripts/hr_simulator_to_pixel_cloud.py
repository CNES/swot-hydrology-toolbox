'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


from netCDF4 import Dataset
from os.path import join, abspath
import argparse
import numpy as np
import os

def write_annotation_file(ann_file, 
                          pixc_file):
    """
    write the river-annotation.rdf file so that lake processor can run
    """
    f = open(ann_file,'w')
    f.write("l2pixc file = %s\n"%pixc_file)
    f.close()
    
    

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('pixel_cloud_file', type=str)
    parser.add_argument('sensor_file', type=str)
    parser.add_argument('output_final_pixel_cloud', default= None, type=str)

    args = vars(parser.parse_args())


    pixc = Dataset(abspath(args['pixel_cloud_file']), 'r')
    sensor = Dataset(abspath(args['sensor_file']), 'r')
    final_pixc = Dataset(abspath(args['output_final_pixel_cloud']), 'w')


    nb_pixels = pixc.dimensions["water_record"].size
    nb_tvp_record = sensor.dimensions["record"].size

    var_pixc_names = [['azimuth_index', int, 'azimuth_index'], ['range_index', int, 'range_index'], ['interferogram', np.float32, None], ['power_plus_y', np.float32, None],\
    ['power_minus_y', np.float32, None], ['coherent_power', np.float64, None],['x_factor_plus_y', np.float64, None],['x_factor_minus_y', np.float64, None],['water_frac', np.float64, None],\
    ['water_frac_uncert', np.float64, None], ['classification', np.int, 'classification'], ['false_detection_rate', np.float64, None], ['missed_detection_rate', np.float64, None],\
    ['prior_water_prob', np.float64, None], ['bright_land_flag', np.float64, None],['layover_impact', np.float32, None], ['num_rare_looks', np.int, None],['latitude', np.float64, 'latitude_medium'],\
    ['longitude', np.float64, 'longitude_medium'], ['height', np.float64, 'height_medium'], ['cross_track', np.float64, 'cross_track_medium'], ['pixel_area', np.float64, 'pixel_area'],\
    ['inc', np.float64, None], ['phase_noise_std', np.float32, 'dphase_medium'], ['dlatitude_dphase', np.float64, None], ['dlongitude_dphase', np.float64, None],\
    ['dheight_dphase', np.float64, 'dheight_dphase_medium'], ['dheight_droll', np.float64, None], ['dheight_dbaseline', np.float64, None], ['dheight_drange', np.float64, None],\
    ['darea_dheight', np.float64, None], ['illumination_time', np.float64, None], ['illumination_time_tai', np.float64, None],['num_med_looks', np.float32, None], ['sig0', np.float32, None],\
    ['phase_unwrapping_region', np.float64, None], ['instrument_range_cor', np.float64, None], ['instrument_phase_cor', np.float64, None], ['instrument_baseline_cor', np.float64, None],\
    ['instrument_attitude_cor', np.float64, None], ['model_dry_tropo_cor', np.float64, None], ['model_wet_tropo_cor', np.float64, None], ['iono_cor_gim_ka', np.float64, None],\
    ['xover_height_cor', np.float64, None], ['load_tide_sol1', np.float64, None], ['load_tide_sol2', np.float64, None], ['pole_tide', np.float64, None], ['solid_earth_tide', np.float64, None],\
    ['geoid', np.float64, None], ['surface_type_flag', np.float64, None], ['pixc_qual', np.float64, None]]
    
      
    var_tvp_names  = [['time', np.float64, 'time'],['time_tai', np.float64, None], ['latitude', np.float64, 'latitude'],['longitude', np.float64, 'longitude'],['height', np.float64, 'altitude'],\
    ['roll', np.float64, None],['pitch', np.float64, None],['yaw', np.float64, None],['velocity_heading', np.float64, 'heading'],['x', np.float64, 'x'],\
    ['y', np.float64, 'y'],['z', np.float64, 'z'],['vx', np.float64, 'velocity_unit_x'],['vy', np.float64, 'velocity_unit_y'],['vz', np.float64, 'velocity_unit_z'],\
    ['plus_y_antenna_x', np.float64, 'baseline_right_x'],['plus_y_antenna_y', np.float64, 'baseline_right_y'],['plus_y_antenna_z', np.float64, 'baseline_right_z'],\
    ['minus_y_antenna_x', np.float64, 'baseline_left_x'],['minus_y_antenna_y', np.float64, 'baseline_left_y'],['minus_y_antenna_z', np.float64, 'baseline_left_z'],\
    ['record_counter', np.float64, None], ['sc_event_flag', np.float64, None],['tvp_qual', np.float64, None]]
  

  
    pixc_group = final_pixc.createGroup("pixel_cloud") 
    pixc_group.createDimension("record", nb_pixels)
    pixc_group.createDimension("depth", 2)

    latitude = pixc.variables['latitude_medium'][:]
    longitude = pixc.variables['longitude_medium'][:]
    
    ## TO BE CHANGED !!!!
    # Get output filename
    # North / south lat flag
    if(np.mean(latitude) > 0 ):
        nord_or_south= "N"
    else:
        nord_or_south= "S"
    # Left / right swath flag
    print(args['pixel_cloud_file'][-17:-13])
    if args['pixel_cloud_file'][-17:-14] == 'left' or args['pixel_cloud_file'][-17:-14] == 'Left' :
        left_or_right = "L"
    else:
        left_or_right = "R"
    # General tile reference
    tile_ref = str(np.int(np.floor(np.mean(latitude)))).zfill(2) + nord_or_south + "-" + left_or_right
                    
    
    presumming_factor = 1
    fill_value = -9990000000.
    
    final_pixc.setncattr('NCProperties', "pixel cloud from 1051 version")
    final_pixc.setncattr('Conventions', "CF-1.7")
    final_pixc.setncattr('title', "Level 2 Pixel Clould Data Product")
    final_pixc.setncattr('institution', "JPL")
    final_pixc.setncattr('source', "Ka-band radar interferometer")
    final_pixc.setncattr('history', "None")
    final_pixc.setncattr('mission_name', "SWOT")
    final_pixc.setncattr('references', "None")
    final_pixc.setncattr('reference_document', "None")
    final_pixc.setncattr('contact', "None")
    final_pixc.setncattr('pass_number', sensor.getncattr('pass_number'))
    final_pixc.setncattr('cycle_number', sensor.getncattr('cycle_number'))
    #WARNING, TO BE CHANGED
    final_pixc.setncattr('tile_name', 0_0_0)
    final_pixc.setncattr('wavelength', 0.008385803020979)
    final_pixc.setncattr('near_range', np.float(pixc.getncattr('near_range')))
    final_pixc.setncattr('nominal_slant_range_spacing', np.float(pixc.getncattr('range_spacing')))
    final_pixc.setncattr('start_time', pixc.getncattr('time_coverage_start'))
    final_pixc.setncattr('stop_time', pixc.getncattr('time_coverage_end'))
    final_pixc.setncattr('polarization', "None")
    final_pixc.setncattr('transmit_antenna', "plus_y")
    final_pixc.setncattr('processing_beamwidth', "0LL")
    final_pixc.setncattr('ephemeris', "0LL")
    final_pixc.setncattr('yaw_flip', "0LL")
    final_pixc.setncattr('hpa_cold', "0LL")
    final_pixc.setncattr('inner_first_lat', latitude[np.argmin(latitude)])
    final_pixc.setncattr('inner_first_lon', longitude[np.argmin(latitude)])
    final_pixc.setncattr('inner_last_lat', latitude[np.argmin(longitude)])
    final_pixc.setncattr('inner_last_lon', longitude[np.argmin(longitude)])
    final_pixc.setncattr('outer_first_lat', latitude[np.argmax(latitude)])
    final_pixc.setncattr('outer_first_lon', longitude[np.argmax(latitude)])
    final_pixc.setncattr('outer_last_lat', latitude[np.argmax(longitude)])
    final_pixc.setncattr('outer_last_lon', longitude[np.argmax(longitude)])
    final_pixc.setncattr('slc_first_line_index_in_tvp', "None")
    final_pixc.setncattr('slc_last_line_index_in_tvp', "None")
    final_pixc.setncattr('xref_input_l1b_hr_slc', "None")
    final_pixc.setncattr('xref_static_karin_cal_file', "None")
    final_pixc.setncattr('xref_ref_dem_file', "None")
    final_pixc.setncattr('xref_water_mask_file', "None")
    final_pixc.setncattr('xref_static_geophys_file', "None")
    final_pixc.setncattr('xref_dynamic_geophys_file', "None")
    final_pixc.setncattr('xref_int_lr_xover_cal_file', "None")
    final_pixc.setncattr('xref_l2_hr_pixc_config_parameters_file', "None")
    final_pixc.setncattr('ellipsoid_semi_major_axis', "None")
    final_pixc.setncattr('ellipsoid_flattening', "None")
 
    pixc_group.setncattr('description', "cloud of geolocated interferogram pixels")
    pixc_group.setncattr('interferogram_size_range', 3500)
    pixc_group.setncattr('interferogram_size_azimuth', nb_tvp_record)
    pixc_group.setncattr('looks_to_efflooks', 1.75)
 

    for v in var_pixc_names: 
        fill_value = -9990000000.
        
        if v[0] == "interferogram":
            pixc_group.createVariable(v[0], v[1], ["record","depth"], fill_value = fill_value)
            pixc_group.variables[v[0]][:,0] = pixc.variables['ifgram_real'][:]
            pixc_group.variables[v[0]][:,1] = pixc.variables['ifgram_imag'][:]
        else:
            pixc_group.createVariable(v[0], v[1], "record", fill_value = fill_value)        
            if v[2] != None:
                pixc_group.variables[v[0]][:] = pixc.variables[v[2]][:]
            if v[0] == 'illumination_time':
                illumination_time = np.zeros(nb_pixels)
                azimuth_index_tmp = pixc.variables['azimuth_index']
                time_tmp = sensor.variables['time']
                for i in range(nb_pixels):
                    if i%1000 == 0:
                        print(i, "/", nb_pixels)
                        
                    illumination_time[i] = time_tmp[np.int(azimuth_index_tmp[i]*presumming_factor)]
             
                pixc_group.variables[v[0]][:] = illumination_time
                    
    ## Temporary fix for SAM version...
    if np.int(np.amax(pixc_group.variables['classification'][:])) == 1:
        print("Warning, only 1 value for classification, assuming SAM version, refactoring to 4")
        pixc_group.variables['classification'][:] = 4
    
    
    tvp_group = final_pixc.createGroup("tvp") 
    tvp_group.createDimension("record", nb_tvp_record)

    for v in var_tvp_names: 
        fill_value = -9990000000.
        tvp_group.createVariable(v[0], v[1], "record", fill_value = fill_value)        
        if v[2] != None:
            tvp_group.variables[v[0]][:] = sensor.variables[v[2]][:]
        

               
        
    pixc.close()
    sensor.close()
    final_pixc.close()

    dir_path = os.path.split(abspath(args['output_final_pixel_cloud']))[0]
    write_annotation_file(os.path.join(dir_path, "pixc-annotation.rdf"), abspath(args['output_final_pixel_cloud'])) 
    
if __name__ == "__main__":
    main()
