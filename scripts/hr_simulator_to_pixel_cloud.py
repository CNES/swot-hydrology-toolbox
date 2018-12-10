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

    var_pixc_names = [['azimuth_index', int, 'azimuth_index'], ['range_index', int, 'range_index'], ['x_factor_left', np.float32, None],  ['x_factor_right', np.float32, None], \
    ['pixel_area', np.float64, 'pixel_area'], ['incidence_angle', np.float64, None], ['classification', np.int, 'classification'], \
    ['continuous_classification', np.float32, 'continuous_classification'], ['ifgram', np.float32, None], ['power_left', np.float32, 'power_left'], ['power_right', np.float32, 'power_right'], \
    ['coherent_power', np.float64, 'coherent_power'], ['num_rare_looks', np.float32, None], ['latitude', np.float64, 'latitude_medium'], ['longitude', np.float64, 'longitude_medium'], \
    ['height', np.float64, 'height_medium'], ['cross_track', np.float64, 'cross_track_medium'], ['phase_noise_std', np.float32, 'dphase_medium'], ['dlatitude_dphase', np.float64, None], \
    ['dlongitude_dphase', np.float64, None], ['dheight_dphase', np.float64, 'dheight_dphase_medium'], \
    ['illumination_time', np.float64, None], ['num_med_looks', np.float32, None], ['sigma0', np.float32, None], ['regions', np.float64, None]]
  
    var_tvp_names  = [['time', np.float64, 'time'],['latitude', np.float64, 'latitude'],['longitude', np.float64, 'longitude'],['height', np.float64, 'altitude'],['heading', np.float64, 'heading'],['x', np.float64, 'x'],\
    ['y', np.float64, 'y'],['z', np.float64, 'z'],['near_range', np.float64, None],['baseline_left_x', np.float64, 'baseline_left_x'],['baseline_left_y', np.float64, 'baseline_left_y'],\
    ['baseline_left_z', np.float64, 'baseline_left_z'],['baseline_right_x', np.float64, 'baseline_right_x'],['baseline_right_y', np.float64, 'baseline_right_y'],['baseline_right_z', \
    np.float64, 'baseline_right_z'], ['vx', np.float64, 'velocity_unit_x'],['vy', np.float64, 'velocity_unit_y'],['vz', np.float64, 'velocity_unit_z'],['ref_leverarm_x', np.float64, None],\
    ['ref_leverarm_y', np.float64, None],['ref_leverarm_z', np.float64, None],['sec_leverarm_x', np.float64, None],['sec_leverarm_y', np.float64, None],['sec_leverarm_z', np.float64, None]]
  

  
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
    pixc_group.setncattr('wavelength', 0.008385803020979)
    pixc_group.setncattr('near_range', np.float(pixc.getncattr('near_range')))
    pixc_group.setncattr('range_spacing', np.float(pixc.getncattr('range_spacing')))
    pixc_group.setncattr('azimuth_spacing', np.float(pixc.getncattr('azimuth_spacing')))
    pixc_group.setncattr('noise_power_left', -116.8458)
    pixc_group.setncattr('noise_power_right', -116.8458)
    pixc_group.setncattr('start_time', pixc.getncattr('time_coverage_start'))
    pixc_group.setncattr('stop_time', pixc.getncattr('time_coverage_end'))
    pixc_group.setncattr('pass_number', sensor.getncattr('pass_number'))
    pixc_group.setncattr('cycle_number', sensor.getncattr('cycle_number'))
    pixc_group.setncattr('tile_ref', tile_ref)
    pixc_group.setncattr('inner_first_lat', latitude[np.argmin(latitude)])
    pixc_group.setncattr('inner_first_lon', longitude[np.argmin(latitude)])
    pixc_group.setncattr('inner_last_lat', latitude[np.argmin(longitude)])
    pixc_group.setncattr('inner_last_lon', longitude[np.argmin(longitude)])
    pixc_group.setncattr('outer_first_lat', latitude[np.argmax(latitude)])
    pixc_group.setncattr('outer_first_lon', longitude[np.argmax(latitude)])
    pixc_group.setncattr('outer_last_lat', latitude[np.argmax(longitude)])
    pixc_group.setncattr('outer_last_lon', longitude[np.argmax(longitude)])
    pixc_group.setncattr('description', fill_value)
    pixc_group.setncattr('nr_pixels', 3500)
    pixc_group.setncattr('nr_lines', np.int(nb_tvp_record/presumming_factor))
    pixc_group.setncattr('looks_to_efflooks', 0)
    
    
    for v in var_pixc_names: 
        fill_value = -9990000000.
        if v[0] == "ifgram":
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
