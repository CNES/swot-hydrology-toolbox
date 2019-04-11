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
import shutil
import lib.my_api as my_api
import my_rdf

from lib.my_variables import RAD2DEG, DEG2RAD, GEN_APPROX_RAD_EARTH
from lib.tropo_module import Tropo_module
from lib.roll_module import Roll_module
import mathematical_function as math_fct
from write_polygons import orbitAttributes, project_array
import cnes.modules.geoloc.lib.geoloc as my_geoloc

    
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('pixel_cloud_file', type=str)
    parser.add_argument('output_final_pixel_cloud', default= None, type=str)
    parser.add_argument('rdf', default= None, type=str)

    args = vars(parser.parse_args())
    shutil.copyfile(abspath(args['pixel_cloud_file']), abspath(args['output_final_pixel_cloud']))
    
    rdf = my_rdf.myRdfReader(args['rdf'])
    parameters = rdf.parameters
    
    pixc = Dataset(abspath(args['pixel_cloud_file']), 'r')
    final_pixc = Dataset(abspath(args['output_final_pixel_cloud']), 'a')

    az = pixc.groups['pixel_cloud'].variables['azimuth_index'][:]
    col = pixc.groups['pixel_cloud'].variables['range_index'][:]
    lat = pixc.groups['pixel_cloud'].variables['latitude'][:]
    lon = pixc.groups['pixel_cloud'].variables['longitude'][:]
    height = pixc.groups['pixel_cloud'].variables['height'][:]
    pixel_cloud_time = pixc.groups['pixel_cloud'].variables['illumination_time'][:]
    cross_track = pixc.groups['pixel_cloud'].variables['cross_track'][:]
    
    orbit_time = pixc.groups['tvp'].variables['time'][:]
    
    cycle_number = pixc.getncattr('cycle_number')
    pass_number = pixc.getncattr('pass_number')
    swath_side = pixc.getncattr('swath_side')
    if swath_side == 'L':
        swath_side = 'Left'
    if swath_side == 'R':
        swath_side = 'Right'
                
    near_range = pixc.getncattr('near_range')
    nominal_slant_range_spacing = pixc.getncattr('nominal_slant_range_spacing')
    
    ri = near_range + nominal_slant_range_spacing*col

    delta_h = 0.
    
    if parameters['Tropo model'] == 'gaussian':

        my_api.printInfo("Applying wet tropo gaussian model")
        tropo = Tropo_module(parameters['Tropo model'])
        tropo_error = tropo.calculate_tropo_error_gaussian(az, col, parameters['Tropo error stdv'], parameters['Tropo error mean'],  parameters['Tropo error correlation']) 
        delta_h += tropo_error
        
    elif parameters['Tropo model'] == 'map':
        
        my_api.printInfo("Applying wet tropo map gaussian model")
        tropo = Tropo_module(parameters['Tropo model'])
        tropo_error = tropo.calculate_tropo_error_map(np.nanmean(lat), az, col, parameters['Tropo error map file'], parameters['Tropo error correlation'])
        delta_h += tropo_error
        
    else:
        my_api.printInfo("No tropo model applied")
            
    
    if parameters['roll_repository_name'] != None:
            
        my_api.printInfo("Applying roll residual error")

        roll = Roll_module(parameters['roll_repository_name'])
        roll.get_roll_file_associated_to_orbit_and_cycle(pass_number, cycle_number, swap = False, delta_time = -1541.907908)
        roll.interpolate_roll_on_sensor_grid(orbit_time)
        
        # Apply roll for each pixel
        roll.interpolate_roll_on_pixelcloud(orbit_time, pixel_cloud_time, cross_track)
        delta_h_roll = (roll.roll1_err_cloud)
        delta_h += delta_h_roll
        
    else:
        my_api.printInfo("No roll error applied")

    print(delta_h)   


    IN_attributes = orbitAttributes()
    
    IN_attributes.lat = pixc.groups['tvp'].variables['latitude'][:] * DEG2RAD
    IN_attributes.lon = pixc.groups['tvp'].variables['longitude'][:] * DEG2RAD
    IN_attributes.heading_init = pixc.groups['tvp'].variables['velocity_heading'][:] * DEG2RAD
    IN_attributes.alt = pixc.groups['tvp'].variables['altitude'][:] 
    
    IN_attributes.vx = pixc.groups['tvp'].variables['vx'][:] 
    IN_attributes.vy = pixc.groups['tvp'].variables['vy'][:] 
    IN_attributes.vz = pixc.groups['tvp'].variables['vz'][:]
        
    if parameters['noisy_geoloc']:
        
        coordinates_sensor = np.dstack((IN_attributes.lat[az]*RAD2DEG, IN_attributes.lon[az]*RAD2DEG, IN_attributes.alt[az]))[0]
        xyz_sensor = project_array(coordinates_sensor, srcp='latlon', dstp='geocent')
        coordinates_pixc = np.dstack((lat, lon, height))[0]     
        xyz_pixc = project_array(coordinates_pixc, srcp='latlon', dstp='geocent')
        vxyz_sensor = np.dstack((IN_attributes.vx[az], IN_attributes.vy[az], IN_attributes.vz[az]))[0]
        
        ri = np.sqrt((xyz_sensor[:,0]-xyz_pixc[:,0])**2+(xyz_sensor[:,1]-xyz_pixc[:,1])**2+(xyz_sensor[:,2]-xyz_pixc[:,2])**2)
        
                
        p_final, p_final_llh, h_mu, (iter_grad,nfev_minimize_scalar) = my_geoloc.pointcloud_height_geoloc_vect(xyz_pixc, height,
                                                                                                           xyz_sensor,
                                                                                                           vxyz_sensor,
                                                                                                           ri, height + delta_h, 
                                                                                                           recompute_Doppler=True, recompute_R=True, verbose=False, 
                                                                                                           max_iter_grad=1, height_goal=1.e-3, safe_flag=True)
             
        final_pixc.groups['pixel_cloud'].variables['latitude'][:] = p_final_llh[:,0]
        final_pixc.groups['pixel_cloud'].variables['longitude'][:] = p_final_llh[:,1]
        final_pixc.groups['pixel_cloud'].variables['height'][:] = p_final_llh[:,2]
     
    else:
        final_pixc.groups['pixel_cloud'].variables['height'][:] = height + delta_h
                                                                                                               
    pixc.close()
    final_pixc.close()
    
if __name__ == "__main__":
    main()
