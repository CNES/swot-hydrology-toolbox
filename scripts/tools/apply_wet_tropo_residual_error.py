'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


from netCDF4 import Dataset
from os.path import abspath
import argparse
import numpy as np
import shutil
import my_rdf
import pyproj as pyproj


from lib.my_variables import RAD2DEG, DEG2RAD
from lib.tropo_module import Tropo_module
from lib.roll_module import Roll_module
from write_polygons import orbitAttributes
import cnes.modules.geoloc.lib.geoloc as my_geoloc
from datetime import datetime

    
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
    pixel_cloud_time = pixc.groups['pixel_cloud'].variables['illumination_time'][:].filled()
    cross_track = pixc.groups['pixel_cloud'].variables['cross_track'][:]
    
    orbit_time = pixc.groups['tvp'].variables['time'][:].filled()
    
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

    tropo = Tropo_module(parameters['Tropo model'], min(col), max(col), min(az), max(az), \
    float(parameters['Tropo error stdv']), float(parameters['Tropo error mean']), float(parameters['Tropo error correlation']), \
    parameters['Tropo error map file'])

    tropo.generate_tropo_field_over_pass(min(lat))

    tropo.apply_tropo_error_on_pixels(az, col)
    tropo_2d_field = tropo.tropo_2d_field
    delta_h += tropo_2d_field        
   
    utc_start_crossover = datetime(2021, 9, 1)
    utc_end_crossover = datetime(2022, 9, 1)
    utc_ref_simu = datetime(2000, 1, 1)
    utc_date_simu = datetime(2019, 2, 5)
    
    diff_between_date = (utc_start_crossover-utc_ref_simu).total_seconds()
    diff_between_simulation_and_reference_crossover = (utc_date_simu-utc_start_crossover).total_seconds()
    
    repeat_cycle_period = 1802697.1564
    nb_cycle_shift = diff_between_simulation_and_reference_crossover//repeat_cycle_period+1

    if diff_between_simulation_and_reference_crossover<0:
        print("Warning, Simulation date (",utc_date_simu,") is below crossover simulated start time : ", utc_start_crossover)
        
    if diff_between_simulation_and_reference_crossover > (utc_end_crossover-utc_start_crossover).total_seconds():
        print("Warning, Simulation date is (",utc_date_simu,") after crossover simulated start time : ", utc_end_crossover)      
    try:
        if parameters['roll_repository_name'] != None:
            print("Applying roll residual error")

            roll = Roll_module(parameters['roll_repository_name'])
            # ~ roll.get_roll_file_associated_to_orbit_and_cycle(pass_number, cycle_number, delta_time = diff_between_date)
            roll.get_roll_file_associated_to_orbit_and_cycle(pass_number, cycle_number, delta_time = diff_between_date+nb_cycle_shift*repeat_cycle_period)
            
            roll.interpolate_roll_on_sensor_grid(orbit_time)
            
            # Apply roll for each pixel
            roll.interpolate_roll_on_pixelcloud(orbit_time, pixel_cloud_time, cross_track)
            delta_h_roll = (roll.roll1_err_cloud)
            delta_h += delta_h_roll
            
        else:
            print("No roll error applied")
    except:
        print("Error during crossover application, No roll error applied")
        
    IN_attributes = orbitAttributes()
    
    IN_attributes.lat = pixc.groups['tvp'].variables['latitude'][:] * DEG2RAD
    IN_attributes.lon = pixc.groups['tvp'].variables['longitude'][:] * DEG2RAD
    IN_attributes.heading_init = pixc.groups['tvp'].variables['velocity_heading'][:] * DEG2RAD
    IN_attributes.alt = pixc.groups['tvp'].variables['altitude'][:] 
    
    IN_attributes.vx = pixc.groups['tvp'].variables['vx'][:] 
    IN_attributes.vy = pixc.groups['tvp'].variables['vy'][:] 
    IN_attributes.vz = pixc.groups['tvp'].variables['vz'][:]
        
        
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
                                                                                                       recompute_doppler=True, recompute_range=True, verbose=False, 
                                                                                                       max_iter_grad=1, height_goal=1.e-3)
         
    final_pixc.groups['pixel_cloud'].variables['latitude'][:] = p_final_llh[:,0]
    final_pixc.groups['pixel_cloud'].variables['longitude'][:] = p_final_llh[:,1]
    final_pixc.groups['pixel_cloud'].variables['height'][:] = p_final_llh[:,2]
 
                                                                                                               
    pixc.close()
    final_pixc.close()

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

    
if __name__ == "__main__":
    main()
