'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


import numpy as np
import matplotlib.pyplot as plt
from lib.my_variables import RAD2DEG
from scipy.spatial import cKDTree
import lib.my_api as my_api
from copy import deepcopy

def get_tiles_from_orbit(my_attributes, orbit_number):
    
    # Retrieve the tile database (pass_number/tile_number/nadir_lon/nadir_lat/nadir_heading)
    tile_db = my_attributes.tile_database
    
    # Subset the tile DB to the portion related to the orbit number
    #tile_db_orbit = tile_db[np.where(tile_db[:,0] == IN_orbit_number)[0],:]
    
    # ~ tmp_orbit_number = orbit_number - 331  # Pass 1 in tile database file = pass 332 in last KML file (sept2015-v2)
    # ~ if tmp_orbit_number < 1:
        # ~ tmp_orbit_number += 584

    tmp_orbit_number = orbit_number
        
    tile_db_orbit = tile_db[np.where(tile_db[:, 0] == tmp_orbit_number)[0], :]
    
    # Shift the lon/lat coordinates to compute center for each tile
    tile_db_orbit_center = np.copy(tile_db_orbit)
    for i in range(len(tile_db_orbit_center[:, 2])-1):
        tile_db_orbit_center[i, 2] = 0.5*(tile_db_orbit[i,2]+tile_db_orbit[i+1,2])
        tile_db_orbit_center[i, 3] = 0.5*(tile_db_orbit[i,3]+tile_db_orbit[i+1,3])
    # Compute the indices of nadir_lat_min and nadir_lat_max
    
    nadir_lat_argmin = int(np.argmin(my_attributes.lat*RAD2DEG))
    nadir_lat_argmax = int(np.argmax(my_attributes.lat*RAD2DEG))
    # Get long and lat in degrees, associated to nadir_min_lat
    nadir_lat_deg_min = my_attributes.lat[nadir_lat_argmin]*RAD2DEG
    nadir_lon_deg_min = my_attributes.lon[nadir_lat_argmin]*RAD2DEG
    # Get long and lat in degrees, associated to nadir_max_lat
    nadir_lat_deg_max = my_attributes.lat[nadir_lat_argmax]*RAD2DEG
    nadir_lon_deg_max = my_attributes.lon[nadir_lat_argmax]*RAD2DEG
    
    # Construct the kd-tree for quick nearest-neighbor lookup        
    tree = cKDTree(tile_db_orbit_center[:, 2:4])
    
    # Retrieve index of tile_db_orbit the nearest of nadir_min_lat
    ind_min = tree.query([nadir_lat_deg_min, nadir_lon_deg_min])
    # Retrieve index of tile_db_orbit the nearest of nadir_max_lat
    ind_max = tree.query([nadir_lat_deg_max, nadir_lon_deg_max])
    tile_db_orbit_cropped = tile_db_orbit_center[max(0, min(ind_max[1], ind_min[1])-1):min(len(tile_db_orbit_center), max(ind_max[1], ind_min[1])+2), :]
    vect_lat_lon_db_cropped = np.zeros([max(0, tile_db_orbit_cropped.shape[0]-1), 2])
    
    
    for i in range(max(0, tile_db_orbit_cropped.shape[0]-1)):
        vect_lat_lon_db_cropped[i,0] = tile_db_orbit_cropped[i+1, 2]-tile_db_orbit_cropped[i, 2]
        vect_lat_lon_db_cropped[i,1] = tile_db_orbit_cropped[i+1, 3]-tile_db_orbit_cropped[i, 3]
        nb_az_traj = max(nadir_lat_argmax, nadir_lat_argmin)- min(nadir_lat_argmax, nadir_lat_argmin) + 1
        tile_values = np.zeros(nb_az_traj, int)

    for i in range(min(nadir_lat_argmax, nadir_lat_argmin)-1, max(nadir_lat_argmax, nadir_lat_argmin)+1):
        dist = np.abs(((my_attributes.lat[i]*RAD2DEG-tile_db_orbit_cropped[:-1, 2])*vect_lat_lon_db_cropped[:, 0] + (my_attributes.lon[i]*RAD2DEG-tile_db_orbit_cropped[:-1, 3])*vect_lat_lon_db_cropped[:, 1])/np.sqrt(vect_lat_lon_db_cropped[:, 0]**2+vect_lat_lon_db_cropped[:, 1]**2))
        tile_values[i] = tile_db_orbit_cropped[np.argmin(dist),1]

    tile_list = np.unique(tile_values)

    my_api.printInfo("[my_tiling] [get_tiles_from_orbit] Simulation over tiles number: %s" % str(tile_list))
    
    return tile_values, tile_list


def crop_orbit(my_attributes, tile_values, tile_number, tropo_map_rg_az):

    my_new_attributes = deepcopy(my_attributes)

    my_api.printInfo("[my_tiling] [crop_orbit] == Dealing with tile number %03d" % tile_number)
    nadir_az = np.where(tile_values == tile_number)[0]
    my_api.printInfo("[my_tiling] [crop_orbit] nadir az contains %d pixels with %d and %d min and max value " %(nadir_az.size, min(nadir_az), max(nadir_az)))

    nb_pix_overlap_begin = 50
    nb_pix_overlap_end = 50

    if min(nadir_az) > nb_pix_overlap_begin:
        add_nadir = np.arange(min(nadir_az)-1-nb_pix_overlap_begin, min(nadir_az)-1)
        nadir_az = np.concatenate((nadir_az, add_nadir))
    else :
        nb_pix_overlap_begin = min(nadir_az)
        add_nadir = np.arange(0, min(nadir_az) - 1)
        nadir_az = np.concatenate((nadir_az, add_nadir))
    nadir_az = np.sort(nadir_az)

    if max(nadir_az) < len(my_attributes.orbit_time) -1 - nb_pix_overlap_end:
        add_nadir = np.arange(max(nadir_az)+1, max(nadir_az)+1+nb_pix_overlap_end)
        nadir_az = np.concatenate((nadir_az, add_nadir))
    else :
        nb_pix_overlap_end = len(my_attributes.orbit_time) - 1 - max(nadir_az)
        add_nadir = np.arange(max(nadir_az)+1, max(nadir_az)+1+ nb_pix_overlap_end)
        nadir_az = np.concatenate((nadir_az, add_nadir))
    my_api.printInfo("[my_tiling] [crop_orbit] nadir az contains %d pixels " % nadir_az.size)
    my_api.printInfo("[my_tiling] [crop_orbit] Tile contains %d and %d of overlaping azimuth pixel in the begining and the end of the tile" %(nb_pix_overlap_begin, nb_pix_overlap_end))

    nadir_az = np.sort(nadir_az)

    my_new_attributes.nb_pix_overlap_begin = nb_pix_overlap_begin
    my_new_attributes.nb_pix_overlap_end = nb_pix_overlap_end

    # Get azimuth indices corresponding to this integer value of latitude
    az_min = np.sort(nadir_az)[0]  # Min azimuth index, to remove from tile azimuth indices vector
    az_max = np.sort(nadir_az)[-1]  # Min azimuth index, to remove from tile azimuth indices vector
    my_api.printInfo("[my_tiling] [crop_orbit] = %d pixels in azimuth (index %d put to 0)" % (nadir_az.size, az_min))

    # Cropping orbit to only simulate tile area
    my_new_attributes.orbit_time = (my_attributes.orbit_time[nadir_az])
    my_new_attributes.x = my_attributes.x[nadir_az]
    my_new_attributes.y = my_attributes.y[nadir_az]
    my_new_attributes.z = my_attributes.z[nadir_az]

    my_new_attributes.lon = (my_attributes.lon[nadir_az])
    my_new_attributes.lon_init = (my_attributes.lon[nadir_az])
    
    my_new_attributes.lat = (my_attributes.lat[nadir_az])
    my_new_attributes.lat_init = (my_new_attributes.lat_init[nadir_az])

    my_new_attributes.heading = (my_attributes.heading[nadir_az])
    my_new_attributes.heading_init = (my_attributes.heading[nadir_az])

    my_new_attributes.alt = (my_attributes.alt[nadir_az])

    my_new_attributes.cosphi_init = my_attributes.cosphi_init[nadir_az]
    my_new_attributes.sinphi_init = my_attributes.sinphi_init[nadir_az]
    my_new_attributes.costheta_init = my_attributes.costheta_init[nadir_az]
    my_new_attributes.sintheta_init = my_attributes.sintheta_init[nadir_az]
    my_new_attributes.cospsi_init = my_attributes.cospsi_init[nadir_az]
    my_new_attributes.sinpsi_init = my_attributes.sinpsi_init[nadir_az]

    my_new_attributes.tile_number = tile_number

    if  tropo_map_rg_az is None:
        my_api.printInfo("[my_tiling] [crop_orbit] = Tropo field not applied")
        my_new_attributes.tropo_map_rg_az = None
    else:
        my_new_attributes.tropo_map_rg_az = tropo_map_rg_az[az_min:az_max,:]


    return my_new_attributes
