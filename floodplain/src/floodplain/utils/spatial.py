# -*- coding: utf8 -*-
'''
Library provided extra methods

Copyright (c) 2018, CNES
'''

import numpy as np
from shapely.geometry import Point, Polygon
import utm
from functools import partial
import pyproj
from shapely.ops import transform
from typing import Tuple, List


# Global variables
GEN_RAD_EARTH_EQ = 6378137.0 # Radius (in meters) of the Earth model (WGS84 ellipsoid) at the equator
GEN_RAD_EARTH_POLE = 6356752.31425 # Radius (in meters) of the Earth model to the pole 
GEN_APPROX_RAD_EARTH = (2*GEN_RAD_EARTH_EQ + GEN_RAD_EARTH_POLE)/3 # Radius (in meters) of the sphere equivalent to ellipsoid


def convert_to_utm_coords(latitude: float,longitude: float) -> Tuple[float,float]:
    '''
    Convert latitude/longitude to utm coordinates

    :param latitude: Latitude
    :param longitude: Longitude
    :return: Point in utm coordinates
    '''
    return utm.from_latlon(latitude,longitude)[0:2]


def convert_to_utm_point(latitude: float,longitude: float, precision = None) -> Point:
    '''
    Convert latitude/longitude to utm coordinates

    :param latitude: Latitude
    :param longitude: Longitude
    :return: Point in utm coordinates
    '''
    x,y = convert_to_utm_coords(latitude, longitude) 
    if precision is not None:
        x=round(x,precision)
        y=round(y,precision)
    return Point(x,y)


def convert_polygon_utm_to_latlon(polygon: Polygon, zone_number: int, north: bool = True):
    '''
    Convert a polygon in utm coordinates to lat/lon coordinates

    :param polygon: Polygon in utm coodinates
    :return: Polygon in lat/lon coordinates
    '''
    pos = "north" if north else "south"
    proj = f"+proj=utm +zone={zone_number} +{pos} +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    
    project = partial(pyproj.transform,
                      pyproj.Proj(proj), # source coordinate system
                      pyproj.Proj(init='epsg:4326')) # destination coordinate system (lat/lon WGS84)
    return transform(project, polygon)  # apply projection


def convert_polygons_utm_to_latlon(polygons: List[Polygon], zone_number: int, north: bool = True):
    '''
    Convert a list of polygons in utm coordinates to lat/lon coordinates

    :param polygon: Polygon in utm coodinates
    :return: Polygon in lat/lon coordinates
    '''
    return [convert_polygon_utm_to_latlon(polygon,zone_number,north) for polygon in polygons]


def compute_binary_mask(IN_sizeX, IN_sizeY, IN_X, IN_Y):
    '''
    Creates a 2D binary matrix from X and Y 1D vectors
    i.e. for each ind: mat[IN_X[ind],IN_Y[ind]] = 1

    :param IN_sizeX: number of pixels in X dimension
    :type IN_sizeX: int
    :param IN_sizeY: number of pixels in Y dimension
    :type IN_sizeY: int
    :param IN_X: X indices of "1" pixels
    :type IN_X: 1D vector of int
    :param IN_Y: Y indices of "1" pixels
    :type IN_Y: 1D vector of int

    :return: 2D matrix with "1" for each (IN_X_i, IN_Y_i) and 0 elsewhere
    :rtype: 2D binary matrix of int 0/1
    '''
    
    # 0 - Deal with exceptions
    # 0.1 - Input vectors size must be the same
    if IN_X.size != IN_Y.size:
        raise ValueError("computeBinMat(IN_X, IN_Y) : IN_X and IN_Y must be the same size ; currently : IN_X = %d and IN_Y = %d" % ( IN_X.size , IN_Y.size ))
    else:
        nb_pts = IN_X.size
    
    # 0.2 - max(X) < IN_sizeX
    if np.max(IN_X) >= IN_sizeX:
        raise ValueError("computeBinMat(IN_X, IN_Y) : elements of IN_X must be less than IN_sizeX")
    # 0.3 - max(X) < IN_sizeX
    if np.max(IN_Y) >= IN_sizeY:
        raise ValueError("computeBinMat(IN_X, IN_Y) : elements of IN_Y must be less than IN_sizeY")

    # 1 - Init output binary image
    OUT_binIm = np.zeros((IN_sizeX, IN_sizeY))

    # 2 - Put 1 for every pixels defined by the input vectors
    OUT_binIm[IN_X, IN_Y] = 1

    return OUT_binIm


def deg2rad(IN_deg):
    '''
    Convert angles in degrees to radians
    
    :param IN_deg: angles to convert
    :type IN_deg: scalar or 1D-array
    
    :return angles in radians
    :rtype: same as input
    '''
    return IN_deg * np.pi / 180.


def rad2deg(IN_rad):
    '''
    Convert angles in radians to degrees
    
    :param IN_rad: angles to convert
    :type IN_rad: scalar or 1D-array
    
    :return angles in degrees
    :rtype: same as input
    '''
    return IN_rad * 180. / np.pi

def llh2xyz(IN_lon, IN_lat, IN_height, IN_flag_rad=True):
    '''
    Convert geographic coordinates (longitude, latitude, height) to cartesian coordinates (x, y, z) 
        
    :param IN_lon: longitude in degrees east
    :type IN_lon: scalar or 1D-array
    :param IN_lat: latitude in degrees north
    :type IN_lat: scalar or 1D-array
    :param IN_height: height in meters
    :type IN_height: scalar or 1D-array
    :param IN_flag_rad: =False if IN_lan and IN_lat are in degrees, =True if they are in radians (default)
    :type IN_flag_rad: boolean
    
    :return OUT_x, OUT_y, OUT_z: cartesian coordinates
    :type: scalars or 1D-array
    
    .. author: Damien DESROCHES & Alejandro BOHE - CNES DSO/SI/TR
    '''

    # 1 - Convert geodetic coordinates in radians
    if ( IN_flag_rad ):
        lat = IN_lat
        lon = IN_lon
    else:
        lat = deg2rad(IN_lat)
        lon = deg2rad(IN_lon)

    # 2 - Compute earth excentricity
    e = np.sqrt(GEN_RAD_EARTH_EQ ** 2 - GEN_RAD_EARTH_POLE ** 2) / GEN_RAD_EARTH_EQ

    # 3 - Compute earth radius for latitude lat
    Rn = GEN_RAD_EARTH_EQ / np.sqrt(1 - e**2 * np.sin(lat)**2 )

    # 4 - Compute cartesian coordinates
    OUT_x = (Rn + IN_height) * np.cos(lat) * np.cos(lon)
    OUT_y = (Rn + IN_height) * np.cos(lat) * np.sin(lon)
    OUT_z = (Rn * (1.-e**2) + IN_height) * np.sin(lat)

    return OUT_x, OUT_y, OUT_z


def xyz2llh(IN_x, IN_y, IN_z):
    '''
    Convert cartesian coordinates (x, y, z) to geographic coordinates (lat, lon, height) 
    
    :param IN_x: coordinate along x-axis
    :type IN_x: scalar or 1D-array
    :param IN_y: coordinate along y-axis
    :type IN_y: scalar or 1D-array
    :param IN_z: coordinate along z-axis
    :type IN_z: scalar or 1D-array
    
    :return OUT_long, OUT_lat, OUT_height: geographic coordinates (longitude, latitude, height) 
    :rtype: same as input (scalar or 1D-array)
    '''

    # 3.0 - Compute ellipsoide variables
    e = np.sqrt(GEN_RAD_EARTH_EQ ** 2 - GEN_RAD_EARTH_POLE ** 2) / GEN_RAD_EARTH_EQ # Ellipsoid excentricity
    f = 1 - np.sqrt(1 - e**2) # Flattening factor
    r = np.sqrt(IN_x**2 + IN_y**2 + IN_z**2) # Distance between center of earth and point
    d = np.sqrt(IN_x**2 + IN_y**2) # Distance between center of earth and point with latitude 0

    # 3.1 - Compute longitudes
    OUT_lon = np.arctan2(IN_y, IN_x)

    # 3.2 - Compute latitudes
    mu = np.arctan(IN_z/d * ( (1-f) + e**2 * GEN_RAD_EARTH_EQ/r))
    OUT_lat = np.arctan( (IN_z*(1-f)+e**2 * GEN_RAD_EARTH_EQ * np.sin(mu)**3) / ((1-f)*(d-e**2 * GEN_RAD_EARTH_EQ * np.cos(mu)**3)) )

    # 3.3 - Compute heights
    OUT_height = (d*np.cos(OUT_lat)) + (IN_z*np.sin(OUT_lat)) - GEN_RAD_EARTH_EQ*np.sqrt(1-(e**2)*np.sin(OUT_lat)**2)

    return rad2deg(OUT_lon), rad2deg(OUT_lat), OUT_height

