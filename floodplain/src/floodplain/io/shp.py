# -*- coding: utf8 -*-
'''
Shaepile writer

Copyright (c) 2018 CNES. All rights reserved.
'''

import fiona
import fiona.crs
import shapely
from shapely.geometry import Point, Polygon, mapping, shape
from collections import OrderedDict
from typing import List, Tuple
import geopandas as gpd

def from_file(filename: str, wse_flag : bool = False, wse_name : str = 'elevation') -> List[Polygon]:
    '''
    Read to shapefile (lat/lon WGS84)

    :param filename: File path
    :return polygons
    ''' 
    polygons = []
    wse = []
    with fiona.open(filename) as c:
        for element in iter(c):
            #use the shape function of Shapely
            polygons.append(shape(element['geometry']))
            if wse_flag:
                wse.append(element['properties'][wse_name])
    if wse_flag:
        return polygons, wse
    else:
        return polygons


def points_to_file(filename: str, points: List[Tuple]):
    '''
    Write to shapefile (lat/lon WGS84)

    :param filename: File path
    :param points: List of points (each pooint is represented 
                   by the following tuple (lon, lat, height)
    ''' 
    # Driver
    driver = "ESRI Shapefile"
    crs = fiona.crs.from_epsg(4326) # WGS84
    # Schema
    schema = {'properties': OrderedDict([('longitude', 'float:24.15'), 
                                     ('latitude', 'float:24.15'), 
                                     ('elevation', 'float:24.15'),]),
              'geometry': 'Point'}

    with fiona.open(filename,'w', driver=driver, crs=crs, schema=schema) as c:
        for point in points:
            prop = {'longitude': float(point[0]),
                    'latitude': float(point[1]),
                    'elevation': float(point[2]),
                   }
            c.write({'geometry': shapely.geometry.mapping(Point(point[0],point[1])), 'properties': prop})


def gdf_to_file(filename: str, gdf: gpd.GeoDataFrame, index: bool = False):
    '''
    Write to shapefile (lat/lon WGS84)

    :param filename: File path
    :param gdf: DataFrame
    ''' 
    # Driver
    driver = "ESRI Shapefile"
    # Projection
    prj = 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]'
    # Schema
    if index:
        schema = {'properties': OrderedDict([('longitude', 'float:24.15'), 
                                 ('latitude', 'float:24.15'), 
                                 ('elevation', 'float:24.15'),
                                 ('range', 'int'),
                                 ('azimuth', 'int'),]),
                  'geometry': 'Point'}
        gdf[['geometry','longitude','latitude','elevation','range','azimuth']].to_file(filename,driver=driver, schema=schema, crs_wkt=prj)
    else:
        schema = {'properties': OrderedDict([('longitude', 'float:24.15'), 
                                 ('latitude', 'float:24.15'), 
                                 ('elevation', 'float:24.15'),]),
                  'geometry': 'Point'}
        gdf[['geometry','longitude','latitude','elevation']].to_file(filename,driver=driver, schema=schema, crs_wkt=prj)


def polygons_to_file(filename: str, polygons: List[Polygon], wse: List[float] = None):
    '''
    Write to shapefile (lat/lon WGS84)

    :param filename: File path
    :param polygon: List of polygons
    ''' 
    # Driver
    driver = "ESRI Shapefile"
    crs = fiona.crs.from_epsg(4326) # WGS84
    # Schema
    if wse != None:
        schema = {'geometry': 'Polygon','properties': {'id': 'int', 'elevation': 'float'}}
    else:
        schema = {'geometry': 'Polygon','properties': {'id': 'int'}}
        
        
    with fiona.open(filename, 'w', driver=driver, crs=crs, schema=schema) as c:

        for i,polygon in enumerate(polygons):
            if wse != None:
                c.write({'geometry': mapping(polygon),'properties': {'id': i, 'elevation':wse[i]}})
            else:
                c.write({'geometry': mapping(polygon),'properties': {'id': i}})

            

