# -*- coding: utf8 -*-
'''
Library provided methods to manipulate polygon

Copyright (c) 2018, CNES
'''

from typing import TypeVar, List
from shapely.geometry import Polygon,LinearRing, MultiPolygon

def extract_polygon_boundaries(polygons: List[Polygon]) -> List[LinearRing]:
    '''
    Extract boundaries from a list of polygons
   
    :param polygons: List of polygons
    :return: Boundaries
    '''
    boundaries = []

    if polygons.type == 'MultiPolygon':
        for polygon in polygons:
            boundaries.append(polygon.exterior)
            for interior in polygon.interiors:
                boundaries.append(interior)
        return boundaries

    elif polygons.type == 'Polygon':
        boundaries.append(polygons.exterior)
        for interior in polygons.interiors:
            boundaries.append(interior)
        return boundaries
    else:
        print("ERROR : Polygon not in good type")

def filter_polygon(polygon: Polygon, threshold: float = 1.0) -> Polygon:
    '''
    Remove small interior of a polygon
    
    :param polygon: Polygon
    :param threshold: Threshold value
    :return: filtered polygon
    '''
    ext = polygon.exterior
    ints = [interior for interior in polygon.interiors if Polygon(interior).area > threshold]
    return Polygon(ext, ints)


def filter_polygons(polygons: List[Polygon], threshold: float = 1.0) -> List[Polygon]:
    '''
    Remove small polygons (interior/extrerior) in a list of polygons
   
    :param polygons: List of polygons
    :param threshold: Threshold value
    :return: filtered list of polygons
    '''
    filtered_polygons = []
    
    if polygons.type == 'MultiPolygon':
        # Loop over polygon
        for polygon in polygons:
            # Keep only polygon with an aera greater than the threshold.
            if Polygon(polygon.exterior).area > threshold:
                #Filter polygon interiors
                filtered_polygons.append(filter_polygon(polygon,threshold))
        return filtered_polygons
    elif polygons.type == 'Polygon':
        # Keep only polygon with an aera greater than the threshold.
        if Polygon(polygons.exterior).area > threshold:
            #Filter polygon interiors
            filtered_polygons.append(filter_polygon(polygons,threshold))
        return filtered_polygons
    else:
        print("ERROR : Polygon not in good type")

def print_stats(polygons):
    '''
    Print statistics for a list of polygons
   
    :param polygons: List of polygons
    '''
    print(f"Number of polygons = {len(polygons)}")
    for i,polygon in enumerate(polygons):
        print(f"Polygon {i}")
        print(f"  Informations:")
        print(f"    - area = {polygon.area}")
        print(f"    - area (ext) = {Polygon(polygon.exterior).area}")
        print(f"    - nb interiors = {len(polygon.interiors)}")
        for j, interior in enumerate(polygon.interiors):
            print(f"    - area {j} (int) = {Polygon(interior).area}")
    


    




