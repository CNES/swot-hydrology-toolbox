# -*- coding: utf8 -*-
'''
Create ply file from a list of pixel cloud files

Copyright (c) 2018, CNES
'''

import shapely.wkt

import logging
import mahotas as mh
import numpy as np
import geopandas as gpd
import pandas as pd
from typing import List
from scipy import spatial
import utm
from shapely.geometry import Polygon, LinearRing,Point, MultiPolygon
from random import randint
from shapely.ops import unary_union
from osgeo import ogr
from scipy.ndimage.morphology import binary_erosion

from floodplain.utils.spatial import compute_binary_mask, convert_to_utm_point
from floodplain.geom.alpha_shape import alpha_shape_with_cgal, alpha_shape_with_cascaded_union
from floodplain.constants import WATER_LABELS, WATER_LABEL, WATER_NEAR_LAND_LABEL, DARK_WATER_LABEL, WATER_NEAR_LAND_LOW_COH_LABEL, WATER_LOW_COH_LABEL

import math
import numpy as np
from osgeo import ogr
import pandas as pd
from random import randint
from scipy.spatial import distance, Delaunay
from scipy.ndimage.morphology import binary_erosion
import shapely.geometry as geometry
from shapely.geometry import Point, LineString, MultiPoint, MultiLineString, Polygon, MultiPolygon, LinearRing, GeometryCollection
from shapely.ops import cascaded_union
from shapely.ops import unary_union
from skimage.measure import find_contours

def extract_water_points(pixelcloud: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    '''
    Extract pixel cloud with water labels

    :param pixelcloud: Pixel cloud Dataframe

    :return points: Extracted water points 
    '''

    # Extract pixel cloud with water labels
    return pixelcloud.loc[(pixelcloud.classification == WATER_LABEL) |
                          (pixelcloud.classification == WATER_NEAR_LAND_LABEL) |
                          (pixelcloud.classification == DARK_WATER_LABEL) |
                          (pixelcloud.classification == WATER_NEAR_LAND_LOW_COH_LABEL) |
                          (pixelcloud.classification == WATER_LOW_COH_LABEL) ].copy()

def extract_contiguous_water_points(water: gpd.GeoDataFrame, 
                         range_max: int, 
                         azimuth_max: int,
                         threshold: float = 100000.0) -> gpd.GeoDataFrame:
    '''
    Extract water points  
    Steps :
      1) Create water mask
      2) Labelize regions 
      3) Compute regions area and remove small regions

    :param water: Water Dataframe
    :param range_max: Range size
    :param azimuth_max: Azimuth size
    :param threshold: Threshold for the keeping region area

    :return points: Extracted water points 
    '''

    # Create water mask
    water_mask = compute_binary_mask(azimuth_max,
                               range_max,
                               water['azimuth_index'].values,
                               water['range_index'].values)

    # Boundary conditions
    bc = np.ones((3,3))
    # Labelize regions
    labeled, nr_objects = mh.label(water_mask, Bc=bc)
    
    # Region size
    sizes = mh.labeled.labeled_size(labeled)
    df_sizes = pd.DataFrame(data=sizes[1:], columns=["size"], index=list(range(1,len(sizes))))
    # Compute region area
    water['region'] = water.apply(lambda row: labeled[row['azimuth_index'],row['range_index']], axis=1)
    df_sizes['area'] = water.groupby(['region'])['pixel_area'].sum()
    df_sizes.sort_values(['size'], ascending=False)
    # Keep regions with area superior to threshold
    keep_regions = df_sizes.loc[df_sizes.area > threshold].index.tolist()
    # Extract points corresponding to remaining regions
    return water.loc[water['region'].isin(keep_regions)].copy()
    
def remove_isolated_points(water: gpd.GeoDataFrame, slc_along_track_resolution: float, 
                           factor: float = 1.0,
                           nb_neighbors: int = 4) -> gpd.GeoDataFrame:
    '''
    Remove isolated points

    :param water: Water Dataframe

    :return points: Extracted water points 
    '''
    # Add utm coordiantes
    water['geometry_utm'] = water.apply(lambda row: convert_to_utm_point(row['latitude'],row['longitude'], precision=6), axis=1)    
    water['utm_x'] = water.apply(lambda row: row['geometry_utm'].x, axis=1)
    water['utm_y'] = water.apply(lambda row: row['geometry_utm'].y, axis=1)
    tree = spatial.cKDTree(np.array([ [pt.x,pt.y] for pt in water["geometry_utm"].values ])) 
    # Find all points within distance r of point(s) x
    water["neighbors"] = water.apply(lambda row: len(tree.query_ball_point([row['geometry_utm'].x,row['geometry_utm'].y], factor*np.sqrt((row['range_spacing'])**2+(slc_along_track_resolution)**2))), axis=1)
    # Keep points with at least 4 neighbors
    return water.loc[(water.neighbors > nb_neighbors)].copy()






def split_lake_image(in_range, in_azimuth, in_v_long, in_v_lat, nb_pix_max=1e4):
    """
    Split input image into list of sub-images, with 10 000 pixel maximum.

    :param in_range: range of points
    :type in_range: 1D-array of int
    :param in_azimuth: azimuth of points
    :type in_azimuth: 1D-array of int
    :param in_v_long: longitude of points
    :type in_v_long: 1D-array of float
    :param in_v_lat: latitude of points
    :type in_v_lat: 1D-array of float

    :return: list of tuple of splitted image in sub images
    :rtype: list of tulple(in_range, in_azimuth, in_v_long, in_v_lat)
    """

    if len(in_range) > nb_pix_max:
        margin = 10
        if np.max(in_range) - np.min(in_range) > np.max(in_azimuth) - np.min(in_azimuth):
            rg_half = int(np.max(in_range) - np.min(in_range)) / 2 + np.min(in_range)
            idx1 = np.where(in_range < rg_half + margin)
            idx2 = np.where(in_range > rg_half - margin)
        else:
            az_half = int(np.max(in_azimuth) - np.min(in_azimuth)) / 2 + np.min(in_azimuth)
            idx1 = np.where(in_azimuth < az_half + margin)
            idx2 = np.where(in_azimuth > az_half - margin)
        retour = split_lake_image(in_range[idx1], in_azimuth[idx1], in_v_long[idx1], in_v_lat[idx1], nb_pix_max) +\
                 split_lake_image(in_range[idx2], in_azimuth[idx2], in_v_long[idx2], in_v_lat[idx2], nb_pix_max)
    else:
        retour = [(in_range, in_azimuth, in_v_long, in_v_lat)]
    return retour

def get_concave_hull_from_radar_vectorisation(in_range, in_azimuth, in_v_long, in_v_lat, in_hull_method):
    """
    Compute the concave hull of a set of points using radar vectorisation, split image into sub-images if needed

    :param in_range: range of points
    :type in_range: 1D-array of int
    :param in_azimuth: azimuth of points
    :type in_azimuth: 1D-array of int
    :param in_v_long: longitude of points
    :type in_v_long: 1D-array of float
    :param in_v_lat: latitude of points
    :type in_v_lat: 1D-array of float
    :param in_hull_method: hull computation method, can be 2.0 or 2.1
    :type in_hull_method: float

    :return: the hull of the input set of points
    :rtype: OGRMultiPolygon
    """

    logger = logging.getLogger("my_hull")

    nb_pix_max = 5e4
    # Split lake sar image in sub_images
    list_sub_img_to_process = split_lake_image(in_range, in_azimuth, in_v_long, in_v_lat, nb_pix_max)
    if len(list_sub_img_to_process) > 1:
        logger.debug("Lake image with %d pixels is splitted in %d sub images of less than %d pixels" % (len(in_range),
                                                                                                       len(list_sub_img_to_process), nb_pix_max))

    multi_poly = ogr.Geometry(ogr.wkbMultiPolygon)

    for i, (sub_range, sub_azimuth, sub_v_long, sub_v_lat) in enumerate(list_sub_img_to_process):
        if i%100 == 0:
            logger.debug("%d over %d sub-image already processed" % ( i+1, len(list_sub_img_to_process)))
        sub_poly = get_polygon_from_binar_image(sub_range, sub_azimuth, sub_v_long, sub_v_lat, in_hull_method)
        if sub_poly.GetGeometryName() == "MULTIPOLYGON":
            for geom in sub_poly:
                multi_poly.AddGeometry(geom)
        else:
            multi_poly.AddGeometry(sub_poly)

    if not multi_poly.IsEmpty():
        poly = multi_poly.UnionCascaded()
    else:
        poly = multi_poly

    return poly

def get_polygon_from_binar_image( in_range, in_azimuth, in_v_long, in_v_lat, in_hull_method):
    """
    Compute the concave hull of a set of points using radar vectorisation

    :param in_range: range of points
    :type in_range: 1D-array of int
    :param in_azimuth: azimuth of points
    :type in_azimuth: 1D-array of int
    :param in_v_long: longitude of points
    :type in_v_long: 1D-array of float
    :param in_v_lat: latitude of points
    :type in_v_lat: 1D-array of float
    :param in_hull_method: hull computation method, can be 2.0 or 2.1
    :type in_hull_method: float

    :return: the hull of the input set of points
    :rtype: OGRMultiPolygon
    """
    # Get instance of service config file
    logger = logging.getLogger("my_hull")

    # 1 - Get relative range and azimuth (1st pixel = 1)
    lake_x = in_range - np.min(in_range) + 1
    lake_y = in_azimuth - np.min(in_azimuth) + 1

    # 2 - Compute contour using range and azimuth image coordinates
    lake_contour_int = compute_contour_from_range_azimuth(lake_x, lake_y)
    # 3 - Convert (azimuth, range) contour into polygon
    logger.debug("Inital polygon contains 1 external ring and %d holes rings " % (len(lake_contour_int) - 1))

    # 4 - Convert contours in lon lat
    # multi_ring_list contains every ring composing the polygon, the first ring contains exterior coordinates, 
    # all other rings are holes in the polygon
    multi_list_of_points = []
    for contour in lake_contour_int:  # Loop over contours
        if len(contour) < 3:
            logger.debug("Current contour contains only %d points ... skip" % (len(contour)))
            continue
        elif len(contour) > 8000:
            logger.warning("Current contour contains %d points ... can be time consuming" % len(contour))
        else:
            logger.debug("Building new ring with %d points " % len(contour))

        list_of_points = []

        for i, (y, x) in enumerate(contour):  # Look over azimuth and range indices

            # Retrieve lon/lat coordinates from range and azimtu coordinates
            # if current range and azimuth are found in input range and azimuth list

            if np.logical_and(lake_x == x, lake_y == y).any():
                point_idx = np.where(np.logical_and(lake_x == x, lake_y == y))[0][0]
                lon = in_v_long[point_idx]
                lat = in_v_lat[point_idx]

                new_point = (lon, lat)
                if in_hull_method == 2.0:
                    list_of_points.append(new_point)
                elif in_hull_method == 2.1:
                    # Add new point :
                    #     - if new_point not in list
                    #     - if list contains more than 3 points : check if new points create a crossing between segments
                    if new_point not in list_of_points:
                        if len(list_of_points) < 3:
                            list_of_points.append(new_point)
                        else:
                            list_of_points = add_new_point_to_list_of_point(new_point, list_of_points)

            else:
                logger.debug("Point of coordinates %d, %d not found -> Point removed" % (y, x))
        if not list_of_points:
            logger.debug("Ring contains 0 points => Discarded")
            continue

        if in_hull_method == 2.0:
            multi_list_of_points.append(list_of_points)
        else:
            # Add first point to the end of list
            ring = build_ring_from_list_of_points(list_of_points)

            # Check if ring does not intersect multi ring
            ring = build_ring_without_integrity_issues_multi_ring(ring, multi_list_of_points)

            if len(ring) > 3:
                # logger.debug("Adding new ring containing %d points to multi ring" % (len(ring[:-1])))
                multi_list_of_points.append(ring)

    # 5 - Build ogr geometry
    lake_poly = get_ogr_polygon_from_list_of_list_of_points(multi_list_of_points)

    if not lake_poly.IsValid():
        logger.warning("Polygon is invalid -> Polygon is downgraded into a valid geometry")
        lake_poly = lake_poly.Buffer(0)
    return lake_poly
    
def compute_contour_from_range_azimuth(in_x, in_y):
    """
    Compute contour of object with coordinates in_x and in_y

    :param in_x: x coordinates
    :type in_x: 1D-array of int
    :param in_y: y coordinates
    :type in_y: 1D-array of int

    :return: geometry object coordinates
    :rtype: list of list of tuple of int
    """
    # Get instance of service config file
    logger = logging.getLogger("my_hull")

    # 1 - Get image (1 pixel around lake)
    lake_img = compute_bin_mat(np.max(in_x) + 2, np.max(in_y) + 2, in_x, in_y, verbose=False)
    logger.debug("Processing image with size %d x %d, with %d water pixels" % (
    lake_img.shape[0], lake_img.shape[1], np.sum(lake_img)))

    # 2 - Filter isolated pixels
    # Make a erosion to check that a contour can be computed
    lake_erosion = binary_erosion(lake_img, structure=np.ones((2, 2))).astype(lake_img.dtype)

    # 3 - Compute boundaries (there might be more than one if there are islands in lake)
    if np.sum(lake_erosion) == 0:
        logger.warning("Current lake pixc no not contain a valid polygon")
        lake_contours = []
    else:
        lake_contours = find_contours(lake_img, 0.99999999999)

    # 4 - Round contour range and azimuth coordinates to units, since they are indices in input parameters
    # Removing duplicate points from lake contour
    lake_contour_int = []

    for contour in lake_contours:
        contour_points = []

        [contour_points.append(tuple(np.round(coords, 0))) for coords in contour if
         tuple(np.round(coords, 0)) not in contour_points]

        lake_contour_int.append(contour_points)

    return lake_contour_int

def compute_bin_mat(in_size_x, in_size_y, in_x, in_y, verbose=True):
    """
    Creates a 2D binary matrix from Y and Y 1D vectors
    i.e. for each ind: mat[in_x[ind],in_y[ind]] = 1

    :param in_size_x: number of pixels in X dimension
    :type in_size_x: int
    :param in_size_y: number of pixels in Y dimension
    :type in_size_y: int
    :param in_x: X indices of "1" pixels
    :type in_x: 1D vector of int
    :param in_y: Y indices of "1" pixels
    :type in_y: 1D vector of int
    :param verbose: print log debug flag
    :type verbose: bool

    :return: 2D matrix with "1" for each (in_x_i, in_y_i) and 0 elsewhere
    :rtype: 2D binary matrix of int 0/1
    """
    logger = logging.getLogger("my_tools")
    if verbose:
        logger.debug("- start -")

    # 0 - Deal with exceptions
    # 0.1 - Input vectors size must be the same
    if in_x.size != in_y.size:
        message = "compute_bin_mat(in_x, in_y) : in_x and in_y must be the same size ; currently : in_x = %d and in_y = %d" % (in_x.size, in_y.size)
        raise service_error.ProcessingError(message, logger)
    else:
        nb_pts = in_x.size
    if verbose:
        logger.debug("> Nb pixels to deal with = %d" % nb_pts)
    # 0.2 - max(X) < in_size_x
    if np.max(in_x) >= in_size_x:
        message = "compute_bin_mat(in_x, in_y) : elements of in_x must be less than in_size_x"
        raise service_error.ProcessingError(message, logger)
    # 0.3 - max(X) < in_size_x
    if np.max(in_y) >= in_size_y:
        message = "compute_bin_mat(in_x, in_y) : elements of in_y must be less than in_size_y"
        raise service_error.ProcessingError(message, logger)
    # 1 - Init output binary image
    out_bin_im = np.zeros((in_size_y, in_size_x))
    if verbose:
        logger.debug("> Binary matrix size = (X=%d , Y=%d)" % (in_size_x, in_size_y))

    # 2 - Put 1 for every pixels defined by the input vectors
    for ind in range(nb_pts):
        out_bin_im[in_y[ind], in_x[ind]] = 1

    return out_bin_im

def get_ogr_polygon_from_list_of_list_of_points(list_list_of_points):
    """
    return ogr geometry polygon from multi ring list

    :param list_list_of_points: list of rings
    :type list_list_of_points: list of list of tuple of 2 floats

    :return: lake polygon
    :rtype: OGRPolygon
    """
    logger = logging.getLogger("my_hull")

    # 1. build a list of ogr geometry "polygon" and "linear ring"
    list_of_ogr_poly = []
    list_of_ogr_ring = []
    for ring in list_list_of_points:
        lake_ring_geom = ogr.Geometry(ogr.wkbLinearRing)
        for (lon, lat) in ring:
            lake_ring_geom.AddPoint(lon, lat)
        lake_ring_geom.CloseRings()
        
        poly_tmp =  ogr.Geometry(ogr.wkbPolygon)
        poly_tmp.AddGeometry(lake_ring_geom)
        list_of_ogr_ring.append(lake_ring_geom)
        list_of_ogr_poly.append((poly_tmp, poly_tmp.Centroid()))
    # 2. build a list of outer / inner rings, with inner and outer relations
    list_of_outer_ring = []
    list_of_inner_ring = []
    for i, (poly_i, centroid_i) in enumerate(list_of_ogr_poly):
        is_i_inner = False
        for j, (poly_j, centroid_j) in enumerate(list_of_ogr_poly):
            if i == j:
                continue
            if poly_i.Within(poly_j):
                if j not in list_of_outer_ring:
                    list_of_outer_ring.append(j)
                list_of_inner_ring.append((j,i))
                is_i_inner = True
        if is_i_inner == False:
            list_of_outer_ring.append(i)

    logger.debug("Current geoemety contains %d inner rings" %len(list_of_outer_ring))
    # 3. Build out geometry multipolygon, with outer / inner rings in the right order.
    multi_poly_geom = ogr.Geometry(ogr.wkbMultiPolygon)
    for i in list_of_outer_ring:
        poly_geom = list_of_ogr_poly[i][0]
        for (out_ring, in_ring) in list_of_inner_ring:
            if i == out_ring:
                poly_geom.AddGeometry(list_of_ogr_ring[in_ring])
        multi_poly_geom.AddGeometry(poly_geom)

    return multi_poly_geom

def add_new_point_to_list_of_point(new_point, point_list):
    """
    Add new point to the list of point if it doesn't imply a crossing in the line formed by the point_list

    :param new_point: lon lat coordinates of new point
    :type new_point: tuple of floats
    :param point_list: list of point composing the contour
    :type point_list: list of tuple of floats


    :return: list of point with new point
    :rtype: list of tuple of 2 floats
    """
    logger = logging.getLogger("my_hull")

    # 1. Add new point to the list
    point_list.append(new_point)

    # 2. Check if last point implies a crossing in the line
    s1 = LineString([point_list[-2], point_list[-1]])
    s2 = LineString(point_list[:-2])
    inter = compute_segment_intersection(s1, s2)

    # 3. If an intersection in the line is found, the point close to the intersection is removed
    while inter:
        if inter:
            p = get_closest_point_from_list_of_point(inter, point_list[-4:])
            logger.debug("Removing point %f, %f" % (p[0], p[1]))
            point_list.remove(p)

        # list of point cannot be smaller than 3 points
        if len(point_list) < 4:
            break

        s1 = LineString([point_list[-2], point_list[-1]])
        s2 = LineString(point_list[:-2])

        inter = compute_segment_intersection(s1, s2)

    return point_list

def build_ring_without_integrity_issues_multi_ring(new_ring, multi_ring):
    """
    return ring without intersection with multi_ring

    :param new_ring: new ring
    :type new_ring: list of tuple of 2 floats
    :param multi_ring: list of rings
    :type multi_ring: list list of tuple of 2 floats


    :return: ring without intersection with multi_ring
    :rtype: list of tuple of 2 floats
    """

    logger = logging.getLogger("my_hull")

    s1 = LineString(new_ring)
    s2 = MultiLineString(multi_ring)

    inter_list = compute_segment_intersection(s1, s2)

    while inter_list and len(new_ring) > 3:

        p = get_closest_point_from_list_of_point(inter_list, new_ring)
        logger.debug("Removing point %f, %f" % (p[0], p[1]))
        new_ring.remove(p)

        if p == new_ring[-1]:
            new_ring.remove(p)
            new_ring.append(new_ring[0])

        s1 = LineString(new_ring)
        s2 = MultiLineString(multi_ring)

        inter_list = compute_segment_intersection(s1, s2)

    return new_ring

def get_closest_point_from_list_of_point(point, list_of_points):
    """
    Return element of list_of_points the closest of point

    :param point: lon lat coordinates of point
    :type point: tuple of 2 floats
    :param list_of_points: segment formed by the line without last two points
    :type list_of_points: list of tuple of 2 floats

    :return: intersection point coordinates
    :rtype: tuple of floats
    """

    if isinstance(point, list):
        id = distance.cdist([point], list_of_points).argmin(axis=1)[0]
    else:
        id = distance.cdist([point], list_of_points).argmin()
    return list_of_points[id]

def compute_segment_intersection(s1, s2):
    """
    Return the intersection point of s1 and s2

    :param s1: segment formed by the two last points of the line
    :param s1: segment formed by the two last points of the line
    :type s1: Shapely LineString
    :param s2: segment formed by the line without last two points
    :type s2: Shapely LineString

    :return: intersection point coordinates
    :rtype: tuple of floats
    """
    logger = logging.getLogger("my_hull")
    inter = s1.intersection(s2)

    if inter:
        # if the intersection is more than one point, return coordinates of the first intersection element
        if isinstance(inter, Point):
            inter = (inter.x, inter.y)
        elif isinstance(inter, MultiPoint):
            inter = (inter[0].x, inter[0].y)
        elif isinstance(inter, LineString):
            inter = (np.mean([inter.coords[0][0], inter.coords[1][0]]),
                     np.mean([inter.coords[0][1], inter.coords[1][1]]))
        elif isinstance(inter, MultiLineString):
            inter = (np.mean([inter[0].coords[0][0], inter[0].coords[1][0]]),
                     np.mean([inter[0].coords[0][1], inter[0].coords[1][1]]))
        elif isinstance(inter, GeometryCollection):
            for geom in inter:
                if isinstance(geom, Point):
                    inter_out = (geom.x, geom.y)
                elif isinstance(geom, MultiPoint):
                    inter_out = (geom[0].x, geom[0].y)
                elif isinstance(geom, LineString):
                    inter_out = (np.mean([geom.coords[0][0], geom.coords[1][0]]),
                             np.mean([geom.coords[0][1], geom.coords[1][1]]))
                elif isinstance(geom, MultiLineString):
                    inter_out = (np.mean([geom[0].coords[0][0], geom[0].coords[1][0]]),
                             np.mean([geom[0].coords[0][1], geom[0].coords[1][1]]))
                else:
                    logger.warning("Unknown intersection geometry type : %s " % (type(geom)))
                    inter_out = None

                inter = inter_out
                break
        else:
            logger.warning("Unknown intersection geometry type : %s " % (type(inter)))
            inter = None

    return inter
    
def build_ring_from_list_of_points(list_of_points):
    """
    Add first point of the list and the and of the list. Take in account crossings problems

    :param list_of_points: list of point composing the contour
    :type list_of_points: list of tuple of 2 floats

    :return:  list of point composing the contour, with same first and last point
    :rtype: list of tuple of 2 floats
    """
    logger = logging.getLogger("my_hull")

    if len(list_of_points) > 3:

        s1 = LineString([list_of_points[-1], list_of_points[0]])
        s2 = LineString(list_of_points[1:-1])

        inter_list = compute_segment_intersection(s1, s2)

        while inter_list and len(list_of_points) > 3:

            p = get_closest_point_from_list_of_point(inter_list, list_of_points[:2] + list_of_points[-2:])
            logger.debug("Removing point %f, %f" % (p[0], p[1]))

            list_of_points.remove(p)

            if len(list_of_points) > 3:
                s1 = LineString([list_of_points[-1], list_of_points[0]])
                s2 = LineString(list_of_points[1:-1])
                inter_list = compute_segment_intersection(s1, s2)

    list_of_points.append(list_of_points[0])

    return list_of_points

def from_ogr_to_shapely(ogr_geom):
    
    # Creating a copy of the input OGR geometry. This is done in order to 
    # ensure that when we drop the M-values, we are only doing so in a 
    # local copy of the geometry, not in the actual original geometry. 
    
    
    ogr_geom_copy = ogr_geom.Clone()
    ogr_geom_copy.FlattenTo2D() 

    # Dropping the M-values
    ogr_geom_copy.SetMeasured(False)

    # Generating a new shapely geometry
    shapely_geom = shapely.wkt.loads(ogr_geom_copy.ExportToWkt())

    return shapely_geom
              
def compute_alpha_shape(water: gpd.GeoDataFrame, nb_pix_range: int, method: str = "cgal") -> gpd.GeoDataFrame:
    '''
    Copute alpha_shape

    :param water: Water Dataframe

    :return points: Extracted water points 
    '''
    data = np.array(water[['utm_x','utm_y']].values)
    in_v_long = water['longitude'].values
    in_v_lat = water['latitude'].values
    in_range = water['range_index'].values
    in_azimuth = water['azimuth_index'].values
        
    polygons = []
    
    if method == 1.2:
        polygons = alpha_shape_with_cascaded_union(data, alpha) 
        
    if method == 2.0 or method == 2.1:
        coords = np.zeros((in_v_long.size, 2))
        coords[:, 0], coords[:, 1] = water['utm_x'].values, water['utm_y'].values
        
        in_hull_method = 2.0
        polygons = from_ogr_to_shapely(get_concave_hull_from_radar_vectorisation(in_range, in_azimuth, coords[:, 0], coords[:, 1], method))
    
    elif method == 1.0:
        in_nb_pix_range = nb_pix_range
        
        coords = np.zeros((in_v_long.size, 2))
        coords[:, 0], coords[:, 1] = water['utm_x'].values, water['utm_y'].values

        # Separate swath in different zone : near and far range. This step is useful to adjust alpha regarding range sampling.
        nb_range_alpha_sampling = 20
        
        
        for i in range(nb_range_alpha_sampling):
            rg_bound_min = i*(in_nb_pix_range / nb_range_alpha_sampling)-10
            rg_bound_max = (i+1)*(in_nb_pix_range / nb_range_alpha_sampling)+10
            
            rg_idx = np.where(np.logical_and(in_range < rg_bound_max, in_range > rg_bound_min))
            alpha_near_rg, dist_near_rg = evaluate_alpha_from_x_pixc_geolocation(coords[:, 0][rg_idx], in_range[rg_idx], in_azimuth[rg_idx])
            if dist_near_rg :
                concave_hull_rg = alpha_shape_with_cgal(coords[rg_idx], alpha_near_rg)
            else :
                concave_hull_rg = MultiPolygon()
                
            if i==0:
                concave_hull = concave_hull_rg
            else:
                concave_hull = unary_union([concave_hull, concave_hull_rg])
                
       
        if concave_hull.type == MultiPolygon:
            polygons = list(concave_hull)
        else:
            polygons = concave_hull
    
    else: 
        raise Exception("Method not found")

    return polygons 
   

           
def evaluate_alpha_from_x_pixc_geolocation(in_x, in_range, in_azimuth):
    """
    Find the best value of alpha regarding the distance between neighbour pixels.

    :param in_x: X projected coordinates of pixels
    :type in_x: 1D-array of float
    :param in_range: range of pixels
    :type in_range: 1D-array of int
    :param in_azimuth: azimuth of pixels
    :type in_azimuth: 1D-array of int

    :return alpha value from 250 to 5000 and distance between 2 pixels variation max (mean + 2*sdt)
    :rtype: tuple of (int, float)
    """

    if len(in_x) == 0:
        alpha = None
        x_2sigma = None
    else:
        # dist_list contains distances between two neighbour pixel along the X dimension
        # Only 20 measurements are needed
        dist_list = []
        cpt = 0
        while(len(dist_list) < 20):
            cpt += 1
            if cpt > 30:
                break
            idx = randint(0, len(in_x)-1)
            dist = get_dist_if_neighbour(idx, in_x, in_range, in_azimuth)
            if dist:
                dist_list.append(dist[0])
    
        x_mean = np.mean(dist_list)
        x_std = np.std(dist_list)
    
        x_2sigma = x_mean + 2 * x_std
    
        if not dist_list:
            alpha = 5000
        elif x_2sigma < 10 :
            alpha = 250
        elif x_2sigma > 10 and x_2sigma < 30 :
            alpha = 500
        elif x_2sigma > 30 and x_2sigma < 50 :
            alpha = 1000
        elif x_2sigma > 50 and x_2sigma < 100:
            alpha = 2000
        else :
            alpha = 5000
        alpha = 5000
        
        ## TODO : Improve the alpha shape tunning / divide "big" polygon into multiples to have a better alpha shape well suited 
        # ~ alpha = alpha*2.5
    return alpha, x_2sigma

def get_dist_if_neighbour(idx, in_x, in_range, in_azimuth):
    """
    For a given pixel of indice idx :
        - find if the pixel have a direct neighbour in range
        - compute the distance in X dimension between the pixel and its neighbour if neighbouir exists
        - return None if neighbour do not exists


    :param idx: indice of pixel
    :type idx: int
    :param in_x: X projected coordinates of pixels
    :type in_x: 1D-array of float
    :param in_range: range of pixels
    :type in_range: 1D-array of int
    :param in_azimuth: azimuth of pixels
    :type in_azimuth: 1D-array of int

    :return distance between pixel of indice idx and its neighbour
    :rtype: float
    """

    rg = in_range[idx]
    az = in_azimuth[idx]
    # Find pixel neighbour
    x_neighbour = np.where(np.logical_and(in_azimuth == az, in_range == rg+1))
    dist = None
    if x_neighbour:
        # compute distance if neighbour exist
        dist = np.abs(in_x[idx] - in_x[x_neighbour])
    else:
        dist = None
    return dist


def get_borders(data: gpd.GeoDataFrame) -> Polygon:
    '''
    Get borders

    :param data: Point cloud information in DataFrame
    :returns: Square polygon which contains point cloud
    '''
    # Image borders
    y_min = data['latitude'].values.min()
    y_max = data['latitude'].values.max()
    x_min = data['longitude'].values.min()
    x_max = data['longitude'].values.max()
    return Polygon([(x_min,y_min),
                    (x_max,y_min),
                    (x_max,y_max),
                    (x_min,y_max),
                    (x_min,y_min)]).exterior


def get_utm_borders(data: gpd.GeoDataFrame) -> Polygon:
    '''
    Get borders in utm coodinates

    :param data: Point cloud information in DataFrame
    :returns: Square polygon which contains point cloud
    '''
    # Image borders
    y_min = data['utm_y'].values.min()
    y_max = data['utm_y'].values.max()
    x_min = data['utm_x'].values.min()
    x_max = data['utm_x'].values.max()
    return Polygon([(x_min,y_min),
                    (x_max,y_min),
                    (x_max,y_max),
                    (x_min,y_max),
                    (x_min,y_min)]).exterior




def extract_boundary_points(data: gpd.GeoDataFrame, boundaries: List[LinearRing]):
    '''
    Extract points from boundaries

    :param data: Point cloud information in DataFrame
    :param boundaries: Polygon boundaries
    
    :return: Boundary points list
    :rtype: GeoPandas Dataframe

    '''
    # Create index
    data = data.reset_index()
    data = data.set_index(['utm_x','utm_y'])

    # ~ get_point_from_coord = lambda coord: tuple(data.loc[coord,['longitude','latitude','elevation','range_index','azimuth_index','time']].values)  
    
    def get_point_from_coord(coord):
        try:
            return tuple(data.loc[coord,['longitude','latitude','elevation','range_index','azimuth_index','time']].values)  
        except:
            print("intersection not found")
            return None
        return tuple(data.loc[coord,['longitude','latitude','elevation','range_index','azimuth_index','time']].values)  
        
            
    points = [get_point_from_coord(coord) for boundary in boundaries for coord in boundary.coords ]
    Not_none_values = filter(None.__ne__, points)
    points = list(Not_none_values)
    
    # ~ get_mean_height = lambda coord: np.array(data.loc[coord,['wse']]).flatten()
    # ~ mean_wse = []
    # ~ for boundary in boundaries:
        # ~ toto = [get_mean_height(coord) for coord in boundary.coords]
        # ~ mean_wse.append(np.mean(toto))
    
    # Convert to dataframe
    df = pd.DataFrame(points,columns=['longitude','latitude','elevation','range','azimuth','time'])
        
    lambdafunc = lambda row: pd.Series([*utm.from_latlon(row['latitude'],
                                                   row['longitude'])[0:2],
                                  row['elevation'],
                                      ]) 
    df[['x','y','z']] = df.apply(lambdafunc,axis=1)
    df = df.astype(dtype={'longitude':'double','latitude':'double', 'elevation':'double',
                          'range':'int','azimuth':'int',
                         'x':'double','y':'double', 'z':'double', 'time':'int'})
    

    geom = [Point(x,y) for x, y in zip(df['longitude'], df['latitude'])]

    # ~ return gpd.GeoDataFrame(df, geometry=geom), mean_wse
    return gpd.GeoDataFrame(df, geometry=geom)

def compute_mean_wse(data: gpd.GeoDataFrame, polygons: List[Polygon]):
    
    def get_mean_height(coord):
        try:
            return tuple(data.loc[coord,['elevation']])
        except:
            print("intersection not found")
            return None

    mean_wse = []
    data = data.reset_index()
    data = data.set_index(['utm_x','utm_y'])
    for poly in polygons:
        toto = [get_mean_height(coord) for coord in poly.exterior.coords]
        Not_none_values = filter(None.__ne__, toto)
        toto = list(Not_none_values)
        mean_wse.append(np.mean(toto))
    return mean_wse
    

def remove_borders(data: gpd.GeoDataFrame, range_max:int, azimuth_max:int, min_range_indices_to_remove:np.array):
    '''
    Remove borders points

    :param data: Boundary points list
    :param range_max: Maximum range value for a tile
    :param azimuth_max: Maximum azimuth value for a tile
    :param min_range_indices_to_remove: minimum range_indice to remove for each azimtuh line
    
    :return: Boundary points list filtered
    :rtype: GeoPandas Dataframe
    '''
    # ~ data = data.loc[data.range != 0]
    data = data.loc[data.azimuth != 0]
    data = data.loc[data.range != range_max]
    data = data.loc[data.azimuth != azimuth_max]
    for i in range(azimuth_max):
        data = data.loc[np.logical_or((data['azimuth'] != i),(data['range'] != min_range_indices_to_remove[i]))]
        
    return data

def remove_near_range_pixels(water, azimuth_max, cross_track_min = 5000):
    '''
    Filter points whose cross-track distance in below a threshold

    :param water: GeoPandas Dataframe
    :param azimuth_max: Maximum azimuth value for a tile
    :param cross_track_min: minimum cross-track value to keep in data
    
    :return: GeoPandas Dataframe filtered
    :rtype: GeoPandas Dataframe
    :return: Numpy 1D Array with range indices corresponding to the border of the tile, for each azimuth position
    :rtype: numpy array    
    '''
    
    water_fil = water.loc[(np.abs(water['cross_track']) > cross_track_min)]
    water_removed = water.loc[(np.abs(water['cross_track']) <= cross_track_min)]
    
    min_range_indices_to_remove = np.zeros([azimuth_max])
    
    for i in range(azimuth_max):
        try:
            min_after_filtering = np.nanmin(water_fil.loc[water_fil['azimuth_index']==i]['range_index'])
            max_filtered_area = np.nanmax(water_removed.loc[water_removed['azimuth_index']==i]['range_index'])
            if max_filtered_area == min_after_filtering -1 :
                min_range_indices_to_remove[i] = min_after_filtering
        except:
            pass
    return water_fil, min_range_indices_to_remove

# End
