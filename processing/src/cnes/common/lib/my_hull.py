# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2019/05/17:version initiale.
# VERSION:2.0.0:DM:#91:2020/07/03:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: my_hull.py
    :synopsis: deal with hull computation (=lake boundary computation)
        Created on 2018/08/27

.. moduleauthor:: Claire POTTIER (CNES DSO/SI/TR) - Cécile CAZALS (CS) - Emmanuelle Sarrazin (CNES)

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

    Library uses alpha_shape implementation based CGAL See https://www.cgal.org/

"""

import logging
import math
import numpy as np
from osgeo import ogr
import pandas as pd
from random import randint
from scipy.spatial import distance, Delaunay
import shapely.geometry as geometry
from shapely.geometry import Point, LineString, MultiPoint, MultiLineString, Polygon, MultiPolygon, LinearRing, GeometryCollection
from shapely.ops import cascaded_union
from shapely.ops import unary_union
from skimage.measure import find_contours

import cnes.common.service_config_file as service_config_file
import cnes.common.service_error as service_error

import cnes.common.lib.my_tools as my_tools


def compute_lake_boundaries(in_v_long, in_v_lat, in_range, in_azimuth, in_nb_pix_range):
    """
    Compute the hull of a set of points determined by their coordinates given in input parameters

    :param in_v_long: longitudes of points for which to compute the hull
    :type in_v_long: 1D-array of float
    :param in_v_lat: latitudes of points for which to compute the hull
    :type in_v_lat: 1D-array of float
    :param in_range: range of points for which to compute the hull
    :type in_range: 1D-array of int
    :param in_azimuth: azimuth of points for which to compute the hull
    :type in_azimuth: 1D-array of int
    :param in_nb_pix_range: maximal number of pixel in range
    :type in_nb_pix_range: int

    :return the hull of the input set of points
    :rtype: OGRMultiPolygon
    """
    # Get instance of service config file
    cfg = service_config_file.get_instance()
    logger = logging.getLogger("my_hull")
    logger.debug("Computing lake boundaries over %d pixels" %(len(in_v_long)))

    hull_method = cfg.getfloat("CONFIG_PARAMS", "HULL_METHOD")
    if len(in_v_long) < 4:
        logger.debug("Hull computation method : Convex hull because less than 4 pixels")
        retour = get_convex_hull(in_v_long, in_v_lat)

    elif hull_method == 0:  # 0 - CONVEX HULL
        logger.debug("Hull computation method : Convex hull")
        retour = get_convex_hull(in_v_long, in_v_lat)
        
    elif hull_method == 1.0:  # 1.0 : CONCAV HULL - Delaunay triangulation with CGAL
        logger.debug("Hull computation method : Concav hull computed in ground geometry, based on Delaunay triangulation - using CGAL library")
        retour = get_concave_hull_from_cgal_triangulation(in_v_long, in_v_lat, in_range, in_azimuth, in_nb_pix_range)
        
    elif hull_method == 1.1:  # 1.1 - CONCAV HULL - Delaunay triangulation 
        logger.debug("Hull computation method : Delauney triangulation")
        retour = get_concave_hull_from_basic_triangulation(in_v_long, in_v_lat, in_range, in_nb_pix_range)

    elif hull_method == 2:  # 2 - CONCAV HULL - Radar vectorisation method
        logger.debug("Hull computation method : radar vectorization (default)")
        retour = get_concave_hull_from_radar_vectorisation(in_range, in_azimuth, in_v_long, in_v_lat)

    else:
        message = "Concave hull computation method not understood"
        raise service_error.ProcessingError(message, logger)

    if retour.GetGeometryName() not in ["POLYGON", "MULTIPOLYGON"]:
        retour = ogr.Geometry(ogr.wkbPolygon)
    return retour
	

#######################################


def remove_holes_triangles(polygon) :
    """
    Remove holes

    :param polygon: polygon
    :type polygon: geometry.polygon or geometry.MultiPolygon

    :return polygon without triangle
    :rtype: polygon
    """
    polygon_without_triangle = polygon

    if type(polygon) == geometry.Polygon:
        polygon_without_triangle = Polygon(polygon.exterior,
                                        [hole for hole in polygon.interiors if len(hole.coords) > 4])

    if type(polygon) == geometry.MultiPolygon:
        triangle_union_smooth = []
        for poly in polygon:
            triangle_union_smooth.append(Polygon(poly.exterior,
                                                 [hole for hole in poly.interiors if len(hole.coords) > 4]))
            polygon_without_triangle = MultiPolygon(triangle_union_smooth)

    return polygon_without_triangle


#######################################


def get_convex_hull(in_v_long, in_v_lat):
    """
    Compute the convex hull of a set of points

    :param in_v_long: longitudes of points for which to compute the hull
    :type in_v_long: 1D-array of float
    :param in_v_lat: latitudes of points for which to compute the hull
    :type in_v_lat: 1D-array of float

    :return the hull of the input set of points
    :rtype: OGRPolygon
    """

    # 1.1 - Group pixels in a OGRGeometryCollection
    pixc_pts = ogr.Geometry(ogr.wkbGeometryCollection)
    for indp in np.arange(in_v_long.size):
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(in_v_long[indp], in_v_lat[indp])
        pixc_pts.AddGeometry(point)

    # 1.2 - Compute the hull of these points
    return pixc_pts.ConvexHull()


#######################################


def get_concave_hull_from_basic_triangulation(in_v_long, in_v_lat, in_range, in_nb_pix_range):
    """
    Compute the concave hull of a set of points using classifcal Delauney triangulation.

    :param in_v_long: longitudes of points for which to compute the hull
    :type in_v_long: 1D-array of float
    :param in_v_lat: latitudes of points for which to compute the hull
    :type in_v_lat: 1D-array of float
    :param in_range: range of points for which to compute the hull
    :type in_range: 1D-array of int
    :param in_nb_pix_range: maximal number of pixel in range
    :type in_nb_pix_range: int

    :return the hull of the input set of points
    :rtype: OGRMultiPolygon
    """
    # Get instance of service config file
    cfg = service_config_file.get_instance()

    logger = logging.getLogger("my_hull")

    nb_pix = in_v_long.size
    # ============================ HULL COMPUTATION FOR BIG LAKES ===========================
    # ============================  TO CHANGE AS SOON AS POSSIBLE ===========================
    # if the lake contains more than 500 000 pixels, the number of pixel is restricted to 500 000. Pixels are randomly selected.
    nb_pix_max_delauney = cfg.getint("CONFIG_PARAMS", "NB_PIX_MAX_DELAUNEY")
    if nb_pix > nb_pix_max_delauney:
        id_selected = np.random.randint(nb_pix, size=int(nb_pix_max_delauney))
        logger.warning(
            "The number of pixel of the lake is reduced to %d" % (nb_pix_max_delauney))
        in_v_long = in_v_long[id_selected]
        in_v_lat = in_v_lat[id_selected]
        in_range = in_range[id_selected]
        alpha_ratio = nb_pix_max_delauney/nb_pix
        nb_pix = in_v_long.size
        logger.debug("alpha_ratio : " + str(alpha_ratio))
    else:
        alpha_ratio = 1.
        
    # =======================================================================================

    # 1.1 - Gather coordinates in a 2D-array
    coords = np.zeros((nb_pix, 2))

    # 1.2 - Transform geographical coordinates into utm coordinates
    coords[:, 0], coords[:, 1], utm_epsg_code = my_tools.get_utm_coords_from_lonlat(in_v_long, in_v_lat)

    # 1.3 - Compute alpha shape
    hull_method = cfg.getfloat("CONFIG_PARAMS", "HULL_METHOD")
    if hull_method == 1.1:  # Without alpha parameter
        concave_hull_utm = get_concav_hull_bis(coords)
    else:  # With alpha parameter
        alpha = alpha_ratio*( 0.03 + 0.01 * in_range / in_nb_pix_range)  # alpha parameter ranges from 0.03 to 0.04 following the range index
        concave_hull_utm = alpha_shape(coords, alpha)

    # 1.4 - Transform concave hull polygon into geographical coordinates
    concave_hull = my_tools.get_lon_lat_polygon_from_utm(concave_hull_utm, utm_epsg_code)

    # 1.5 - Convert Shapely geometry to OGRPolygon or OGRMultiPolygon
    return ogr.CreateGeometryFromWkb(concave_hull.wkb)


def alpha_shape(in_coords, in_alpha):
    """
    Compute the alpha shape (concave hull) of a set of points.
    
    :param in_coords: set of points coordinates
    :type in_coords: 2D-array of size (nb_pixels, 2=lon/lat)
    :param in_alpha: alpha value to influence the gooeyness of the border. Smaller numbers don't fall inward as much as larger numbers. Too large, and you lose everything!
    :type in_alpha: 1D-array of int
    
    :return: Shapely.MultiPolygons which is the hull of the input set of points
    
    Retrieved from http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
    """
    
    # Number of points
    nb_pts = in_coords.shape[0]

    # 0 - Particular case of nb points <= 3: return convex hull
    if nb_pts <= 3:
        # 0.1 - Init Shapely set of pointsn and Aggregate coordinates to the set of points
        points = [geometry.point.Point(in_coords[indp, 0], in_coords[indp, 1], 0) for indp in np.arange(nb_pts)]
        # 0.2 - Return convex hull
        retour = geometry.MultiPoint(list(points)).convex_hull
    else:
        # 1 - Compute Delaunay triangulation
        tri = Delaunay(in_coords)  
        # tri = specific object with attributes (cf. https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.Delaunay.html)

        # 2 - Select triangles following it's shape
        # 2.1 - Compute CircumRatio for each triangle
        circum_r = [get_circum_ratio(in_coords[ia], in_coords[ib], in_coords[ic]) for ia, ib, ic in tri.vertices]

        # 2.2 - Compute mean alpha parameter for each triangle
        mean_alpha = [1.0 / np.mean([in_alpha[ia], in_alpha[ib], in_alpha[ic]]) for ia, ib, ic in tri.vertices]

        # 2.3 - Select triangles
        triangles_selected = np.where(np.array(circum_r) < mean_alpha)

        # 2.4 - Compute a list of shapely polygon correpsonding to selected triangles
        list_triangle = [geometry.Polygon([(in_coords[ia][0], in_coords[ia][1]), (in_coords[ib][0], in_coords[ib][1]),
                                           (in_coords[ic][0], in_coords[ic][1])]) for ia, ib, ic in tri.vertices[triangles_selected]]

        # 3 - Union of selected triangles
        triangle_union = cascaded_union(list_triangle)

        retour = triangle_union
    return retour


def get_circum_ratio(pa, pb, pc):
    """
    Compute the circumscribing circle radius of a triangle. The triangle is given by its corner coordinates pa, pb, pc.
    
    :param pa: geographical coordinates of pa
    :type pa: tuple of float (lon, lat)
    :param pb: geographical coordinates of pb
    :type pb: tuple of float (lon, lat)
    :param pc: geographical coordinates of pc
    :type pc: tuple of float (lon, lat)
    
    :return: circumscribing circle radius of the triangle
    :rtype: float
    """
    
    # Lengths of sides of triangle
    a = math.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
    b = math.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
    c = math.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)

    # Semiperimeter of triangle
    s = (a + b + c) / 2.0

    # Area of triangle by Heron's formula
    try:
        area = math.sqrt(s * (s - a) * (s - b) * (s - c))
    except ValueError:
        circum_r = 1
    else:
        # Circumscribing circle radius
        if area != 0:
            circum_r = a * b * c / (4.0 * area)
        else:
            circum_r = 0

    return circum_r


#######################################


def get_concav_hull_bis(in_coords):
    """
    Compute the Delaunay triangulation and return the concave hull of cloud point in_coords

    :param in_coords: set of points coordinates
    :type in_coords: 2D-array of float ; size = (nb_pixels, 2=lon/lat)

    :return the hull of the input set of points
    :rtype: Shapely.MultiPolygons

    Partially retrieved from http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
    """

    # Number of points
    nb_pts = in_coords.shape[0]

    # 0 - Particular case of nb points <= 3: return convex hull
    if nb_pts <= 3:
        # 0.1 - Init Shapely set of points and aggregate its coordinates
        points = [geometry.point.Point(in_coords[indp, 0], in_coords[indp, 1], 0) for indp in np.arange(nb_pts)]
        # 0.2 - Return convex hull
        retour = geometry.MultiPoint(list(points)).convex_hull
    else:
        # 1 - Compute Delaunay triangulation
        tri = Delaunay(in_coords)
        # tri = specific object with attributes (cf. https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.Delaunay.html)

        # 2 - Select triangles following their shape

        # 2.1 - Get maximal length of each triangle
        maxseg = np.array([get_max_segment(in_coords[ia], in_coords[ib], in_coords[ic]) for ia, ib, ic in tri.vertices])

        # 2.2 - Get median length of maximal segments
        maxseg_median = np.median(maxseg)

        # 2.3 - Select triangles
        triangles_selected_idx = np.where(maxseg < maxseg_median * 3)

        # 2.4 - Compute a list of shapely polygon correpsonding to selected triangles
        list_triangle = [geometry.Polygon([(in_coords[ia][0], in_coords[ia][1]), (in_coords[ib][0], in_coords[ib][1]),
                                           (in_coords[ic][0], in_coords[ic][1])]) for ia, ib, ic in
                         tri.vertices[triangles_selected_idx]]

        # 3 - Union of selected triangles
        triangle_union = cascaded_union(list_triangle)

        # 4 - Deletion of some artifacts related to Delaunay triangulation: holes may appears in some cases
        if triangle_union.type == "Polygon":

            # Get holes in concerned polygon
            holes_in_poly = [list(poly.coords) for poly in list(triangle_union.interiors)]

            for poly_holes in list(triangle_union.interiors):
                # If the hole contains only 3 points (triangle) remove the hole
                if len(poly_holes.coords) == 4:
                    triangle_union = geometry.Polygon(triangle_union.exterior.coords,
                                                      holes_in_poly.remove(list(poly_holes.coords)))

        elif triangle_union.type == "MultiPolygon":

            multipolygon = []

            for polygon in triangle_union:

                # Get holes in concerned polygon
                holes_in_poly = [list(poly.coords) for poly in list(polygon.interiors)]

                for poly_holes in list(polygon.interiors):
                    # If the hole contains only 3 points (triangle) remove the hole
                    if len(poly_holes.coords) == 4:
                        multipolygon.append(
                            geometry.Polygon(polygon.exterior.coords, holes_in_poly.remove(list(poly_holes.coords))))

            triangle_union = cascaded_union(multipolygon)

        retour = triangle_union

    return retour


def get_max_segment(pa, pb, pc):
    """
    Compute the length of the longuest triangle edge. The triangle is given by its corner coordinates pa, pb, pc.

    :param pa: geographical coordinates of pa
    :type pa: tuple of float (lon, lat)
    :param pb: geographical coordinates of pb
    :type pb: tuple of float (lon, lat)
    :param pc: geographical coordinates of pc
    :type pc: tuple of float (lon, lat)

    :return: length of the longuest triangle edge
    :rtype: float
    """

    # Lengths of triangle edges
    a = math.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
    b = math.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
    c = math.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)

    # Get and return length of the longuest triangle edge
    return np.max([a, b, c])


#######################################

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

    retour = [(in_range, in_azimuth, in_v_long, in_v_lat)]
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

def get_concave_hull_from_radar_vectorisation(in_range, in_azimuth, in_v_long, in_v_lat):
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

    :return: the hull of the input set of points
    :rtype: OGRMultiPolygon
    """

    logger = logging.getLogger("my_hull")

    nb_pix_max = 5e4
    # Split lake sar image in sub_images
    list_sub_img_to_process = split_lake_image(in_range, in_azimuth, in_v_long, in_v_lat, nb_pix_max)
    if len(list_sub_img_to_process) > 1:
        logger.info("Lake image with %d pixels is splitted in %d sub images of less than %d pixels" % (len(in_range), 
                    len(list_sub_img_to_process), nb_pix_max))

    multi_poly = ogr.Geometry(ogr.wkbMultiPolygon)

    for i, (sub_range, sub_azimuth, sub_v_long, sub_v_lat) in enumerate(list_sub_img_to_process) :
        if i%100 == 0 :
            logger.info("%d over %d sub-image already processed" % ( i+1, len(list_sub_img_to_process)))
        sub_poly = get_polygon_from_binar_image(sub_range, sub_azimuth, sub_v_long, sub_v_lat)
        if sub_poly.GetGeometryName() == "MULTIPOLYGON":
            for geom in sub_poly:
                multi_poly.AddGeometry(geom)
        else :
            multi_poly.AddGeometry(sub_poly)

    poly = multi_poly.UnionCascaded()

    return poly

def get_polygon_from_binar_image( in_range, in_azimuth, in_v_long, in_v_lat):
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

    :return: the hull of the input set of points
    :rtype: OGRMultiPolygon
    """
    # Get instance of service config file
    cfg = service_config_file.get_instance()
    # get param nb_pix_max_contour
    nb_pix_max_contour = cfg.getint("CONFIG_PARAMS", "NB_PIX_MAX_CONTOUR")

    logger = logging.getLogger("my_hull")

    # 1 - Get relative range and azimuth (1st pixel = 1)
    lake_x = in_range - np.min(in_range) + 1
    lake_y = in_azimuth - np.min(in_azimuth) + 1

    # 2 - Get image (1 pixel around lake)
    lake_img = my_tools.compute_bin_mat(np.max(lake_x) + 2, np.max(lake_y) + 2, lake_x, lake_y, verbose = False)

    logger.info("Processing image with size %d x %d, with %d water pixels" %(lake_img.shape[0], lake_img.shape[1], np.sum(lake_img)))


    # 3 - Compute boundaries (there might be more than one if there are islands in lake)
    lake_contours = find_contours(lake_img, 0.99999999999)

    # 4 - Round contour range and azimuth coordinates to units, since they are indices in input parameters
    # Removing duplicate points from lake contour
    lake_contour_int = []

    for contour in lake_contours:
        contour_points = []

        [contour_points.append(tuple(np.round(coords, 0))) for coords in contour if
         tuple(np.round(coords, 0)) not in contour_points]

        lake_contour_int.append(contour_points)

    # 5 - Convert (azimuth, range) contour into polygon
    logger.debug("Inital polygon contains 1 external ring and %d holes rings " % (len(lake_contour_int) - 1))

    # multi_ring_list contains every ring composing the polygon, the first ring contains exterior coordinates, 
    # all other rings are holes in the polygon
    multi_ring_list = []

    for contour in lake_contour_int:  # Loop over contours
        logger.debug("Building new ring with %d points " % (len(contour)))
        # ============================ HULL COMPUTATION FOR BIG LAKES ===========================
        # ============================  TO CHANGE AS SOON AS POSSIBLE ===========================
        # if the lake contour contains more than 8 000 pixels, the number of pixel is restricted to 8 000..
        if len(contour) > nb_pix_max_contour:
            logger.warning("Number of contour points is reduced to %d points" % (nb_pix_max_contour))
            id_selected_points = np.round(np.linspace(0, len(contour) - 1, nb_pix_max_contour), 0).astype('int')
            contour = [contour[idx] for idx in id_selected_points]
        # =======================================================================================

        list_of_points = []

        for (y, x) in contour:  # Look over azimuth and range indices

            # Retrieve lon/lat coordinates from range and azimtu coordinates
            # if current range and azimuth are found in input range and azimuth list
            if np.where(np.logical_and(lake_x == x, lake_y == y))[0].any():
                point_idx = np.where(np.logical_and(lake_x == x, lake_y == y))[0][0]
                lon = in_v_long[point_idx]
                lat = in_v_lat[point_idx]

                new_point = (lon, lat)

                # Add new point :
                #     - if new_point not in list
                #     - if list contains more than 3 points : check if new points create a crossing between segments
                if new_point not in list_of_points:
                    if len(list_of_points) < 3:
                        list_of_points.append(new_point)
                    else:
                        list_of_points = add_new_point_to_list_of_point(new_point, list_of_points)
                # code commenté : list_of_points = addNewPoint(new_point, list_of_points)
            else:
                logger.debug("Point of coordinates %d, %d not found -> Point removed" % (y, x))

        if not list_of_points:
            logger.debug("Ring contains 0 points => Discarded")
            continue

        # Add first point to the end of list
        ring = build_ring_from_list_of_points(list_of_points)

        # Check if ring does not intersect multi ring
        ring = build_ring_without_integrity_issues_multi_ring(ring, multi_ring_list)

        if len(ring) > 3:
            # logger.debug("Adding new ring containing %d points to multi ring" % (len(ring[:-1])))
            multi_ring_list.append(ring)

        else:
            logger.debug("Ring contains less than 2 points => Discarded")

    lake_poly = get_ogr_polygon_from_ring_list(multi_ring_list)

    if not lake_poly.IsValid():
        logger.debug("Polygon is invalid -> Polygon is downgraded into a valid geometry")
        lake_poly = lake_poly.Buffer(0)

    return lake_poly


def compute_segment_intersection(s1, s2):
    """
    Return the intersection point of s1 and s2

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
            for geom in inter :
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
                else :
                    logger.warning("Unknown intersection geometry type : %s " % (type(geom)))
                    inter_out = None

                inter = inter_out
                break
        else:
            logger.warning("Unknown intersection geometry type : %s " % (type(inter)))
            inter = None

    return inter


def get_closest_point_from_list_of_ring(point, list_of_ring):
    """
    Return element of list_of_ring the closest of point

    :param point: lon lat coordinates of point
    :type point: tuple of 2 floats
    :param list_of_ring: TODO
    :type list_of_ring: TODO

    :return: intersection point coordinates
    :rtype: tuple of floats
    """
    logger = logging.getLogger("my_hull")
    min_dist = 100000
    for i, ring in enumerate(list_of_ring):
        if np.min(distance.cdist(point, ring)) < min_dist and len(ring) > 4:
            min_dist = np.min(distance.cdist(point, ring))
            min_point_id = distance.cdist(point, ring).argmin(axis=1)[0]
            min_ring_id = i

    min_point_coord = list_of_ring[min_ring_id][min_point_id]
    logger.debug("Removing point %f, %f" % (min_point_coord[0], min_point_coord[1]))
    list_of_ring[min_ring_id].remove(min_point_coord)

    return list_of_ring


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


def build_ring_from_list_of_points(list_of_points):
    """
    Add first point of the list and the and of the list. Take in account crossings problems

    :param list_of_points: list of point composing the contour
    :type list_of_points: list of tuple of 2 floats

    :return:  list of point composing the contour, with same first and last point
    :rtype: list of tuple of 2 floats
    """
    logger = logging.getLogger("my_hull")

    if len(list_of_points) > 3 :

        s1 = LineString([list_of_points[-1], list_of_points[0]])
        s2 = LineString(list_of_points[1:-1])

        inter_list = compute_segment_intersection(s1, s2)

        while inter_list and len(list_of_points) > 3 :

            p = get_closest_point_from_list_of_point(inter_list, list_of_points[:2] + list_of_points[-2:])
            logger.debug("Removing point %f, %f" % (p[0], p[1]))

            list_of_points.remove(p)

            if len(list_of_points) > 3:
                s1 = LineString([list_of_points[-1], list_of_points[0]])
                s2 = LineString(list_of_points[1:-1])
                inter_list = compute_segment_intersection(s1, s2)

    list_of_points.append(list_of_points[0])

    return list_of_points


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


def get_ogr_polygon_from_ring_list(multi_ring_list) :
    """
    return ogr geometry polygon from multi ring list

    :param multi_ring_list: list of rings
    :type multi_ring_list: list of list of tuple of 2 floats

    :return: lake polygon
    :rtype: OGRPolygon
    """
    # 1. build a list of ogr geometry "polygon" and "linear ring"
    list_of_ogr_poly = []
    list_of_ogr_ring = []
    for ring in multi_ring_list :
        lake_ring_geom = ogr.Geometry(ogr.wkbLinearRing)
        for (lon, lat) in ring :
            lake_ring_geom.AddPoint(lon, lat)
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
            if i == j :
                continue
            if poly_i.Within(poly_j):
                if j not in list_of_outer_ring:
                    list_of_outer_ring.append(j)
                list_of_inner_ring.append((j,i))
                is_i_inner = True
        if is_i_inner == False :
            list_of_outer_ring.append(i)

    # 3. Build out geometry multipolygon, with outer / inner rings in the right order.
    multi_poly_geom = ogr.Geometry(ogr.wkbMultiPolygon)
    for i in list_of_outer_ring:
        poly_geom = list_of_ogr_poly[i][0]
        for (out_ring, in_ring) in list_of_inner_ring:
            if i == out_ring:
                poly_geom.AddGeometry(list_of_ogr_ring[in_ring])
        multi_poly_geom.AddGeometry(poly_geom)

    return multi_poly_geom


#######################################


def get_concave_hull_from_cgal_triangulation(in_v_long, in_v_lat, in_range, in_azimuth, in_nb_pix_range):
    """
    Compute the concave hull of a set of points using CGAL Delauney triangulation.

    :param in_v_long: longitudes of points for which to compute the hull
    :type in_v_long: 1D-array of float
    :param in_v_lat: latitudes of points for which to compute the hull
    :type in_v_lat: 1D-array of float
    :param in_range: range of points for which to compute the hull
    :type in_range: 1D-array of int
    :param in_azimuth: azimuth of points for which to compute the hull
    :type in_azimuth: 1D-array of int
    :param in_nb_pix_range: maximal number of pixel in range
    :type in_nb_pix_range: int

    :return the hull of the input set of points
    :rtype: OGRMultiPolygon
    """
    logger = logging.getLogger("my_hull")

    try:
        from CGAL import CGAL_Alpha_shape_2
        from CGAL.CGAL_Kernel import Point_2, Polygon_2, Vector_2
    except:
        logger.error("CGAL library must be install to use this option, use an other HULL_METHOD")
        raise

    # Transform lonlat coordinates to UTM coordinates
    coords = np.zeros((in_v_long.size, 2))
    coords[:, 0], coords[:, 1], utm_epsg_code = my_tools.get_utm_coords_from_lonlat(in_v_long, in_v_lat)

    # Separate swath in two zone : near and far range. This step is useful to adjust alpha regarding range sampling.
    near_rg_idx = (np.where(in_range < (in_nb_pix_range / 2 + 10)))
    far_rg_idx = (np.where(in_range > (in_nb_pix_range / 2 - 10)))

    # Near range processing
    alpha_near_rg, dist_near_rg = evaluate_alpha_from_x_pixc_geolocation(coords[:, 0][near_rg_idx], in_range[near_rg_idx], in_azimuth[near_rg_idx])
    if dist_near_rg :
        concave_hull_utm_near_rg = alpha_shape_with_cgal(coords[near_rg_idx], alpha_near_rg)
        concave_hull_near_rg = my_tools.get_lon_lat_polygon_from_utm(concave_hull_utm_near_rg, utm_epsg_code)
        logger.info("Alpha far rg value is %d and X var (mean + 2std) = %f " % (alpha_near_rg, dist_near_rg))
    else :
        concave_hull_near_rg = MultiPolygon()

    # Far range processing
    alpha_far_rg, dist_far_rg = evaluate_alpha_from_x_pixc_geolocation(coords[:, 0][far_rg_idx], in_range[far_rg_idx], in_azimuth[far_rg_idx])
    if dist_far_rg :
        concave_hull_utm_far_rg = alpha_shape_with_cgal(coords[far_rg_idx], alpha_far_rg)
        concave_hull_far_rg = my_tools.get_lon_lat_polygon_from_utm(concave_hull_utm_far_rg, utm_epsg_code)
        logger.info("Alpha far rg value is %d and X var (mean + 2std) = %f " % (alpha_far_rg, dist_far_rg))
    else :
        concave_hull_far_rg = MultiPolygon()

    concave_hull = unary_union([concave_hull_far_rg, concave_hull_near_rg])

    retour = ogr.CreateGeometryFromWkt(concave_hull.wkt)

    logger.info("Polygon is valid : %r" %(retour.IsValid()))

    if retour.IsEmpty():
        retour = get_convex_hull(in_v_long, in_v_lat)

    return retour


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

    return alpha, x_2sigma


def get_angle(p, q, r):
    """
    Compute angle between 3 points
    Function coming from https://github.com/pgRouting/pgrouting/blob/master/src/alpha_shape/alpha.c

    :param p: Point
    :param q: Point
    :param r: Point
    :return: angle
    """
    logger = logging.getLogger("my_hull")

    try:
        from CGAL import CGAL_Alpha_shape_2
        from CGAL.CGAL_Kernel import Point_2, Polygon_2, Vector_2
    except:
        logger.error("CGAL library must be install to use this option, use an other HULL_METHOD")
        raise
        
    v1 = Vector_2(q, p)
    v2 = Vector_2(q, r)

    crossproduct = v1.x() * v2.y() - v1.y() * v2.x();
    dot = v1.x() * v2.x() + v1.y() * v2.y();
    angle = math.atan2(crossproduct, dot)
    if angle < 0.0:
        angle += 2 * np.pi
    return angle


def find_rings(vertices, edges):
    """
    Identify rings from a list of vertices and edges

    :param vertices: List of vertices
    :param edges: List of edges
    :return: rings
    """
    logger = logging.getLogger("my_hull")

    try:
        from CGAL import CGAL_Alpha_shape_2
        from CGAL.CGAL_Kernel import Point_2, Polygon_2, Vector_2
    except:
        logger.error("CGAL library must be install to use this option, use an other HULL_METHOD")
        raise
        
    # Create a dataframe for segments
    df_edges = pd.DataFrame(data=[[edge.source(),
                                   edge.target(),
                                   True] for edge in edges],
                            columns=['source', 'target', 'unused'])
    # Create a dataframe for vertices
    index = pd.MultiIndex.from_tuples(vertices, names=['x', 'y'])
    df_vertices = pd.DataFrame(columns=['segments'],
                               index=index)
    # Remove duplicates index
    df_vertices = df_vertices[~df_vertices.index.duplicated(keep='first')]
    # Initialize segment list
    df_vertices['segments'] = ''
    df_vertices['segments'] = df_vertices['segments'].apply(list)
    # Add segment list
    for row in df_edges.itertuples():
        index = getattr(row, "Index")
        source = getattr(row, "source")
        df_vertices.loc[(source.x(), source.y())]['segments'].append(index)
    # Create polygons
    rings = []
    nb = df_edges.shape[0]
    while not df_edges[df_edges.unused].empty:
        # New polygon
        ring = Polygon_2()
        # Get an edge
        ind = 0
        current_edge = df_edges[df_edges.unused].iloc[0]
        start = current_edge['source']
        current_start = start
        current_end = current_edge['target']
        df_edges.at[current_edge.name, 'unused'] = False
        ring.push_back(start)
        ring.push_back(current_end)
        while current_end != start and not df_edges[df_edges.unused].empty:
            ind += 1
            # Find next possible edges
            next_indexes = df_vertices.loc[(current_end.x(), current_end.y())]['segments']
            if len(next_indexes) == 1:
                next_index = next_indexes[0]
            elif len(next_indexes) > 1:
                next_angles = []
                for j in next_indexes:
                    target = df_edges.loc[j]["target"]
                    angle = get_angle(current_start, current_end, target)
                    next_angles.append((angle, j))
                sorted_next_angles = sorted(next_angles, key=lambda tup: tup[0], reverse=False)
                next_index = sorted_next_angles[0][1]
                next_indexes.remove(next_index)
                df_vertices.at[(current_end.x(), current_end.y()), 'segments'] = next_indexes
            current_edge = df_edges.loc[next_index]
            current_start = current_edge['source']
            current_end = current_edge['target']
            df_edges.at[next_index, 'unused'] = False
            ring.push_back(current_end)
            if ind == nb:
                raise Exception("Too many")
        rings.append(ring)
    return rings


def find_polygons(rings):
    """
    Identify polygons from rings

    :param rings: List of rings
    :return: polygons
    """

    logger = logging.getLogger("my_hull")

    try:
        from CGAL import CGAL_Alpha_shape_2
        from CGAL.CGAL_Kernel import Point_2, Polygon_2, Vector_2
    except:
        logger.error("CGAL library must be install to use this option, use an other HULL_METHOD")
        raise
        
    polygons = []
    # Create dataframe
    df = pd.DataFrame(data=rings, columns=['ring'])
    # Convert to shapely linearring and compute area
    df['ring'] = df.apply(lambda row: LinearRing([[vertex.x(), vertex.y()] for vertex in row.ring.vertices()]), axis=1)
    df['area'] = df.apply(lambda row: Polygon(row.ring).area, axis=1)
    # Sort by area
    df = df.sort_values(by=['area'], ascending=False)

    # Find interiors
    df['unused'] = True
    while not df[df.unused].empty:
        # Get an ring (with the larger surface)
        current = df[df.unused].iloc[0]
        ext = current.ring
        polygon = Polygon(ext)
        df.at[current.name, 'unused'] = False
        indexes = []
        ints = []
        for index, row in df[df.unused].iterrows():
            if polygon.contains(row.ring):
                ints.append(row.ring)
                indexes.append(index)
        for index in indexes:
            df.at[index, 'unused'] = False
        if len(ints) > 0:
            polygons.append(Polygon(ext, ints))
        else:
            polygons.append(Polygon(ext))
    return polygons


def alpha_shape_with_cgal(coords, alpha):
    """
    Compute the alpha shape of a set of points.
    Retrieved from http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/

    :param coords : Coordinates of points
    :param alpha: List of alpha values to influence the gooeyness of the border. Smaller numbers don't fall inward as much as larger numbers. 
    Too large, and you lose everything!
    :return: Shapely.MultiPolygons which is the hull of the input set of points
    """

    logger = logging.getLogger("my_hull")

    try:
        from CGAL import CGAL_Alpha_shape_2
        from CGAL.CGAL_Kernel import Point_2, Polygon_2, Vector_2
    except:
        logger.error("CGAL library must be install to use this option, use an other HULL_METHOD")
        raise
        
    alpha_value = np.mean(alpha)
    # Convert to CGAL point
    points = [Point_2(pt[0], pt[1]) for pt in coords]
    # Compute alpha shape
    a = CGAL_Alpha_shape_2.Alpha_shape_2()
    a.make_alpha_shape(points)
    a.set_alpha(alpha_value)
    a.set_mode(CGAL_Alpha_shape_2.REGULARIZED)
    alpha_shape_edges = [a.segment(it) for it in a.alpha_shape_edges()]

    alpha_shape_vertices = [(vertex.point().x(), vertex.point().y()) for vertex in a.alpha_shape_vertices()]

    if len(alpha_shape_vertices) == 0:
        poly_shp = Polygon()
    else:
        rings = find_rings(alpha_shape_vertices, alpha_shape_edges)
        polygons_list = find_polygons(rings)
    
        if not polygons_list:
            poly_shp = Polygon()
        elif len(polygons_list) > 1 :
            poly_shp = MultiPolygon(polygons_list)
        else :
            poly_shp = Polygon(polygons_list[0])

    return poly_shp

