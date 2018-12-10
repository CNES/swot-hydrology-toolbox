# -*- coding: utf8 -*-
"""
.. module:: my_hull.py
    :synopsis: deal with hull computation (=lake boundary computation)
    Created on 08/27/2018

.. moduleauthor:: Claire POTTIER (CNES DSO/SI/TR) - Cécile CAZALS (CS)

Copyright (c) 2018 CNES. All rights reserved.
"""

import math
import numpy as np
from osgeo import ogr
from scipy.spatial import Delaunay
import shapely.geometry as geometry
from shapely.ops import cascaded_union
from skimage.measure import find_contours
import logging

import cnes.common.lib.my_tools as my_tools
import cnes.common.lib_lake.locnes_variables as my_var
import cnes.common.serviceError as serviceError


def computeLakeBoundaries(IN_vLong, IN_vLat, IN_range, IN_azimuth, IN_nb_pix_range):
    """
    Compute the hull of a set of points determined by their coordinates given in input parameters

    :param IN_vLong: longitudes of points for which to compute the hull
    :type IN_vLong: 1D-array of float
    :param IN_vLat: latitudes of points for which to compute the hull
    :type IN_vLat: 1D-array of float
    :param IN_range: range of points for which to compute the hull
    :type IN_range: 1D-array of int
    :param IN_azimuth: azimuth of points for which to compute the hull
    :type IN_azimuth: 1D-array of int
    :param IN_nb_pix_range: maximal number of pixel in range
    :type IN_nb_pix_range: int

    :return the hull of the input set of points
    :rtype: OGRMultiPolygon
    """
    logger = logging.getLogger("my_hull")
    logger.debug("[my_hull] == computeLakeBoundaries ==")

    # Nb pixels
    nb_pix = IN_vLong.size

    if my_var.HULL_METHOD == 0:  # 1 - CONVEX HULL
        
        # 1.1 - Group pixels in a OGRGeometryCollection
        pixc_pts = ogr.Geometry(ogr.wkbGeometryCollection)
        for indp in np.arange(nb_pix):
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(IN_vLong[indp], IN_vLat[indp])
            pixc_pts.AddGeometry(point)

        # 1.2 - Compute the hull of these points
        return pixc_pts.ConvexHull()

    elif math.floor(my_var.HULL_METHOD) == 1:  # 1 - CONCAV HULL - Delaunay triangulation method

        # 1.1 - Together coordinates in a 2D-array
        coords = np.zeros((nb_pix, 2))
        coords[:, 0] = IN_vLong
        coords[:, 1] = IN_vLat

        # 1.2 - Compute alpha shape
        
        if my_var.HULL_METHOD == 1.1:  # Without alpha parameter
            concave_hull = getConcavHull_bis(coords)
        else:  # With alpha parameter
            alpha = (2000 + 2000 * IN_range / IN_nb_pix_range).astype('int')  # alpha parameter ranges from 2000 to 4000 following the range index
            concave_hull = alpha_shape(coords, alpha)

        # 1.3 - Convert Shapely geometry to OGRPolygon
        return ogr.CreateGeometryFromWkb(concave_hull.wkb)

    elif my_var.HULL_METHOD == 2:  # 2 - CONCAV HULL - Radar vectorisation method

        radar_vect_concav_hull = getConcaveHullFromRadarVectorisation(IN_range, IN_azimuth, IN_vLong, IN_vLat)

        return radar_vect_concav_hull

    else:
        message = "Concave hull computation method not understood"
        raise serviceError.ProcessingError(message, logger)


#######################################


def alpha_shape(IN_coords, IN_alpha):
    """
    Compute the alpha shape (concave hull) of a set of points.
    
    :param IN_coords: set of points coordinates
    :type IN_coords: 2D-array of size (nb_pixels, 2=lon/lat)
    :param IN_alpha: alpha value to influence the gooeyness of the border. Smaller numbers don't fall inward as much as larger numbers. Too large, and you lose everything!
    :type IN_alpha: 1D-array of int
    
    :return Shapely.MultiPolygons which is the hull of the input set of points
    
    Retrieved from http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
    """
    
    # Number of points
    nb_pts = IN_coords.shape[0]

    # 0 - Particular case of nb points <= 3: return convex hull
    if nb_pts <= 3:
        # 0.1 - Init Shapely set of pointsn and Aggregate coordinates to the set of points
        points = [geometry.point.Point(IN_coords[indp, 0], IN_coords[indp, 1], 0) for indp in np.arange(nb_pts)]
        # 0.2 - Return convex hull
        return geometry.MultiPoint(list(points)).convex_hull

    # 1 - Compute Delaunay triangulation
    tri = Delaunay(IN_coords)  # tri = specific object with attributes (cf. https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.Delaunay.html)

    # 2 - Select triangles following it's shape
    # 2.1 - Compute CircumRatio for each triangle
    circum_r = [getCircumRatio(IN_coords[ia], IN_coords[ib], IN_coords[ic]) for ia, ib, ic in tri.vertices]

    # 2.2 - Compute mean alpha parameter for each triangle
    mean_alpha = [1.0 / np.mean([IN_alpha[ia], IN_alpha[ib], IN_alpha[ic]]) for ia, ib, ic in tri.vertices]

    # 2.3 - Select triangles
    triangles_selected = np.where(np.array(circum_r) < mean_alpha)

    # 2.4 - Compute a list of shapely polygon correpsonding to selected triangles
    list_triangle = [geometry.Polygon([(IN_coords[ia][0], IN_coords[ia][1]), (IN_coords[ib][0], IN_coords[ib][1]), (IN_coords[ic][0], IN_coords[ic][1])]) for ia, ib, ic in tri.vertices[triangles_selected]]

    # 3 - Union of selected triangles
    triangle_union = cascaded_union(list_triangle)

    return triangle_union


def getCircumRatio(pa, pb, pc):
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
    except:
        return 1

    # Circumscribing circle radius
    circum_r = 0
    if area != 0:
        circum_r = a * b * c / (4.0 * area)

    return circum_r


#######################################


def getConcavHull_bis(IN_coords):
    """
    Compute the Delaunay triangulation and return the concave hull of cloud point IN_coords

    :param IN_coords: set of points coordinates
    :type IN_coords: 2D-array of float ; size = (nb_pixels, 2=lon/lat)

    :return the hull of the input set of points
    :rtype: Shapely.MultiPolygons

    Partially retrieved from http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
    """

    # Number of points
    nb_pts = IN_coords.shape[0]

    # 0 - Particular case of nb points <= 3: return convex hull
    if nb_pts <= 3:
        # 0.1 - Init Shapely set of points and aggregate its coordinates
        points = [geometry.point.Point(IN_coords[indp, 0], IN_coords[indp, 1], 0) for indp in np.arange(nb_pts)]
        # 0.2 - Return convex hull
        return geometry.MultiPoint(list(points)).convex_hull

    # 1 - Compute Delaunay triangulation
    tri = Delaunay(IN_coords)  # tri = specific object with attributes (cf. https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.Delaunay.html)

    # 2 - Select triangles following their shape
    
    # 2.1 - Get maximal length of each triangle
    maxseg = np.array([getMaxSegment(IN_coords[ia], IN_coords[ib], IN_coords[ic]) for ia, ib, ic in tri.vertices])

    # 2.2 - Get median length of maximal segments
    maxseg_median = np.median(maxseg)

    # 2.3 - Select triangles
    triangles_selected_idx = np.where(maxseg < maxseg_median * 3)

    # 2.4 - Compute a list of shapely polygon correpsonding to selected triangles
    list_triangle = [geometry.Polygon([(IN_coords[ia][0], IN_coords[ia][1]), (IN_coords[ib][0], IN_coords[ib][1]), (IN_coords[ic][0], IN_coords[ic][1])]) for ia, ib, ic in tri.vertices[triangles_selected_idx]]

    # 3 - Union of selected triangles
    triangle_union = cascaded_union(list_triangle)

    # 4 - Deletion of some artifacts related to Delaunay triangulation: holes may appears in some cases
    if triangle_union.type == "Polygon":

        # Get holes in concerned polygon
        holes_in_poly = [list(poly.coords) for poly in list(triangle_union.interiors)]

        for poly_holes in list(triangle_union.interiors):
            # If the hole contains only 3 points (triangle) remove the hole
            if len(poly_holes.coords) == 4:
                triangle_union = geometry.Polygon(triangle_union.exterior.coords, holes_in_poly.remove(list(poly_holes.coords)))

    elif triangle_union.type == "MultiPolygon":

        multipolygon = []

        for polygon in triangle_union:

            # Get holes in concerned polygon
            holes_in_poly = [list(poly.coords) for poly in list(polygon.interiors)]

            for poly_holes in list(polygon.interiors):
                # If the hole contains only 3 points (triangle) remove the hole
                if len(poly_holes.coords) == 4:
                    multipolygon.append(geometry.Polygon(polygon.exterior.coords, holes_in_poly.remove(list(poly_holes.coords))))

        triangle_union = cascaded_union(multipolygon)

    return triangle_union


def getMaxSegment(pa, pb, pc):
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


def getConcaveHullFromRadarVectorisation(IN_range, IN_azimuth, IN_vLong, IN_vLat):
    """
    Compute the concave hull of a set of points using radar vectorisation
    
    :param IN_range: range of points
    :type IN_range: 1D-array of int
    :param IN_azimuth: azimuth of points
    :type IN_azimuth: 1D-array of int
    :param IN_vLong: longitude of points
    :type IN_vLong: 1D-array of float
    :param IN_vLat: latitude of points
    :type IN_vLat: 1D-array of float
    
    :return: the hull of the input set of points
    :rtype: Shapely.MultiPolygons
    """
    logger = logging.getLogger("my_hull")
    logger.debug("[my_hull] == getConcaveHullFromRadarVectorisation ==")

    # 1 - Get relative range and azimuth (1st pixel = 1)
    lake_x = IN_range - np.min(IN_range) + 1
    lake_y = IN_azimuth - np.min(IN_azimuth) + 1

    # 2 - Get image (1 pixel around lake)
    lake_img = my_tools.computeBinMat(np.max(lake_x) + 2, np.max(lake_y) + 2, lake_x, lake_y)

    # 3 - Compute boundaries (there might be more than one if there are islands in lake)
    lake_contour = find_contours(lake_img, 0.99999999999)

    # 4 - Round contour range and azimuth coordinates to units

    # 4.1 - Get boundary as a list of integers (since they are indices)
    lake_contour_int = []
    for contour in lake_contour:
        lake_contour_int.append(np.around(contour, 0))

    # 4.2 - Convert (azimuth, range) into polygon
    lake_poly = ogr.Geometry(ogr.wkbPolygon)
    lake_ring = None
    for contour in lake_contour_int:  # Loop over contours
        lake_ring = ogr.Geometry(ogr.wkbLinearRing)
        for (y, x) in contour[::2]:  # Look over azimuth and range indices
            if np.where(np.logical_and(lake_x == x, lake_y == y))[0]:
                point_idx = np.where(np.logical_and(lake_x == x, lake_y == y))[0][0]
                lon = IN_vLong[point_idx]
                lat = IN_vLat[point_idx]
                lake_ring.AddPoint(lon, lat)
        lake_ring.AddPoint(lake_ring.GetPoint(0)[0], lake_ring.GetPoint(0)[1])
        lake_poly.AddGeometry(lake_ring)

    return lake_poly
