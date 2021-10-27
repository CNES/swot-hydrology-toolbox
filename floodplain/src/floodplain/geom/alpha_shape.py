# -*- coding: utf8 -*-
'''
Library with alpha_shape implementation
based on cascaded_union or CGAL

See http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
See https://www.cgal.org/

Copyright (c) 2018, CNES
'''

import sys
import numpy as np
from shapely.geometry import Point, MultiPoint, Polygon, MultiPolygon, LinearRing
from shapely.ops import cascaded_union
from scipy.spatial import Delaunay
import math
from CGAL import CGAL_Alpha_shape_2
from CGAL.CGAL_Kernel import Point_2, Segment_2, Polygon_2, Vector_2
from typing import List  , Tuple
import pandas as pd

def alpha_shape_with_cascaded_union(coords: np.ndarray, 
                                    alpha: np.ndarray) -> List[Polygon]:
    '''
    Compute the alpha shape of a set of points.
    Retrieved from http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
    
    :param coords : Coordinates of points 
    :param alpha: List of alpha values to influence the
        gooeyness of the border. Smaller numbers
        don't fall inward as much as larger numbers.
        Too large, and you lose everything!
    :return Shapely.MultiPolygons which is the hull of the input set of points
    '''
    # Number of points
    nb_pts = coords.shape[0]
    # 0 - Particular case of nb points <= 3: return convex hull
    if ( nb_pts <= 3 ):
        # When you have a triangle, there is no sense
        # in computing an alpha shape.
        # 0.1 - Init Shapely set of points
        points = []
        # 0.2 - Aggregate coordinates to the set of points
        for indp in np.arange(nb_pts):
            points.append(Point(coords[indp,0],coords[indp,1],0))
        # 0.3 - Return convex hull
        return MultiPoint(list(points)).convex_hull

    # Compute Delaunay's triangulation
    # tri = specific object with attributes 
    # (cf. https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.Delaunay.html)
    tri = Delaunay(coords)

    list_triangle = []
    # loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.vertices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]
        # Lengths of sides of triangle
        a = math.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
        b = math.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
        c = math.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)
        # Semiperimeter of triangle
        s = (a + b + c)/2.0
        # Area of triangle by Heron's formula
        area = math.sqrt(s*(s-a)*(s-b)*(s-c))
        circum_r = 0.0
        if area > 1.e-8:
            circum_r = a*b*c/(4.0*area)
        # Compute alpha criteria
        mean_alpha = np.mean([alpha[ia], alpha[ib], alpha[ic]])
        # Here's the radius filter.
        if circum_r < mean_alpha:
            gg = [(pa[0], pa[1]),(pb[0], pb[1]),(pc[0], pc[1])]
            tt = Polygon(gg)
            list_triangle.append(tt)  
    result = cascaded_union(list_triangle)
    polygons = []
    if result.geom_type == "Polygon":
        polygons = [result]
    elif result.geom_type == "MultiPolygon":
        polygons = [polygon for polygon in result] 
    return polygons


def get_angle(p: Point_2, q: Point_2, r: Point_2) -> float:
    '''
    Compute angle between 3 points
    tan Theta = || cross_product(v1,v2) || / dot_product(v1,v2)

    :param p: Point
    :param q: Point
    :param r: Point
    :return angle
    '''
    v1 = Vector_2(q, p)
    v2 = Vector_2(q, r)
    magnitude_cross_product = v1.x() * v2.y() - v1.y() * v2.x()
    dot_product = v1.x() * v2.x() + v1.y() * v2.y()
    angle = math.atan2(magnitude_cross_product, dot_product)
    if angle < 0.0:
        angle += 2 * np.pi
    return angle


def find_rings(vertices: List[Tuple],edges: List[Segment_2]) -> List[Polygon_2]:
    '''
    Identify rings from a list of vertices and edges

    :param vertices: List of vertices
    :param edges: List of edges
    :return rings
    '''
    # Create a dataframe for edges (list of segments)
    # For each edge the column unused spêcifies if the edge is still unused
    df_edges = pd.DataFrame(data=[[edge.source(),
                                   edge.target(),
                                  True] for edge in edges],
                            columns=['source','target','unused'])
    # Create a dataframe for vertices 
    # For each vertex, it specifies the position of the vertex and the list of 
    # segments in which the vertex is used as source point 
    index = pd.MultiIndex.from_tuples(vertices, names=['x', 'y'])
    df_vertices = pd.DataFrame(columns=['segments'], 
                               index=index)
    # Remove duplicates index
    df_vertices = df_vertices[~df_vertices.index.duplicated(keep='first')]
    # Initialize segment list
    df_vertices['segments'] = ''
    df_vertices['segments'] = df_vertices['segments'].apply(list)
    # Add segment list
    for row in  df_edges.itertuples():
        index = getattr(row,"Index")
        source = getattr(row, "source")
        df_vertices.loc[(source.x(),source.y())]['segments'].append(index)
    # Create polygons
    rings = []
    nb = df_edges.shape[0]
    nb_poly = 0
    while not df_edges[df_edges.unused].empty:
        # Create a counter to avoid infinite loop
        ind = 0
        # Initialize the creation of a new polygon
        ring = Polygon_2() 
        # Get an edge
        current_edge = df_edges[df_edges.unused].iloc[0]
        start = current_edge['source'] # Start point of the polygon
        current_start = start
        current_end = current_edge['target']
        # Tag the edge as used
        df_edges.at[current_edge.name,'unused'] = False
        # From the edge get the two first points of the polygon
        ring.push_back(start)
        ring.push_back(current_end)
        # Remove the current edge from the segment list
        mask = df_vertices.segments.apply(lambda x: int(current_edge.name) in x)
        df_vertices[mask] = df_vertices.segments.apply(lambda L: [x for x in L if x!=current_edge.name])
        # Loop to create the polygon
        while current_end != start and not df_edges[df_edges.unused].empty:
            ind += 1
            # Find next possible edges
            next_indexes = df_vertices.loc[(current_end.x(),current_end.y())]['segments']
            # If one possibity
            if len(next_indexes) == 1:
                next_index = next_indexes[0]
            # If there are several possiblties compute the angle to find the netx edges
            elif len(next_indexes) > 1:
                next_angles = []
                for j in next_indexes:
                    target = df_edges.loc[j]["target"]
                    angle = get_angle(current_start, current_end, target)
                    next_angles.append((angle, j))
                sorted_next_angles = sorted(next_angles, key=lambda tup: tup[0],reverse=False)
                next_index = sorted_next_angles[0][1]
                next_indexes.remove(next_index)
                df_vertices.at[(current_end.x(),current_end.y()),'segments'] = next_indexes
            # When the next edge is found, add the target point to the polygon
            current_edge = df_edges.loc[next_index]
            current_start = current_edge['source']
            current_end = current_edge['target']
            df_edges.at[next_index,'unused'] = False
            ring.push_back(current_end)
            if ind == nb:
                raise Exception("Too many")
        # End of polygon creation
        # Add the polygon to the polygon list
        rings.append(ring)
    return rings

def find_polygons(rings: List[Polygon_2]) -> List[Polygon]:
    '''
    Identify polygons from rings

    :param rings: List of rings
    :return polygons
    '''
    polygons = []
    # Create dataframe
    df = pd.DataFrame(data=rings,columns=['ring'])
    # Convert to shapely linearring and compute area 
    df['ring'] = df.apply(lambda row: LinearRing([[vertex.x(),vertex.y()] for vertex in row.ring.vertices()]), axis = 1)
    df['area'] = df.apply(lambda row: Polygon(row.ring).area, axis = 1)
    # Sort by area
    df = df.sort_values(by=['area'], ascending=False)

    # Find interiors
    df['unused'] = True
    while not df[df.unused].empty:
        # Get an ring (with the larger surface)
        current = df[df.unused].iloc[0]
        ext = current.ring
        polygon = Polygon(ext)
        df.at[current.name,'unused'] = False
        indexes = []
        ints = []
        for index, row in df[df.unused].iterrows():
            if polygon.contains(row.ring):
                ints.append(row.ring)
                indexes.append(index)
        for index in indexes:
            df.at[index,'unused'] = False
        if len(ints) > 0:
            polygons.append(Polygon(ext, ints))
        else:
            polygons.append(Polygon(ext))
    return polygons

# ~ def alpha_shape_with_cgal(coords: np.ndarray,
                                    # ~ alpha: np.array) -> List[Polygon]:
    # ~ '''
    # ~ Compute the alpha shape of a set of points.
    # ~ Retrieved from http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
    
    # ~ :param coords : Coordinates of points 
    # ~ :param alpha: List of alpha values to influence the
        # ~ gooeyness of the border. Smaller numbers
        # ~ don't fall inward as much as larger numbers.
        # ~ Too large, and you lose everything!
    # ~ :return Shapely.MultiPolygons which is the hull of the input set of points
    # ~ '''
    # ~ alpha_value = np.mean(alpha)
    # ~ # Convert to CGAL point
    # ~ points = [Point_2(pt[0],pt[1]) for pt in coords]
    # ~ # Compute alpha shape
    # ~ a = CGAL_Alpha_shape_2.Alpha_shape_2()
    # ~ a.make_alpha_shape(points)
    # ~ a.set_alpha(alpha_value)
    # ~ a.set_mode(CGAL_Alpha_shape_2.REGULARIZED)
    # ~ alpha_shape_edges = [ a.segment(it) for it in a.alpha_shape_edges() ]

    # ~ alpha_shape_vertices = [ (vertex.point().x(),vertex.point().y()) for vertex in a.alpha_shape_vertices() ]

    # ~ if len(alpha_shape_vertices) == 0:
        # ~ return []
    # ~ rings = find_rings(alpha_shape_vertices, alpha_shape_edges)
    # ~ polygons = find_polygons(rings)
    # ~ return polygons


def alpha_shape_with_cgal(coords, alpha):
    """
    Compute the alpha shape of a set of points.
    Retrieved from http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/

    :param coords : Coordinates of points
    :param alpha: List of alpha values to influence the gooeyness of the border. Smaller numbers don't fall inward as much as larger numbers. 
    Too large, and you lose everything!
    :return: Shapely.MultiPolygons which is the hull of the input set of points
    """
        
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
