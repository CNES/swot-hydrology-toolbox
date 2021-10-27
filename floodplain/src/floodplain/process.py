# -*- coding: utf8 -*-
'''
Create ply file from a list of pixel cloud files

Copyright (c) 2018, CNES
'''

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


from floodplain.utils.spatial import compute_binary_mask, convert_to_utm_point
from floodplain.geom.alpha_shape import alpha_shape_with_cgal, alpha_shape_with_cascaded_union

WATER_LABEL = 3
WATER_NEAR_LAND_LABEL = 4
DARK_WATER_LABEL = 23
DARK_WATER_NEAR_LAND_LABEL = 24

WATER_LABELS = [WATER_LABEL,
                WATER_NEAR_LAND_LABEL,
                DARK_WATER_LABEL,
                DARK_WATER_NEAR_LAND_LABEL]


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
                          (pixelcloud.classification == DARK_WATER_NEAR_LAND_LABEL) ].copy()


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
    water['geometry_utm'] = water.apply(lambda row: convert_to_utm_point(row['latitude'],row['longitude']), axis=1)    
    water['utm_x'] = water.apply(lambda row: row['geometry_utm'].x, axis=1)
    water['utm_y'] = water.apply(lambda row: row['geometry_utm'].y, axis=1)
    tree = spatial.cKDTree(np.array([ [pt.x,pt.y] for pt in water["geometry_utm"].values ])) 
    # Find all points within distance r of point(s) x
    water["neighbors"] = water.apply(lambda row: len(tree.query_ball_point([row['geometry_utm'].x,row['geometry_utm'].y], factor*np.sqrt((row['range_spacing'])**2+(slc_along_track_resolution)**2))), axis=1)
    # Keep points with at least 4 neighbors
    return water.loc[(water.neighbors > nb_neighbors)].copy()


# ~ def compute_alpha_shape(water: gpd.GeoDataFrame, nb_pix_range: int, method: str = "cgal") -> gpd.GeoDataFrame:
    # ~ '''
    # ~ Copute alpha_shape

    # ~ :param water: Water Dataframe

    # ~ :return points: Extracted water points 
    # ~ '''
    # ~ data = np.array(water[['utm_x','utm_y']].values)

    # ~ polygons = []
    # ~ if method == "cascaded_union":
        # ~ polygons = alpha_shape_with_cascaded_union(data, alpha) 
    # ~ elif method == "cgal":
                
        # ~ in_v_long = water['longitude'].values
        # ~ in_v_lat = water['latitude'].values
        # ~ in_range = water['range_index'].values
        # ~ in_azimuth = water['azimuth_index'].values
        # ~ in_nb_pix_range = nb_pix_range
        
        # ~ coords = np.zeros((in_v_long.size, 2))
        # ~ coords[:, 0], coords[:, 1] = water['utm_x'].values, water['utm_y'].values

        # ~ # Separate swath in two zone : near and far range. This step is useful to adjust alpha regarding range sampling.
        # ~ near_rg_idx = (np.where(in_range < (in_nb_pix_range / 2 + 10)))
        # ~ far_rg_idx = (np.where(in_range > (in_nb_pix_range / 2 - 10)))
        
        # ~ # Near range processing
        # ~ alpha_near_rg, dist_near_rg = evaluate_alpha_from_x_pixc_geolocation(coords[:, 0][near_rg_idx], in_range[near_rg_idx], in_azimuth[near_rg_idx])
        # ~ if dist_near_rg :
            # ~ concave_hull_near_rg = alpha_shape_with_cgal(coords[near_rg_idx], alpha_near_rg)
        # ~ else :
            # ~ concave_hull_near_rg = MultiPolygon()

        # ~ # Far range processing
        # ~ alpha_far_rg, dist_far_rg = evaluate_alpha_from_x_pixc_geolocation(coords[:, 0][far_rg_idx], in_range[far_rg_idx], in_azimuth[far_rg_idx])
        # ~ if dist_far_rg :
            # ~ concave_hull_far_rg = alpha_shape_with_cgal(coords[far_rg_idx], alpha_far_rg)
        # ~ else :
            # ~ concave_hull_far_rg = MultiPolygon()

        # ~ concave_hull = unary_union([concave_hull_far_rg, concave_hull_near_rg])
        
        # ~ if concave_hull.type == MultiPolygon:
            # ~ polygons = list(concave_hull)
        # ~ else:
            # ~ polygons = concave_hull
    
    # ~ else: 
        # ~ raise Exception("Method not found")
    # ~ return polygons 


       
def compute_alpha_shape(water: gpd.GeoDataFrame, nb_pix_range: int, method: str = "cgal") -> gpd.GeoDataFrame:
    '''
    Copute alpha_shape

    :param water: Water Dataframe

    :return points: Extracted water points 
    '''
    data = np.array(water[['utm_x','utm_y']].values)

    polygons = []
    if method == "cascaded_union":
        polygons = alpha_shape_with_cascaded_union(data, alpha) 
    elif method == "cgal":
                
        in_v_long = water['longitude'].values
        in_v_lat = water['latitude'].values
        in_range = water['range_index'].values
        in_azimuth = water['azimuth_index'].values
        in_nb_pix_range = nb_pix_range
        
        coords = np.zeros((in_v_long.size, 2))
        coords[:, 0], coords[:, 1] = water['utm_x'].values, water['utm_y'].values

        # Separate swath in two zone : near and far range. This step is useful to adjust alpha regarding range sampling.
        
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

    get_point_from_coord = lambda coord: tuple(data.loc[coord,['longitude','latitude','elevation','range_index','azimuth_index']].values)    
    
    points = [get_point_from_coord(coord) for boundary in boundaries for coord in boundary.coords ]

    # ~ get_mean_height = lambda coord: np.array(data.loc[coord,['wse']]).flatten()
    # ~ mean_wse = []
    # ~ for boundary in boundaries:
        # ~ toto = [get_mean_height(coord) for coord in boundary.coords]
        # ~ mean_wse.append(np.mean(toto))
    
    # Convert to dataframe
    df = pd.DataFrame(points,columns=['longitude','latitude','elevation','range','azimuth'])
    lambdafunc = lambda row: pd.Series([*utm.from_latlon(row['latitude'],
                                                   row['longitude'])[0:2],
                                  row['elevation'],
                                      ]) 
    df[['x','y','z']] = df.apply(lambdafunc,axis=1)
    df = df.astype(dtype={'longitude':'double','latitude':'double', 'elevation':'double',
                          'range':'int','azimuth':'int',
                         'x':'double','y':'double', 'z':'double'})
    

    geom = [Point(x,y) for x, y in zip(df['longitude'], df['latitude'])]

    # ~ return gpd.GeoDataFrame(df, geometry=geom), mean_wse
    return gpd.GeoDataFrame(df, geometry=geom)

def compute_mean_wse(data: gpd.GeoDataFrame, polygons: List[Polygon]):
    mean_wse = []
    data = data.reset_index()
    data = data.set_index(['utm_x','utm_y'])
    for poly in polygons:
        get_mean_height = lambda coord: tuple(data.loc[coord,['elevation']])
        toto = [get_mean_height(coord) for coord in poly.exterior.coords]
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
