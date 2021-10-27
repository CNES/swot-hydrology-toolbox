# -*- coding: utf8 -*-
'''
Create raster from floodplain dem pixel cloud

Copyright (c) 2018, CNES
'''
import os
import logging
import numpy as np
from floodplain.geom.alpha_shape import alpha_shape_with_cgal
from floodplain.geom.tools import filter_polygons
from floodplain.utils.spatial import convert_polygons_utm_to_latlon
import utm
import floodplain.io.shp as shp
import argparse
import lib.my_rdf_file as my_rdf
import xarray as xr
from shapely.ops import unary_union


class Extract_Area(object):
    """
    Class Floodplain
    Main class to extract area of interest from set of tiles
    """
    
    def __init__(self, param, input_file = None, output_file = None):
        """
        Constructor: initialize variables
        
        :param in_params: input parameters to run the processor
        :type in_params: dictionary
        """
        
        if input_file == None:
            self.input_file = param.getValue("input_file")
        else:
            self.input_file = input_file 

        if input_file == None:
            self.output_file = param.getValue("output_file")
        else:
            self.output_file = output_file 
                    
        self.method = param.getValue("method_extract")
        
        if self.method == 'alpha_shape':
            self.alpha = int(param.getValue("alpha"))
        if self.method == 'intersection':
            self.smoothing_factor= float(param.getValue("smoothing_factor"))

    def compute_extract_area(self):
        if self.method == 'alpha_shape':
            self.compute_area_of_interest_alpha_shape()
        if self.method == 'intersection':
            self.compute_area_of_interest_intersection()

    def load_input_extract_area(self):
        if self.method == 'alpha_shape':
            cloud_xr = xr.open_dataset(self.input_file) 
            self.input_cloud_df = cloud_xr.to_dataframe()
        if self.method == 'intersection':
            self.input_polygons, self.wse = shp.from_file(self.input_file, wse_flag = True)

    def compute_area_of_interest_alpha_shape(self):
 

        coords = np.zeros((self.input_cloud_df['x'].values.size, 2))
        coords[:, 0], coords[:, 1] = self.input_cloud_df['x'].values, self.input_cloud_df['y'].values
        
        polygons = alpha_shape_with_cgal(coords, self.alpha)
        
        filtered_polygons = filter_polygons(polygons)
        zone_number = utm.latlon_to_zone_number(self.input_cloud_df.iloc[0]['latitude'],
                                                        self.input_cloud_df.iloc[0]['longitude'])
        north = True if self.input_cloud_df.iloc[0]['latitude'] > 0 else False
        self.filtered_polygons_latlon = convert_polygons_utm_to_latlon(filtered_polygons,zone_number,north)
        
            
            
    def compute_area_of_interest_intersection(self):
        
        level = getattr(logging, "INFO")
        logging.basicConfig(filename=None,format='%(asctime)s [%(levelname)s] %(message)s', level=level)
        logging.info("Compute area between min and max water level from intersection method")     
        
        smoothing_factor = self.smoothing_factor
        polygons_dilated = []
        for i in self.input_polygons:
            polygons_dilated.append(i.buffer(smoothing_factor))
        
        logging.info("Compute exterior polygon (max water level)")  
        # Get exterior polygon (max water level)
        polygons_exterior = unary_union(polygons_dilated)
        polygons_interior = []
        
        if not isinstance(polygons_exterior, list): 
            polygons_exterior=[polygons_exterior]
        
        logging.info("Compute interior polygon (max water level)") 
        list_poly = []
        for poly in self.input_polygons:
            list_poly.append(poly)
            
        # Get interior polygon (min water level)
        
            
        # old method sort by area
        
        list_poly.sort(key = get_area)
        good_poly = []
        for poly in list_poly:
            if not good_poly:
                good_poly.append(poly)
            else:
                condition = True
                for little_poly in good_poly:
                    if poly.intersection(little_poly).area > 0.001*poly.area:
                        condition = False
                if condition:
                    good_poly.append(poly)
        for poly in good_poly:
            polygons_interior.append(poly)


            # Sort by water level   

        # ~ poly_wse_list = list(zip(self.input_polygons, self.wse))
        # ~ poly_wse_list.sort(key = takeSecond)
        # ~ list_poly, list_wse = list(zip(*poly_wse_list))
     
        # ~ good_poly = []
        # ~ good_wse = []
        # ~ for i, poly in enumerate(list_poly):
            # ~ if not good_poly:
                # ~ good_poly.append(poly)
                # ~ good_wse.append(list_wse[i])
            # ~ else:
                # ~ condition = True
                # ~ for little_poly in good_poly:
                    # ~ if poly.intersection(little_poly).area > 0.:
                        # ~ condition = False
                # ~ if condition:
                    # ~ good_poly.append(poly)
                    # ~ good_wse.append(list_wse[i])
        
        # ~ to_remove=[]
        # ~ for i, poly_i in enumerate(good_poly):
            # ~ continue_flag = True
            # ~ for k, poly_k in enumerate(list_poly): 
                # ~ if good_wse[i] - list_wse[k]> 0. and continue_flag == True:
                    # ~ if poly_i.intersection(poly_k).area > 0.: #0.01*poly_i.area:
                        # ~ to_remove.append(i)
                        # ~ continue_flag = False
                                
        # ~ wse_interior=[]
        # ~ for i, poly in enumerate(good_poly):
            # ~ if i not in to_remove:
                # ~ polygons_interior.append(poly)
                # ~ wse_interior.append(good_wse[i])

        
        polygons_interior_cleaned = []
        for i in polygons_interior:
            poly_cleaned = i.buffer(smoothing_factor)   
            polygons_interior_cleaned.append(poly_cleaned) 
        polygons_interior_cleaned = [unary_union(polygons_interior_cleaned).buffer(-smoothing_factor)]   
                 
        # Get area between and max water level
        polygons_between = []
        for poly_ext in polygons_exterior:
            poly_difference = poly_ext
            for poly_int in polygons_interior_cleaned:
                poly_difference = poly_difference.difference(poly_int)
            polygons_between.append(poly_difference)
        
        # ~ #Cleaning (erode + dissolve)
        self.polygons_between_cleaned = polygons_between
        

            
        # ~ shp.polygons_to_file("/work/ALT/swot/swotdev/desrochesd/floodplain/run/biglake150/output_dask/output_extract_area_interior.shp", list_poly, wse = list_wse)
        # ~ shp.polygons_to_file("/work/ALT/swot/swotdev/desrochesd/floodplain/run/biglake150/output_dask/output_extract_area_interior2.shp", self.input_polygons[0::50], wse = self.wse[0::50])
        # ~ shp.polygons_to_file("/work/ALT/swot/swotdev/desrochesd/floodplain/run/LacOrient/output/output_extract_area_exterior.shp", polygons_exterior)
        

    def write_extract_area(self):
        if self.method == 'alpha_shape':
            shp.polygons_to_file(self.output_file, self.filtered_polygons_latlon)
        if self.method == 'intersection':
            shp.polygons_to_file(self.output_file, self.polygons_between_cleaned)


def takeSecond(elem):
    return elem[1]
    
def get_area(elem):
    return elem.area

# Main program
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Compute extent mask between min and max water level")
    parser.add_argument("parameter_file", help="parameter_file (*.rdf)")
    args = parser.parse_args()
    parameters = my_rdf.myRdfReader(args.parameter_file)

    if parameters.getValue("method_extract") == 'alpha_shape':
        extract_area = Extract_Area(parameters)
    if parameters.getValue("method_extract") == 'intersection':
        extract_area = Extract_Area(parameters)
        
    extract_area.load_input_extract_area()
    extract_area.compute_extract_area()
    extract_area.write_extract_area()
