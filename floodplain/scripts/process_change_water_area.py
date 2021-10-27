# -*- coding: utf8 -*-
'''
Create shapefile expanding water elevation and surface proportionnaly

Copyright (c) 2021, CNES
'''

import argparse

import lib.my_rdf_file as my_rdf
import floodplain.io.shp as shp

class Expand_Water(object):
    def __init__(self, param, input_file = None, output_file = None):
        
        if input_file == None:
            self.input_file = param.getValue("input_file")
        else:
            self.input_file = input_file

        if output_file == None:
            self.output_file = param.getValue("output_filename")
        else:
            self.output_file = output_file

        self.expansion_factor = param.getValueList("water_expansion_factor")
        self.water_change = param.getValueList("height_change")
        self.height_attribute_name = param.getValue("height_attribute_name")
        print(self.expansion_factor, len(self.expansion_factor))

    def get_polygons(self):
        self.input_polygons, self.ref_height = shp.from_file(self.input_file, wse_flag=True, wse_name=self.height_attribute_name)

    def expand_water_area(self):
        for i in range(len(self.expansion_factor)):
            self.new_polygons=[]
            self.wse=[]
            for polygon, height in zip(self.input_polygons, self.ref_height):
                self.new_polygons.append(polygon.buffer(float(self.expansion_factor[i])))
                self.wse.append(float(height) + float(self.water_change[i]))
            shp.polygons_to_file('_'.join([self.output_file, str(float(height)+float(self.water_change[i])), 'm'])+'.shp', self.new_polygons, wse=self.wse)
    

# Main program
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Expand water mask and change water_height")
    parser.add_argument("parameter_file", help="parameter_file (*.rdf)")
    args = parser.parse_args()

    parameters = my_rdf.myRdfReader(args.parameter_file)

    my_new_water = Expand_Water(parameters)

    my_new_water.get_polygons()
    my_new_water.expand_water_area()
