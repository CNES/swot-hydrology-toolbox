#!/usr/bin/env python
import sys
import fiona
import shapely.geometry
import rasterio.features
import xarray
import affine
import numpy as np
from collections import OrderedDict

def main():
    input_name = sys.argv[1]
    output_name = sys.argv[2]
    source_driver = 'ESRI Shapefile'
    source_crs = {'init': 'epsg:4326'}
    dest_schema = {'properties': OrderedDict([
        ('code', 'int:10'), ('basin_id', 'int:10'), ('name', 'str:254'),
        ('other_name', 'str:254'), ('reach_lake', 'str:254'),
        ('ref_height', 'str:254'), ('ref_date', 'int:10'),
        ('geoid_map', 'str:254'), ('ice_probab', 'str:254'),
        ('id', 'int:10')]), 'geometry': 'Polygon'}

    water_label = 1

    ds = xarray.open_dataset(input_name)
    water_mask = ds.variables['landtype']
    latitude = ds.variables['latitude']
    longitude = ds.variables['longitude']
    ds.close()

    water_pixel = np.where(water_mask.values == water_label)
    water_pixel_latitude, water_pixel_longitude = \
        latitude.values[water_pixel[0]], longitude.values[water_pixel[1]]

    dlon = longitude[1]-longitude[0]
    dlat = latitude[1]-latitude[0]

    tr = affine.Affine(
        dlon, 0, longitude[0]-0.5*dlon, 0, dlat, latitude[0]-0.5*dlat)

    shapes = rasterio.features.shapes(water_mask.values, transform=tr)
    polygons = [shapely.geometry.Polygon(
        shape[0]["coordinates"][0]) for shape in shapes if shape[1] == 1]

    with fiona.open(output_name,'w', driver=source_driver, crs=source_crs,
                    schema=dest_schema) as c:
        for p in polygons:
            prop = {'code': None, 'basin_id': None, 'name': None,
                    'other_name': None, 'reach_lake': None, 'ref_height': None,
                    'ref_date': None, 'geoid_map': None, 'ice_probab': None,
                    'id': None}
            c.write({'geometry': shapely.geometry.mapping(p), 
                     'properties': prop})

if __name__ == "__main__":
    main()
