import fiona
import shapely.geometry
import numpy as np
from netCDF4 import Dataset
import xarray as xr
from collections import OrderedDict
import rasterio.features
import shapely.geometry as geometry
from affine import Affine

#~ input_name="/work/ALT/swot/swotpub/SWOT_Simulator_data/input_data/garonne_aval/gdem-dem-truth/gdem-dem-garonne_aval-0.nc"
#~ input_name="/work/ALT/swot/swotpub/SWOT_Simulator_data/input_data/camargue/gdem-dem-truth/gdem-dem-camargue-0.nc"
input_name="/work/ALT/swot/swotpub/SWOT_Simulator_data/input_data/sacramento/gdem-dem-truth/gdem-dem-sacramento-1.nc"

#~ output_name="/work/ALT/swot/swotdev/desrochesd/swot-hydrology-toolbox/scripts/tools/camargue.shp"
output_name="/work/ALT/swot/swotdev/desrochesd/swot-hydrology-toolbox/scripts/tools/sacramento.shp"

source_driver = 'ESRI Shapefile'
source_crs = {'init': 'epsg:4326'}
dest_schema = {'properties': OrderedDict([('code', 'int:10'), ('basin_id', 'int:10'), ('name', 'str:254'), ('other_name', 'str:254'), ('reach_lake', 'str:254'), ('ref_height', 'str:254'), ('ref_date', 'int:10'), ('geoid_map', 'str:254'), ('ice_probab', 'str:254'), ('id', 'int:10')]), 'geometry': 'Polygon'}

water_label = 1

ds = xr.open_dataset(input_name)
water_mask = ds.variables['landtype']
latitude = ds.variables['latitude']         
longitude = ds.variables['longitude']          

water_pixel = np.where(water_mask.values == water_label)
water_pixel_latitude, water_pixel_longitude = latitude.values[water_pixel[0]], longitude.values[water_pixel[1]]

mypoly=[]

dlon = longitude[1]-longitude[0]
dlat = latitude[1]-latitude[0]

tr = Affine(dlon,0,longitude[0]-0.5*dlon,0,dlat,latitude[0]-0.5*dlat)

shapes = rasterio.features.shapes(water_mask.values, transform=tr)
polygons =[geometry.Polygon(shape[0]["coordinates"][0]) for shape in shapes if shape[1] == 1]
 
with fiona.open(output_name,'w', driver=source_driver, crs=source_crs, schema = dest_schema) as c:
    
    for p in polygons:
            #~ point = geometry.Point(float(water_pixel_longitude[i]), float(water_pixel_latitude[i]))
            prop = {'code':None, 'basin_id': None, 'name':None,'other_name':None, 'reach_lake':None,'ref_height':None,'ref_date':None,'geoid_map':None,'ice_probab':None,'id':None}
            c.write({'geometry': geometry.mapping(p), 'properties': prop})
