
FPDEM_BASENAME = 'SWOT_L2_HR_FPDEM'
FPDEM_POINTCLOUD_BASENAME = 'SWOT_L2_HR_FPDEM_Ungridded'
FPDEM_RASTER_BASENAME = 'SWOT_L2_HR_FPDEM_Gridded'

MASK_SUFFIX = '_extract_area.shp'
POLYGON_SUFFIX = '_polygon.shp'

import os

def compute_name(output_directory, output_basename, tile_name, first_date_name, last_date_name):
    return os.path.join(output_directory, output_basename+'_'+tile_name+'_'+first_date_name+'_'+last_date_name+"_PGA2"+"_01.nc")

