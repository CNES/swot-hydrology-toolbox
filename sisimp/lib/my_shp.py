from osgeo import ogr
import os

def open_shp(path_file, in_poly=None):


    # 1 - Open shapefile in read-only access
    shp_driver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Shapefile driver
    shp_data_source = shp_driver.Open(path_file, 0)

    # 2 - Get the lake_layer
    layer = shp_data_source.GetLayer()

    # 3 - Select some lakes among BD using in_poly
    if in_poly is not None:
        layer.SetSpatialFilter(in_poly)

    # 4 - Create an output datasource in memory
    mem_driver = ogr.GetDriverByName('MEMORY')  # Memory driver
    data_source = mem_driver.CreateDataSource('memData')

    # 5 - Open the memory datasource with write access
    mem_driver.Open('memData', 1)

    # 6 - Copy the lake_layer to memory
    data_source.CopyLayer(layer, 'lake_db')

    # 7 - Get memory lake_layer
    out_layer = data_source.GetLayer()
    out_layer.ResetReading()

    layer.SetSpatialFilter(None)

    # 8 - Close shapefile
    shp_data_source.Destroy()

    return out_layer, data_source


def overwrite_shapefile(ds_name, ds_format, geom_type, srs, overwrite=False):
    drv = ogr.GetDriverByName(ds_format)
    if os.path.exists(ds_name) and overwrite is True:
        drv.DeleteDataSource(ds_name)
    ds = drv.CreateDataSource(ds_name)
    lyr_name = os.path.splitext(os.path.basename(ds_name))[0]
    lyr = ds.CreateLayer(lyr_name, srs, geom_type)
    return ds, lyr