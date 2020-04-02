from osgeo import ogr
import os

def open_shp(path_file, in_poly=None):

    # 1 - Open shapefile in read-only access
    shp_driver = ogr.GetDriverByName(str('ESRI Shapefile'))  # Shapefile driver
    shp_data_source = shp_driver.Open(path_file, 0)

    # 2 - Get the lake_layer
    layer = shp_data_source.GetLayer()
    layer_multi_lake = shp_data_source.GetLayer()

    # 3 - Select some lakes among BD using in_poly
    id_lake=[]
    code_data=False
    if in_poly is not None:
        layer.SetSpatialFilter(in_poly)

        for feature in layer:
            try:
                if feature.GetField("code") not in id_lake:
                    id_lake.append(feature.GetField("code"))
                    code_data=True
            except ValueError:
                if feature.GetField("id") not in id_lake:
                    id_lake.append(feature.GetField("id"))

        if len(id_lake) == 1:
            id_lake.append(id_lake[0])
        if code_data:
            layer_multi_lake.SetAttributeFilter("code IN {}".format(tuple(id_lake)))
        else:
            layer_multi_lake.SetAttributeFilter("id IN {}".format(tuple(id_lake)))

    # 4 - Create an output datasource in memory
    mem_driver = ogr.GetDriverByName('MEMORY')  # Memory driver
    data_source = mem_driver.CreateDataSource('memData')

    # 5 - Open the memory datasource with write access
    mem_driver.Open('memData', 1)

    # 6 - Copy the lake_layer to memory
    data_source.CopyLayer(layer_multi_lake, 'lake_db')

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
