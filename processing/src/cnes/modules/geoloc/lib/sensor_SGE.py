import cnes.modules.geoloc.lib.my_netcdf_file as myCdf


class SensorSGE(object):
    
    def __init__(self, IN_pixc_main_file):
    
        # 1 - Retrieve needed information from pixel cloud main file
        sensor_main = myCdf.myNcReader(IN_pixc_main_file)
        # 1.1 - Number of pixels in range dimension
        self.altitude = sensor_main.getVarValue("altitude")
        # 1.1 - Number of pixels in range dimension
        self.azimuth_spacing = sensor_main.getAttValue("azimuth_spacing")
