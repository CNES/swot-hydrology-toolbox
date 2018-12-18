'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


from netCDF4 import Dataset
import numpy as np
import os.path
import collections
import fiona
import copy


class PixelCloud(object):
    @classmethod
    def from_file(cls, IN_pixc_file):
        self = cls()
        pixc = Dataset(IN_pixc_file, 'r')

        data_dict = pixc.groups['pixel_cloud']
        attr_id = pixc.groups['pixel_cloud']
      
        #self.sensor_s = np.array(pixc.variables["sensor_s"])
        self.illumination_time = np.array(data_dict["illumination_time"])
        self.azimuth_index = np.array(data_dict["azimuth_index"])
        self.range_index = np.array(data_dict["range_index"])
        self.num_rare_looks = np.array(data_dict["num_rare_looks"])
        # get some of the attributes too
        self.wavelength = pixc.wavelength
        self.near_range = pixc.near_range
        self.range_spacing = pixc.nominal_slant_range_spacing
        self.azimuth_spacing = 21.875
        self.tile_ref = pixc.tile_number
        
        shape = attr_id.interferogram_shape
        self.nr_lines = int(shape.split(",")[0])
        self.nr_pixels = int((shape.split(",")[1]).split("(")[0])
                
        pixc.close()
        return self

    @classmethod
    def from_dict(cls, pixc):
        self = cls()
        #self.sensor_s = pixc['sensor_s']
        self.sensor_s = pixc['illumination_time']
        return self  
    
class Sensor(object):
    @classmethod
    def from_pixc(cls, pixc_file):
        self = cls()
        with Dataset(pixc_file, "r") as ifp:
            self.time = np.array(ifp.groups['tvp']['time'])
            self.nadir_x = np.array(ifp.groups['tvp']['x'])
            self.nadir_y = np.array(ifp.groups['tvp']['y'])
            self.nadir_z = np.array(ifp.groups['tvp']['z'])
            self.nadir_vx = np.array(ifp.groups['tvp']['vx'])
            self.nadir_vy = np.array(ifp.groups['tvp']['vy'])
            self.nadir_vz = np.array(ifp.groups['tvp']['vz'])
            self.ref_leverarm_x = np.array(ifp.groups['tvp']['plus_y_antenna_x'])
            self.ref_leverarm_y = np.array(ifp.groups['tvp']['plus_y_antenna_y'])
            self.ref_leverarm_z = np.array(ifp.groups['tvp']['plus_y_antenna_z'])
            self.sec_leverarm_x = np.array(ifp.groups['tvp']['minus_y_antenna_x'])
            self.sec_leverarm_y = np.array(ifp.groups['tvp']['minus_y_antenna_y'])
            self.sec_leverarm_z = np.array(ifp.groups['tvp']['minus_y_antenna_z'])
        return self

    @classmethod
    def from_file(cls, IN_pixc_sensor_file):
        self = cls()
        pixc_sensor = Dataset(IN_pixc_sensor_file, 'r')
        self.time = np.array(pixc_sensor.variables["time"])
        self.nadir_x = np.array(pixc_sensor.variables["x"])
        self.nadir_y = np.array(pixc_sensor.variables["y"])
        self.nadir_z = np.array(pixc_sensor.variables["z"])
        self.nadir_vx = np.array(pixc_sensor.variables["vx"])
        self.nadir_vy = np.array(pixc_sensor.variables["vy"])
        self.nadir_vz = np.array(pixc_sensor.variables["vz"])
        # get also the baseline info
        self.ref_leverarm_x = np.array(
            pixc_sensor.variables["plus_y_antenna_x"])
        self.ref_leverarm_y = np.array(
            pixc_sensor.variables["plus_y_antenna_y"])
        self.ref_leverarm_z = np.array(
            pixc_sensor.variables["plus_y_antenna_z"])
        self.sec_leverarm_x = np.array(
            pixc_sensor.variables["minus_y_antenna_x"])
        self.sec_leverarm_y = np.array(
            pixc_sensor.variables["minus_y_antenna_y"])
        self.sec_leverarm_z = np.array(
            pixc_sensor.variables["minus_y_antenna_z"])
            
        pixc_sensor.close()
        return self

    @classmethod
    def from_dict(cls, sensor):
        self = cls()
        self.time = sensor['time']
        self.nadir_x = sensor['x']
        self.nadir_y = sensor['y']
        self.nadir_z = sensor['z']
        self.nadir_vx = sensor['vx']
        self.nadir_vy = sensor['vy']
        self.nadir_vz = sensor['vz']
        # get also the baseline info
        self.sec_leverarm_x = sensor['sec_leverarm_x']
        self.sec_leverarm_y = sensor['sec_leverarm_y']
        self.sec_leverarm_z = sensor['sec_leverarm_z']
        self.ref_leverarm_x = sensor['ref_leverarm_x']
        self.ref_leverarm_y = sensor['ref_leverarm_y']
        self.ref_leverarm_z = sensor['ref_leverarm_z']
        return self

class RiverTile(object):
    @classmethod
    def from_file(cls, IN_river_main_file):
        self = cls()

        root, ext = os.path.splitext(IN_river_main_file)

        if ext == ".nc":
            self.read_nc(IN_river_main_file)
        elif ext == ".shp":
            self.read_shp(IN_river_main_file)

        return self

    @classmethod
    def from_node_outputs(cls, node_outputs):
        self = cls()
        self.reach_indx = node_outputs['reach_indx']
        self.node_indx = node_outputs['node_indx']
        self.s = node_outputs['s']
        self.h_n_ave = node_outputs['h_n_ave']
        self.h_n_ave_fit = copy.deepcopy(node_outputs['h_n_ave'])
        return self

    def read_nc(self, IN_river_main_file):
        # 1 - Retrieve needed information from pixel cloud main file
        riverfile = Dataset(IN_river_main_file, 'r')

        nodes = riverfile.groups["nodes"].variables

        self.reach_indx = np.array(nodes["reach_indx"])
        self.node_indx = np.array(nodes["node_indx"])
        self.s = np.array(nodes["s"])
        self.h_n_ave = np.array(nodes["h_n_ave"])
        self.h_n_ave_fit = copy.deepcopy(np.array(nodes["h_n_ave"]))

        riverfile.close()

    def read_shp(self, IN_river_main_file):
        with fiona.open(IN_river_main_file) as source:
            n = len(source)
            self.reach_indx = np.zeros((n))
            self.node_indx = np.zeros((n))
            self.s = np.zeros((n))
            self.h_n_ave = np.zeros((n))
            for i, f in enumerate(source):
                self.reach_indx[i] = f['properties']['reach_indx']
                self.node_indx[i] = f['properties']['node_indx']
                self.s[i] = f['properties']['s']
                self.h_n_ave[i] = f['properties']['h_n_ave']

class PixcvecRiver(object):
    def __init__(self, pixcvec_file):

        self.pixc_main_file = pixcvec_file
        # 1 - Retrieve needed information from pixel cloud main file
        pixc_main = Dataset(pixcvec_file, 'r')

        # 1.3 - Cycle number
        self.cycle_num = pixc_main.getncattr("cycle_number")
        # 1.4 - Pass number
        self.pass_num = pixc_main.getncattr("pass_number")
        # 1.5 - Tile reference
        self.tile_ref = pixc_main.getncattr("tile_name")
        
        # 1.7 - Range indices of water pixels
        self.range_idx = np.array(pixc_main.variables["range_index"])
        # 1.8 - Azimuth indices of water pixels
        self.azimuth_idx = np.array(pixc_main.variables["azimuth_index"])
        # this is to handle partial reading of the pixc as well 
        try:
            # 1.14 - Longitude
            self.longitude = np.array(pixc_main.variables["longitude_vectorproc"])
            # 1.15 - Latitude
            self.latitude = np.array(pixc_main.variables["latitude_vectorproc"])
            # 1.16 - Height
            self.height = np.array(pixc_main.variables["height_vectorproc"])
            # 1.16 - Node index
            self.node_index = np.array(pixc_main.variables["node_index"])
            # 1.16 - reach index
            self.reach_index = np.array(pixc_main.variables["reach_index"])
            # 1.16 - river tag
            self.river_tag = np.array(pixc_main.variables["river_tag"])
            # 1.16 - along_reach
            self.along_reach = np.array(pixc_main.variables["along_reach"])
            # 1.16 - cross_reach
            self.cross_reach = np.array(pixc_main.variables["cross_reach"])
            # 1.16 - distance_to_node
            self.distance_to_node = np.array(pixc_main.variables["distance_to_node"])

            #~ self.range_spacing = pixc_main.getncattr("nominal_slant_range_spacing")
            self.range_spacing = 0.74948114156723
            
            near_range = pixc_main.getncattr("near_range")
            if isinstance(near_range, collections.Iterable):
                print("Warning, near_range is a list.")
                self.near_range = 0
            else:
                self.near_range = near_range
        except KeyError:
            pass

        pixc_main.close()

    def add_new_variable(self, varname, datatype, array, size):

        pixc_main = Dataset(self.pixc_main_file, 'a')

        #pixc_main.createDimension('record2', 10)

        varlist = pixc_main.createVariable(varname, datatype, ('record'))
        varlist[:] = array

        pixc_main.close()

    def set_variable(self, varname, array):
        pixc_main = Dataset(self.pixc_main_file, 'a')
        pixc_main.variables[varname][:] = array
        pixc_main.close()

    def set_index_file(self, lat_corr, lon_corr, height_corr, pixc_index):
        pixc_main = Dataset(self.pixc_main_file, 'a')
        var = pixc_main.createVariable(
            'pixc_index', pixc_index.dtype, ('record',))
        var[:] = pixc_index[:]
        pixc_main.variables["latitude_vectorproc"][:] = lat_corr
        pixc_main.variables["longitude_vectorproc"][:] = lon_corr
        pixc_main.variables["height_vectorproc"][:] = height_corr
        pixc_main.close()
