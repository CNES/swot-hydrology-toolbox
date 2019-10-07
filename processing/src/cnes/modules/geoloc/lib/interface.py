# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2019/05/17:version initiale.
# FIN-HISTORIQUE
# ======================================================
'''
.. module:: interface.py

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''

import numpy as np
import os.path
import collections
import fiona
import copy
import logging
import cnes.common.lib.my_netcdf_file as my_nc
import cnes.common.service_error as service_error


class PixelCloud(object):
    """
        class PixelCloud
    """

    @classmethod
    def from_file(cls, in_pixc_file):
        """
        Constructor of pixel cloud
            :param in_pixc_file: full path of L2_HR_PIXC file
            :type in_pixc_file: string

        Variables of the object:
        - illumination_time / 1D-array of float:Time of illumination of each pixel
        - azimuth_index / 1D-array of int: azimuth indices of water pixels
        - range_index / 1D-array of int: range indices of water pixels
        - num_rare_looks / 1D-array of float Number of rare looks
        - wavelength / wavelength of the satellite measuring instrument
        - near_range / 1D-array of float: near range distance for each nadir point
        - range_spacing / float: nominal slant range spacing
        - azimuth_spacing / float: nominal slant azimuth spacing
        - tile_ref / list of string : list of tile references to process
        - nr_lines / int : number of tines to process
        - nb_pixels / int : number of pixels to process
        """

        self = cls()
        # 1 - Initiate logging service
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("PixelCloud init with %s file", in_pixc_file)

        try:
            # 2 - Open pixel cloud file in reading mode
            pixc_reader = my_nc.MyNcReader(in_pixc_file)
            data_dict = pixc_reader.content.groups['pixel_cloud']
            attr_id = pixc_reader.content.groups['pixel_cloud']

            # 3 - Initialization of the pixel cloud variables
            self.illumination_time = np.array(data_dict["illumination_time"])
            self.azimuth_index = np.array(data_dict["azimuth_index"])
            self.range_index = np.array(data_dict["range_index"])
            self.num_rare_looks = np.array(data_dict["num_rare_looks"])
            self.wavelength = pixc_reader.get_att_value("wavelength")
            self.near_range = pixc_reader.get_att_value("near_range")
            self.range_spacing = pixc_reader.get_att_value("nominal_slant_range_spacing")
            self.azimuth_spacing = 21.875
            self.tile_ref = pixc_reader.get_att_value("tile_name")
            self.nr_lines = int(attr_id.interferogram_size_azimuth)
            self.nr_pixels = int(attr_id.interferogram_size_range)

            # 4 - Close pixel cloud file
            pixc_reader.close()

        except Exception as exc:
            logger.debug(exc)
            raise

        return self

    @classmethod
    def from_dict(cls, pixc):
        """
        Constructor of pixel cloud
            :param pixc: Pixel dictionary
            :type pixc: dictionary

        Variables of the object:
        - illumination_time / 1D-array of float:Time of illumination of each pixel
        """
        self = cls()
        # 1 - Initiate logging service
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("PixelCloud init with a dictionary")

        try:
            # 2 - Initialization of the Pixel cloud variables
            self.sensor_s = pixc['illumination_time']
        except Exception as exc:
            logger.debug(exc)
            raise

        return self


class Sensor(object):
    """
        class Sensor
    """
    @classmethod
    def from_pixc(cls, pixc_file):
        """
        Constructor of Sensor
            :param pixc_file: full path of Sensor file
            :type pixc_file: string

        Variables of the object:
        - nadir_[x|y|z] / 1D-array of float: [x|y|z] cartesian coordinates of each nadir pixel (= variables named [x|y|z] in L2_HR_PIXC file)
        - nadir_[vx|vy|vz] / 1D-array of float: velocity vector of each nadir pixel in cartesian coordinates (= variables named velocity_unit_[x|y|z]
        in L2_HR_PIXC file)
        - ref_leverarm_[x|y|z] / 1D-array of float: [x|y|z]: reference of level arm 
        - sef_leverarm_[x|y|z] / 1D-array of float: [x|y|z]: level arm value
        """

        self = cls()
        # 1 - Initiate logging service
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Sensor init with %s pixc_file", pixc_file)

        try:
            # 2 - Open sensor file in reading mode
            sensor_reader = my_nc.MyNcReader(pixc_file)
            ifp = sensor_reader.content

            # 3 - Initialization of the Sensor variables
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

            # 4 - Close Sensor file
            sensor_reader.close()

        except Exception as exc:
            message = str(exc.__class__) + str(exc)
            logger.debug(message)
            raise

        return self

    @classmethod
    def from_file(cls, in_pixc_sensor_file):
        """
            :param in_pixc_sensor_file: full path of Sensor file
            :type in_pixc_sensor_file: string
        """

        self = cls()
        # 1 - Initiate logging service
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Sensor init with %s file sensor", in_pixc_sensor_file)
        try:
            # 2 - Open Sensor file in reading mode
            pixc_reader = my_nc.MyNcReader(in_pixc_sensor_file)
            pixc_sensor = pixc_reader.content

            # 3 - Initialization of the Sensor variables
            self.time = np.array(pixc_sensor.variables["time"])
            self.nadir_x = np.array(pixc_sensor.variables["x"])
            self.nadir_y = np.array(pixc_sensor.variables["y"])
            self.nadir_z = np.array(pixc_sensor.variables["z"])
            self.nadir_vx = np.array(pixc_sensor.variables["vx"])
            self.nadir_vy = np.array(pixc_sensor.variables["vy"])
            self.nadir_vz = np.array(pixc_sensor.variables["vz"])
            self.ref_leverarm_x = np.array(pixc_sensor.variables["plus_y_antenna_x"])
            self.ref_leverarm_y = np.array(pixc_sensor.variables["plus_y_antenna_y"])
            self.ref_leverarm_z = np.array(pixc_sensor.variables["plus_y_antenna_z"])
            self.sec_leverarm_x = np.array(pixc_sensor.variables["minus_y_antenna_x"])
            self.sec_leverarm_y = np.array(pixc_sensor.variables["minus_y_antenna_y"])
            self.sec_leverarm_z = np.array(pixc_sensor.variables["minus_y_antenna_z"])

            # 4 - Close Sensor file
            pixc_sensor.close()

        except Exception as exc:
            message = str(exc.__class__) + str(exc)
            logger.debug(message)
            raise

        return self

    @classmethod
    def from_dict(cls, sensor):
        """
        Constructor of Sensor
            :param sensor: Sensor dictionary
            :type sensor: dictionary
        """

        self = cls()
        # 1 - Initiate logging service
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("Sensor init with a dictionary")

        try:
            # 2 - Initialization of the Sensor variables
            self.time = sensor['time']
            self.nadir_x = sensor['x']
            self.nadir_y = sensor['y']
            self.nadir_z = sensor['z']
            self.nadir_vx = sensor['vx']
            self.nadir_vy = sensor['vy']
            self.nadir_vz = sensor['vz']
            self.sec_leverarm_x = sensor['sec_leverarm_x']
            self.sec_leverarm_y = sensor['sec_leverarm_y']
            self.sec_leverarm_z = sensor['sec_leverarm_z']
            self.ref_leverarm_x = sensor['ref_leverarm_x']
            self.ref_leverarm_y = sensor['ref_leverarm_y']
            self.ref_leverarm_z = sensor['ref_leverarm_z']

        except Exception as exc:
            message = str(exc.__class__) + str(exc)
            logger.debug(message)
            raise

        return self


class RiverTile(object):
    """
        class RiverTile
    """
    @classmethod
    def from_file(cls, in_river_main_file):
        """
        Constructor of RiverTile
            :param in_river_main_file: full path of L2_HR_PIXC file
            :type in_river_main_file: string

        Variables of the object:

        - From netcdf file (.nc)
            - [reach/node]_indx / [reach/node] 1D-array of float: tag associated to river node database, corresponding to the(reach index, node index)
             tuple from the pixel cloud
            - s 2D-array of float: curvilinear coordinate of distance along reach
            - h_n_ave / 1D-array of float:  node heights for each node previous node
            - h_n_ave_fit / 1D-array of float:  node heights for each node previous node

        - From netcdf file (.shp)
            - [reach/node]_indx / [reach/node] 1D-array of float: tag associated to river node database, corresponding to the(reach index, node index)
             tuple from the pixel cloud
            - s 2D-array of float: curvilinear coordinate of distance along reach
            - h_n_ave / 1D-array of float:  node heights for each node previous node

        """
        self = cls()
        root, ext = os.path.splitext(in_river_main_file)

        if ext == ".nc":
            self.read_nc(in_river_main_file)
        elif ext == ".shp":
            self.read_shp(in_river_main_file)
        else:
            # 1 - Initiate logging service
            logger = logging.getLogger(self.__class__.__name__)
            message = str(in_river_main_file) + " wrong extension expected *.shp or *.nc"
            logger.debug(message)
            raise service_error.FileError(message)
        return self

    @classmethod
    def from_node_outputs(cls, node_outputs):
        self = cls()
        # 1 - Initiate logging service
        logger = logging.getLogger(self.__class__.__name__)
        try:
            # 2 - Initialization of the RiverTile variables
            self.reach_indx = node_outputs['reach_indx']
            self.node_indx = node_outputs['node_indx']
            self.s = node_outputs['s']
            self.h_n_ave = node_outputs['h_n_ave']
            self.h_n_ave_fit = copy.deepcopy(node_outputs['h_n_ave'])
            self.fit_height = node_outputs['fit_height']
            # self.fit_height = []

        except Exception as exc:
            message = str(exc.__class__) + str(exc)
            logger.debug(message)
            raise

        return self

    def read_nc(self, in_river_main_file):
        # 1 - Initiate logging service
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("RiverTile init with %s netcdf file", in_river_main_file)

        # 2 - Open pixel RiverTile file in reading mode
        riverfile_reader = my_nc.MyNcReader(in_river_main_file)
        riverfile = riverfile_reader.content
        nodes = riverfile.groups["nodes"].variables
        try:
            # 2 - Initialization of the RiverTile variables
            self.reach_indx = np.array(nodes["reach_indx"])
            self.node_indx = np.array(nodes["node_indx"])
            self.s = np.array(nodes["s"])
            self.h_n_ave = np.array(nodes["h_n_ave"])
            self.h_n_ave_fit = copy.deepcopy(np.array(nodes["h_n_ave"]))

            # 3 - Close rivertile file
            riverfile_reader.close()

        except Exception as exc:
            message = str(exc.__class__) + str(exc)
            logger.debug(message)
            raise

    def read_shp(self, in_river_main_file):
        # 1 - Initiate logging service
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("RiverTile init with %s shape file", in_river_main_file)
        with fiona.open(in_river_main_file) as source:
            try:
                # 2 - Initialization of the RiverTile variables
                n = len(source)
                self.reach_indx = np.zeros(n)
                self.node_indx = np.zeros(n)
                self.s = np.zeros(n)
                self.h_n_ave = np.zeros(n)
                for i, f in enumerate(source):
                    self.reach_indx[i] = f['properties']['reach_indx']
                    self.node_indx[i] = f['properties']['node_indx']
                    self.s[i] = f['properties']['s']
                    self.h_n_ave[i] = f['properties']['h_n_ave']

                # 3 - Close rivertile file
                source.close()

            except Exception as exc:
                message = str(exc.__class__) + str(exc)
                logger.debug(message)
                raise


class PixcvecRiver(object):
    """
        class PixcvecRiver
    """
    def __init__(self, pixcvec_file):
        """
        Constructor of PixcvecRiver
          :param pixcvec_file: full_path_to_PIXCVecRiver_file
          :type pixcvec_file: string
        Variables of the object:
           - pixc_main_file / full_path_to_PIXCVecRiver_file
           - cycle_num / int: cycle number
           - pass_num / int: pass number
           - tile_ref / list of string : list of tile references to process ex: ['42N-R', '43N-R', '44N-R', '45N-R']
           - range_idx / 1D-array of int: range indices of water pixels
           - azimuth_idx / 1D-array of int: azimuth indices of water pixels
           - longitude / 1D-array of float: longitude of water pixels
           - latitude / 1D-array of float: latitude of water pixels
           - height / 1D-array of float: height of water pixels
           - node_index / 1D-array of float: tag associated to river node database
           - range_idx / 1D-array of int: range indices of water pixels
           - river tag / 1D-array of str: tag associated to river database
           - along_reach / 1D-array of float: along reach water for each pixel
           - cross_reach / 1D-array of float: cross reach water for each pixel
           - distance_to_node / 1D-array of float: distance for each pixel to node
           - range_spacing / int:nominal slant range spacing
           - near_range / 1D-array of float: near range distance for each nadir point
        """

        # 0 - Initiate logging service
        logger = logging.getLogger(self.__class__.__name__)
        logger.info(" PixcvecRiver init with %s file", pixcvec_file)
        try:
            # 1 - Open PixcvecRiver cloud file in reading mode
            pixc_reader = my_nc.MyNcReader(pixcvec_file)
            pixc_main = pixc_reader.content

            # 2 - Initialization of the PixcvecRiver variables
            self.pixc_main_file = pixcvec_file
            self.cycle_num = pixc_main.getncattr("cycle_number")
            self.pass_num = pixc_main.getncattr("pass_number")
            self.tile_ref = pixc_main.getncattr("tile_name")
            self.range_idx = np.array(pixc_main.variables["range_index"])
            self.azimuth_idx = np.array(pixc_main.variables["azimuth_index"])
            self.longitude = np.array(pixc_main.variables["longitude_vectorproc"])
            self.latitude = np.array(pixc_main.variables["latitude_vectorproc"])
            self.height = np.array(pixc_main.variables["height_vectorproc"])
            self.node_index = np.array(pixc_main.variables["node_index"])
            self.reach_index = np.array(pixc_main.variables["reach_index"])
            self.along_reach = np.array(pixc_main.variables["along_reach"])
            self.cross_reach = np.array(pixc_main.variables["cross_reach"])
            self.distance_to_node = np.array(pixc_main.variables["distance_to_node"])
            self.range_spacing = pixc_main.getncattr("nominal_slant_range_spacing")

            near_range = pixc_main.getncattr("near_range")
            if isinstance(near_range, collections.Iterable):
                print("Warning, near_range is a list.")
                self.near_range = 0
            else:
                self.near_range = near_range

            # 3 - Close PixcvecRiver file
            pixc_main.close()

        except Exception as exc:
            message = str(exc.__class__) + str(exc)
            logger.debug(message)
            raise



    def set_variable(self, varname, array):
        # 1 - Open PixcvecRiver cloud file in adding mode
        pixc_reader = my_nc.MyNcReader(self.pixc_main_file, mode='a')
        pixc_main = pixc_reader.content
        pixc_main.variables[varname][:] = array
        pixc_main.close()

    def set_index_file(self, lat_corr, lon_corr, height_corr, pixc_index):
        # 1 - Open PixcvecRiver cloud file in adding mode
        pixc_reader = my_nc.MyNcReader(self.pixc_main_file, mode='a')
        pixc_main = pixc_reader.content
        var = pixc_main.createVariable(
            'pixc_index', pixc_index.dtype, ('record',))
        var[:] = pixc_index[:]
        pixc_main.variables["latitude_vectorproc"][:] = lat_corr
        pixc_main.variables["longitude_vectorproc"][:] = lon_corr
        pixc_main.variables["height_vectorproc"][:] = height_corr
        pixc_main.close()

