# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2019/05/17:version initiale.
# VERSION:2.0.0:DM:#91:2020/07/03:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: proc_pixc_vec_sp.py
   :synopsis: Deals with SWOT pixel cloud complementary files (L2_HR_PIXCVec product) related to 1 single-pass and 1 swath (Left or Right)
    Created on 27/09/2017

.. moduleauthor:: Cécile Cazals - CS

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import absolute_import, division, print_function, unicode_literals 

import os
import numpy as np

import cnes.common.lib_lake.locnes_filenames as my_names
import cnes.common.lib_lake.proc_pixc_vec as proc_pixc_vec
import cnes.common.service_config_file as service_config_file

import logging

class PixCVecSP(object):
    """
    class PixCVecSP
    """
    def __init__(self, in_lake_tile_pixcvec_file_list, in_obj_pixc_edge_sp, in_lake_sp_dir, in_continent=None):
        """
        This class is designed to process L2_HR_PIXCVec product over 2 swath.

        :param in_lake_tile_pixcvec_file_list: list of LakeTile_pixcvec files full path concerning current swath
        :type in_lake_tile_pixcvec_file_list: list of string
        :param in_obj_pixc_edge_sp: list of subset of PixC for edge objects current swath
        :type in_obj_pixc_edge_sp: proc_pixc_sp.PixCEdgeSP
        :param in_lake_sp_dir: output LakeSP directory
        :type in_lake_sp_dir: string
        :param in_continent: continent covered by the tile (if global var CONTINENT_FILE exists)
        :type in_continent: string
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        lake_tile_pixcvec_file_path_r = [ file for file in in_lake_tile_pixcvec_file_list if "R_20" in file]
        lake_tile_pixcvec_file_path_r.sort()
        lake_tile_pixcvec_file_path_l = [ file for file in in_lake_tile_pixcvec_file_list if "L_20" in file]
        lake_tile_pixcvec_file_path_l.sort()

        self.pixcvec_r = PixCVecSwath(lake_tile_pixcvec_file_path_r, in_obj_pixc_edge_sp.pixc_edge_r, in_lake_sp_dir, in_continent)
        self.pixcvec_l = PixCVecSwath(lake_tile_pixcvec_file_path_l, in_obj_pixc_edge_sp.pixc_edge_l, in_lake_sp_dir, in_continent)


class PixCVecSwath(object):
    """
    class PixCVecSwath
    """
    def __init__(self, in_lake_tile_pixcvec_file_list, in_obj_pixc_edge_sp, in_lake_sp_dir, in_continent):
        """
        This class is designed to process L2_HR_PIXCVec product over a swath .
        All pixels involved in entities covering more than one tile are processed here.
        The geolocation of pixels is improved for those pixels and PIXCVec NetCDF file is updated.
        
        The initialization of a PixCVecSP consists in:
         - set class attributes
         - copy input LakeTile_pixcvec files into PIXCVec files

        NP: this object is related to one swath
        
        :param in_lake_tile_pixcvec_file_list: list of LakeTile_pixcvec files full path concerning current swath
        :type in_lake_tile_pixcvec_file_list: list of string
        :param in_obj_pixc_edge_sp: list of subset of PixC for edge objects current swath
        :type in_obj_pixc_edge_sp: proc_pixc_sp.PixCEdgeSP
        :param in_lake_sp_dir: output LakeSP directory
        :type in_lake_sp_dir: string
        :param in_continent: continent covered by the tile (if global var CONTINENT_FILE exists)
        :type in_continent: string


        Variables of the object:
        
            - From LakeTile_pixcvec file
                - longitude_vectorproc / 1D-array of float: improved longitude of water pixels (= variable named longitude_vectorproc in LakeTile_pixcvec file)
                - latitude_vectorproc / 1D-array of float: improved latitude of water pixels (= variable named latitude_vectorproc in LakeTile_pixcvec file)
                - height_vectorproc / 1D-array of float: improved height of water pixels (= variable named height_vectorproc in LakeTile_pixcvec file)
                - river_lake_other_tag / 1D-array of float: tag associated to river and lake databases (= variable named river_lake_other_tag in LakeTile_pixcvec file)

            - From process:
                - lake_tile_pixcvec_file_list / list of str : list of input LakeTile_pixcvec files
                - pixc_vec_file_list / list of str : list of output PIXCVec files
                - obj_pixc_edge_sp / proc_pixc_sp.PixCEdgeSP object : subset of PixC related to pixels of objects at top/bottom edge of a PixC tile (output of PGE_LakeTile)
                - nb_water_pix / int : number of pixels to process
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Init variables
        # Get instance of service config file
        self.cfg = service_config_file.get_instance()
        # List of LakeTile_pixcvec files concerning current swath
        in_lake_tile_pixcvec_file_list.sort()
        self.lake_tile_pixcvec_file_list = in_lake_tile_pixcvec_file_list
        # List of output PIXCVec files
        self.pixc_vec_file_list = []
        # List of PixC_SP objects of current swath
        self.obj_pixc_edge_sp = in_obj_pixc_edge_sp
        # Continent processed
        self.continent = in_continent

        # Initialize PIXCVec variables to 0 or ""
        self.longitude_vectorproc = np.zeros(self.obj_pixc_edge_sp.nb_pixels, dtype=np.double)
        self.latitude_vectorproc = np.zeros(self.obj_pixc_edge_sp.nb_pixels, dtype=np.double)
        self.height_vectorproc = np.zeros(self.obj_pixc_edge_sp.nb_pixels, dtype=np.float)
        
        self.reach_id = np.chararray(self.obj_pixc_edge_sp.nb_pixels, itemsize=11)
        self.reach_id[:] = ""
        self.node_id = np.chararray(self.obj_pixc_edge_sp.nb_pixels, itemsize=14)
        self.node_id[:] = ""
        self.lake_id = np.chararray(self.obj_pixc_edge_sp.nb_pixels, itemsize=10)
        self.lake_id[:] = ""
        self.obs_id = np.chararray(self.obj_pixc_edge_sp.nb_pixels, itemsize=13)
        self.obs_id[:] = ""
        
        self.ice_clim_f = np.zeros(self.obj_pixc_edge_sp.nb_pixels, dtype=np.uint8)
        self.ice_dyn_f = np.zeros(self.obj_pixc_edge_sp.nb_pixels, dtype=np.uint8)

        # Init a list of tiles ref processed in thios class
        self.tile_number_list = []

        # 2 - List of output files computation
        for lake_tile_pixcvec_file in self.lake_tile_pixcvec_file_list:
            
            # 2.1 - Compute output PIXCVec file full path
            pixc_vec_file = my_names.compute_pixcvec_filename(lake_tile_pixcvec_file, in_lake_sp_dir)
            
            # 2.2 - Remove if exists
            if os.path.exists(pixc_vec_file):
                os.remove(pixc_vec_file)
            
            # 2.3 - Copy LakeTile_pixcvec file full path to PIXCVec file list
            self.pixc_vec_file_list.append(pixc_vec_file)

            # 2.3 - Extact tile number from PixC Vec file
            tile_number = int(my_names.get_info_from_filename(lake_tile_pixcvec_file, "LakeTile")["tile_ref"][:-1])
            self.tile_number_list.append(tile_number)

    # ----------------------------------------

    def update_pixc_vec(self, in_write_to_shp=False, in_proc_metadata=None):
        """
        This function updates PIXCVec netcdf files obtained in output of PGE_LakeTile with improved longitude, latitude, 
        height and lake_tile_id or prior_id if exists.
        
        :param in_write_to_shp: =True to also write the shapefile of the PIXCVec product (default=False)
        :type in_write_to_shp: boolean
        :param in_proc_metadata: processing metadata
        :type in_proc_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)

        for tile_idx in range(len(self.tile_number_list)):  # Loop over tiles

            tile_num = self.tile_number_list[tile_idx]

            # 1 - Get input and output files info
            # 1.1 - LakeTile_pixcvec path and filename
            lake_tile_pixcvec_file = self.lake_tile_pixcvec_file_list[tile_idx]
            logger.info("Updating file %s" %(lake_tile_pixcvec_file))

            # 1.2 - Get PIXCVec associated path and filename
            pixc_vec_file = self.pixc_vec_file_list[tile_idx]

            # 2 - Init proc_pixc_vec.PixelCloudVec object
            obj_pixc_vec = proc_pixc_vec.PixelCloudVec("SP")
            obj_pixc_vec.set_from_pixcvec_file(lake_tile_pixcvec_file)
            obj_pixc_vec.pixcvec_metadata["xref_l2_hr_lake_sp_config_parameter_file"] = in_proc_metadata["xref_l2_hr_lake_sp_config_parameter_file"]

            # obj_pixc_vec.setContinent(self.continent)

            # 3 - Get corresponding obj_pixc_edge_sp tile_idx
            pixc_sp_idx = np.where(self.obj_pixc_edge_sp.tile_index == tile_num)[0]

            # 4 - Update PIXCVec info
            if pixc_sp_idx.any():  # Only when pixels need to be updated
                logger.debug("Updating %d pixels of pixc vec" %(pixc_sp_idx.size))
                # 4.1 - Retrieve corresponding indices in PixC_SP object
                # pixc_sp_idx = np.where(self.obj_pixc_edge_sp.tile_index == obj_pixc_edge_sp_tile_idx)[0]

                # 4.2 - Retrieve corresponding indices in original PIXC
                pixc_tile_idx = self.obj_pixc_edge_sp.edge_index[pixc_sp_idx]

                # 4.3 - Update geolocation information if computed
                if self.cfg.getboolean('CONFIG_PARAMS', 'IMP_GEOLOC'):
                    obj_pixc_vec.longitude_vectorproc[pixc_tile_idx] = self.longitude_vectorproc[pixc_sp_idx]
                    obj_pixc_vec.latitude_vectorproc[pixc_tile_idx] = self.latitude_vectorproc[pixc_sp_idx]
                    obj_pixc_vec.height_vectorproc[pixc_tile_idx] = self.height_vectorproc[pixc_sp_idx]

                # 4.4 - Update identifiers
                obj_pixc_vec.lake_id[pixc_tile_idx] = self.lake_id[pixc_sp_idx]
                obj_pixc_vec.obs_id[pixc_tile_idx] = self.obs_id[pixc_sp_idx]
                
                # 4.5 - Update ice flags
                obj_pixc_vec.ice_clim_f[pixc_tile_idx] = self.ice_clim_f[pixc_sp_idx]
                obj_pixc_vec.ice_dyn_f[pixc_tile_idx] = self.ice_dyn_f[pixc_sp_idx]

            else :
                logger.debug("Updating 0 pixels of pixc vec")

            # 5 - Write PIXCVec file
            obj_pixc_vec.write_file(pixc_vec_file, None)

            # 6 - Write associated shapefile if asked
            if in_write_to_shp:
                obj_pixc_vec.write_file_as_shp(pixc_vec_file.replace('.nc', '.shp'), self.obj_pixc_edge_sp)
