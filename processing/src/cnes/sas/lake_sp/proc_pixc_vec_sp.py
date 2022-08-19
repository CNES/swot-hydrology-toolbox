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
# VERSION:3.0.0:DM:#91:2021/03/12:Poursuite industrialisation
# VERSION:3.1.0:DM:#91:2021/05/21:Poursuite industrialisation
# VERSION:4.0.0:DM:#91:2022/05/05:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: proc_pixc_vec_sp.py
   :synopsis: Deals with SWOT pixel cloud complementary files (L2_HR_PIXCVec product) related to 1 single-pass and 1 swath (Left or Right)
    Created on 27/09/2017

.. moduleauthor:: Cécile Cazals - CS
                    Claire POTTIER - CNES DSO/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import absolute_import, division, print_function, unicode_literals 

import logging
import numpy as np
import os

import cnes.common.lib.my_tools as my_tools
import cnes.common.lib_lake.locnes_filenames as my_names
import cnes.common.lib_lake.proc_pixc_vec as proc_pixc_vec

import cnes.common.service_config_file as service_config_file


class PixCVecSwath(object):
    """
    class PixCVecSwath
    Manage L2_HR_PIXCVec product over a single swath
    All pixels involved in water regions covering more than one tile are processed here.
    The geolocation of pixels is improved for those pixels and PIXCVec NetCDF file is updated.
    """
    def __init__(self, in_list_laketile_pixcvec, in_obj_pixc_edge_sp, in_continent_id):
        """
        Constructor
        
        :param in_list_laketile_pixcvec: list of LakeTile_PIXCVec files to process
        :type in_list_laketile_pixcvec: list
        :param in_obj_pixc_edge_sp: subset of PixC for along-track edge objects
        :type in_obj_pixc_edge_sp: proc_pixc_sp.PixCEdgeSP
        :param in_continent_id: 2-letter identifier of the continent covered by the swath (if global var CONTINENT_FILE exists)
        :type in_continent_id: string

        Variables of the object:
        
            - All variables from LakeTile_PIXCVec file 
            + From process:
                - list_laketile_pixcvec / list of str: list of input LakeTile_PIXCVec files
                - obj_pixc_edge_sp / proc_pixc_sp.PixCEdgeSP: subset of PixC related to pixels of objects at top/bottom edge of a PixC tile
                  (output of PGE_LakeTile)
                - continent_id / string : 2-letter identifier of the continent covered by the swath
        """
        # Get instance of service config file
        self.cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # List of LakeTile_PIXCVec files concerning current swath
        in_list_laketile_pixcvec.sort()
        self.list_laketile_pixcvec = in_list_laketile_pixcvec
        # List of PixC_SP objects of current swath
        self.obj_pixc_edge_sp = in_obj_pixc_edge_sp
        # Continent identifier processed
        self.continent_id = in_continent_id

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

    # ----------------------------------------

    def update_pixcvec(self, in_output_dir, in_write_to_shp=False, in_proc_metadata=None, in_compress=True):
        """
        This function updates PIXCVec netcdf files obtained in output of PGE_LakeTile with improved longitude, latitude, 
        height and lake_tile_id or prior_id if exists.
        
        :param in_output_dir: full path of the output directory
        :type in_output_dir: string
        :param in_write_to_shp: =True to also write the shapefile of the PIXCVec product (default=False)
        :type in_write_to_shp: boolean
        :param in_proc_metadata: processing metadata
        :type in_proc_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        tmp_pixcvec_metadata = {}

        # Loop over LakeTile_PIXCVec input files
        for lake_tile_pixcvec_file in self.list_laketile_pixcvec:
            logger.debug("Updating file %s" % lake_tile_pixcvec_file)
            
            # 1 - Compute output PIXCVec file full path
            pixcvec_file = my_names.compute_pixcvec_filename(lake_tile_pixcvec_file, in_output_dir)
            
            # 2 - Extract tile number from PixC Vec file
            tile_number = int(my_names.get_info_from_filename(lake_tile_pixcvec_file, "LakeTile")["tile_ref"][:-1])

            # 3 - Init proc_pixc_vec.PixelCloudVec object
            obj_pixcvec = proc_pixc_vec.PixelCloudVec("SP")
            obj_pixcvec.set_from_pixcvec_file(lake_tile_pixcvec_file)

            # 4 - Get list of continent identifiers of current LakeTile_PIXCVec
            list_continent_id = obj_pixcvec.pixcvec_metadata["continent_id"].split(";")

            # 5 - Get corresponding obj_pixc_edge_sp tile_idx
            pixc_sp_idx = np.where(self.obj_pixc_edge_sp.tile_index == tile_number)[0]

            # 6 - Update PIXCVec info
            if pixc_sp_idx.any():  # Only when pixels need to be updated
                logger.debug("Updating %d pixels of LakeTile_PIXCVec file" % pixc_sp_idx.size)

                # 6.1 - Retrieve corresponding indices in original PIXC
                pixc_tile_idx = self.obj_pixc_edge_sp.edge_index[pixc_sp_idx]

                # 6.2 - Update geolocation information if computed
                if self.cfg.getboolean('CONFIG_PARAMS', 'IMP_GEOLOC'):
                    obj_pixcvec.longitude_vectorproc[pixc_tile_idx] = self.longitude_vectorproc[pixc_sp_idx]
                    obj_pixcvec.latitude_vectorproc[pixc_tile_idx] = self.latitude_vectorproc[pixc_sp_idx]
                    obj_pixcvec.height_vectorproc[pixc_tile_idx] = self.height_vectorproc[pixc_sp_idx]

                # 6.3 - Update identifiers
                obj_pixcvec.lake_id[pixc_tile_idx] = self.lake_id[pixc_sp_idx]
                obj_pixcvec.obs_id[pixc_tile_idx] = self.obs_id[pixc_sp_idx]
                
                # 6.4 - Update ice flags
                obj_pixcvec.ice_clim_f[pixc_tile_idx] = self.ice_clim_f[pixc_sp_idx]
                obj_pixcvec.ice_dyn_f[pixc_tile_idx] = self.ice_dyn_f[pixc_sp_idx]
                
                # 6.5 - Update time_coverage_start and _end global attributes
                tmp_nadir_time = self.obj_pixc_edge_sp.nadir_time[pixc_sp_idx]
                tmp_min = my_tools.convert_utc_to_str(min(tmp_nadir_time))
                if tmp_min < obj_pixcvec.pixcvec_metadata["time_coverage_start"]:
                    obj_pixcvec.pixcvec_metadata["time_coverage_start"] = tmp_min
                tmp_max = my_tools.convert_utc_to_str(max(tmp_nadir_time))
                if tmp_max < obj_pixcvec.pixcvec_metadata["time_coverage_end"]:
                    obj_pixcvec.pixcvec_metadata["time_coverage_end"] = tmp_max

            else:
                logger.debug("NO pixel of LakeTile_PIXCVec file need to be updated")

            if in_proc_metadata is not None:
                for key in in_proc_metadata.keys():
                    obj_pixcvec.pixcvec_metadata[key] = in_proc_metadata[key]

            # 7 - Write PIXCVec file only if major continent of tile corresponds to current continent processed
            if list_continent_id[0] == self.continent_id:
                
                # 7.1 - Write PIXCVec file
                obj_pixcvec.write_file(pixcvec_file, None, in_compress=in_compress)

                # 7.2 - Write associated shapefile if asked
                if in_write_to_shp:
                    obj_pixcvec.write_file_as_shp(pixcvec_file.replace('.nc', '.shp'), self.obj_pixc_edge_sp)

            tmp_pixcvec_metadata[obj_pixcvec.pixcvec_metadata["tile_name"]] = obj_pixcvec.pixcvec_metadata

        return tmp_pixcvec_metadata
