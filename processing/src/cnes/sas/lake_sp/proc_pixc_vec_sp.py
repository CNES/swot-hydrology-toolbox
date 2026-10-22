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
# VERSION:5.0.0:DM:#91:2022/12/08:Poursuite industrialisation
# VERSION:6.1.0:DM:#91:2024/09/20:Poursuite industrialisation
# VERSION:6.3.2:DM:#900:2026/07/06: [PGE] [LAKES] Ajout d'une interface niveau moyen des lacs
# VERSION:6.3.2:DM:#900:2026/07/06:SSALTO-565 [PGE] [LAKES] Ajout d'une interface coeff NN pour calcul ice_dyn_f
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

import cnes.common.service_config_file as service_config_file

import cnes.common.lib.my_tools as my_tools
import cnes.common.lib_lake.locnes_filenames as my_names
import cnes.common.lib_lake.proc_pixc_vec as proc_pixc_vec


class PixCVecSwath(object):
    """
    class PixCVecSwath
    Manage L2_HR_PIXCVec product over a single swath
    All pixels involved in water regions covering more than one tile are processed here.
    The geolocation of pixels is improved for those pixels and PIXCVec NetCDF file is updated.
    
    !!! IMPORTANT NOTE: If in_obj_pixc_edge_sp is None, it means that the error mode is activated: 
        PIXCVec products are just a copy of LakeTile PIXCVec files !!!
    """
    def __init__(self, in_list_laketile_pixcvec, in_obj_pixc_edge_sp, in_continent_id):
        """
        Constructor
        
        :param in_list_laketile_pixcvec: list of LakeTile_PIXCVec files to process
        :type in_list_laketile_pixcvec: list
        :param in_obj_pixc_edge_sp: subset of PixC for along-track edge objects; is None when the error mode is activated
        :type in_obj_pixc_edge_sp: proc_pixc_sp.PixCEdgeSP
        :param in_continent_id: 2-letter identifier of the continent covered by the swath (if global var CONTINENT_FILE exists)
        :type in_continent_id: string

        Variables of the object:
        
            - All variables from LakeTile_PIXCVec file 
            - From processing:
                - list_laketile_pixcvec / list of str: list of input LakeTile_PIXCVec files
                - obj_pixc_edge_sp / proc_pixc_sp.PixCEdgeSP: subset of PixC related to pixels of objects at top/bottom edge of a PixC tile \
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
        
        # List of PixC_SP objects of current swath or None in case of activation of error mode
        self.obj_pixc_edge_sp = in_obj_pixc_edge_sp
        
        # Continent identifier processed
        self.continent_id = in_continent_id
        
        if self.obj_pixc_edge_sp is not None:

            # Initialize PIXCVec variables to 0 or ""
            self.longitude_vectorproc = np.zeros(self.obj_pixc_edge_sp.nb_pixels, dtype=np.float64)
            self.latitude_vectorproc = np.zeros(self.obj_pixc_edge_sp.nb_pixels, dtype=np.float64)
            self.height_vectorproc = np.zeros(self.obj_pixc_edge_sp.nb_pixels, dtype=np.float32)
            
            self.delta_mean_lake_surface = np.zeros(self.obj_pixc_edge_sp.nb_pixels, dtype=np.float64)
            
            self.used_for_lake_wse = np.zeros(self.obj_pixc_edge_sp.nb_pixels, dtype=np.int8)
            self.used_for_lake_wse2 = np.zeros(self.obj_pixc_edge_sp.nb_pixels, dtype=np.int8)
            
            self.reach_id = np.chararray(self.obj_pixc_edge_sp.nb_pixels, itemsize=11)
            self.reach_id[:] = ""
            self.node_id = np.chararray(self.obj_pixc_edge_sp.nb_pixels, itemsize=14)
            self.node_id[:] = ""
            self.lake_id = np.chararray(self.obj_pixc_edge_sp.nb_pixels, itemsize=10)
            self.lake_id[:] = ""
            self.obs_id = np.chararray(self.obj_pixc_edge_sp.nb_pixels, itemsize=13)
            self.obs_id[:] = ""
            
            self.ice_clim_f = np.zeros(self.obj_pixc_edge_sp.nb_pixels, dtype=np.int8)
            self.ice_dyn_f = np.zeros(self.obj_pixc_edge_sp.nb_pixels, dtype=np.int8)

    # ----------------------------------------

    def update_pixcvec(self, in_output_dir, in_write_to_shp=False, in_proc_metadata=None, in_compress='zlib'):
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

        tmp_pixcvec_metadata = {} # WARNING: mandatory for the PGE

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

            # 5 - Remove if exists and if PIXCVec is updated in current SP processing
            if os.path.exists(pixcvec_file) and list_continent_id[0] == self.continent_id:
                logger.debug("Remove file %s" % pixcvec_file)
                os.remove(pixcvec_file)
                
            if self.obj_pixc_edge_sp is not None:

                # 6 - Get corresponding obj_pixc_edge_sp tile_idx
                pixc_sp_idx = np.where((self.obj_pixc_edge_sp.tile_index == tile_number) & \
                                       (self.obj_pixc_edge_sp.edge_loc < 10))[0]
    
                # 7 - Update PIXCVec info
                if pixc_sp_idx.any():  # Only when pixels need to be updated
                    logger.debug("Updating %d pixels of LakeTile_PIXCVec file" % pixc_sp_idx.size)
    
                    # 7.1 - Retrieve corresponding indices in original PIXC
                    pixc_tile_idx = self.obj_pixc_edge_sp.edge_index[pixc_sp_idx]
    
                    # 7.2 - Update geolocation information if computed
                    if self.cfg.getboolean('CONFIG_PARAMS', 'IMP_GEOLOC'):
                        obj_pixcvec.longitude_vectorproc[pixc_tile_idx] = self.longitude_vectorproc[pixc_sp_idx]
                        obj_pixcvec.latitude_vectorproc[pixc_tile_idx] = self.latitude_vectorproc[pixc_sp_idx]
                        obj_pixcvec.height_vectorproc[pixc_tile_idx] = self.height_vectorproc[pixc_sp_idx]
                        
                    # 7.3 - Update difference between Mean Lake Surface and EGM2008 if exists
                    obj_pixcvec.delta_mean_lake_surface[pixc_tile_idx] = self.delta_mean_lake_surface[pixc_sp_idx]
                    
                    # 7.4 - Update pixel assignation for wse and wse2 computation
                    obj_pixcvec.used_for_lake_wse[pixc_tile_idx] = self.used_for_lake_wse[pixc_sp_idx]
                    obj_pixcvec.used_for_lake_wse2[pixc_tile_idx] = self.used_for_lake_wse2[pixc_sp_idx]
    
                    # 7.5 - Update identifiers
                    obj_pixcvec.lake_id[pixc_tile_idx] = self.lake_id[pixc_sp_idx]
                    obj_pixcvec.obs_id[pixc_tile_idx] = self.obs_id[pixc_sp_idx]
                    
                    # 7.6 - Update ice flags
                    obj_pixcvec.ice_clim_f[pixc_tile_idx] = self.ice_clim_f[pixc_sp_idx]
                    obj_pixcvec.ice_dyn_f[pixc_tile_idx] = self.ice_dyn_f[pixc_sp_idx]
                    
                    # 7.7 - Update time_coverage_start and _end global attributes
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
                for key, value in in_proc_metadata.items():
                    obj_pixcvec.pixcvec_metadata[key] = value
            if self.obj_pixc_edge_sp is None:
                obj_pixcvec.pixcvec_metadata["is_complete"] = "False"
            else:
                obj_pixcvec.pixcvec_metadata["is_complete"] = "True"

            # 8 - Write PIXCVec file only if major continent of tile corresponds to current continent processed
            if list_continent_id[0] == self.continent_id:
                
                # 8.1 - Write PIXCVec file
                obj_pixcvec.write_file(pixcvec_file, None, in_compress=in_compress)

                # 8.2 - Write associated shapefile if asked
                if (self.obj_pixc_edge_sp is not None) and in_write_to_shp:
                    obj_pixcvec.write_file_as_shp(pixcvec_file.replace('.nc', '.shp'), self.obj_pixc_edge_sp)
            else :
                logger.debug("PIXCVec not generated because tile intersects more continent %s than current continent %s" \
                             % (list_continent_id[0], self.continent_id))

            # WARNING: mandatory for the PGE
            tmp_pixcvec_metadata[obj_pixcvec.pixcvec_metadata["tile_name"]] = obj_pixcvec.pixcvec_metadata

        return tmp_pixcvec_metadata # WARNING: mandatory for the PGE
