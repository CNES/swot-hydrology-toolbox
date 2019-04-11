# -*- coding: utf8 -*-
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

import cnes.common.lib.my_api as my_api
import cnes.common.lib_lake.locnes_filenames as my_names
import cnes.common.lib_lake.locnes_variables as my_var
import cnes.common.lib_lake.proc_pixc_vec as proc_pixc_vec


class PixC_Vec_SP(object):

    def __init__(self, IN_lake_tile_pixcvec_file_list, IN_objPixcEdgeSP, IN_lake_sp_dir, IN_continent):
        """
        This class is designed to process L2_HR_PIXCVec product .
        All pixels involved in entities covering more than one tile are processed here.
        The geolocation of pixels is improved for those pixels and PIXCVec NetCDF file is updated.
        
        The initialization of a PixC_Vec_SP consists in:
         - set class attributes
         - copy input LakeTile_pixcvec files into PIXCVec files

        NP: this object is related to one swath
        
        :param IN_lake_tile_pixcvec_file_list: list of LakeTile_pixcvec files full path concerning current swath
        :type IN_lake_tile_pixcvec_file_list: list of string
        :param IN_objPixcEdgeSP: list of subset of PixC for edge objects current swath
        :type IN_objPixcEdgeSP: proc_pixc_sp.PixC_Edge_SP
        :param IN_lake_sp_dir: output LakeSP directory
        :type IN_lake_sp_dir: string
        :param IN_continent: continent covered by the tile (if global var CONTINENT_FILE exists)
        :type IN_continent: string


        Variables of the object:
        
            - From LakeTile_pixcvec file
                - longitude_vectorproc / 1D-array of float: improved longitude of water pixels (= variable named longitude_vectorproc in LakeTile_pixcvec file)
                - latitude_vectorproc / 1D-array of float: improved latitude of water pixels (= variable named latitude_vectorproc in LakeTile_pixcvec file)
                - height_vectorproc / 1D-array of float: improved height of water pixels (= variable named height_vectorproc in LakeTile_pixcvec file)
                - river_lake_other_tag / 1D-array of float: tag associated to river and lake databases (= variable named river_lake_other_tag in LakeTile_pixcvec file)

            - From process:
                - lake_tile_pixcvec_file_list / list of str : list of input LakeTile_pixcvec files
                - pixc_vec_file_list / list of str : list of output PIXCVec files
                - objPixcEdgeSP / PixC_Edge_SP.PixC_Edge_SP object : subset of PixC related to pixels of objects at top/bottom edge of a PixC tile (output of PGE_LakeTile)
                - nb_water_pix / int : number of pixels to process
        """
        my_api.printInfo("[PixelCloudVecEdge] == INIT ==")
        
        # 1 - Init variables
        # List of LakeTile_pixcvec files concerning current swath
        self.lake_tile_pixcvec_file_list = IN_lake_tile_pixcvec_file_list
        # List of output PIXCVec files
        self.pixc_vec_file_list = []
        # List of PixC_SP objects of current swath
        self.objPixcEdgeSP = IN_objPixcEdgeSP
        # Continent processed
        self.continent = IN_continent
        # Initialize PIXCVec variables to 0
        self.longitude_vectorproc = np.zeros(self.objPixcEdgeSP.nb_pixels)
        self.latitude_vectorproc = np.zeros(self.objPixcEdgeSP.nb_pixels)
        self.height_vectorproc = np.zeros(self.objPixcEdgeSP.nb_pixels)
        self.tag = np.empty(self.objPixcEdgeSP.nb_pixels, dtype=object)
        self.tag [:] = ""

        # Init a list of tiles ref processed in thios class
        self.tile_ref_list = []

        # 2 - List of output files computation
        for lake_tile_pixcvec_file in self.lake_tile_pixcvec_file_list:
            
            # 2.1 - Compute output PIXCVec file full path
            pixc_vec_file = my_names.computePixcvecFilename(lake_tile_pixcvec_file, IN_lake_sp_dir)
            
            # 2.2 - Remove if exists
            if os.path.exists(pixc_vec_file):
                os.remove(pixc_vec_file)
            
            # 2.3 - Copy LakeTile_pixcvec file full path to PIXCVec file list
            self.pixc_vec_file_list.append(pixc_vec_file)

            # 2.3 - Extact tile ref from PixC Vec file name
            tile_ref = my_names.getInfoFromFilename(lake_tile_pixcvec_file, "LakeTile")["tile_ref"]
            self.tile_ref_list.append(tile_ref)

    # ----------------------------------------

    def updatePixcVec(self, IN_write_to_shp=False):
        """
        This function updates PIXCVec netcdf files obtained in output of PGE_LakeTile with improved longitude, latitude, height and lake_tile_id or prior_id if exists.
        
        :param IN_write_to_shp: =True to also write the shapefile of the PIXCVec product (default=False)
        :type IN_write_to_shp: boolean
        """
        my_api.printInfo("[PixelCloudVecEdge] == updatePixcVec ==")

        for tile_idx in range(len(self.tile_ref_list)):  # Loop over tiles

            tile_ref = self.tile_ref_list[tile_idx]
            
            # 1 - Get input and output files info
            # 1.1 - LakeTile_pixcvec path and filename
            lake_tile_pixcvec_file = self.lake_tile_pixcvec_file_list[tile_idx]
            # 1.2 - Get PIXCVec associated path and filename
            pixc_vec_file = self.pixc_vec_file_list[tile_idx]

            # 2 - Init proc_pixc_vec.PixelCloudVec object
            objPixCVec = proc_pixc_vec.PixelCloudVec("SP", lake_tile_pixcvec_file)
            objPixCVec.setContinent(self.continent)

            # 3 - Get corresponding objPixcEdgeSP tile_idx
            objPixcEdgeSP_tile_idx = self.objPixcEdgeSP.tile_ref.index(str(tile_ref))

            # 4 - Update PIXCVec info
            if (self.objPixcEdgeSP.tile_idx == objPixcEdgeSP_tile_idx).any():  # Only when pixels need to be updated

                # 4.1 - Retrieve corresponding indices in PixC_SP object
                pixc_sp_idx = np.where(self.objPixcEdgeSP.tile_idx == objPixcEdgeSP_tile_idx)[0]
                
                # 4.2 - Retrieve corresponding indices in original PIXC
                pixc_tile_idx = self.objPixcEdgeSP.edge_idx[pixc_sp_idx]

                # 4.3 - Update geolocation information if computed
                if my_var.IMP_GEOLOC:
                    objPixCVec.longitude_vectorproc[pixc_tile_idx] = self.longitude_vectorproc[pixc_sp_idx]
                    objPixCVec.latitude_vectorproc[pixc_tile_idx] = self.latitude_vectorproc[pixc_sp_idx]
                    objPixCVec.height_vectorproc[pixc_tile_idx] = self.height_vectorproc[pixc_sp_idx]
                    
                # 4.4 - Update tag
                objPixCVec.tag[pixc_tile_idx] = str(self.tag[pixc_sp_idx])

            # 5 - Write PIXCVec file
            objPixCVec.write_file(pixc_vec_file, None)

            # 6 - Write associated shapefile if asked
            if IN_write_to_shp:
                objPixCVec.write_file_asShp(pixc_vec_file.replace('.nc', '.shp'))
            
            my_api.printInfo("")
