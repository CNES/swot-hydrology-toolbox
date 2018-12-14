#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
.. module:: sas_lake_tile.py
    :synopsis: Process PGE_L2_HR_LakeTile, i.e. generate L2_HR_LakeTile product from one tile of L2_HR_PIXC product and associated L2_HR_PIXCVec product
    Created on 02/27/2017

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

Copyright (c) 2017 CNES. All rights reserved.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import os
import logging

import cnes.sas.lake_tile.proc_pixc as proc_pixc
import cnes.common.lib.my_api as my_api
import cnes.common.lib.my_shp_file as my_shp
import cnes.common.lib.my_timer as my_timer
import cnes.common.lib.my_tools as my_tools
import cnes.common.lib_lake.lake_db as lake_db
import cnes.common.lib_lake.locnes_filenames as my_names
import cnes.common.lib_lake.locnes_variables as my_var
import cnes.common.lib_lake.proc_lake as proc_lake
import cnes.common.lib_lake.proc_pixc_vec as proc_pixc_vec


class SASLakeTile(object):
    """
        Class SASLakeTile
        SAS lake tile class
    """
    def __init__(self, IN_pixc_file, IN_pixc_vec_river_file, IN_output_dir, IN_shp_option=False):
        """
        Constructor: initialize variables

        :param IN_pixc_file: PIXC file full path
        :type IN_pixc_file: string
        :param IN_pixc_vec_river_file: PIXCVecRiver file full path
        :type IN_pixc_vec_river_file: string
        :param IN_output_dir: output directory full path
        :type IN_output_dir: string
        :param IN_shp_option: to also produce LakeTile_edge and LakeTile_pixcvec as shapefiles (=True); else=False (default)
        :type IN_shp_option: boolean
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("[lakeTileProcessing] == INIT ==")

        # Input paths
        self.pixc_file = IN_pixc_file  # PixC file
        self.pixc_vec_river_file = IN_pixc_vec_river_file  # Associated PIXCVecRiver file
        self.output_dir = IN_output_dir  # Output directory

        # Flag to produce LakeTile_edge and LakeTile_pixcvec shapefiles
        self.flag_prod_shp = False
        if IN_shp_option:
            self.flag_prod_shp = True

        # LakeTile filenames
        self.lake_tile_filenames = None

        # Objects
        self.objLakeDb = None  # Lake DB object
        self.objPixc = None  # PIXC object
        self.objPixcVec = None  # PIXCVecRiver object
        self.objLake = None  # LakeTile object

    def run_preprocessing(self):

        """
        Process PGE_LakeTile IN, i.e. test input paths, retrieve orbit infos, open lake database and init objects
        """
        logger = logging.getLogger(self.__class__.__name__)

        logger.info("")
        logger.info("")
        logger.info("[lakeTileProcessing] PRE-PROCESSING...")
        logger.info("")

        # 1 - Test existance and file format of input paths
        logger.info("[lakeTileProcessing] > 1 - Testing existence of input paths...")
        # 1.1 - PIXC file
        message = "[lakeTileProcessing]   INPUT PIXC file = %s" % self.pixc_file
        logger.info(message)
        my_tools.testFile(self.pixc_file, IN_extent=".nc")
        # 1.2 - PIXCVecRiver file
        message = "[lakeTileProcessing]   INPUT PIXCVecRiver file = %s" % self.pixc_vec_river_file
        logger.info(message)
        my_tools.testFile(self.pixc_vec_river_file, IN_extent=".nc")
        # 1.3 - Output directory
        message = "[lakeTileProcessing]   OUTPUT DIR = %s" % self.output_dir
        logger.info(message)
        my_tools.testDir(self.output_dir)
        logger.info("")

        # 2 - Retrieve orbit info from PIXC filename and compute output filenames
        logger.info("[lakeTileProcessing] > 2 - Retrieving tile infos from PIXC filename...")
        self.lake_tile_filenames = my_names.lakeTileFilenames(self.pixc_file, self.pixc_vec_river_file, self.output_dir)
        logger.info("")

        # 3 - Objects initialisation
        logger.info("[lakeTileProcessing] > 3 - Init and format intput objects...")
        logger.info("")

        # 3.1 - Init PIXCVec product by retrieving data from the pixel cloud complementary file after river processing
        logger.info("[lakeTileProcessing] > 3a - Init pixel cloud complementary file...")
        self.objPixcVec = proc_pixc_vec.PixelCloudVec("TILE", self.pixc_vec_river_file)
        logger.info("")

        # 3.2 - Retrieve needed data from the pixel cloud
        logger.info("[lakeTileProcessing] > 3b - Retrieving needed data from the pixel cloud...")
        self.objPixc = proc_pixc.PixelCloud(self.pixc_file, self.objPixcVec.reject_idx)
        logger.info("")

        # 3.3 - Reshape PIXCVec arrays
        logger.info("[lakeTileProcessing] > 3c - Reshape PIXCVecRiver arrays...")
        self.objPixcVec.reshape(self.objPixc)
        logger.info("")

        # 2 - Retrieve lake Db layer
        logger.info("[lakeTileProcessing] > 2 - Retrieving lake database layer...")
        if my_var.LAKE_DB == "":
            logger.info("[lakeTileProcessing] NO database specified -> NO link of SWOT obs with a priori lake")
        else:
            if os.path.exists(my_var.LAKE_DB):
                type_db = my_var.LAKE_DB.split('.')[-1]  # Type of database
                if type_db == "shp":  # Shapefile format
                    self.objLakeDb = lake_db.LakeDb_shp(my_var.LAKE_DB, self.objPixc.tile_poly)
                elif type_db == "sqlite":  # SQLite format
                    self.objLakeDb = lake_db.LakeDb_sqlite(my_var.LAKE_DB, self.objPixc.tile_poly)
                else:
                    my_api.exitWithError("[lakeTileProcessing] Lake a priori database format (%s) is unknown: must be .shp or .sqlite" % type_db)
            else:
                my_api.exitWithError("[lakeTileProcessing]   ERROR = %s doesn't exist" % my_var.LAKE_DB)
        logger.info("")

        # 4 - Initialize lake product
        logger.info("[lakeTileProcessing] > 4c - Reshape PIXCVecRiver arrays...")
        self.objLake = proc_lake.LakeProduct("TILE",
                                             self.objPixc,
                                             self.objPixcVec,
                                             self.objLakeDb,
                                             os.path.basename(self.lake_tile_filenames.lake_tile_shp_file).split(".")[0],
                                             IN_id_prefix=self.lake_tile_filenames.lake_id_prefix)
        logger.info("")


    def run_processing(self):
        """
        Process SAS_L2_HR_LakeTile
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("")
        logger.info("")
        logger.info("[lakeTileProcessing] PROCESSING...")
        logger.info("")

        timer_proc = my_timer.Timer()
        timer_proc.start()

        # Processing only if PixC pixels are selected
        if self.objPixc.nb_selected != 0:

            # 2 - F2-F3-F3b = Identify all separate entities in the water mask
            logger.info("[lakeTileProcessing] 1 - Identifying all separate entities in the water mask...")
            self.objPixc.computeSeparateEntities()
            logger.info("[lakeTileProcessing] " + timer_proc.info(0))
            logger.info("")

            # 3 - F4 = Retrieve pixels corresponding to lakes and new objects entirely inside the tile
            logger.info("[lakeTileProcessing] 2 - Getting pixels corresponding to lakes and new objects entirely inside the tile...")
            self.objPixc.computeObjInsideTile()
            logger.info("[lakeTileProcessing] " + timer_proc.info(0))
            logger.info("")

            # 4 - F6 = Fill lake product
            logger.info("[lakeTileProcessing] 3 - Filling LakeTile product...")
            self.objLake.computeLakeProducts(self.objPixc.labels_inside)
            logger.info("[lakeTileProcessing] " + timer_proc.info(0))
            logger.info("")

        else:
            logger.info("[lakeTileProcessing] NO selected PixC => empty lake tile product generated")
            logger.info("")

    def run_postprocessing(self):
        """
        Process PGE_L2_HR_LakeTile OUT, i.e. convert output data in L2_HR_LakeTile product and close files
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("")
        logger.info("")
        logger.info("[lakeTileProcessing] POST-PROCESSING...")
        logger.info("")

        # 1 - Write LakeTile shapefile
        logger.info("[lakeTileProcessing] 1 - Writing LakeTile memory layer to shapefile...")
        my_shp.write_mem_layer_as_shp(self.objLake.layer, self.lake_tile_filenames.lake_tile_shp_file)
        self.objLake.dataSource.Destroy()  # Close memory layer
        logger.info("")
        # Write XML metadatafile for shapefile
        self.objLake.writeMetadataFile("%s.xml" % self.lake_tile_filenames.lake_tile_shp_file)

        # 2 - Write PIXCVec for objects entirely inside tile
        logger.info("[lakeTileProcessing] 2 - Writing LakeTile_pixcvec file...")
        self.objPixcVec.write_file(self.lake_tile_filenames.lake_tile_pixcvec_file, None, True)
        if self.flag_prod_shp and (self.objPixcVec.nb_water_pix != 0):
            self.objPixcVec.write_file_asShp(self.lake_tile_filenames.lake_tile_pixcvec_file.replace(".nc", ".shp"), IN_classif=self.objPixc.origin_classif)
        logger.info("")

        # 3 - Write intermediate NetCDF file with indices of pixels (and associated label) related to objects at the top/bottom edges of the tile
        logger.info("[lakeTileProcessing] 3 - Writing LakeTile_edge file...")
        self.objPixc.write_edge_file(self.lake_tile_filenames.lake_tile_edge_file, None, True)
        if self.flag_prod_shp and (self.objPixc.nb_edge_pix != 0):
            self.objPixc.write_edge_file_asShp(self.lake_tile_filenames.lake_tile_edge_file.replace(".nc", ".shp"))
        logger.info("")

        # 4 - Close lake database
        if self.objLakeDb is not None:
            logger.info("[lakeTileProcessing] 4 - Closing lake database...")
            self.objLakeDb.close_db()
            logger.info("")


#######################################
