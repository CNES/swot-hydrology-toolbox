# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:6.2.6:DM:#900:2026/07/06: [PGE] [LAKES] Ajout d'une interface niveau moyen des lacs
# VERSION:6.2.6:DM:#900:2026/07/06:SSALTO-565 [PGE] [LAKES] Ajout d'une interface coeff NN pour calcul ice_dyn_f
# VERSION:6.2.7:DM:#906:2026/07/28:SSALTO-872 [PGE][LAKES] Gestion des erreurs/warning GDAL
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: proc_mls.py
   :synopsis: Deals with Mean Lake Surfaces from SWOT Observations
    Created on 2026/06/11

.. moduleauthor:: Claire POTTIER - CNES DTN/TPI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import logging
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import LinearNDInterpolator
import tarfile

import cnes.common.service_config_file as service_config_file


class MeanLakeSurface(object):
    """
    class MeanLakeSurfaces
    Deal with Mean Lake Surfaces from SWOT Observations
    """
    
    def __init__(self):
        """
        Constructor: init mean lake surface object
        
        Variables of the object:
            - mls_data / dict: to store MLS; key1=lake_id ; key2=variable or attribute name ; value=value
        """
        # Get instance of service config file
        self.cfg = service_config_file.get_instance()
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")

        # Init MLS variables
        self.mls_data = dict()  # key1=lake_id ; key2=variable or attribute name ; value=value

    # ----------------------------------------
        
    def load_data(self, in_mls_tar, list_lakeid=None):
        """
        Load data in a dictionary
        
        :param in_mls_tar: full path of MLS TAR file
        :type in_mls_tar: string
        :param list_lakeid: list of needed lake_id
        :type list_lakeid: list of string
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("MLS TAR file = %s" % in_mls_tar)
        
        with tarfile.open(in_mls_tar, "r:*") as tar:
            
            for cur_fic_mls in tar.getmembers():
                
                # 1.1 - Retrieve lake_id
                lake_id = cur_fic_mls.name.split("_")[0]
                # 1.2 - Stop if lake_id is not needed
                if (list_lakeid is not None) and (lake_id not in list_lakeid):
                    continue
                # 1.3 - Create key for dictionary
                self.mls_data[lake_id] = dict()
                
                logger.debug("> Loading %s" % cur_fic_mls.name)
                
                # 2 - Load NetCDF file in memory
                mls_reader = tar.extractfile(cur_fic_mls)
                if mls_reader is None:
                    logger.warning("UNABLE to read %s" % cur_fic_mls.name)
                    continue
                
                # 3.1 - Extract NetCDF from tar in a bytes object
                raw_bytes = mls_reader.read()
                # 3.2 - Close NetCDF
                mls_reader.close()
                
                with Dataset(cur_fic_mls, mode="r", memory=raw_bytes) as nc:
                
                    # 4.1 - Extract and store needed variables
                    self.mls_data[lake_id]["longitude"] = nc.variables["longitude"][:]
                    self.mls_data[lake_id]["latitude"] = nc.variables["latitude"][:]
                    self.mls_data[lake_id]["height"] = nc.variables["height"][:]
                    # 4.2 - Extract and store needed global attributes
                    self.mls_data[lake_id]["geospatial_lat_min"] = getattr(nc, "geospatial_lat_min", None)
                    self.mls_data[lake_id]["geospatial_lat_max"] = getattr(nc, "geospatial_lat_max", None)
                    self.mls_data[lake_id]["geospatial_lon_min"] = getattr(nc, "geospatial_lon_min", None)
                    self.mls_data[lake_id]["geospatial_lon_max"] = getattr(nc, "geospatial_lon_max", None)
                    self.mls_data[lake_id]["spatial_resolution"] = getattr(nc, "spatial_resolution", None)

    # ----------------------------------------

    def interpolate(self, in_lakeid, in_longitude, in_latitude):
        """
        Interpolate MLS for PLD lake_id at (in_longitude, in_latitude) coordinates
        
        :param in_lakeid: identifier of the selected PLD lake
        :type in_lakeid: string
        :param in_longitude: longitude of the point where the interpolation of the MLS is performed
        :type in_longitude: numpy.array of float
        :param in_latitude: latitude of the point where the interpolation of the MLS is performed
        :type in_latitude: numpy.array of float
        
        :return: out_mls_interp = interpolated MLS; =None if there is no MLS for input the lake_id
        :rtype: out_mls_interp = numpy.array of float or NoneType
        """
        logger = logging.getLogger(self.__class__.__name__)
        
        if in_lakeid not in self.mls_data.keys():
            logger.debug("lake_id %s not in Mean Lake Surface database => DON'T correct the geoid" % in_lakeid)
            return None
        
        logger.debug("Interpolate MLS for lake_id %s" % in_lakeid)
            
        # 1.1 - Build the grid
        points = np.column_stack((self.mls_data[in_lakeid]["longitude"].ravel(), self.mls_data[in_lakeid]["latitude"].ravel()))
        values = self.mls_data[in_lakeid]["height"].ravel()
        # 1.2 - Keep only not NaN points
        if np.isnan(values).any():
            mask = ~np.isnan(values)
            points = points[mask]
            values = values[mask]
            
        # 2 - Make the interpolator
        interp = LinearNDInterpolator(points, values, fill_value=np.nan)

        # 3.1 - Make a table from the coordinates of points to interpolate
        target_points = np.column_stack((in_longitude, in_latitude))
        # 3.2 - Compute the interpolation
        out_mls_interp = interp(target_points)
        
        return out_mls_interp
            
