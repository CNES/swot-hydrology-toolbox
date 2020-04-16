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
"""
.. module:: locnes_products_netcdf.py
    :synopsis: Deals with SWOT NetCDF products
     Created on 02/18/2019

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2019 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
"""

from collections import OrderedDict
import datetime
import logging
import numpy as np

import cnes.common.lib.my_netcdf_file as my_nc
import cnes.common.lib.my_tools as my_tools


class NetCDFProduct(object):
    """
    Deal with some of SWOT NetCDF products: LakeTile_pixcvec, LakeTile_edge and PIXCVec
    """
    
    def __init__(self):
        """
        Constructor: set the general values
        
        Variables of the object:
        - dim / dict: key=dimension name ; value=number of elements in the 1-D array variables
        - is_complex / boolean: =True if there is at least 1 complex variable in the NetCDF (ie 2 dimensions); =False if not (default)
        - variables / OrderedDict: dictionary having key=variables names and value=OrderedDict of NetCDF variable attributes
        - metadata / OrderedDict: dictionary having key=global attributes and value=their values
        """
        
        # 1 - Init dim
        self.dim = {}
        self.dim["name"] = ""
        self.dim["value"] = 0
        self.is_complex = False
        
        # 2 - Init variables
        self.variables = OrderedDict()
        
        # 3 - Init metadata
        self.metadata = OrderedDict()
        self.metadata["conventions"] = "CF-1.7"
        self.metadata["title"] = ""
        self.metadata["short_name"] = ""
        self.metadata["institution"] = "CNES"
        self.metadata["source"] = ""
        self.metadata["history"] = "%sZ: Creation" % my_tools.swot_timeformat(datetime.datetime.utcnow(), in_format=1)
        self.metadata["mission_name"] = "SWOT"
        self.metadata["references"] = ""
        self.metadata["reference_document"] = ""
        self.metadata["contact"] = "claire.pottier@cnes.fr"
        
    #----------------------------------------
        
    def set_metadata_val(self, in_metadata):
        """
        Setter of metadata value
        
        :param in_metadata: metadata stored as a dictionary with key=metadata key and value=metadata value
        :type in_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        metadata_keys = self.metadata.keys()
        
        for key, value in in_metadata.items():
            if key in metadata_keys:
                self.metadata[key] = value
            else:
                logger.debug("%s key is not used as metadata" % key)

    #----------------------------------------
        
    def write_product(self, in_out_file, in_size, in_variables):
        """
        Write NetCDF product
        
        :param in_out_file: output file full path
        :type in_out_file: string
        :param in_size: number of values in 1-D arrays stored in in_variables
        :type in_size: int
        :param in_variables: dictionary with key=variable name and value=variable value
        :type in_param: OrderedDict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        # Open file in writing mode
        nc_writer = my_nc.MyNcWriter(in_out_file)
        
        # 1 - Write dimension(s)
        self.dim["value"] = in_size
        nc_writer.add_dimension(self.dim["name"], self.dim["value"])
        # Add specific dimension if NetCDF with complex variables
        if self.is_complex:
            nc_writer.add_dimension("complex_depth", 2)
        # Add specific dimensions for char variables
        for key, value in in_variables.items():
            if np.dtype(self.variables[key]['dtype']).char == 'S':
                dim_name = "nchar_" + key
                dim_size = int(self.variables[key]['dtype'][1:])
                nc_writer.add_dimension(dim_name, dim_size)
        
        # 2 - Write variables
        if in_size == 0:
            logger.info("Empty NetCDF file generated")
        else:
            self.write_variables(nc_writer, in_variables)
        
        # 3 - Write global attributes
        self.write_metadata(nc_writer)
        
        # Close file
        nc_writer.close()
        
    def write_variables(self, in_nc_writer, in_variables):
        """
        Write NetCDF variables listed in input
        
        :param in_nc_writer: NetCDF file writer
        :type in_nc_writer: my_netcdf_file.MyNcWriter
        :param in_variables: dictionary with key=variable name and value=variable value
        :type in_param: OrderedDict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        var_keys = self.variables.keys()
        
        for key, value in in_variables.items():
            if key in var_keys:
                # Add variable depending on its type
                if np.dtype(self.variables[key]['dtype']).char == 'S':  # Table of char
                    dim_nchar = "nchar_" + key
                    in_nc_writer.add_variable(key, 'c', (self.dim["name"], dim_nchar), in_attributes=self.variables[key])
                elif self.is_complex and (len(value.shape) == 2):  # Complex variable
                    in_nc_writer.add_variable(key, self.variables[key]['dtype'], (self.dim["name"], "complex_depth"), in_attributes=self.variables[key])
                else:
                    in_nc_writer.add_variable(key, self.variables[key]['dtype'], self.dim["name"], in_attributes=self.variables[key])
                # Fill variable
                in_nc_writer.fill_variable(key, value)
            else:
                logger.debug("Variable %s key is not known" % key)
    
    def write_metadata(self, in_nc_writer):
        """
        Write global attributes
        
        :param in_nc_writer: NetCDF file writer
        :type in_nc_writer: my_netcdf_file.MyNcWriter
        """          
        logger = logging.getLogger(self.__class__.__name__)
        logger.debug("- start -")
        
        for key, value in self.metadata.items():
            in_nc_writer.add_global_attribute(key, value)


#######################################


class CommonPixcvec(NetCDFProduct):
    """
    Gather items common to LakeTile_pixcvec and PIXCVec NetCDF files
    """
    
    def __init__(self):
        """
        Constructor
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Init NetCDF_product
        super().__init__()
        
        # 2 - Init dimension name
        self.dim["name"] = "points"
        
        # 3 - Init variables info; they are common to LakeTile_pixcvec and PIXCVec NetCDF files
        # All variables are defined with the same scheme :
        # type of variable
        # list of variable attribute name: value

        # azimuth_index
        self.variables["azimuth_index"] = {'dtype': np.int32,   
                                              'long_name': "rare interferogram azimuth index",
                                              'units': 1,
                                              'valid_min': 0,
                                              'valid_max': 999999,
                                              'comment': "Rare interferogram azimuth index"}
        # range_index
        self.variables["range_index"] = {'dtype': np.int32,
                                          'long_name': "rare interferogram range index",
                                          'units': 1,
                                          'valid_min': 0,
                                          'valid_max': 999999,
                                          'comment': "Rare interferogram range index"}
        # latitude_vectorproc
        self.variables["latitude_vectorproc"] = {'dtype': np.double,
                                                  'long_name': "height-constrained geolocation latitude",
                                                  'standard_name': "latitude",
                                                  'units': "degrees_north",
                                                  'valid_min': -80.0,
                                                  'valid_max': 80.0,
                                                  'comment': "Height-constrained geodetic latitude of the pixel. Units are in degrees north of the equator."}
        # longitude_vectorproc
        self.variables["longitude_vectorproc"] = {'dtype': np.double,
                                                  'long_name': "height-constrained geolocation longitude",
                                                  'standard_name': "longitude",
                                                  'units': "degrees_east",
                                                  'valid_min': -180.0,
                                                  'valid_max': 180.0,
                                                  'comment': "Height-constrained geodetic longitude of the pixel. Positive=degrees east of the Greenwich meridian. Negative=degrees west of the Greenwich meridian."}
        # height_vectorproc
        self.variables["height_vectorproc"] = {'dtype': np.float,
                                                  'long_name': "height above reference ellipsoid",
                                                  'units': "m",
                                                  'valid_min': -1500,
                                                  'valid_max': 1500,
                                                  'comment': "Height-constrained height of the pixel above the reference ellipsoid."}
        # reach_id
        self.variables["reach_id"] = {'dtype': 'S11',
                                      'long_name': "identifier of the associated prior river reach",
                                      'comment': "Unique reach identifier from the prior river database. \
                                                  The format of the identifier is CBBBBBRRRRT, where C=continent, B=basin, R=reach, T=type."}
        # node_id
        self.variables["node_id"] = {'dtype': 'S14',
                                      'long_name': "identifier of the associated prior river node",
                                      'comment': "Unique node identifier from the prior river database. \
                                                  The format of the identifier is CBBBBBRRRRNNNT, where C=continent, B=basin, R=reach, N=node, T=type of water body."}
        # lake_id
        self.variables["lake_id"] = {'dtype': 'S10',
                                      'long_name': "identifier of the associated prior lake",
                                      'comment': "Identifier of the lake from the lake prior database) associated to the pixel. \
                                                    The format of the identifier is CBBNNNNNNT, where C=continent, B=basin, N=counter within the basin, T=type of water body."}
        # obs_id
        self.variables["obs_id"] = {'dtype': 'S13',
                                      'long_name': "identifier of the observed object",
                                      'comment': "Tile-specific identifier of the observed object associated to the pixel. \
                                                      The format of the identifier is CBBTTTSNNNNNN, where C=continent, B=basin, T=tile number, S=swath side, N=lake counter within the PIXC tile."}
        # ice_clim_f
        self.variables["ice_clim_f"] = {'dtype': np.uint8,
                                          'long_name': "climatological ice cover flag",
                                          'flag_meanings': "no_ice_cover partial_ice_cover full_ice_cover",
                                          'flag_values': "0 1 2",
                                          'institution': "University of North Carolina",
                                          'comment': "Climatological ice cover flag indicating whether the pixel is ice-covered on the day of the observation \
                                                      based on external climatological information (not the SWOT measurement). \
                                                      Values of 0, 1, and 2 indicate that the surface is not ice covered, partially ice covered, \
                                                      and fully ice covered, respectively. A value of 255 indicates that this flag is not available."}
        # ice_dyn_f
        self.variables["ice_dyn_f"] = {'dtype': np.uint8,
                                          'long_name': "dynamical ice cover flag",
                                          'flag_meanings': "no_ice_cover partial_ice_cover full_ice_cover",
                                          'flag_values': "0 1 2",
                                          'institution': "University of North Carolina",
                                          'comment': "Dynamic ice cover flag indicating whether the pixel is ice-covered on the day of the observation \
                                                      based on analysis of external satellite optical data. Values of 0, 1, and 2 indicate \
                                                      that the surface is not ice covered, partially ice covered, and fully ice covered, respectively. \
                                                      A value of 255 indicates that this flag is not available."}
                
        # 4 - Init metadata common to LakeTile_pixcvec and PIXCVec files
        # 4.1 - Metadata retrieved from L2_HR_PIXC product
        # The values are defined later in TBD function
        self.metadata["cycle_number"] = -9999
        self.metadata["pass_number"] = -9999
        self.metadata["tile_number"] = -9999
        self.metadata["swath_side"] = ""
        self.metadata["tile_name"] = ""
        self.metadata["time_coverage_start"] = ""
        self.metadata["time_coverage_end"] = ""
        self.metadata["inner_first_latitude"] = -9999.0
        self.metadata["inner_first_longitude"] = -9999.0
        self.metadata["inner_last_latitude"] = -9999.0
        self.metadata["inner_last_longitude"] = -9999.0
        self.metadata["outer_first_latitude"] = -9999.0
        self.metadata["outer_first_longitude"] = -9999.0
        self.metadata["outer_last_latitude"] = -9999.0
        self.metadata["outer_last_longitude"] = -9999.0
        self.metadata["continent"] = ""
        self.metadata["ellipsoid_semi_major_axis"] = ""
        self.metadata["ellipsoid_flattening"] = ""
        # 4.2 - Processing metadata
        self.metadata["xref_input_l2_hr_pixc_file"] = ""
        self.metadata["xref_static_river_db_file"] = ""
        self.metadata["xref_static_lake_db_file"] = ""


class LakeTilePixcvecProduct(CommonPixcvec):
    """
    Deal with LakeTile_pixcvec NetCDF file
    """
    
    def __init__(self, in_pixcvecriver_metadata=None, in_proc_metadata=None):
        """
        Constructor; specific to LakeTile_pixcvec NetCDF file
        
        :param in_pixcvecriver_metadata: metadata specific to L2_HR_PIXCVecRiver product
        :type in_pixcvecriver_metadata: dict
        :param in_proc_metadata: metadata specific to LakeTile processing
        :type in_proc_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Init CommonPixcvec
        super().__init__()
                
        # 2 - Init metadata specific to LakeTile_pixcvec file
        # 2.1 - Update general metadata
        self.metadata["title"] = "Level 2 KaRIn high rate lake tile vector product"
        self.metadata["short_name"] = "SWOT_L2_HR_LakeTile_Pixcvec"
        # code commenté self.metadata["references"] = ""
        self.metadata["reference_document"] = "SWOT-TN-CDM-0677-CNES"
        if in_pixcvecriver_metadata is not None:
            self.set_metadata_val(in_pixcvecriver_metadata)
        # 2.3 - Processing metadata
        self.metadata["xref_input_l2_hr_pixc_vec_river_file"] = ""
        self.metadata["xref_l2_hr_lake_tile_param_file"] = ""
        if in_proc_metadata is not None:
            self.set_metadata_val(in_proc_metadata)


class PixcvecProduct(CommonPixcvec):
    """
    Deal with PIXCVec NetCDF file
    """
    
    def __init__(self, in_laketile_pixcvec_metadata=None, in_proc_metadata=None):
        """
        Constructor; specific to PIXCVec NetCDF file
        
        :param in_laketile_pixcvec_metadata: metadata specific to L2_HR_LakeTile_Pixcvec product
        :type in_laketile_pixcvec_metadata: dict
        :param in_proc_metadata: metadata specific to LakeTile processing
        :type in_proc_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Init CommonPixcvec
        super().__init__()
                
        # 4 - Init metadata specific to PIXCVec file
        # 4.1 - Update general metadata
        self.metadata["title"] = "Level 2 KaRIn high rate pixel cloud vector attribute product"
        self.metadata.pop("short_name", None)  # Metadata short_name has to be removed
        # code commenté self.metadata["references"] = ""
        self.metadata["reference_document"] = "SWOT-TN-CDM-0677-CNES"
        # 4.2 - Metadata retrieved from L2_HR_LakeTile_Pixcvec product
        if in_laketile_pixcvec_metadata is not None:
            self.set_metadata_val(in_laketile_pixcvec_metadata)
        # 4.3 - Processing metadata
        self.metadata["xref_l2_hr_lake_sp_param_file"] = ""
        if in_proc_metadata is not None:
            self.set_metadata_val(in_proc_metadata)


#######################################


class LakeTileEdgeProduct(NetCDFProduct):
    """
    Deal with LakeTile_edge NetCDF file
    """
    
    def __init__(self, in_pixc_metadata=None, in_proc_metadata=None):
        """
        Constructor; specific to LakeTile_edge NetCDF file
        
        :param in_pixc_metadata: metadata specific to L2_HR_PIXC product
        :type in_pixc_metadata: dict
        :param in_proc_metadata: metadata specific to LakeTile processing
        :type in_proc_metadata: dict
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.info("- start -")
        
        # 1 - Init NetCDF_product
        super().__init__()
        
        # 2 - Init dimension name
        self.dim["name"] = "points"
        self.is_complex = True
        
        # 3 - Init variables info specific to LakeTile_edge file
        # All variables are defined with the same scheme :
        # type of variable
        # list of variable attribute name: value
        
        # edge_index
        self.variables["edge_index"] = {'dtype': np.int32,
                      "long_name":	"index of pixel in pixel cloud",
                      "units": 1,
                      "valid_min": 0,
                      "valid_max": 999999,
                      "comment": "Index of pixel in the associated L2_HR_PIXC 1-D array"}
        # edge_label
        self.variables["edge_label"] = {'dtype': np.int32,
                      "long_name": "object label",
                      "units": 1,
                      "valid_min": 0,
                      "valid_max": 999999,
                      "comment": "Object label for each pixel contained in objects at top/bottom edges"}
        # edge_loc
        self.variables["edge_loc"] = {'dtype': np.int8,
                      "long_name": "edge location",
                      "flag_meanings": "bottom top both",
                      "flag_values": "0 1 2",
                      "valid_min": 0,
                      "valid_max": 2,
                      "comment": "Object edge location (0=bottom 1=top 2=both) for each pixel contained in objects at top/bottom edges"}
        # azimuth_index
        self.variables["azimuth_index"] = {'dtype': np.int32,
           'long_name': "rare interferogram azimuth index",
           'units': 1,
           'valid_min': 0,
           'valid_max': 999999,
           'comment': "Rare interferogram azimuth index"}
        # range_index
        self.variables["range_index"] = {'dtype': np.int32,
            'long_name': "rare interferogram range index",
            'units': 1,
            'valid_min': 0,
            'valid_max': 999999,
            'comment': "Rare interferogram range index"}
        # interferogram
        self.variables["interferogram"] = {'dtype': np.float,
                      "long_name": "rare interferogram",
                      "units": 1,
                      "valid_min": -999999,
                      "valid_max": 999999,
                      "coordinates": "longitude latitude",
                      "comment": "Complex unflattened rare interferogram."}
        # power_plus_y
        self.variables["power_plus_y"] = {'dtype': np.float,
                      "long_name": "power for plus_y channel",
                      "units": 1,
                      "valid_min": 0,
                      "valid_max": 999999,
                      "coordinates": "longitude latitude",
                      "comment": "Power for the plus_y channel (arbitrary units that give sigma0 when noise subtracted and normalized by the X factor)."}
        # power_minus_y
        self.variables["power_minus_y"] = {'dtype': np.float,
                      "long_name": "power for minus_y channel",
                      "units": 1,
                      "valid_min": 0,
                      "valid_max": 999999,
                      "coordinates": "longitude latitude",
                      "comment": "Power for the minus_y channel (arbitrary units that give sigma0 when noise subtracted and normalized by the X factor)."}
        # water_frac
        self.variables["water_frac"] = {'dtype': np.float,
                      "long_name": "water fraction",
                      "units": 1,
                      "valid_min": -999999,
                      "valid_max": 999999,
                      "comment": "Noisy estimate of the fraction of the pixel that is water."}
        # water_frac_uncert
        self.variables["water_frac_uncert"] = {'dtype': np.float,
                      "long_name":	"water fraction uncertainty",
                      "units": 1,
                      "valid_min": 0,
                      "valid_max": 999999,
                      "comment": "Uncertainty estimate of the water fraction estimate (width of noisy water frac estimate distribution)."}
        # classification
        self.variables["classification"] = {'dtype': np.int8,
                      "long_name": "classification",
                      "flag_meanings": "land land_near_water water_near_land open_water land_near_dark_water dark_water_edge dark_water",
                      "flag_values": "1 2 3 4 22 23 24",
                      "valid_min": 0,
                      "valid_max": 99,
                      "comment": "Flags indicating water detection results."}
        # false_detection_rate
        self.variables["false_detection_rate"] = {'dtype': np.float,
                      "long_name": "false detection rate",
                      "units": 1,
                      "valid_min": 0,
                      "valid_max": 1,
                      "comment": "Probability of falsely detecting water when there is none."}
        # missed_detection_rate
        self.variables["missed_detection_rate"] = {'dtype': np.float,
                      "long_name": "missed detection rate",
                      "units": 1,
                      "valid_min": 0,
                      "valid_max": 1,
                      "comment": "Probability of falsely detecting no water when there is water."}
        # bright_land_flag
        self.variables["bright_land_flag"] = {'dtype': np.int8, 
                      'long_name': 'bright land flag',
                      'flag_meanings': "not_bright_land bright_land",
                      'flag_values': "0 1",
                      'valid_min': 0,
                      'valid_max': 1,
                      'comment': "Flag indicating areas that are not typically water but are expected to be bright (e.g., urban areas, ice).\
                      This flag can be used to exclude detected water pixels in downstream processing."}
        # layover_impact
        self.variables["layover_impact"] = {'dtype': np.float,
                      'long_name': "layover impact",
                      'units': "m",
                      'valid_min': -999999,
                      'valid_max': 999999,
                      'comment': "Estimate of the height error caused by layover, which may not be reliable on a pixel by pixel basis,\
                      but may be useful to augment aggregated height uncertainties."}
        # eff_num_rare_looks
        self.variables["eff_num_rare_looks"] = {'dtype': np.int8,
                      'long_name': "number of rare looks",
                      'units': 1,
                      'valid_min': 0,
                      'valid_max': 999999,
                      'comment': "Number of rare looks taken (real looks or number of azimuth pixels multilooked, not effective looks)."}
        # latitude
        self.variables["latitude"] = {'dtype': np.double,
            'long_name': "latitude (positive N, negative S)",
            'standard_name': "latitude",
            'units': "degrees_north",
            'valid_min': -90.0,
            'valid_max': 90.0,
            'comment': "Geodetic latitude [-90,90] of the pixel."}
        # longitude
        self.variables["longitude"] = {'dtype': np.double,
            'long_name': "longitude (degrees East)",
            'standard_name': "longitude",
            'units': "degrees_east",
            'valid_min': -180.0,
            'valid_max': 180.0,
            'comment': "Longitude [-180,180] of the pixel."}
        # height
        self.variables["height"] = {'dtype': np.float,
            'long_name': "height above reference ellipsoid",
            'units': "m",
            'valid_min': -999999,
            'valid_max': 999999,
            'comment': "Height of the pixel above the reference ellipsoid."}
        # cross_track
        self.variables["cross_track"] = {'dtype': np.float,
                      "long_name": "approximate cross-track location",
                      "units": "m",
                      "valid_min": -999999,
                      "valid_max": 999999,
                      "comment": "Approximate cross-track location of the pixel."}
        # pixel_area
        self.variables["pixel_area"] = {'dtype': np.float,
                      "long_name": "pixel area",
                      "units": "m^2",
                      "valid_min": 0,
                      "valid_max": 999999,
                      "comment": "Pixel area."}
        # inc
        self.variables["inc"] = {'dtype': np.float,
                      "long_name": "incidence angle",
                      "units": "degrees",
                      "valid_min": 0,
                      "valid_max": 999999,
                      "comment": "Incidence angle."}
        # phase_noise_std
        self.variables["phase_noise_std"] = {'dtype': np.float,
                      "long_name": "phase noise standard deviation",
                      "units": "radians",
                      "valid_min": -999999,
                      "valid_max": 999999,
                      "comment": "Estimate of the phase noise standard deviation."}
        # dlatitude_dphase
        self.variables["dlatitude_dphase"] = {'dtype': np.float,
                      "long_name": "sensitivity of latitude estimate to interferogram phase",
                      "units": "degrees/radian",
                      "valid_min": -999999,
                      "valid_max": 999999,
                      "comment": "Sensitivity of the latitude estimate to the interferogram phase."}
        # dlongitude_dphase
        self.variables["dlongitude_dphase"] = {'dtype': np.float,
                      "long_name": "sensitivity of longitude estimate to interferogram phase",
                      "units": "degrees/radian",
                      "valid_min": -999999,
                      "valid_max": 999999,
                      "comment": "Sensitivity of the longitude estimate to the interferogram phase."}
        # dheight_dphase
        self.variables["dheight_dphase"] = {'dtype': np.float,
                      "long_name": "sensititvity of height estimate to interferogram phase",
                      "units": "m/radian",
                      "valid_min": -999999,
                      "valid_max": 999999,
                      "comment": "Sensititvity of the height estimate to the interferogram phase."}
        # dheight_drange
        self.variables["dheight_drange"] = {'dtype': np.float,
            'long_name': "sensititvity of height estimate to range (delay)", 
            'units': "m/m",
            'valid_min': -999999,
            'valid_max': 999999,
            'comment': "Sensititvity of the height estimate to the range (delay)."}
        # darea_dheight
        self.variables["darea_dheight"] = {'dtype': np.float,
            'long_name': "sensititvity of pixel area to reference height", 
            'units': "m^2/m",
            'valid_min': -999999,
            'valid_max': 999999,
            'comment': "Sensititvity of the pixel area to the reference height."}
        # eff_num_medium_looks
        self.variables["eff_num_medium_looks"] = {'dtype': np.int32,
            'long_name': "number of medium looks", 
            'units': "1",
            'valid_min': 0,
            'valid_max': 999999,
            'comment': "Number of medium looks taken (number of pixels of the same class as this pixel in the adaptive averaging window)."}
        # model_dry_tropo_cor
        self.variables["model_dry_tropo_cor"] = {'dtype': np.float,
            'long_name': "dry troposphere vertical correction", 
            'source': "European Centre for Medium-Range Weather Forecasts",
            'units': "m",
            'valid_min': -2.5,
            'valid_max': -2.0,
            'comment': "Equivalent vertical correction due to dry troposphere delay. The reported pixel height,\
            latitude and longitude are computed after adding negative media corrections to uncorrected range along slant-range paths,\
            accounting for the differential delay between the two KaRIn antennas. The equivalent vertical correction is computed by \
            applying obliquity factors to the slant-path correction. Adding the reported correction to the reported pixel height \
            results in the uncorrected pixel height."}
        # model_wet_tropo_cor
        self.variables["model_wet_tropo_cor"] = {'dtype': np.float,
            'long_name': "wet troposphere vertical correction", 
            'source': "European Centre for Medium-Range Weather Forecasts",
            'units': "m",
            'valid_min': -0.5,
            'valid_max': 0,
            'comment': "Equivalent vertical correction due to wet troposphere delay. The reported pixel height, \
            latitude and longitude are computed after adding negative media corrections to uncorrected range along slant-range paths,\
            accounting for the differential delay between the two KaRIn antennas. The equivalent vertical correction is computed by \
            applying obliquity factors to the slant-path correction. Adding the reported correction to the reported pixel height \
            results in the uncorrected pixel height."}
        # iono_cor_gim_ka
        self.variables["iono_cor_gim_ka"] = {'dtype': np.float,
            'long_name': "ionosphere vertical correction", 
            'source': "NASA/JPL Global Ionosphere Map",
            'units': "m",
            'valid_min': -0.1,
            'valid_max': 0,
            'comment': "Equivalent vertical correction due to ionosphere delay. The reported pixel height, \
            latitude and longitude are computed after adding negative media corrections to uncorrected range along slant-range paths,\
            accounting for the differential delay between the two KaRIn antennas. The equivalent vertical correction is computed by \
            applying obliquity factors to the slant-path correction. Adding the reported correction to the reported pixel height \
            results in the uncorrected pixel height."}
        # height_cor_xover
        self.variables["height_cor_xover"] = {'dtype': np.float,
            'long_name': "corssover calibration height correction", 
            'units': "m",
            'valid_min': -1000,
            'valid_max': 1000,
            'comment': "Equivalent height correction estimated from crossover calibration.  The correction is applied before\
            geolocation in terms of roll, baseline dilation, etc., but reported as an equivalent height correction."}
        # geoid
        self.variables["geoid"] = {'dtype': np.float,
            'long_name': "geoid height", 
            'standard_name': "geoid_height_above_reference_ellipsoid",
            'source': "EGM2008",
            'institution': "GSFC",
            'units': "m",
            'valid_min': -999999,
            'valid_max': 999999,
            'comment': "Geoid height above the reference ellipsoid.  The value is computed from the EGM2008\
            geoid model with a correction to refer the value to the mean tide system (i.e., includes the zero-frequency permanent tide)"}
        # solid_earth_tide
        self.variables["solid_earth_tide"] = {'dtype': np.float,
            'long_name': "solid earth tide", 
            'source': "Cartwright and Edden [1973] Corrected tables of tidal harmonics - J. Geophys. J. R. Astr. Soc., 33, 253-264",
            'units': "m",
            'valid_min': -999999,
            'valid_max': 999999,
            'comment': "Cartwright/Taylor solid earth tide; units are meters above the reference ellipsoid. \
            The zero-frequency permanent tide is not included."}
        # load_tide_sol1
        self.variables["load_tide_sol1"] = {'dtype': np.float,
            'long_name': "geocentric load tide height (solution 1)", 
            'institution': "GSFC",
            'units': "m",
            'valid_min': -999999,
            'valid_max': 999999,
            'comment': "GOT4.10; units are meters above the reference ellipsoid. This term is reported for reference but\
            is not removed from the pixel geolocation."}
        # load_tide_sol2
        self.variables["load_tide_sol2"] = {'dtype': np.float,
            'long_name': "geocentric load tide height (solution 2)", 
            'institution': "LEGOS/CNES",
            'units': "m",
            'valid_min': -999999,
            'valid_max': 999999,
            'comment': "FES2014; units are meters above the reference ellipsoid. This term is reported for reference but\
            is not removed from the pixel geolocation."}
        # pole_tide
        self.variables["pole_tide"] = {'dtype': np.float,
            'long_name': "geocentric pole tide height", 
            'units': "m",
            'valid_min': -999999,
            'valid_max': 999999,
            'comment': "Geocentric pole tide height; units are meters above the reference ellipsoid."}
        # pixc_qual
        self.variables["pixc_qual"] = {'dtype': np.int8,
            'flag_meanings': "good bad",
            'flag_values': "0 1",
            'valid_min': 0,
            'valid_max': 1,
            'comment': "Quality flag for pixel cloud data"}
        # nadir_time
        self.variables["nadir_time"] = {'dtype': np.double,
            'long_name': "time in UTC",
            'standard_name': "time",
            'calendar': "gregorian",
            'tai_utc_difference': "[Value of TAI-UTC at time of first record]",
            'leap_second': "YYYY-MM-DD hh:mm:ss",
            'units': "seconds since 2000-01-01 00:00:00.000",
            'comment': "Time of measurement in seconds in the UTC time scale since 1 Jan 2000 00:00:00 UTC.\
            [tai_utc_difference] is the difference between TAI and UTC reference time (seconds) for the first\
            measurement of the data set. If a leap second occurs within the data set, the attribute leap_second\
            is set to the UTC time at which the leap second occurs."}
        # nadir_time_tai
        self.variables["nadir_time_tai"] = {'dtype': np.double,
            'long_name': "time in TAI",
            'standard_name': "time",
            'calendar': "gregorian",
            'units': "seconds since 2000-01-01 00:00:00.000",
            'comment': "Time of measurement in seconds in the TAI time scale since 1 Jan 2000 00:00:00 TAI.\
            This time scale contains no leap seconds. The difference (in seconds) with time in UTC is given by \
            the attribute [time:tai_utc_difference]."}
        # nadir_latitude
        self.variables["nadir_latitude"] = {'dtype': np.double,
            'long_name': "latitude (positive N, negative S) of the spacecraft",
            'units': "degrees north",
            'valid_min': -90.0,
            'valid_max': 90.0,
            'comment': "Geodetic latitude of the KMSF origin with respect to the reference ellipsoid."}
        # nadir_longitude
        self.variables["nadir_longitude"] = {'dtype': np.double,
            'long_name': "longitude (degrees East) of the spacecraft",
            'units': "degrees east",
            'valid_min': -180.0,
            'valid_max': 180.0,
            'comment': "Longitude of the KMSF origin."}
        # nadir_x
        self.variables["nadir_x"] = {'dtype': np.double,
            'long_name': "x coordinate of the spacecraft in the ECEF frame",
            'units': "m",
            'valid_min': -10000000.0,
            'valid_max': 10000000.0,
            'comment': "x coordinate of the KMSF origin in the ECEF frame."}
        # nadir_y
        self.variables["nadir_y"] = {'dtype': np.double,
            'long_name': "y coordinate of the spacecraft in the ECEF frame",
            'units': "m",
            'valid_min': -10000000.0,
            'valid_max': 10000000.0,
            'comment': "y coordinate of the KMSF origin in the ECEF frame."}
        # nadir_z
        self.variables["nadir_z"] = {'dtype': np.double,
            'long_name': "z coordinate of the spacecraft in the ECEF frame",
            'units': "m",
            'valid_min': -10000000.0,
            'valid_max': 10000000.0,
            'comment': "z coordinate of the KMSF origin in the ECEF frame."}
        # nadir_vx
        self.variables["nadir_vx"] = {'dtype': np.double,
            'long_name': "x component of the spacecraft velocity in the ECEF frame",
            'units': "m/s",
            'valid_min': -10000.0,
            'valid_max': 10000.0,
            'comment': "KMSF velocity component in x direction in the ECEF frame."}
        # nadir_vy
        self.variables["nadir_vy"] = {'dtype': np.double,
            'long_name': "y component of the spacecraft velocity in the ECEF frame",
            'units': "m/s",
            'valid_min': -10000.0,
            'valid_max': 10000.0,
            'comment': "KMSF velocity component in y direction in the ECEF frame."}
        # nadir_vz
        self.variables["nadir_vz"] = {'dtype': np.double,
            'long_name': "z component of the spacecraft velocity in the ECEF frame",
            'units': "m/s",
            'valid_min': -10000.0,
            'valid_max': 10000.0,
            'comment': "KMSF velocity component in z direction in the ECEF frame."}
        # nadir_plus_y_antenna_x
        self.variables["nadir_plus_y_antenna_x"] = {'dtype': np.double,
            'long_name': "x coordinate of plus_y antenna phase center in the ECEF frame",
            'units': "m",
            'valid_min': -10000000.0,
            'valid_max': 10000000.0,
            'comment': "x coordinate of plus_y antenna phase center in the ECEF frame."}
        # nadir_plus_y_antenna_y
        self.variables["nadir_plus_y_antenna_y"] = {'dtype': np.double,
            'long_name': "y coordinate of plus_y antenna phase center in the ECEF frame",
            'units': "m",
            'valid_min': -10000000.0,
            'valid_max': 10000000.0,
            'comment': "y coordinate of plus_y antenna phase center in the ECEF frame."}
        # nadir_plus_y_antenna_z
        self.variables["nadir_plus_y_antenna_z"] = {'dtype': np.double,
            'long_name': "z coordinate of plus_y antenna phase center in the ECEF frame",
            'units': "m",
            'valid_min': -10000000.0,
            'valid_max': 10000000.0,
            'comment': "z coordinate of plus_y antenna phase center in the ECEF frame."}
        # nadir_minus_y_antenna_x
        self.variables["nadir_minus_y_antenna_x"] = {'dtype': np.double,
            'long_name': "x coordinate of minus_y antenna phase center in the ECEF frame",
            'units': "m",
            'valid_min': -10000000.0,
            'valid_max': 10000000.0,
            'comment': "x coordinate of minus_y antenna phase center in the ECEF frame."}
        # nadir_minus_y_antenna_y
        self.variables["nadir_minus_y_antenna_y"] = {'dtype': np.double,
            'long_name': "y coordinate of minus_y antenna phase center in the ECEF frame",
            'units': "m",
            'valid_min': -10000000.0,
            'valid_max': 10000000.0,
            'comment': "y coordinate of minus_y antenna phase center in the ECEF frame."}
        # nadir_minus_y_antenna_z
        self.variables["nadir_minus_y_antenna_z"] = {'dtype': np.double,
            'long_name': "z coordinate of minus_y antenna phase center in the ECEF frame",
            'units': "m",
            'valid_min': -10000000.0,
            'valid_max': 10000000.0,
            'comment': "z coordinate of minus_y antenna phase center in the ECEF frame."}
        # nadir_sc_event_flag
        self.variables["nadir_sc_event_flag"] = {'dtype': np.int8,
            'flag_meanings': "nominal not_nominal",
            'flag_values': "0 1",
            'valid_min': 0,
            'valid_max': 1,
            'comment': "Spacecraft event flag"}
        # nadir_tvp_qual
        self.variables["nadir_tvp_qual"] = {'dtype': np.int8,
            'flag_meanings': "good bad",
            'flag_values': "0 1",
            'valid_min': 0,
            'valid_max': 1,
            'comment': "Quality flag for TVP data"}
        
        # 4 - Init metadata specific to LakeTile_pixcvec file
        # 4.1 - Update general metadata
        self.metadata["title"] = "Level 2 KaRIn high rate lake tile vector product"
        self.metadata["short_name"] = "SWOT_L2_HR_LakeTile_Edge"
        # code commenté self.metadata["references"] = ""
        self.metadata["reference_document"] = "SWOT-TN-CDM-0673-CNES"
        # 4.2 - Metadata retrieved from L2_HR_PIXC product
        # The values are defined later in TBD function
        self.metadata["cycle_number"] = -9999
        self.metadata["pass_number"] = -9999
        self.metadata["tile_number"] = -9999
        self.metadata["swath_side"] = ""
        self.metadata["tile_name"] = ""
        self.metadata["time_coverage_start"] = ""
        self.metadata["time_coverage_end"] = ""
        self.metadata["inner_first_latitude"] = -9999.0
        self.metadata["inner_first_longitude"] = -9999.0
        self.metadata["inner_last_latitude"] = -9999.0
        self.metadata["inner_last_longitude"] = -9999.0
        self.metadata["outer_first_latitude"] = -9999.0
        self.metadata["outer_first_longitude"] = -9999.0
        self.metadata["outer_last_latitude"] = -9999.0
        self.metadata["outer_last_longitude"] = -9999.0
        self.metadata["continent"] = ""
        self.metadata["wavelength"] = -9999.0
        self.metadata["near_range"] = -9999.0
        self.metadata["nominal_slant_range_spacing"] = -9999.0
        self.metadata["interferogram_size_range"] = -9999
        self.metadata["interferogram_size_azimuth"] = -9999
        self.metadata["looks_to_efflooks"] = -9999.0
        self.metadata["ellipsoid_semi_major_axis"] = ""
        self.metadata["ellipsoid_flattening"] = ""
        if in_pixc_metadata is not None:
            self.set_metadata_val(in_pixc_metadata)
        # 4.3 - Processing metadata
        self.metadata["xref_input_l2_hr_pixc_file"] = ""
        self.metadata["xref_l2_hr_lake_tile_param_file"] = ""
        if in_proc_metadata is not None:
            self.set_metadata_val(in_proc_metadata)
