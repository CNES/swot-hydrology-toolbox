import netCDF4 as nc
import numpy as np
import os
import ogr, osr
from scipy.interpolate import LinearNDInterpolator

__author__ = """Kevin Larnier """ \
             """<kevin.larnier@c-s.fr>"""
__version__ = """0.1"""


class TrueHeightModel():
    """ Class for application of true height to a sisimp pointcloud
    """
  
    def __init__(self, trueheight_filename, pt_lat, pt_lon, pt_height, **kwargs):
        """ Initialize the class
        
            Arguments:
                pointcloud_filename(str) : path to the pointcloud file
                trueheight_filename(str) : path to the true height file
                output_filename(str) : path to the output file
        """
        
        # Set-up
        self.th_filename = trueheight_filename
        self.th_attr = {'height' : "height",
                        'lat' : "latitude",
                        'lon' : "longitude"}
        self.pt_lon = pt_lon
        self.pt_lat = pt_lat
        self.pt_height = pt_height

        self.verbose = False
        
        # Apply keyword arguments
        for key in kwargs:

          if key == "trueheight_latitude":
            if not isinstance(kwargs[key], str):
              raise ValueError("'trueheight_latitude' argument must be a string")
            self.th_attr['lat'] = kwargs[key]
          elif key == "trueheight_longitude":
            if not isinstance(kwargs[key], str):
              raise ValueError("'trueheight_longitude' argument must be a string")
            self.th_attr['lon'] = kwargs[key]
          elif key == "trueheight_height":
            if not isinstance(kwargs[key], str):
              raise ValueError("'trueheight_height' argument must be a string")
            self.th_attr['height'] = kwargs[key]
          elif key == "verbose":
            if not isinstance(kwargs[key], bool):
              raise ValueError("'verbose' argument must be True or False")
            self.verbose = kwargs[key]
          else:
            raise ValueError("Unknown argument '%s'" % key)

    def read_netcdf_trueheight(self, filename, lon_varname, lat_varname, height_varname):
        """ Read true height from a NetCDF file
        
            Arguments:
                filename(str) : path to the pointcloud file
                lon_varname(str) : variable name for the longitude
                lat_varname(str) : variable name for the latitude
                height_varname(str) : variable name for the height
        """
      
        if self.verbose:
          print(" - Read true height from NetCDF file: %s" % filename)

        # Open NetCDF dataset
        dataset = nc.Dataset(filename)
        if dataset is None:
            raise IOError("Unable to load pointcloud file: %s" % filename)
        
        # Extract values
        self.th_lat = dataset.variables[lat_varname][:]
        self.th_lon = dataset.variables[lon_varname][:]
        self.th_height = dataset.variables[height_varname][:]
        
    def apply_trueheight(self):
        """ Compute final heights
        """
      
        if self.verbose:
            print(" - Apply true height")
          
        th_points = np.stack((self.th_lon, self.th_lat), axis=1)
        interpolator = LinearNDInterpolator(th_points, self.th_height)
        
        pt_points = np.stack((self.pt_lon, self.pt_lat), axis=1)
        interpolated_height = interpolator(pt_points)
        
        self.final_height = self.pt_height + interpolated_height
        
    def apply_model(self):
        
        # Read true height file
        if os.path.splitext(self.th_filename)[1] == ".nc":
            self.read_netcdf_trueheight(self.th_filename, self.th_attr['lon'], 
                                        self.th_attr['lat'], self.th_attr['height'])
        elif os.path.splitext(self.th_filename)[1] == ".shp":
            self.read_shapefile_trueheight(self.th_filename, self.th_attr['lon'], 
                                          self.th_attr['lat'], self.th_attr['height'])
        else:
            raise RuntimeError("True height file must be in NetCDF-4 or Shapefile format")
        
        # Apply true height
        self.apply_trueheight()    
