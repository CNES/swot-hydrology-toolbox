from netCDF4 import Dataset
import numpy
import xarray as xr
import sys
import argparse
import shutil

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('gdem', help='gdem file', type=str)
    parser.add_argument('output', help='subsampled output gdem', type=str)
    args = parser.parse_args()


    gdem = Dataset(args.gdem, "r")
    output = Dataset(args.output, "w")

    subsample_factor = 10

    latitude = gdem.variables["latitude"][::subsample_factor]
    longitude = gdem.variables["longitude"][::subsample_factor]
    elevation = gdem.variables["elevation"][::subsample_factor, ::subsample_factor]
    landtype = gdem.variables["landtype"][::subsample_factor, ::subsample_factor]
    
    output.createDimension("latitude", latitude.size)
    output.createDimension("longitude", longitude.size)

    print(latitude.shape)
    print(longitude.shape)
    print(elevation.shape)
    print(landtype.shape)
    
    
    output.createVariable("latitude", "double", dimensions=("latitude"), fill_value=-9990000000.)
    output.createVariable("longitude", "double", dimensions=("longitude"), fill_value=-9990000000.)
    output.createVariable("elevation", "double", dimensions=("latitude", "longitude"), fill_value=-128)
    output.createVariable("landtype", "byte", dimensions=("latitude", "longitude"), fill_value=-9990000000.)

    output.variables["latitude"][:] = latitude
    output.variables["longitude"][:] = longitude
    output.variables["elevation"][:] = elevation
    output.variables["landtype"][:] = landtype

    output.close()
