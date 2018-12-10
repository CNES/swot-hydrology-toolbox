import sys,os,shutil
import numpy
import scipy.interpolate
import matplotlib.pylab as plt
import netCDF4 as nc
from scipy.interpolate import griddata
import shutil


deltat=0. # time lag to be applied if reference time or pass number is different from CNES reference. 


class Roll_module(object):

    def __init__(self, In_filename):
        try:
            fid = nc.Dataset(In_filename, 'r')
        except:
            raise Exception("Roll file not found")
            
        self.time=numpy.array(fid.variables['time'])+deltat
        self.lon_nadir=numpy.array(fid.variables['lon_nadir'])
        self.lat_nadir=numpy.array(fid.variables['lat_nadir'])
        self.roll1_err=numpy.array(fid.variables['roll1_err'])
        self.roll2_err=numpy.array(fid.variables['roll2_err'])
        self.roll1_cor=numpy.array(fid.variables['roll1_cor'])
        self.roll2_cor=numpy.array(fid.variables['roll2_cor'])
        fid.close()


    def interpolate_roll_on_sensor_grid(self, sensor_time):

        # Interpolate simulated roll variables on sensor time grid
        ensind=numpy.where(((self.time>=sensor_time[0]-2.)&(self.time<=sensor_time[-1]+2.)))
        f=scipy.interpolate.interp1d(self.time[ensind],self.roll1_err[ensind],kind='linear')
        self.roll1_err_sens=f(sensor_time)
        f=scipy.interpolate.interp1d(self.time[ensind],self.roll2_err[ensind],kind='linear')
        self.roll2_err_sens=f(sensor_time)
        f=scipy.interpolate.interp1d(self.time[ensind],self.roll1_cor[ensind],kind='linear')
        self.roll1_cor_sens=f(sensor_time)
        f=scipy.interpolate.interp1d(self.time[ensind],self.roll2_cor[ensind],kind='linear')
        self.roll2_cor_sens=f(sensor_time)
        
    def interpolate_roll_on_pixelcloud(self, sensor_time, cloud_time):

        # Interpolate roll variables on pixel cloud
        f=scipy.interpolate.interp1d(sensor_time,self.roll1_err_sens,kind='linear')
        self.roll1_err_cloud=f(cloud_time)
        f=scipy.interpolate.interp1d(sensor_time,self.roll2_err_sens,kind='linear')
        self.roll2_err_cloud=f(cloud_time)
        f=scipy.interpolate.interp1d(sensor_time,self.roll1_cor_sens,kind='linear')
        self.roll1_cor_cloud=f(cloud_time)
        f=scipy.interpolate.interp1d(sensor_time,self.roll2_cor_sens,kind='linear')
        self.roll2_cor_cloud=f(cloud_time)        
