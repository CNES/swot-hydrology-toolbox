'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


import netCDF4 as nc
import numpy
import scipy.interpolate
from matplotlib import pyplot 
import matplotlib.pyplot as plt


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


if __name__ == "__main__":
    
    roll_file = '/work/ALT/swot/swototar/XOVERCAL/roll_phase/roll_phase_c016_p574.nc'
    fid = nc.Dataset(roll_file, 'r')
    
    x = numpy.array(fid.variables['time_karin'])
        
    sl1 = numpy.array(fid.variables['sl1'])
    sl2 = numpy.array(fid.variables['sl2'])
    
    slopeleft = (sl1-sl2)
    sloperight = (sl1+sl2)
    #~ slopeleft_p = (sl1p-sl2p)
    #~ sloperight_p = (sl1p+sl2p)
        
    error_left = slopeleft*10e-3*10.
    #~ error_left_p = slopeleft_p*10e-3*10.
    error_right = sloperight*10e-3*10.
    #~ error_right_p = sloperight_p*10e-3*10.
    
    plt.figure()
    plt.plot(x, error_right)
    #~ plt.plot(xp, error_right_p)
    plt.show()   
