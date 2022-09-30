'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''

import os
import netCDF4 as nc
import numpy
import scipy.interpolate
from matplotlib import pyplot
import matplotlib.pyplot as plt
import numpy as np
import re
import lib.my_api as my_api

# time lag to be applied if reference time or pass number is different from CNES
# reference.
deltat=0.

class Roll_module(object):

    def __init__(self, In_repository):
        self.In_repository = In_repository

    def get_roll_file_associated_to_orbit_and_cycle(
            self, orbit_number, cycle_number, delta_time=0):
        root_name = "GLOBAL_swot292_c"

        ############### VERY DIRTY FIX #################
        # ~ dt = 0
        # ~ orbit_number_shifted = orbit_number
        orbit_number_shifted = orbit_number-330
        dt = 1018647.36577032
        if orbit_number_shifted < 1:
            orbit_number_shifted+= 584
            dt = -784049.790623216
        ############### VERY DIRTY FIX #################

        file_name = os.path.join(self.In_repository, root_name \
                                 + str(cycle_number).zfill(2) + "_p" \
                                 + str(orbit_number_shifted).zfill(3)) +".nc"
        my_api.printInfo("Roll file used : %s " % file_name)

        self.read_roll_file(file_name, delta_time = delta_time+dt)

    def read_roll_file(self, In_filename, delta_time = 0.):
        try:
            fid = nc.Dataset(In_filename, 'r')
        except:
            raise Exception("Roll file not found")

        self.time=numpy.array(fid.variables['time'])+delta_time
        self.lon_nadir=numpy.array(fid.variables['lon_nadir'])
        self.lat_nadir=numpy.array(fid.variables['lat_nadir'])
        self.roll1_err=numpy.array(fid.variables['roll_err'])
        self.xtrack=numpy.array(fid.variables['x_ac'])*1e3

        fid.close()

    def interpolate_roll_on_sensor_grid(self, sensor_time):
        self.roll1_err_sens = []
        self.lon_sens = []
        self.lat_sens = []

        for i in range(self.roll1_err.shape[1]):
            # Interpolate simulated roll variables on sensor time grid
            ensind=numpy.where(((self.time>=sensor_time[0]-2.) \
                                & (self.time<=sensor_time[-1]+2.)))
            try:
                f=scipy.interpolate.interp1d(
                    self.time[ensind], (self.roll1_err[:,i])[ensind], kind='linear')
            except:
                print("time between crossover files and simulation time are not consistent")
            self.roll1_err_sens.append(f(sensor_time))

        f2=scipy.interpolate.interp1d(
            self.time[ensind], (self.lon_nadir[:])[ensind], kind='linear')
        self.lon_sens.append(f2(sensor_time))

        f3=scipy.interpolate.interp1d(
            self.time[ensind], (self.lat_nadir[:])[ensind], kind='linear')
        self.lat_sens.append(f3(sensor_time))

    def interpolate_roll_on_pixelcloud(self, sensor_time, cloud_time, y,
                                       bounds_error=False, fill_value=0):
        # Interpolate roll variables on pixel cloud
        self.roll1_err_cloud = scipy.interpolate.interpn(
            (self.xtrack, sensor_time), np.array(self.roll1_err_sens),
            (y, cloud_time), method='linear', bounds_error=bounds_error,
            fill_value=fill_value)

if __name__ == "__main__":

    #~ roll_file_err = '/work/ALT/swot/swototar/XOVERCAL/karin/GLOBAL_swot292_c16_p114.nc'
    #~ roll_file_cor = '/work/ALT/swot/swototar/XOVERCAL/roll_phase/roll_phase_c016_p114.nc'
    roll_file_err = '/work/ALT/swot/swototar/XOVERCAL/karin/GLOBAL_swot292_c16_p244.nc'
    roll_file_cor = '/work/ALT/swot/swototar/XOVERCAL/roll_phase/roll_phase_c016_p244.nc'

    fid_err = nc.Dataset(roll_file_err, 'r')
    fid_cor = nc.Dataset(roll_file_cor, 'r')

    time_err = numpy.array(fid_err.variables['time'])
    roll_err = numpy.array(fid_err.variables['roll_err'])*0
    phase_err = numpy.array(fid_err.variables['phase_err'])
    x_ac = (numpy.array(fid_err.variables['x_ac'])[:])*10e3

    time_cor = numpy.array(fid_cor.variables['time_karin'])
    sl1 = (numpy.array(fid_cor.variables['sl1']))*10e-6
    sl2 = (numpy.array(fid_cor.variables['sl2']))*10e-6

    error_left = roll_err[:,0] + phase_err[:,0]
    error_right = roll_err[:,-1] + phase_err[:,-1]

    corr_left = -(sl1-sl2)*x_ac[0]
    corr_right = -(sl1+sl2)*x_ac[-1]

    plt.figure()
    plt.plot(phase_err[:,0])
    plt.figure()
    plt.plot(roll_err[:,0])

    #~ plt.figure()
    #~ plt.plot(error_left)
    #~ plt.figure()
    #~ plt.plot(error_right)
    #~ plt.figure()
    #~ plt.plot(corr_left)
    #~ plt.figure()
    #~ plt.plot(corr_right)

    #~ plt.plot(x, error_left)
    #~ plt.plot(xp, error_right_p)
    plt.show()
