#!/usr/bin/env python

import sys
import os
import os.path

import cnes.modules.geoloc.lib.tools as lib
import cnes.modules.geoloc.lib.my_netcdf_file as myCdf
import cnes.modules.geoloc.lib.gdem as gdem
import cnes.modules.geoloc.lib.function_mean_h as function_mean_h
import cnes.modules.geoloc.lib.pixel_cloud as pixel_cloud
import cnes.modules.geoloc.lib.geoloc as geoloc

import numpy as np

from scipy import spatial
from timeit import default_timer as timer
import itertools as itertools

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from math import atan2
import timeit

import cProfile 
import profile 

def normalize(vect):
    norm = np.sqrt(np.dot(vect,vect))
    return vect/norm #don t use me with np.zeroes(#)

def dist(a,b):
    return np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)
## last JPL version
#~ input_intf_path = "/work/ALT/swot/swotdev/SIMU_JPL/simu_US_land_-5dB/simu1/output/0001/cycle_0001_pass_0107_RightSwath_nlcd50_darkWater_25m_CBE/intf_truth.RightSwath.nc"
#~ input_sensor_path = "/work/ALT/swot/swotdev/SIMU_JPL/simu_US_land_-5dB/simu1/output/0001/cycle_0001_pass_0107_RightSwath_nlcd50_darkWater_25m_CBE/sensor.nc"
## november binaries, SAM simulator
#~ input_intf_path = "/work/ALT/swot/swothr/users/desrochesd/runs/old/test_liv12_06.ref/makeInterferogramJPL_1.01/outputs/intf_test.LeftSwath.nc"
#~ input_sensor_path = "/work/ALT/swot/swothr/users/desrochesd/runs/old/test_liv12_06.ref/makeSensor_2.00/outputs/sensorfile.nc"


#~ intf_main = myCdf.myNcReader(input_intf_path, dim=2)
#~ df, nb_az, nb_rg = intf_main.getVarValue2d('record', 'nr_pixels')
#~ dataframe = df
        
#~ near_range = intf_main.content.attrs['near_range']
#~ range_spacing = intf_main.content.attrs['range_spacing']

        
#~ no_layover_x = dataframe.no_layover_x.values.reshape(nb_rg,  nb_az)
#~ no_layover_y = dataframe.no_layover_y.values.reshape(nb_rg,  nb_az)
#~ no_layover_z = dataframe.no_layover_z.values.reshape(nb_rg,  nb_az)
#~ height = dataframe.height.values.reshape(nb_rg,  nb_az)


#~ sensor_main = myCdf.myNcReader(input_sensor_path)
#~ x_sensor = sensor_main.getVarValue("x")
#~ y_sensor = sensor_main.getVarValue("y")
#~ z_sensor = sensor_main.getVarValue("z")

#~ vx_sensor = sensor_main.getVarValue("velocity_unit_x")
#~ vy_sensor = sensor_main.getVarValue("velocity_unit_y")
#~ vz_sensor = sensor_main.getVarValue("velocity_unit_z")

#~ ind=1000
#~ col=1000

#~ s=np.array([x_sensor[ind], y_sensor[ind], z_sensor[ind]])
#~ vs=np.array([vx_sensor[ind], vy_sensor[ind], vz_sensor[ind]])
#~ p=np.array([no_layover_x[col,ind], no_layover_y[col,ind], no_layover_z[col,ind]])
#~ h=height[col,ind]


#~ R_target = float(near_range + col*range_spacing)
#~ h_target = h

p = np.array([  604036.2612028 ,-5151738.12882137 ,3700899.70950969])

s = np.array([  691911.1260945 ,-5874186.70691969 ,4224481.28190375])

vs = np.array([ 0.27379512 ,-0.54234502 ,-0.79429094])



R_target = 896544.255928 
#~ h_target = 1279.56469727

R_e = 6378137.0 
R_p = 6356752.31425 

# Used height obtained from ecef2llh conversion using improved geolocation ellipsoid to be consistent
h_true  = geoloc.convert_ecef2llh(p[0], p[1], p[2], R_e, R_p)[2]

h_target = h_true + 2.

#~ p=[p,p,p]
#~ s=[s,s,s]
#~ vs=[vs,vs,vs]
#~ h_target=[h_target,h_target,h_target]


start = timer()
p_corr = geoloc.pointcloud_ellipsoidgeoloc_improved(p, s, vs, R_target, h_target, keep_Doppler=False, keep_R=False, verbose = False, numerical_proj = False)
print(timer()-start)

start = timer()
p_corr_true = geoloc.pointcloud_ellipsoidgeoloc_improved(p, s, vs, R_target, h_target, keep_Doppler=True, keep_R=True, verbose = False, numerical_proj = False)
print(timer()-start)

start = timer()
p_corr_taylor = geoloc.pointcloud_ellipsoidgeoloc_improved_Taylortheta2mu1(p, s, vs, R_target, h_target, h_true, keep_Doppler=False, keep_R=False, verbose = False)
print(timer()-start)


print('')
print('p : ', p)
print(geoloc.convert_ecef2llh(p[0], p[1], p[2], R_e, R_p)[2])

print('')
print('p_corr : ', p_corr, dist(p, p_corr))
print(geoloc.convert_ecef2llh(p_corr[0], p_corr[1], p_corr[2], R_e, R_p)[2])

print('')
print('p_corr_true : ', p_corr_true, dist(p, p_corr_true))
print(geoloc.convert_ecef2llh(p_corr_true[0], p_corr_true[1], p_corr_true[2], R_e, R_p)[2])

print('')
print('p_corr_taylor : ', p_corr_taylor, dist(p, p_corr_taylor))
print(geoloc.convert_ecef2llh(p_corr_taylor[0], p_corr_taylor[1], p_corr_taylor[2], R_e, R_p)[2])

print('')


#~ print('')
#~ print(p)
#~ print('')
#~ print(p_corr)
#~ print('')
#~ print(p_corr_true)
