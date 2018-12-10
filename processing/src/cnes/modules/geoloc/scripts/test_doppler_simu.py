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
import scipy.optimize as opt
from math import atan2
import timeit

import cProfile 
import profile 

def normalize(vect):
    norm = np.sqrt(np.dot(vect,vect))
    return vect/norm #don t use me with np.zeroes(#)
    
    
## last JPL version
input_intf_path = "/work/ALT/swot/swotdev/SIMU_JPL/simu_US_land_-5dB/simu1/output/0001/cycle_0001_pass_0107_RightSwath_nlcd50_darkWater_25m_CBE/intf_truth.RightSwath.nc"
input_sensor_path = "/work/ALT/swot/swotdev/SIMU_JPL/simu_US_land_-5dB/simu1/output/0001/cycle_0001_pass_0107_RightSwath_nlcd50_darkWater_25m_CBE/sensor.nc"
## november binaries, SAM simulator
#~ input_intf_path = "/work/ALT/swot/swothr/users/desrochesd/runs/old/test_liv12_06.ref/makeInterferogramJPL_1.01/outputs/intf_test.LeftSwath.nc"
#~ input_sensor_path = "/work/ALT/swot/swothr/users/desrochesd/runs/old/test_liv12_06.ref/makeSensor_2.00/outputs/sensorfile.nc"


intf_main = myCdf.myNcReader(input_intf_path, dim=2)
df, nb_az, nb_rg = intf_main.getVarValue2d('record', 'nr_pixels')
dataframe = df
        
near_range = intf_main.content.attrs['near_range']
range_spacing = intf_main.content.attrs['range_spacing']

        
no_layover_x = dataframe.no_layover_x.values.reshape(nb_rg,  nb_az)
no_layover_y = dataframe.no_layover_y.values.reshape(nb_rg,  nb_az)
no_layover_z = dataframe.no_layover_z.values.reshape(nb_rg,  nb_az)
height = dataframe.height.values.reshape(nb_rg,  nb_az)


sensor_main = myCdf.myNcReader(input_sensor_path)
x_sensor = sensor_main.getVarValue("x")
y_sensor = sensor_main.getVarValue("y")
z_sensor = sensor_main.getVarValue("z")

vx_sensor = sensor_main.getVarValue("velocity_unit_x")
vy_sensor = sensor_main.getVarValue("velocity_unit_y")
vz_sensor = sensor_main.getVarValue("velocity_unit_z")

ind=10000
col=2000

s=np.array([x_sensor[ind], y_sensor[ind], z_sensor[ind]])
vs=np.array([vx_sensor[ind], vy_sensor[ind], vz_sensor[ind]])
p=np.array([no_layover_x[col,ind], no_layover_y[col,ind], no_layover_z[col,ind]])
h=height[col,ind]

#~ ## Verification Doppler 0 

print(s, vs, p)
print(h)
# sensor position:
s_hat = normalize(s)
# An adapted orthonormal basis:
v_hat = normalize(vs) #points along track
u_hat = normalize(np.cross(s_hat, v_hat)) #s and v are almost ortho no danger of colinear (points roughly cross track)
w_hat = np.cross(u_hat, v_hat) # points roughly to nadir
r_s = np.linalg.norm(s)
delta = np.dot(p-s,v_hat)
print('doppler = ', delta)

