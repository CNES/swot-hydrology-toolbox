import sys,os,shutil
import numpy
import scipy.interpolate
import matplotlib.pylab as plt
import netCDF4 as nc
from scipy.interpolate import griddata
import shutil
#from mpl_toolkits.basemap import Basemap

deltat=0. # time lag to be applied if reference time or pass number is different from CNES reference. 

sim = type('', (), {})()
sens = type('', (), {})()
cloud = type('', (), {})()

# Read roll error and roll correction data base (simulated with the Ocean Science simulator and cross-calibation Algorithm). Roll unit is micro-radiant. 
# roll1_err: roll error for negative cross-track distance
# roll2_err: roll error for positive cross-track distance 
# roll1_cor: roll correction for negative cross-track distance
# roll2_cor: roll correction for positive cross-track distance
# This file is in the module package, add correct destination
filename='/work/ALT/swot/swotdev/desrochesd/swot-sds/sisimp/lib/data_sim_roll_v1.nc'
fid = nc.Dataset(filename, 'r')
sim.time=numpy.array(fid.variables['time'])+deltat
sim.lon_nadir=numpy.array(fid.variables['lon_nadir'])
sim.lat_nadir=numpy.array(fid.variables['lat_nadir'])
sim.roll1_err=numpy.array(fid.variables['roll1_err'])
sim.roll2_err=numpy.array(fid.variables['roll2_err'])
sim.roll1_cor=numpy.array(fid.variables['roll1_cor'])
sim.roll2_cor=numpy.array(fid.variables['roll2_cor'])
fid.close()

# Read sensor file correponding associated with the cloud point simulation
# Exemple of Seine river (From D. Desroches) is taken here
filehydro='/home/cubelmann/Data/image_point_cloud/sensor_seine.nc'
fid = nc.Dataset(filehydro, 'r')
sens.time=fid.variables['time'][:]
sens.lat=fid.variables['latitude'][:]
sens.lon=fid.variables['longitude'][:]
fid.close()


# Interpolate simulated roll variables on sensor time grid
ensind=numpy.where(((sim.time>=sens.time[0]-2.)&(sim.time<=sens.time[-1]+2.)))
f=scipy.interpolate.interp1d(sim.time[ensind],sim.roll1_err[ensind],kind='linear')
sens.roll1_err=f(sens.time)
f=scipy.interpolate.interp1d(sim.time[ensind],sim.roll2_err[ensind],kind='linear')
sens.roll2_err=f(sens.time)
f=scipy.interpolate.interp1d(sim.time[ensind],sim.roll1_cor[ensind],kind='linear')
sens.roll1_cor=f(sens.time)
f=scipy.interpolate.interp1d(sim.time[ensind],sim.roll2_cor[ensind],kind='linear')
sens.roll2_cor=f(sens.time)


# Read pixel-cloud file
# Exemple of Seine river (From D. Desroches) is taken here
filehydrocloud='/home/cubelmann/Data/image_point_cloud/intf_testLeftSwath.L2PC_seine_wind.nc'
iswath=-1 # index to specify if left swath (-1) or right swath (+1) as current format does not specify.
fid = nc.Dataset(filehydrocloud, 'r')
azimuth_index=fid.variables['azimuth_index'][:]
cloud.range=fid.variables['range'][:]
cloud.lon=fid.variables['longitude_medium'][:]
cloud.lat=fid.variables['latitude_medium'][:]
cloud.height=fid.variables['height_medium'][:]
cloud.xac=iswath*fid.variables['cross_track_welldone'][:]
fid.close()
cloud.time=sens.time[azimuth_index]


# Interpolate roll variables on pixel cloud
f=scipy.interpolate.interp1d(sens.time,sens.roll1_err,kind='linear')
cloud.roll1_err=f(cloud.time)
f=scipy.interpolate.interp1d(sens.time,sens.roll2_err,kind='linear')
cloud.roll2_err=f(cloud.time)
f=scipy.interpolate.interp1d(sens.time,sens.roll1_cor,kind='linear')
cloud.roll1_cor=f(cloud.time)
f=scipy.interpolate.interp1d(sens.time,sens.roll2_cor,kind='linear')
cloud.roll2_cor=f(cloud.time)

# Conversion in equivalent height: *xac*1e-6 for microrad to meter correction
cloud.hroll_err=numpy.zeros(len(cloud.xac))
cloud.hroll_err[cloud.xac<0] = cloud.roll1_err[cloud.xac<0]*cloud.xac[cloud.xac<0]*1e-6
cloud.hroll_err[cloud.xac>0] = cloud.roll2_err[cloud.xac>0]*cloud.xac[cloud.xac>0]*1e-6
cloud.hroll_cor=numpy.zeros(len(cloud.xac))
cloud.hroll_cor[cloud.xac<0] = cloud.roll1_cor[cloud.xac<0]*cloud.xac[cloud.xac<0]*1e-6
cloud.hroll_cor[cloud.xac>0] = cloud.roll2_cor[cloud.xac>0]*cloud.xac[cloud.xac>0]*1e-6



# Write roll error height signature in pixel-cloud file
fid = nc.Dataset(filehydrocloud, 'a')
rec_dim=fid.dimensions['water_record']

newvar = fid.createVariable('height_error','f4','water_record')
newvar[:]=cloud.hroll_err
newvar.units= "m"
newvar.long_name = 'height error induced by uncorrected roll error'
newvar.FillValue=-9.99e+09

newvar = fid.createVariable('height_res_error','f4','water_record')
newvar[:]=cloud.hroll_err
newvar.units= "m"
newvar.long_name = 'height  residual error after roll correction'
newvar.FillValue=-9.99e+09
fid.close()




