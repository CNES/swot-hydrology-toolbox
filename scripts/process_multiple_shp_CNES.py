#!/usr/bin/env python

'''
 This is a python wrapper for processing multiple shapefiles through the
 CNES simulator. The wrapper is useful if you have a time series of model
 output. Specify inputs below.

 Author: Nicholas J. Elmer
         NASA Marshall Space Flight Center
         Huntsville, Alabama, USA
         nicholas.j.elmer@nasa.gov

 Date: November 26, 2019

 This wrapper follows these steps:
 1) Create the SWOT orbit passplan.txt file (if it has not already
    been created). The passplan.txt file specifies the time of each SWOT 
    overpass over the domain extent specified in parameter_orbit.rdf.

 2) Read the passplan.txt cycle number, orbit number, and overpass time.

 3) Find the model output shapefile corresponding to the correct overpass
    time. If model output is hourly, you should set 'roundhour = True' to
    ensure the correct model output file is selected.

 4) Update the paramater_sisimp.rdf file with the correct shapefile path,
    cycle number, and orbit number.

 5) Run the simulator, repeating until all valid (i.e., corresponding to
    all passplan.txt cycle and orbits) shapefiles have been processed.
'''

import os
import numpy as np
import datetime

### INPUTS ###
#-------------------------------------------------------------------------
# name of test case with height shapefiles
# (e.g., river_and_lake, afrique, etc.)
testcase = 'mytestcase'

# shapefile filename date in datetime convention.
# See https://docs.python.org/2/library/datetime.html for more information.
modeltimefmt = '%Y-%m-%d_%H:%M:%S'

# shapefile filenames should follow modeltime+shapefile_base convention.
# If not, you will have to modify Line 124.
# Do not include file extension .shp in shapefile_base. 
# Example: If shapefile filename is 2019-11-26_12:00:00_myheights.shp, then:
#   modeltimefmt = '%Y-%m-%d_%H:%M:%S'
#   shapefile_base = '_myheights'
shapefile_base = '_myheights' 

# Round to nearest hour. Useful if model output is hourly.
roundhour = True

#Specify filepaths
#file path to parameter_sisimp.rdf file
#  (typically ./rdf/parameter_sisimp.rdf)
rdfpath = './rdf/parameter_sisimp.rdf.test'

#file path to passplan.txt
#  (typically ./output/orbit/passplan.txt)
planpath = './output/orbit/passplan.txt'

#file path to mask/height shapefiles
#  (typically $SWOT_HYDROLOGY_TOOLBOX/test/%s/data)
maskpath = '$SWOT_HYDROLOGY_TOOLBOX/test/%s/data' %(testcase)












### NO EDITS BELOW THIS LINE ###
#-------------------------------------------------------------------------

# Check if passplan.txt exists.
# If not, run select_orbit_cnes.py using parameter_orbit.rdf
print('Checking passplan...')
if not os.path.exists(planpath):
  print('Passplan does not exist. Running parameter_orbit.rdf...')
  print('...Stay tuned...')
  os.system('python $SWOT_HYDROLOGY_TOOLBOX/select_orbit_cnes/select_orbit_cnes.py rdf/parameter_orbit.rdf output/orbit')

# open passplan.txt to read orbit information
cycle = np.genfromtxt(planpath, dtype=int, usecols=(1))
orbit = np.genfromtxt(planpath, dtype=int, usecols=(2))
date = np.genfromtxt(planpath, dtype=str, usecols=(6))
time = np.genfromtxt(planpath, dtype=str, usecols=(7))
ntimes = len(time)

# process each shapefile using the CNES simulator
for i in range(ntimes):
  
  # convert passplan time to datetime tuple
  proctime = datetime.datetime.strptime(date[i]+' '+time[i], '%Y-%m-%d %H:%M:%S')
  # round to the nearest hour, if requested, and then convert to model output
  #     shapefile timestamp convention
  if roundhour:
   roundtime = proctime.replace(second=0, minute=0, hour=proctime.hour) + \
     datetime.timedelta(hours=proctime.minute//30)
   modeltime = roundtime.strftime(modeltimefmt)
  else:
   modeltime = proctime.strftime(modeltimefmt)
  print(modeltime, orbit[i], cycle[i])

  # Create backup of rdf, named parameter_sisimp.rdf.orig
  if not os.path.exists(rdfpath+'.orig'):
    os.system('mv %s %s' %(rdfpath,rdfpath+'.orig'))
  else:
    os.system('rm %s' %(rdfpath))
  # Then, create modified rdf, named parameter_sisimp.rdf
  lines = open(rdfpath+'.orig').read().splitlines()
  for j, line in zip(range(len(lines)),lines):
    if line[0:16] == 'Shapefile path =':
      lines[j] = 'Shapefile path = %s/%s%s' %(maskpath,modeltime,shapefile_base)
    elif line[0:16] == 'Multiple orbit =':
      lines[j] = 'Multiple orbit = no'
    elif line[0:7] == 'Orbit =':
      lines[j] = 'Orbit = %i' %(orbit[i])
    elif line[0:14] == 'Cycle number =':
      lines[j] = 'Cycle number = %i' %(cycle[i])
  open(rdfpath, 'w').write('\n'.join(lines))

  # run CNES simulator using modified parameter_sisimp.rdf
  print('...Processing with simulator...')
  os.system('python $SWOT_HYDROLOGY_TOOLBOX/sisimp/proc_sisimp.py rdf/parameter_sisimp.rdf')

print('All shapefiles have been processed!')









