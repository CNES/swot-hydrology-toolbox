#!/usr/bin/env python

'''
 This is a python wrapper for processing multiple shapefiles through the
 CNES simulator. The wrapper is useful if you have a time series of model
 output. Specify inputs below.
 
 Script adapted by Manon Quellec C-S (manon.quellec@c-s.fr) 
 from that written by Nicholas J. Elmer NASA Marshall Space Flight Center 
 (nicholas.j.elmer@nasa.gov) 
 
 Date: January 23, 2020
 
 This wrapper follows these steps:
 1) Create the SWOT orbit passplan.txt file (if it has not already
    been created). The passplan.txt file specifies the time of each SWOT 
    overpass over the domain extent specified in parameter_orbit.rdf.
 2) Read the passplan.txt cycle number and orbit number.
 3) Find the water mask shapefile whose date corresponds, or is closest 
    to each date in the passplan file.
 4) Update the paramater_sisimp.rdf file with the correct shapefile path,
    cycle number, and orbit number.
 5) Run the simulator (SISIMP), repeating until all valid (i.e., 
    corresponding to all passplan.txt cycle and orbits) shapefiles have 
    been processed.
 6) Run the simulator (LakeTile) from the PixelClouds generated in the 
    previous step.
'''

import os
import numpy as np
import datetime


#### INPUTS ####
##-------------------------------------------------------------------------
## Name of test case folder (containing the watermask shapefiles folder)
testcase = 'LT6_Cas_Der'

## Shapefile filename date in datetime convention.
## See https://docs.python.org/2/library/datetime.html for more information.
## Example: If shapefile filename is "20170923_Der_S2_SVM_SERTIT.shp", then:
modeltimefmt = '%Y%m%d'
shapefile_basename = '_Der_S2_SVM_SERTIT_WGS84' 

## Specify filepaths :
## file path to parameter_orbit.rdf file and output orbit 
orbit_rdf_path = '%s/config/parameter_orbit.rdf' %(testcase)
output_orbit_path = '%s/1_select_orbit' %(testcase)

## file path to parameter_sisimp.rdf file
sisimp_rdf_path = '%s/config/parameter_sisimp_light.rdf' %(testcase)

## file path to multi_lake_tile_command.cfg file
lake_tile_cfg_path = '%s/config/multi_lake_tile_command.cfg' %(testcase)

## file path to passplan.txt
planpath = '%s/1_select_orbit/passplan.txt' %(testcase)

## path to watermask shapefiles
maskpath = '%s/0_water_mask/' %(testcase)






#### NO EDITS BELOW THIS LINE ####
##-------------------------------------------------------------------------

def is_date_in_list_of_string(date, list_of_string):
	for string in list_of_string:
		if date in string:
			return string
	return None


def find_mask(date, list_of_mask):
	find_condition = True
	ind = 0
	l_delta = []
	for i in range(200):
		l_delta.append(int(-i))
		l_delta.append(int(i))

	date_of_good_mask = None
	date = datetime.datetime.strptime(date, modeltimefmt)

	while (find_condition and ind < 200):
		test_date = date + datetime.timedelta(days=l_delta[ind])
		test_date = test_date.strftime(modeltimefmt)
		date_of_good_mask = is_date_in_list_of_string(test_date, list_of_mask)
		if date_of_good_mask:
			find_condition = False
		ind += 1
	return date_of_good_mask


######
## Check if passplan.txt exists.
## If not, run select_orbit_cnes.py using parameter_orbit.rdf
print('Checking passplan...')
if not os.path.exists(planpath):
	print('Passplan does not exist. Running parameter_orbit.rdf...')
	print('...Stay tuned...')
	os.system('python /work/ALT/swot/swotpub/modules/install/swot-hydrology-toolbox/select_orbit_cnes/select_orbit_cnes.py %s %s', %(orbit_rdf_path, output_orbit_path))

## Open passplan.txt to read orbit information
cycle = np.genfromtxt(planpath, dtype=int, usecols=(1))
orbit = np.genfromtxt(planpath, dtype=int, usecols=(2))
date = np.genfromtxt(planpath, dtype=str, usecols=(6))
ndates = len(date)

proctime = datetime.datetime.strptime(date[0], '%Y-%m-%d')
modeltime = proctime.strftime(modeltimefmt)

list_of_mask = []
for file in os.listdir(maskpath):
	if file.endswith("shp"):
		list_of_mask.append(file)
list_of_mask.sort()

## Process each shapefile using the CNES simulator
for i in range(ndates):
	## Convert passplan time to datetime tuple
	proctime = datetime.datetime.strptime(date[i], '%Y-%m-%d')
	modeltime = proctime.strftime(modeltimefmt)
	
	mask = find_mask(modeltime, list_of_mask)
	print("===========================================================================================================")
	print("Acquisition date (passplan) : %s ; Corresponding Water Mask : %s" %(modeltime, mask))
	print("===========================================================================================================")

	## Create backup of rdf, named parameter_sisimp.rdf.orig
	if not os.path.exists(sisimp_rdf_path + '.orig'):
		os.system('cp %s %s' %(sisimp_rdf_path, sisimp_rdf_path + '.orig'))
	
	## Then, create modified rdf, named parameter_sisimp.rdf
	lines = open(sisimp_rdf_path).read().splitlines()
	mask = '.'.join(mask.split('.')[:-1])
	for j, line in zip(range(len(lines)),lines):
		if line[0:26] == 'Run directory for orbits =':
			lines[j] = 'Run directory for orbits = ./%s/1_select_orbit' %(testcase)
		elif line[0:16] == 'Shapefile path =':
			lines[j] = 'Shapefile path = %s/%s' %(maskpath, mask)
		elif line[0:18] == 'Output directory =':
			lines[j] = 'Output directory = ./%s/2_sisimp' %(testcase)
		elif line[0:16] == 'Multiple orbit =':
			lines[j] = 'Multiple orbit = no'
		elif "Orbit = " in line:
			lines[j] = 'Orbit = %d' %(orbit[i])
		elif "Cycle number = " in line:
			lines[j] = 'Cycle number = %d' %(cycle[i])
	open(sisimp_rdf_path, 'w').write('\n'.join(lines))

	## Run CNES simulator (SISIMP) using modified parameter_sisimp.rdf
	print('.....Processing with SISIMP.....')
	os.system("python /work/ALT/swot/swotpub/modules/install/swot-hydrology-toolbox/sisimp/proc_sisimp.py %s" %(sisimp_rdf_path))
	
## Run CNES simulator (Lake_Tile)
print('.....Processing with Lake_Tile.....')
os.system("python /work/ALT/swot/swotpub/modules/install/swot-hydrology-toolbox/processing/PGE/lake_tile/multi_lake_tile.py %s" %(lake_tile_cfg_path))

print('\n ALL SHAPEFILES HAVE BEEN PROCESSED ! ')
