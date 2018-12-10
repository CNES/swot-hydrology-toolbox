# select_orbit_cnes
## Purpose
The ```select_orbit_cnes``` software computes orbit files specific to the studied area from a set of theoretical orbit files.

## Run the software
```
usage: select_orbit_cnes.py [-h] param_file output_dir

Compute orbit files specific to the studied area

positional arguments:
  param_file  full path to the parameter file (*.rdf)
  output_dir  full path to the output directory

optional arguments:
  -h, --help  show this help message and exit
```

### Parameter file
The parameter file is in RDF format and must contain the following:
```
!== Mission specific parameters ==
Mission start time = 2014-01-01                                                 
Cycle duration (sec) = 1802645.8059698 


!== Orbit parameters ==

!Directory of theoretical orbit files
Orbit repository = <full_path_to_orbit_files>

!Studied area bounding box
Area south latitude (deg) = <min_lat>
Area north latitude (deg) = <max_lat>
Area west longitude (deg) = <min_long>
Area east longitude (deg) = <max_long>

!PixC mapping parameters
Azimuth spacing (m) = 21.875000 
Swath width (m) = 120000.000000
NR cross track (m) = 10000.000


!== Pass plan parameters ==
passplan = yes
simulation_start_time = 2015-04-01
simulation_stop_time = 2015-05-01


!== Output parameters ==
GDEM Orbit prefix = test.gdem_orbit 
```
In particular:
* ___Orbit repository___ is the full path for the directory storing the theoretical orbit files
* ___Area south|north latitude / west|east longitude___ is the bounding box of the studied area
* ___passplan___ is the flag to generate the pass plan or not (see hereafter)

## Output files
### Orbit files
The orbit file is a NetCDF file,corresponding to a subset of the theoretical orbit over the studied area.

Its pattern is ```<prefix>_cycle_<cccc>_pass_<pppp>.nc``` where:
* ___prefix___ is a string specified in the parameter file (cf. hereafter)
* ___cccc___ is the cycle number, on 4 digits
* ___pppp___ is the pass number, on 4 digits

In this file, each ```record``` corresponds to a point along the nadir track. For each record, the following values are given:
* ```time```
* ```heading```
* geodetic coordinates: ```longitude```, ```latitude```, ```altitude```
* cartesian coordinates: ```x```, ```y```, ```z```

### Pass plan
If ```passplan = yes``` in the parameter file, an ASCII file named ```passplan.txt``` is generated in the output directory. This file contains the sequence of passes overflying the studied area between the ```simulation_start_time``` and the ```simulation_stop_time``` specified in the parameter file.

It looks like:
```
# Mission start:    2014-01-01 00:00:00
# Simulation start: 2015-04-01 00:00:00
# Simulation stop:  2015-05-01 00:00:00
#
#run      cycle orbit MissionTime year DayOfYear       date     time
c021_t573     21  573    39620418 2015  94.56965 2015-04-04 13:40:18
c022_t267     22  267    40478498 2015 104.50113 2015-04-14 12:01:38
c022_t338     22  338    40699259 2015 107.05624 2015-04-17 01:20:59
```
