# Large scale simulator of pixel cloud

## Purpose
This module simulates pixel cloud with realistic errors from a shapefile of water mask, and over large scale.

### Input files

#### Shapefile of water mask
This is a polygon shapefile. Each object corresponds to a water body. Specific attributes are required depending on the configuration of the simulator:
* ```RIV_FLAG = 1``` for river objects, ```=0``` otherwise; to use with ```Create dummy pixc vec river file = yes``` to distinguish river pixels from other pixels; if RIV_FLAG not set, all water bodies are considered as lakes
*  An attribute named ```HEIGHT``` (or other to fill in the parameter file, see hereafter) gives the height associated to each water body. If the value is set to NULL, then the height is considered to be 0.

#### Orbit files
They are the result of the select_orbit_cnes module, i.e:
* ```<prefix>_cycle_<cccc>_pass_<pppp>.nc``` (cccc = cycle number on 4 digits; pppp = pass number on 4 digits) corresponds to a subset of the theoretical orbit over the studied area
* ```passplan.txt``` is the pass plan, i.e. the sequence of passes overflying the studied area; ```passplan_t<pppp>.txt``` contains the information for pass ```pppp```.


### Output files
#### The pixel cloud PIXC file
The filename pattern is ```SWOT_L2_HR_PIXC_I_<ccc>_<ppp>_<ttt><s>_<yyyyMMddThhmmss>_<yyyyMMddThhmmss>_<CRID>_<nn>``` where:
* __ccc__ is the cycle number, on 3 digits
* __ppp__ is the pass number, on 3 digits
* __ttt__ is the tile reference, on 3 digits; for ex: "45N"
* __s__ is the swath: =L for Left swath =R for right swath
* First __yyyyMMddThhmmss__ is the start date of the tile
* Second __yyyyMMddThhmmss__ is the end date of the tile
* __CRID__ is the Composite Release IDentifier (set to Dx0000)
* __nn__ is a product counter with the same CRID

NB: a shapefile version of the PixC is produced if asked

#### The PIXCVecRiver file (optionnal)

This is the file complementary to the PIXC file, containing improved geolocation for river objects (cf. above).

The filename pattern is ```SWOT_L2_HR_PIXCVecRiver_<ccc>_<ppp>_<ttt><s>_<yyyyMMddThhmmss>_<yyyyMMddThhmmss>_<CRID>_<nn>``` where:
* __ccc__ is the cycle number, on 3 digits
* __ppp__ is the pass number, on 3 digits
* __ttt__ is the tile reference, on 3 digits; for ex: "45N"
* __s__ is the swath: =L for Left swath =R for right swath
* First __yyyyMMddThhmmss__ is the start date of the tile
* Second __yyyyMMddThhmmss__ is the end date of the tile
* __CRID__ is the Composite Release IDentifier (set to Dx0000)
* __nn__ is a product counter with the same CRID

NB: a shapefile version of the PixC is produced if asked

## Run simulator

```
usage: proc_sisimp.py [-h] [-v [VERBOSE]] [-l [LOGFILE]] parameter_file

positional arguments:
  parameter_file

optional arguments:
  -h, --help            show this help message and exit
  -v [VERBOSE], --verbose [VERBOSE]
                        Verbose level (DEBUG of INFO=default)
  -l [LOGFILE], --logfile [LOGFILE]
                        Write prints to a logfile in addition to the console
```

### Parameter file

#### With min options
The parameter file is in RDF format and must contain at least the following (```parameter_sisimp_light.rdf```):
```
Run directory for orbits = <full_path_to_orbit_files>
Shapefile path = <full_path_to_shp_without_.shp>
Output directory = <full_path_to_output_directory>

!### Noise and error files 
Noise file path = $SWOT_HYDROLOGY_TOOLBOX/sisimp/data/height_noise_presum2.txt

!### Orbit parameters: 
Multiple orbit = <no|yes|passplan>
Orbit = <cycle_number> 

!### Files in output
Create shapefile = <yes|no>
Create dummy pixc vec river file = <yes|no>
```
In particular:
* ___Orbit___ is needed only if ```Multiple orbit = no```
* ___Create shapefile___ is to generate PixC files also in shapefile format
* ___Create dummy pixc vec river file___ is to generate arbitrary PIXCVecRiver file (normally processed by RiverObs)

#### All options
The full parameter file can contain the following attributes (```parameter_sisimp_full.rdf```):
```
Run directory for orbits = $SWOT_HYDROLOGY_TOOLBOX/test/river_and_lake/output/orbit
Shapefile path = $SWOT_HYDROLOGY_TOOLBOX/test/river_and_lake/data/river_and_lake
Output directory = $SWOT_HYDROLOGY_TOOLBOX/test/river_and_lake/output/simu

!### Noise and error files 
Noise file path = $SWOT_HYDROLOGY_TOOLBOX/sisimp/data/height_noise_presum2.txt
roll_file_name = /work/ALT/swot/swotpub/SWOT_Simulator_data/input_data/large_scale_sim_scene/data_sim_roll_v1.nc

!### Orbit parameters
!3 options =
!Multiple orbit = yes (default) => all orbit files in Orbit directory will be processed
!Multiple orbit = no => set the Orbit to a correct number found in the Orbit directory; only this orbit file will be processed
!Multiple orbit = passplan => orbit files will be processed according to passplan.txt file (generated if passplan = yes in select_orbit.rdf)
Multiple orbit = no !=no/yes/passplan
Orbit = 267

!### Simulation parameters
Swath width = 120000.
NR cross track = 10000.
Sensor wavelength = 0.0084
Range sampling = 0.75 !by default
Number of pixels in range = 3117 !3117 (10-60km) or 3500 (extended swath 5-65km)
Orbit jitter = 1000
Water flag = 4

!### Error parameters
Height bias std = 0.1
Noise multiplier factor = 0.5 !1/sqrt(Nl) where Nl is 4
Geolocalisation improvement = no !yes (no cross-track noise applied) or no

!### Dark water parameters
Dark water = no !yes (dark water simulation is executed) or no
!Scale factor dw = 2.0  !correlation length of simulated dark water areas (float)
!Dark water percentage = 10  !probability of total dark water simulated (int between 0 and 100)
!Dark water flag = 24  !classification flag for detected dark water (usually 24)
!Dark water seed = 1234567891  !int seed to be used for reproducible dark water simulations

!### Files in output
Create shapefile = yes !Produce output files also as shapefiles
Create dummy pixc vec river file = yes !Produce L2_HR_PIXCVecRiver product associated to PixC files


!== Height ==

!### Option 1 - Constant height model, uniform over each water body; height varies sinusoidally with time
!Specific option: [Height model A = 0.] & [!Height model] => no height applied, height in output onlyt contains errors
!### Constant height model parameters (same height for each water body) 
Constant height model A = 0.  !=0 to disable
Constant height model t0 = 47076
Constant height model period (days) = 365.25
!### Complex 2D height model parameters (2D variations added for lakes > [Height model min area]) 
!Height model = gaussian    !=polynomial or gaussian; if disabled, only constant height model
!Height 2d model min area = 100.     !(ha) min area of water bodies on which to add complex 2D height model (default=100.)
!Height 2d model stdv = 0.2        !stdv for gaussian model (ie Height model = gaussian)

!### Option 2 - Height is given from a specific attribute in the shapefile of water bodies
!Height model = reference_height
!Height shp attribute name = ref_height     !Name of the attribute (default=HEIGHT)

!### Option 3 - Height is given in a NetCDF file
!Height model = reference_file
!True height file = $SWOT_HYDROLOGY_TOOLBOX/test/river_and_lake/data/true_height_model_river_and_lake.nc
```
In particular:
* There are 3 options to set the height of the water bodies (cf. above)
* ___Dark water___ can be set to simulate dark water over portions of water bodies (cf. above)