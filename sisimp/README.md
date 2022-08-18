# Large scale simulator of pixel cloud

This module simulates pixel cloud with realistic errors from a shapefile of water mask, and over large scale.

## Input files

### Shapefile of water mask
This is a polygon shapefile. Each object corresponds to a water body. Specific attributes are required depending on the configuration of the simulator:
*  An attribute named ```HEIGHT``` (or other to fill in the parameter file, see hereafter) gives the height associated to each water body. If the value is set to NULL, then the height is considered to be 0.

### Orbit files
They are the result of the select_orbit_cnes module, i.e:
* ```<prefix>_cycle_<cccc>_pass_<pppp>.nc``` (cccc = cycle number on 4 digits; pppp = pass number on 4 digits) corresponds to a subset of the theoretical orbit over the studied area
* ```passplan.txt``` is the pass plan, i.e. the sequence of passes overflying the studied area; ```passplan_t<pppp>.txt``` contains the information for pass ```pppp```.


## Output files
### The pixel cloud PIXC file
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

### The PIXCVecRiver file (optionnal)

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

# Usage
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

## Parameter file

### With min options
The parameter file is in RDF format and must contain at least the following (```parameter_sisimp_light.rdf```):
```
!### PATHS
Run directory for orbits = <full_path_to_orbit_files>
Shapefile path = <full_path_to_shp_without_.shp>
Output directory = <full_path_to_output_directory>

!### OPTIONAL FILES IN OUTPUT
Create shapefile = <yes|no>
Create dummy PIXCVecRiver file = <yes|no>


!###################
!# ORBIT SELECTION #
!###################

Multiple orbit = <no/yes/passplan>
Orbit = <cycle_number> 
```
In particular:
* ___Orbit___ is needed only if ```Multiple orbit = no```
* ___Create shapefile___ is to generate PixC files also in shapefile format
* ___Create dummy pixc vec river file___ is to generate arbitrary PIXCVecRiver file (normally processed by RiverObs; useful when you only simulate lakes and you don't want to run RiverObs)
* <span style="color: red;">__Please, note that only instrument error is applied with this light configuration__</span> (cf. below for dark water error, wet tropo error, and cross-over residual roll error)

### All options
The full parameter file can contain the following attributes (```parameter_sisimp_full.rdf```):
```
!### PATHS
Run directory for orbits = <full_path_to_orbit_files>
Shapefile path = <full_path_to_shp_without_.shp>
Output directory = <full_path_to_output_directory>

!### AUXILIARY DATA
Geoid = <full_path_to_geoid_map>

!### PIXC TILING (must be coherent with orbit used in select_orbit_cnes processing)
Tile database path = <full_path_to_tile_files>
Longitude delay = <value_in_seconds>

!### OPTIONAL FILES IN OUTPUT
Create shapefile = <yes|no>
Create dummy PIXCVecRiver file = <yes|no>


!###################
!# ORBIT SELECTION #
!###################

!== 3 options: Multiple orbit = yes/no/passplan (default=yes)

!### Option 1 (default) - All orbit files in Orbit directory will be processed
Multiple orbit = yes

!### Option 2 - Set the Orbit to a correct number found in the Orbit directory; only this orbit file will be processed
!Multiple orbit = no
!Orbit = 29
!Cycle number = 1  !Optional; default=1

!### Option 3 - Orbit files will be processed according to a passplan file 
!Multiple orbit = passplan
!Passplan path = <full_path_to_passplan>  !Optional; default=passplan generated by previous select_orbit_cnes module


!#######################
!# HEIGHT MODELISATION #
!#######################

!== 5 options: Height model = None/polynomial/gaussian/reference_height/reference_file (default=None => constant height model)

!### Option 1 (default) - Constant height model, uniform over each water body; height varies sinusoidally with time
!Specific option: [Height model A = 0.] & [!Height model] => no height applied, height in output only contains errors
!### Constant height model parameters (same height for each water body) 
Constant height model A = 10.     !=0 to disable
Constant height model t0 = 47076
Constant height model period (days) = 365.25
!### Complex 2D height model parameters (2D variations added for lakes > [Height model min area]) 
!Height model = polynomial        !=polynomial or gaussian; if disabled, only constant height model
!Height 2d model min area = 100.  !(ha) min area of water bodies on which to add complex 2D height model (default=100.)
!Height 2d model stdv = 1.        !IF Height model = gaussian: stdv for gaussian model 

!### Option 2 - Height is given from a specific attribute in the shapefile of water bodies (height is a single value per water body)
!Height model = reference_height
!Height shp attribute name = ref_height  !Name of the attribute (default=HEIGHT)
!Height ref multitemp = True             !IF True, set constant height model parameters (see option 1 above)

!### Option 3 - Height is given in a NetCDF file (height is a 2D map over each water body)
!Height model = reference_file
!True height file = <full_path_to_NetCDF_height_file>


!#########################
!# SIMULATION PARAMETERS #
!#########################

!== Instrument caracteristics ==
Swath width = 120000.
NR cross track = 5000.
NR cross track filter = False  !=True for select only PixC beyond [NR cross track] distance
Sensor wavelength = 0.008385803
Baseline = 10.
Range sampling = 0.75
Number of pixels in range = 4575

!== Orbit ==
Orbit jitter = 1000

!== Classification flags ==
Land flag = 1
Land water flag = 2
Water land flag = 3
Water flag = 4
Dark water flag = 5


###########################
!# NOISE AND ERROR CONFIG #
###########################

!== Noise parameters ==
Noise file path = $SWOT_HYDROLOGY_TOOLBOX/sisimp/data/height_noise_presum2.txt
Noise file path for land= $SWOT_HYDROLOGY_TOOLBOX/sisimp/data/height_noise_presum2_land.txt
Noise file path for dark water= $SWOT_HYDROLOGY_TOOLBOX/sisimp/data/height_noise_presum2_land.txt
Noise multiplier factor = 0.5     !1/sqrt(Nl) where Nl is 4
Height bias std = 0.  !Deprecated
Geolocalisation improvement = no  !yes (no cross-track noise applied) or no (default)

!== Dark water ==
Dark water = yes                     !If yes: dark water simulation is executed (default=no)
Dark water percentage = 10           !int (between 0 and 100) probability of total dark water simulated
Dark water detected percentage = 90  !int (between 0 and 100) percentage of detected dark water
Scale factor non detected dw = 0.5   !float to parameterize the correlation length of simulated non detected dark water areas
Dark water correlation length = 50
!Dark water seed = 12345678          !int seed to be used for reproducible dark water simulations

!== Tropo error model ==
!== 3 options: Tropo model = None/gaussian/map (default=None)
!Tropo model = gaussian
!Tropo error correlation = 5000
!Tropo error stdv = 0.03
!Tropo error mean = 0.01

!////////// ON CNES CLUSTER ONLY //////////

!Tropo model = map
!Tropo error correlation = 5000
!Tropo error map file = <full_path_to_wet_tropo_map>

!== Cross-over residual roll error ==
!roll_repository_name = <full_path_to_roll_dir>

!//////////////////////////////////////////
```
In particular:
* There are 3 options to set the height of the water bodies (cf. above)
* ___Dark water___ can be set to simulate dark water over portions of water bodies (cf. above)
* __Tropo residual error model__ can be set to include wet tropo error 
* <span style="color: red;">__Please, note that tropo error map and roll error can only be applied when run on CNES machines__</span>
