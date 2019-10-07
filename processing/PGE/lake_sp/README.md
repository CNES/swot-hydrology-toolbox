# LakeSP processor

## Purpose
The LakeSP processor computes a lake single-pass product from tiles of LakeTile files related to specific cycle, pass and continent:

![LakeSP overview](overview_lake_sp.png)

### Input files
The LakeTile pattern is ```SWOT_L2_HR_LakeTile_<ccc>_<ppp>_<ttt><s>_<yyyyMMddThhmmss>_<yyyyMMddThhmmss>_<CRID>_<nn>``` where:
* __ccc__ is the cycle number, on 3 digits
* __ppp__ is the pass number, on 3 digits
* __ttt__ is the tile number, on 3 digits
* __s__ is the swath: =L for Left swath =R for right swath
* First __yyyyMMddThhmmss__ is the start date of the tile
* Second __yyyyMMddThhmmss__ is the end date of the tile
* __CRID__ is the Composite Release IDentifier
* __nn__ is a product counter with the same CRID

The LakeTile product is composed of 3 files:
* "[pattern].shp" contains the lake product for each lake entirely covered by the swath
* "[pattern]_pixcvec.nc" contains the PIXCVec product (ie. improved geolocation and link to a priori data for each pixel already processed by the RiverTile processor and pixels of lake entirely covered by the swath)
* "[pattern]_edge.nc" contains the subset of PIXC file for pixels belonging to lakes cut at the top or bottom of the tile edge

### Output files

#### The lake single-pass file
The LakeSP pattern is ```SWOT_L2_HR_LakeSP_<ccc>_<ppp>_<CC>_<yyyyMMddThhmmss>_<yyyyMMddThhmmss>_<CRID>_<nn>``` where:
* __ccc__ is the cycle number, on 3 digits
* __ppp__ is the pass number, on 3 digits
* __CC__ is the continent reference, on 2 letters; for ex: "EU" for Europe; the continental splitting has been retrieved from HydroBASINS splitting
* First __yyyyMMddThhmmss__ is the start date of the tile
* Second __yyyyMMddThhmmss__ is the end date of the tile
* __CRID__ is the Composite Release IDentifier
* __nn__ is a product counter with the same CRID

#### The PIXCVec file
This is the file complementary to the PIXC file, containing improved geolocation and link to a priori data for each pixel.
The PIXCVec pattern is ```SWOT_L2_HR_PIXCVec_<ccc>_<ppp>_<ttt><s>_<yyyyMMddThhmmss>_<yyyyMMddThhmmss>_<CRID>_<nn>``` where:
* __ccc__ is the cycle number, on 3 digits
* __ppp__ is the pass number, on 3 digits
* __ttt__ is the tile number, on 3 digits
* __s__ is the swath: =L for Left swath =R for right swath
* First __yyyyMMddThhmmss__ is the start date of the tile
* Second __yyyyMMddThhmmss__ is the end date of the tile
* __CRID__ is the Composite Release IDentifier
* __nn__ is a product counter with the same CRID

## Run the processor
```
usage: pge_lake_sp.py [-h] [-shp] [-l] [-v {DEBUG,INFO}]
                      indir_or_param_file [output_dir] [cycle_num] [pass_num]

Compute SWOT LakeSP product from LakeTile products corresponding to a specific (cycle, pass). If indir_or_param_file is a parameter file (*.cfg), all input parameters are only read in the parameter file.

positional arguments:
  indir_or_param_file   LakeTile directory or parameter file (*.cfg)
  output_dir            output directory
  cycle_num             cycle number
  pass_num              pass number

optional arguments:
  -h, --help            show this help message and exit
  -shp                  convert output NetCDF file as shapefile
  -l, --logfile         write prints to a logfile
  -v {DEBUG,INFO}, --verbose {DEBUG,INFO}     
  						verbose level
```

### Option 1 = in-line command
The first way to run the LakeSP processor is through the following in-line command:
```
python pge_lake_sp.py [-shp] [-l] [-v {DEBUG,INFO}] input_dir output_dir cycle_num pass_num
```
where:
* __input_dir__ *(mandatory)* is the directory containing LakeTile files
* __output_dir__ *(mandatory)* is the output directory
* __cycle_num__ *(mandatory)* is the cycle number
* __pass_num__ *(mandatory)* is the pass number
* __-shp__ *(optional)*: if set, PIXCVec file (cf. above) is produced not only in NetCDF (nominal) but also in shapefile format (optional)
* __-l__ *(optional)*: if set, logs are printed in a dedicated file named ```pge_lake_sp_<run_date_in_yyyyMMdd-hhmmss>.log``` located in the output directory
* __-v__ *(optional)* is the verbose level; it may be "DEBUG" or "INFO" (default); if not set, the default value is used

### Option 2 = parameter file
The second way to run the LakeSP processor is through a parameter file:
```
python pge_lake_sp.py [-l] [-v {DEBUG,INFO}] param_file
```
where:
* __param_file__ *(mandatory)* is the full path of the parameter file; the file extension has to be *.cfg (Python ConfigParser format; see [here](https://docs.python.org/3/library/configparser.html) for more details)
* __-l__ *(optional)*: if set, logs are printed in a dedicated file named ```pge_lake_sp_<run_date_in_yyyyMMdd-hhmmss>.log``` located in the output directory
* __-v__ *(optional)* is the verbose level; it may be "DEBUG" or "INFO" (default); if not set, the default value is used

NB: if shp parameter is set as in the in-line command list, its values is automatically replaced by the one in the parameter file.

The parameter file must contain the following:
```
[PATHS]
LakeTile directory = <input_directory>
Output directory = <output_directory>

[DATABASES]
# OPTIION 1 : SQLITE lake database containing  : lake_table, lake_influence_table, basin_table
LAKE_DB = /work/ALT/swot/swotpub/BD/BD_lakes/20190624_EU/EU_lakedb.sqlite

# OPTION 2 : SHP lake database
# Prior lake database
# LAKE_DB = /work/ALT/swot/swotpub/BD/BD_lakes/20190624_EU/EU_lakedb.shp
# Lake identifier attribute name in the prior lake database and influence_lake_db
# LAKE_DB_ID = lakedb_id

[TILES_INFOS]
# Format = int; if empty, deal with all LakeTile files in LakeTile directory
Cycle number = 0
# Format = int; if empty, deal with all LakeTile files of cycle "Cycle" in LakeTile directory
Pass number = 17

[OPTIONS]
# To also produce LakeTile_edge and LakeTile_pixcvec as shapefiles (=True); else=False (default)
Produce shp = <True|False>

[LOGGING]
# Log filename
logFile = <full_path_to_log_file>
# Log level put inside the file
logfilelevel = <DEBUG|INFO>
# Is log console output ?
logConsole = <True|False>
# Log level print in console
logconsolelevel = <DEBUG|INFO>

```

## Configuration parameters
The configuration parameters used by LakeSP are retrieved from [LakeTile].shp.xml metadata file. If they differ between one tile and another, the software stops with an error.


## Multi-tile processing
```
usage: multi_lake_sp.py [-h] [-cycle CYCLE] [-pass PASS] [-shp] [-l] [-v {DEBUG,INFO}]
                        indir_or_param_file [output_dir]

Compute multiple SWOT LakeSP products from LakeTile products corresponding to
one or more specific (cycle, pass). If indir_or_param_file is a parameter file
(*.cfg), all input parameters are only read in the parameter file.

positional arguments:
  indir_or_param_file   LakeTile directory or parameter file (*.cfg)
  output_dir            output directory

optional arguments:
  -h, --help            show this help message and exit
  -cycle CYCLE          cycle number
  -pass PASS            pass number
  -shp                  convert output NetCDF file as shapefile
  -l, --logfile         write prints to a logfile
  -v {DEBUG,INFO}, --verbose {DEBUG,INFO}
                        verbose level
```

The parameter file must contain the following:
```
[PATHS]
LakeTile directory = <input_directory>
Output directory = <output_directory>

[TILES_INFOS]
# Format = int; if empty, deal with all LakeTile files in LakeTile directory
#Cycle number = 0
# Format = int; if empty, deal with all LakeTile files of cycle "Cycle" in LakeTile directory
#Pass number = 17

[OPTIONS]
# To also produce LakeTile_edge and LakeTile_pixcvec as shapefiles (=True); else=False (default)
Produce shp = <True|False>
```

NB: if not used, parameters HAVE TO be removed or in comment (#)

## Algorithm main steps

![Alt text](workflowGitlab_lake_sp.png?raw=true "Workflow diagram")

