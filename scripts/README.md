# l2pixc_to_rivertile
## Purpose
This program converts one or more tiles of pixel cloud into as many as river product(s).

## Run the software
```
usage: l2pixc_to_rivertile.py [-h] [--riverobs_path RIVEROBS_PATH]
                              [--nogdem [NOGDEM]] [--noshp [NOSHP]] [-f]
                              l2pixc_annotation_file output_dir
                              parameter_riverobs

positional arguments:
  l2pixc_annotation_file
                        PixC annotation file
  output_dir            output directory
  parameter_riverobs    param file for RiverObs

optional arguments:
  -h, --help            show this help message and exit
  --riverobs_path RIVEROBS_PATH
                        full path to RiverObs software (if not in os.environ)
  --nogdem [NOGDEM]     if true, don't call riverobs with gdem
  --noshp [NOSHP]       if true, don't produce shapefiles
  -f, --force           force overwrite existing outputs; default is to quit
```
In particular:
* ___l2pixc_annotation_file___ is an output of the large scale simulator (sisimp) or the hr_simulator_to_pixel_cloud.py tool; its pattern is like: ```pixc_annotation_<ccc>_<ppp>_<ttt-t>``` where ___ccc___ is the cycle number over 3 digits, ___ppp___ is the pass number over 3 digits and ___ttt-t___ is the tile reference (ex: 45N-R for the tile of 1° of latitude at nadir, starting at 45° North, Right swath); it is a key-value file, containing at least: ```l2pixc file = <full_path_to_pixc_file>```; __l2pixc_annotation_file__ may be a directory: in this case, the script takes in input all the annotation files in this directory
* ___riverobs_path___ is the full path to RiverObs software; if not filled, the software uses the environment variable named "RIVEROBS"

### Parameter file
The ```parameter_riverobs``` parameter file is in RDF format and must contain the following:
```
width_db_file             (-) = None
shape_file_root           (-) = <full_path_to_river_database_shp_without_ext>
l2_file                   (-) = 'REPLACE_ME'
fout_reach                (-) = 'REPLACE_ME'
fout_node                 (-) = 'REPLACE_ME'
fout_index                (-) = 'REPLACE_ME'
lonmin                    (-) = 'REPLACE_ME'
lonmax                    (-) = 'REPLACE_ME'
latmin                    (-) = 'REPLACE_ME'
latmax                    (-) = 'REPLACE_ME'
bounding_box              (-) = 'lonmin, latmin, lonmax, latmax'
lat_kwd                   (-) = 'latitude'
lon_kwd                   (-) = 'longitude'
class_kwd                 (-) = 'classification'
height_kwd                (-) = 'height'
xtrack_kwd                (-) = 'cross_track'
class_list                (-) = [2, 3, 4, 22, 23, 24]
fractional_inundation_kwd (-) = 'water_frac'
use_fractional_inundation (-) = [True, True, False, False, False, False]
use_segmentation          (-) = [False, True, True, False, True, True]
use_heights               (-) = [False, False, True, False, False, False]
min_points                (-) = 100
clip_buffer               (-) = 20.0
use_width_db              (-) = False
ds                        (-) = 300.0
refine_centerline         (-) = False
smooth                    (-) = 0.01
alpha                     (-) = 1
max_iter                  (-) = 1
scalar_max_width          (-) = 600.0
minobs                    (-) = 10
trim_ends                 (-) = False
fit_types                 (-) = ['OLS', 'WLS', 'RLM']
min_fit_points            (-) = 3
do_improved_geolocation   (-) = True
geolocation_method        (-) = taylor
```
Every parameter can remain the same, except:
* ___shape_file_root___: full path to the shapefile corresponding to the river database; the extension of the file (i.e. ".shp") has to be removed

## Output files
### River output
The river output is a NetCDF file, with 2 groups inside, one dedicated to nodes, the other to reaches.
In this file:
* each ```reaches``` corresponds to a river reach defined in the river database,
* each ```nodes``` corresponds to a river node (ie a point every 200m along the centerline of a reach) defined in the river database.

Its pattern is ```rivertile_<ccc>_<ppp>_<ttt-t>.nc``` where:
* ___ccc___ is the cycle number, on 3 digits
* ___ppp___ is the pass number, on 3 digits
* ___ttt-t___ is the tile reference; for example: 45N-R for the tile of 1° of latitude at nadir, starting at 45° North, Right swath

The river output may be produced in shapefile format, except if the ```noshp``` parameter is set.

### PIXCVec file
The PIXCVec file is a NetCDF file containing the improved geolocation as well as link to river a priori database for pxeils of pixel cloud processed as river pixels.

Its pattern is ```pixcvec_<ccc>_<ppp>_<ttt-t>.nc``` where:
* ___ccc___ is the cycle number, on 3 digits
* ___ppp___ is the pass number, on 3 digits
* ___ttt-t___ is the tile reference; for example: 45N-R for the tile of 1° of latitude at nadir, starting at 45° North, Right swath

The PIXCVec file may be produced in shapefile format, except if the ```noshp``` parameter is set.

### River annotation file
The river annotation file is a working text file, which will be used later in rivertile_to_laketile. Its pattern is ```river-annotation_<ccc>_<ppp>_<ttt-t>.nc``` where:
* ___ccc___ is the cycle number, on 3 digits
* ___ppp___ is the pass number, on 3 digits
* ___ttt-t___ is the tile reference; for example: 45N-R for the tile of 1° of latitude at nadir, starting at 45° North, Right swath

Its content is like:
```
pixc file = <full_path_to_pixc_file>
gdem pixc file = None
pixcvec file = <full_path_to_pixcvec_file>
```

# rivertile_to_laketile
## Purpose
This program computes a LakeTile product from a PIXC (pixel cloud) file (output from large scale simulator sisimp or hr_simulator_to_pixel_cloud.py program) and its associated PIXCVecRiver (output from l2pixc_to_rivertile.py program) on the same tile.

## Run the software
```
usage: rivertile_to_laketile.py [-h] [--parameter_laketile PARAMETER_LAKETILE]
                                [--nogdem [NOGDEM]]
                                river_annotation_file output_dir

positional arguments:
  river_annotation_file
                        river annotation file (output from
                        l2pixc_to_rivertile) or directory
  output_dir            output directory

optional arguments:
  -h, --help            show this help message and exit
  --parameter_laketile PARAMETER_LAKETILE
                        param file, if variables diff from default
  --nogdem [NOGDEM]     if true, don't call riverobs with gdem
```
In particular:
* ___river_annotation_file___ is an output of the rivertile_to_laketile.py program; its pattern is like: ```river-annotation_<ccc>_<ppp>_<ttt-t>``` where ___ccc___ is the cycle number over 3 digits, ___ppp___ is the pass number over 3 digits and ___ttt-t___ is the tile reference (ex: 45N-R for the tile of 1° of latitude at nadir, starting at 45° North, Right swath); it is a key-value file, containing: 

```
pixc file = <full_path_to_pixc_file>
gdem pixc file = None
pixcvec file = <full_path_to_pixcvec_file>
```
It may also be a directory containing several river annotation files.


### Parameter file
The ```parameter_laketile``` parameter file is in ConfigParser format and must contain the following:
```
[PATHS]
PIXC file = REPLACE_ME
PIXCVecRiver file = REPLACE_ME
Output directory = REPLACE_ME

[OPTIONS]
# To also produce LakeTile_edge and LakeTile_pixcvec as shapefiles (=True); else=False (default)
Produce shp = True
# Verbose level; =DEBUG or =INFO
Verbose level = DEBUG
# Generate log file in the output directory (=True) or print on screen (=False; default)
Log file = False

[CONFIG_OVERWRITE]

# Lake a priori database
LAKE_DB = /work/ALT/swot/swotpub/BD/BD_lakes/20170905_France/priordb_lakes_france.shp
# Lake identifier attribute name in the database
LAKE_DB_ID = basin_id
# Shapefile with polygons of continents
CONTINENT_FILE = /work/ALT/swot/swotpub/BD/major_basins/FAO/major_hydrobasins.shp

# Water flag: 3=water near land edge  4=interior water
FLAG_WATER = 3;4
# Dark water flag: 23=darkwater near land  24=interior dark water
FLAG_DARK = 23;24
# Layover flag
FLAG_LAYOVER = 12;13;14

# Min size for a lake to generate a lake product (=polygone + attributes) for it
MIN_SIZE = 1.0  # In ha

# To improve PixC golocation (=True) or not (=False)
IMP_GEOLOC = True

# Method to compute lake boundary or polygon hull
# 0=convex hull 1=concav hull (1.0=with alpha param (default) 1.1=without) 2=concav hull radar vectorisation
HULL_METHOD = 1.0

# Maximal standard deviation of height inside a lake
STD_HEIGHT_MAX = 10

# Big lakes parameters for improved geoloc; used only if imp geoloc=1
BIGLAKE_MODEL = polynomial
BIGLAKE_MIN_SIZE = 5000
BIGLAKE_GRID_SPACING = 4000
BIGLAKE_GRID_RES = 8000

```
[CONFIG_OVERWRITE] section is optional and may not contain all these variables (or they can be commented). Therefore, the default values will be used (see [LOCNES README file](https://gitlab.cnes.fr/swotdev/swot-hydro-toolbox/tree/develop/processing/src/cnes/sas/lake_tile) for more info).

## Output files
See dedicated [LOCNES README file](https://gitlab.cnes.fr/swotdev/swot-hydro-toolbox/tree/develop/processing/src/cnes/sas/lake_tile).