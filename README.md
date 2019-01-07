# Installation procedure

## Get the repository
1. Clone __swot_hydrology_toolbox__ repo

```bash
% git clone https://github.com/cnes/swot-hydrology-toolbox.git
```
The repository __swot_hydrology_toolbox__ should be assignated to the __SWOT_HYDROLOGY_TOOLBOX__ variable.

```bash
% export SWOT_HYDROLOGY_TOOLBOX=your_installation_path/swot-hydrology-toolbox
```

2. Clone __RiverObs__ repo

```bash
% git clone https://github.com/SWOTAlgorithms/RiverObs.git
```

The repository __RiverObs__ be assignated to the __RIVEROBS__ variable.

```bash
% export RIVEROBS=your_installation_path/RiverObs
```

## Dependancies

The dependencies of swot-hydrology-toolbox are:
* GDAL
* netcdf
* proj4
* libspatialindex
* and the following Python modules:
  o numpy
  o scipy
  o matplotlib
  o scikit-learn
  o scikit-image
  o lxml
  o netCDF4
  o xarray
  o dask
  o distributed
  o pyproj
  o jupyter
  o notebook
  o statsmodels
  o pysal
  o pandas
  o pytables
  o Shapely
  o Fiona
  o sphinx
  o numpydoc
  o rtree
  o mahotas
  o utm

## Python environment installation

### Setting up an conda environment

To create an conda environment, execute :

```bash
cd $SWOT_HYDROLOGY_TOOLBOX
conda env create -f environment.yml
```

To activate this environment, if the first option was used, type
```bash
conda activate swot-env
```

To deactivate this environment, type
```bash
conda deactivate
```

## Execute the toolbox

After activating your Python environment, you have to set your PYTHONPATH variables:
```bash
export PYTHONPATH=$SWOT_HYDROLOGY_TOOLBOX/processing/src/:$RIVEROBS/src:$PYTHONPATH
```


