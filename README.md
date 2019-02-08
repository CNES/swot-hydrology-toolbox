# SWOT Hydrology Toolbox.

Copyright (C) 2018 Centre National d'Etudes Spatiales

This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

## Background 

SWOT (Surface Water and Ocean Topography) is an innovative radar altimetry satellite mission projected for launch in 2021, prepared jointly by NASA’s Jet Propulsion Laboratory (JPL) and Centre National d’Etudes Spatiales (CNES), with contributions from the UK Space Agency (UKSA) and the Canadian Space Agency (CSA). SWOT features a high-rate (HR) data mode for continental hydrology, and a low-rate (LR) data mode dedicated mainly to oceanography. For more information, refer to https://swot.cnes.fr/en/ or https://swot.jpl.nasa.gov/. 

## Objectives 
* Provide open-source tools that, together with JPL’s RiverObs tool (https://github.com/SWOTAlgorithms/RiverObs.git), enable end-users to generate virtually all SWOT HR level-2 (L2) products with fairly (but not fully) representative characteristics (see section on caveats below)
  * Get familiar with product content and formats, use the data to conduct studies...
* Give end-users access to the L2+ processing prototypes 
  * Validate methodology, propose improvements...
* As far as possible let the toolbox correspond directly to the processing prototypes that are evolving towards operational processing chains 
  * Coded in high-level programming language (Python 3), for simple distribution, installation and use

```
Note that both algorithms and products are still under development and will be updated regularly.
```

## Content 
* SISIMP: Large scale simulator of L2_HR_PIXC products (with orbit selection tool)
* LOCNES: Generation of L2_HR_LAKE* and L2_HR_PIXC_VEC products 
* Improved geolocation library (used by RiverObs and LOCNES)
* Module to generate L2_HR_RASTER products (under development, not yet included)
* Overall script permitting to run all steps consecutively (with example dataset)
* Tools for file format conversion etc.
* Potentially other modules in the future (for ex. floodplain DEM generation)

## Associated RiverObs version
The RiverObs version compatible with the master branch of the swot-hydrology-toolbox is :

* commit 10199d47cce881f70bd46178a05d132d3a0534ab
* Merge: 8b9a57e ac1c0e4
* Author: Alex Fore <Alexander.Fore@jpl.nasa.gov>
* Date:   Wed Jan 16 07:30:51 2019 -0800


## Caveats
Although the large-scale simulator included in the toolbox provides fairly representative statistical errors, several simplifications or approximations are made:
* Spherical Earth approximate geolocation equations (loss of accuracy at high latitudes, >60°)
* No topography is taken into account 
  * Radar geometry grid constructed on sphere
  * No geometric distortions, no layover 
* Simplified representation of water height (several options) 
  * Spatially constant height for each water body (but possibility to vary height over time, cycle)
  * Spatially random correlated heights, and 2D polynomial model (with synthetic slopes)
  * Also possible to inject “true” heights from models (after simulation)
  * Random effective instrument noise added to height (and propagated to geolocation)
* Idealized pixel cloud processing (i.e. water detection, phase unwrapping etc. are implicitly assumed to be perfect)

```
If a higher degree of realism is necessary to conduct a study, lower-level simulators and processors need to be employed. 
These are not publicly available, but SWOT science team members can contact the SWOT algorithm development team for support. 
```

Products formats and algorithms:
* The product formats are not yet in their final version, but are regularly updated to reflect the current working definitions, and will evolve into the official product formats as and when these are defined. Some data fields are at this stage void (various flags, uncertainty indicators…).
* The processing algorithms will also continue to evolve, as today's prototypes are progessively refined into operational software. 

## Installation procedure

### Get the repository
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

The repository __RiverObs__ should be assignated to the __RIVEROBS__ variable.

```bash
% export RIVEROBS=your_installation_path/RiverObs
```

### Dependencies

The dependencies of swot-hydrology-toolbox are:
* GDAL
* netcdf
* proj4
* libspatialindex
* and the following Python modules:
  * numpy
  * scipy
  * matplotlib
  * scikit-learn
  * scikit-image
  * lxml
  * netCDF4
  * xarray
  * dask
  * distributed
  * pyproj
  * jupyter
  * notebook
  * statsmodels
  * pysal
  * pandas
  * pytables
  * Shapely
  * Fiona
  * sphinx
  * numpydoc
  * rtree
  * mahotas
  * utm

### Python environment installation

#### Setting up a conda environment

To create a conda environment, execute

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

An example dataset showing how to configure and run simulation and processing is available under /test.

The needed input for the overall chain include:
* An orbit file (provided)
* A water mask in shapefile format covering the area you want so simulate (see example under /test)
* A river database in shapefile format (e.g. GRWL)
* A lake database in shapefile format
* Various configuration files (examples provided)
