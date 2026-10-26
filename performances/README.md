# SWOT LAKE performance assessment computations

Software and data are distributed WITHOUT ANY WARRANTY.
Refer to LICENSE.txt under https://github.com/CNES/swot-hydrology-toolbox/.

## Objective 
Provide the scripts and data used to generate the results of the SWOT missions global evaluation of lake elevation and area accuracy based on L2_HR_LakeSP_Prior products. Several filtering combinations are tested, based on flags embedded in the products (partial_f, ice_f, quality_f). 

## Directory structure and content 
* **perfos_lakesp_with_filters.ipynb**: the notebook used to calculate the statistics for absolute elevation errors and relative surface area errors for lakes, and the input data.
* **input**: directory with CSV files containing the differences between in situ and SWOT observations: `diff_area.csv` for area comparisons and `diff_wse.csv` for water surface elevation (WSE) comparisons over lakes.

## Associated SHT version
develop branch

## Installation procedure

### Get the repository
1. Clone only __performances__ directory from __swot_hydrology_toolbox__ repository

```
% git clone --filter=blob:none --no-checkout --branch develop https://github.com/CNES/swot-hydrology-toolbox.git
% cd swot-hydrology-toolbox
% git sparse-checkout init --cone
% git sparse-checkout set performances
% git checkout
```
2. The repository __swot_hydrology_toolbox__ should be assignated to the __SWOT_HYDROLOGY_TOOLBOX__ variable.

```bash
% export SWOT_HYDROLOGY_TOOLBOX=your_installation_path/swot-hydrology-toolbox
```

### Dependencies
The dependencies of performances are:
- pandas
- matplotlib
- numpy
- notebook
- ipykernel
- shapely
- geopandas

### Example of environment.yml
name: base
channels:
  - defaults
  - conda-forge
dependencies:
  - python=3.13
  - conda=26.1.1
  - setuptools
  - pip
  - wheel
  - shapely
  - geopandas
  - ipykernel
  - jupyterlab
  - conda-anaconda-telemetry[version='>=0.3.0']
  - conda-anaconda-tos[version='>=0.2.2']
  - conda-content-trust
  - conda-libmamba-solver
  - menuinst
  - anaconda-anon-usage
  - anaconda-auth
  
### Python environment installation

#### Creating a conda environment (default method)

To create a conda environment, execute:

```
cd $SWOT_HYDROLOGY_TOOLBOX/performances
conda create -n perfos_env python=3.10 pandas matplotlib numpy shapely geopandas jupyterlab notebook ipykernel -y
conda activate perfos_env
python -c "import pandas, matplotlib, numpy, shapely, geopandas, ipykernel, notebook, jupyterlab; print('OK')"
```

#### Creating a conda environment from a .yaml file

To create a conda environment, execute:

```
cd $SWOT_HYDROLOGY_TOOLBOX/performances
conda env create -f environment.yml -n perfos_env
```

### Python environment activation

To activate this environment, type:
```
conda activate perfos_env
```

To deactivate this environment, type
```
conda deactivate
```

## Setting up Python path
After activating your Python environment, you have to set your PYTHONPATH variables:
```
export PYTHONPATH=$SWOT_HYDROLOGY_TOOLBOX/performances/:$PYTHONPATH
```

## To run a lake performance assessment by using the notebook:
```
jupyter notebook perfos_lakesp_with_filters.ipynb 
```
