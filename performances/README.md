# SWOT LAKE performance assessment computations

Software and data are distributed WITHOUT ANY WARRANTY.
Refer to LICENSE.txt under https://github.com/CNES/swot-hydrology-toolbox/.

## Objective 
Provide the scripts and data used to generate the results of the SWOT missions global evaluation of lake elevation and area accuracy based on L2_HR_LakeSP_Prior products. Several filtering combinations are tested, based on flags embedded in the products (partial_f, ice_f, quality_f). 

## Directory structure and content 
* **examples**: the notebook used to calculate the statistics for absolute elevation errors and relative surface area errors for lakes, and the input data
* **src**:  libraries and the main Python script called from the notebook.

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
* pandas
* matplotlib
* numpy
* notebook
* shapely
* geopandas


### Example of environment.yml
name: base
channels:
  - defaults
dependencies:
  - python=3.13
  - conda=26.1.1
  - setuptools
  - pip
  - wheel
  - shapely
  - geopandas
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
export PYTHONPATH=$SWOT_HYDROLOGY_TOOLBOX/performances/src/:$PYTHONPATH
```

## To run a lake performance assessment

### By using the notebook:
```
jupyter notebook examples/perfos_lakesp_with_filters.ipynb 
```
### (deprecated) By executing python code:

#### Import libraries
```
import sys
from datetime import date
import os
import pandas as pd
import src.stats_lakesp as my_perfos
```
#### Update Python path
```
swot_path = os.environ.get('SWOT_HYDROLOGY_TOOLBOX')
sys.path.append(os.path.join(swot_path, "performances"))
```
#### Define input/output directories
```
input_dir = os.path.join(swot_path, "performances/examples/input")
output_dir = os.path.join(swot_path, "performances/examples/output")
out_dir = os.path.join(output_dir, "lakes_perfos_area_20260217")
```
#### Compute the 1-sigma error statistics, example for the area_total variable:

##### Load input data
```
df_diff = pd.read_csv(os.path.join(input_dir, "diff_area.csv"), sep=";")
```
##### Filter the data to compare the SWOT data with reference data from Pléiades and RCM water masks
```
df_filtered = df_diff[df_diff['_fuid'].str.startswith(('SLUM_TP:', 'RCM_TP:'))]
```
##### Filter the data based on the quality indicator partial_f 
```
df_subset = df_filtered[df_filtered["partial_f"] == 0]
```
##### Run the statistical computation: The algorithm's performance is compared to the targets set by the scientific specifications ('goal'): 
```
suffix_dir = "M_partial_f_eq0"
my_perfos.plot_stats(os.path.join(out_dir, suffix_dir), df_subset, which_stats="area_total", which_perfos="goal")
```
##### The statistics images are saved to disk, under:
```
Statistics images are saved to the disk: $out_dir/M_partial_f_eq0
```


