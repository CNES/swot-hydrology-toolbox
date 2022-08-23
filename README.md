# SWOT Hydrology Toolbox.

Copyright (C) 2018-2020 Centre National d'Etudes Spatiales

This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

## Background 

SWOT (Surface Water and Ocean Topography) is an innovative radar altimetry satellite mission projected for launch in 2021, prepared jointly by NASA’s Jet Propulsion Laboratory (JPL) and Centre National d’Etudes Spatiales (CNES), with contributions from the UK Space Agency (UKSA) and the Canadian Space Agency (CSA). SWOT features a high-rate (HR) data mode for continental hydrology, and a low-rate (LR) data mode dedicated mainly to oceanography. For more information, refer to https://swot.cnes.fr/en/ or https://swot.jpl.nasa.gov/. 

## Objectives 
* Provide open-source tools that, together with JPL’s RiverObs tool (https://github.com/SWOTAlgorithms/RiverObs.git), enable end-users to generate virtually all SWOT HR level-2 (L2) products with fairly (but not fully) representative characteristics (see section on caveats below)
  * Get familiar with product content and formats, use the data to conduct studies...
* Give end-users access to the L2+ HR processing prototypes 
  * Validate methodology, propose improvements...
* As far as possible let the toolbox correspond directly to the processing prototypes that are evolving towards operational processing chains 
  * Coded in high-level programming language (Python 3), for simple distribution, installation and use

```
Note that both algorithms and products are still under development and will be updated regularly.
```

## Content 
* SISIMP: Large scale simulator of L2_HR_PIXC products (with orbit selection tool)
* LOCNES: Generation of L2_HR_LakeTile, L2_HR_LakeSP and L2_HR_PIXCVec products 
* Improved geolocation library (used by RiverObs and LOCNES)
* Module to generate L2_HR_Raster products (under development, not yet included)
* Flooplain DEM prototype
* Overall script permitting to run all steps consecutively (with example dataset)
* Tools for file format conversion etc.
* Potentially other modules in the future 

## Associated SHT version
release_version_08_23_2022 branch

## Associated RiverObs version
develop branch

commit 3ddd70202765298f6ddf7091edfb75c0a7b5ff1c 
Author: Alex Fore <Alexander.Fore@jpl.nasa.gov>
Date:   Wed May 18 05:49:40 2022 -0700

    add validate_inputs to rivertile SAS and check for pixc with no pixels


updated batch validation tool for compatibility with plot reach

River database SWORD v12 available here:
http://gaia.geosci.unc.edu/SWORD/SWORD_v08.zip

Don't forget to modify parameter_river.rdf
reach_db_path (-) = /work/ALT/swot/swotpub/BD/BD_rivers/SWORD_v08/Reaches_Nodes/netcdf

You can also try to use a more recent RiverObs version, but don't forget to use the associated SWORD version.

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
* Idealized pixel cloud processing 
* Synthetic "dark water" model (correlated random fields used to simulate low reflectivity areas)
* Geoid (mean tide corrected EGM-2008), tropospheric and cross-over residual errors simulated

When to use the simulator:
* To familiarize with SWOT HR products, if needed over large areas and over time (multitemporal series)
* To study the inpact of the geometrical shapes of the water bodies on River and Lake processing 
* When a simplified representation of phenomenology, hydrological characteristics and errors is sufficient (e.g. no layover, artifical water slope, basic error models)  

```
If a higher degree of realism is necessary to conduct a study, lower-level simulators and processors need to be employed. 
These are not publicly available, but SWOT Science Team members can contact the SWOT Algorithm Development Team for support. 
```

Product formats and algorithms:
* The product formats correspond to the current official versions, but are likely to evolve. Some data fields are at this stage void (various flags, some uncertainty indicators…).
* The processing algorithms will also continue to evolve, as today's prototypes are progessively refined into operational software. 

Last modifications:
In the large scale simulator:
* Uncertainties of geolocated heights added
* Multilook adaptive averaging implemented
* Land pixels around water bodies added (label 1 and 2 in classification field)
* Near-range computaton improved
* Some fields added or made more realistic

In the processing chain
* Lake tile processing improved, new product format (three shapefiles)
* Lake single pass processing added (cf leman and france_pekel new dataset to test it, can't be pushed on github, but can't be share through CNES cluster if needed)
```bash
% python ../../scripts/laketile_to_lakesp.py output/lakesp rdf/multi_lake_sp_command.cfg output/lake/lake-annotation_*
```


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
* CGAL (optional, if using HULL_METHOD=1.0)
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
  * pygeodesy
  
### Example of environment.yml

name: sht_env_test
channels:
  - conda-forge
  - defaults
dependencies:
  - cgal=4.12=py36h8634a1c_1
  - _libgcc_mutex=0.1=main
  - alabaster=0.7.12=py36_0
  - asn1crypto=0.24.0=py36_0
  - atomicwrites=1.3.0=py_0
  - attrs=18.2.0=py36h28b3542_0
  - babel=2.6.0=py36_0
  - backcall=0.1.0=py36_0
  - blas=1.0=mkl
  - bleach=3.0.2=py36_0
  - blosc=1.14.4=hdbcaa40_0
  - bokeh=1.0.3=py36_0
  - boost-cpp=1.67.0=h14c3975_4
  - bzip2=1.0.6=h14c3975_5
  - ca-certificates=2019.1.23=0
  - cairo=1.14.12=h8948797_3
  - certifi=2018.11.29=py36_0
  - cffi=1.11.5=py36he75722e_1
  - cftime=1.0.3.4=py36hdd07704_0
  - chardet=3.0.4=py36_1
  - click=7.0=py36_0
  - click-plugins=1.0.4=py36_0
  - cligj=0.5.0=py36_0
  - cloudpickle=0.6.1=py36_0
  - cryptography=2.4.2=py36h1ba5d50_0
  - curl=7.63.0=hbc83047_1000
  - cycler=0.10.0=py36_0
  - cytoolz=0.9.0.1=py36h14c3975_1
  - dask=1.0.0=py36_0
  - dask-core=1.0.0=py36_0
  - dbus=1.13.2=h714fa37_1
  - decorator=4.3.0=py36_0
  - descartes=1.1.0=py36_0
  - distributed=1.25.2=py36_0
  - docutils=0.14=py36_0
  - entrypoints=0.2.3=py36_2
  - expat=2.2.6=he6710b0_0
  - fiona=1.8.4=py36hc38cc03_0
  - fontconfig=2.13.0=h9420a91_0
  - freetype=2.9.1=h8a8886c_1
  - freexl=1.0.5=h14c3975_0
  - gdal=2.3.3=py36hbb2a789_0
  - geopandas=0.4.0=py36_1
  - geos=3.7.1=he6710b0_0
  - giflib=5.1.4=h14c3975_1
  - git=2.20.1=pl526hacde149_0
  - glib=2.56.2=hd408876_0
  - gmp=6.1.2=h6c8ec71_1
  - gst-plugins-base=1.14.0=hbbd80ab_1
  - gstreamer=1.14.0=hb453b48_1
  - hdf4=4.2.13=h3ca952b_2
  - hdf5=1.10.4=hb1b8bf9_0
  - heapdict=1.0.0=py36_2
  - icu=58.2=h9c2bf20_1
  - idna=2.8=py36_0
  - imageio=2.4.1=py36_0
  - imagesize=1.1.0=py36_0
  - intel-openmp=2019.1=144
  - ipykernel=5.1.0=py36h39e3cac_0
  - ipython=7.2.0=py36h39e3cac_0
  - ipython_genutils=0.2.0=py36_0
  - ipywidgets=7.4.2=py36_0
  - jedi=0.13.2=py36_0
  - jinja2=2.10=py36_0
  - jpeg=9b=h024ee3a_2
  - json-c=0.13.1=h1bed415_0
  - jsonschema=2.6.0=py36_0
  - jupyter=1.0.0=py36_7
  - jupyter_client=5.2.4=py36_0
  - jupyter_console=6.0.0=py36_0
  - jupyter_core=4.4.0=py36_0
  - kealib=1.4.7=hd0c454d_6
  - kiwisolver=1.0.1=py36hf484d3e_0
  - krb5=1.16.1=h173b8e3_7
  - libboost=1.67.0=h46d08c1_4
  - libcurl=7.63.0=h20c2e04_1000
  - libdap4=3.19.1=h6ec2957_0
  - libedit=3.1.20170329=h6b74fdf_2
  - libffi=3.2.1=hd88cf55_4
  - libgcc-ng=9.1.0=hdf63c60_0
  - libgdal=2.3.3=h2e7e64b_0
  - libgfortran-ng=7.3.0=hdf63c60_0
  - libkml=1.3.0=h590aaf7_4
  - libllvm10=10.0.1=hbcb73fb_5
  - libnetcdf=4.6.1=h11d0813_2
  - libpng=1.6.36=hbc83047_0
  - libpq=11.1=h20c2e04_0
  - libsodium=1.0.16=h1bed415_0
  - libspatialindex=1.8.5=h20b78c2_2
  - libspatialite=4.3.0a=hb08deb6_19
  - libssh2=1.8.0=h1ba5d50_4
  - libstdcxx-ng=8.2.0=hdf63c60_1
  - libtiff=4.0.9=he85c1e1_2
  - libuuid=1.0.3=h1bed415_2
  - libxcb=1.13=h1bed415_1
  - libxml2=2.9.8=h26e45fe_1
  - libxslt=1.1.32=h1312cb7_0
  - line_profiler=2.1.2=py36h14c3975_0
  - llvmlite=0.34.0=py36h269e1b5_4
  - locket=0.2.0=py36_1
  - lxml=4.3.0=py36hefd8a0e_0
  - lzo=2.10=h49e0be7_2
  - markupsafe=1.1.0=py36h7b6447c_0
  - matplotlib=3.0.2=py36h5429711_0
  - mistune=0.8.4=py36h7b6447c_0
  - mkl=2019.1=144
  - mkl_fft=1.0.10=py36ha843d7b_0
  - mkl_random=1.0.2=py36hd81dba3_0
  - more-itertools=5.0.0=py36_0
  - mpfr=4.0.1=hdf1c602_3
  - msgpack-python=0.5.6=py36h6bb024c_1
  - munch=2.3.2=py36_0
  - mypy=0.660=py36_0
  - mypy_extensions=0.4.1=py36_0
  - nbconvert=5.3.1=py36_0
  - nbformat=4.4.0=py36_0
  - ncurses=6.1=he6710b0_1
  - netcdf4=1.4.2=py36h808af73_0
  - networkx=2.2=py36_1
  - notebook=5.7.4=py36_0
  - numba=0.51.2=py36h0573a6f_1
  - numexpr=2.6.9=py36h9e4a6bb_0
  - numpy=1.15.4=py36h7e9f1db_0
  - numpy-base=1.15.4=py36hde5b4d6_0
  - numpydoc=0.8.0=py36_0
  - olefile=0.46=py36_0
  - openjpeg=2.3.0=h05c96fa_1
  - openssl=1.1.1b=h7b6447c_0
  - packaging=18.0=py36_0
  - pandas=0.23.4=py36h04863e7_0
  - pandoc=2.2.3.2=0
  - pandocfilters=1.4.2=py36_1
  - parso=0.3.1=py36_0
  - partd=0.3.9=py36_0
  - patsy=0.5.1=py36_0
  - pcre=8.42=h439df22_0
  - perl=5.26.2=h14c3975_0
  - pexpect=4.6.0=py36_0
  - pickleshare=0.7.5=py36_0
  - pillow=5.4.1=py36h34e0f95_0
  - pip=18.1=py36_0
  - pixman=0.36.0=h7b6447c_0
  - pluggy=0.8.1=py36_0
  - poppler=0.65.0=h581218d_1
  - poppler-data=0.4.9=0
  - proj4=5.2.0=he6710b0_1
  - prometheus_client=0.5.0=py36_0
  - prompt_toolkit=2.0.7=py36_0
  - psutil=5.4.8=py36h7b6447c_0
  - psycopg2=2.7.6.1=py36h1ba5d50_0
  - ptyprocess=0.6.0=py36_0
  - py=1.7.0=py36_0
  - pycparser=2.19=py36_0
  - pygments=2.3.1=py36_0
  - pyopenssl=18.0.0=py36_0
  - pyparsing=2.3.0=py36_0
  - pyproj=1.9.5.1=py36h14380d9_1
  - pyqt=5.9.2=py36h05f1152_2
  - pysal=1.14.4.post1=py36_1
  - pyshp=2.0.1=py36_0
  - pysocks=1.6.8=py36_0
  - pytables=3.4.4=py36h71ec239_0
  - pytest=4.2.0=py36_0
  - python=3.6.7=h0371630_0
  - python-dateutil=2.7.5=py36_0
  - pytz=2018.7=py36_0
  - pywavelets=1.0.1=py36hdd07704_0
  - pyyaml=3.13=py36h14c3975_0
  - pyzmq=17.1.2=py36h14c3975_0
  - qt=5.9.7=h5867ecd_1
  - qtconsole=4.4.3=py36_0
  - readline=7.0=h7b6447c_5
  - requests=2.21.0=py36_0
  - rtree=0.8.3=py36_0
  - scikit-image=0.14.1=py36he6710b0_0
  - scikit-learn=0.20.2=py36hd81dba3_0
  - scipy=1.1.0=py36h7c811a0_2
  - send2trash=1.5.0=py36_0
  - setuptools=40.6.3=py36_0
  - shapely=1.6.4=py36h86c5351_0
  - sip=4.19.8=py36hf484d3e_0
  - six=1.12.0=py36_0
  - snappy=1.1.7=hbae5bb6_3
  - snowballstemmer=1.2.1=py36_0
  - sortedcontainers=2.1.0=py36_0
  - sphinx=1.8.2=py36_0
  - sphinxcontrib=1.0=py36_1
  - sphinxcontrib-websupport=1.1.0=py36_1
  - sqlalchemy=1.2.16=py36h7b6447c_0
  - sqlite=3.26.0=h7b6447c_0
  - statsmodels=0.9.0=py36h035aef0_0
  - tbb=2020.3=hfd86e86_0
  - tbb4py=2020.3=py36hfd86e86_0
  - tblib=1.3.2=py36_0
  - terminado=0.8.1=py36_1
  - testpath=0.4.2=py36_0
  - tk=8.6.8=hbc83047_0
  - toolz=0.9.0=py36_0
  - tornado=5.1.1=py36h7b6447c_0
  - traitlets=4.3.2=py36_0
  - typed-ast=1.1.0=py36h14c3975_0
  - urllib3=1.24.1=py36_0
  - wcwidth=0.1.7=py36_0
  - webencodings=0.5.1=py36_1
  - wheel=0.32.3=py36_0
  - widgetsnbextension=3.4.2=py36_0
  - xarray=0.11.0=py36_0
  - xerces-c=3.2.2=h780794e_0
  - xz=5.2.4=h14c3975_4
  - yaml=0.1.7=had09818_2
  - zeromq=4.2.5=hf484d3e_1
  - zict=0.1.3=py36_0
  - zlib=1.2.11=h7b6447c_3
  
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

Then, you will need to install manually utm and pygeodesy packages with 'pip' or 'easy_install':
```
pip install utm
pip install PyGeodesy=='20.05.20
```

## Execute the toolbox

After activating your Python environment, you have to set your PYTHONPATH variables:
```bash
export PYTHONPATH=$SWOT_HYDROLOGY_TOOLBOX/processing/:$PYTHONPATH
export PYTHONPATH=$SWOT_HYDROLOGY_TOOLBOX/processing/src/:$PYTHONPATH
export PYTHONPATH=$SWOT_HYDROLOGY_TOOLBOX/processing/src/cnes/sas:$PYTHONPATH
export PYTHONPATH=$SWOT_HYDROLOGY_TOOLBOX/sisimp/:$PYTHONPATH
export PYTHONPATH=$RIVEROBS/src:$PYTHONPATH


```

An example dataset showing how to configure and run simulation and processing is available under /test.

The needed input for the overall chain include:
* An orbit file (provided)
* A water mask in shapefile format covering the area you want so simulate (see example under /test)
* A river database in shapefile format (e.g. GRWL)
* A lake database in shapefile format
* Various configuration files (examples provided)


The difference step are described below : 

The WIKI section contains a more details description of the simulation steps and the associated parameters.

### PIXC simulation by Large Scale Simulator
```bash
% python $SWOT_HYDROLOGY_TOOLBOX/select_orbit_cnes/select_orbit_cnes.py rdf/parameter_orbit.rdf output/orbit
% python $SWOT_HYDROLOGY_TOOLBOX/sisimp/proc_sisimp.py rdf/parameter_sisimp.rdf
```

### River processing
```bash
% python $SWOT_HYDROLOGY_TOOLBOX/scripts/l2pixc_to_rivertile.py output/simu output/river rdf/parameter_river.rdf
```

### LakeTile and LakeSP processing
```bash
% python $SWOT_HYDROLOGY_TOOLBOX/scripts/rivertile_to_laketile.py output/river output/laketile rdf/parameter_lak
etile.cfg -output_dir_lakesp output/lakesp
```

