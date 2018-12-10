# Installation procedure

## Get the repository
1. Clone __swot_hydrology_toolbox__ repo

```bash
% git clone https://gitlab.cnes.fr/swotdev/swot-hydrology-toolbox.git
```
The repository __swot_hydrology_toolbox__ should be assignated to the __SWOT_HYDROLOGY_TOOLBOX__ variable.

```bash
% export SWOT_HYDROLOGY_TOOLBOX=your_installation_path/swot-hydrology-toolbox
```
2. Clone __RiverObs__ repo

If you are on the CNES network, you should need to bypass the proxi by using these commands :
```bash
% export http_proxy='http://CNESNET\login:password@proxy-http2.cnes.fr:8050'
% export https_proxy='http://CNESNET\login:password@proxy-http2.cnes.fr:8050'
% export ftp_proxy='http://CNESNET\login:password@proxy-http2.cnes.fr:8050'
```

```bash
% git clone https://github.com/SWOTAlgorithms/RiverObs.git
```


## Dependancies

The specific dependancies for swot are installed here on the HAL cluster: /work/ALT/swot/swotpub/modules.
To load it on HAL, you need to set the variable MODULEPATH.

```bash
export MODULEPATH=/work/ALT/swot/swotpub/modules/modulefiles/:$MODULEPATH
```

* Module hdf/5.1.8.12 (specific module SWOT)
* Module libspatialindex/1.8.5 (specific module SWOT)
* Module gdal/2.1.1-py3.5
* Module netcdf/4.4.1
* Module pyenv/20180411 (specific SWOT python environment)

You can also make your own [propre virtualenv](creation_virtualenv).

## Creation of modulefile for swot-hydrology-toolbox

if you don't have, you need a repository to stock the modulesfiles files
For example
```bash
mkdir $HOME/modulefiles
```
You need to add this repository to the variable __MODULEPATH__.
```bash
% export MODULEPATH=<repository_modulesfiles>:$MODULEPATH
```
In the examples :
```bash
% export MODULEPATH=$HOME/modulefiles:$MODULEPATH
```

In this repository, create a repository swot-hydrology-toolbox.
```bash
mkdir -p $HOME/modulefiles/swot-hydrology-toolbox
```

In the repository swot-hydrology-toolbox, you need to create the file develop.lua which correspond to a version of swot-hydrology-toolbox

```bash
% cat swot-hydrology-toolbox/develop.lua
-- -*- lua -*-
-- Aide du module accessible avec la commande module help
help(
[[
Module permettant de charger l'environnement pour exécuter la chaine du JPL (swot-sds).
]])

-- Information du modulefile
local os_disponible = "rh7"
local nom           = "swot-hydrology-toolbox"
local version       = "develop"
local informations  = system_information()
local rhos          = informations.os

-- Information du module accessible avec la commande module whatis
whatis("Nom     : "..nom)
whatis("Version : "..version)
whatis("Os      : "..os_disponible)

-- Verification du modulefile
check_os(os_disponible)

-- Variable du modulefile  
local home=REPLACE_ME (example: "/work/ALT/swot/swotdev/desrochesd/swot-hydrology-toolbox")
local home_riverobs=REPLACE_ME (example: "/work/ALT/swot/swotdev/desrochesd/RiverObs")

-- Dépendance du modulefile
depend("cmake")
depend("pyenv/20180411")
depend("libspatialindex")
depend("gdal/2.1.1-py3.5")
depend("netcdf/4.4.1")
depend("hdf/5.1.8.12")

-- Action du modulefile 
setenv("SWOT_HYDROLOGY_TOOLBOX",home)
setenv("RIVEROBS", home_riverobs)
setenv("CURRENTARCH","LINUX")
prepend_path("PYTHONPATH",pathJoin(home_riverobs,"src"))
prepend_path("PYTHONPATH",pathJoin(home,"processing/src"))
prepend_path("PYTHONPATH",pathJoin(home,"processing/src/cnes/sas"))

prepend_path("PATH",pathJoin(home_riverobs,"src/bin"))

```

You need to replace __REPLACE_ME__ by the path corresponding to swot_hydrology_toolbox  and by the version of the python module to load