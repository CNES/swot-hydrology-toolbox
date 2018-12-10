# pixc_simulator
Large-scale pixelc cloud simulator.

CNES sisimp code delivered by Claire 14Feb2018

This module is dedicated to simulation of pixel clouds from water mask shapefiles at large scale.

# Options
Processing options are set in the RDF file, an example is given in the current project in the file parameter_sisimp.rdf.

Options of parameter_file are :
* Noise file path : path to noise file. An example of noise file is given in file height_noise_presum2.txt
* Run directory for orbits : path to the orbit folder. Folling files and directories should be avaliable :
    * makeGdemOrbit_2.01 directory
    * passplan.txt
    * passplan_tXXX.txt, XXX being the orbit number
* Shapefile path : path to the water mask shapefile without extention ".shp"
* Geolocalisation improvement : yes (no cross-track noise applied) or no
* Real pixc file format : yes (PixC as defined in PDD) or no(PixC as provided by classic simulator)
* Create shapefile : yes or no (works only with real pixel cloud option)
* Create dummy pixc vec river file : yes (L2_HR_PIXC_VEC_RIVER product associated to PixC files) or no (works only with real pixel cloud option)
* Water flag : classification flag for water (usually 4)
* Output directory : path to the output directory
* Cycle : cycle number
* Orbit : orbit number

<<<<<<< HEAD
# Options for height calculation

Options of parameter_file :
* Height model : name of the model to apply (polynomial, gaussian, height_file, height_shapefile)
* Height name : field where find height value
=======
Options for Dark Water simulation are : 
* Dark water =  yes or no If yes dark water simulation is executed
* Scale factor dw = float to parameterize the correlation length of simulated dark water areas
* Dark water percentage = int (between 0 and 100) probability of total dark water simulated
* Dark water flag = classification flag for detected dark water (usually 24)
* Dark water seed = int seed to be used for reproducible dark water simulations
>>>>>>> origin/dark_water

# Run

To run sisimp, use this command in a terminal :
```
python sisimp.py [-h] [-v [VERBOSE]] [-l [LOGFILE]] parameter_file
```

optional arguments:
* -h, --help            show this help message and exit
* -v [VERBOSE], --verbose [VERBOSE] Verbose level : DEBUG INFO
* -l Write prints to a logfile, if not -l write prints in the console

# Workflow
On comming ...
