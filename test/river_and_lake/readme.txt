How to run :

## Large Scale simulator
python $SWOT_HYDROLOGY_TOOLBOX/select_orbit_cnes/select_orbit_cnes.py rdf/parameter_orbit.rdf output/orbit
python $SWOT_HYDROLOGY_TOOLBOX/sisimp/proc_sisimp.py rdf/parameter_sisimp.rdf


## River processing
python $SWOT_HYDROLOGY_TOOLBOX/scripts/l2pixc_to_rivertile.py output/simu output/river rdf/parameter_river.rdf --nogdem

## Lake Processing
python $SWOT_HYDROLOGY_TOOLBOX/scripts/rivertile_to_laketile.py output/river output/lake -f --parameter_laketile rdf/parameter_laketile.cfg --nogdem

