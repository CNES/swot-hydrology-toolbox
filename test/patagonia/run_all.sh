source clean.sh
cd /work/ALT/swot/swothr/users/loussertp/simu_opensource/swot-hydrology-toolbox/test/patagonia/
python $SWOT_HYDROLOGY_TOOLBOX/sisimp/proc_sisimp.py rdf/parameter_sisimp.rdf

## River processing
python $SWOT_HYDROLOGY_TOOLBOX/scripts/l2pixc_to_rivertile.py output/simu output/river rdf/parameter_river.rdf 

## Lake Processing
python $SWOT_HYDROLOGY_TOOLBOX/scripts/rivertile_to_laketile.py output/river output/lake rdf/parameter_laketile.cfg --nogdem
