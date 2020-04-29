source clean.sh
python $SWOT_HYDROLOGY_TOOLBOX/sisimp/proc_sisimp.py rdf/parameter_sisimp.rdf

## River processing
python $SWOT_HYDROLOGY_TOOLBOX/scripts/l2pixc_to_rivertile.py output/simu output/river rdf/parameter_river.rdf

## Lake Processing
python $SWOT_HYDROLOGY_TOOLBOX/scripts/rivertile_to_laketile.py output/river output/lake rdf/parameter_laketile.cfg 

## LakeSP Processing
python $SWOT_HYDROLOGY_TOOLBOX/scripts/laketile_to_lakesp.py output/lakesp rdf/multi_lake_sp_command.cfg output/lake/lake-annotation_*
