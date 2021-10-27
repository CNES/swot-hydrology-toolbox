# 1. Generation des différentes dates
'''
Paramètres a modifier rdf/parameter_change_water.rdf
input_file - height_change - water_expansion_factor - height_attribute_name
'''
python process_change_water_area.py rdf/parameter_change_water_area.rdf

# 2. Generation des fichiers PIXC/PIXCVEC

cd ../../swot-hydrology-toolbox/test/<nouveau_dossier>
source clean.sh

# 2.1 Trouver les orbites correspondantes
'''
Paramètres a modifier rdf/parameter_orbit.rdf
DEM south latitude - DEM north latitude - DEM east longitude - DEM west longitude
'''
python $SWOT_HYDROLOGY_TOOLBOX/sisimp/proc_sisimp.py rdf/parameter_sisimp.rdf

# 2.2 Générer les fichiers pixc
'''
Paramètres a modifier rdf/parameter_sisimp.rdf
Shapefile path - Output directory - Multiple_orbit - Orbit - Cycle number
'''
python $SWOT_HYDROLOGY_TOOLBOX/sisimp/proc_sisimp.py rdf/parameter_sisimp.rdf

# 2.3 Lancer les scripts river & lake
python $SWOT_HYDROLOGY_TOOLBOX/scripts/l2pixc_to_rivertile.py output/simu output/river /work/ALT/swot/swotdev/desrochesd/swot-hydrology-toolbox/test/river_and_lake/rdf/parameter_river.rdf
python $SWOT_HYDROLOGY_TOOLBOX/scripts/rivertile_to_laketile.py output/river output/lake /work/ALT/swot/swotdev/desrochesd/swot-hydrology-toolbox/test/river_and_lake/rdf/parameter_laketile.cfg -output_dir_lakesp output/lakesp

# 3. Exécution des scripts floodplain
cd ../../../floodplain/scripts

# 3.1 Lancer le floodplain_dem
'''
Paramètres a modifier rdf/parameter_floodplain.rdf
latitude & longitude (min & max) - pixc file list - pixvec file list - output_directory
'''
python process_floodplain.py rdf/parameter_floodplain.rdf

# 3.2 Lancer le process_extract_area
'''
Paramètres a modifier rdf/parameter_extract_area.rdf
input_file
'''
python process_extract_area.py rdf/parameter_extract_area.rdf

# 3.3 Lancer le process_raster
'''
Paramètres a modifier rdf/parameter_raster.rdf
input_file
'''
python process_aster.py rdf/parameter_raster.rdf
