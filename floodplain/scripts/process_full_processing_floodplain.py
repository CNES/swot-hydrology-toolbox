import argparse
import os
import logging
import itertools
import lib.my_rdf_file as my_rdf
import lib.my_timer as my_timer
import multiprocessing as mp
import pandas as pd

from floodplain.io.names import FPDEM_BASENAME, POLYGON_SUFFIX, FPDEM_POINTCLOUD_BASENAME, MASK_SUFFIX, compute_name
from functools import partial

from process_floodplain import Floodplain
from process_extract_area import Extract_Area
from process_raster import FPDEM_Raster

def process_fpdem_boundaries(fpdem):
    res_data = pd.DataFrame()
    res_poly = []
    res_wse=[]
            
    #### Multiprocessing option ######
    if fpdem.option == "multiprocessing":
        # Create multiprocessing Pool
        logging.info("multiprocessing")
        cpu_number = 12
        print("mp.cpu_count() = ", cpu_number)
        # ~ pool = mp.Pool(mp.cpu_count())
        pool = mp.Pool(cpu_number)
        # Input/Output files list
        res_poly=[]
        res_data=pd.DataFrame()
        args=[]
        
        for pixc_file, vec_file in zip(fpdem.inputpixcfiles,fpdem.inputvecfiles):
            args.append([pixc_file, vec_file])            
        
        res = pool.starmap(partial(fpdem.process_data), args)

        pool.close()
        pool.join()

        for i in range(len(res)):
            res_data=pd.concat([res_data, res[i][0]])
            res_poly += res[i][1]
            res_wse.append(res[i][2])

    elif fpdem.option == "dask":
        logging.info("dask option") 
        #### Dask option ####
        
        dask.config.set(scheduler='threads')
        cluster = PBSCluster(cores=fpdem.n_cores_per_worker, memory=fpdem.mem_per_worker, processes=fpdem.processes,
                             local_directory='$TMPDIR',
                             project='floodplain_dem',
                             name='floodplain_dem_worker',
                             walltime=fpdem.walltime,
                             interface='ib0',
                             log_directory=fpdem.output_path,
                             death_timeout=7200,
                             python='python -m tbb --ipc -a',
                             env_extra=[
                                 'module load conda',
                                 'conda activate {}'.format(fpdem.floodplain_dem_env),
                                 'export PYTHONPATH=$PYTHONPATH:{}'
                                 .format(fpdem.floodplain_src_path), 
                                 'export PYTHONPATH=$PYTHONPATH:{}'
                                 .format(fpdem.floodplain_scripts_path)
                             ])
        dask_client = Client(cluster)
        cluster.scale(jobs = fpdem.n_workers)  
        # important, actually creates workers
        with open(os.path.join(fpdem.output_path, 'bokeh_addr.txt'), 'w') as bokeh_file:
            bokeh_file.write(cluster.dashboard_link + '\n')

        input_data = list(zip(fpdem.inputpixcfiles,fpdem.inputvecfiles))
    
             
        ## Version delayed
        for i in range(len(input_data)):
            results = dask.delayed(fpdem.process_data)(input_data[i][0],input_data[i][1])
            res_data = dask.delayed(concat)(res_data, results[0])
            res_poly = dask.delayed(add)(res_poly, results[1])
            res_wse = dask.delayed(add)(res_wse, results[2])
            
        res = dask.delayed(return_output)(res_data, res_poly, res_wse)
        output = dask.compute(res)
        
        res_data, res_poly, res_wse = output[0][0], output[0][1], output[0][2]
        dask_client.close()
    
    # No multiprocessing
    else:
        for pixc_file, vec_file in zip(fpdem.inputpixcfiles, fpdem.inputvecfiles):
            result_data, result_poly, wse = fpdem.process_data(pixc_file, vec_file)
            res_data = pd.concat([res_data,result_data])
            res_poly += result_poly
            res_wse.append(wse)
         
    if fpdem.option != "dask":
        flatten = itertools.chain.from_iterable
        res_wse = list(flatten(res_wse))
                    
    fpdem.res_pointcloud = res_data
    fpdem.res_poly = res_poly
    fpdem.mean_wse = res_wse

# Main program
if __name__ == "__main__":

    print("===== Floodplain DEM processing = BEGIN =====")
    print("")
    timer = my_timer.Timer()
    timer.start()
     
    # Get parameters 
    parser = argparse.ArgumentParser(description="Compute fpdem intermediate product from multiple tiles of PIXC products and their associated PIXCVecRiver products.")
    parser.add_argument("parameter_file", help="parameter_file (*.rdf)")
    args = parser.parse_args()
    
    level = getattr(logging, "INFO")
    logging.basicConfig(filename=None,format='%(asctime)s [%(levelname)s] %(message)s', level=level)
    parameters = my_rdf.myRdfReader(args.parameter_file)


    # Prepare output names
    output_directory = parameters.getValue("output_directory")
    output_fpdem_pointcloud_name = compute_name(output_directory, FPDEM_POINTCLOUD_BASENAME, parameters.getValue("tile_name") \
                                ,parameters.getValue("first_date_name"), parameters.getValue("last_date_name"))
    output_poly_name = os.path.join(output_directory, FPDEM_BASENAME+POLYGON_SUFFIX)                            
    output_mask_name = os.path.join(output_directory, FPDEM_BASENAME+MASK_SUFFIX)
    output_fpdem_raster_name = compute_name(output_directory, FPDEM_BASENAME, parameters.getValue("tile_name") \
                                ,parameters.getValue("first_date_name"), parameters.getValue("last_date_name"))
                                
                                
    # Compute floodplain dem pixels boundaries extraction     
    fpdem = Floodplain(parameters)
    process_fpdem_boundaries(fpdem)
    # ~ fpdem.get_input_structure()
    #fpdem.compute_fpdem_pointcloud_boundaries()
    fpdem.write_fpdem_pointcloud_output()
    
    # Compute region of intereset where raster floodplain dem product in computed (between min and max water extent)
    if parameters.getValue("method_extract") == 'alpha_shape':
        extract_area = Extract_Area(parameters, input_file=output_fpdem_pointcloud_name, output_file=output_mask_name)
    if parameters.getValue("method_extract") == 'intersection':
        extract_area = Extract_Area(parameters, input_file=output_poly_name, output_file=output_mask_name)
        
    extract_area.load_input_extract_area()
    extract_area.compute_extract_area()
    extract_area.write_extract_area()
    
     # Compute rasterization of pixel cloud boundaries and produce final FPDEM product 
     
    fpdem_raster = FPDEM_Raster(parameters, input_file=output_fpdem_pointcloud_name \
                                          , output_file=output_fpdem_raster_name \
                                          , mask=output_mask_name)
    fpdem_raster.load_input_fpdem_raster()
    fpdem_raster.compute_fpdem_raster()
    fpdem_raster.write_fpdem_raster()
    
   
    print("")
    print(timer.stop())
    print("===== Floodplain DEM processing  = END =====")        
