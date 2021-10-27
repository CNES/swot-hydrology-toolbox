# -*- coding: utf8 -*-
'''
Create ply or shp from pixel cloud

Copyright (c) 2018, CNES
'''

import argparse
import os
import sys
import logging
import datetime
import traceback
import pandas as pd
import numpy as np
import multiprocessing as mp
import utm
import numpy as np

import psutil

from typing import List
from functools import partial

import itertools  

from floodplain.process import compute_mean_wse, extract_water_points, extract_contiguous_water_points, remove_isolated_points, compute_alpha_shape, extract_boundary_points, remove_borders, remove_near_range_pixels
from floodplain.io.nc import PixcReader
import floodplain.io.ply as ply
import floodplain.io.shp as shp
from floodplain.geom.tools import filter_polygons, extract_polygon_boundaries
from floodplain.utils.spatial import convert_polygons_utm_to_latlon
from floodplain.io.names import FPDEM_BASENAME, POLYGON_SUFFIX, FPDEM_POINTCLOUD_BASENAME, compute_name
from floodplain.io.nc import write_raster_ungridded

import lib.my_rdf_file as my_rdf

import dask.bag as db
from dask.distributed import Client
from dask_jobqueue import PBSCluster
import dask


class Floodplain(object):
    """
    Class Floodplain
    Main class to run many Floodplain processes
    """
    
    def __init__(self, param):
        """
        Constructor: initialize variables
        
        :param in_params: input parameters to run the processor
        :type in_params: dictionary
        """
        self.lat_max = float(param.getValue("latitude max"))
        self.lon_max = float(param.getValue("longitude max"))
        self.lat_min = float(param.getValue("latitude min"))
        self.lon_min = float(param.getValue("longitude min"))
        self.cross_track_min = float(param.getValue("cross track min value"))
        self.threshold = 10000.

        self.pixcfiles = param.getValue("pixc file list").split(" ")
        self.vecfiles = param.getValue("pixcvec file list").split(" ")
        
        
        
        self.option = param.getValue("threading option")
        self.output_path = param.getValue("output_directory")
        
        self.tile_name = param.getValue("tile_name") 
        self.first_date_name = param.getValue("first_date_name")
        self.last_date_name = param.getValue("last_date_name")
        
        ## For dask multiprocessing
        self.n_workers = 12
        self.n_cores_per_worker = 1
        self.processes = 1
        self.mem_per_worker = '{}GiB'.format(self.n_cores_per_worker * 60)   
        self.floodplain_src_path = '/work/ALT/swot/swotdev/desrochesd/floodplain/src'
        self.floodplain_scripts_path = '/work/ALT/swot/swotdev/desrochesd/floodplain/scripts'
        self.floodplain_dem_env = '/work/ALT/swot/swotpub/modules/conda-envs/swot-floodplain-env-20190327'
        self.walltime = '12:00:00'
        

        inputpixcfiles=os.popen("ls "+self.pixcfiles[0]).readlines()
        inputvecfiles=os.popen("ls "+self.vecfiles[0]).readlines()
        
        self.inputpixcfiles = []
        self.inputvecfiles = []
        
        for i in range(0,len(inputpixcfiles)):
            for k in range(0,len(inputvecfiles)): 
                if (inputpixcfiles[i].rstrip('\n')).split("_PIXC_")[1] == (inputvecfiles[k].rstrip('\n')).split("_PIXCVec_")[1] :
                    self.inputpixcfiles.append(inputpixcfiles[i].rstrip('\n'))
                    self.inputvecfiles.append(inputvecfiles[k].rstrip('\n'))

        if len(self.inputpixcfiles) != len(inputpixcfiles):
            logging.warning(f"{len(inputpixcfiles)-len(self.inputpixcfiles)} pixvec file is missing")
            

        
        logging.info(f"{len(self.inputvecfiles)} PIXC/PIXCVec couple files will be processed")
        

    def compute_fpdem_pointcloud_boundaries(self):
        '''
        Create cloud or shp file from input files (netcdf).

        :param parameters: parameters file 
     
        '''
        
        res_data = pd.DataFrame()
        res_poly = []
        res_wse=[]
            
        #### Multiprocessing option ######
        if self.option == "multiprocessing":
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
            
            for pixc_file, vec_file in zip(self.inputpixcfiles,self.inputvecfiles):
                args.append([pixc_file, vec_file])            
            
            res = pool.starmap(partial(self.process_data), args)

            pool.close()
            pool.join()

            for i in range(len(res)):
                res_data=pd.concat([res_data, res[i][0]])
                res_poly += res[i][1]
                res_wse.append(res[i][2])

        elif self.option == "dask":
            logging.info("dask option") 
            #### Dask option ####
            
            dask.config.set(scheduler='threads')
            cluster = PBSCluster(cores=self.n_cores_per_worker, memory=self.mem_per_worker, processes=self.processes,
                                 local_directory='$TMPDIR',
                                 project='floodplain_dem',
                                 name='floodplain_dem_worker',
                                 walltime=self.walltime,
                                 interface='ib0',
                                 log_directory=self.output_path,
                                 death_timeout=7200,
                                 python='python -m tbb --ipc -a',
                                 env_extra=[
                                     'module load conda',
                                     'conda activate {}'.format(self.floodplain_dem_env),
                                     'export PYTHONPATH=$PYTHONPATH:{}'
                                     .format(self.floodplain_src_path), 
                                     'export PYTHONPATH=$PYTHONPATH:{}'
                                     .format(self.floodplain_scripts_path)
                                 ])
            dask_client = Client(cluster)
            cluster.scale(jobs = self.n_workers)  
            # important, actually creates workers
            with open(os.path.join(self.output_path, 'bokeh_addr.txt'), 'w') as bokeh_file:
                bokeh_file.write(cluster.dashboard_link + '\n')

            input_data = list(zip(self.inputpixcfiles,self.inputvecfiles))
        
                 
            ## Version delayed
            for i in range(len(input_data)):
                results = dask.delayed(self.process_data)(input_data[i][0],input_data[i][1])
                res_data = dask.delayed(concat)(res_data, results[0])
                res_poly = dask.delayed(add)(res_poly, results[1])
                res_wse = dask.delayed(add)(res_wse, results[2])
                
            res = dask.delayed(return_output)(res_data, res_poly, res_wse)
            output = dask.compute(res)
            
            res_data, res_poly, res_wse = output[0][0], output[0][1], output[0][2]
            dask_client.close()
        
        # No multiprocessing
        else:
            for pixc_file, vec_file in zip(self.inputpixcfiles, self.inputvecfiles):
                result_data, result_poly, wse = self.process_data(pixc_file, vec_file)
                res_data = pd.concat([res_data,result_data])
                res_poly += result_poly
                res_wse.append(wse)
             
        if self.option != "dask":
            flatten = itertools.chain.from_iterable
            res_wse = list(flatten(res_wse))
                        
        self.res_pointcloud = res_data
        self.res_poly = res_poly
        self.mean_wse = res_wse
        

    def process_data(self, pixc_file, vec_file):
        res_poly_process=[]
        try:
            pixc_reader = PixcReader(pixc_file, vec_file)
            range_max = pixc_reader.range_max
            azimuth_max = pixc_reader.azimuth_max
            
            # Extract points
            water = extract_water_points(pixc_reader.get_data())
            # Remove pixels to close to near_range
            water, min_range_indices_to_remove = remove_near_range_pixels(water, azimuth_max, cross_track_min = self.cross_track_min)
            
            # Extract contiguous water points by keeping only water bodies whose areas is  greater than threshold
            water = extract_contiguous_water_points(water,
                                                    range_max,
                                                    azimuth_max,
                                                    threshold=self.threshold)
            ## TODO : NEIGHBOOR POINTS REMOVING DESACTIVATED, need to improve the criteria
            water = remove_isolated_points(water, pixc_reader.along_track_sampling,
                                           factor=1.2,
                                           nb_neighbors=0)
                                           
            # Compute alpha shape polygonization to produce boundaries shapefile
            polygons = compute_alpha_shape(water, pixc_reader.range_max, "cgal")
            
            # Post treatment
            
            boundaries = extract_polygon_boundaries(polygons)
            data = extract_boundary_points(water, boundaries)
            filtered_polygons = filter_polygons(polygons)
            mean_wse = compute_mean_wse(water, filtered_polygons)
            
            filtered_data = remove_borders(data, pixc_reader.range_max-1, pixc_reader.azimuth_max-1, min_range_indices_to_remove)
            filtered_data = filtered_data.loc[(filtered_data['longitude'] > self.lon_min) & (filtered_data['longitude'] < self.lon_max)
                                & (filtered_data['latitude'] > self.lat_min) & (filtered_data['latitude'] < self.lat_max)]

            
            # Convert polygon
            if len(filtered_data)>0:
                zone_number = utm.latlon_to_zone_number(filtered_data.iloc[0]['latitude'],
                                                        filtered_data.iloc[0]['longitude'])
                north = True if data.iloc[0]['latitude'] > 0 else False
                filtered_polygons_latlon = convert_polygons_utm_to_latlon(filtered_polygons,zone_number,north)
                res_poly_process = filtered_polygons_latlon
                logging.info(f"Process files: {pixc_file} {vec_file}")
                logging.info("Processing: OK")
            else:
                logging.warning(f"Output is empty for {pixc_file} and {vec_file}")
        except:
            traceback.print_exc()
            logging.warning(f"Error in files {pixc_file} and {vec_file}")    
        


        try:
            logging.info(f"Process files: {pixc_file} {vec_file}")

            return filtered_data, res_poly_process, mean_wse
        except:
            logging.info(f"Process files: {pixc_file} {vec_file}")
            filtered_data = pd.DataFrame()
            res_poly_process = []
            return filtered_data, res_poly_process, []


    def write_fpdem_pointcloud_output(self):
                # Combine outputs into on product        
        if len(self.res_pointcloud)>0:
            try:
               
                # Add valid_date and epsg attributes
                res_pointcloud = valid_date(self.res_pointcloud)
                
                outputfile_root = os.path.join(self.output_path, FPDEM_BASENAME)
                
                ply.gdf_to_file(outputfile_root+'.ply',res_pointcloud,mode="text")
                
                # Also write shp file
                shp.gdf_to_file(outputfile_root + ".shp",res_pointcloud, index=True)
                shp.polygons_to_file(outputfile_root + POLYGON_SUFFIX,self.res_poly, wse=self.mean_wse)

                #TODO : write netcdf file
                res_pointcloud_ds = res_pointcloud[['longitude', 'latitude', 'elevation', 'x', 'y', 'z', 'data_valid']].to_xarray()
                
                zone_number = utm.latlon_to_zone_number(res_pointcloud.iloc[0]['latitude'], res_pointcloud.iloc[0]['longitude'])
                if np.mean(res_pointcloud.iloc[0]['latitude']) >= 0:
                    res_pointcloud_ds.attrs['espg']  = int(32600 + zone_number) 
                else:
                    res_pointcloud_ds.attrs['espg']  = int(32700 + zone_number)
                
                output_fpdem_pointcloud_name = compute_name(self.output_path, FPDEM_POINTCLOUD_BASENAME, self.tile_name \
                                    ,self.first_date_name, self.last_date_name)
                
                write_raster_ungridded(res_pointcloud_ds, output_fpdem_pointcloud_name)


            except:
                traceback.print_exc()
                logging.error(f"Writing output file")
            else:
                logging.info("Write output file: OK")
        else:
            logging.warning("Output is empty for this dataset")



def return_output(a,b,c):
    return [a, b, c]

def concat(a, b):
    return pd.concat([a, b])

def add(a,b):
    if b != 0:
        return a+b
    else:
        return a
        
def valid_date(dataFrame):
    '''
    Create valid_date attribute to create raster file
    '''
    dataFrame['data_valid']=pd.Series([1 for i in range(dataFrame.size)])
    return dataFrame
    
# Main program
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Compute fpdem intermediate product from multiple tiles of PIXC products and their associated PIXCVecRiver products.")
    parser.add_argument("parameter_file", help="parameter_file (*.rdf)")
    args = parser.parse_args()

    level = getattr(logging, "INFO")
    logging.basicConfig(filename=None,format='%(asctime)s [%(levelname)s] %(message)s', level=level)

    parameters = my_rdf.myRdfReader(args.parameter_file)
    fpdem = Floodplain(parameters)
    fpdem.compute_fpdem_pointcloud_boundaries()
    fpdem.write_fpdem_pointcloud_output()
