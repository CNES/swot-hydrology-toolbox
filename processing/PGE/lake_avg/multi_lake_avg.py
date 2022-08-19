#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:3.2.0:DM:#91:2021/10/27:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: multi_lake_avg.py
    :synopsis: Process PGE_L2_HR_LakeAvg, i.e. generate L2_HR_LakeAvg
    product from multiple cycles of L2_HR_LakeSP products over multiple basins 
    Created on 05/03/2021

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
import configparser as cfg
import datetime
from lxml import etree as ET
import multiprocessing as mp
import os
import sys

import cnes.common.lib.my_tools as my_tools
import cnes.common.lib.my_timer as my_timer
import cnes.common.lib_lake.lake_db as lake_db
import cnes.common.lib_lake.locnes_filenames as locnes_filenames

import pge_lake_avg


class MultiLakeAvg(object):
    """
    Class MultiLakeAvg
    Main class to run many SAS LakeAvg processes
    """
    
    def __init__(self, in_params):
        """
        Constructor: initialize variables
        
        :param in_params: input parameters to run the processor
        :type in_params: dictionary
        """
        print("[multiLakeAvgProcessing] == INIT ==")

        # Directories
        self.param_file = in_params["param_file"]  # Parameter file
        self.lakesp_dir = in_params["lakesp_dir"]  # LakeSP_Prior files directory
        self.output_dir = in_params["output_dir"]  # Output directory

        # PLD infos
        self.lake_db = in_params["LAKE_DB"]
        self.lake_db_id = in_params["LAKE_DB_ID"]

        # Cycle number and basin code
        self.cycle_num = my_params["cycle_num"]  # Cycle number
        self.basin_code = my_params["basin_code"]  # Basin code
        
        # Options
        self.flag_inc_file_counter = in_params["flag_inc_file_counter"]
        self.flag_write_full_path = in_params["flag_write_full_path"]
        self.nb_proc = in_params["nb_proc"]
        
        # Log level
        self.error_file = in_params["errorFile"]
        self.log_file = in_params["logFile"]
        self.log_file_level = in_params["logfilelevel"]
        self.log_console = in_params["logConsole"]
        self.log_console_level = in_params["logconsolelevel"]
        
        # File information
        self.institution = in_params["institution"]
        self.product_version = in_params["product_version"]
        self.crid_lakesp = in_params["crid_lakesp"]
        self.crid_lakeavg = in_params["crid_lakeavg"]
        self.pge_version = in_params["pge_version"]
        self.contact = in_params["contact"]

        # List of input files
        self.list_cycle_basin = dict()  # List of cycles and basin_code to handle
        self.cmd_file_path_list = []  # List of command files
        
        # Nb_proc parameter
        self.nb_proc = in_params["nb_proc"]
    # -------------------------------------------

    def run_preprocessing(self):
        """
        Retrieve the list of input files, i.e. L2_HR_LakeSP_Prior files
        """

        print("")
        print("")
        print("[multiLakeAvgProcessing] PRE-PROCESSING...")
        print("")

        # 1 - Test existence of directories
        print("[multiLakeAvgProcessing] > 1 - Testing existence of working directories...")
        # 1.1 - LakeSP_Prior directory
        print("[multiLakeAvgProcessing]   INPUT DIR for LakeSP_Prior files = %s" % self.lakesp_dir)
        my_tools.test_dir(self.lakesp_dir)
        # 1.2 - Output directory
        print("[multiLakeAvgProcessing]   OUTPUT DIR = %s" % self.output_dir)
        my_tools.test_dir(self.output_dir)
        print("")

        # 2 - Get input files
        print("[multiLakeAvgProcessing] > 2 - Retrieving input files...")

        # 2.1 - Compute file prefix regarding cycle condition and basin_code in LakeSP_Prior.shp.xml files
        lakesp_prefix = locnes_filenames.LAKE_SP_PREFIX["prior"]
        lakesp_suffix = locnes_filenames.LAKE_SP_SHP_SUFFIX
        cond_prefix = lakesp_prefix  # Deal with all LakeSP_Prior files in self.lakesp_dir

        flag_list = False
        cond_continent = None
        if (self.cycle_num is None) or (self.cycle_num == "-1"):
            flag_list = True
            
            if (self.basin_code is None) or (self.basin_code == "-1"):
                print("[multiLakeAvgProcessing]   All LakeSP_Prior files in the input directory")
                
            else:  # Deal with LakeSP_Prior files with basin_code = self.basin_code
                continent = lake_db.compute_continent_id_from_basin_code(str(self.basin_code)[0])
                cond_continent = "_%s_" % continent
                if self.basin_code < 10:  # Manage all basins of continent
                    print("[multiLakeAvgProcessing]   LakeSP_Prior files with continent=%s (digit=%s)" % (continent, str(self.basin_code)[0]))
                    
                else:
                    print("[multiLakeAvgProcessing]   LakeSP_Prior files with basin_code=%02d (continent=%s)" % (self.basin_code, continent))
            
        else:  # Deal with LakeSP_Prior files with cycle = self.cycle_num
            cond_prefix += "%03d" % self.cycle_num
            if (self.basin_code is None) or (self.basin_code == "-1"):
                print("[multiLakeAvgProcessing]   LakeSP_Prior files with cycle=%03d" % self.cycle_num)
                flag_list = True
            else:  # Deal with LakeSP_Prior files with cycle = self.cycle_num and basin_code = self.basin_code
                continent = lake_db.compute_continent_id_from_basin_code(str(self.basin_code)[0])
                cond_continent = "_%s_" % continent
                if self.basin_code < 10:  # Manage all basins of continent
                    print("[multiLakeAvgProcessing]   LakeSP_Prior files with cycle=%03d and continent=%s (digit=%s)" % (self.cycle_num, continent, str(self.basin_code)[0]))
                    flag_list = True
                else:
                    self.list_cycle_basin[self.cycle_num] = {self.basin_code}
                    print("[multiLakeAvgProcessing]   LakeSP_Prior files with cycle=%03d and basin_code=%02d (continent=%s)" % (self.cycle_num, self.basin_code, continent))

        if flag_list:

            # 2.2 - List all files in self.lakesp_dir
            tmp_list_lakesp_files = os.listdir(self.lakesp_dir)
    
            # 2.3 - For each listed file, select LakeSP_Prior shapefiles corresponding to specified cycle number and basin code 
            for cur_file in tmp_list_lakesp_files:
    
                # Test if file meets the conditions
                flag_ok = False
                if cond_continent is None:
                    if cur_file.startswith(cond_prefix) and cur_file.endswith(lakesp_suffix):
                        flag_ok = True
                else:
                    if cur_file.startswith(cond_prefix) and cur_file.endswith(lakesp_suffix) and (cond_continent in cur_file):
                        flag_ok = True
                        
                # Process if yes
                if flag_ok:
                        
                    # 2.4 - Retrieve cycle number
                    cycle_num = int(cur_file.split("_")[5])
                    if cycle_num not in self.list_cycle_basin.keys():
                        self.list_cycle_basin[cycle_num] = set()
                            
                    # 2.5 - Get basin_code list from LakeSP_Prior metadata file (ie .shp.xml)
                    metadata = ET.parse(os.path.join(self.lakesp_dir, cur_file + ".xml"))
                    # Retrieve basin_code 
                    cur_basin_code_txt = metadata.xpath("//swot_product/global_metadata/basin_code")[0].text
                    if not cur_basin_code_txt:
                        cur_basin_code_list = [""]
                    else:
                        cur_basin_code_list = cur_basin_code_txt.split(";")
                        
                    # 2.6 - Add CB basins to list of selected basins for the cycle
                    for cur_basin_code in cur_basin_code_list:
                        if cur_basin_code != "":
                            self.list_cycle_basin[cycle_num].add(int(cur_basin_code[0:2]))
        
        # Print some numbers
        tmp_print = "[multiLakeAvgProcessing]   * Cycle numbers to process = "
        for cycle_num in self.list_cycle_basin.keys():
            tmp_print += "%03d - " % cycle_num
        print(tmp_print)
        for cycle_num, cycle_basins in self.list_cycle_basin.items():
            tmp_print = "[multiLakeAvgProcessing]   * Cycle number %03d => basin IDs = " % cycle_num
            for basin_code in cycle_basins:
                tmp_print += "%02d - " % basin_code
            print(tmp_print)
        print("")
        
        # 3 - Writing command files
        print("[multiLakeAvgProcessing] > 3 - Writing command files...")
        for cycle_num, cycle_basins in self.list_cycle_basin.items():
            for basin_code in cycle_basins:
                print("[multiLakeAvgProcessing]   * Writing command file cycle_num=%03d / basin_code=%02d" %(cycle_num, basin_code))
                cmd_file_path = self.create_cmd_file(cycle_num, basin_code)
                self.cmd_file_path_list.append(cmd_file_path)
        print("")

    def create_cmd_file(self, in_cycle_num, in_basin_code):
        """
        Create command file for PGE_L2_HR_LakeAvg with CYCLE_NUM=in_cycle_num and BASIN_CODE=in_basin_code
        
        :param in_cycle_num: cycle number to write in command file
        :type in_cycle_num: int
        :param in_basin_code: basin code to write in command file
        :type in_basin_code: int
        
        :return: out_cmd_file = command file full path
        :rtype: string
        """

        # 0.1 - Init error log filename (without date because it is computed in pge_lake_tile.py)
        if self.error_file is None:
            error_file = os.path.join(self.output_dir, 
                                      "ErrorFile_%03d_%02d.log" % (in_cycle_num, in_basin_code))
        else:
            error_file = "%s_%03d_%02d%s" % (os.path.splitext(self.error_file)[0], 
                                             in_cycle_num, in_basin_code, 
                                             os.path.splitext(self.error_file)[1])
        # 0.2 - Init log filename (without date because it is computed in pge_lake_tile.py)
        if self.log_file is None:
            log_file = os.path.join(self.output_dir, 
                                    "LogFile_%03d_%02d.log" % (in_cycle_num, in_basin_code))
        else:
            log_file = "%s_%03d_%02d%s" % (os.path.splitext(self.log_file)[0], 
                                           in_cycle_num, in_basin_code, 
                                           os.path.splitext(self.log_file)[1])
        
        # 1 - Init command filename
        cmd_filename = "lake_avg_command_%03d_%02d.cfg" % (in_cycle_num, in_basin_code)
        out_cmd_file = os.path.join(self.output_dir, cmd_filename)
        
        # 2 - Open command file 
        writer_cmd_file = open(out_cmd_file, "w")  # Open file in writing mode
        
        # 3 - Write command file
        
        # 3.1 - Fill PATHS section
        writer_cmd_file.write("[PATHS]\n")
        if self.param_file is None:
            # Needed by Jenkins script
            writer_cmd_file.write("param_file = %s\n" % os.path.join(sys.path[0], "lake_avg_param.cfg"))
        else:
            writer_cmd_file.write("param_file = %s\n" % self.param_file) 
        writer_cmd_file.write("LakeSP_Prior directory = %s\n" % self.lakesp_dir)
        writer_cmd_file.write("Output directory = %s\n\n" % self.output_dir)
        
        # 3.2 - Fill DATABASES section
        if self.lake_db is not None:
            writer_cmd_file.write("[DATABASES]\n")
            writer_cmd_file.write("LAKE_DB = " + self.lake_db + "\n")
            if self.lake_db.endswith(".shp"):
                if self.lake_db_id is not None:
                    writer_cmd_file.write("LAKE_DB_ID = " + self.lake_db_id + "\n")
            writer_cmd_file.write("\n")
        
        # 3.3 - Fill PASS_INFOS section
        writer_cmd_file.write("[PASS_INFOS]\n")
        writer_cmd_file.write("Cycle number = %03d\n" % in_cycle_num)
        writer_cmd_file.write("Basin code = %02d\n\n" % in_basin_code)
        
        # 3.4 - Fill OPTIONS section
        writer_cmd_file.write("[OPTIONS]\n")
        writer_cmd_file.write("Increment file counter = " + str(self.flag_inc_file_counter) + "\n")
        writer_cmd_file.write("Write full path = " + str(self.flag_write_full_path) + "\n")
        writer_cmd_file.write("Nb_proc = " + str(self.nb_proc) + "\n\n")
        
        # 3.5 - Fill LOGGING section
        writer_cmd_file.write("[LOGGING]\n")
        writer_cmd_file.write("errorFile = " + error_file + "\n")
        writer_cmd_file.write("logFile = " + log_file + "\n")
        writer_cmd_file.write("logfilelevel = " + self.log_file_level + "\n")
        writer_cmd_file.write("logConsole = " + str(self.log_console) + "\n")
        writer_cmd_file.write("logconsolelevel = " + self.log_console_level + "\n\n")

        # 3.6 - Fill FILE_INFORMATION section
        writer_cmd_file.write("[FILE_INFORMATION]\n")    
        writer_cmd_file.write("INSTITUTION = " + self.institution + "\n")
        writer_cmd_file.write("PRODUCT_VERSION = " + self.product_version + "\n")
        writer_cmd_file.write("CRID_LAKESP = " + self.crid_lakesp + "\n")
        writer_cmd_file.write("CRID_LAKEAVG = " + self.crid_lakeavg + "\n")
        writer_cmd_file.write("PGE_VERSION = " + self.pge_version + "\n")
        writer_cmd_file.write("CONTACT = " + self.contact + "\n")
        
        # 4 - Close command file
        writer_cmd_file.close() 
        
        return out_cmd_file
        
    # -------------------------------------------

    def run_processing(self, cmd_file=None):
        """
        Process SAS_L2_HR_LakeAvg for each input command file corresponding to cycle number / basin code pairs
        """

        print("")
        print("")
        print("[multiLakeAvgProcessing] PROCESSING...")
        print("")
        print("")
        
        timer_proc = my_timer.Timer()
        timer_proc.start()  # Init timer
        
        if not cmd_file:
            for indf, cmd_file in enumerate(self.cmd_file_path_list):
                print("")
                print("")
                print("***********************************************")
                print("***** Dealing with command file %d / %d *****" % (indf+1, len(self.cmd_file_path_list)))
                print("***********************************************")
                print("")
                print("")

                # 1 - Init
                my_lake_avg = pge_lake_avg.PGELakeAvg(cmd_file)

                # 2 - Run
                my_lake_avg.start()

                # 3 - Stop
                my_lake_avg.stop()

        else:
            print("")
            print("")
            print("***********************************************")
            print("***** Dealing with command file %s *****" % cmd_file)
            print("***********************************************")
            print("")
            print("")

            # 1 - Init
            my_lake_avg = pge_lake_avg.PGELakeAvg(cmd_file)

            # 2 - Run
            my_lake_avg.start()

            # 3 - Stop
            my_lake_avg.stop()

        print("")
        print("")
        print(timer_proc.stop())  # Print tile process duration
        print("")
        print("")

    def run_multiprocessing(self):
        """
        Multiprocess SAS_L2_HR_LakeAvg
        """
        n_cores = int(mp.cpu_count())
        pool = mp.Pool(n_cores)
        with pool:
            print('Running map')
            print("[multiLakeAvgProcessing] PROCESSING...")
            tmp_result = pool.map(self.run_processing, self.cmd_file_path_list)
            pool.close()
            pool.join()
    

#######################################


def read_command_file(in_filename):
    """
    Read the command file in input and store parameters in a dictionary
    
    :param in_filename: parameter file full path
    :type in_filename: string
    
    :return: out_params = dictionary containing parameters
    :rtype: dict
    """
    print("[multiLakeAvgProcessing] == read_command_file = %s ==" % in_filename)

    # 0 - Init output dictionary
    out_params = {}
    # Default values
    out_params["param_file"] = None
    out_params["LAKE_DB"] = None
    out_params["LAKE_DB_ID"] = None
    out_params["cycle_num"] = None
    out_params["basin_code"] = None
    out_params["flag_inc_file_counter"] = True
    out_params["flag_write_full_path"] = False
    out_params["nb_proc"] = 1
    out_params["errorFile"] = None
    out_params["logFile"] = None
    out_params["logfilelevel"] = "DEBUG"
    out_params["logConsole"] = True
    out_params["logconsolelevel"] = "DEBUG"
    out_params["institution"] = "CNES"
    out_params["product_version"] = "0.0"
    out_params["crid_lakesp"] = "Dx0000"
    out_params["crid_lakeavg"] = "Dx0000"
    out_params["pge_version"] = "0.0"
    out_params["contact"] = "test@cnes.fr"

    # 1 - Read parameter file
    config = cfg.ConfigParser()
    config.read(in_filename)

    # 2 - Retrieve PATHS
    list_paths = config.options("PATHS")
    if "param_file" in list_paths:
        out_params["param_file"] = config.get("PATHS", "param_file")
    out_params["lakesp_dir"] = config.get("PATHS", "LakeSP_Prior directory")
    out_params["output_dir"] = config.get("PATHS", "Output directory")
    
    # 3 - Retrieve DATABASES
    if "DATABASES" in config.sections():
        list_db = config.options("DATABASES")
        # Prior lake database
        if "lake_db" in list_db:
            out_params["LAKE_DB"] = config.get("DATABASES", "LAKE_DB")
        if "lake_db_id" in list_db:
            out_params["LAKE_DB_ID"] = config.get("DATABASES", "LAKE_DB_ID")
            
    # 4 - Retrieve PASS_INFOS
    if "PASS_INFOS" in config.sections():
        list_options = config.options("PASS_INFOS")
        # Cycle number
        if "cycle number" in list_options:
            out_params["cycle_num"] = config.getint("PASS_INFOS", "Cycle number")
        # Basin code
        if "basin code" in list_options:
            out_params["basin_code"] = config.getint("PASS_INFOS", "Basin code")

    # 5 - Retrieve OPTIONS
    if "OPTIONS" in config.sections():
        list_options = config.options("OPTIONS")
        # Flag to increment the file counter in the output filenames (=True, default); else=False
        if "increment file counter" in list_options:
            out_params["flag_inc_file_counter"] = config.get("OPTIONS", "Increment file counter")
        # To write full path in global attributes (=True); to write only basename=False (default)
        if "write full path" in list_options:
            out_params["flag_write_full_path"] = config.get("OPTIONS", "Write full path")
        # Number of processors to use (default=1)
        if "nb_proc" in list_options:
            out_params["nb_proc"] = config.get("OPTIONS", "Nb_proc")
        else:
            out_params["nb_proc"] = 1

    # 6 - Retrieve LOGGING
    if "LOGGING" in config.sections():
        list_options = config.options("LOGGING")
        if "errorFile" in list_options:
            out_params["errorFile"] = config.get("LOGGING", "errorFile")
            if out_params["errorFile"] is not None:
                out_params["errorFile"] = out_params["errorFile"].replace("<date>", datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
        if "logFile" in list_options:
            out_params["logFile"] = config.get("LOGGING", "logFile")
            if out_params["logFile"] is not None:
                out_params["logFile"] = out_params["logFile"].replace("<date>", datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
        if "logfilelevel" in list_options:
            out_params["logfilelevel"] = config.get("LOGGING", "logfilelevel")
        if "logConsole" in list_options:
            out_params["logConsole"] = config.get("LOGGING", "logConsole")
        if "logconsolelevel" in list_options:
            out_params["logconsolelevel"] = config.get("LOGGING", "logconsolelevel")
            
    # 7 - Retrieve FILE_INFORMATION
    if "FILE_INFORMATION" in config.sections():
        list_options = config.options("FILE_INFORMATION")
        if "institution" in list_options:
            out_params["institution"] = config.get("FILE_INFORMATION", "INSTITUTION")
        if "product_version" in list_options:
            out_params["product_version"] = config.get("FILE_INFORMATION", "PRODUCT_VERSION")
        if "crid_lakesp" in list_options:
            out_params["crid_lakesp"] = config.get("FILE_INFORMATION", "CRID_LAKESP")
        if "crid_lakeavg" in list_options:
            out_params["crid_lakeavg"] = config.get("FILE_INFORMATION", "CRID_LAKEAVG")
        if "pge_version" in list_options:
            out_params["pge_version"] = config.get("FILE_INFORMATION", "PGE_VERSION")
        if "contact" in list_options:
            out_params["contact"] = config.get("FILE_INFORMATION", "CONTACT")

    return out_params
    

#######################################


if __name__ == '__main__':

    # 0 - Parse inline parameters
    parser = argparse.ArgumentParser(description="Compute SWOT LakeAvg products from multiple LakeSP_Prior files.")
    parser.add_argument("command_file", help="command file (*.cfg)")
    parser.add_argument("-mp", "--multiproc", help="if true, CB basins will be computed in parallel", nargs='?', type=bool,
                        default=False, const=True)
    args = parser.parse_args()

    print("===== multiLakeAvgProcessing = BEGIN =====")
    print("")
    timer = my_timer.Timer()
    timer.start()

    # 1 - Read command file
    print("WORKING VARIABLES")
    print()
    my_tools.test_file(args.command_file, in_extent=".cfg")  # Test existance and extension
    my_params = read_command_file(args.command_file)  # Read variables in command file

    # 2 - Initialization
    multi_lake_avg = MultiLakeAvg(my_params)
    print(timer.info(0))

    # 3 - Run pre-processing
    multi_lake_avg.run_preprocessing()
    print(timer.info(0))
    
    # 4 - Run processing
    if args.multiproc:
        multi_lake_avg.run_multiprocessing()
    else:
        multi_lake_avg.run_processing()
    
    print(timer.info(0))

    print("")
    print(timer.stop())
    print("===== multiLakeAvgProcessing = END =====")
