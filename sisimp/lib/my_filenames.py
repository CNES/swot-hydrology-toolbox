# -*- coding: utf-8 -*-
"""
.. module:: my_filenames.py
    :synopsis: Deal with filenames

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

Copyright (c) 2018 CNES. All rights reserved.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

from datetime import datetime, timedelta
import os

from lib.my_variables import PATTERN_FOOTPRINT, PATTERN_PIXC, PATTERN_FILE_ANNOT, PATTERN_PIXC_VEC_RIVER

class sisimpFilenames(object):

    def __init__(self, IN_out_dir, IN_mission_start_time, IN_cycle_duration, IN_cycle_num, IN_pass_num):
        """
        Constructor of SISIMP filenames
        
        :param IN_out_dir: output directory
        :type IN_out_dir: string
        :param IN_mission_start_time: mission start time
        :type IN_mission_start_time: string (yyyy-mm-dd)
        :param IN_cycle_duration: number of seconds in a cycle
        :type IN_cycle_duration: int
        :param IN_cycle_num: cycle number
        :type IN_cycle_num: int
        :param IN_pass_num: pass number
        :type IN_pass_num: int
        """
        
        # Init variables
        self.out_dir = IN_out_dir
        self.mission_start_time = datetime.strptime(IN_mission_start_time, '%Y-%m-%d')
        self.cycle_duration = IN_cycle_duration
        self.cycle_num = IN_cycle_num
        self.pass_num = IN_pass_num
        self.tile_ref = None
        self.begin_date = None
        self.end_date = None
        
        # Init filenames
        self.computeFootprintFilename()
        self.pixc_file = None
        self.pixc_vec_river_file = None
        self.file_annot_file = None
        
    def updateWithTileRef(self, IN_tile_ref, IN_sec_in_cycle_start, IN_sec_in_cycle_end):
        """
        Update filenames with tile ref info
        
        :param IN_tile_ref: tile reference
        :type IN_tile_ref: string
        :param IN_sec_in_cycle_start: start date of tile, as number of seconds from begin of current cycle
        :type IN_sec_in_cycle_start: int
        :param IN_sec_in_cycle_end: end date of tile, as number of seconds from begin of current cycle
        :type IN_sec_in_cycle_end: int
        """
        
        self.tile_ref = IN_tile_ref
        
        # Compute tile start date
        self.begin_date = self.computeDate(IN_sec_in_cycle_start)
        # Compute tile end date
        self.end_date = self.computeDate(IN_sec_in_cycle_end)
        
        # Set filenames
        self.computePixcFilename()
        self.computePixcVecRiverFilename()
        self.computeFileAnnotFilename()
    
    #----------------------------------
    
    def computeDate(self, IN_sec_in_cycle):
        """
        Compute date
        
        :param IN_sec_in_cycle: number of seconds from begin of current cycle
        :type IN_sec_in_cycle: int
        
        :return: date in UTC
        :rtype: string YYYYMMDDThhmmss
        """
        
        # Computation
        if self.cycle_num == 0:
            date_in_sec = self.mission_start_time + timedelta(seconds=self.cycle_duration*self.cycle_num+IN_sec_in_cycle)
        else:
            date_in_sec = self.mission_start_time + timedelta(seconds=self.cycle_duration*(self.cycle_num-1)+IN_sec_in_cycle)
        
        # Format
        return datetime.strftime(date_in_sec, '%Y%m%dT%H%M%S')
    
    #----------------------------------
    
    def computeFootprintFilename(self):
        """
        Compute footprint full path
        """
        filename = PATTERN_FOOTPRINT % (self.cycle_num, self.pass_num)
        self.footprint_file = os.path.join(self.out_dir, filename)
    
    def computePixcFilename(self):
        """
        Compute pixel cloud full path
        """
        filename = PATTERN_PIXC % (self.cycle_num, self.pass_num, self.tile_ref, self.begin_date, self.end_date)
        self.pixc_file = os.path.join(self.out_dir, filename) 
        
    def computeFileAnnotFilename(self):
        """
        Compute file annotation full path
        """
        filename = PATTERN_FILE_ANNOT % (self.cycle_num, self.pass_num, self.tile_ref)
        self.file_annot_file = os.path.join(self.out_dir, filename)
        
    def computePixcVecRiverFilename(self):
        """
        Compute PIXCVecRiver full path
        """
        filename = PATTERN_PIXC_VEC_RIVER % (self.cycle_num, self.pass_num, self.tile_ref, self.begin_date, self.end_date)
        self.pixc_vec_river_file = os.path.join(self.out_dir, filename)
         