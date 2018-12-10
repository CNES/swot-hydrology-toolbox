# -*- coding: utf-8 -*-
"""
.. module:: my_passplan.py
    :synopsis: Deal with the orbit file named passplan.txt (= list of cycle/orbits occuring on the studed area during the simulation time period), resulting from SAM select_orbit 
    
    passplan.txt file is similar to:
        
        # Mission start:    2014-01-01 00:00:00
        # Simulation start: 2015-04-01 00:00:00
        # Simulation stop:  2015-05-01 00:00:00
        #
        #run      cycle orbit MissionTime year DayOfYear       date     time autoclean
        c022_t017     22   17    39706831 2015  95.56980 2015-04-05 13:40:31     False

.. moduleauthor:: Claire POTTIER - CNES DSO/SI/TR

Copyright (c) 2018 CNES. All rights reserved.
"""
from __future__ import absolute_import, division, print_function, unicode_literals


class orbitPassplan(object):
    
    def __init__(self, IN_file):
        """
        Constructor
        
        :param IN_file: full path of the passplan.txt file
        :type IN_file: string
        """
        
        # 1 - Set filename attribute
        self.filename = IN_file
        
        # 2 - List of cycle number / orbit number pairs
        self.cycle_orbit_pairs = []  # Init the list
        self.readFile()  # Read file and fill the list
        
        # 3 - Dictionnary giving correspondance between orbit numbers and their indices in self.cycle_orbit_pairs
        self.plan_orbit = {}  # Init the dictionnary
        self.computeOrbitPlan()  # Compute plan according to orbit numbers
    
    #----------------------------------
        
    def readFile(self):
        """
        Read the file and store variables
        """
        
        # 1 - Read the file
        plan_reader = open(self.filename, 'r') # Open file in reading mode
        plan_lines = plan_reader.read().splitlines() # Read file and store line / line
        
        # 2 - Skip comment lines and store cycle number / orbit number pairs in a table
        for line in plan_lines:
            
            # Consider line only if contain "=" (not empty, not a comment, ...)
            if not line.startswith("#"):
                
                # 2.1 - Cut regarding blank 
                TMP_line = line.split(" ")
                
                # 2.2 - Remove empty elements
                TMP_line = [x for x in TMP_line if x != '']
                
                # 2.3 - Add pair
                self.cycle_orbit_pairs.append([int(TMP_line[1]), int(TMP_line[2])])
                
        # 3 - Close file
        plan_reader.close()
    
    #----------------------------------
    
    def getOrbitList(self):
        """
        Retrieve the list of orbits, i.e. the 2nd column of self.cycle_orbit_pairs
        NB: some may be repeated several times
        
        :return: the list of orbits
        :rtype: list of int
        """
        
        OUT_orbit_list = []
        
        for pair in self.cycle_orbit_pairs:
            OUT_orbit_list.append(pair[1])
            
        return OUT_orbit_list
    
    def computeOrbitPlan(self):
        """
        Find all indices corresponding to each orbit number in the list of pairs
        """
        
        orbit_list = self.getOrbitList()  # List of all orbit values according to the plan
        orbit_list_uniq = set(orbit_list)  # Set of uniq orbit values
        
        for orbit_number in orbit_list_uniq:
            
            ind_list = []
            for ind, elem in enumerate(orbit_list):
                if elem == orbit_number:
                    ind_list.append(ind)
            
            self.plan_orbit[orbit_number] = ind_list
    
    #----------------------------------
    
    def printPairs(self):
        """
        Print list of cycle number / orbit number pairs  on screen
        """
        
        print("Cycle\tOrbit")
        for pair in self.cycle_orbit_pairs:
            print("%03d\t%03d" % (pair[0], pair[1]))
            
            
#######################################


if __name__ == "__main__":
    
    plan = orbitPassplan("C:\\Users\\pottierc\\Documents\\workspace_qgis\\sisimp_tests_unitaires\\FT2_Cas2\\1_select_orbit\\passplan.txt")
    
    print("Number of pairs = %d" % len(plan.cycle_orbit_pairs))
    print("")
    plan.printPairs()
    print("")
    print("Orbit plan =")
    print(plan.getOrbitList())
    print("")
    print("Orbit indices in plan =")
    print(plan.plan_orbit)