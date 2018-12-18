'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


import cnes.modules.geoloc.lib.my_netcdf_file as myCdf
import cnes.modules.geoloc.lib.tools as lib
import numpy as np
import itertools as itertools
import copy as copy


class PixelCloud(object):
    
    def __init__(self, IN_pixc_main_file, list_classif_flags):
    
        # 1 - Retrieve needed information from pixel cloud main file
        pixc_main = myCdf.myNcReader(IN_pixc_main_file)
        # 1.1 - Number of pixels in range dimension
        self.nb_pix_range = pixc_main.getAttValue("nr_pixels")
        # 1.2 - Number of pixels in azimuth dimension
        self.nb_pix_azimuth = pixc_main.getAttValue("nr_lines") 
        #~ # 1.3 - Cycle number
        #~ self.cycle_num = pixc_main.getAttValue("cycle_number")
        #~ # 1.4 - Pass number
        #~ self.pass_num = pixc_main.getAttValue("pass_number")
        #~ # 1.5 - Tile reference
        #~ self.tile_ref = pixc_main.getAttValue("tile_ref")
        # 1.6 - Classification flag
        self.classif = pixc_main.getVarValue("classification")       
        
        # Get indices corresponding to input classification flags
        TMP_classif = self.classif # Init temporary classif vector
        
        self.selected_idx = None # Init wanted indices vector
        self.nb_selected = 0

        for classif_flag in list_classif_flags:
            vInd = np.where(TMP_classif == int(classif_flag))[0]
            print("[PixelCloud] %d pixels with classification flag = %d" % ( vInd.size, int(classif_flag)))
            if ( vInd.size != 0 ):
                if ( self.nb_selected == 0 ):
                    self.selected_idx = vInd
                    self.nb_selected += vInd.size
                else:
                    self.selected_idx = np.array(list(itertools.chain(self.selected_idx,vInd)))
                    self.nb_selected += vInd.size
        self.nb_selected = self.selected_idx.size
        print("[PixelCloud] => %d pixels to keep" % ( self.nb_selected ))


        # Keep PixC data only for selected pixels
        if ( self.nb_selected != 0 ):
            # 1.7 - Range indices of water pixels
            self.range_idx = pixc_main.getVarValue("range_index")[self.selected_idx]
            # Number of water pixels
            self.nb_water_pix = self.range_idx.size
            # 1.8 - Azimuth indices of water pixels
            self.azimuth_idx = pixc_main.getVarValue("azimuth_index")[self.selected_idx]
            # 1.12 - Pixel area
            self.pixel_area = pixc_main.getVarValue("pixel_area")[self.selected_idx]
            # 1.14 - Longitude
            self.longitude = pixc_main.getVarValue("longitude")[self.selected_idx]
            # 1.15 - Latitude
            self.latitude = pixc_main.getVarValue("latitude")[self.selected_idx]
            # 1.16 - Height
            self.height = pixc_main.getVarValue("height")[self.selected_idx]

        pixc_main.close()


        
    def computeWaterMask(self):
        '''
        Create the water mask (i.e. a 2D binary matrix) in radar geometry,
        from the pixel cloud (1D-array layers of azimuth_index, range_index, classification and continuous classification)
        
        :return: water mask in radar geometry, i.e. a 2D matrix with "1" for each (IN_X_i, IN_Y_i)  having classification=0 and 0 elsewhere
        :rtype: 2D binary matrix of int 0/1
        '''
        print("[PixelCloud] == computeWaterMask ==")
        
        return lib.computeBinMat(self.nb_pix_range, self.nb_pix_azimuth, self.range_idx, self.azimuth_idx)
    
    def computeSeparateEntities(self):
        '''
        Identify all separate entities in the water mask
        '''
        print("[PixelCloud] == computeSeparateEntities ==")
        
        # 1 - Create the water mask
        water_mask = self.computeWaterMask()
        
        # 2 - Identify all separate entities in a 2D binary mask
        self.sepEntities, self.nb_obj = lib.labelRegion(water_mask)
        
        # 3 - Convert 2D labelled mask in 1D-array layer of the same size of the L2_HR_PIXC layers
        self.labels = lib.convert2dMatIn1dVec(self.range_idx, self.azimuth_idx, self.sepEntities) 
        
    def filterShorterEntities(self, area_min_entities):
        '''
        Filter the shortest entities in radar geometry from the pixel cloud using the pixel area
        Conpute the different labels and the separated entities (cleaned mask)
        '''
        # Dict with labels, area computed for each water body
        compter_area = lib.compter_area(self.labels, self.pixel_area)
        
        self.kept_labels=[]
        for cle,valeur in compter_area.items():
            if valeur > area_min_entities:
                self.kept_labels.append(cle)
        print("Number of large enough regions = ", len(self.kept_labels))

        sepEntities_tmp = copy.copy(self.sepEntities)
        sepEntities_tmp[:,:]=0
        labels_tmp = copy.copy(self.labels)
        labels_tmp[:]=0
                    
        # We keep only areas greater than area_min_entities in self.sepEntities
        for i in self.kept_labels:
            sepEntities_tmp += np.where(self.sepEntities == i, self.sepEntities, 0)
            labels_tmp += np.where(self.labels == i, self.labels, 0)
            
        print("Nb pixels to deal with after filtering = ", np.where(sepEntities_tmp!=0.)[0].size)

        self.sepEntities = sepEntities_tmp
        self.labels = labels_tmp
