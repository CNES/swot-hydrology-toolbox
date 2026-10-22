#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2019/05/17:version initiale.
# VERSION:3.3.0:DM:#91:2024/09/20:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''

import cnes.modules.geoloc.lib.geoloc as geoloc
import numpy as np



p = np.load('/work/SWOT/swotdev/swotdev/desrochesd/swot-sds/swotCNES/src/cnes/modules/geoloc/scripts/data_test/p.save.npy')
h_noisy = np.load('/work/SWOT/swotdev/swotdev/desrochesd/swot-sds/swotCNES/src/cnes/modules/geoloc/scripts/data_test/h_noisy.save.npy')
s = np.load('/work/SWOT/swotdev/swotdev/desrochesd/swot-sds/swotCNES/src/cnes/modules/geoloc/scripts/data_test/s.save.npy')
vs = np.load('/work/SWOT/swotdev/swotdev/desrochesd/swot-sds/swotCNES/src/cnes/modules/geoloc/scripts/data_test/vs.save.npy')
R_target = np.load('/work/SWOT/swotdev/swotdev/desrochesd/swot-sds/swotCNES/src/cnes/modules/geoloc/scripts/data_test/r.save.npy')
h_target = np.load('/work/SWOT/swotdev/swotdev/desrochesd/swot-sds/swotCNES/src/cnes/modules/geoloc/scripts/data_test/h_target.save.npy')
    
         
                                        
p_corr = geoloc.pointcloud_height_geoloc_vect(p, h_noisy, s, vs, R_target, h_target, recompute_Doppler=True, 
                                        recompute_R=True,
                                        verbose = False,
                                        max_iter_grad=1, height_goal = 1.e-3, safe_flag=True)
                                        

np.save('/work/SWOT/swotdev/swotdev/desrochesd/swot-sds/swotCNES/src/cnes/modules/geoloc/scripts/data_test/p_corr.save', p_corr[0])
np.save('/work/SWOT/swotdev/swotdev/desrochesd/swot-sds/swotCNES/src/cnes/modules/geoloc/scripts/data_test/p_corr_llh.save', p_corr[1])

