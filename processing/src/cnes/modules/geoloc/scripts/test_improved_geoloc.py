#!/usr/bin/env python
'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''

import cnes.modules.geoloc.lib.geoloc as geoloc
import numpy as np

p = np.array([  (604036.2612028 ,-5151738.12882137 ,3700899.70950969), ( 604036.2612028 ,-5151738.12882137 ,3700899.70950969)])
s = np.array([  (691911.1260945 ,-5874186.70691969 ,4224481.28190375), ( 691911.1260945 ,-5874186.70691969 ,4224481.28190375)])
vs = np.array([ (0.27379512 ,-0.54234502 ,-0.79429094), ( 0.27379512 ,-0.54234502 ,-0.79429094)])
R_target = np.array([896544.255928, 896544.255928])

h_noisy = geoloc.project_array(p)[:,2]
h_target = h_noisy + 2.

p_corr = np.copy(p)
p_corr[:,:] = 0.
p_corr_llh = np.copy(p)
p_corr_llh[:,:] = 0.


for i in range(p.shape[0]):
    toto = geoloc.pointcloud_height_geoloc(p[i], h_noisy[i], s[i], vs[i], R_target[i], h_target[i], recompute_Doppler=True, 
                                            recompute_R=True,
                                            verbose = False,
                                            max_iter_grad=1, height_goal = 1.e-3, safe_flag=True)
    p_corr[i], p_corr_llh[i] = toto[0], toto[1] 
         
                                        
titi = geoloc.pointcloud_height_geoloc_vect(p, h_noisy, s, vs, R_target, h_target, recompute_Doppler=True, 
                                        recompute_R=True,
                                        verbose = False,
                                        max_iter_grad=1, height_goal = 1.e-3, safe_flag=True)
                                        
p_corr_vect, p_corr_vect_llh = titi[0], titi[1]
                                        
                                    
                                        
print('')
print('p =  ', p)
print('p_llh =  ', geoloc.project_array(p))

print('p_corr = ', p_corr)
print('p_corr_llh = ', p_corr_llh)
print('p_corr_vect = ', p_corr_vect)
print('p_corr_vect_llh = ', p_corr_vect_llh)

print('')
print('delta_h = ',  p_corr_llh[:,2]-h_target)
print('delta_h_vect = ',  p_corr_vect_llh[:,2]-h_target)


