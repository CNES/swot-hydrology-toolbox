'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


import numpy as np
from scipy import spatial

dmax = 5000.


class Function_mean_h(object):
    ### Different functions using dask parallelism to estimate DEM height on a grid using pixel cloud height. Description on the calcul_h_array function. Other are derivated
    def __init__(self, nx, ny, dem, nb_points):
        # x size of recomputed DEM
        self.nx = nx
        # y size of recomputed DEM
        self.ny = ny
        # x,y coordinates of recomputed DEM
        self.dem = dem
        # number of points used to estimate height averaging (TO BE CHANGED TO CORRESPOND TO A 1x1 KM2 AVERAGING)
        self.nb_points = nb_points

    
    def calcul_h(self, h, dist_max):
        # h : pixel heights
        # dist_max : distance minimal to have between a pixel from pix. cloud of the grid pixel on the DEM to estimate an height on the DEM
        # KD Tree of h values (to improve speed computation)
        tree = spatial.KDTree(h[0])
        # Array with h coordinates on the DEM grid
        dem = np.zeros([self.nx,self.ny])

        for m in range(self.nx):
            for n in range(self.ny):
                 # xi, yi coordinates of the DEM ( TO BE CHECKED : BOUNDARIES OR MIDDLE OF THE PIXEL)   
                 xi, yi = self.dem[m, n, 0], self.dem[m, n, 1]
                 # Criteria to calculate a mean(h) on the DEM pixel (TO BE CHANGED, UPDATED : MORE PIXELS, RADAR COORDINATES, INTERSECTION BETWEEN CONCAVE HULL ESTIMATED FROM PIXC AND DEM ??)
                 if tree.query([xi, yi],k=1)[0] < dist_max:
                     # h averaging of the nb_points closest to the grid center
                    dist, ind = tree.query([xi, yi],k = self.nb_points)
                    # Weighting using the distance between points 
                    weights = (dmax-dist)/dmax
                    weights = (np.where(weights < 0., 0., weights))
                    weights = weights*weights
                    dem[m, n] = np.average(h[1][ind], weights = weights)

        return dem


    def calcul_h_array(self, array, h, dist_max):
        # Go to calcul_h
        dem = np.zeros([len(array),len(array[0])])
        for m in range(len(array)):
            for n in range(len(array[0])):
                xi, yi = array[m, n, 0], array[m, n, 1]
                if self.tree.query([xi, yi],k=1)[0] < dist_max:
                    dist, ind = self.tree.query([xi, yi],k = self.nb_points)
                    weights = (dmax-dist)/dmax
                    weights = (np.where(weights < 0., 0., weights))
                    weights = weights*weights
                    dem[m, n] = np.average(h[ind], weights = weights)

        return dem

    def calcul_h_map(self, array, h, tree, dist_max):
        # Go to calcul_h
        xi, yi = array[0], array[1]
        if tree.query([xi, yi],k=1)[0] < dist_max:
            dist, ind = tree.query([xi, yi],k = self.nb_points)
            #~ dem = np.mean(h[ind])
            weights = (dmax-dist)/dmax
            weights = (np.where(weights < 0., 0., weights))
            weights = weights*weights
            dem[m, n] = np.average(h[ind], weights = weights)


    def calcul_h_map_tuple(self, data, dist_max, no_close_point_value=0):
        # Go to calcul_h

        id_i, id_j, label = data[0]  # /!\ label not used
        coords = data[1]
        pts = [(x, y) for x, y, _ in coords]
        h = np.array([c[2] for c in coords])
        tree = spatial.KDTree(pts)
        centroid = self.dem[id_i, id_j, :2]
        if tree.query(centroid ,k=1)[0] < dist_max:
            dist, ind = tree.query(centroid, k = min(self.nb_points, h.size))
            weights = (dmax-dist)/dmax
            weights = (np.where(weights < 0., 0., weights))
            weights = weights*weights
            #~ print(weights)
            print(h.size)
            #~ print(ind.size, weights.size)
            return (data[0], np.average(h[ind], weights = weights))
        else:
            return (data[0], 0.)
