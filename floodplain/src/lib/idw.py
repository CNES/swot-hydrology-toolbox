#!/usr/bin/env python

# -*- coding: utf-8 -*-



# Copyright (C) 2016 Paul Brodersen <paulbrodersen+idw@gmail.com>



# Author: Paul Brodersen <paulbrodersen+idw@gmail.com>



# This program is free software; you can redistribute it and/or

# modify it under the terms of the GNU General Public License

# as published by the Free Software Foundation; either version 3

# of the License, or (at your option) any later version.



# This program is distributed in the hope that it will be useful,

# but WITHOUT ANY WARRANTY; without even the implied warranty of

# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

# GNU General Public License for more details.



# You should have received a copy of the GNU General Public License

# along with this program. If not, see <http://www.gnu.org/licenses/>.

# Modification from Damien Desroches, 05/10/2021


"""

Inverse distance weighting (IDW)

--------------------------------



Compute the score of query points based on the scores of their k-nearest neighbours,

weighted by the inverse of their distances.



Example:

--------



# import idw



# 'train'

idw_tree = idw.tree(X1, z1)



# 'test'

spacing = np.linspace(-5., 5., 100)

X2 = np.meshgrid(spacing, spacing)

grid_shape = X2[0].shape

X2 = np.reshape(X2, (2, -1)).T

z2 = idw_tree(X2)





For a more complete example see demo().



"""

from __future__ import absolute_import, division, print_function, unicode_literals 


import numpy as np

from scipy.spatial import cKDTree



class tree(object):

    """

    Compute the score of query points based on the scores of their k-nearest neighbours,

    weighted by the inverse of their distances.



    @reference:

    https://en.wikipedia.org/wiki/Inverse_distance_weighting



    Arguments:

    ----------

        X: (N, d) ndarray

            Coordinates of N sample points in a d-dimensional space.

        z: (N,) ndarray

            Corresponding scores.

        leafsize: int (default 10)

            Leafsize of KD-tree data structure;

            should be less than 20.



    Returns:

    --------

        tree instance: object



    Example:

    --------



    # 'train'

    idw_tree = tree(X1, z1)



    # 'test'

    spacing = np.linspace(-5., 5., 100)

    X2 = np.meshgrid(spacing, spacing)

    X2 = np.reshape(X2, (2, -1)).T

    z2 = idw_tree(X2)



    See also:

    ---------

    demo()



    """

    def __init__(self, X=None, z=None, leafsize=10):

        if not X is None:

            self.tree = cKDTree(X, leafsize=leafsize )

        if not z is None:

            self.z = np.array(z)



    def fit(self, X=None, z=None, leafsize=10):

        """

        Instantiate KDtree for fast query of k-nearest neighbour distances.



        Arguments:

        ----------

            X: (N, d) ndarray

                Coordinates of N sample points in a d-dimensional space.

            z: (N,) ndarray

                Corresponding scores.

            leafsize: int (default 10)

                Leafsize of KD-tree data structure;

                should be less than 20.



        Returns:

        --------

            idw_tree instance: object



        Notes:

        -------

        Wrapper around __init__().



        """

        return self.__init__(X, z, leafsize)



    def __call__(self, X, k=6, eps=1e-6, p=2, regularize_by=1e-9):

        """

        Compute the score of query points based on the scores of their k-nearest neighbours,

        weighted by the inverse of their distances.



        Arguments:

        ----------

            X: (N, d) ndarray

                Coordinates of N query points in a d-dimensional space.



            k: int (default 6)

                Number of nearest neighbours to use.



            p: int or inf

                Which Minkowski p-norm to use.

                1 is the sum-of-absolute-values "Manhattan" distance

                2 is the usual Euclidean distance

                infinity is the maximum-coordinate-difference distance



            eps: float (default 1e-6)

                Return approximate nearest neighbors; the k-th returned value

                is guaranteed to be no further than (1+eps) times the

                distance to the real k-th nearest neighbor.



            regularise_by: float (default 1e-9)

                Regularise distances to prevent division by zero

                for sample points with the same location as query points.



        Returns:

        --------

            z: (N,) ndarray

                Corresponding scores.

        """

        self.distances, self.idx = self.tree.query(X, k, eps=eps, p=p)

        self.distances += regularize_by

        weights = self.z[self.idx.ravel()].reshape(self.idx.shape)

        mw = np.sum(weights/self.distances, axis=1) / np.sum(1./self.distances, axis=1)
        
        mean_dist = np.mean(self.distances, axis=1)
        mean_z_rel = np.abs(np.std(weights, axis=1)/np.mean(weights, axis=1))

        return mw, np.min(self.distances, axis=1), mean_dist, mean_z_rel


    def transform(self, X, k=6, p=2, eps=1e-6, regularize_by=1e-9):

        """

        Compute the score of query points based on the scores of their k-nearest neighbours,

        weighted by the inverse of their distances.



        Arguments:

        ----------

            X: (N, d) ndarray

                Coordinates of N query points in a d-dimensional space.



            k: int (default 6)

                Number of nearest neighbours to use.



            p: int or inf

                Which Minkowski p-norm to use.

                1 is the sum-of-absolute-values "Manhattan" distance

                2 is the usual Euclidean distance

                infinity is the maximum-coordinate-difference distance



            eps: float (default 1e-6)

                Return approximate nearest neighbors; the k-th returned value

                is guaranteed to be no further than (1+eps) times the

                distance to the real k-th nearest neighbor.



            regularise_by: float (default 1e-9)

                Regularise distances to prevent division by zero

                for sample points with the same location as query points.



        Returns:

        --------

            z: (N,) ndarray

                Corresponding scores.



        Notes:

        ------



        Wrapper around __call__().

        """

        return self.__call__(X, k, eps, p, regularize_by)



def demo():

    import matplotlib.pyplot as plt



    # create sample points with structured scores

    X1 = 10 * np.random.rand(1000, 2) -5



    def func(x, y):

        return np.sin(x**2 + y**2) / (x**2 + y**2)



    z1 = func(X1[:,0], X1[:,1])



    # 'train'

    idw_tree = tree(X1, z1)



    # 'test'

    spacing = np.linspace(-5., 5., 100)

    X2 = np.meshgrid(spacing, spacing)

    grid_shape = X2[0].shape

    X2 = np.reshape(X2, (2, -1)).T

    z2 = idw_tree(X2)



    # plot

    fig, (ax1, ax2, ax3) = plt.subplots(1,3, sharex=True, sharey=True, figsize=(10,3))

    ax1.contourf(spacing, spacing, func(*np.meshgrid(spacing, spacing)))

    ax1.set_title('Ground truth')

    ax2.scatter(X1[:,0], X1[:,1], c=z1, linewidths=0)

    ax2.set_title('Samples')

    ax3.contourf(spacing, spacing, z2.reshape(grid_shape))

    ax3.set_title('Reconstruction')

    plt.show()

    return
