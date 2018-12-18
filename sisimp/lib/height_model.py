'''
 This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.
'''


import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


## Scale between lac_size and correlation window applied to create the random field

def generate_1d_profile(taille_1D, hmean, hdev, l, plot=False):

    h = np.random.randn(taille_1D) + hmean
    a = np.arange(taille_1D,dtype=float)
    b = np.ones((1,taille_1D))

    xi=np.kron(a,b.transpose())
    xj=np.kron(a,b.transpose()).transpose()

    c = (hdev**2)*np.exp(-((xi-xj)**2)/l**2)

    u,s,v = np.linalg.svd(c, full_matrices=False)
    d = np.dot(np.dot(u,np.diag(np.sqrt(s))),v)
    #cc=np.dot(d,d.transpose())

    hh = hmean + np.dot(d,h)

    if plot:
        plt.plot(hh)
        plt.show()


def generate_2d_profile_gaussian(taille_2D, hmean, size_filter, hdev, fact_echelle, plot=False):

    if size_filter == "Default":
        a = int(taille_2D[0] * fact_echelle)
        b = int(taille_2D[1] * fact_echelle)
        size_filter=[a,b]
    h = (np.random.rand(taille_2D[0],taille_2D[1])-0.5)

    filter = gauss_filter(size_filter)
    
    h_corr = signal.convolve(h,filter,mode='same') + hmean
    h_corr = h_corr * hdev/((abs(h_corr.min()) + h_corr.max())/2.)
    
    if plot:
        plt.figure()
        plt.imshow(h_corr)
        plt.colorbar()
        plt.show()
        
    return h_corr


def gauss_filter(size_filter, plot=False):

    sizex = size_filter[0]
    sizey = size_filter[1]

    x, y = np.mgrid[-sizex:sizex+1, -sizey:sizey+1]
    g = np.exp(-0.333*(x**2/float(sizex)+y**2/float(sizey)))
    return g/g.sum()

'''
def generate_2d_profile_2nd_order_list(size_x, size_y, x0, y0, a0, b0, c0, d0, e0, f0, plot=False):
     
    a = np.arange(size_x,dtype=float)
    b = np.ones((1,size_x))
    x=np.kron(a,b.transpose())
    
    a = np.arange(size_x,dtype=float)
    b = np.ones((1,size_x))
    y=np.kron(a,b.transpose()).transpose()
    
    X = x-x0
    Y = y-y0
    
    h = a0*X**2 + b0*Y**2 + c0*X + d0*Y + e0*X*Y + f0 
    
    if plot:
        print(h)
        plt.figure()
        plt.imshow(h)
        plt.colorbar()
        plt.show()
         
    return h
'''

def generate_2d_profile_2nd_order_list(x0, y0, x, y, a0, b0, c0, d0, e0, f0):

    X = x-x0
    Y = y-y0
    
    h = a0*X**2 + b0*Y**2 + c0*X + d0*Y + e0*X*Y + f0 
    
    return h


#############################
    

if __name__ == '__main__':
   
    h = generate_2d_profile_2nd_order_list(100,100,20.,20.,1.e-3,1.e-3,1.e-3,1.e-3,0.,5.)
    print(h)
