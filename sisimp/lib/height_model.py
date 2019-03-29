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

    hh = hmean + np.dot(d,h)

    if plot:
        plt.plot(hh)
        plt.show()

def generate_2d_profile_gaussian(dlat, latmin, latmax, dlon, lonmin, lonmax, height_model_stdv, plot=False, lcorr = 500, seed = None):
        
    Nx = int((latmax-latmin)/dlat)
    Ny = int((lonmax-lonmin)/dlon)
    lcorr = int(lcorr)
    
    lx = lcorr*dlat
    ly = lcorr*dlon
    
    Nx0= Nx
    Ny0= Ny
    
    # Need to increase image size due to larger lcorr compare to image shape
    if Nx < lx+10:
        Nx = lx+10
    if Ny < ly+10:
        Ny = ly+10           
 
    ### zero-padding not seems to improve computation time, TBC
    #~ nx_pad = 2**(int(np.log2(Nx - 1)) + 1)
    #~ ny_pad = 2**(int(np.log2(Ny - 1)) + 1)
    ###    
    
    nx_pad = Nx
    ny_pad = Ny
         
    if seed is not None:
        np.random.seed(int(seed))
    kx = np.fft.fftfreq(nx_pad, d=dlat)
    ky = np.fft.rfftfreq(ny_pad, d=dlon)
    
    ## TBD : Can be optimize : can be calculated only on non-zero value (after filtering by filterkx, filterky to avoid useless gaussian computation)
    hij_real = np.random.normal(0., height_model_stdv/np.sqrt(2), (len(kx),len(ky)))
    hij_imag = np.random.normal(0., height_model_stdv/np.sqrt(2), (len(kx),len(ky)))

    hij = hij_real + 1j*hij_imag
    
    filterkx = np.where(np.abs(kx)>  1/lx, 0, 1)
    filterky = np.where(np.abs(ky)>  1/ly, 0, 1)
            
    kxv, kyv = np.meshgrid(filterky,filterkx)
        
    hij = hij*kxv*kyv

    ## TBD : Can be optimize : check if "manual fft" only on non-zero values and only on 0:Nx0, 0:Ny0 interval is faster than that direct approach
    h_corr = np.sqrt((nx_pad*ny_pad/(4*dlon*dlat/lx/ly)))*np.fft.irfft2(hij, s=(nx_pad,ny_pad))

    h_corr = h_corr[0:Nx0, 0:Ny0]
        
    if plot:
        plt.figure()
        plt.imshow(h_corr)
        plt.show()
    
    return h_corr






def generate_2d_profile_gaussian_old(taille_2D, hmean, size_filter, hdev, fact_echelle, plot=False, seed = None):

    if size_filter == "Default":
        a = int(taille_2D[0] * fact_echelle)
        b = int(taille_2D[1] * fact_echelle)
        size_filter=[a,b]
    np.random.seed(seed)
    h = (np.random.rand(taille_2D[0],taille_2D[1])-0.5)

    filter = gauss_filter(size_filter)
    
    h_corr = signal.convolve(h,filter,mode='same') + hmean
    h_corr = 2*h_corr/(np.abs(np.amin(h_corr))+np.abs(np.amax(h_corr))) * hdev
    
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



def generate_2d_profile_2nd_order_list(x0, y0, x, y, a0, b0, c0, d0, e0, f0):

    X = x-x0
    Y = y-y0
    
    h = a0*X**2 + b0*Y**2 + c0*X + d0*Y + e0*X*Y + f0 
    
    return h


#############################
    

if __name__ == '__main__':
   
    h = generate_2d_profile_2nd_order_list(100,100,20.,20.,1.e-3,1.e-3,1.e-3,1.e-3,0.,5.)
    print(h)
