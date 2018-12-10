'''
Copyright (c) 2018 CNES. All rights reserved.
'''

from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import numpy as np
import scipy.optimize as opt
from math import atan2
import pyproj


GEN_RAD_EARTH_EQ = 6378137.0 # Radius of the Earth model (WGS84 ellipsoid) at the equator
GEN_RAD_EARTH_POLE = 6356752.31425 # Radius of the Earth model at the pole
GEN_RANGE_SPACING = 0.75 # Range spacing of SWOT
R_e = GEN_RAD_EARTH_EQ
R_p = GEN_RAD_EARTH_POLE

def normalize(vect):
    norm = np.sqrt(np.dot(vect,vect))
    return vect/norm #don t use me with np.zeroes(#)

def normalize_vect(a, axis=-1, order=2):
    l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
    l2[l2==0] = 1
    return a / np.expand_dims(l2, axis)


def height_fast(p):
    # minimal and fast height-only function using the analytical method
    # there is code duplication from the more generic function project_on_ellipsoid
    # but this avoids all if statements as well as computing alpha_0
    # The gain in time is only a few percent, but since this function will be used a lot
    # let's make it as fast as possible...
    
    x, y, z = p
    d = np.sqrt(x**2+y**2)
    alpha = np.sqrt(d**2/R_e**2 + z**2/R_p**2) #homotethy factor for an ellipsoid going through p_mu
    
    c_th = z/(alpha*R_p)
    s_th = d/(alpha*R_e)
    
    terme_xy = R_e*s_th-d
    terme_z = R_p*c_th-z
    
    G = terme_xy**2+terme_z**2
    dG_dth = 2*(R_e*c_th*terme_xy - R_p*s_th*terme_z)
    d2G_dth2 = 2*((R_e*c_th)**2+(R_p*s_th)**2-R_e*s_th*terme_xy-R_p*c_th*terme_z) # this is ~R_e**2>>1, so never 0
    delta_theta = -dG_dth/d2G_dth2
    A = G - dG_dth**2/d2G_dth2/2
    h_mu = np.sign(alpha-1.)*np.sqrt(A)
    return h_mu

def pointcloud_height_geoloc(p_noisy, h_noisy,
                                               s, vs, 
                                        R_target, h_target, 
                                        recompute_Doppler=False, 
                                        #if False uses 0 Doppler, else the value computed from p, s, vs
                                        recompute_R=False,
                                        verbose = False,
                                        max_iter_grad=1, height_goal = 1.e-3, safe_flag=True):
    
    # - for gradient only, set safe_flag to False
    # - for bounded only, set max_iter_grad to 0
    # - otherwise, will do gradient first (with up to max_iter_grad) and if height_goal is not met, will do bounded
    # Note that max_iter_grad=0 and safe_flag=False will not just return p_noisy. An initial guess is computed which
    # typically brings the height close to h_target by a couple of tens of centimeters at most (usually much less!)
    
    s_hat = normalize(s)

    # let's define first an adapted orthonormal basis
    v_hat = normalize(vs) #points along track
    u_hat = normalize(np.cross(s_hat, v_hat)) #s and v are almost ortho; no danger of colinear (points roughly cross track)
    w_hat = np.cross(u_hat, v_hat) # points roughly to nadir
    
    r_s = np.linalg.norm(s)
    delta = np.dot(p_noisy-s,v_hat) if recompute_Doppler else 0.
    R = np.linalg.norm(p_noisy-s) if recompute_R else R_target
    r_eq = np.sqrt(R**2-delta**2)
    def p_of_mu(mu, return_trigo=False):
        sin_mu = np.sin(mu)
        cos_mu = np.sqrt(1-sin_mu**2) #with our conventions, mu is small and cos(mu)>0
        p_mu = s + r_eq*(cos_mu*w_hat + sin_mu*u_hat) + delta*v_hat
        if return_trigo:
            return p_mu, cos_mu, sin_mu
        else:
            return p_mu
    
    #initial guess for mu uses the input point as a first start
    mu_0 = atan2(np.dot(p_noisy-s,u_hat),np.dot(p_noisy-s,w_hat))

    #tweak the intial guess: move by dh in the w_hat direction, taking the spherical earth correction factor into account
    R_sph_earth = (R_e+R_p)/2.#crude approx but ok since we are computing an order 1 correction
    Alt_sph = r_s-R_sph_earth
    fact_sph = R_sph_earth/(R_sph_earth+Alt_sph)
    mu_0 += fact_sph*(h_target-h_noisy)/r_eq/np.sin(mu_0)

    # findroot h**2(mu)==h_target using the gradient method (and analytical mu derivatives)
    # note that this can converge poorly if we aim too close to h=0
    # see the notes attached. (or even converge to the wrong solution; we have taken a square...)
    # After a chosen nb of iterations, we check that we got our solution to the prescribed accuracy
    # and if it is not the case, we resort to a safer method (if safe flag on)
    iter_grad = 0
    mu = mu_0
    while iter_grad<=max_iter_grad:
        p_mu, cos_mu, sin_mu = p_of_mu(mu, return_trigo=True)
        
        #computation of A (= height**2 at mu)
        # this is a bit of code duplication wrt height_fast but we will reuse most of the internal variables later...
        z_mu, d_mu = p_mu[2], np.sqrt(p_mu[0]**2+p_mu[1]**2)
        alpha = np.sqrt(d_mu**2/R_e**2 + z_mu**2/R_p**2)
        c_th, s_th = z_mu/(alpha*R_p), d_mu/(alpha*R_e)
        terme_xy, terme_z = R_e*s_th-d_mu, R_p*c_th-z_mu
        G = terme_xy**2+terme_z**2
        dG_dth = 2*(R_e*c_th*terme_xy - R_p*s_th*terme_z)
        d2G_dth2 = 2*((R_e*c_th)**2+(R_p*s_th)**2-R_e*s_th*terme_xy-R_p*c_th*terme_z) # this is ~R_e**2>>1, so never 0
        A = G - dG_dth**2/d2G_dth2/2
        h_mu = np.sign(alpha-1.)*np.sqrt(A)

        h_error_mu = np.abs(h_mu-h_target)
        if iter_grad==0:
            h_error_mu_0 = h_error_mu # save for later
              
        # the second condition is mostly cosmetic and also saves us from the last unused computation of B
        if (h_error_mu<height_goal)|(iter_grad==max_iter_grad):
            break
        
        #computation of B (=derivative wrt mu of the height**2 at mu)
        dp_mu_dmu = r_eq*(-sin_mu*w_hat + cos_mu*u_hat)
        dz_dmu = dp_mu_dmu[2]
        dd_dmu = (dp_mu_dmu[0]*p_mu[0]+dp_mu_dmu[1]*p_mu[1])/d_mu
        dG_dmu = -2*(dd_dmu*terme_xy+dz_dmu*terme_z)
        d2G_dthdmu = 2*(-dd_dmu*R_e*c_th+dz_dmu*R_p*s_th)
        d3G_dth2dmu = 2*(dd_dmu*R_e*s_th+dz_dmu*R_p*c_th)
        B = dG_dmu - d2G_dthdmu*dG_dth/d2G_dth2 + dG_dth**2*d3G_dth2dmu/d2G_dth2**2/2
        
        delta_mu = (h_target**2-A)/B

        # just for debugging purposes
        if verbose:
            print('gradient method, iteration ',iter_grad)
            print('error before updating mu:', h_error_mu)
            print('h2, A, B', h_target**2, A, B)
            print('termes de B',dG_dth**2*d3G_dth2dmu/d2G_dth2**2/2, d2G_dthdmu*dG_dth/d2G_dth2, dG_dmu)
            delta_mu_tab = np.linspace(-3*np.abs(delta_mu),3*np.abs(delta_mu),100)
            Ghat = lambda mu: height(p_of_mu(mu))**2
            Gtilde_lin = lambda dmu: A + B*dmu
            Ghat_tab = np.array([Ghat(mu+ddd) for ddd in delta_mu_tab])  
            Gtilde_lin_tab = np.array([Gtilde_lin(ddd) for ddd in delta_mu_tab])
            plt.figure(figsize=(10,5))
            plt.plot(delta_mu_tab, Ghat_tab, c='r', label='h2(mu)')
            plt.plot(delta_mu_tab, Gtilde_lin_tab, c='b', label='h2(mu) linearized', ls='--')
            plt.axvline(0, c='k')
            plt.axvline(delta_mu, c='gray', ls='--')
            plt.axhline(h_target**2, c='k', ls='--')
            plt.legend()
            #plt.ylim([0,.2**2])
        
        mu += delta_mu # note that if this is the last iteration and we got that far (error criterion not met), then this is never applied
        iter_grad += 1
       
    
    # Now that we are outside of the loop, we check whether we got the right (signed!) height. If not, safe method
    # namely a numerical bounded minimization of (h(mu)-h_target)**2 (which uses numerical derivatives)
    # in a resticted segment around the intial guess mu_0 whose length is chosen so as to cover approximately 
    # a height variation of dh_max
    nfev_minimize_scalar = 0 #to store the number of function evaluations during the numerical optimization
    if (h_error_mu>height_goal)&safe_flag:
        dh_max = 2*h_error_mu_0
        dmu_max = np.abs(fact_sph*dh_max/r_eq/np.sin(mu_0))
        
        def hofmu(mu):
            p_mu = p_of_mu(mu)
            h_mu = height_fast(p_mu)
            return h_mu, p_mu
        
        def herror2ofx(x): # cost function, rescale mu
            mu = mu_0 + x*dmu_max #moving by 1 in x moves by dmu_max in mu so *approximately* by dh_max in h
            h_error_square_mu = (hofmu(mu)[0]-h_target)**2
            return h_error_square_mu
    
        i_verbose = 3 if verbose else 0 
        sol = opt.minimize_scalar(herror2ofx, bounds=(-1.,1.),method='bounded',
                                  options={'xatol':height_goal/dh_max, 'disp':i_verbose})
        mu_sol = mu_0 + dmu_max*sol.x
        h_mu, p_mu = hofmu(mu_sol) #a bit of a waste of computation here for the sake of readability
        # we could just recompute p_mu and extract h from sol (but we then need to be careful about the sign etc...)
        nfev_minimize_scalar = sol.nfev
    
        # just for debugging purposes
        if verbose:
            print('The error was '+str(h_error_mu)+'after the gradient method. Switch to minimize_scalar method.')
            delta_mu_tab = np.linspace(-dmu_max,dmu_max,100)
            h_tab = np.array([hofmu(mu_0+ddd)[0] for ddd in delta_mu_tab])
            plt.figure(figsize=(10,5))
            plt.plot(delta_mu_tab, h_tab, c='r', label='h(mu)')
            plt.axvline(0, c='k')
            plt.axhline(h_target, c='k', ls='--')
            plt.axhline(hofmu(mu_0)[0], c='g', ls='--')
            plt.legend()
            #plt.ylim([0,.2**2])
            
    p_final_llh = convert_ecef2llh(p_mu[0], p_mu[1], p_mu[2], GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE)
    #~ p_noisy_llh = convert_ecef2llh(p_noisy[0], p_noisy[1], p_noisy[2], GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE)
    #~ print(p_noisy_llh[0], p_final_llh[0], p_noisy_llh[1], p_final_llh[1], p_noisy_llh[2], p_final_llh[2])
    
    return p_mu, p_final_llh, h_mu, (iter_grad,nfev_minimize_scalar)










def pointcloud_height_geoloc_vect(p_noisy, h_noisy,
                                               s, vs, 
                                        R_target, h_target, 
                                        recompute_Doppler=False, 
                                        #if False uses 0 Doppler, else the value computed from p, s, vs
                                        recompute_R=False,
                                        verbose = False,
                                        max_iter_grad=1, height_goal = 1.e-3, safe_flag=True):
    
    # - for gradient only, set safe_flag to False
    # - for bounded only, set max_iter_grad to 0
    # - otherwise, will do gradient first (with up to max_iter_grad) and if height_goal is not met, will do bounded
    # Note that max_iter_grad=0 and safe_flag=False will not just return p_noisy. An initial guess is computed which
    # typically brings the height close to h_target by a couple of tens of centimeters at most (usually much less!)
    
    s_hat = normalize_vect(s)

    # let's define first an adapted orthonormal basis
    v_hat = normalize_vect(vs) #points along track


    u_hat = normalize_vect(np.cross(s_hat, v_hat)) #s and v are almost ortho; no danger of colinear (points roughly cross track)
    w_hat = np.cross(u_hat, v_hat) # points roughly to nadir
    
    r_s = np.linalg.norm(s, axis=1)
    
    delta = np.einsum('ij,ij->i',p_noisy-s,v_hat) if recompute_Doppler else 0.
    
    R = np.linalg.norm(p_noisy-s, axis = 1) if recompute_R else R_target
    r_eq = np.sqrt(R**2-delta**2)
        
    def p_of_mu_vect(mu, return_trigo=False):
        sin_mu = np.sin(mu)
        cos_mu = np.sqrt(1-sin_mu**2) #with our conventions, mu is small and cos(mu)>0
        p_mu = s + np.einsum('i,ij->ij',r_eq,np.einsum('i,ij->ij',cos_mu,w_hat) + np.einsum('i,ij->ij',sin_mu,u_hat)) + np.einsum('i,ij->ij',delta,v_hat)
        if return_trigo:
            return p_mu, cos_mu, sin_mu
        else:
            return p_mu

    def p_of_mu(mu, indice, return_trigo=False):
        sin_mu = np.sin(mu)
        cos_mu = np.sqrt(1-sin_mu**2) #with our conventions, mu is small and cos(mu)>0
        p_mu = s[indice] + r_eq[indice]*(cos_mu*w_hat[indice] + sin_mu*u_hat[indice]) + delta[indice]*v_hat[indice]
        if return_trigo:
            return p_mu, cos_mu, sin_mu
        else:
            return p_mu
                
    #initial guess for mu uses the input point as a first start
    mu_0 = np.arctan2(np.einsum('ij,ij->i',p_noisy-s,u_hat),np.einsum('ij,ij->i',p_noisy-s,w_hat))
    #tweak the intial guess: move by dh in the w_hat direction, taking the spherical earth correction factor into account
    R_sph_earth = (R_e+R_p)/2.#crude approx but ok since we are computing an order 1 correction
    Alt_sph = r_s-R_sph_earth
    fact_sph = R_sph_earth/(R_sph_earth+Alt_sph)
    mu_0 += fact_sph*(h_target-h_noisy)/r_eq/np.sin(mu_0)

    # findroot h**2(mu)==h_target using the gradient method (and analytical mu derivatives)
    # note that this can converge poorly if we aim too close to h=0
    # see the notes attached. (or even converge to the wrong solution; we have taken a square...)
    # After a chosen nb of iterations, we check that we got our solution to the prescribed accuracy
    # and if it is not the case, we resort to a safer method (if safe flag on)
    iter_grad = 0
    mu = np.copy(mu_0)
    while iter_grad<=max_iter_grad:
        p_mu, cos_mu, sin_mu = p_of_mu_vect(mu, return_trigo=True)
        #computation of A (= height**2 at mu)
        # this is a bit of code duplication wrt height_fast but we will reuse most of the internal variables later...
        z_mu, d_mu = p_mu[:,2], np.sqrt(p_mu[:,0]**2+p_mu[:,1]**2)
        alpha = np.sqrt(d_mu**2/R_e**2 + z_mu**2/R_p**2)
        c_th, s_th = z_mu/(alpha*R_p), d_mu/(alpha*R_e)
        terme_xy, terme_z = R_e*s_th-d_mu, R_p*c_th-z_mu
        G = terme_xy**2+terme_z**2
        dG_dth = 2*(R_e*c_th*terme_xy - R_p*s_th*terme_z)
        d2G_dth2 = 2*((R_e*c_th)**2+(R_p*s_th)**2-R_e*s_th*terme_xy-R_p*c_th*terme_z) # this is ~R_e**2>>1, so never 0
        A = G - dG_dth**2/d2G_dth2/2
        h_mu = np.sign(alpha-1.)*np.sqrt(A)
        
        h_error_mu = np.abs(h_mu-h_target)
        if iter_grad==0:
            h_error_mu_0 = h_error_mu # save for later
              
        #~ # the second condition is mostly cosmetic and also saves us from the last unused computation of B
        #~ if (h_error_mu<height_goal)|(iter_grad==max_iter_grad):
            #~ break

        #computation of B (=derivative wrt mu of the height**2 at mu)
        dp_mu_dmu = np.einsum('i,ij->ij',r_eq,np.einsum('i,ij->ij',-sin_mu,w_hat) + np.einsum('i,ij->ij',cos_mu,u_hat))
        dz_dmu = dp_mu_dmu[:,2]
        dd_dmu = (dp_mu_dmu[:,0]*p_mu[:,0]+dp_mu_dmu[:,1]*p_mu[:,1])/d_mu
        dG_dmu = -2*(dd_dmu*terme_xy+dz_dmu*terme_z)
        d2G_dthdmu = 2*(-dd_dmu*R_e*c_th+dz_dmu*R_p*s_th)
        d3G_dth2dmu = 2*(dd_dmu*R_e*s_th+dz_dmu*R_p*c_th)
        B = dG_dmu - d2G_dthdmu*dG_dth/d2G_dth2 + dG_dth**2*d3G_dth2dmu/d2G_dth2**2/2
        
        delta_mu = (h_target**2-A)/B

        
        mu += delta_mu # note that if this is the last iteration and we got that far (error criterion not met), then this is never applied
        iter_grad += 1
       
    
    # Now that we are outside of the loop, we check whether we got the right (signed!) height. If not, safe method
    # namely a numerical bounded minimization of (h(mu)-h_target)**2 (which uses numerical derivatives)
    # in a resticted segment around the intial guess mu_0 whose length is chosen so as to cover approximately 
    # a height variation of dh_max
       
    
    nfev_minimize_scalar = 0 #to store the number of function evaluations during the numerical optimization
    

    safe_flag = True
    if safe_flag:
        indices = np.where(h_error_mu > height_goal)
        for i in indices[0]:
            
            dh_max = 2*h_error_mu_0[i]
            dmu_max = np.abs(fact_sph[i]*dh_max/r_eq[i]/np.sin(mu_0[i]))
            #~ print(mu_0[15])

            def hofmu(mu):
                p_mu = p_of_mu(mu, i)
                h_mu = height_fast(p_mu)
                #~ print(mu_0[15])

                return h_mu, p_mu
            
            def herror2ofx(x): # cost function, rescale mu
                mu = mu_0[i] + x*dmu_max #moving by 1 in x moves by dmu_max in mu so *approximately* by dh_max in h
                h_error_square_mu = (hofmu(mu)[0]-h_target[i])**2
                return h_error_square_mu
        
            i_verbose = 3 if verbose else 0 
            sol = opt.minimize_scalar(herror2ofx, bounds=(-1.,1.),method='bounded',
                                      options={'xatol':height_goal/dh_max, 'disp':i_verbose})
            mu_sol = mu_0[i] + dmu_max*sol.x
            h_mu[i], p_mu[i] = hofmu(mu_sol) #a bit of a waste of computation here for the sake of readability
            # we could just recompute p_mu and extract h from sol (but we then need to be careful about the sign etc...)
            nfev_minimize_scalar = sol.nfev
        

    
    p_final_llh = np.zeros([len(p_mu),3],np.double)
    

    #~ for i in range(len(p_mu)):
        #~ p_final_llh[i] = convert_ecef2llh(p_mu[i,0], p_mu[i,1], p_mu[i,2], GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE)
    
    

    def project_array(coordinates, srcp='geocent', dstp='latlon'):
        """
        Project a numpy (n,2) array in projection srcp to projection dstp
        Returns a numpy (n,2) array.
        """
        p1 = pyproj.Proj(proj=srcp, datum='WGS84')
        p2 = pyproj.Proj(proj=dstp, datum='WGS84')
        fx, fy, fz = pyproj.transform(p1, p2, coordinates[:,0], coordinates[:,1], coordinates[:,2])
        # Re-create (n,2) coordinates
        # Inversion of lat and lon !
        return np.dstack([fy, fx, fz])[0]

    p_final_llh = project_array(p_mu)
    
    
    return p_mu, p_final_llh, h_mu, (iter_grad,nfev_minimize_scalar)

















########################## Old  version   ####################################




def convert_llh2ecef(lat, lon, height, rad_e=R_e, rad_p=R_p):
    """
    Converting from geodetic coordinates (lat, lon, height) to cartesian coordinates (x, y, z)

    Args:
        lat(float): latitude of the record

        lon(float): longitude of the record

        height(float): height of the record

        rad_e(float): radius of the Earth model (WGS84 ellipsoid) at the equator

        rad_p(float): radius of the Earth model at the pole

    Returns:
        tuple. Contain : base axes (x, y and z)
    """
    lat = lat*np.pi/180.
    lon = lon*np.pi/180.

    e = np.sqrt(rad_e*rad_e - rad_p*rad_p)/rad_e
    Rn = rad_e/np.sqrt(1-e*e*np.sin(lat)*np.sin(lat))

    axe_x = (Rn + height) * np.cos(lat) * np.cos(lon)
    axe_y = (Rn + height) * np.cos(lat) * np.sin(lon)
    axe_z = (Rn * (1.-e*e) + height) * np.sin(lat)

    return axe_x, axe_y, axe_z

def convert_ecef2llh(axe_x, axe_y, axe_z, rad_e=R_e, rad_p=R_p):
    """
    Converting from cartesian coordinates (x, y, z) to geodetic coordinates (lat, lon, height)

    Args:
        axe_x(float): base axe 'x' of the record

        axe_y(float): base axe 'y' of the record

        axe_z(float): base axe 'z' of the record

        rad_e(float): radius of the Earth model (WGS84 ellipsoid) at the equator

        rad_p(float): radius of the Earth model at the pole

    Returns:
        tuple. Contain : coordinates (latitude, longitude, height)
    """
    e = np.sqrt(rad_e*rad_e - rad_p*rad_p)/rad_e

    lat1 = 0.
    lon = np.arctan2(axe_y, axe_x)


    P = np.sqrt(axe_x*axe_x + axe_y*axe_y)
    lat0 = np.arctan(axe_z/(P*(1.-e*e)))
    N0 = rad_e/np.sqrt(1.-e*e*np.sin(lat0)*np.sin(lat0))

    val_acc = 1.
    acc = 1.e-20

    while val_acc > acc:
        if axe_z == 0.:
            lat1 = 0.
            N1 = rad_e/np.sqrt(1.-e*e*np.sin(lat1)*np.sin(lat1))
            val_acc = np.abs(np.abs(lat0) - np.abs(lat1))
            lat0 = lat1
            N0 = N1
        else:
            lat1 = np.arctan((axe_z/P)*(1.+(e*e*N0*np.sin(lat0)/axe_z)))
            N1 = rad_e/np.sqrt(1.-e*e*np.sin(lat1)*np.sin(lat1))
            val_acc = np.abs(np.abs(lat0) - np.abs(lat1))
            lat0 = lat1
            N0 = N1

    height = P/np.cos(lat1) - N1

    lat = lat1 * 180./np.pi
    lon = lon * 180./np.pi

    return lat, lon, height
    
def proj_ellipsoid(R_e, R_p, p, verbose=False, numerical_min=True):
    x, y, z = p
    phi = atan2(y, x)
    d = np.sqrt(x**2+y**2)
    #~ r = np.linalg.norm(p)
    #~ R_e2 = R_e**2
    #~ R_p2 = R_p**2
    alpha = np.sqrt(d**2/R_e**2 + z**2/R_p**2)

    #homothetie between ellipsoide and ellipse on p
    #Development around theta_0 (parametrization of an homotetic ellipse (x2+yp2)/(alpha*Re)2 + (zp2)/(alpha*Rp)2 = 1
    theta_0 = np.arccos(z/(alpha*R_p))

    if numerical_min:
        ## Minimization of (xp-x)2+(yp-y)2+(zp-z)2 with xp, yp, zp are the p coordinates and z,y,z are the coordinates of the ellipsoide define by
        ## x=Re.sin(theta)cos(phi)
        ## y=Re.sin(theta)sin(phi)
        ## x=Rp.cos(theta)
        ### Minimal for phi = atan(yp,xp)
        ### After development :

        r = np.linalg.norm(p)

        def G(theta):
            terme_cst = r**2+R_p**2
            terme_sin2 = (R_e**2-R_p**2)*np.sin(theta)**2
            terme_sin = -2*R_e*d*np.sin(theta)
            terme_cos = -2*R_p*z*np.cos(theta)
            return terme_cst+terme_sin2+terme_cos+terme_sin

        #numerical minimization
        ### to be parametrized
        solve_bisec = opt.minimize_scalar(G, bounds=(.9999*theta_0,1.0001*theta_0), method='bounded')
        theta_sol = solve_bisec.x

    else:
        #analytical minimization (linearization use me only with points close to the ellipsoid)
        #marginally faster
        ct = np.cos(theta_0)
        st = np.sin(theta_0)
        delta_theta_num = -1*(ct*st*(R_e**2-R_p**2) - R_e*d*ct + R_p*z*st)
        delta_theta_den = ((ct*ct-st*st)*(R_e**2-R_p**2) + R_e*d*st + R_p*z*ct)
        delta_theta = delta_theta_num/delta_theta_den
        theta_sol = theta_0+delta_theta

    if verbose:
        import matplotlib.pyplot as plt
        check_dtetha_vals = np.linspace(-1e-7,1e-7,51)
        check_vals = np.array([np.sqrt(G(theta_0+dtheta)) for dtheta in check_dtetha_vals])
        plt.plot(check_dtetha_vals, check_vals)
        plt.axvline(theta_sol-theta_0, c='g')

    x_sol = R_e*np.sin(theta_sol)*np.cos(phi)
    y_sol = R_e*np.sin(theta_sol)*np.sin(phi)
    z_sol = R_p*np.cos(theta_sol)
    return np.array([x_sol, y_sol, z_sol])

def height_above_ellipsoid(p, compute_proj=False, numerical_min=False, verbose=False):
    x, y, z = p
    d = np.sqrt(x**2+y**2)
    alpha = np.sqrt(d**2/R_e**2 + z**2/R_p**2) #homotethy factor for an ellipsoid going through p_mu
    theta_0 = np.arccos(z/(alpha*R_p)) #search/expand around this value
    if numerical_min:
        #numerical minimization
        def G(theta):
            terme_xy = R_e*np.sin(theta)-d #precision 1e-9 m
            terme_z = R_p*np.cos(theta)-z #precision 1e-9 m
            return terme_xy**2+terme_z**2
        solve_bisec = opt.minimize_scalar(G, bounds=(.9999*theta_0,1.0001*theta_0), method='bounded',
                                          options={'disp':0, 'xatol':1e-12})
        theta_sol = solve_bisec.x
        A = solve_bisec.fun
    else:
        #analytical minimization by linearizing in theta
        #c_th, s_th = np.cos(theta_0), np.sin(theta_0) #in fact, we know them as functions of z, d and alpha
        c_th, s_th = z/(alpha*R_p), d/(alpha*R_e)
        # also note that in fact theta_0 is only necessary to get the projection but not the height in this approach
        # we could do with only its cos and sin. We would avoid one arccos. 
        # We keep it for readability and implement a separate minimal and fast version for height only
        terme_xy = R_e*s_th-d
        terme_z = R_p*c_th-z
        G = terme_xy**2+terme_z**2
        dG_dth = 2*(R_e*c_th*terme_xy - R_p*s_th*terme_z)
        d2G_dth2 = 2*((R_e*c_th)**2+(R_p*s_th)**2-R_e*s_th*terme_xy-R_p*c_th*terme_z) # this is ~R_e**2>>1, so never 0
        delta_theta = -dG_dth/d2G_dth2
        theta_sol = theta_0+delta_theta
        A = G - dG_dth**2/d2G_dth2/2
        # note: now sqrt(G(theta_0 + delta_theta)) and the height computed taking the norm of the difference of coords should agree up to 1e-9
        # the linear approx is extremely good to determine delta_theta (i.e. our ini guess is good enough that G is very quadratic )
        # G(theta_0)-delta_theta_num**2/delta_theta_den and G(theta_sol) are also extremely close

    # once we computed theta_sol, h is almost for free so we compute it no matter what.
    # by contrast, the coordinates of the projection require a atan and 4 sines/cosines
    # we only compute it on demand
    h_mu = np.sign(alpha-1.)*np.sqrt(A)
    if compute_proj:
        phi = atan2(y, x)
        cos_theta = np.cos(theta_sol)
        sin_theta = np.sqrt(1-cos_theta**2) #always positive
        x_sol = R_e*sin_theta*np.cos(phi)
        y_sol = R_e*sin_theta*np.sin(phi)
        z_sol = R_p*cos_theta
        return h_mu, np.array([x_sol, y_sol, z_sol])
    else:
        return h_mu

def pointcloud_height_geoloc_old(p_noisy, s, vs, R_target, h_target, h_brut,
                             recompute_Doppler = True,
                             recompute_R = True,
                             verbose = False,
                             maxiter = 2,
                             tolerance_height = 0.01):
    """
    Height constrained improved geolocation

    p_noisy:           initial position of the noisy pixel in x,y,z
    s:                 sensor position corresponding to the pixel
    vs:                sensor velocity corresponding the the pixel
    R_target:          target range used for geolocation
    h_target:          target height used for geolocation
    h_brut:            initial (noisy) height
    recompute_Doppler: if True recompute doppler from p, s, vs, else use 0 doppler
    recompute_R:       if True recompute R_target using point position, else use R_target
    verbose:           display debug info
    maxiter:           maximum number of iterations
    tolerance_height:  stopping condition on height error (m)

    Returns:
        p_final:       Pixel position in x,y,z
        p_final_llh:   Pixel position in lat, lon, height
    """

    # Height threshold to fallback to the numerical solver for precision problem when target point is too close to the ellipsoid
    # hardcoded because this is a temporary fix
    threshold_target_ellipsoid = 0.5
    if abs(h_target) < threshold_target_ellipsoid:
                                                                   
        p_final_a, p_final_llh_a = pointcloud_height_geoloc(p_noisy, s, vs,
                                                            R_target, h_target+2., h_brut,
                                                            recompute_Doppler=recompute_Doppler,
                                                            recompute_R=recompute_R,
                                                            verbose=verbose)
                                                                                                                                                                  
        p_final_b, p_final_llh_b = pointcloud_height_geoloc(p_noisy, s, vs,
                                                            R_target, h_target-2., h_brut,
                                                            recompute_Doppler=recompute_Doppler,
                                                            recompute_R=recompute_R,
                                                            verbose=verbose) 

        
        s_hat = normalize(s)
        v_hat = normalize(vs) #points along track
        u_hat = normalize(np.cross(s_hat, v_hat)) #s and v are almost ortho; no danger of colinear (points roughly cross track)
        w_hat = np.cross(u_hat, v_hat) # points roughly to nadir

        r_s = np.linalg.norm(s)

        delta = np.dot(p_noisy-s,v_hat) if recompute_Doppler else 0.
        R = np.linalg.norm(p_noisy-s) if recompute_R else R_target

        mu_a = atan2(np.dot(p_final_a-s,u_hat),np.dot(p_final_a-s,w_hat))
        mu_b = atan2(np.dot(p_final_b-s,u_hat),np.dot(p_final_b-s,w_hat))
        
        delta_mu = (mu_a+mu_b)/2.
        
        r_eq = np.sqrt(R**2-delta**2)

        p_final = s + r_eq*(np.cos(delta_mu)*w_hat + np.sin(delta_mu)*u_hat) + delta*v_hat
        p_final_llh = convert_ecef2llh(p_final[0], p_final[1], p_final[2], GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE)
        return p_final, p_final_llh

    # Iterate
    i = 0
    iterate = True
    while iterate:
        i += 1
        p_final, p_final_llh = pointcloud_ellipsoidgeoloc_improved_Taylortheta2mu1(p_noisy, s, vs,
                                                                                   R_target, h_target, h_brut,
                                                                                   recompute_Doppler=recompute_Doppler,
                                                                                   recompute_R=recompute_R,
                                                                                   verbose=verbose)

        # Update for next iteration
        p_noisy = p_final
        h_brut = p_final_llh[2]

        height_error = np.abs(p_final_llh[2] - h_target)
        iterate = (height_error > tolerance_height and i < maxiter)

    return p_final, p_final_llh

def pointcloud_ellipsoidgeoloc_improved_Taylortheta2mu1(p_noisy, s, vs,
                                        R_target, h_target, h_brut,
                                        recompute_Doppler=True,
                                        recompute_R=True,
                                        verbose = False):

    s_hat = normalize(s)
    # let's define first an adapted orthonormal basis
    v_hat = normalize(vs) #points along track
    u_hat = normalize(np.cross(s_hat, v_hat)) #s and v are almost ortho; no danger of colinear (points roughly cross track)
    w_hat = np.cross(u_hat, v_hat) # points roughly to nadir

    r_s = np.linalg.norm(s)

    delta = np.dot(p_noisy-s,v_hat) if recompute_Doppler else 0.
    R = np.linalg.norm(p_noisy-s) if recompute_R else R_target

    #OLDlet's choose the sign of mu (discriminates between both swaths) based on the input point
    # mu_input_sign = np.sign(np.dot(p_noisy-s,u_hat))

    #initial guess for mu uses the input point for now
    mu = atan2(np.dot(p_noisy-s,u_hat),np.dot(p_noisy-s,w_hat))

    h_noisy = height_above_ellispoid(p_noisy, numerical_min=False)
    r_eq = np.sqrt(R**2-delta**2)
    mu += (h_target-h_brut)/r_eq/np.sin(mu) #tweak the intial guess: move by dh in the w_hat direction

    sinmu = np.sin(mu)
    cosmu = np.cos(mu)

    p_mu = s + r_eq*(cosmu*w_hat + sinmu*u_hat) + delta*v_hat
    dp_mu_dmu = r_eq*(-sinmu*w_hat + cosmu*u_hat)

    r2_mu = np.dot(p_mu,p_mu)
    z_mu = p_mu[2]
    d_mu = np.sqrt(r2_mu-z_mu*z_mu)

    dr2_dmu = 2*np.dot(p_mu,dp_mu_dmu)
    dz_dmu = dp_mu_dmu[2]
    dd_dmu = (dr2_dmu/2.-z_mu*dz_dmu)/d_mu


    #initial guess for theta
    alpha = np.sqrt(d_mu**2/R_e**2 + z_mu**2/R_p**2) #homotethy factor for an ellipsoid going through p_noisy
    theta_0 = np.arccos(z_mu/(alpha*R_p))

    DR2 = R_e**2-R_p**2
    s_th, c_th = np.sin(theta_0), np.cos(theta_0)
    s_th2 , c_th2 = s_th**2, c_th**2
    G = r2_mu+R_p**2 + DR2*s_th2 - 2*R_e*d_mu*s_th - 2*R_p*z_mu*c_th
    dG_dth = 2*DR2*s_th*c_th - 2*R_e*d_mu*c_th + 2*R_p*z_mu*s_th
    d2G_dth2 = 2*DR2*(c_th2-s_th2) + 2*R_e*d_mu*s_th + 2*R_p*z_mu*c_th
    dG_dmu = dr2_dmu - 2*R_e*dd_dmu*s_th -2*R_p*dz_dmu*c_th
    d2G_dthdmu = -2*R_e*dd_dmu*c_th + 2*R_p*dz_dmu*s_th
    d3G_dth2dmu = 2*R_e*dd_dmu*s_th + 2*R_p*dz_dmu*c_th

    A = G - dG_dth**2/d2G_dth2/2
    B = dG_dmu - d2G_dthdmu*dG_dth/d2G_dth2 + dG_dth**2*d3G_dth2dmu/d2G_dth2**2/2

    delta_mu = (h_target**2-A)/B

    p_final = s + r_eq*(np.cos(mu+delta_mu)*w_hat + np.sin(mu+delta_mu)*u_hat) + delta*v_hat

    # result in lat, lon, height
    p_final_llh = convert_ecef2llh(p_final[0], p_final[1], p_final[2], GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE)
    
    return p_final, np.array(p_final_llh)

