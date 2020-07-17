# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2019/05/17:version initiale.
# VERSION:2.0.0:DM:#91:2020/07/03:Poursuite industrialisation
# FIN-HISTORIQUE
# ======================================================
'''
.. module:: geoloc.py

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


'''

from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import numpy 
import scipy.optimize as opt
import pyproj as pyproj
import logging
from cnes.common.lib.my_variables import GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE


def normalize_vect(a, axis=-1, order=2):
    """
    Normalize n vectors in multi-dimensional arrays

    :param a: vector array
    :type a: numpy 2D-array of size (n, 3=x/y/z) with n the number of vector ([[x1,y1,z1],...[xn,yn,zn]])
             or numpy 2D-array of size (3=x/y/z, n) with n the number of vector ([[x1,...,xn],[y1,...,Yn][z1,...,zn]])
    :param axis: axis of normalization
    :type axis: integer (-1 for a numpy 2D-array of size (n, 3=x/y/z) or 0 for a numpy 2D-array of size (3=x/y/z, n) )
    :param order: define the norm (order=2 for the Euclidean norm)
    :type axis: integer

    :return normalized vector array
    :rtype: numpy 2D-array of float
    """

    # 0 -Compute the norm
    l2 = numpy.atleast_1d(numpy.linalg.norm(a, order, axis))

    # 1 -replace the 0 by 1 to avoid div by zero
    l2[l2 == 0] = 1

    # 3 - normalize a
    return a / numpy.expand_dims(l2, axis)


def height_fast(p):
    """
    Return an approximate height of p point using the analytical method

    :param p: point
    :type p: 1D-array of float: [x|y|z] cartesian coordinates

    :return: approximate height
    :rtype: float
    """
    # minimal and fast height-only function using the analytical method
    # there is code duplication from the more generic function project_on_ellipsoid
    # but this avoids all if statements as well as computing alpha_0
    # The gain in time is only a few percent, but since this function will be used a lot


    # 0 - resolution of G_tilde(mu) = A + d_mu*B = h_mu**2 considering B is equal to 0

    # 1 - variables initialization
    x, y, z = p
    d = numpy.sqrt(x ** 2 + y ** 2)
    # compute the homothety factor for an ellipsoid going through p_mu
    alpha = numpy.sqrt(d ** 2 / GEN_RAD_EARTH_EQ ** 2 + z ** 2 / GEN_RAD_EARTH_POLE ** 2)
    # compute cosinus theta
    c_th = z / (alpha * GEN_RAD_EARTH_POLE)
    # compute sinus theta
    s_th = d / (alpha * GEN_RAD_EARTH_EQ)
    terme_xy = GEN_RAD_EARTH_EQ * s_th - d
    terme_z = GEN_RAD_EARTH_POLE * c_th - z

    # 2 -computation of G
    G = terme_xy ** 2 + terme_z ** 2
    # derivative of G
    dG_dth = 2 * (GEN_RAD_EARTH_EQ * c_th * terme_xy - GEN_RAD_EARTH_POLE * s_th * terme_z)
    # second derivative of G
    d2G_dth2 = 2 * ((GEN_RAD_EARTH_EQ * c_th) ** 2 + (
        GEN_RAD_EARTH_POLE * s_th) ** 2 - GEN_RAD_EARTH_EQ * s_th * terme_xy - GEN_RAD_EARTH_POLE * c_th * terme_z)  # this is ~GEN_RAD_EARTH_EQ**2>>1, so never 0

    # 3 - computation of A
    A = G - dG_dth ** 2 / d2G_dth2 / 2

    # 4 - computation of h_mu
    h_mu = numpy.sign(alpha - 1.) * numpy.sqrt(A)
    return h_mu


def project_array(coordinates, srcp='geocent', dstp='latlon'):
    """
    Project a numpy (n,3) array in projection srcp to projection dstp
    Returns a numpy (n,3) array.

    :param coordinates: coordinates of points
    :type coordinates: 2D-array of float; size (3=x/y/z, nb_points)
    :param srcp: input system of coordinates
    :type srcp: string
    :param dstp: output system of coordinates
    :type dstp: string

    :return: dstp projection
    :rtype: float

    """
    input_proj = pyproj.Proj(proj=srcp, datum='WGS84')
    output_proj = pyproj.Proj(proj=dstp, datum='WGS84')
    fx, fy, fz = pyproj.transform(input_proj, output_proj, coordinates[:, 0], coordinates[:, 1], coordinates[:, 2])
    # Re-create (n,3) coordinates
    # Inversion of lat and lon !
    return numpy.dstack([fy, fx, fz])[0]


def p_of_mu_vect(mu, s, r_eq, delta, u_hat, v_hat, w_hat):
    """
    Compute the p points,the sinus associated and the cosinus associated for analytic problem
    Return the p points,the sinus associated and the cosinus associated

    :param mu: mu angles of analytic problem
    :type mu: numpy 1D-array of float;size (nb_points)
    :param s: position of the sensor in cartesian coordinates
    :type s: numpy 2D-array of float; size (3=x/y/z, nb_points)
    :param r_eq: range of problem
    :type r_eq: numpy 1D-array of float;size (nb_points)
    :param delta: correction factor of v in (u,v,w)basis problem
    :type delta: numpy 1D-array of float;size (nb_points)
    :param u_hat: u basis vector of (u,v,w)basis problem
    :type u_hat: numpy 2D-array of float; size (3=x/y/z, nb_points)
    :param v_hat: v basis vector of (u,v,w)basis problem
    :type u_hat: numpy 2D-array of float; size (3=x/y/z, nb_points)
    :param w_hat: w basis vector of (u,v,w)basis problem
    :type u_hat: numpy 2D-array of float; size (3=x/y/z, nb_points)
    :return: p points,associate sinus,associate cosinus
    :rtype:  numpy 2D-array of float;size (3=x/y/z, nb_points)
             numpy 1D-array of float;size (nb_points)
             numpy 1D-array of float;size (nb_points)
    """

    sin_mu = numpy.sin(mu)
    cos_mu = numpy.sqrt(1 - sin_mu ** 2)  # with our conventions, mu is small and cos(mu)>0
    p_mu = s + numpy.einsum('i,ij->ij', r_eq, numpy.einsum('i,ij->ij', cos_mu, w_hat) + numpy.einsum('i,ij->ij', sin_mu, u_hat)) \
           + numpy.einsum('i,ij->ij', delta, v_hat)

    return p_mu, cos_mu, sin_mu


def p_of_mu(mu, s, r_eq, delta, u_hat, v_hat, w_hat):
    """
    Compute the p point,the sinus associated and the cosinus associated for analytic problem
    Return the p point,the sinus associated and the cosinus associated

    :param mu: mu angle of analytic problem
    :type mu: float
    :param s: position of the sensor in cartesian coordinates
    :type s: numpy 1D-array of float; size (3=x/y/z)
    :param r_eq: range of problem
    :type r_eq: float
    :param delta: correction factor of v in (u,v,w)basis problem
    :type delta: float
    :param u_hat: u basis vector of (u,v,w)basis problem
    :type u_hat: numpy 1D-array of float; size (3=x/y/z)
    :param v_hat: v basis vector of (u,v,w)basis problem
    :type u_hat: numpy 1D-array of float; size (3=x/y/z)
    :param w_hat: w basis vector of (u,v,w)basis problem
    :type u_hat: numpy 1D-array of float; size (3=x/y/z)
    :return: p point,associate sinus,associate cosinus
    :rtype:  numpy 1D-array of float; size (3=x/y/z)
             float
             float
    """
    sin_mu = numpy.sin(mu)
    cos_mu = numpy.sqrt(1 - sin_mu ** 2)  # with our conventions, mu is small and cos(mu)>0
    p_mu = s + r_eq * (cos_mu * w_hat + sin_mu * u_hat) + delta * v_hat

    return p_mu


def h_of_mu(mu, s, r_eq, delta, u_hat, v_hat, w_hat):
    """
    Compute the p point, and the height associated
    Return the height of computed point and this point

    :param mu: mu angle of analytic problem
    :type mu: float
    :param s: position of the sensor in cartesian coordinates
    :type s: numpy 1D-array of float; size (3=x/y/z)
    :param r_eq: range of problem
    :type r_eq: float
    :param delta: correction factor of v in (u,v,w)basis problem
    :type delta: float
    :param u_hat: u basis vector of (u,v,w)basis problem
    :type u_hat: numpy 1D-array of float; size (3=x/y/z)
    :param v_hat: v basis vector of (u,v,w)basis problem
    :type u_hat: numpy 1D-array of float; size (3=x/y/z)
    :param w_hat: w basis vector of (u,v,w)basis problem
    :type u_hat: numpy 1D-array of float; size (3=x/y/z)
    :return: height of p point, p point
    :rtype:  float; numpy 1D-array of float; size (3=x/y/z)

    """
    p_mu = p_of_mu(mu, s, r_eq, delta, u_hat, v_hat, w_hat)
    h_mu = height_fast(p_mu)
    return h_mu, p_mu


def herror2ofx(x, dmu_max, mu_0, s, w_hat, u_hat, v_hat, r_eq, delta, h_target):
    """
     cost function for rescale mu
     return the cost

    :param x: infix variable for cost function
    :type x: float
    :param dmu_max: mu derivative
    :type dmu_max: float
    :param mu_0: initial guess for mu angle
    :type mu_0: float
    :param s: position of the sensor in cartesian coordinates
    :type s: numpy 1D-array of float; size (3=x/y/z)
    :param r_eq: range of problem
    :type r_eq: float
    :param delta: correction factor of v in (u,v,w)basis problem
    :type delta: float
    :param u_hat: u basis vector of (u,v,w)basis problem
    :type u_hat: numpy 1D-array of float; size (3=x/y/z)
    :param v_hat: v basis vector of (u,v,w)basis problem
    :type u_hat: numpy 1D-array of float; size (3=x/y/z)
    :param w_hat: w basis vector of (u,v,w)basis problem
    :type u_hat: numpy 1D-array of float; size (3=x/y/z)
    :param h_target: height targeted to re-geolocate pixel without noise
    :type h_target: float

    :return: the cost
    :rtype: float
    """
    mu = mu_0 + x * dmu_max  # moving by 1 in x moves by dmu_max in mu so *approximately* by dh_max in h
    h_error_square_mu = (h_of_mu(mu, s, r_eq, delta, u_hat, v_hat, w_hat)[0] - h_target) ** 2
    return h_error_square_mu


def pointcloud_height_geoloc_vect(p_noisy, h_noisy, s, vs, range_target, h_target,
                                  recompute_doppler=False,
                                  recompute_range=False,
                                  verbose=False,
                                  max_iter_grad=1, height_goal=1.e-3):
    """
    Compute the height of noisy point

    :param p_noisy: position of the noisy pixel in cartesian coordinates
    :type p_noisy: 2D-array of float; size (3=x/y/z, nb_points)
    :param h_noisy: height of the noisy pixel in geographic coordinates
    :type h_noisy: 1D-array of float
    :param s: position of the sensor in cartesian coordinates
    :type s: 2D-array of float; size (3=x/y/z, nb_points)
    :param vs: velocity of the sensor in cartesian coordinates
    :type vs: 2D-array of float; size (3=x/y/z, nb_points)
    :param range_target: Range targeted (can be overwritted by recomputing true range of the pixel using recomputeR = True)
    :type range_target: 1D-array of float
    :param h_target: height targeted to re-geolocate pixel without noise
    :type h_target: 1D-array of float
    :param recompute_doppler: option to recomputed Doppler (Zdop) instead of using nul vector
    :type recompute_doppler: boolean
    :param max_iter_grad: number of iteration for gradient method
    :type max_iter_grad: integer
    :param recompute_range: option to recomputed range instead of using Rtarget
    :type recompute_range: boolean
    :param verbose: display progress of optimize function
    :type verbose: boolean
    :param max_iter_grad: number of iteration
    :type max_iter_grad: integer
    :param height_goal: precision criteria for using numerical method
    :type height_goal: float

    :return:  position of pixel in geocentric coordinates, position of pixel in geographic coordinates,pixel height
    :rtype: numpy 2D-array of float; size (3=x/y/z, nb_points)
            numpy 2D-array of float; size (3=x/y/z, nb_points)
            numpy 1D-array of float; size (3=x/y/z)
    """
    # 0 - Initiate logging service
    logger = logging.getLogger("pointcloud_height_geoloc_vect()")
    logger.info("PROCESSING")

    # 1 define variable of problem

    # 1.1 - normalize the position of the sensor in cartesian coordinates
    s_hat = normalize_vect(s)

    # 1.2 - define an adapted orthonormal basis
    v_hat = normalize_vect(vs)  # points along track
    u_hat = normalize_vect(
        numpy.cross(s_hat, v_hat))  # s and v are almost ortho; no danger of colinear (points roughly cross track)
    w_hat = numpy.cross(u_hat, v_hat)  # points roughly to nadir

    # 1.3 - define sensor vector position
    r_s = numpy.linalg.norm(s, axis=1)

    # 1.4 - compute delta
    if recompute_doppler:
        delta = numpy.einsum('ij,ij->i', p_noisy - s, v_hat)
    else:
        delta = numpy.zeros(len(v_hat))

    # 1.5 - first computation of range
    if recompute_range:
        r_eq = numpy.linalg.norm(p_noisy - s, axis=1)
    else:
        r_eq = range_target

    # 1.6 - correction of range
    r_eq = numpy.sqrt(r_eq ** 2 - delta ** 2)

    # 1.7 - initial guess for mu uses the input point as a first start
    mu_0 = numpy.arctan2(numpy.einsum('ij,ij->i', p_noisy - s, u_hat), numpy.einsum('ij,ij->i', p_noisy - s, w_hat))

    # 1.8 - tweak the intial guess: move by dh in the w_hat direction, taking the spherical earth correction factor into account
    # 1.81 - crude approximation  (acceptable but we are computing an order 1 correction)
    R_sph_earth = (GEN_RAD_EARTH_EQ + GEN_RAD_EARTH_POLE) / 2.
    # 1.82 - altitude approximation for earth spherical model
    Alt_sph = r_s - R_sph_earth
    #  1.83 - correction factor
    fact_sph = R_sph_earth / (R_sph_earth + Alt_sph)
    # 1.84 - mu correction for earth spherical model
    mu_0 += fact_sph * (h_target - h_noisy) / r_eq / numpy.sin(mu_0)

    # findroot h**2(mu)==h_target using the gradient method (and analytical mu derivatives)
    # note that this can converge poorly if we aim too close to h=0
    iter_grad = 0
    mu = numpy.copy(mu_0)
    while iter_grad <= max_iter_grad:
        # 2 - resolution of G_tilde(mu) = A + d_mu*B = h_mu**2

        # redefinition of problem variable due to mu redefinition for each iteration
        p_mu, cos_mu, sin_mu = p_of_mu_vect(mu, s, r_eq, delta, u_hat, v_hat, w_hat)
        # 2.1 - resolution of G_tilde(mu) = A + d_mu*B = h_mu**2 considering B is equal to 0
        # this is a bit of code duplication wrt height_fast but we will reuse most of the internal variables later...
        x_mu, y_mu, z_mu = p_mu[:, 0], p_mu[:, 1], p_mu[:, 2]
        d_mu = numpy.sqrt(x_mu ** 2 + y_mu ** 2)
        # compute the homothety factor for an ellipsoid going through p_mu
        alpha = numpy.sqrt(d_mu ** 2 / GEN_RAD_EARTH_EQ ** 2 + z_mu ** 2 / GEN_RAD_EARTH_POLE ** 2)
        # compute cosinus theta
        c_th = z_mu / (alpha * GEN_RAD_EARTH_POLE)
        # compute sinus theta
        s_th = d_mu / (alpha * GEN_RAD_EARTH_EQ)
        terme_xy = GEN_RAD_EARTH_EQ * s_th - d_mu
        terme_z = GEN_RAD_EARTH_POLE * c_th - z_mu

        # 2 -computation of G
        G = terme_xy ** 2 + terme_z ** 2
        # derivative of G
        dG_dth = 2 * (GEN_RAD_EARTH_EQ * c_th * terme_xy - GEN_RAD_EARTH_POLE * s_th * terme_z)
        # second derivative of G
        d2G_dth2 = 2 * ((GEN_RAD_EARTH_EQ * c_th) ** 2 + (
            GEN_RAD_EARTH_POLE * s_th) ** 2 - GEN_RAD_EARTH_EQ * s_th * terme_xy - GEN_RAD_EARTH_POLE * c_th * terme_z)  # this is ~GEN_RAD_EARTH_EQ**2>>1, so never 0

        # 3 - computation of A
        A = G - dG_dth ** 2 / d2G_dth2 / 2

        # 4 - computation of h_mu
        h_mu = numpy.sign(alpha - 1.) * numpy.sqrt(A)

        h_error_mu = numpy.abs(h_mu - h_target)
        if iter_grad == 0:
            h_error_mu_0 = h_error_mu  # save for later

        # 5 - compute B for G_tilde(mu) = A + d_mu*B = h_mu**2
        dp_mu_dmu = numpy.einsum('i,ij->ij', r_eq,
                              numpy.einsum('i,ij->ij', -sin_mu, w_hat) + numpy.einsum('i,ij->ij', cos_mu, u_hat))
        dz_dmu = dp_mu_dmu[:, 2]
        dd_dmu = (dp_mu_dmu[:, 0] * p_mu[:, 0] + dp_mu_dmu[:, 1] * p_mu[:, 1]) / d_mu
        dG_dmu = -2 * (dd_dmu * terme_xy + dz_dmu * terme_z)
        d2G_dthdmu = 2 * (-dd_dmu * GEN_RAD_EARTH_EQ * c_th + dz_dmu * GEN_RAD_EARTH_POLE * s_th)
        d3G_dth2dmu = 2 * (dd_dmu * GEN_RAD_EARTH_EQ * s_th + dz_dmu * GEN_RAD_EARTH_POLE * c_th)
        B = dG_dmu - d2G_dthdmu * dG_dth / d2G_dth2 + dG_dth ** 2 * d3G_dth2dmu / d2G_dth2 ** 2 / 2

        # 7 - compute delta_mu
        delta_mu = (h_target ** 2 - A) / B

        # 6 - rescale mu before iteration
        mu += delta_mu  # note that if this is the last iteration and we got that far (error criterion not met), then this is never applied
        iter_grad += 1

    nfev_minimize_scalar = 0  # to store the number of function evaluations during the numerical optimization
    indices = numpy.where(h_error_mu > height_goal)
    for i in indices[0]:
        dh_max = 2 * h_error_mu_0[i]
        dmu_max = numpy.abs(fact_sph[i] * dh_max / r_eq[i] / numpy.sin(mu_0[i]))

        if verbose:
            i_verbose = 3
        else:
            i_verbose = 0

        sol = opt.minimize_scalar(
            lambda x: herror2ofx(x, dmu_max, mu_0[i], s[i], w_hat[i], u_hat[i], v_hat[i], r_eq[i], delta[i],
                                 h_target[i]), bounds=(-1., 1.), method='bounded',
            options={'xatol': height_goal / dh_max, 'disp': i_verbose})
        mu_sol = mu_0[i] + dmu_max * sol.x
        h_mu[i], p_mu[i] = h_of_mu(mu_sol, s[i], r_eq[i], delta[i], u_hat[i], v_hat[i],
                                   w_hat[i])  # a bit of a waste of computation here for the sake of readability

        # we could just recompute p_mu and extract h from sol (but we then need to be careful about the sign etc...)
        nfev_minimize_scalar = sol.nfev

    # projection of geocentric coordinates into geographic coordinates
    p_final_llh = project_array(p_mu)
    logger.info("ENDING")

    return p_mu, p_final_llh, h_mu, (iter_grad, nfev_minimize_scalar)



def convert_llh2ecef(lat, lon, height, rad_e=GEN_RAD_EARTH_EQ, rad_p=GEN_RAD_EARTH_POLE):
    """
    Converting from geodetic coordinates (lat, lon, height) to cartesian coordinates (x, y, z)

    :arg lat: latitude of the record
    :type lat: float
    :arg lon: longitude of the record
    :type lon: float
    :arg height: height of the record
    :type height: float
    :arg rad_e: radius of the Earth model (WGS84 ellipsoid) at the equator
    :type rad_e: float
    :arg rad_p: radius of the Earth model at the pole
    :type rad_p: float
    :returns: tuple. Contain : base axes (x, y and z)
    :rtype: float,float,float

    """
    lat = lat * numpy.pi / 180.
    lon = lon * numpy.pi / 180.

    e = numpy.sqrt(rad_e * rad_e - rad_p * rad_p) / rad_e
    Rn = rad_e / numpy.sqrt(1 - e * e * numpy.sin(lat) * numpy.sin(lat))

    axe_x = (Rn + height) * numpy.cos(lat) * numpy.cos(lon)
    axe_y = (Rn + height) * numpy.cos(lat) * numpy.sin(lon)
    axe_z = (Rn * (1. - e * e) + height) * numpy.sin(lat)

    return axe_x, axe_y, axe_z



