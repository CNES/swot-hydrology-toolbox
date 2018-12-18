"""
.. module swot_hr.common.utils.vincenty_direct_formula.py
    :synopsis: Vincenty Direct and Inverse Solution of Geodesics on the Ellipsoid (c) Chris Veness 2005-2012              
                                                                                               
 from: Vincenty direct formula - T Vincenty, "Direct and Inverse Solutions of Geodesics on the  
       Ellipsoid with application of nested equations", Survey Review, vol XXII no 176, 1975    
       http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf 
    Created on March 2014.

.. moduleauthor: Capgemini

    $Id: vincenty_direct_formula.py 1465 2016-07-01 10:05:12Z nestival $
This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""

from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import math
import numpy

from ressources.utils.common_exception import CommonException

GEN_RAD_EARTH_EQ = 6378137
GEN_RAD_EARTH_POLE = 6356752.314245179497563967
ERROR_GEN_DIST_VINCENTY = "Error during distance vincenty process : "

def dest_vincenty(lat1, lon1, brng, dist, rad_e = GEN_RAD_EARTH_EQ, rad_p = GEN_RAD_EARTH_POLE):
    """
    Calculates destination point given start point lat/long, bearing & distance, 
    using Vincenty inverse formula for ellipsoids

    Args:
        lat1(double): latitude of the first point in decimal degrees

        lon1(double): longitude of the first point in decimal degrees

        brng(double): initial bearing in decimal degrees

        dist(double): distance along bearing in meters

        rad_e(float): radius of the Earth model (WGS84 ellipsoid) at the equator (optional)

        rad_p(float): radius of the Earth model at the pole (optional)
    
        Returns:
        Lat, Lon. destination point
    """
    f = 1/298.257223563  # WGS-84 ellipsiod
    s = dist
    alpha1 = brng*math.pi/180 
    sinAlpha1 = math.sin(alpha1)
    cosAlpha1 = math.cos(alpha1)
    
    tanU1 = (1-f) * math.tan(lat1*math.pi/180)
    cosU1 = 1 / math.sqrt((1 + tanU1*tanU1))
    sinU1 = tanU1*cosU1
    sigma1 = math.atan2(tanU1, cosAlpha1)
    sinAlpha = cosU1 * sinAlpha1
    cosSqAlpha = 1 - sinAlpha*sinAlpha
    uSq = cosSqAlpha * (rad_e*rad_e - rad_p*rad_p) / (rad_p*rad_p)
    A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
    B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))
    
    sigma = s / (rad_p*A)
    sigmaP = 2*math.pi
    while (abs(sigma-sigmaP) > 1e-12):
        cos2SigmaM = math.cos(2*sigma1 + sigma)
        sinSigma = math.sin(sigma)
        cosSigma = math.cos(sigma)
        deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-\
          B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)))
        sigmaP = sigma
        sigma = s / (rad_p*A) + deltaSigma
    
    
    tmp = sinU1*sinSigma - cosU1*cosSigma*cosAlpha1
    lat2 = math.atan2(sinU1*cosSigma + cosU1*sinSigma*cosAlpha1, (1-f)*math.sqrt(sinAlpha*sinAlpha + tmp*tmp))
    lambd = math.atan2(sinSigma*sinAlpha1, cosU1*cosSigma - sinU1*sinSigma*cosAlpha1)
    C = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha))
    L = lambd - (1-C) * f * sinAlpha * (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)))
    lon2 = (lon1*math.pi/180+L+3*math.pi)%(2*math.pi) - math.pi # normalise to -180...+180
    
    revAz = math.atan2(sinAlpha, -tmp) # final bearing, if required
    finalBearing = revAz * 180/math.pi
    lat2 = lat2 * 180/math.pi
    lon2 = lon2 * 180/math.pi
    return (lat2, lon2, finalBearing)

def dist_vincenty(lat1, lon1, lat2, lon2, rad_e = GEN_RAD_EARTH_EQ, rad_p = GEN_RAD_EARTH_POLE):
    """
    Calculate the distance between two points on the surface of ellipsoids by using Vincenty formula.

    Args:
        lat1(float): latitude of the record A

        lon1(float): longitude of the record A

        lat2(float): latitude of the record B

        lon2(float): longitude of the record B

        rad_e(float): radius of the Earth model (WGS84 ellipsoid) at the equator (optional)

        rad_p(float): radius of the Earth model at the pole (optional)

    Raises:
        CommonException

    Returns:
        float. distance
    """
    if lat1 == 0.:
        lat1 = 0.000000001
    if lat2 == 0.:
        lat2 = 0.000000001

    f = 1/298.257223563  # WGS-84 ellipsoid params
    L = (lon2-lon1) * math.pi/180
    U1 = math.atan((1-f) * math.tan(lat1*math.pi/180))
    U2 = math.atan((1-f) * math.tan(lat2*math.pi/180))
    sinU1 = math.sin(U1)
    cosU1 = math.cos(U1)
    sinU2 = math.sin(U2)
    cosU2 = math.cos(U2)

    if (lat1 == lat2) and (lon1 == lon2):
        return 0.
    else:
        lambd = L
        lambda_p = 1000
        iter_limit = 100
        while numpy.abs(lambd-lambda_p) > 1e-12 and iter_limit > 0:
            iter_limit = iter_limit - 1
            sin_lambda = math.sin(lambd)
            cos_lambda = math.cos(lambd)
            sin_sigma = math.sqrt((cosU2*sin_lambda) * (cosU2*sin_lambda) + (cosU1*sinU2-sinU1*cosU2*cos_lambda) * (cosU1*sinU2-sinU1*cosU2*cos_lambda))

            if sin_sigma == 0:
                exception = CommonException(ERROR_GEN_DIST_VINCENTY + "sin_sigma = 0: co-incident points")
                raise exception

            cos_sigma = sinU1*sinU2 + cosU1*cosU2*cos_lambda
            sigma = math.atan2(sin_sigma, cos_sigma)
            sin_alpha = cosU1 * cosU2 * sin_lambda / sin_sigma
            cos_sq_alpha = 1 - sin_alpha *  sin_alpha
            cos_2_sigma_m = cos_sigma - 2 * sinU1 * sinU2 / cos_sq_alpha

            if numpy.isnan(cos_2_sigma_m):
                cos_2_sigma_m = 0 #  // equatorial line: cos_sq_alpha=0 (6)

            C = f / 16 * cos_sq_alpha * (4 + f*(4 - 3 * cos_sq_alpha))
            lambda_p = lambd
            lambd = L + (1-C) * f * sin_alpha * (sigma + C*sin_sigma*(cos_2_sigma_m+C*cos_sigma*(-1+2*cos_2_sigma_m*cos_2_sigma_m)))


        if iter_limit == 0:
            exception = CommonException(ERROR_GEN_DIST_VINCENTY + "formula failed to converge, lambda diff = " + str(numpy.abs(lambd-lambda_p)))
            raise exception

        uSq = cos_sq_alpha * (rad_e*rad_e - rad_p*rad_p) / (rad_p*rad_p)
        A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
        B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))
        delta_sigma = B*sin_sigma*(cos_2_sigma_m+B/4*(cos_sigma*(-1+2*cos_2_sigma_m*cos_2_sigma_m)- B/6*cos_2_sigma_m*(-3+4*sin_sigma*sin_sigma)*(-3+4*cos_2_sigma_m*cos_2_sigma_m)))
        result = rad_p*A*(sigma-delta_sigma)

        return result
