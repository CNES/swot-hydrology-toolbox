"""
.. module inversion_algo.py
    :synopsis: Algorithm classes for the SAR equation
    Created on 14 sept. 2015

.. moduleauthor: Capgemini

    $Id: inversion_algo.py 1520 2017-09-27 14:36:11Z chruiz $
This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""

from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import numpy as np


class inversionCore(object):
    """
    This class permit to resolving geolocation equations.

    This class contain only static methods
    """


    @staticmethod
    def newton_raphson(axe_x, axe_y, axe_z, sensor_x, sensor_y, sensor_z, vsx, vsy, vsz, rad_e, rad_p, height_mean, range_rec, accuracy, fdop_cen):
        """
        Resolving a system of 3 non linear equations: 
            
                - SAR Range equation
                - SAR Doppler equation
                - Earth model equation (WGS-84 ellipsoide) 

        Args:
            axe_x(float): base axe 'x' of the record

            axe_y(float): base axe 'y' of the record

            axe_z(float): base axe 'z' of the record

            sensor_x(float): 'x' sensor coordinates in the Earth centered Cartesian frame at the time of the line acquisition

            sensor_y(float): 'y' sensor coordinates in the Earth centered Cartesian frame at the time of the line acquisition

            sensor_z(float): 'z' sensor coordinates in the Earth centered Cartesian frame at the time of the line acquisition

            vsx(float): 'x' velocity components expressed in the Earth centered Cartesian frame at the time of the line acquisition

            vsy(float): 'y' velocity components expressed in the Earth centered Cartesian frame at the time of the line acquisition

            vsz(float): 'z' velocity components expressed in the Earth centered Cartesian frame at the time of the line acquisition

            rad_e(float): radius of the Earth model (WGS84 ellipsoid) at the equator

            rad_p(float): radius of the Earth model at the pole

            height_mean(float): mean height of the water

            range_rec(float): range between the sensor and the pixel in concern
            
            fdop_cen(float) : Doppler centroid frequency (Hz) for the pixel

            accuracy(float): accuracy for the calculation

        Returns:
            tuple. Contain : modified base axes (x, y and z)
        """
        # equ 2,3 et 4
        val_accuracy = 1.
        
        landa = 3e8/35e9

        iteration = 0
        converged = True
        
#        epsilon = 1.e-1
#        axe_x_org = axe_x + epsilon
#        axe_y_org = axe_y + epsilon
#        axe_z_org = axe_z + epsilon
        
        
        
        while (val_accuracy > accuracy):
            a10 = (sensor_x - axe_x)*(sensor_x - axe_x) + (sensor_y - axe_y)*(sensor_y - axe_y) + (sensor_z - axe_z)*(sensor_z - axe_z) - range_rec*range_rec
            a11 = 2*axe_x - 2*sensor_x
            a12 = 2*axe_y - 2*sensor_y
            a13 = 2*axe_z - 2*sensor_z

            a20 = -2*vsx*(sensor_x - axe_x)/(landa*range_rec) - 2*vsy*(sensor_y - axe_y)/(landa*range_rec) - 2*vsz*(sensor_z - axe_z)/(landa*range_rec) - fdop_cen
            a21 = 2*vsx/(range_rec*landa)
            a22 = 2*vsy/(range_rec*landa)
            a23 = 2*vsz/(range_rec*landa)

            a30 = (axe_x*axe_x)/((rad_e+height_mean)*(rad_e+height_mean)) + (axe_y*axe_y)/((rad_e+height_mean)*(rad_e+height_mean)) + (axe_z*axe_z)/(rad_p*rad_p) -1
            a31 = 2*axe_x/((rad_e+height_mean)*(rad_e+height_mean))
            a32 = 2*axe_y/((rad_e+height_mean)*(rad_e+height_mean))
            a33 = 2*axe_z/(rad_p*rad_p)

            A = np.array([[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]])
            b = np.array([-a10, -a20, -a30])

            x = np.linalg.solve(A, b)
            axe_x = axe_x + x[0]
            axe_y = axe_y + x[1]
            axe_z = axe_z + x[2]

            val_accuracy = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
	          
            iteration = iteration + 1
            if (iteration > 20):
                #print("ITERATION MAX REACHED!!") 
                converged = False
#                axe_x = axe_x_org
#                axe_y = axe_y_org
#                axe_z = axe_z_org
                break
             		   

        return axe_x, axe_y, axe_z, converged


    @staticmethod
    def convert_llh2ecef(lat, lon, height, rad_e, rad_p):
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


    @staticmethod
    def convert_ecef2llh(axe_x, axe_y, axe_z, rad_e, rad_p):
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
