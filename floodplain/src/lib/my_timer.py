# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2019/05/17:version initiale.
# FIN-HISTORIQUE
# ======================================================
"""
.. module:: my_timer.py
    :synopsis: Gestion du temps d'execution d'un code
     Version 2.0 = Modification de la fonction "info" (prise en compte d'un parametre 
     pour indiquer si on donne la duree depuis la date start ou depuis le dernier info)
     Version 1.0 = Created on 2016/03/06

.. moduleauthor:: Claire POTTIER - CNES DCT/SI/TR

..
   This file is part of the SWOT Hydrology Toolbox
   Copyright (C) 2018 Centre National d’Etudes Spatiales
   This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""
from __future__ import absolute_import, division, print_function, unicode_literals

import time
import math


class Timer(object):
    """
        class Timer
    """
    def __init__(self):
        """
        Constructeur du timer


        Variables of the object:
         - start_time : debut d'execution du code
         - tmp_time : temps intermediaire
        """
        self.start_time = 0.0
        self.tmp_time = 0.0

    # ----------------------------------------

    def start(self):
        """
        Initialise le compteur de depart
        """
        self.start_time = time.time()

    # ----------------------------------------

    def info(self, in_flag_debut=1):
        """
        Calcule la duree depuis une date de reference:
         - la derniere info demandee (ou debut d'execution pour la 1e info depuis la date d'execution) si in_flag_debut=0
         - la date de debut d'execution si in_flag_debut=1
         
        La date de l'execution de cette fonction est sauvegardee en tant que date intermediaire pour une utilisation ulterieure.

        :param in_flag_debut: =0 si la date de reference correspond a la derniere info demandee, =1 si elle correspond a la date de debut d'execution
        :type in_flag_debut: int

        :return: la duree depuis la date de reference
        :rtype: str
        """

        # 1 - On calcule la duree
        if (self.tmp_time == 0) or (in_flag_debut == 1):
            cur_time = time.time() - self.start_time
        else:
            cur_time = time.time() - self.tmp_time

        # 2 - On met a jour le temps intermediaire
        self.tmp_time = time.time()

        # 3 - Retourne la duree au format hh:mm:ss
        return "[Timer-INFO] Step executed in %s" % self.print_delay(int(cur_time))

    def stop(self):
        """
        Calcule la duree totale d'execution

        :return: la duree totale d'execution
        :rtype: str
        """

        # 1 - On calcule la duree
        total_time = time.time() - self.start_time

        # 2 - Retourne la duree au format hh:mm:ss
        return "[Timer-INFO] Total execution time in %s" % self.print_delay(int(total_time))

    # ----------------------------------------

    def print_delay(self, in_time):
        """
        Convertit la duree en parametre d'entree en chaine de caracteres

        :param in_time: duree en secondes
        :type in_time: int

        :return: la duree correspondant a in_time au format "HH:MM:SS" ou "DD days HH:MM:SS" si besoin
        :rtype: str
        """

        # Calcul des secondes
        tmp_mm = math.floor(in_time / 60.0)
        ss = in_time - 60.0 * tmp_mm

        # Calcul des minutes
        hh = math.floor(tmp_mm / 60.0)
        mm = tmp_mm - 60.0 * hh

        # Calcul des heures si > 24h
        if hh >= 24.0:
            dd = math.floor(hh / 24.0)
            hh -= 24.0 * dd
            retour = "%d days %02d:%02d:%02d" % (dd, hh, mm, ss)
        else:
            retour = "%02d:%02d:%02d" % (hh, mm, ss)
        return retour

#######################################


if __name__ == '__main__':

    my_time = Timer()
    my_time.start()
    time.sleep(5)
    print(my_time.info())
    time.sleep(5)
    print(my_time.info())
    print(my_time.stop())
