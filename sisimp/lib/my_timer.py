# -*- coding: utf-8 -*-
'''
.. module my_timer.py
    :synopsis: Gestion du temps d'execution d'un code
    Version 2.0 = Modification de la fonction "info" (prise en compte d'un parametre pour indique si on donne la duree depuis la date start ou depuis le dernier info)
    Version 1.0 = Created on 03/06/2016

.. module author: Claire POTTIER - CNES DCT/SI/TR

This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


'''
from __future__ import absolute_import, division, print_function, unicode_literals
 
import time
import math


class Timer(object):
    
    def __init__(self):
        '''
        Constructeur du timer
            
        ** Attributes **
        - start_time : debut d'execution du code
        - tmp_time : temps intermediaire
        ''' 
        self.start_time = 0.0
        self.tmp_time = 0.0
    
    #----------------------------------------
    
    def start(self):
        '''
        Initialise le compteur de depart
        '''
        self.start_time = time.time()
    
    #----------------------------------------
        
    def info(self, IN_flagDebut=1):
        '''
        Calcule la duree depuis une date de reference:
        - la derniere info demandee (ou debut d'execution pour la 1e info depuis la date d'execution) si IN_flagDebut=0
        - la date de debut d'execution si IN_flagDebut=1
        La date de l'execution de cette fonction est sauvegardee en tant que date intermediaire pour une utilisation ulterieure.
            
        :param IN_flagDebut: =0 si la date de reference correspond a la derniere info demandee, =1 si elle correspond a la date de debut d'execution
        :return: la duree depuis la date de reference
        :rtype: str
        '''
        
        # 1 - On calcule la duree
        if ( self.tmp_time == 0 ) or ( IN_flagDebut == 1) :
            curTime = time.time() - self.start_time
        else:
            curTime = time.time() - self.tmp_time
        
        # 2 - On met a jour le temps intermediaire
        self.tmp_time = time.time()
            
        # 3 - Retourne la duree au format hh:mm:ss
        return "[Timer-INFO] Step executed in %s" % ( self.printDelay(curTime) )
  
    def stop(self):  
        '''
        Calcule la duree totale d'execution
            
        :return: la duree totale d'execution
        :rtype: str
        '''
        
        # 1 - On calcule la duree
        totalTime = time.time() - self.start_time
        
        # 2 - Retourne la duree au format hh:mm:ss
        return "[Timer-INFO] Total execution time in %s" % ( self.printDelay(totalTime) )
    
    #----------------------------------------
    
    def printDelay(self, IN_time):
        '''
        Convertit la duree en parametre d'entree en chaine de caracteres
            
        :param IN_time: duree en secondes
        :type IN_time: int
            
        :return: la duree correspondant a IN_time au format "HH:MM:SS" ou "DD days HH:MM:SS" si besoin
        :rtype: str
        '''
        
        # Calcul des secondes
        tmp_MM = math.floor(IN_time / 60.0)
        SS = IN_time - 60.0 * tmp_MM
        
        # Calcul des minutes
        HH = math.floor(tmp_MM / 60.0)
        MM = tmp_MM - 60.0 * HH
        
        # Calcul des heures si > 24h
        if HH >= 24.0:
            DD = math.floor(HH / 24.0)
            HH -= 24.0 * DD
            return "%d days %02d:%02d:%02d" % ( DD , HH , MM , SS )
        else:
            return "%02d:%02d:%02d" % ( HH , MM , SS )
            
#######################################

if __name__ == '__main__':
    
    myTime = Timer()
    myTime.start()
    time.sleep(5)
    print(myTime.info())
    time.sleep(5)
    print(myTime.info())
    print(myTime.stop())
