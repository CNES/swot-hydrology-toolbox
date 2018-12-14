#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# $Id:
#
# ======================================================
#
# Project : SWOT
# Program : common python modules
# Produit par Capgemini.
#
# ======================================================
# HISTORIQUE
#
# FIN-HISTORIQUE
# ======================================================

import unittest
import time
import cnes.common.lib.my_timer as my_timer


class TestMyTimer(unittest.TestCase):
    """
        class for unitary test of my_timer module
    """
    def test_Timer(self):
        """
            unitary test for Timer
        """
        # Init
        temps = my_timer.Timer()
        # start
        temps.start()
        # info
        time.sleep(1)
        info1 = temps.info(IN_flagDebut=1)
        ref1 = "[Timer-INFO] Step executed in 00:00:01"
        print(info1)
        self.assertEqual(info1, ref1)
        time.sleep(2)
        info2 = temps.info(IN_flagDebut=1)
        ref2 = "[Timer-INFO] Step executed in 00:00:03"
        print(info2)
        self.assertEqual(info2, ref2)
        time.sleep(1)
        info3 = temps.info(IN_flagDebut=0)
        ref3 = "[Timer-INFO] Step executed in 00:00:01"
        print(info3)
        self.assertEqual(info3, ref3)
        # change start time
        temps.start_time = temps.start_time - 100000
        # stop
        info4 = temps.stop()
        ref4 = "[Timer-INFO] Total execution time in 1 days 03:46:44"
        print(info4)
        self.assertEqual(info4, ref4)
