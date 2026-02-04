#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : SWOT KARIN
#
# ======================================================
# HISTORIQUE
# VERSION:1.0.0:::2020/09/11: initial version
# FIN-HISTORIQUE
# ======================================================
"""
.. module run_rivers_estimation.py
    :synopsis: Metrix to compute relative errors

.. moduleauthor: CNES (L. Labat-Allee)

Copyright (c) 2020 CNES. All rights reserved.
"""
import logging
import numpy as np


class Scores(object):
    """
    Class for computing scores from Confusion Matrix results
    """

    def __init__(self, fp, fn, tp, tn):
        """
        Constructor
        """
        # FP False Positive
        self.fp = fp
        # FN False Negative
        self.fn = fn
        # TP True Positive
        self.tp = tp
        # TN True Negative
        self.tn = tn

    def recall(self):
        """
        Compute recall (TPR)
        """

        # Recall: True Positive Rate
        return float(self.tp) / (float(self.tp) + float(self.fn))

    def false_positive_rate(self):
        """
        Compute False Positive Rate (FPR)
        """

        # Recall: True Positive Rate
        return float(self.fp) / (float(self.fp) + float(self.tn))

    def error_rate(self):
        """
        Compute Error Rate (ER)
        """

        # ER: Error Rate
        if self.tp + self.fn > 0:
            er = (float(self.fp) + float(self.fn)) / (float(self.tp) + float(self.fn))
        else:
            er = -1

        return er

    def mcc_score(self):
        """
        Compute Matthews Correlation Coefficent (MCC)
        """

        # MCC: Matthews Correlation Coefficent
        if self.tp > 0 and self.tn > 0 and self.fp > 0 and self.fn > 0:
            mcc = float(self.tp * self.tn - self.fp * self.fn) / \
                  float(np.sqrt(float((self.tp + self.fp) *
                                      (self.tp + self.fn) * (self.tn + self.fp) * (self.tn + self.fn))))
        else:
            mcc = -1

        return mcc

    def precision(self):
        """
        Compute Precision
        """

        return self.tp / (self.tp + self.fp)

    def accuracy(self):
        """
        Compute Accuracy
        """

        return (self.tp + self.tn) / (self.tp + self.fn + self.tn + self.fp)

    def f1_score(self):
        """
        Compute F1-score
        """

        precision = self.precision()
        recall = self.recall()
        return 2 * (precision * recall) / (precision + recall)

    def nb_pixels(self):
        """
        Compute total number of pixels
        :return: total number of pixels
        """

        return (self.tp + self.fp + self.tn + self.fn)

#----------------------------------------

def compute_surface_relative_error(s_ref_area, s_compared_area):
    """
    For one reach (one water body), compute the relative error in surface area

    Args:
            s_ref_area (float): reference value (ideal or truth)
            s_compared_area (float): measured value
    Returns:
            (numpy.float64) : surface area relative error
    """

    relative_err = abs(s_compared_area - s_ref_area) / s_ref_area

    return relative_err

#----------------------------------------

class Metrix(object):
    """
    Class for initializing a metrix structure
    """
    def __init__(self, name, **kwargs):

        # Metrix name: 'reference_water_area'/'reference_wse' or 'confusion_matrix'
        self._name = name
        self.__dict__.update(**kwargs)

    def tojson(self):
        return self.__dict__


    def compute_relative_errors(self):
        """
        Compute relative error from metrix:
        - if name is 'reference_water_area': from a list of water surface (compared, reference)
        - if name is 'confusion_matrix': from confusion matrix metrix (tp, fp, fn)

        Returns:
            (numpy.float64) : tab of relative errors
        """
        logger = logging.getLogger(self.__class__.__name__)

        if self._name in ['reference_surface_area']:

            # Get compared and reference values
            try:
                reference = np.asarray(getattr(self, 'reference'))   # reference: Ideal or Truth
                compared = np.asarray(getattr(self, 'compared'))     # compared
                assert len(compared) == len(reference)
                
            except Exception:
                logger.error("Error when getting metrix for {0} case.".format(self._name))
                raise

            # Relative error
            #er = abs(reference - compared) / reference
            er = (compared - reference) / reference

        elif self._name == 'confusion_matrix':

            # Get confusion matrix metrix
            try:
                # True positive
                tp = np.asarray(getattr(self, 'tp'))
                # False positive
                fp = np.asarray(getattr(self, 'fp'))
                # False negative
                fn = np.asarray(getattr(self, 'fn'))
            except:
                logger.error("Error when getting metrix for {0} case".format(self._name))

            # ER: Error Rate
            er = np.zeros_like(tp)
            ind_neg = np.where((tp + fn) <= 0)
            ind_pos = np.where((tp + fn) > 0)
            er[ind_neg] = -1
            er[ind_pos] = (fp[ind_pos] + fn[ind_pos]) / (tp[ind_pos] + fn[ind_pos])
            
        else:
            logger.error("{0} case is not defined!".format(self._name))
            er = -1

        return er


    def compute_absolute_errors(self):
        """
        Compute absolute error from metrix:
        - if name is in ['reference_wse']: from a list of water surface elevations (compared, reference)

        Returns:
            (numpy.float64) : tab of absolute errors
        """
        logger = logging.getLogger(self.__class__.__name__)
        if self._name in ['reference_wse']:

            # Get compared and reference values
            try:
                reference = np.asarray(getattr(self, 'reference'))        # reference: Ideal or Truth
                compared = np.asarray(getattr(self, 'compared'))     # compared
                assert len(compared) == len(reference)
                
            except Exception:
                logger.error("Error when getting metrix for {0} case.".format(self._name))
                raise

            # Absolute error
            #er = abs(reference - compared)
  
            er = compared - reference
            
        else:
            logger.error("{0} case is not defined!".format(self._name))
            er = -1

        return er
