"""
.. module swot_hr.common.rdf.rdf_exception.py
    :synopsis: Exceptions for JPL RDF Files
    Created on 24 avr. 2013

.. moduleauthor: Capgemini

    $Id: rdf_exception.py 1465 2016-07-01 10:05:12Z nestival $
This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""

from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

from ressources.utils.common_exception import CommonException

class RdfException(CommonException):
    """
    Exceptions for JPL RDF Files
    """


    def __init__(self, message):
        """
        Constructor
        
        Args;
            message(str) : Exception message
        """
        CommonException.__init__(self, message)
