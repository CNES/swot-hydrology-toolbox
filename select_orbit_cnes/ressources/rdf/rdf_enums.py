"""
.. module swot_hr.common.rdf.rdf_enums.py
    :synopsis:
    Created on 24 avr. 2013

.. moduleauthor: Capgemini

    $Id: rdf_enums.py 1465 2016-07-01 10:05:12Z nestival $
This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""

from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

# Error messages
RDF_WRONG_PATH = "The RDF file does not exist"
RDF_WRONG_SYNTAX = "Syntax error, the line could not be parsed"
RDF_WRONG_SECTION = "Section does not exist"
RDF_WRONG_PARAM = "Parameter does not exist in this section"
RDF_NULL_PARAM = "Parameter must be indicated"

# Default section name
RDF_DEFAULT = "DEFAULT"

# Special characters
RDF_EOL = ";"
RDF_UNIT = "("
RDF_EQUAL = "="
RDF_COMMENT = "!"