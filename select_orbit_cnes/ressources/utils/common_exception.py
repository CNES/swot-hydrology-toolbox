"""
.. module swot_hr.common.utils.common_exception.py
    :synopsis: Base exception for the SWOT HR programs
    Created on 12 juin 2013

.. moduleauthor: Capgemini

    $Id: common_exception.py 1465 2016-07-01 10:05:12Z nestival $
This file is part of the SWOT Hydrology Toolbox
 Copyright (C) 2018 Centre National d’Etudes Spatiales
 This software is released under open source license LGPL v.3 and is distributed WITHOUT ANY WARRANTY, read LICENSE.txt for further details.


"""

from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

class CommonException(BaseException):
    """
    Generic exception to describe errors that could be raised during the 
    SWOT HR elements' execution.
    
    This class provides two added functions that make abstratction of error
    message formatting:
        
        * *add_message*, just append a message to the already existing message
        * *add_data*, append the message and the data to the already existing message
    """

    # Message seperator
    MESSAGE_SEPERATOR = " :: "

    def __init__(self, message):
        """
        Constructor
        
        Args;
            message(str) : Exception message
        """
        BaseException.__init__(self)
        self.message = message
        
    def add_message(self, message):
        """
        Adds a message into the exception
        
        Args;
            message(str) : Exception message
        """
        self.message += CommonException.MESSAGE_SEPERATOR + message
            
    
    def add_data(self, message, data):
        """
        Adds information and data to the previously registered message
        
        Args;
            message(str) : Exception message
            data(str) : The data value associated to the message
        """
        self.message += CommonException.MESSAGE_SEPERATOR + message +  CommonException.MESSAGE_SEPERATOR +str(data)
        
    
    def __repr__(self, *args, **kwargs):
        return BaseException.__repr__(self, *args, **kwargs) + CommonException.MESSAGE_SEPERATOR + self.message
        
    def __str__(self, *args, **kwargs):
        return BaseException.__str__(self, *args, **kwargs) + CommonException.MESSAGE_SEPERATOR + self.message
        
        