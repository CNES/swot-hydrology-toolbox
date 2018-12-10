"""
Package used to read RDF files used in the JPL simulator.
These files look like ConfigParser files but they are slightly
differerent.
This is the reason why a specific reader has been developped.
The authorized syntax is as follows:

    SECTION NAME
    param name        (unit) = value;    ! comment
    
The reader present here allows to get the section names, parameter
names and parameter values in a given file.
"""