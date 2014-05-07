"""
Module to handle Astronomical data
"""

import numpy as np
import pdb
from astropy.io import fits
from astropy.table import Table

# Convert a list of classes to a Table
#   Useful for generating a binary FITS table (eventually)
#   Assumed that every element in the list is an identical instance
def lclass_to_table(arr):

    from astropy.table import Column
    # Error check
    #if type(self) != instance:

    # Find the attributes
    attr = (arr[0]).__dict__.keys()
    if len(attr) == 0:
        return -1

    # Loop on attributes
    narr = len(arr)
    for iattr in attr:
        # Type
        
        col = Column(name=iattr, data=np.zeros(arr))
        # Loop on the list
        #for iclass in arr:
    
