"""
#;+ 
#; NAME:
#; arrays
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for array utilities
#;   10-Dec-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

from astropy import units as u

from xastropy.xutils import xdebug as xdb

#
def lst_to_array(lst,mask=None):
    ''' Convert a list into an array
    -- includes quantities
    '''
    # Set mask
    if mask is None:
        mask = (lst == lst)
    
    if isinstance(lst[0],u.quantity.Quantity):
        #Mask and create an array
        unit = lst[0].unit
        newlst= []
        for ii in range(len(mask)):
            if mask[ii]:
                newlst.append(lst[ii].value)
        return np.array(newlst) * unit
    else:
        return np.array(lst)[mask]
