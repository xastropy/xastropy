"""
#;+ 
#; NAME:
#; stats.basic
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for basic stat calculations
#;   04-Dec-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
from scipy.interpolate import interp1d

from xastropy.xutils import xdebug as xdb

# def perc
# def lin_to_log

def lin_to_log(x, sig):
    """ Convert linear value+error to log 

    Parameters:
      x: float
      sig: float 

    Returns:
      logx, sig_logx
        Log value and error in log

    JXP 26 Mar 2015
    """
    logx = np.log10( x ) 
    lgvar = ((1. / (np.log(10.0)*x))**2) * sig**2
    sig_logx = np.sqrt(lgvar)

    return logx, sig_logx

def perc(x, per=0.68):
    """ Calculate the percentile bounds of a distribution

    Parameters:
      x: float
        numpy array of values
      per: float (0.68)
          Percentile for the calulation

    Returns:
      xper: array
        Value at lower, value at upper

    JXP 04 Dec 2014
    """
    #
    npt = len(x)

    # Sort
    xsort = np.sort(x)
    perx = (np.arange(npt)+1) / npt

    f = interp1d(perx,xsort)

    frac = (1.-per) / 2.

    # Fill
    xper = np.zeros(2)
    try:
        xper[0] = f( frac )
    except ValueError:
        xper[0] = np.min(x)

    try:
        xper[1] = f( 1.-frac )
    except ValueError:
        xper[1] = np.max(x)

    #xdb.set_trace()

    # Return
    return xper
