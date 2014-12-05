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

from xastropy.xutils import xdebug as xdb

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

    # Median first
    #medv = np.median(x)

    # Sort
    xsort = np.sort(x)
    perx = (np.arange(npt)+1) / npt

    from scipy.interpolate import interp1d
    f = interp1d(perx,xsort)
    frac = (1.-per) / 2.
    xper = f( [frac, 1.-frac]) 
    #xdb.set_trace()

    # Return
    return xper
