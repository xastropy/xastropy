"""
#;+ 
#; NAME:
#; igm.tau_eff
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for tau effective
#;   07-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function

import numpy as np
import os

from xastropy.xutils import xdebug as xdb


#    Calculate tau_effective for the Lyman series using the EW
#    approximation (e.g. Zuo 93)
def ew_teff_lyman():
    """ tau effective (follows ew_teff_lyman.pro from XIDL)

    Parameters:
      ilambda: float
        Observed wavelength 
      zem: float 
        Emission redshift of the source [sets which Lyman lines are included]

    JXP 07 Nov 2014
    """
