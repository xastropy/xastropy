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

from astropy.io import fits

from xastropy.xutils import xdebug as xdb


#    Calculate tau_effective for the Lyman series using the EW
#    approximation (e.g. Zuo 93)
def ew_teff_lyman(ilambda, zem, fN_model, NHI_MIN=11.5, NHI_MAX=22.0, N_eval=5000,
                  EW_SPLINE=None, bval=24.):
    """ tau effective (follows ew_teff_lyman.pro from XIDL)

    Parameters:
      ilambda: float
        Observed wavelength 
      zem: float 
        Emission redshift of the source [sets which Lyman lines are included]
      bva: float
         -- Characteristics Doppler parameter for the Lya forest
         -- [Options: 24, 35 km/s]
      NHI_MIN: float
         -- Minimum log HI column for integration [default = 11.5]
      NHI_MAX: float
         -- Maximum log HI column for integration [default = 22.0]

    Returns:
      teff: 
        Total effective opacity of all lines contributing

    JXP 07 Nov 2014
    """
    # Lambda
    if not isinstance(ilambda,float):
        raise ValueError('igm.tau_eff: ilambda must be a flaat for now')
    Lambda = ilambda[0]

    # Read in EW spline (if needed)
    if EW_SPLINE == None:
        if int(bval) == 24: EW_FIL = os.environ.get('XIDL_DIR')+'/IGM/EW_SPLINE_b24.fits'
        elif int(bval) == 35: EW_FIL = os.environ.get('XIDL_DIR')+'/IGM/EW_SPLINE_b35.fits'
        else: 
            raise ValueError('igm.tau_eff: Not ready for this bvalue %g' % bval)
        hdu = fits.open(EW_FIL)
        EW_SPLINE = hdu[1].data

    # Assumed Line list
    wrest = ([1215.6701d,  1025.7223,       972.53680,       949.74310,       937.80350, $
            930.74830,  926.22570,       923.15040,       920.96310,       919.35140, $
            918.12940,  917.18060,       916.42900,       915.82400,       915.32900, $
            914.91900,  914.57600,       914.28600,       914.03900,       913.82600, $
            913.64100,  913.48000,       913.33900,       913.21500,       913.10400, $
            913.00600,  912.91800,       912.83900,       912.76800,       912.70300, $
            912.64500] )


    # Find the lines
    gd_Lyman = np.where(Lambda/(1+zem) LT wrest)[0]
    nlyman = len(gd_Lyman)
    if nlyman == 0: return, -1

    # N_HI grid
    lgNval = NHI_MIN + (NHI_MAX-NHI_MIN)*np.arange(N_eval)/(N_eval-1) # Base 10 
    dlgN = lgNval[1]-lgNval[0]
    Nval = 10.**lgNval
    teff_lyman = np.zeros(nlyman)
