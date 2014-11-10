"""
#;+ 
#; NAME:
#; igm.igm_utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for miscellanious IGM items
#;   07-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function

import numpy as np
import os

from astropy.io import fits

from xastropy.xutils import xdebug as xdb

# cosm_xz -- Calculates X(z), the absorption path length

# cosm_xz -- Calculates X(z), the absorption path length of dxdz
def cosm_xz(z, cosmo=None, zmin=0., flg=0):
    """ Calculates X(z) -- absorption path length or dXdz

    Parameters:
      z: float or ndarray
        Redshift to evaluate at.  May be an array
      cosmo: astropy.cosmology
        Cosmological model to adopt 
      zmin: float (0.)
        Minimum redshift to evaluate at 
      flg: int (0)
        Flag controlling the output
          0 = X(z)
          1 = dX/dz at z

    Returns:
      Xz or dXdz: 
        Absorption path or its derivative

    ToDo:

    JXP 08 Nov 2014
    """

    # Cosmology
    if cosmo == None:
        from astropy.cosmology import Planck13
        cosmo = Planck13

    # Flat?
    if cosmo.Ok(0.) == 0:
        if flg == 0:  # X(z)
            t1 = np.sqrt(cosmo.Om0*(1+z)**3 + cosmo.Ode0)  
            t2 = np.sqrt(cosmo.Om0*(1+zmin)**3 + cosmo.Ode0)  
            # X(z)
            rslt = 2 * (t1-t2) / (3. * cosmo.Om0)
        elif flg == 1:  # dX/dz
            rslt = (cosmo.H0 / cosmo.H(z)) * (1+z)**2
        else: raise ValueError('igm_utils.cosm_xz: Bad flg %d' % flg)
    else:
        raise ValueError('igm_utils.cosm_xz: Not prepared for non-flat cosmology')

    #
    return rslt
