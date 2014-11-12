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

from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.utils.misc import isiterable

from xastropy.xutils import xdebug as xdb

# cosm_xz -- Calculates X(z), the absorption path length
# X_Cosmo -- Class that inherits astropy.cosmology FlatLambdaCDM class

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

####
class X_Cosmo(FlatLambdaCDM):
    """A class for extending the astropy Class

    """
    # Initialize with a .dat file
    def __init__(self, H0=70, Om0=0.3):
        # Generate 
        FlatLambdaCDM.__init__(self,H0,Om0)

    #
    def physical_distance(self,z):
        """ Physical line-of-sight distance in Mpc at a given redshift.

        Parameters
        ----------
        z : array_like
          Input redshifts.  Must be 1D or scalar.

        Returns
        -------
        d : ndarray, or float if input scalar
          Physical distance in Mpc to each input redshift.
        """

        from scipy.integrate import quad
        if not isiterable(z):
            #xdb.set_trace()
            return self._hubble_distance * quad(self.phys_inv_efunc, 0, z)[0]

        out = [quad(self.phys_inv_efunc, 0, redshift)[0] for redshift in z]
        return self._hubble_distance * np.array(out)

    def phys_inv_efunc(self, z):
        r""" Function used to calculate :math:`\frac{1}{H_z (1+z)}`.

        Parameters
        ----------
        z : array_like
          Input redshifts.

        Returns
        -------
        E/(1+z) : ndarray, or float if input scalar
          The inverse redshift scaling of the Hubble constant.

        Notes
        -----
        The return value, E/(1+z), is defined such that :math:`H_z = H_0 /
        E`.
        """

        if isiterable(z):
            z = np.asarray(z)
        Om0, Ode0, Ok0 = self._Om0, self._Ode0, self._Ok0
        if self._massivenu:
            Or = self._Ogamma0 * (1 + self.nu_relative_density(z))
        else:
            Or = self._Ogamma0 + self._Onu0
        zp1 = 1.0 + z

        return 1.0 / np.sqrt(zp1 ** 2 * ((Or * zp1 + Om0) * zp1 + Ok0) + Ode0) / zp1






## #################################    
## #################################    
## TESTING
## #################################    
if __name__ == '__main__':

    # dXdz
    z=3
    print('dX/dz at z=%g is %g' % (z, cosm_xz(z, flg=1)) )

    # Physical distance
    cosmo = X_Cosmo(H0=75.) # Vanilla
    dist = cosmo.physical_distance(z)
    print('dist to z=%g is %g Mpc' % (z, dist.value))
