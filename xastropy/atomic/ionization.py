"""
#;+ 
#; NAME:
#; ionization
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for ionization of atoms
#;   03-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function

import numpy as np
import os, pdb
from astropy.io import fits, ascii
import xastropy
#from astropy import constants as const

def photo_cross(Z, ion, E, datfil=None):

    """ Estimate photo-ionization cross-section using Fit parameters
    from Verner et al. 1996, ApJ, 465, 487

    Parameters
    ----------
    Z: Atomic number
    ion : Ionization state (1=Neutral)
    E : Energy to calcualte at [eV]

    Returns
    -------
    sigma : Cross-section (cm^2)
    """

    # Read data
    if datfil == None:
        datfil = xastropy.__path__[0]+'/data/atomic/verner96_photoion_table1.dat'
    dat = ascii.read(datfil)

    # Match
    #pdb.set_trace()
    mt = np.where((Z == dat['Z']) & (ion == dat['N']))[0]
    nmt = len(mt)
    if nmt == 0:
        raise ValueError('photo_cross: %d,%d pair not in our table' % (Z,ion))
    idx = mt[0]
    #
    x = E/dat['E0'][idx] - dat['y0'][idx]
    y = np.sqrt(x**2 + dat['y1'][idx]**2)

    F = (((x-1.)**2 + dat['yw'][idx]**2) * y**(0.5*dat['P'][idx] - 5.5) * 
            (1 + np.sqrt(y/dat['ya'][idx]) )**(-1.*dat['P'][idx]))

    sigma = dat['s0'][idx] * F * 1e-18 # cm^2

    # Energy threshold
    low = np.where(E < dat['Eth'][idx])[0]
    if len(low) > 0: sigma[low] = 0.

    return sigma

# Testing
if __name__ == '__main__':
    # Hydrogen
    print(photo_cross(1,1,13.6))
    print(photo_cross(1,1,13.6*np.arange(1,11)))
