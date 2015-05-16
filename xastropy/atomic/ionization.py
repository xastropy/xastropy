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
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp
from astropy.io import fits, ascii
from astropy import units as u 
#from astropy import constants as const

from xastropy.atomic.elements import ELEMENTS
from xastropy.outils import roman
from xastropy.xutils import xdebug as xdb

from astropy.utils.misc import isiterable

# Path for xastropy
xa_path = imp.find_module('xastropy')[1]

#def ion_name(ion):
#def photo_cross(Z, ion, E, datfil=None, silent=False):

########################## ##########################
########################## ##########################
def ion_name(ion,flg=0,nspace=None):
    """ Convert ion into a string
    JXP on 16 Nov 2014

    Parameters
    ----------
    ion: tuple (Z,ion)
         dict with tags of 'Z' and 'ion'
    flg: int (0)
      0: Roman numeral
      1: Latex with ion notation (e.g C^+)
    nspace: int  (0)
      Number of spaces to insert

    Returns
    -------
    name : string
      e.g. Si II, {\rm Si}^{+}
    """
    if type(ion) is tuple:
        elm = ELEMENTS[ion[0]]
        str_elm = elm.symbol
    else: 
        return ion_name( (ion['Z'], ion['ion']) )
        #raise ValueError('ionization.ion_name: Not ready for this input yet.')

    # Ion state
    if flg == 0: # Roman
        if nspace is None: nspace = 0
        str_ion = roman.toRoman(ion[1]) 
        spc = ' '*nspace
        outp = str_elm+spc+str_ion
    elif flg == 1: # LaTeX
        if ion[1] == 0:
            raise ValueError('ionization.ion_name: Not ready for this input yet.')
        elif ion[1] == 1:
            str_ion = '^0'
        elif ion[1] == 2:
            str_ion = '^{+}'
        elif ion[1] == 3:
            str_ion = '^{++}'
        else:
            str_ion = '^{+'+str(ion[1]-1)+'}'
        outp = '{\\rm '+str_elm+'}'+str_ion
    else:
        raise ValueError('ionization.ion_name: Not ready for this flg.')

    return outp


########################## ##########################
########################## ##########################
def photo_cross(Z, ion, E, datfil=None, silent=False):
    """ Estimate photo-ionization cross-section using Fit parameters
    from Verner et al. 1996, ApJ, 465, 487
    JXP on 04 Nov 2014

    Parameters
    ----------
    Z: Atomic number
    ion : Ionization state (1=Neutral)
    E : Energy to calculate at [eV]

    Returns
    -------
    sigma : Cross-section (cm^2)
    """

    # Read data
    if datfil == None:
        datfil = xa_path+'/data/atomic/verner96_photoion_table1.dat'
    dat = ascii.read(datfil)

    # Deal with Units
    if not isinstance(E,u.quantity.Quantity):
        if silent is False: print('photo_cross: Assuming eV for input energy')
        E = E * u.eV

    # Match
    #pdb.set_trace()
    mt = np.where((Z == dat['Z']) & (ion == dat['N']))[0]
    nmt = len(mt)
    if nmt == 0:
        raise ValueError('photo_cross: %d,%d pair not in our table' % (Z,ion))
    idx = mt[0]
    #
    x = E/(dat['E0'][idx]*u.eV) - dat['y0'][idx]
    y = np.sqrt(x**2 + dat['y1'][idx]**2)

    F = (((x-1.)**2 + dat['yw'][idx]**2) * y**(0.5*dat['P'][idx] - 5.5) * 
            (1 + np.sqrt(y/dat['ya'][idx]) )**(-1.*dat['P'][idx]))

    sigma = dat['s0'][idx] * F * 1e-18 * u.cm**2 

    # Energy threshold
    low = np.where(E < dat['Eth'][idx]*u.eV)[0]
    if len(low) > 0: sigma[low] = 0.

    return sigma

# Testing
if __name__ == '__main__':
    # Hydrogen

    print(photo_cross(1,1,13.6*u.eV))
    print(photo_cross(1,1,13.6*np.arange(1,11)*u.eV))
