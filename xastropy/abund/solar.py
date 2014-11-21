"""
#;+ 
#; NAME:
#; solar
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for solar abundances
#;   21-Nov-2014 by JXP
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
from xastropy.xutils import xdebug as xdb

from astropy.utils.misc import isiterable

# Path for xastropy
xa_path = imp.find_module('xastropy')[1]

#def ion_name(ion):
#def photo_cross(Z, ion, E, datfil=None, silent=False):

########################## ##########################
########################## ##########################
def abund(Z,dat_file=None,table=None):
    """ Report back the solar abundance

    Parameters
    ----------
    Z: int or string (can be an array of either)
      Atomic number or name

    Returns
    -------
    out_abnd : float (scalar or an array)
      Solar abundance.  Meteoritic if available.

    JXP on 21 Nov 2014
    """
    if table is None:
        # Data file
        if dat_file is None:
            dat_file = xa_path+'/data/abund/solar_Apslund09.dat'
    
        # Read table
        names=('name', 'abund', 'Z')
        table = ascii.read(dat_file, format='no_header', names=names) 


    # Iterate?
    if isiterable(Z): 
        out_abnd = []
        for iZ in Z:
            out_abnd.append(abund(iZ,table=table))
        out_abnd = np.array(out_abnd)
    else: 
        if isinstance(Z,int):
            idx = np.where( table['Z'] == Z )[0]
            if len(idx) == 0:
                raise ValueError('abund.solar.abund: Z={:d} not in {:s}'.format(Z,dat_file))
            out_abnd = table['abund'][idx[0]]
        else:
            raise ValueError('abund.solar.abund: Not ready for this input yet.')

    return out_abnd



# Testing
if __name__ == '__main__':

    # Abundances
    print('C = {:g}'.format(abund(6)))
    print('O = {:g}'.format(abund(8)))
    tmp = abund([6,8])
    print('C = {:g} and O = {:g}'.format(tmp[0], tmp[1]))
