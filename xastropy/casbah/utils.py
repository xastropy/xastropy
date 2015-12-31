"""
#;+ 
#; NAME:
#; utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for CASBAH utilities
#;   02-Jan-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp, glob, copy
from astropy.io import fits, ascii
from astropy import units as u 
from astropy.coordinates import SkyCoord

from xastropy.obs import radec as xra

from xastropy.xutils import xdebug as xdb

def get_filename(field, ftype):
    '''Generate a CASBAH file given field and type

    Parameters:
    -----------
    field: tuple
      (Name, ra, dec)
    ftype: str
      Type of filename desired
    '''
    path = os.getenv('CASBAH_GALAXIES')
    if ftype == 'MULTI_OBJ':
        filename = path+'/'+field[0]+'/'+field[0]+'_MULTI_OBJ.ascii'
    elif ftype == 'FIELD_PATH':
        filename = path+'/'+field[0]
    elif ftype == 'TARGETS':
        filename = path+'/'+field[0]+'/'+field[0]+'_targets.ascii'
    elif ftype == 'SDSS':
        filename = path+'/'+field[0]+'/'+field[0]+'_SDSS.fits'
    elif ftype == 'HECTOSPEC':
        filename = path+'/'+field[0]+'/'+field[0]+'_HECTOSPEC.fits'
    elif ftype == 'DEIMOS_TARG_FIG':
        filename = path+'/'+field[0]+'/'+field[0]+'_deimostarg.pdf'
    elif ftype == 'HECTO_TARG_FIG':
        filename = path+'/'+field[0]+'/'+field[0]+'_hectotarg.pdf'
    else:
        raise ValueError('Not ready for this ftype: {:s}'.format(ftype))
    # Return
    return filename

