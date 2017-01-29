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
import pdb
from astropy.io import fits, ascii
from astropy import units as u 
from astropy.table import Table

#from xastropy.xutils import xdebug as xdb


def get_filename(field, ftype, **kwargs):
    '''Generate a CASBAH file given field and type

    Parameters:
    -----------
    field: tuple
      (Name, ra, dec)
    ftype: str
      Type of filename desired
    '''
    path = os.getenv('CASBAH_DIR')
    subdir = 'Galaxies/'
    if ftype == 'MULTI_OBJ':
        filename = path+'/'+subdir+field[0]+'/'+field[0]+'_MULTI_OBJ.ascii'
    elif ftype == 'FIELD_PATH':
        filename = path+'/'+subdir+field[0]
    elif ftype == 'TARGETS':
        filename = path+'/'+subdir+field[0]+'/'+field[0]+'_targets.ascii'
    elif ftype == 'SDSS':
        filename = path+'/'+subdir+field[0]+'/'+field[0]+'_SDSS.fits'
    elif ftype == 'HECTOSPEC':
        filename = path+'/'+subdir+field[0]+'/'+field[0]+'_HECTOSPEC.fits'
    elif ftype == 'DEIMOS_TARG_FIG':
        filename = path+'/'+subdir+field[0]+'/'+field[0]+'_deimostarg.pdf'
    elif ftype == 'HECTO_TARG_FIG':
        filename = path+'/'+subdir+field[0]+'/'+field[0]+'_hectotarg.pdf'
    elif ftype == 'IMAGING':
        if 'rename_deimos' in kwargs.keys():
            img_fil = kwargs['orig_file']
            ipos = img_fil.find('.mos')
            filter=img_fil[ipos-1:ipos]
            filename = path+'/'+subdir+field[0]+'/'+field[0]+'_IMG_{:s}_{:s}.fits'.format('LBC', filter)
    else:
        raise ValueError('Not ready for this ftype: {:s}'.format(ftype))
    # Return
    return filename


def load(field, ftype, **kwargs):
    """ Load input file type
    Parameters
    ----------
    field : tuple
      name,ra,dec
    ftype
    kwargs

    Returns
    -------
    data : Table or ndarray

    """
    # Get filename first
    filename = get_filename(field, ftype, **kwargs)
    # Load
    if ftype == 'TARGETS':
        data = Table.read(filename, format='ascii.fixed_width')
    else:
        raise ValueError('Not ready for this ftype: {:s}'.format(ftype))
    # Return
    return data
