"""
#;+ 
#; NAME:
#; fits
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for FITS utilities
#;   10-Dec-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import gzip
import subprocess, glob

from astropy import units as u
from astropy.io import ascii, fits
from astropy.table import QTable, Table, Column

from xastropy.xutils import xdebug as xdb

# def bintab_to_table(fits_fil,exten=1, silent=False):
# def table_to_fits(table, outfil, compress=False, comment=None):

#
def bintab_to_table(fits_fil,exten=1, silent=False):
    ''' Read a binary FITS table into an astropy QTable

    Parameters
    ---------
    fits_fil: sting
      File
    exten: int (1)
      Extension of the FITS file
    ''' 
    # Look for file
    fil = glob.glob(fits_fil+'*')
    if len(fil) == 1:
        if not silent:
            print('x.fits.bintab_to_table: Reading {:s}'.format(fil[0]))
        infil = fil[0]
    elif len(fil) == 0:
        raise IOError('x.fits.bintab_to_table: No file {:s}'.format(fits_fil))
    else:
        # Taking specified one
        infil = fits_fil

    return QTable( fits.open(infil)[exten].data )

#
def table_to_fits(itable, outfil, compress=False, comment=None):
    ''' Write an astropy Table as a FITS binary table
    ---
    Parameters
    ---------
    table: astropy.Table
    outfil: string
    compress: bool (False)
       gzip compress?
    comment: string 
       Comment to insert into the header
    '''

    # Comments
    table = Table(itable)
    if not comment is None:
        table.meta['comments'] = [comment]
    table.write(outfil, format='fits', overwrite=True)

    # Compress?
    if compress:
        subprocess.call(["gzip", "-f", outfil])

def write_quick_fits(arr_list, outfil, clobber=True):
    ''' Write a list of arrays to a FITS file
    Parameters
    ---------
    arr_list: list of ndarray
    outfil: str
    '''
    hdulist = None
    for arr in arr_list:
        if hdulist is None:
            hdulist = fits.HDUList([fits.PrimaryHDU(arr)])
        else:
            hdulist.append(fits.ImageHDU(arr))
    # Write
    hdulist.writeto(outfil, clobber=clobber)
    print('Wrote: {:s}'.format(outfil))


