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
import subprocess

from astropy import units as u
from astropy.io import ascii, fits
from astropy.table import Table, Column

from xastropy.xutils import xdebug as xdb

# def bintab_to_table(fits_fil):
# def table_to_fits(table, outfil, compress=False, comment=None):

#
def bintab_to_table(fits_fil):
    ''' Read a binary FITS table into an astropy Table
    ''' 
    return Table( fits.open(fits_fil)[1].data )

def table_to_fits(table, outfil, compress=False, comment=None):
    ''' Write an astropy Table as a FITS binary table
    ---
    Parameters
    table: astropy.Table
    outfil: string
    compress: bool (False)
       gzip compress?
    comment: string 
       Comment to insert into the header
    '''

    # Comments
    if not comment is None:
        table.meta['comments'] = [comment]
    table.write(outfil, format='fits', overwrite=True)

    # Compress?
    if compress:
        subprocess.call(["gzip", "-f", outfil])
    '''
    # Uses the astropy.fits approach
    # Generate the header
    prihdr = fits.Header()
    if not comment is None:
        prihdr['COMMENT'] = comment
    prihdu = fits.PrimaryHDU(header=prihdr)

    # Table HDU
    table_hdu = fits.BinTableHDU.from_columns(np.array(table.filled()))

    # Write
    thdulist = fits.HDUList([prihdu, table_hdu])
    thdulist.writeto(outfil,clobber=True)

    # Compress?
    if compress:
        subprocess.call(["gzip", "-f", outfil])
    '''
