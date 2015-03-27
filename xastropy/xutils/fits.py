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

#
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
