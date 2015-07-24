"""
#;+ 
#; NAME:
#; ginga
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for Ginga stuff
#;   Jul-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import gzip, os
import subprocess, glob

from astropy import units as u
from astropy.io import ascii, fits
from astropy.table import QTable, Table, Column

from ginga.util import grc

#from xastropy.xutils import xdebug as xdb


#
def show_fits(fits_file, host='localhost', port=9000, **kwargs):
    ''' Read a binary FITS file into an active Ginga window

    Parameters
    ---------
    fits_file: str
      Filename  
    host: str, optional
      Name of the host; Default='localhost'
    port: int, optional
      Value of the port; Default=9000
    ''' 
    # Full path?
    if fits_file[0] != '/':
        fil = os.getcwd()+'/'+fits_file
    else:
        fil = fits_file
    # Check for file
    if not os.path.isfile(fil):
        raise ValueError('File={:s} not found!'.format(fil))
    # Connect to ginga RC (should have error checking here)
    client = grc.RemoteClient(host, port)
    # Connect to ginga widget
    method = getattr(client, 'ginga')
    # Set up args 
    command = 'load_file'
    args = [command,fil]
    # Execute
    res = method(*args, **kwargs)

def show_img(img,tmp_dir='./',tmp_file='tmp_ginga.fits', **kwargs):
    '''Push an image to an Ginga window
    Currently uses the kludge of writing a temp file
    '''
    # Generate filename with full path
    if tmp_file[0] != '/':
        kludge_fil = os.getcwd()+'/'+tmp_dir+tmp_file
    else:
        kludge_fil = tmp_file
    # Error checking
    if len(img.shape) != 2:
        raise ValueError('Image is wrong shape!')
    # Generate hdu
    hdu = fits.PrimaryHDU(img)
    hdulist = fits.HDUList([hdu])
    # Write
    print('WARNING: Writing kludge file to {:s}'.format(kludge_fil))
    hdulist.writeto(kludge_fil,clobber=True)
    # Show
    show_fits(kludge_fil, **kwargs)
