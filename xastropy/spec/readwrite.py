"""
#;+ 
#; NAME:
#; readwrite
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for read/write of Spectra
#;   07-Sep-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""

# Import libraries
import numpy as np
from astropy.io import fits
from astropy.io import ascii 
from astropy.table import Table
from astropy.table import Column
import pdb

#### ###############################
#  Read Spectrum from FITS file
#  Return Barak-favored Table
#
def readspec(file, inflg=None):
    from xastropy.spec import readwrite as rw

    # Initialize
    dat = None
    if inflg == None:
        inflg = 0

    # Read header
    hdulist = fits.open(file)

    ## #################
    # Binary FITS table?
    if hdulist[0].header['NAXIS'] == 0:
        # Flux 
        flux_tags = ['SPEC','FLUX','FLAM','FX']
        fx = rw.get_table_column(flux_tags, hdulist)
        if fx == None:
            print 'spec.readwrite: Binary FITS Table but no Flux tag'
            return
        # Error
        sig_tags = ['ERROR','ERR','SIGMA_FLUX','FLAM_SIG']
        sig = rw.get_table_column(sig_tags, hdulist)
        if sig == None:
            ivar_tags = ['IVAR']
            if ivar == None:
                print 'spec.readwrite: Binary FITS Table but no error tags'
                return
            else: 
                sig = fltarr(ivar.size)
                gdi = np.where( ivar > 0.)
                sig[gdi] = sqrt(1./strct.ivar[gdi])
        # Wavelength
        wave_tags = ['WAVE','WAVELENGTH','LAMBDA']
        wave = rw.get_table_column(wave_tags, hdulist)
        if wave == None:
            print 'spec.readwrite: Binary FITS Table but no wavelength tag'
            return
    elif hdulist[0].header['NAXIS'] == 1: # Data in the zero extension
        # Look for wavelength info
        if 'CRVAL1' in hdulist[0].header.keys():
            print 'spec.readwrite: Not ready for this yet!'
            return
        else:  # ASSUMING MULTI-EXTENSION
            if len(hdulist) <= 2:
                print 'spec.readwrite: No wavelength info but only 2 extensions!'
                return
            fx = hdulist[0].data.flatten()
            sig = hdulist[1].data.flatten()
            wave = hdulist[2].data.flatten()
    else:  # Should not be here
        print 'spec.readwrite: Looks like an image'
        return dat

    # Generate output

    # Return 
    return fx, sig, wave


#### ###############################
#  Grab values from the Binary FITS Table
def get_table_column(tags, hdulist):
    dat = None
    ii = 0
    #pdb.set_trace()
    while(ii < len(tags)):
        if tags[ii] in hdulist[1].columns.names: 
            dat = hdulist[1].data[tags[ii]]
            break  # Break with first hit
        else:
            ii = ii + 1
    # Return
    return dat.flatten()
