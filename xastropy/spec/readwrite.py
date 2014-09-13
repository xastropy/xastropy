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
import barak
import numpy as np
from astropy.io import fits
from astropy.io import ascii 
import os, pdb

#### ###############################
#  Read Spectrum from FITS file
#  Return Barak-favored Table
#  from xastropy.spec import readwrite as xsr
#  sp = xsr.readspec('SDSSJ114435.54+095921.7_F.fits',outfil='SDSSJ114435.54+095921.7.fits')
#
def readspec(specfil, inflg=None, efil=None, outfil=None):
    from xastropy.spec import readwrite as rw
    from xastropy.files import general as xfg
    #from xastropy.plotting import x_guis as xpxg
    from astropy.table import Table
    from astropy.table import Column
    from barak import spec as bs
    #reload(bs)

    # Initialize
    dat = None
    chk = None
    if inflg == None:
        inflg = 0

    # Read header
    datfil = xfg.chk_for_gz(specfil,chk=chk)
    if chk == 0:
        print 'xastropy.spec.readwrite: File does not exist ', specfil
        return -1
    hdulist = fits.open(datfil)

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
            # Flux
            fx = hdulist[0].data.flatten()
            # Wavelength
            wave = rw.setwave(hdulist[0].header)
            # Error
            if efil == None:
                ipos = max(specfil.find('F.fits'),specfil.find('f.fits'))
                if ipos < 0: # No error array
                    sig = np.zeros(fx.size)
                else:
                    if specfil.find('F.fits') > 0:
                        efil = xfg.chk_for_gz(specfil[0:ipos]+'E.fits')
                    else:
                        efil = xfg.chk_for_gz(specfil[0:ipos]+'e.fits')
            if efil != None:
                sig=fits.getdata(efil)
            #pdb.set_trace()
            #xpxg.plot_1d_arrays(sig)
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

    # Generate Barak output
    try:
        co = fits.getdata(name+'_c.fits')
    except:
        co = np.nan*np.ones(len(fx))
    # Plot?
    if show_plot:
            xpxg.plot_1d_arrays(wave,fx,sig,co)
    sp = bs.Spectrum(wa=wave, fl=fx, er=sig, co=co, filename=specfil)

    # Write to disk?
    if outfil != None:
        sp.fits_write(outfil,overwrite=True)

    # Return 
    return sp

#### ###############################
#### ###############################
#  Set wavelength array using Header cards
def setwave(hdr):

    # Initialize
    SCL = 1.
    
    # Parse the header
    npix = hdr['NAXIS1'] 
    crpix1 = hdr['CRPIX1'] if 'CRPIX1' in hdr else 1.
    crval1 = hdr['CRVAL1'] if 'CRVAL1' in hdr else 1.
    cdelt1 = hdr['CDELT1'] if 'CDELT1' in hdr else 1.
    ctype1 = hdr['CTYPE1'] if 'CTYPE1' in hdr else None
    dcflag = hdr['DC-FLAG'] if 'DC-FLAG' in hdr else None

    # Generate
    if dcflag == 1:
        wave = SCL * 10.**(crval1 + ( cdelt1 * np.arange(npix) + 1. - crpix1) ) # Log

    # Return
    return wave

#### ###############################
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
