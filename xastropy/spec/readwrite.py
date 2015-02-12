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

from __future__ import print_function, absolute_import, division, unicode_literals

# Import libraries
import numpy as np
from astropy.io import fits
from astropy.io import ascii 
from astropy.nddata import StdDevUncertainty
import os

from specutils.io import read_fits as spec_read_fits
from specutils.spectrum1d import Spectrum1D

from xastropy.xutils import xdebug as xdb

#### ###############################
#  Read Spectrum1D from FITS file
#  from xastropy.spec import readwrite as xsr
#  sp = xsr.readspec('SDSSJ114435.54+095921.7_F.fits',outfil='SDSSJ114435.54+095921.7.fits')
#
def readspec(specfil, inflg=None, efil=None, outfil=None, show_plot=0,
             use_barak=False, verbose=False, flux_tags=None, sig_tags=None, multi_ivar=False):
    ''' 
    multi_ivar: Bool (False)
      BOSS format of  flux, ivar, log10(wave) in multi-extension FITS
    '''
    from xastropy.spec import readwrite as rw
    from xastropy.files import general as xfg
    #from xastropy.plotting import x_guis as xpxg
    from astropy.table import Table
    from astropy.table import Column
    #reload(bs)

    # Initialize
    dat = None
    chk = None
    if inflg == None:
        inflg = 0

    # Read header
    datfil = xfg.chk_for_gz(specfil,chk=chk)
    if chk == 0:
        print('xastropy.spec.readwrite: File does not exist ', specfil)
        return -1
    hdulist = fits.open(os.path.expanduser(datfil))

    ## #################
    # Binary FITS table?
    if hdulist[0].header['NAXIS'] == 0:
        # Flux 
        if flux_tags is None:
            flux_tags = ['SPEC', 'FLUX','FLAM','FX', 'FLUXSTIS']
        fx, fx_tag = rw.get_table_column(flux_tags, hdulist)
        #xdb.set_trace()
        if fx is None:
            print('spec.readwrite: Binary FITS Table but no Flux tag')
            return
        # Error
        if sig_tags is None:
            sig_tags = ['ERROR','ERR','SIGMA_FLUX','FLAM_SIG', 'SIGMA_UP', 'ERRSTIS']
        sig, sig_tag = rw.get_table_column(sig_tags, hdulist)
        if sig is None:
            ivar_tags = ['IVAR']
            ivar, ivar_tag = rw.get_table_column(ivar_tags, hdulist)
            if ivar is None:
                print('spec.readwrite: Binary FITS Table but no error tags')
                return
            else: 
                sig = np.zeros(ivar.size)
                gdi = np.where( ivar > 0.)[0]
                sig[gdi] = np.sqrt(1./ivar[gdi])
        # Wavelength
        wave_tags = ['WAVE','WAVELENGTH','LAMBDA','LOGLAM', 'WAVESTIS']
        wave, wave_tag = rw.get_table_column(wave_tags, hdulist)
        if wave_tag == 'LOGLAM':
            wave = 10.**wave
        if wave is None:
            print('spec.readwrite: Binary FITS Table but no wavelength tag')
            return
    elif hdulist[0].header['NAXIS'] == 1: # Data in the zero extension
        # Look for wavelength info
        try:
            ctype1 = hdulist[0].header['CTYPE1']
        except KeyError:
            ctype1 = ''
        if ('CRVAL1' in hdulist[0].header.keys()) and (multi_ivar is False) and (not ctype1 in ['RA---TAN']):
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
                    efil=os.path.expanduser(efil)

            # Generate Spectrum1D
            spec1d = spec_read_fits.read_fits_spectrum1d(os.path.expanduser(datfil), efil=efil)
            #spec1d = spec_read_fits.read_fits_spectrum1d(os.path.expanduser(datfil))

        else:  # ASSUMING MULTI-EXTENSION
            if len(hdulist) <= 2:
                print('spec.readwrite: No wavelength info but only 2 extensions!')
                return
            fx = hdulist[0].data.flatten()
            sig = hdulist[1].data.flatten()
            wave = hdulist[2].data.flatten()
            if multi_ivar is True:
                ivar = np.zeros(len(sig))
                gdp = np.where(sig > 0.)[0]
                ivar[gdp] = np.sqrt(1./sig[gdp])
                wave = 10.**wave
    else:  # Should not be here
        print('spec.readwrite: Looks like an image')
        return dat

    # Generate, as needed
    if 'spec1d' not in locals():
        # Generate
        spec1d = Spectrum1D.from_array(wave, fx, uncertainty=StdDevUncertainty(sig))

    spec1d.filename = specfil

    # Continuum?
    try:
        co = fits.getdata(name+'_c.fits')
    except:
        try:
            npix = len(fx)
        except UnboundLocalError:
            npix = len(spec1d.flux)
        co = np.nan*np.ones(npix)
    # Generate a Barak Spectrum Class?
    hd = hdulist[0].header
    if use_barak is True:
        # Barak
        from barak import spec as bs
        spec1d = bs.Spectrum(wa=wave, fl=fx, er=sig, co=co, filename=specfil)
        spec1d.header = hd

    # Plot?
    if show_plot:
            xpxg.plot_1d_arrays(wave,fx,sig,co)

    # Write to disk? Unlikely
    if outfil != None:
        if use_barak is True:
            spec1d.fits_write(outfil,overwrite=True)
        else:
            xdb.set_trace() # Not ready

    # Return 
    return spec1d


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
    if (dcflag == 1) or (cdelt1 < 1e-4):
        wave = SCL * 10.**(crval1 + ( cdelt1 * np.arange(npix) + 1. - crpix1) ) # Log
    xdb.set_trace()

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
    if dat is not None:
        return dat.flatten(), tags[ii]
    else: 
        return dat, 'NONE'








#### ###############################
# Testing
if __name__ == '__main__':
    flg_test = 0
    flg_test += 1 # MagE
    flg_test += 2**1 # LRIS LowRedux

    # Standard log-linear read (MagE)
    if (flg_test % 2**1) >= 2**0:
        fil = '~/PROGETTI/LLSZ3/data/normalize/UM669_nF.fits'
        #efil = '~ers/xavier/PROGETTI/LLSZ3/data/normalize/UM669_nE.fits'
        myspec = readspec(fil)
        #xdb.xplot(myspec.dispersion, myspec.flux)
        xdb.xplot(myspec.dispersion, myspec.flux, myspec.sig)
        #xdb.xplot(myspec.dispersion, myspec.flux, myspec.uncertainty.array)

    # LowRedux
    if (flg_test % 2**2) >= 2**1:
        fil = '/Users/xavier/Dropbox/QSOPairs/data/LRIS_redux/SDSSJ234704.25+150146.4_F.fits'
        myspec = readspec(fil)
        xdb.xplot(myspec.dispersion, myspec.flux, myspec.uncertainty.array)
        #xdb.set_trace()
