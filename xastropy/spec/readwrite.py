"""
#;+ 
#; NAME:
#; io  (used to be readwrite)
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
import os

from astropy.io import fits, ascii
from astropy.nddata import StdDevUncertainty
from astropy import units as u
from astropy.io.fits.fitsrec import FITS_rec
from astropy.io.fits.hdu.table import BinTableHDU

from specutils.io import read_fits as spec_read_fits

from xastropy.xutils import xdebug as xdb
from xastropy.spec.utils import XSpectrum1D

#### ###############################
#  Generate Spectrum1D from FITS file
#
def readspec(specfil, inflg=None, efil=None, outfil=None, show_plot=0,
             use_barak=False, verbose=False, flux_tags=None, sig_tags=None, multi_ivar=False):
    ''' 
    specfil: string or Table
    multi_ivar: Bool (False)
      BOSS format of  flux, ivar, log10(wave) in multi-extension FITS
    '''
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

    # Check specfil type
    if type(specfil) is Table:
        datfil = 'None'
        # Dummy hdulist
        hdulist = [fits.PrimaryHDU(), specfil]
    else:
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
            flux_tags = ['SPEC', 'FLUX','FLAM','FX', 'FLUXSTIS', 'FLUX_OPT', 'fl']
        fx, fx_tag = get_table_column(flux_tags, hdulist)
        #xdb.set_trace()
        if fx is None:
            print('spec.readwrite: Binary FITS Table but no Flux tag')
            return
        # Error
        if sig_tags is None:
            sig_tags = ['ERROR','ERR','SIGMA_FLUX','FLAM_SIG', 'SIGMA_UP', 'ERRSTIS', 'FLUXERR', 'er']
        sig, sig_tag = get_table_column(sig_tags, hdulist)
        if sig is None:
            ivar_tags = ['IVAR', 'IVAR_OPT']
            ivar, ivar_tag = get_table_column(ivar_tags, hdulist)
            if ivar is None:
                print('spec.readwrite: Binary FITS Table but no error tags')
                return
            else: 
                sig = np.zeros(ivar.size)
                gdi = np.where( ivar > 0.)[0]
                sig[gdi] = np.sqrt(1./ivar[gdi])
        # Wavelength
        wave_tags = ['WAVE','WAVELENGTH','LAMBDA','LOGLAM', 'WAVESTIS', 'WAVE_OPT', 'wa']
        wave, wave_tag = get_table_column(wave_tags, hdulist)
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
            spec1d = spec_read_fits.read_fits_spectrum1d(os.path.expanduser(datfil),
                                                         dispersion_unit='AA',
                                                         efil=efil)
            xspec1d = XSpectrum1D.from_spec1d(spec1d)

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
    if 'xspec1d' not in locals():
        # Give Ang as default
        if not hasattr(wave, 'unit'):
            dispersion_unit = u.AA
        else:
            dispersion_unit = None
        #xdb.set_trace()
        xspec1d = XSpectrum1D.from_array(u.Quantity(wave), u.Quantity(fx),
                                         uncertainty=StdDevUncertainty(sig),
                                         dispersion_unit=dispersion_unit)

    xspec1d.filename = specfil

    # Continuum?
    try:
        co = fits.getdata(name+'_c.fits')
    except:
        try:
            npix = len(fx)
        except UnboundLocalError:
            npix = len(xspec1d.flux)
        co = np.nan*np.ones(npix)
    # Generate a Barak Spectrum Class?
    hd = hdulist[0].header
    if use_barak is True:
        # Barak
        raise ValueError('Avoid!')
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
    return xspec1d


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
#  Grab values from the Binary FITS Table or Table
def get_table_column(tags, hdulist):
    dat = None
    ii = 0
    # Use Table
    if type(hdulist[1]) is BinTableHDU:
        tab = Table(hdulist[1])
    else:
        tab = hdulist[1]

    # Grab
    for tag in tags:
        if tag in tab.dtype.names: 
            dat = tab[tag]
            break  # Break with first hit

    '''
    For BinTableHDU (deprecated)
    while(ii < len(tags)):
        if tags[ii] in hdulist[1].columns.names: 
            dat = hdulist[1].data[tags[ii]]
            break  # Break with first hit
        else:
            ii = ii + 1
    '''
    # Return
    if dat is not None:
        return dat.flatten(), tag
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
        fil = '/Users/xavier/Dropbox/QSOPairs/data/LRIS_redux/SDSSJ231254.65-025403.1_b400_F.fits.gz'
        myspec = readspec(fil)
        xdb.xplot(myspec.dispersion, myspec.flux, myspec.uncertainty.array)
        #xdb.set_trace()
