"""
#;+ 
#; NAME:
#; casbah.galaxy_data
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for generating data files for galaxies in CASBAH
#;   01-Jan-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp
from astropy.io import fits, ascii
from astropy import units as u 
from astropy.coordinates import SkyCoord
#from astropy import constants as const

from xastropy.xutils import xdebug as xdb

########################## ##########################
########################## ##########################
def galaxy_attrib():
    """ List of properties expected to be stored for CASBAH galaxies

    JXP on 01 Dec 2015
    """
    attrib = [ (str('RA'), float),                 # RA (J2000)
               (str('DEC'), float),                # DEC (J2000)
               (str('Z'), float),                  # Redshift
               (str('Z_ERR'), float),              # Redshift uncertainty
               (str('SDSS_MAG'), float, (5,)),     # ugriz photometry from SDSS
               (str('SDSS_MAGERR'), float, (5,)),    # ugriz photometry uncertainties
               (str('TELESCOPE'), '|S80'),            # Telescope(s) used
               (str('INSTRUMENT'), '|S80')            # Instrument(s) used
               ]

    #tmp = np.recarray( (1,), dtype=attrib)
    #tmp['TELESCOPE'] = 'BLAHHH'
    #xdb.set_trace()

    return attrib



########################## ##########################
########################## ##########################
def grab_sdss_galaxies(radec, radius=0.1*u.deg,outfil='tmp.fits'):
    """ Grab SDSS galaxies

    radius: float (0.1)
      Search radius in deg
    outfil: str ('tmp.fits')
      Name of output file for FITS table

    JXP on 01 Jan 2015
    """
    from astroquery.sdss import SDSS
    from astropy import coordinates as coords


    cC = coords.SkyCoord(ra=radec[0]*u.degree, dec=radec[1]*u.degree)

    photoobj_fs = ['ra', 'dec', 'objid', 'run', 'rerun', 'camcol', 'field']
    mags = ['petroMag_u', 'petroMag_g', 'petroMag_r', 'petroMag_i', 'petroMag_z'] 
    magsErr = ['petroMagErr_u', 'petroMagErr_g', 'petroMagErr_r', 'petroMagErr_i', 'petroMagErr_z']

    phot_catalog = SDSS.query_region(cC,spectro=True,radius=radius,
                                     photoobj_fields=photoobj_fs+mags+magsErr) # Unique
    spec_catalog = SDSS.query_region(cC,spectro=True, radius=radius) # Duplicates exist

    cgal = SkyCoord(ra=phot_catalog['ra']*u.degree, dec=phot_catalog['dec']*u.degree)
    cC   = SkyCoord(ra=radec[0]*u.degree, dec=radec[1]*u.degree)
    sepgal = cgal.separation(cC) #in degrees

    for ii,phot in enumerate(phot_catalog):
        # Spectro match?
        mt = np.where( np.abs(np.array(spec_catalog['ra']) - phot['ra']) < 1e-3)[0]
        if len(mt) == 0:
            print('No spectral match!')
            xdb.set_trace()
        # Duplicate?
        mt = np.where( cgal.separation(cgal[ii]).to('arcsec') < 1.*u.Unit('arcsec'))[0] 
        if len(mt) != 1:
            print('Two photometric with same RA/DEC')
            xdb.set_trace()

    nobj = len(phot_catalog)
    spec_hdus = SDSS.get_spectra(matches=spec_catalog)

    # Generate output table
    attribs = galaxy_attrib()
    npix = 5000 #len( spec_hdus[0][1].data.flux )
    spec_attrib = [(str('FLUX'), np.float32, (npix,)),
                   (str('SIG'), np.float32, (npix,)),
                   (str('WAVE'), np.float64, (npix,))]
    tbl = np.recarray( (nobj,), dtype=attribs+spec_attrib)

    tbl['RA'] = phot_catalog['ra']
    tbl['DEC'] = phot_catalog['dec']
    tbl['TELESCOPE'] = str('SDSS 2.5-M')

    # Deal with spectra separately (for now)
    npix = 5000 #len( spec_hdus[0][1].data.flux )
    
    for idx,obj in enumerate(phot_catalog): 
        #print('idx = {:d}'.format(idx))

        # Grab spectra (there may be duplicates)
        mt = np.where( np.abs(np.array(spec_catalog['ra']) - obj['ra']) < 1e-3)[0]
        if len(mt) > 1:
            # Use BOSS if you have it
            mmt = np.where( spec_catalog[mt]['instrument'] == 'BOSS')[0]
            if len(mmt) > 0:
                mt = mmt[0]
            else:
                mt = mt[0]
        else:
            mt = mt[0]

        tbl[idx]['INSTRUMENT'] = spec_catalog[mt]['instrument']
        spec = spec_hdus[mt][1].data
        npp = len(spec.flux)
        tbl[idx]['FLUX'][0:npp] = spec.flux
        sig = np.zeros(npp)
        gdi = np.where(spec.ivar > 0.)[0]
        if len(gdi) > 0:
            sig[gdi] = np.sqrt( 1./spec.ivar[gdi] )
        tbl[idx]['SIG'][0:npp] = sig
        tbl[idx]['WAVE'][0:npp] = 10.**spec.loglam

        # Redshifts
        meta = spec_hdus[mt][2].data
        for attrib in ['Z','Z_ERR']:
            tbl[idx][attrib] = meta[attrib]

        # Fill in rest
        tbl[idx].SDSS_MAG = np.array( [obj[phot] for phot in mags])
        tbl[idx].SDSS_MAGERR = np.array( [obj[phot] for phot in magsErr])

    # Could clip on redshift to excise stars/quasars

    # Write to FITS file
                
    prihdr = fits.Header()
    prihdr['COMMENT'] = 'SDSS Spectra'
    prihdu = fits.PrimaryHDU(header=prihdr)

    tbhdu = fits.BinTableHDU(tbl)

    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(outfil,clobber=True)

    print('Wrote SDSS table to {:s}'.format(outfil))

    
# ################
if __name__ == "__main__":

    flg_fig = 0 
    flg_fig += 2**0  # SDSS search
    
    # XSpec
    if (flg_fig % 2**1) >= 2**0:
        radec = (212.34957,26.30585)
        grab_sdss_galaxies(radec, radius=0.1*u.degree) 
