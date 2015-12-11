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
from astropy.io import fits
from astropy import units as u 
from astropy.table.table import Table

from xastropy.xutils import xdebug as xdb


def galaxy_attrib():
    """ List of properties expected to be stored for CASBAH galaxies

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
    # Return
    return attrib


def grab_sdss_spectra(radec, radius=0.1*u.deg, outfil=None,
                      debug=False, maxsep=None, timeout=600., zmin=None):
    """ Grab SDSS spectra

    Parameters
    ----------
    radec : tuple
      RA, DEC in deg
    radius : float, optional (0.1*u.deg)
      Search radius -- Astroquery actually makes a box, not a circle
    timeout : float, optional
      Timeout limit for connection with SDSS
    outfil : str ('tmp.fits')
      Name of output file for FITS table
    maxsep : float (None) :: Mpc
      Maximum separation to include 
    zmin : float (None)
      Minimum redshift to include

    Returns
    -------
    tbl : Table

    """
    from astroquery.sdss import SDSS
    from astropy import coordinates as coords
    from astropy.cosmology import Planck13 as cosmo 
    from astropy.coordinates import SkyCoord

    cC = coords.SkyCoord(ra=radec[0], dec=radec[1])

    # Query
    photoobj_fs = ['ra', 'dec', 'objid', 'run', 'rerun', 'camcol', 'field']
    mags = ['petroMag_u', 'petroMag_g', 'petroMag_r', 'petroMag_i', 'petroMag_z'] 
    magsErr = ['petroMagErr_u', 'petroMagErr_g', 'petroMagErr_r', 'petroMagErr_i', 'petroMagErr_z']

    phot_catalog = SDSS.query_region(cC,spectro=True,radius=radius, timeout=timeout,
                                     photoobj_fields=photoobj_fs+mags+magsErr) # Unique
    spec_catalog = SDSS.query_region(cC,spectro=True, radius=radius, timeout=timeout) # Duplicates exist
    nobj = len(phot_catalog)

    #
    print('grab_sdss_spectra: Found {:d} sources in the search box.'.format(nobj))

    # Coordinates
    cgal = SkyCoord(ra=phot_catalog['ra']*u.degree, dec=phot_catalog['dec']*u.degree)
    sgal = SkyCoord(ra=spec_catalog['ra']*u.degree, dec=spec_catalog['dec']*u.degree)
    sepgal = cgal.separation(cC) #in degrees

    # Check for problems and parse z
    zobj = np.zeros(nobj)
    idx, d2d, d3d = coords.match_coordinates_sky(cgal, sgal, nthneighbor=1)
    if np.max(d2d)*3600. > 1.*u.Unit('arcsec'):
        print('No spectral match!')
        xdb.set_trace()
    else:
        zobj = spec_catalog['z'][idx]

    idx, d2d, d3d = coords.match_coordinates_sky(cgal, cgal, nthneighbor=2)
    if np.min(d2d.to('arcsec')) < 1.*u.Unit('arcsec'):
        print('Two photometric sources with same RA/DEC')
        xdb.set_trace()

    #xdb.set_trace()


    # Cut on Separation
    if not maxsep is None:
        print('grab_sdss_spectra: Restricting to {:g} Mpc separation.'.format(maxsep))
        sepgal_kpc = cosmo.kpc_comoving_per_arcmin(zobj) * sepgal.to('arcmin')
        sepgal_mpc = sepgal_kpc.to('Mpc')
        gdg = np.where( sepgal_mpc < (maxsep * u.Unit('Mpc')))[0]
        phot_catalog = phot_catalog[gdg]
        #xdb.set_trace()

    nobj = len(phot_catalog)
    print('grab_sdss_spectra: Grabbing data for {:d} sources.'.format(nobj))

    # Grab Spectra from SDSS

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
        mt = np.where( sgal.separation(cgal[idx]).to('arcsec') < 1.*u.Unit('arcsec'))[0] 
        if len(mt) > 1:
            # Use BOSS if you have it
            mmt = np.where( spec_catalog[mt]['instrument'] == 'BOSS')[0]
            if len(mmt) > 0:
                mt = mt[mmt[0]]
            else:
                mt = mt[0]
        elif len(mt) == 0:
            xdb.set_trace()
        else:
            mt = mt[0]

        # Grab spectra
        spec_hdus = SDSS.get_spectra(matches=Table(spec_catalog[mt]))

        tbl[idx]['INSTRUMENT'] = spec_catalog[mt]['instrument']
        spec = spec_hdus[0][1].data
        npp = len(spec.flux)
        tbl[idx]['FLUX'][0:npp] = spec.flux
        sig = np.zeros(npp)
        gdi = np.where(spec.ivar > 0.)[0]
        if len(gdi) > 0:
            sig[gdi] = np.sqrt( 1./spec.ivar[gdi] )
        tbl[idx]['SIG'][0:npp] = sig
        tbl[idx]['WAVE'][0:npp] = 10.**spec.loglam

        # Redshifts
        meta = spec_hdus[0][2].data
        for attrib in ['Z','Z_ERR']:
            tbl[idx][attrib] = meta[attrib]

        if debug:
            sep_to_qso = cgal[idx].separation(cC).to('arcmin')
            print('z = {:g}, Separation = {:g}'.format(tbl[idx].Z, sep_to_qso))
            xdb.set_trace()

        # Fill in rest
        tbl[idx].SDSS_MAG = np.array( [obj[phot] for phot in mags])
        tbl[idx].SDSS_MAGERR = np.array( [obj[phot] for phot in magsErr])

    # Clip on redshift to excise stars/quasars
    if zmin is not None:
        gd = np.where(tbl['Z'] > zmin)[0]
        tbl = tbl[gd]

    # Write to FITS file
    if outfil is not None:
        prihdr = fits.Header()
        prihdr['COMMENT'] = 'SDSS Spectra'
        prihdu = fits.PrimaryHDU(header=prihdr)

        tbhdu = fits.BinTableHDU(tbl)

        thdulist = fits.HDUList([prihdu, tbhdu])
        thdulist.writeto(outfil,clobber=True)

    print('Wrote SDSS table to {:s}'.format(outfil))
    return tbl

    
# ################
if __name__ == "__main__":

    flg_fig = 0 
    flg_fig += 2**0  # SDSS search
    
    # XSpec
    if (flg_fig % 2**1) >= 2**0:
        radec = (212.34957*u.deg,26.30585*u.deg)
        grab_sdss_spectra(radec, radius=1.*u.degree/12.) 
