'''
#;+ 
#; NAME:
#; sdss.qso
#;    Version 1.1
#;
#; PURPOSE:
#;   Class for SDSS QSO
#;     2015 Written by JXP
#;-
#;------------------------------------------------------------------------------
'''

# Import libraries
import numpy as np
import os

from astropy.table import QTable, Column
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.units import Quantity

from xastropy.obs import radec as xor
from xastropy.xutils import xdebug as xdb

class SdssQso(object):
    '''Class to handle a single SDSS Quasar

    Parameters:
    ----------
    coord: SkyCoord, optional
      RA/Dec of the sightline
    z: float, optional
      Emission redshift
    database: SdssQuasars class, optional
      Required for grabbing data, etc.
    '''
    # Init
    def __init__(self, coord=None, z=0., database=None, verbose=True):
        # Init
        if coord is None:
            radec = (0.*u.deg, 0.*u.deg)
            self.coord = SkyCoord(ra=radec[0], dec=radec[0])
        else:
            self.coord = coord
        self.z = z 
        self.verbose = verbose
        self.database = database
        # None init
        self._specfil = None

    def get_specfil(self):
        '''Parse the SDSS spectrum file
        Requires a link to the database Class
        '''
        if self.database is None:
            raise IOError('SdssQso: Need to be linked to an SDSS Database')
        # Generate file name (DR4 is different)
        pnm = '{0:04d}'.format(
            self.database._data[self.database.index]['PLATE'])
        #fnm = '{0:04d}'.format(
        #    self.database._data[self.database.index]['FIBERID'])
        fnm = '{0:03d}'.format(
            self.database._data[self.database.index]['FIBERID'])
        mjd = str(self.database._data[self.database.index]['MJD'])
        sfil = self.database._datdir+pnm+'/1d/'+'spSpec-' 
        # Finish
        self._specfil = sfil+mjd+'-'+pnm+'-'+fnm+'.fit'  # Is usually gzipped

    def load_spec(self):
        '''Input the Spectrum
        '''
        from linetools.spectra.xspectrum1d import XSpectrum1D
        if self._specfil is None:
            self.get_specfil()
        #
        if self.verbose:
            print('SdssQso: Loading spectrum from {:s}'.format(self._specfil))
        self.spec = XSpectrum1D.from_file(self._specfil)

    def __repr__(self):
        ''' For printing
        '''
        return '[{:s}: {:s} {:s}, z={:g}]'.format(self.__class__.__name__, 
            self.coord.ra.to_string(unit=u.hour,sep=':',pad=True), 
            self.coord.dec.to_string(sep=':',pad=True,alwayssign=True), self.z)
