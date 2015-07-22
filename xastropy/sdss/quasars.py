'''
#;+ 
#; NAME:
#; sdss.quasars
#;    Version 1.1
#;
#; PURPOSE:
#;   Basic routines for dealing with SDSS quasar database
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

class SdssQuasars(object):
    '''Class to handle a data release of SDSS quasars

    Parameters:
    ----------
    version: str, optional
       'DR7'   :: JXP version of DR7 [DEFAULT]
    '''
    # Init
    def __init__(self, version=None, verbose=True):
        # Set version
        if version is None:
            self._version = 'DR7'
        else:
            self._version = version
        self.verbose = verbose
        # Paths
        self._path = os.getenv('SDSSPATH')+'/'+self._version+'_QSO/' 
        self._datdir = self._path+'spectro/1d_26/'

        # Grab Summary FITS file
        self._summf = path+version.lower()+'_qso.fits.gz'
        if self.verbose:
            print('SDSS_QUASAR: Using summary file {:s}'.format(summf))
        self._data = QTable.read(summf)


    #### ###############################
    def get_qso(self,inp):
        '''Grab a QSO from the SDSS database
        Similar to sdss_objinf in XIDL

        Parameters:
        ----------
        inp: tuple or str
          tuple: (PLATE,FIBER)
          string (JXXXXXX.X+XXXXXX.X format)

        Returns:
        ----------
        row: QTable
          Single line from SDSS DR7 table
        '''
        # Branch on inp
        if isinstance(inp,tuple):
            mt = np.where( (self._data['PLATE']==inp[0]) & 
                (self._data['FIBERID']==inp[1]))[0]
        elif isinstance(inp,basestring):
            # Get RA/DEC
            radec = xor.stod1(inp)
            # Match
            mt = np.where( (np.abs(self._data['RAOBJ']-radec[0].value)<5e-3) & 
                (np.abs(self._data['DECOBJ']-radec[1].value)<5e-3))[0]
        else:
            raise ValueError('SDSS_QUASAR: Bad input type')

        # Parse and return
        if len(mt) == 0:
            if isinstance(inp,tuple):
                print('Plate={:d}, FIBERID={:d} not found in SDSS-{:s}'.format(
                    inp[0],inp[1],self._version))
            elif isinstance(inp,basestring):
                print('Quasar {:s} not found in SDSS-{:s}'.format(inp,self._version))
            return
        elif len(mt) == 1:
            return self._data[mt]
        else:
            raise ValueError('Not expecting this')

    def __repr__(self):
        ''' For printing
        '''
        # Generate 
        return '[SDSS Quasars: Version={:s}]'.format(self._version)

