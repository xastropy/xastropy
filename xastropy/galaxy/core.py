"""
#;+ 
#; NAME:
#; galaxy.core
#;    Version 1.0
#;
#; PURPOSE:
#;    Core routines for galaxy analysis
#;   29-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import os, copy, sys

import numpy as np

from astropy import units as u
from astropy.io import ascii 
from astropy.coordinates import SkyCoord

from xastropy.obs import radec as xra

from xastropy.xutils import xdebug as xdb

# Class for LLS Absorption Lines 
class Galaxy(object):
    """A Galaxy Class

    Inputs:
    ra, dec

    Attributes:
      name: string
        Name(s)
      z: float
        Adopted redshift
      coord: Coordinates
      mstar: float
        Stellar mass (MsolMass)
      
    """
    # Initialize with a .dat file
    def __init__(self, ra, dec, z=None):

        self.z = z
        
        # Coord
        self.coord = xra.to_coord( (ra,dec) ) 

        # Name
        self.name = ('J'+
                    self.coord.ra.to_string(unit=u.hour,sep='',pad=True)+
                    self.coord.dec.to_string(sep='',pad=True,alwayssign=True))

    # #############
    def __repr__(self):
        return ('[Galaxy: {:s} {:s} {:s}, z={:g}]'.format(
                 self.name, 
                 self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True), 
                 self.z) )













## #################################    
## #################################    
## TESTING
## #################################    
if __name__ == '__main__':

    # Instantiate
    gal = Galaxy()
    print(gal)
