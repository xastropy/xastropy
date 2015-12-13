"""
#;+ 
#; NAME:
#; cgm.core 
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for core routines of CGM analysis
#;   29-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
from astropy import units as u

from xastropy.igm.abs_sys.abssys_utils import AbslineSystem
from xastropy.galaxy.core import Galaxy
from xastropy.obs import radec as xra

from xastropy.xutils import xdebug as xdb
from xastropy.xutils import arrays as xu_array

# Path for xastropy
#xa_path = imp.find_module('xastropy')[1]

########################## ##########################
########################## ##########################
class CGMSys(object):
    """ Class for a CGM system
    Combines absorption lines with a Galaxy

    Inputs:
    ----------
    gal_ra: str, float, Quantity
      RA for galaxy
    gal_dec: str, float, Quantity
      DEC for galaxy
    gal_z: float
      Galaxy redshift
    bg_ra: str, float, Quantity
      RA for background source
    bg_dec: str, float, Quantity
      DEC for background source
    bg_z: float
      Redshift of background source

    Attributes
    ----------
    rho: float
      Impact parameter (u.kpc)

    JXP on 29 Nov 2014
    """

    # Initialize 
    def __init__(self, gal_ra, gal_dec, gal_z, bg_ra, bg_dec, bg_z,
        cosmo=None, verbose=False):

        # Galaxy
        self.galaxy = Galaxy(gal_ra, gal_dec, z=gal_z)

        # Name
        self.name = ('CGM'+
                    self.galaxy.coord.ra.to_string(unit=u.hour,sep='',pad=True)+
                    self.galaxy.coord.dec.to_string(sep='',pad=True,alwayssign=True))

        # Absorption system
        self.abs_sys = CGMAbs()
        self.abs_sys.coord = xra.to_coord( (bg_ra,bg_dec) ) # Background source
        self.abs_sys.zem = bg_z

        # Calcualte rho
        if cosmo is None:
            from astropy.cosmology import WMAP9 as cosmo
            if verbose is True:
                print('cgm.core: Using WMAP9 cosmology')
        ang_sep = self.abs_sys.coord.separation(self.galaxy.coord).to('arcmin')
        kpc_amin = cosmo.kpc_comoving_per_arcmin( self.galaxy.z ) # kpc per arcmin
        self.rho = ang_sep * kpc_amin / (1+self.galaxy.z) # Physical
        #xdb.set_trace()


    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'CGM'

    # Output
    def __repr__(self):
        return ('[{:s}: Galaxy RA/DEC={:s}{:s}, zgal={:g}, rho={:g}]'.format(
                self.__class__.__name__,
                 self.abs_sys.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.abs_sys.coord.dec.to_string(sep=':',pad=True,alwayssign=True),
                 self.galaxy.z, self.rho))

# Class for CGM Absorption
class CGMAbs(AbslineSystem):
    """A CGM absorption system

    Attributes:
    """
    def __init__(self): 
        # Generate with type
        AbslineSystem.__init__(self,'CGM')

        # Init
        self.ions = None

    # Output
    def __repr__(self):
        return ('[{:s}: {:s} {:s}, {:g}]'.format(
                self.__class__.__name__,
                 self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True),
                 self.zabs))

    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'CGM'




# ###################### #######################
# Testing
if __name__ == '__main__':

    # Initialize
    tmp = CGMAbs()
    print(tmp)

    tmp2 = CGMAbsSurvey()
    print(tmp2)
