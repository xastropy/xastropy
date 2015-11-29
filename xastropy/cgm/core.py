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


# Class for CGM Survey
class CGMAbsSurvey(object):
    """A CGM Survey class in absorption

    Attributes:
    """
    # Initialize with a .dat file
    def __init__(self, tree=None, survey=''):

        from xastropy.igm.abs_sys.abs_survey import AbslineSurvey

        # Name of survey
        self.survey = ''
        self.ref = ''
        self.nsys = 0

        # Generate with type
        #Absline_Survey.__init__(self, '', abs_type='CGM', tree=tree)

        self.cgm_abs = []
        self.mask = None

    # Extend attributes
    def __getattr__(self, k):
        # Try Self first
        try:
            lst = [getattr(cgm_abs,k) for cgm_abs in self.cgm_abs]
        except AttributeError:
            # Try AbsLine_Sys next
            try:
                lst = [getattr(cgm_abs.abs_sys,k) for cgm_abs in self.cgm_abs] 
            except AttributeError:
                # Galaxy?
                try:
                    lst = [getattr(cgm_abs.galaxy,k) for cgm_abs in self.cgm_abs] 
                except AttributeError:
                    print('cgm.core: Attribute not found!')
                    xdb.set_trace()
        # Return array
        return xu_array.lst_to_array(lst,mask=self.mask)

    # Kinematics (for convenience)
    def abs_kin(self, lbl):
        """  Create a Table of the Kinematic info

        Parameters
        ----------
        lbl : string
          Label for the Kinematics dict
        TODO: 
          Add wrest!!
        """
        from astropy.table import Table

        keys = self.cgm_abs[0].abs_sys.kin[lbl].keys
        t = Table(names=keys,
                  dtype=self.cgm_abs[0].abs_sys.kin[lbl].key_dtype)

        for cgm_abs in self.cgm_abs:
            try:
                kdict = cgm_abs.abs_sys.kin[lbl]
            except KeyError:
                # No dict.  Filling in zeros
                row =  [0 for key in keys]
                t.add_row( row )   
                continue
            # Filling
            row = [kdict[key] for key in keys]
            t.add_row( row )   
        return t

    # Printing
    def __repr__(self):
        str1 = '[CGM_Survey: {:s} nsys={:d}, ref={:s}]\n'.format(self.survey, self.nsys, self.ref) 
        for ii in range(self.nsys):
            str1 = str1+self.cgm_abs[ii].abs_sys.__repr__()+'\n'
        return str1


# ###################### #######################
# Testing
if __name__ == '__main__':

    # Initialize
    tmp = CGMAbs()
    print(tmp)

    tmp2 = CGMAbsSurvey()
    print(tmp2)
