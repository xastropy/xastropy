"""
#;+ 
#; NAME:
#; lines_utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for SpectralLine class
#;   03-Mar-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import barak
import xastropy
import numpy as np
import matplotlib.pyplot as plt
import pdb
from astropy import constants as const
from astropy import units as u

import xastropy.atomic as xatom

# Class for Spectral line
class SpectralLine(object):
    """Class for a spectral line.  Emission or absorption 

    Attributes:
        ltype: string
          type of line, e.g.  Abs, Emiss
        wrest: float
          Rest wavelength of the spectral feature
    """
    __metaclass__ = ABCMeta

    # Initialize with wavelength
    def __init__(self, ltype, wrest):
        """  Initiator

        Parameters
        ----------
        ltype : string
          Type of Spectral line, (Abs
        wrest : float
          Rest wavelength
        """

        # Required
        self.ltype = ltype
        if ltype not in ['Abs']:
            raise ValueError('spec/lines: Not ready for type {:s}'.format(ltype))

        if type(wrest) is not float:
            raise ValueError('spec/lines: Rest wavelength must be a float!')
        self.wrest = wrest

        # Other
        self.atomic = {} # Atomic Data
        self.analy = {} # Analysis inputs (e.g. spectrum; from .clm file or AbsID)
        self.measure = {} # Measured quantities (e.g. column, EW, centroid)

        # Fill atomic data
        self.fill_data()

    # Output
    def __repr__(self):
        return ('[{:s}: wrest={:g}]'.format(
                self.__class__.__name__, self.wrest))

# Class for Generic Absorption Line System
class AbsLine(SpectralLine):
    """Spectral absorption line
    """
    def print_specline_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'AbsLine'

    # Fill Analy
    def fill_data(self):
        import xastropy.spec.abs_line as xspa

        # Data
        self.atomic = xspa.abs_line_data(self.wrest)

        #
        self.analy['VLIM'] = [0., 0.]*u.km/u.s # Velocity limit of line
        self.analy['FLG_ANLY'] = 1 # Analyze
        self.analy['FLG_EYE'] = 0
        self.analy['FLG_LIMIT'] = 0 # No limit
        self.analy['DATFIL'] = '' 
        self.analy['IONNM'] = self.atomic['name']

