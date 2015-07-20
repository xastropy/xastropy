"""
#;+
#; NAME:
#; utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for spectral utilities
#;       Primarily overloads of Spectrum1D
#;   07-Sep-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
import astropy as apy

from astropy import units as u
from astropy import constants as const
from astropy.io import fits

#from specutils import Spectrum1D
#from specutils.wcs import BaseSpectrum1DWCS, Spectrum1DLookupWCS
#from specutils.wcs.specwcs import Spectrum1DPolynomialWCS

from xastropy.xutils import xdebug as xdb

'''
DEPRECATED
'''