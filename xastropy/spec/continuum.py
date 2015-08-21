"""
#;+
#; NAME:
#; continuum 
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for continuum code 
#;   20-Aug-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp
import astropy as apy

from astropy import units as u
from astropy import constants as const
from astropy.io import fits, ascii

from linetools.spectra.xspectrum1d import XSpectrum1D

from xastropy.xutils import xdebug as xdb

xa_path = imp.find_module('xastropy')[1]

def init_conti_dict(Norm=0., tilt=0., piv_wv=0.):
    '''Initialize a continuum conti_dict
    Parameters:
    ----------
    Norm: float, optional
      Normaliztion
    tilt: float, optional
      Power-law tilt to continuum
    piv_wv: float, optional
      Pivot wave for tilt.  Best kept *without* units

    Returns:
    ---------
    conti_dict: dict 
      Useful for simple modeling
    '''
    conti_dict = {'Norm': 0., 'tilt': 0., 'piv_wv': np.median(spec.dispersion.value)}
    #
    return conti_dict

def get_telfer_spec(zqso=0.):
    '''Generate a Telfer QSO composite spectrum
    Paraemters:
    ----------
    zqso: float, optional
      Redshift of the QSO

    Returns:
    --------
    telfer_spec: XSpectrum1D
      Spectrum
    '''
    telfer = ascii.read(
        xa_path+'/data/quasar/telfer_hst_comp01_rq.ascii', comment='#')
    scale = telfer['flux'][(telfer['wrest'] == 1450.)]
    telfer_spec = XSpectrum1D.from_tuple((telfer['wrest']*(1+zqso),
        telfer['flux']/scale[0])) # Observer frame
    # Return
    return telfer_spec
 
