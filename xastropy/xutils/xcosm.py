"""
#;+ 
#; NAME:
#; xcosm
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for cosmology things
#;    Mainly a wrapper on astropy.cosmology
#;   10-Dec-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import gzip
import subprocess, glob

from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

from xastropy.xutils import xdebug as xdb

# def bintab_to_table(fits_fil,exten=1, silent=False):
# def table_to_fits(table, outfil, compress=False, comment=None):

#
def quick_cosm(H0=70., Om0=0.3):
    """Return a FlatLCDM cosmology from astropy
    Parameters
    ----------
    H0 : float
        Hubble's parameter in km/s/Mpc
    0m0 : float
        Matter density parameter
    """
    print('LCDM with H0={:g}, Om0={:g}'.format(H0,Om0))
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    #
    return cosmo

