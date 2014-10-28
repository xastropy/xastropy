"""
#;+ 
#; NAME:
#; lls_utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for Absorption Systems
#;   27-Oct-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""

from __future__ import print_function

import numpy as np
import pdb
from astropy.io import ascii 
from xastropy.abs_sys.abssys_util import Absline_System

# Class for LLS Absorption Lines 
class LLS_System(Absline_System):
    """An LLS absorption system

    Attributes:
        tau_ll: Opacity at the Lyman limit
    """

    # Initialize with a .dat file
    def __init__(self, dat_file):
        data = ascii.read('UM669.z2927.dat', data_start=0, guess=False,format='no_header',
                          delimeter='!')
