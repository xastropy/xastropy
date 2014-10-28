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
from xastropy.abs_sys.abssys_utils import Absline_System
from astropy import units as u

# Class for LLS Absorption Lines 
class LLS_System(Absline_System):
    """An LLS absorption system

    Attributes:
        tau_ll: Opacity at the Lyman limit
    """

    # Initialize with a .dat file
    def __init__(self, dat_file=None):
        # Generate with type
        Absline_System.__init__(self,'LLS')
        #
        if dat_file != None:
            self.parse_dat_file(dat_file)
        # Set tau_LL
        self.tau_LL = (10.**self.NHI)*6.3391597e-18 # Should replace with photocross

    # Modify standard dat parsing
    def parse_dat_file(self,dat_file):
        # Standard Call
        out_list = Absline_System.parse_dat_file(self,dat_file,flg_out=1)
        datdic = out_list[0]
        # LLS keys
        self.MH = float(datdic['[M/H]ave'])

    # Output
    def __repr__(self):
        return ('[LLS_System: %s %s, %g, NHI=%g, M/H=%g]' %
                (self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True),
                 self.zabs, self.NHI, self.MH))

if __name__ == '__main__':
    # Test Absorption System
    tmp1 = LLS_System(dat_file='/Users/xavier/LLS/Data/UM669.z2927.dat')
    print(tmp1)
