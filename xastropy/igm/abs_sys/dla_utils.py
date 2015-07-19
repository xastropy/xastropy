"""
#;+ 
#; NAME:
#; dla_utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for DLAs
#;   06-Dec-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import os, copy, sys
import numpy as np
import yaml

from astropy import units as u
from astropy.io import ascii 

from xastropy.igm.abs_sys.abssys_utils import AbslineSystem, Abs_Sub_System
from xastropy.igm.abs_sys.abs_survey import AbslineSurvey
from xastropy.spec import abs_line, voigt
from xastropy.atomic import ionization as xatomi
from xastropy.xutils import xdebug as xdb

#class DLA_System(Absline_System):
#class DLA_Survey(Absline_Survey):

# Class for DLA Absorption Lines 
class DLASystem(AbslineSystem):
    """A DLA absorption system

    Attributes:
        tau_ll: Opacity at the Lyman limit
    """
    # Initialize with a .dat file
    def __init__(self, dat_file=None, tree=None):
        # Generate with type
        AbslineSystem.__init__(self,'DLA')
        # Over-ride tree?
        if tree != None:
            self.tree = tree
        else:
            self.tree = ''

        # Parse .dat file
        if dat_file != None:
            self.dat_file = self.tree+dat_file
            print('dla_utils: Reading {:s}'.format(self.dat_file))
            self.parse_dat_file(self.dat_file)
            # QSO keys
            self.qso = self.datdict['QSO name']
            self.zqso = float(self.datdict['QSO zem'])
            # Abund
            self.flg_MH = float(self.datdict['flg_mtl'])
            self.MH = float(self.datdict['[M/H]'])
            self.sigMH = float(self.datdict['sig([M/H])'])

        # Init
        self.ions = None
        self.zpeak = None

    # Output
    def __repr__(self):
        return ('[{:s}: {:s} {:s}, {:g}, NHI={:g}, M/H={:g}]'.format(
                self.__class__.__name__,
                 self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True),
                 self.zabs, self.NHI, self.MH))

    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'DLA'

# #######################################################################
# #######################################################################
# #######################################################################
# Class for DLA Survey
class DLASurvey(AbslineSurvey):
    """An DLA Survey class

    Attributes:
        
    """
    # Initialize with a .dat file
    def __init__(self, dat_file, tree=None):
        # Generate with type
        AbslineSurvey.__init__(self,dat_file,abs_type='DLA', tree=tree)


    # Default sample of DLA:  Neeleman
    @classmethod
    def default_sample(cls):
        # Local
        dla = cls('Lists/Neeleman13.lst', tree=os.environ.get('DLA'))
        dla.ref = 'Neeleman+13'

        # Return
        return dla












## #################################    
## #################################    
## TESTING
## #################################    
if __name__ == '__main__':

    flg_test = 0
    #flg_test = 1  # ions
    #
    #flg_test += 2**9 # DLA Survey NHI
    flg_test += 2**10 # DLA Survey ions

    # Test Absorption System
    print('-------------------------')
    tmp1 = DLA_System(dat_file='Data/PH957.z2309.dat',
                      tree=os.environ.get('DLA'))
    print(tmp1)

    # Test ions
    if (flg_test % 2**1) >= 2**0:
        print('-------------------------')
        clm_fil = tmp1.tree+tmp1.datdict['Abund file']
        tmp1.get_ions()
        # Print
        print('Si II: ')
        print(tmp1.ions[(14,2)])
        print(tmp1.ions.trans[15]) # CIV 1550
    


    # #############################
    # DLA Survey
    if (flg_test % 2**10) >= 2**9:
        print('-------------------------')
        dla = DLA_Survey('Lists/metal_MAR_all.lst', tree=os.environ.get('DLA'))
        dla.fill_ions()
        xdb.xhist(dla.NHI, binsz=0.10)

    # DLA Survey ions
    if (flg_test % 2**11) >= 2**10:
        dla = DLA_Survey('Lists/metal_MAR_all.lst', tree=os.environ.get('DLA'))
        dla.fill_ions()
        xdb.xhist(dla.ions((6,4),skip_null=True)['clm'], binsz=0.3,
                  xlabel=r'$\log_{10} N({\rm C}^{+3})$')
        xdb.xhist(dla.ions((14,2),skip_null=True)['clm'], binsz=0.3,
                  xlabel=r'$\log_{10} N({\rm Si}^{+})$')
