"""
#;+ 
#; NAME:
#; abssys_utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for Absorption Systems
#;   23-Oct-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

from astropy.io import ascii 
from astropy import units as u
from astropy.coordinates import SkyCoord

from xastropy.igm.abs_sys.ionic_clm import Ions_Clm, Ionic_Clm_File
from xastropy.xutils import xdebug as xdb

###################### ######################
###################### ######################
# Class for Absorption Line Survey
class Absline_Survey(object):
    """A survey of absorption line systems. Each system may be a
    collection of Absline_System's

    Attributes:
        nsys: An integer representing the number of absorption systems
        abs_type: Type of Absorption system (DLA, LLS)
        ref: Reference to the Survey
    """

    # Init
    def __init__(self, flist, abs_type=None, tree='', ref=None):
        # Expecting a list of files describing the absorption systems
        """  Initiator

        Parameters
        ----------
        flist : string
          ASCII file giving a list of systems (usually .dat files)
        abs_type : string
          Type of Abs Line System, e.g.  MgII, DLA, LLS
        ref : string
          Reference(s) for the survey
        """
        data = ascii.read(tree+flist, data_start=0, guess=False,format='no_header')

        self.flist = flist
        self.dat_files = list(data['col1'])
        self.nsys = len(self.dat_files)
        self.tree = tree
        print('Read %d files from %s in the tree %s' % (self.nsys, self.flist, self.tree))

        # Generate AbsSys list
        self.abs_sys = []
        for dat_file in self.dat_files:
            if abs_type == 'LLS':
                from xastropy.igm.abs_sys.lls_utils import LLS_System
                #xdb.set_trace()
                self.abs_sys.append(LLS_System(dat_file=dat_file,tree=tree))
            else: # Generic
                from xastropy.igm.abs_sys.abssys_utils import Generic_System
                self.abs_sys.append(Generic_System(abs_type,dat_file=tree+dat_file))

        # Other
        self.abs_type = abs_type
        self.ref = ref
        #

    # Get attributes
    def __getattr__(self, k):
        return np.array( [getattr(abs_sys,k) for abs_sys in self.abs_sys] )

    # Get ions
    def fill_ions(self):
        '''
        Loop on systems to fill in ions
        '''
        for abs_sys in self.abs_sys:
            abs_sys.fill_ions()  # This may be overloaded!

    # Get ions
    def ions(self,iZion, skip_null=False):
        '''
        Generate a Table of columns and so on
        Restrict to those systems where flg_clm > 0

        Parameters
        ----------
        iZion : tuple
           Z, ion   e.g. (6,4) for CIV
        skip_null : boolean (False)
           Skip systems without an entry, else pad with zeros 

        Returns
        ----------
        Table of values for the Survey
        '''
        from astropy.table import Table
        keys = (u'name',) + self.abs_sys[0].ions.keys
        key_dtype= ('<U32',) + self.abs_sys[0].ions.key_dtype
        t = Table(names=keys, dtype=key_dtype)

        # Loop on systems
        for abs_sys in self.abs_sys:
            # Grab
            try:
                idict = abs_sys.ions[iZion]
            except KeyError:
                if skip_null is False:
                    row = [abs_sys.name] + [0 for key in keys[1:]]
                    t.add_row( row )   
                continue
            # Cut on flg_clm
            if idict['flg_clm'] > 0:
                row = [abs_sys.name] + [idict[key] for key in keys[1:]]
                t.add_row( row )   # This could be slow
        # Return
        return t

    # Printing
    def __repr__(self):
        return '[Absline_Survey: {s} {s}, nsys={d}, type={s}, ref={s}]'.format(self.tree, self.flist,
                                                            self.nsys, self.abs_type, self.ref)







    
###################### ###################### ######################
###################### ###################### ######################
###################### ###################### ######################
# Testing
###################### ###################### ######################
if __name__ == '__main__':

    # Test Absorption System
    tmp1 = Absline_System('LLS')
    tmp1.parse_dat_file('/Users/xavier/LLS/Data/UM669.z2927.dat')
    print(tmp1)

    #pdb.set_trace()

    # Test the Survey
    tmp = Absline_Survey('Lists/lls_metals.lst',abs_type='LLS',
                         tree='/Users/xavier/LLS/')
    print(tmp)
    print('z  NHI')
    xdb.xpcol(tmp.zabs, tmp.NHI)
    
    #xdb.set_trace()
    
    print('abssys_utils: All done testing..')
        
