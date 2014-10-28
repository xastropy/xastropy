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

from __future__ import print_function

import numpy as np
import pdb
from astropy.io import ascii 
from astropy import units as u
from astropy.coordinates import SkyCoord

# Class for Absorption Line Survey
class Absline_Survey(object):
    """A survey of absorption line systems

    Attributes:
        nsys: An integer representing the number of absorption systems
        abs_type: Type of Absorption system (DLA, LLS)
        ref: Reference to the Survey
    """

    # Init
    def __init__(self, flist, abs_type=None, tree='', ref=None):
        # Expecting a list of files describing the absorption systems
        data = ascii.read(tree+flist, data_start=0, guess=False,format='no_header')

        self.flist = flist
        self.dat_files = list(data['col1'])
        self.nsys = len(self.dat_files)
        self.tree = tree
        print('Read %d files from %s in the tree %s' % (self.nsys, self.flist, self.tree))

        # Generate AbsSys list
        self.abs_sys = []
        for dat_file in self.dat_files:
            self.abs_sys.append(Absline_System(tree+dat_file))

        # Other
        self.abs_type = abs_type
        self.ref = ref
        #

    # Printing
    def __repr__(self):
        return '[Absline_Survey: %s %s, %d, %s, %s]' % (self.tree, self.flist,
                                                        self.nsys, self.abs_type, self.ref)

###################### ######################
###################### ######################
###################### ######################
# Class for Absorption Line System
class Absline_System(object):
    """An absorption line system

    Attributes:
        coord: Coordinates
        dec: Declination
        epoch: Epoch (e.g. 2000.0)
    """

    # Init
    def __init__(self, zabs=0., NHI=0., epoch=2000.):
        self.zabs = zabs
        self.NHI = NHI
        self.epoch = epoch

    def parse_dat_file(self,dat_file,verbose=False):
        # Define
        datdic = {}
        # Open
        f=open(dat_file,'r')
        for line in f:
            tmp=line.split(' ! ')
            tkey=tmp[1].strip()
            key=tkey.replace(' ','')
            val=tmp[0].strip()
            datdic[key]=val
        f.close()
        #pdb.set_trace()

        # Pull attributes
        try:  # RA/DEC
            self.coord = SkyCoord(datdic['RA(2000)'], datdic['DEC(2000)'],
                         'icrs', unit=(u.hour, u.deg))
            #print(datdic['RA(2000)'], datdic['DEC(2000)'])
            #pdb.set_trace()
            if self.epoch != 2000.:
                assert False, 'Not ready for non-2000'
        except:
            print('blah exception')
        # Return
        if verbose: print(datdic)

    def __repr__(self, zabs=0., NHI=0., epoch=2000.):
        return ('[Absline_System: %s %s, z=%g, NHI=%g]' %
                (self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True),
                 self.zabs, self.NHI)

######################
# Testing
if __name__ == '__main__':
    # Test Absorption System
    tmp1 = Absline_System()
    tmp1.parse_dat_file('/Users/xavier/LLS/Data/UM669.z2927.dat')
    print(tmp1)

    # Test the Survey
    tmp = Absline_Survey('Lists/lls_metals.lst',abs_type='LLS',
                         tree='/Users/xavier/LLS/')
    print(tmp)
    
        
