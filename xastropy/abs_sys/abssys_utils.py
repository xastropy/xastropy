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
    def __init__(self, abs_type, zabs=0., NHI=0., epoch=2000.):
        self.zabs = zabs
        self.NHI = NHI
        self.epoch = epoch
        self.abs_type = abs_type

    def parse_dat_file(self,dat_file,verbose=False,flg_out=None):
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

        #  #########
        # Pull attributes

        # RA/DEC
        try:
            ras,decs = (datdic['RA(2000)'], datdic['DEC(2000)'])
            #print(datdic['RA(2000)'], datdic['DEC(2000)'])
            #pdb.set_trace()
        except:
            ras, decs = ('00 00 00', '+00 00 00')
        self.coord = SkyCoord(ras, decs, 'icrs', unit=(u.hour, u.deg))

        # zabs
        try: 
            self.zabs = float(datdic['zabs'])
        except: self.zabs=0.

        # NHI
        try: 
            self.NHI = float(datdic['NHI']) # DLA format
        except:
            try:
                self.NHI = float(datdic['NHItot']) # LLS format
            except: self.NHI=0.

        # NHIsig
        try: 
            key_sigNHI = datdic['sig(NHI)'] # DLA format
        except:
            try:
                key_sigNHI = datdic['NHIsig'] # LLS format
            except: key_sigNHI='0.0'
        #pdb.set_trace()
        self.sigNHI = np.array(map(float,key_sigNHI.split()))
        #print('sigNHI: ', self.sigNHI)

        # Finish
        if verbose: print(datdic)
        if flg_out != None:
            if (flg_out % 2) == 1: ret_val = [datdic]
            else: ret_val = [0]
            return ret_val
        
                        

    def __repr__(self):
        return ('[Absline_System: %s %s %s, %g, NHI=%g]' %
                (self.abs_type,
                 self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True),
                 self.zabs, self.NHI))

######################
# Testing
if __name__ == '__main__':

    # Test Absorption System
    tmp1 = Absline_System('LLS')
    tmp1.parse_dat_file('/Users/xavier/LLS/Data/UM669.z2927.dat')
    print(tmp1)

    # Test the Survey
    tmp = Absline_Survey('Lists/lls_metals.lst',abs_type='LLS',
                         tree='/Users/xavier/LLS/')
    print(tmp)
    
        
