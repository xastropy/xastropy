"""
#;+ 
#; NAME:
#; lls_utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for Lyman Limit Systems
#;   27-Oct-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""

from __future__ import print_function

import numpy as np
import pdb
from astropy.io import ascii 
from xastropy.abs_sys.abssys_utils import Absline_System
from xastropy.abs_sys.ionic_clm import Ionic_Clm
from xastropy.abs_sys.ionic_clm import Ionic_Clm_File
from astropy import units as u

# Class for LLS Absorption Lines 
class LLS_System(Absline_System):
    """An LLS absorption system

    Attributes:
        tau_ll: Opacity at the Lyman limit
    """

    # Initialize with a .dat file
    def __init__(self, dat_file=None, tree=None):
        # Generate with type
        Absline_System.__init__(self,'LLS')
        # Over-ride tree?
        if tree != None: self.tree = tree
        # Parse .dat file
        if dat_file != None:
            self.parse_dat_file(tree+dat_file)

        # Set tau_LL
        self.tau_LL = (10.**self.NHI)*6.3391597e-18 # Should replace with photocross
        self.ionic = {}

        # Name

    # Modify standard dat parsing
    def parse_dat_file(self,dat_file):
        # Standard Call
        out_list = Absline_System.parse_dat_file(self,dat_file,flg_out=1)
        datdic = out_list[0]

        # LLS keys
        self.MH = float(datdic['[M/H]ave'])
        self.nsub = int(datdic['Nsubsys'])
        self.cldyfil = datdic['CloudyGridFile']

        # LLS Subsystems
        if self.nsub > 0:
            subsys = {}
            lbls= map(chr, range(65, 91))
            keys = (['zabs','NHI','NHIsig','NH','NHsig','logx','sigx','b','bsig','Abundfile',
                     'U','Usig','flg_low','flg_alpha','[alpha/H]','sig[a/H]',
                     'flg_Fe','[Fe/H]','sig[Fe/H]','VPFITfile'])
            for i in range(self.nsub):
                # Generate
                subsys[lbls[i]] = self.subsys()
                # Fill in
                for key in list(subsys[lbls[i]].keys()):
                    try:
                        tmpc = datdic[lbls[i]+key]
                    except:
                        print('lls_utils: Key %s not found', lbls[i]+key)
                    else:  # Convert
                        val = subsys[lbls[i]][key]
                        #pdb.set_trace()
                        if val.__class__ == np.ndarray:  
                            subsys[lbls[i]][key] = np.array(map(float,tmpc.split()))
                        else: # Single value
                            subsys[lbls[i]][key] = (map(type(val),[tmpc]))[0]
            # Encode
            self.subsys = subsys

    # Generate the Ionic dictionary
    def get_ions(self):
        lbls= map(chr, range(65, 91))
        # Loop on Sub-Systems
        for kk in range(self.nsub):
            # Read .clm file
            tmp = Ionic_Clm_File(self.tree+self.subsys[lbls[kk]]['Abundfile'])
            # Fill it up
            print('lls_utils.get_ions: The next line needs to be changed!')
            Dumb_Class = type('Dummy_Object', (object,), {})
            self.subsys[lbls[kk]]['Ionic'] = Dumb_Class()
            self.subsys[lbls[kk]]['Ionic'].analy = tmp # THIS NEEDS TO BE CHANGED
            
            #pdb.set_trace()
            #self.ionic[lbls[kk],
        #pdb.set_trace()

    # Subsystem Dict
    def subsys(self):
        keys = (['zabs','NHI','NHIsig','NH','NHsig','logx','sigx','b','bsig','Abundfile',
                'U','Usig','flg_low','flg_alpha','[alpha/H]','sig[a/H]',
                'flg_Fe','[Fe/H]','sig[Fe/H]','VPFITfile'])
        values = ([0., 0., np.zeros(2), 0., np.zeros(2), 0., np.zeros(2), 0., 0.,
                   '', 0., np.zeros(2), 0, 0, 0., 0., 0, 0., 0., ''])
        return dict(zip(keys,values))

    # Output
    def __repr__(self):
        return ('[LLS_System: %s %s, %g, NHI=%g, M/H=%g]' %
                (self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True),
                 self.zabs, self.NHI, self.MH))

if __name__ == '__main__':
    # Test Absorption System
    tmp1 = LLS_System(dat_file='Data/UM669.z2927.dat',
                      tree='/Users/xavier/LLS/')
    print(tmp1)
    print(tmp1.subsys)
    tmp1.get_ions()
    #tmp1.ionic['1215.6701'] = Ionic_Clm(1215.6701)
    #print(tmp1.ionic)
    #print(tmp1.ionic['1215.6701'].wave)
    
