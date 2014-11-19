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
from abc import ABCMeta, abstractmethod

from astropy.io import ascii 
from astropy import units as u
from astropy.coordinates import SkyCoord

from xastropy.igm.abs_sys.ionic_clm import Ions_Clm, Ionic_Clm_File
from xastropy.xutils import xdebug as xdb

###################### ######################
###################### ######################
###################### ######################
# Class for Absorption Line System
class Absline_System(object):
    """An absorption line system

    Attributes:
        name: Coordinates
        coord: Coordinates
        epoch: Epoch (e.g. 2000.0)
        zabs : float
          Absorption redshift
        NHI:  float
          Log10 of the HI column density
        sigNHI:  np.array(2)
          Log10 error of the HI column density (-/+)
        ions:  Ions_Clm Class
    """

    __metaclass__ = ABCMeta

    # Init
    def __init__(self, abs_type, zabs=0., NHI=0., epoch=2000., dat_file=None, tree=None):
        """  Initiator

        Parameters
        ----------
        abs_type : string
          Type of Abs Line System, e.g.  MgII, DLA, LLS
        dat_file : string
          ASCII .dat file summarizing the system
        """
        self.zabs = zabs
        self.NHI = NHI
        self.epoch = epoch
        # Abs type
        if abs_type == None:
            self.abs_type = 'NONE'
        else: self.abs_type = abs_type
        # Tree
        if tree == None: tree = ''
        self.tree = tree
        # Fill in
        if dat_file != None:
            print('absys_utils: Reading %s file' % dat_file)
            self.parse_dat_file(dat_file)

    # Read a .dat file
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

        # Name
        self.name = ('J'+
                    self.coord.ra.to_string(unit=u.hour,sep='',pad=True)+
                    self.coord.dec.to_string(sep='',pad=True,alwayssign=True))

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
            except: key_sigNHI='0.0 0.0'
        #pdb.set_trace()
        self.sigNHI = np.array(map(float,key_sigNHI.split()))
        #print('sigNHI: ', self.sigNHI)

        # Finish
        if verbose: print(datdic)
        if flg_out != None:
            if (flg_out % 2) == 1: ret_val = [datdic]
            else: ret_val = [0]
            return ret_val

    # Parse the ion files
    def get_ions(self,clm_fil):
        # Read .clm file
        self.clm_analy = Ionic_Clm_File(clm_fil)

        # Read .all file
        #xdb.set_trace()
        all_fil = self.tree+self.clm_analy.ion_fil.split('.ion')[0]+'.all'
        self.ions = Ions_Clm(all_fil)
            
    @abstractmethod
    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        pass

    # #############
    def __repr__(self):
        return ('[Absline_System: %s %s %s %s, %g, NHI=%g]' %
                (self.name, self.abs_type,
                 self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True),
                 self.zabs, self.NHI))


# Class for Generic Absorption Line System
class Generic_System(Absline_System):
    """A simple absorption system

    """
    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'Generic'

# Class for Generic Absorption Line System
class Abs_Sub_System(Absline_System):
    """A simple absorption system

    """
    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'SubSystem'






    
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
        
