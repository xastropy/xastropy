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
import imp
from abc import ABCMeta, abstractmethod

from astropy.io import ascii 
from astropy import units as u

from xastropy.igm.abs_sys.ionic_clm import Ions_Clm, Ionic_Clm_File
from xastropy.xutils import xdebug as xdb
from xastropy.xutils import arrays as xarray

xa_path = imp.find_module('xastropy')[1]

###################### ######################
###################### ######################
# Class for Absorption Line Survey
class AbslineSurvey(object):
    """A survey of absorption line systems. Each system may be a
    collection of Absline_System's

    Attributes:
        nsys: An integer representing the number of absorption systems
        abs_type: Type of Absorption system (DLA, LLS)
        ref: Reference to the Survey
    """

    __metaclass__ = ABCMeta

    # Init
    def __init__(self, abs_type, flist=None, tree='', ref=''):
        # Expecting a list of files describing the absorption systems
        """  Initiator

        Parameters
        ----------
        flist : string
          ASCII file giving a list of systems (usually .dat files)
        abs_type : string
          Type of Abs Line System, e.g.  MgII, DLA, LLS
        mask : Boolean array
          Mask for the systems
        ref : string
          Reference(s) for the survey
        """
        self.flist = flist
        self.abs_type = abs_type
        self.ref = ref
        self._abs_sys = []

        # Load up (if possible)
        if flist is not None:
            data = ascii.read(tree+flist, data_start=0, guess=False,format='no_header')

            self.dat_files = list(data['col1'])
            self.nsys = len(self.dat_files)
            self.tree = tree
            print('Read {:d} files from {:s} in the tree {:s}'.format(
                self.nsys, self.flist, self.tree))
            # Generate AbsSys list
            for dat_file in self.dat_files:
                if abs_type == 'LLS':
                    from xastropy.igm.abs_sys.lls_utils import LLSSystem
                    self._abs_sys.append(LLSSystem(dat_file=dat_file,tree=tree))
                elif abs_type == 'DLA':
                    from xastropy.igm.abs_sys.dla_utils import DLA_System
                    self._abs_sys.append(DLA_System(dat_file=dat_file,tree=tree))
                else: # Generic
                    from xastropy.igm.abs_sys.abssys_utils import GenericAbsSystem
                    self._abs_sys.append(GenericAbsSystem(abs_type,dat_file=tree+dat_file))
            # Mask
            self.mask = (np.ones(self.nsys) == np.ones(self.nsys))
        else:
            self.nsys = 0 
            self.mask = None

    # Get abs_sys (with mask)
    def abs_sys(self):
        lst = self._abs_sys
        return xarray.lst_to_array(lst,mask=self.mask)

    # Get attributes
    def __getattr__(self, k):
        try:
            lst = [getattr(abs_sys,k) for abs_sys in self._abs_sys]
        except ValueError:
            raise ValueError

        return xarray.lst_to_array(lst,mask=self.mask)


    # Get ions
    def fill_ions(self): # This may be overloaded!
        '''
        Loop on systems to fill in ions
        '''
        for abs_sys in self._abs_sys:
            abs_sys.get_ions()

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
        keys = (u'name',) + self.abs_sys()[0].ions.keys
        key_dtype= ('<U32',) + self.abs_sys()[0].ions.key_dtype
        t = Table(names=keys, dtype=key_dtype)

        # Loop on systems (Masked)
        for abs_sys in self.abs_sys():
            ## Mask?
            #if not self.mask[self.abs_sys.index(abs_sys)]: 
            #    continue
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
            else:
                if skip_null is False:
                    row = [abs_sys.name] + [0 for key in keys[1:]]
                    t.add_row( row )   
        # Return
        return t

    # Mask
    def upd_mask(self, msk, increment=False):
        '''
        Update the Mask for the abs_sys

        Parameters
        ----------
        msk : array (usually Boolean)
           Mask of systems 
        increment : Boolean (False)
           Increment the mask (i.e. keep False as False)
        '''
        if len(msk) == len(self._abs_sys):  # Boolean mask
            if increment is False:
                self.mask = msk
            else:
                self.mask = (self.mask == True) & (msk == True)
        else:
            raise ValueError('abs_survey: Needs developing!')

    # Printing
    def __repr__(self):
        if self.flist is not None:
            return '[AbslineSurvey: {:s} {:s}, nsys={:d}, type={:s}, ref={:s}]'.format(
                self.tree, self.flist, self.nsys, self.abs_type, self.ref)
        else:
            return '[AbslineSurvey: nsys={:d}, type={:s}, ref={:s}]'.format(
                self.nsys, self.abs_type, self.ref)

# Class for Generic Absorption Line System
class GenericAbsSurvey(AbslineSurvey):
    """A simple absorption line survey
    """
    def __init__(self, **kwargs):
        AbslineSurvey.__init__(self,'Generic',**kwargs)


###################### ###################### ######################
###################### ###################### ######################
###################### ###################### ######################
# Testing
###################### ###################### ######################
if __name__ == '__main__':

    # Test the Survey
    tmp = AbslineSurvey('Lists/lls_metals.lst',abs_type='LLS',
                         tree='/Users/xavier/LLS/')
    print(tmp)
    print('z  NHI')
    xdb.xpcol(tmp.zabs, tmp.NHI)
    
    #xdb.set_trace()
    
    print('abssys_utils: All done testing..')
        
