"""
#;+ 
#; NAME:
#; lls_literature
#;   Ordered by publication date
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for loading up literature data on Lyman Limit Systems
#;   29-Jun-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import os, copy, sys, imp, glob
import numpy as np
import yaml, time

from astropy import units as u
from astropy.io import ascii 

from xastropy.igm.abs_sys.lls_utils import LLSSystem
from xastropy.igm.abs_sys.ionclms import IonClms
from xastropy.obs import radec as xor 
from xastropy.xutils import xdebug as xdb

xa_path = imp.find_module('xastropy')[1]

#class LLSSystem(AbslineSystem):
#class LLS_Survey(Absline_Survey):

def zonak2004():
    '''Zoank, S. et al. 2004, ApJ, 2004, 606, 196
    PG1634+706
    HST+Keck spectra
    MgII, SiIV, SiIII from Table 2.  Summing Subsystems A (Model 2) and B
       Errors estimated by JXP (not reported)
       SiIII in A may be a model
       SiIV in B may be a model
    Total NHI from LL. Taken from Fig 3 caption.  
       Error estimated by JXP 
    Not all EWs in Table 1 included
    Adopting their M/H
    '''
    # Setup
    radec = xor.stod1('J163428.9897+703132.422') # SIMBAD
    lls = LLSSystem(name='PG1634+706_z1.041', RA=radec[0], Dec=radec[1], zem=1.337,
        zabs=1.0414, vlim=[-200., 30.]*u.km/u.s, NHI=17.23, MH=-1.4,
        sigNHI=np.array([0.15,0.15])) 
    # SubSystems
    lls.mk_subsys(2) 
    # Abundances
    adict = dict(MgII={'clm': log_sum([11.45,11.90,12.02,11.68]), 'sig_clm': 0.05, 'flg_clm': 1},
        SiIII={'clm': log_sum([12.5,12.5,12.8,12.7]), 'sig_clm': 0.25, 'flg_clm': 1},
        SiIV={'clm': log_sum([10.9,10.8,11.2,11.1]), 'sig_clm': 0.15, 'flg_clm': 1} )
    lls.subsys['A']._ionclms = IonClms(idict=adict)
    bdict = dict(SiIII={'clm': log_sum([11.8,12.8,12.4]), 'sig_clm': 0.15, 'flg_clm': 1},
        SiIV={'clm': log_sum([11.2,12.2,11.8]), 'sig_clm': 0.15, 'flg_clm': 1} )
    lls.subsys['B']._ionclms = IonClms(idict=bdict)
    # Total
    lls._ionclms = lls.subsys['A']._ionclms.sum(lls.subsys['B']._ionclms)
    lls.Refs.append('Zon04')
    # Return
    return lls

#####
def log_sum(logN):
    '''Sum up logN values return the log
    '''
    Nsum = np.sum(10.**np.array(logN))
    return np.log10(Nsum)

######
if __name__ == '__main__':
     
    flg_test = 0
    flg_test += 2**0  # Zonak2004

    # Test ions
    if (flg_test % 2**1) >= 2**0:
        lls = zonak2004()
        print(lls)
        #xdb.set_trace()
    
    # Plot the LLS
