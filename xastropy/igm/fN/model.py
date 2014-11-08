"""
#;+ 
#; NAME:
#; fN.model
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for calculating fN models
#;   07-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""

from __future__ import print_function

import numpy as np
import os
from xastropy.xutils import xdebug as xdb
from xastropy.spec import abs_line, voigt

#from astropy.io import fits, ascii

class fN_Model(object):
    """A Class for fN models

    Attributes:
       fN_mtype: string
          Model type for the fN
           'plaw' -- Power Laws
           'Hspline' -- Hermite monotonic spline 
       zmnx: tuple
          Redshift range where this model applies (zmin,zmax)
       npivot: int 
          Number f pivots
    """

    # Initialize with type
    def __init__(self, fN_mtype, zmnx=(0.,0.), npivot=0):
        self.fN_mtype = fN_mtype  # Should probably check the choice
        self.zmnx = zmnx  
        self.npivot = npivot  

    # Output
    def __repr__(self):
        return ('[%s: %s_%s zmnx=(%g,%g)]' %
                (self.__class__.__name__,
                 self.fN_mtype, self.zmnx[0], self.zmnx[1] ) )


## #################################    
## #################################    
## TESTING
## #################################    
if __name__ == '__main__':

    # Reproduce the main figure from P14

    # MCMC Analysis
    chain_file = os.environ.get('DROPBOX_DIR')+'IGM/fN/MCMC/mcmc_spline_k13r13o13n12_8.fits'
    #strct = x_mcmc_chain_stats('../Analysis/'+chain_file)

    # Data
