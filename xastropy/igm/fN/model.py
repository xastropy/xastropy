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
from scipy import interpolate as scii

from xastropy.xutils import xdebug as xdb
from xastropy.spec import abs_line, voigt
from xastropy.stats import mcmc

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
       pivots: array 
          log NHI values for the pivots
       zpivot: float (2.4)
          Pivot for redshift evolution
       gamma: float (1.5)
          Power law for dN/dX
    """

    # Initialize with type
    def __init__(self, fN_mtype, zmnx=(0.,0.), pivots=None,
                 param=None, zpivot=2.4):
        self.fN_mtype = fN_mtype  # Should probably check the choice
        self.zmnx = zmnx  

        # Pivots
        if pivots == None: self.pivots = np.zeros(2)
        else: self.pivots = pivots
        self.npivot = len(pivots)

        # Param
        if param == None: self.param = np.zeros(self.npivot)
        else:
            self.param = param
            #if np.amax(self.pivots) < 99.:
            #    self.pivots.append(99.)
            #    self.param = np.append(self.param,-30.)
            # Init
            if fN_mtype == 'Hspline':
                self.model = scii.PchipInterpolator(self.pivots, self.param)

        # Redshift (needs updating)
        self.zpivot = 2.4
        self.gamma = 1.5
    ##
    # Evaluate
    def eval(self, z, NHI_values):
        """ Evaluate the model at a set of NHI values

        Parameters:
        z: float
          Redshift for evaluation
        NHI_values: array
          NHI values

        Returns:
        fN: array
          Array of f(NHI) values

        JXP 07 Nov 2014
        """
        # Exception checking?

        # Evaluate
        if self.fN_mtype == 'Hspline': 
            log_fNHI = self.model.__call__(NHI_values)
        else: 
            raise ValueError('fN.model: Not ready for this model type %s' % self.fN_mtype)

        # Redshift
        log_fNHI += self.gamma * np.log10((1+z)/(1+self.zpivot))
        return log_fNHI
    ##
    # Output
    def __repr__(self):
        return ('[%s: %s zmnx=(%g,%g)]' %
                (self.__class__.__name__,
                 self.fN_mtype, self.zmnx[0], self.zmnx[1] ) )


## #################################    
## #################################    
## TESTING
## #################################    
if __name__ == '__main__':

    # Reproduce the main figure from P14
    
    from xastropy.igm.fN import data as fN_data
    from xastropy.igm.fN import model as xifm

    # MCMC Analysis
    chain_file = os.environ.get('DROPBOX_DIR')+'IGM/fN/MCMC/mcmc_spline_k13r13o13n12_8.fits.gz'
    outp = mcmc.chain_stats(chain_file)

    # Build a model
    NHI_pivots = [12., 15., 17.0, 18.0, 20.0, 21., 21.5, 22.]
    fN_model = xifm.fN_Model('Hspline', zmnx=(0.5,3.0),
                        pivots=NHI_pivots, param=outp['best_p'])
    #xdb.set_trace()
    print(fN_model)

    # Plot with Data
    fN_data.tst_fn_data(fN_model=fN_model)
