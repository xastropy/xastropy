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
import os, pickle, imp
from scipy import interpolate as scii

from xastropy.xutils import xdebug as xdb
from xastropy.spec import abs_line, voigt
from xastropy.stats import mcmc

from astropy.io import fits

# Path for xastropy
xa_path = imp.find_module('xastropy')[1]

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
                 param=None, zpivot=2.4, gamma=1.5):
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
        self.zpivot = zpivot
        self.gamma = gamma
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

#########
def default_model(recalc=False, pckl_fil=None, use_mcmc=False):
    """
    Pass back a default fN_model from Prochaska+13
      Tested against XIDL code by JXP on 09 Nov 2014

    Parameters:
    recalc: boolean (False)
      Recalucate the default model
    use_mcmc: boolean (False)
      Use the MCMC chain to generate the model
    """
    if pckl_fil==None:
        pckl_fil = xa_path+'/igm/fN/fN_model_P13.p'

    if recalc is True:
        
        if use_mcmc == True:
            # MCMC Analysis
            chain_file = (os.environ.get('DROPBOX_DIR')+
                        'IGM/fN/MCMC/mcmc_spline_k13r13o13n12_8.fits.gz')
            outp = mcmc.chain_stats(chain_file)
    
            # Build a model
            NHI_pivots = [12., 15., 17.0, 18.0, 20.0, 21., 21.5, 22.]
            fN_model = fN_Model('Hspline', zmnx=(0.5,3.0),
                            pivots=NHI_pivots, param=outp['best_p'])
        else:
            # Input the f(N) at z=2.4
            fN_file = (os.environ.get('DROPBOX_DIR')+
                        'IGM/fN/fN_spline_z24.fits.gz')
            hdu = fits.open(fN_file)
            fN_data = hdu[1].data
            #xdb.set_trace()
            # Instantiate
            fN_model = fN_Model('Hspline', zmnx=(0.5,3.0),
                            pivots=np.array(fN_data['LGN']).flatten(),
                            param=np.array(fN_data['FN']).flatten())
        # Write
        print('default_model: Writing %s' % pckl_fil)
        pickle.dump( fN_model, open( pckl_fil, "wb" ), -1)
    else: fN_model = pickle.load( open( pckl_fil, "rb" ) )
        
    # Return
    return fN_model







## #################################    
## #################################    
## TESTING
## #################################    
if __name__ == '__main__':

    from matplotlib import pyplot as plt
    
    from xastropy.igm.fN import data as fN_data
    from xastropy.igm.fN import model as xifm

    flg_test = 6
    
    if (flg_test % 2) == 1:
        # MCMC Analysis
        chain_file = os.environ.get('DROPBOX_DIR')+'IGM/fN/MCMC/mcmc_spline_k13r13o13n12_8.fits.gz'
        outp = mcmc.chain_stats(chain_file)

        # Build a model
        NHI_pivots = [12., 15., 17.0, 18.0, 20.0, 21., 21.5, 22.]
        fN_model = xifm.fN_Model('Hspline', zmnx=(0.5,3.0),
                            pivots=NHI_pivots, param=outp['best_p'])
        #xdb.set_trace()
        print(fN_model)

    # Compare default against P+13
    if (flg_test % 4) >= 2:
        fN_model = xifm.default_model()
        #p13_pivots = ([12.000000,       15.000000,       17.000000,      18.000000,
        #        20.000000,       21.000000,       21.500000,       22.000000])
        #p13_fN = ([-9.7233782,      -14.411575,      -17.941498,      -19.392904,
        #        -21.282812,      -22.820860,      -23.945236,      -25.502331])
        p13_file = (os.environ.get('DROPBOX_DIR')+'IGM/fN/fN_spline_z24.fits.gz')
        hdu = fits.open(p13_file)
        p13_data = hdu[1].data
        
        plt.clf()
        plt.scatter(p13_data['LGN'],p13_data['FN'])
        #plt.plot(p13_data['LGN'],p13_data['FN'], '-')
        xplt = np.linspace(12., 22, 10000)
        yplt = fN_model.eval(2.4,xplt)
        plt.plot(xplt, yplt, '-', color='red')
        #plt.plot(xplt, yplt, '-')
        plt.show()
        #xdb.set_trace()


    # Reproduce the main figure from P14
    # Plot with Data
    if (flg_test % 8) >= 4:
        fN_data.tst_fn_data(fN_model=fN_model)
