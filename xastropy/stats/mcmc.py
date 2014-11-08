"""
#;+ 
#; NAME:
#; stats.mcmc
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for MCMC calcualtions
#;   07-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""

from __future__ import print_function

import numpy as np
import os

from xastropy.xutils import xdebug as xdb
from xastropy.spec import abs_line, voigt

from astropy.io import fits


def chain_stats(chain_file, burn_frac=0.3, cl=0.683):
    """ Turn an MCMC chain into stats
    Port of x_mcmc_chain_stats from XIDL

    Parameters:
      chain_file: string
          Name of MCMC file
      burn_frac: float (0.3)
          Fraction of chain to burn
      cl: float (0.683)
          Confidence interval

    Returns:
      A dictionary with the key outputs

    JXP 07 Nov 2014
    """

    # Read
    hdu = fits.open(chain_file)
    chain = hdu[0].data
    like = hdu[1].data

    # Param
    nparm = chain.shape[-1]
    if len(chain.shape) < 3: nchain = 1
    else: nchain = chain.shape[0]

    # Output
    outp = {}

    # Burn
    burn = int( np.round(chain.shape[1] * burn_frac ) )
    chain = chain[:,burn:,:]    # Burn
    like = like[:,burn:]    # Burn
    sz = chain.shape

    # Reshape
    if len(sz) == 3:
        chain = chain.reshape(sz[0]*sz[1],sz[2])
        chain = chain.transpose()
        sz = chain.shape
        like = like.flatten()

    # Maximize
    imx = np.argmax(like)
    outp['best_p'] = chain[:,imx]
    outp['sig'] = np.zeros((sz[0],2))

    # Confidence limits (68%)
    cnt = 0
    for row in chain: #qq=0L,nparm-1 do begin
        #; Sort
        srt = np.sort(row)
        # Simple CDF
        lowv = srt[np.round(sz[1]*(1-cl)/2.)]
        outp['sig'][cnt,0] = np.fabs(outp['best_p'][cnt] - lowv)
        #
        hiv = srt[np.round(sz[1]*(1.-((1-cl)/2.)))]
        outp['sig'][cnt,1] = np.fabs(hiv-outp['best_p'][cnt])
        cnt += 1

    return outp




    
## #################################    
## #################################    
## TESTING
## #################################    
if __name__ == '__main__':

    # Test mcmc_chain_stats
    mcmc_file = os.environ.get('DROPBOX_DIR')+'IGM/fN/MCMC/mcmc_spline_k13r13o13n12_8.fits.gz'
    outp = chain_stats(mcmc_file)
    print(outp)

