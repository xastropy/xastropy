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
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os

from xastropy.xutils import xdebug as xdb

from astropy.io import fits





# For Alix
def test():
    import time
    return time.time()


    
## #################################    
## #################################    
## TESTING
## #################################    
if __name__ == '__main__':

    # Test mcmc_chain_stats
    mcmc_file = os.environ.get('DROPBOX_DIR')+'IGM/fN/MCMC/mcmc_spline_k13r13o13n12_8.fits.gz'
    outp = chain_stats(mcmc_file)
    print(outp)

