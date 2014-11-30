"""
#;+ 
#; NAME:
#; fN.mcmc
#;    Version 1.0
#;
#; PURPOSE:
This module will fit f(N) data
using a Markov chain Monte Carlo (MCMC) approach.
#;   27-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os, pickle, imp
import numpy as np
import pymc
#import MCMC_errors

from xastropy.xutils import xdebug as xdb
from xastropy.igm.fN import model as xifm
from xastropy.igm.fN import data as xifd

xa_path = imp.find_module('xastropy')[1]

#######################################
#   DEFINE THE MODEL
#######################################

def set_fn_model(flg=0):
    '''
	Load up f(N) data

    Parameters
    ----------

    Returns
    -------
    fN_data :: List of fN_Constraint Classes

    JXP on 27 Nov 2014
    '''
    if flg==0: # I may choose to pickle a few of these
        sfN_model = xifm.default_model(recalc=True,use_mcmc=True) # Hermite Spline
    elif flg==1:
        sfN_model = fNmodel.fN_Model('Gamma')
    else: 
        raise ValueError('mcmc.set_model: Not ready for this type of fN model {:d}'.format(flg))
    #
    #print(sfN_model)
    return sfN_model


#######################################
#          READ IN THE DATA
#######################################

def set_fn_data():
    '''
    Load up f(N) data

    Parameters
    ----------

    Returns
    -------
    fN_data :: List of fN_Constraint Classes

    JXP on 27 Nov 2014
    '''
    fn_file = xa_path+'/igm/fN/fn_constraints_z2.5_vanilla.fits'
    k13r13_file = xa_path+'/igm/fN/fn_constraints_K13R13_vanilla.fits'
    n12_file = xa_path+'/igm/fN/fn_constraints_N12_vanilla.fits'
    all_fN_cs = xifd.fn_data_from_fits([fn_file,k13r13_file,n12_file])

    return all_fN_cs

##########################################
#   Prepare the variables and their limits
##########################################
def set_pymc_var(fN_model,lim=2.):
    '''
    Generate pymc variables

    Parameters
    ----------
    fN_model
    lim: float (2.)
      Range limits for Uniform stochastic value

    Returns
    -------
    Array of pymc Stochastic variables

    JXP on 27 Nov 2014
    '''
    if fN_model.fN_mtype == 'Hspline': 
        iparm=np.array([])
        for ii in range(len(fN_model.param)):
            nm = str('p')+str(ii)
            doc = str('SplinePointNHI_')+str(fN_model.pivots[ii])
            iparm = np.append(iparm, pymc.Uniform(nm, lower=fN_model.param[ii]-lim,
                                                upper=fN_model.param[ii]+lim, doc=doc))
    else:
        raise ValueError('mcmc: Not ready for this type of fN model {:s}'.format(fN_model.fN_mtype))
    # Return
    return iparm

##########################################
#   PyMC model for f(N) data
##########################################
@pymc.deterministic(plot=False)
def pymc_fn_model(parm=parm):
    #
    fN_model.param = np.array([iparm.value for iparm in parm])
	out = func_gauss_oned(waves, p)
	return out





#####
if __name__ == '__main__':

    # Generate the f(N) model
    fN_model = set_fn_model()
    print(fN_model)

    # Set Data
    fN_data = set_fn_data()

    # Set variables
    parm = set_pymc_var(fN_model)

    # Check plot
    if False:
        fN_model.param = np.array([iparm.value for iparm in parm])
        xifd.tst_fn_data(fN_model=fN_model)

    # Build PyMC model for standard f(N) data

    # Set model
    xdb.set_trace()
    print('mcmc: All done')
