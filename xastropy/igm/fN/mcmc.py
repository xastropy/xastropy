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
from xastropy.igm.fN import data as fN_data

xa_path = imp.find_module('xastropy')[1]

#######################################
#   DEFINE THE MODEL
#######################################

def set_model(flg=0):
	"""
	Generate the f(N) model

    Parameters
    ----------
    flg: int (0)
      Sets the type of model
      0 = Hermite Spline
      1 = Inoue model

    Returns
    -------
    fN_model Class

    JXP on 27 Nov 2014
	"""
    if flg==0: # I may choose to pickle a few of these
        fN_model = xifm.default_model(recalc=True,use_mcmc=True) # Hermite Spline
    elif flg==1:
        fN_model = fNmodel.fN_Model('Gamma')
    else: 
        raise ValueError('mcmc.set_model: Not ready for this type of fN model {:d}'.format(flg))
	
	return fN_model

#######################################
#          READ IN THE DATA
#######################################

def set_fNdata():
	"""
	Load up f(N) data

    Parameters
    ----------

    Returns
    -------
    fN_data :: List of fN_Constraint Classes

    JXP on 27 Nov 2014
	"""

    fn_file = xa_path+'/igm/fN/fn_constraints_z2.5_vanilla.fits'
    k13r13_file = xa_path+'/igm/fN/fn_constraints_K13R13_vanilla.fits'
    n12_file = xa_path+'/igm/fN/fn_constraints_N12_vanilla.fits'
    all_fN_cs = fN_data.fn_data_from_fits([fn_file,k13r13_file,n12_file])

    return all_fN_cs

##########################################
#   Prepare the variables and their limits
##########################################
def set_pymc_var(fN_model,lim=2.):
	"""
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
	"""
    print "Preparing variables and their limits"
    if fN_model.fN_mtype == 'Hspline': 
        parm=np.array([])
        for ii in range(len(fN_model.parm)):
            nm = 'p'+str(ii)
            doc = 'SplinePointNHI_'+str(fN_model.pivots[ii])
            parm = np.append(parm, pymc.Uniform(nm, lower=fN_model.parm[ii]-lim,
                                                upper=fN_model.parm[ii]+lim, doc=doc))
    else:
        raise ValueError('mcmc: Not ready for this type of fN model {:s}'.format(fN_model.fN_mtype))
    # Return
    return parm


##########################################
#   The wrapper for the MCMC.
##########################################
'''
@pymc.deterministic(plot=False)
def model(parm=parm):
	""" Gaussian Model """
	p = np.array([parm[0],parm[1],6563.0,parm[2]])
	out = func_gauss_oned(waves, p)
	return out

##########################################
#   Setup the data array to be fitted
##########################################

data = pymc.Normal('data', mu=model, tau=1.0/error**2, value=flux, observed=True)

#######################################
#   RUN THE MCMC
#######################################

MC = pymc.MCMC([parm, model, data])
# Run a total of 40000 samples, but ignore the first 10000.
# Verbose just prints some details to screen.
MC.sample(40000, 10000, verbose=2)

#######################################
#   PRINT THE RESULTS
#######################################

# Print the best values and their errors
MCMC_errors.print_errors(MC)

# Draw a contour plot with 1 & 2 sigma errors
MCMC_errors.draw_contours(MC, 'p0', 'p1')

# Save the individual distributions to a file to check convergence
pymc.Matplot.plot(MC)
'''

if __name__ == '__main__':

    # Generate the model
    fN_model = set_model()
    xdb.set_trace()
