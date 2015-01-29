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

    # Cut
    fN_cs = [fN_c for fN_c in all_fN_cs
             if ((fN_c.ref != 'K02') & (fN_c.ref != 'PW09'))]

    return fN_cs

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
# Main run call
##########################################
def run(fN_cs, fN_model, parm):

    #
    pymc_list = [parm]

    # Combine f(N) data
    all_NHI = []
    all_fN = []
    all_sigfN = []
    all_z = []
    for fN_c in fN_cs: 
        if fN_c.fN_dtype != 'fN':
            continue
        ip = range(fN_c.data['NPT'])
        val = np.where(fN_c.data['FN'][ip] > -90)[0] # Deal with limits later
        ipv = np.array(ip)[val]
        # Append the NHI
        NHI = np.median(fN_c.data['BINS'][:,ipv],0)
        all_NHI += list(NHI)
        # Append the f(N)
        all_fN += list(fN_c.data['FN'][ipv])
        # Append the Error
        fNerror = np.median(fN_c.data['SIG_FN'][:,ipv],0)
        all_sigfN += list(fNerror)
        # Append zeval
        for ii in range(len(ipv)):
            all_z.append(fN_c.zeval)
    fN_input = (np.array(all_NHI), np.array(all_z))

    # TEST
    xdb.set_trace()
    if False:
        log_fNX = fN_model.eval( fN_input, 0. )
        xdb.set_trace()

    # Define f(N) model for PyMC
    @pymc.deterministic(plot=False)
    def pymc_fn_model(parm=parm):
        # Set parameters
        fN_model.param = parm
        #
        log_fNX = fN_model.eval( fN_input, 0. )
        #
        return log_fNX
    pymc_list.append(pymc_fn_model)

    # Define f(N) data for PyMC
    fNvalue=np.array(all_fN)
    pymc_fN_data = pymc.Normal(str('fNdata'), mu=pymc_fn_model, tau=1.0/np.array(all_sigfN)**2,
                               value=fNvalue, observed=True)
    pymc_list.append(pymc_fN_data)


    #######################################
    #   RUN THE MCMC
    #######################################

    MC = pymc.MCMC(pymc_list)
    # Run a total of 40000 samples, but ignore the first 10000.
    # Verbose just prints some details to screen.
    MC.sample(40000, 10000, verbose=2)

    #######################################
    #   PRINT THE RESULTS
    #######################################

    # Print the best values and their errors
    print_errors(MC)
    
    # Draw a contour plot with 1 & 2 sigma errors
    #MCMC_errors.draw_contours(MC, 'p0', 'p1')

    # Save the individual distributions to a file to check convergence
    #pymc.Matplot.plot(MC)

def geterrors(array):
	arrsort = np.sort(array)
	arrsize = np.size(array)
	value = arrsort[int(round(0.5*arrsize))]
	err1 = np.array([value-arrsort[int(round(0.15866*arrsize))],arrsort[int(round(0.84134*arrsize))]-value])
	err2 = np.array([value-arrsort[int(round(0.02275*arrsize))],arrsort[int(round(0.97725*arrsize))]-value])
	return value, err1, err2

def print_errors(MC):
    keys = MC.stats().keys()
    keys_size = len(keys)
    ival=0
    xdb.set_trace()

    for i in range(0,keys_size):
        teststr = 'p'+str(ival)
        try: pval, perr1, perr2 = geterrors(MC.trace(teststr)[:])
        except: continue
        print('{:5.4f +{:5.4f}-{:5.4f} (1sig) +{:5.4f}-{:5.4f} (2sig)'.format(
            pval, perr1[0], perr1[1], perr2[0], perr2[1]))
        ival += 1


#####
if __name__ == '__main__':

    # Set Data
    fN_data = set_fn_data()

    # Set f(N) functional form 
    fN_model = set_fn_model()

    # Set variables
    parm = set_pymc_var(fN_model)
    fN_model.param = np.array([iparm.value for iparm in parm])

    # Check plot
    if False:
        xifd.tst_fn_data(fN_model=fN_model)

    # Run
    run(fN_data, fN_model, parm)

    # Set model
    xdb.set_trace()
    print('mcmc: All done')
