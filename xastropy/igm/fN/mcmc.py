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
from scipy import interpolate as scii

from xastropy.xutils import xdebug as xdb
from xastropy.igm.fN import model as xifm
from xastropy.igm.fN import data as xifd
from xastropy.igm import tau_eff

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
        tmp = [13., 15., 17., 21.5, 22.]
        val = sfN_model.eval(tmp, 2.5)
        # Not using this! -- JXP 29 Jan 2015
        #xdb.set_trace()
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

def set_fn_data(sources=None, extra_fNc=None):
    '''
    Load up f(N) data

    Parameters
    ----------

    Returns
    -------
    fN_data :: List of fN_Constraint Classes

    JXP on 27 Nov 2014
    '''
    if sources is None:
        sources = ['OPB07', 'OPW12', 'OPW13', 'K05', 'K13R13', 'N12']
    else:
        sources = [src.upper() for src in sources]

    fn_file = xa_path+'/igm/fN/fn_constraints_z2.5_vanilla.fits'
    k13r13_file = xa_path+'/igm/fN/fn_constraints_K13R13_vanilla.fits'
    n12_file = xa_path+'/igm/fN/fn_constraints_N12_vanilla.fits'
    all_fN_cs = xifd.fn_data_from_fits([fn_file,k13r13_file,n12_file])

    # Add on, e.g. user-supplied
    if not extra_fNc is None:
        all_fN_cs.append(extra_fNc)

    # Cut
    #fN_cs = [fN_c for fN_c in all_fN_cs
    #         if ((fN_c.ref != 'K02') & (fN_c.ref != 'PW09'))]

    # Include good data sources
    fN_cs = []
    for fN_c in all_fN_cs:
        # In list?
        if fN_c.ref in sources:
            print('Using {:s} as a constraint'.format(fN_c.ref))
            # Append
            fN_cs.append(fN_c)
            # Pop
            idx = sources.index(fN_c.ref)
            sources.pop(idx)
    
    # Check that all the desired sources were used
    if len(sources) > 0:
        xdb.set_trace()

    #xdb.set_trace()

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
    # Deal with model Type
    if fN_model.fN_mtype == 'Hspline': 
        iparm=np.array([])
        # Loop on parameters to create an array of pymc Stochatsic variable objects
        for ii in range(len(fN_model.param)):
            nm = str('p')+str(ii)
            doc = str('SplinePointNHI_')+str(fN_model.pivots[ii])
            #iparm = np.append(iparm, pymc.Uniform(nm, lower=fN_model.param[ii]-lim,
            #                                    upper=fN_model.param[ii]+lim, doc=doc))
            #xdb.set_trace()
            iparm = np.append(iparm, pymc.Normal(nm, mu=fN_model.param[ii], tau=1./0.025, doc=doc))
    else:
        raise ValueError('mcmc: Not ready for this type of fN model {:s}'.format(fN_model.fN_mtype))
    # Return
    return iparm




##########################################
# Main run call
##########################################
def run(fN_cs, fN_model, parm, debug=0):

    #
    pymc_list = [parm]

    # Parse data and combine as warranted
    all_NHI = []
    all_fN = []
    all_sigfN = []
    all_z = []
    flg_teff = 0
    flg_LLS = 0
    for fN_c in fN_cs: 
        # Standard f(N)
        if fN_c.fN_dtype == 'fN':
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
        elif fN_c.fN_dtype == 'teff': # teff_Lya
            if flg_teff:
                raise ValueError('Only one teff allowed for now!')
            else:
                flg_teff = 1
            teff=float(fN_c.data['TEFF'])
            D_A = 1. - np.exp(-1. * teff)
            SIGDA_LIMIT = 0.1  # Allows for systemtics and b-value uncertainty
            sig_teff = np.max([fN_c.data['SIG_TEFF'], (SIGDA_LIMIT*teff)])
            teff_zeval = float(fN_c.data['Z_TEFF'])

            # Save input for later usage
            teff_input = (teff_zeval, fN_c.data['NHI_MNX'][0], fN_c.data['NHI_MNX'][1])
        elif fN_c.fN_dtype == 'l(X)': # teff_Lya
            if flg_LLS:
                raise ValueError('Only one teff allowed for now!')
            else:
                flg_LLS = 1
            LLS_lx = fN_c.data['LX']
            LLS_siglx = fN_c.data['SIG_LX']
            LLS_input = (fN_c.zeval, fN_c.data['TAU_LIM'])
            
    # 
    fN_input = (np.array(all_NHI), np.array(all_z))
    #flg_teff = 0

    #######################################
    #   Generate the Models
    #######################################

    # Define f(N) model for PyMC
    @pymc.deterministic(plot=False)
    def pymc_fn_model(parm=parm):
        # Set parameters
        fN_model.param = parm
        fN_model.model = scii.PchipInterpolator(fN_model.pivots, fN_model.param)
        #
        log_fNX = fN_model.eval( fN_input, 0. )
        #
        return log_fNX
    pymc_list.append(pymc_fn_model)

    # Define teff model for PyMC
    if flg_teff:
        @pymc.deterministic(plot=False)
        def pymc_teff_model(parm=parm):
            # Set parameters
            fN_model.param = parm
            fN_model.model = scii.PchipInterpolator(fN_model.pivots, fN_model.param)
            # Calculate teff
            model_teff = tau_eff.ew_teff_lyman(1215.6701*(1+teff_input[0]), teff_input[0]+0.1,
                                               fN_model, NHI_MIN=teff_input[1], NHI_MAX=teff_input[2])
            return model_teff
        pymc_list.append(pymc_teff_model)

    # Define l(X)_LLS model for PyMC
    if flg_LLS:
        @pymc.deterministic(plot=False)
        def pymc_lls_model(parm=parm): 
            # Set parameters 
            fN_model.param = parm
            fN_model.model = scii.PchipInterpolator(fN_model.pivots, fN_model.param)
            # Calculate l(X)
            lX = fN_model.calc_lox(LLS_input[0], 
                                    17.19+np.log10(LLS_input[1]), 22.) 
            return lX
        pymc_list.append(pymc_lls_model)

    #######################################
    #   Generate the Data
    #######################################

    # Define f(N) data for PyMC
    fNvalue=np.array(all_fN)
    pymc_fN_data = pymc.Normal(str('fNdata'), mu=pymc_fn_model, tau=1.0/np.array(all_sigfN)**2,
                               value=fNvalue, observed=True)
    pymc_list.append(pymc_fN_data)

    # Define teff data for PyMC
    if flg_teff:
        pymc_teff_data = pymc.Normal(str('teffdata'), mu=pymc_teff_model, tau=1.0/np.array(sig_teff)**2,
                                value=teff, observed=True)
        pymc_list.append(pymc_teff_data)

    # Define l(X)_LLS model for PyMC
    if flg_LLS:
        pymc_lls_data = pymc.Normal(str('LLSdata'), mu=pymc_lls_model, tau=1.0/np.array(LLS_siglx)**2,
                                value=LLS_lx, observed=True)
        pymc_list.append(pymc_lls_data)


    #######################################
    #   RUN THE MCMC
    #######################################


    MC = pymc.MCMC(pymc_list)#,verbose=2)
    # Force step method to be Metropolis!
    for ss in MC.stochastics-MC.observed_stochastics:
        MC.use_step_method(pymc.Metropolis, ss, proposal_sd=0.025, proposal_distribution='Normal')
    #xdb.set_trace()

    # Run a total of 40000 samples, but ignore the first 10000.
    # Verbose just prints some details to screen.
    #xdb.set_trace()
    #MC.sample(20000, 3000, verbose=2, tune_interval=500)
    MC.sample(2000, 300, verbose=2, tune_interval=200)
    #MC.isample(10000, 1000, verbose=2)

    #######################################
    #   PRINT THE RESULTS
    #######################################

    # Print the best values and their errors
    best_pval = print_errors(MC)

    fN_model.param = best_pval
    fN_model.model = scii.PchipInterpolator(fN_model.pivots, fN_model.param)
    if debug:
        xifd.tst_fn_data(fN_model=fN_model)
        xdb.xhist(MC.trace(str('p0'))[:])
        xdb.set_trace()

    
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

    all_pval = []
    for ival in range(keys_size):
        teststr = str('p'+str(ival))
        #xdb.set_trace()
        if not teststr in keys:
            continue
        try: pval, perr1, perr2 = geterrors(MC.trace(teststr)[:])
        except: continue
        print('{:s} {:5.4f} +{:5.4f}-{:5.4f} (1sig) +{:5.4f}-{:5.4f} (2sig)'.format(
            teststr, pval, perr1[0], perr1[1], perr2[0], perr2[1]))
        all_pval.append(pval)
        #ival += 1
    return all_pval

def mcmc_main(flg_model=0, flg_plot=0):
    '''
    flg_model = Flag controlling the f(N) model fitted
       0: JXP spline
       1: Inoue+14 functional form
    '''

    import argparse

    # PARSE 
    parser = argparse.ArgumentParser(description='MCMC for f(N)')
    parser.add_argument("model", type=int, help="Model flag (0=JXP, 1=Inoue+14)")
    parser.add_argument("--plot", help="Plot the final model", action="store_true")
    #parser.add_argument("-id_dir", type=str,
    #                    help="Directory for ID files (ID_LINES is default)")
    
    args = parser.parse_args()

    flg_model = args.model
    if args.plot:
        flg_plot = 1

    # ##########################
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

    # Plot?
    if flg_plot:
        xifd.tst_fn_data(fN_model=fN_model)

#####
if __name__ == '__main__':
    import sys

    if len(sys.argv) == 1: # TESTING
        sys.argv.append('0')
        mcmc_main(flg_plot=1)
    else:
        mcmc_main()

    # Set model
    #xdb.set_trace()
    print('mcmc: All done')
