"""
#;+ 
#; NAME:
#; igm.tau_eff
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for tau effective
#;   07-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function

import os, imp, pickle

import numpy as np
from scipy import interpolate

from xastropy.igm import igm_utils as xigmu
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u

from xastropy.xutils import xdebug as xdb

# Path for xastropy
xa_path = imp.find_module('xastropy')[1]


# def ew_teff_lyman -- Calcualte tau_effective for the HI Lyman series
# def mk_ew_lyman_spline -- Generates a Pickle file for EW splines
# def teff_obs(z)

#    Calculate tau_effective for the Lyman series using the EW
#    approximation (e.g. Zuo 93)


    

# ###
# Generate a pickle file of a Spline of EW vs NHI for the Lyman series
def mk_ew_lyman_spline(bval,ew_fil=None):
    """ Generate a pickle file of a Spline of EW vs NHI for the Lyman series

    Parameters:
      bval: float
        Doppler parameter (km/s)
      ew_fil: string ('EW_SPLINE_b##.p')
        Name of output pickle file
    
    """

    
    from astropy import constants as const
    from xastropy.spec import abs_line as xsab
    from xastropy.spec import voigt as xsv

    # Outfil
    if ew_fil == None:
        ew_fil = 'EW_SPLINE_b'+str(int(bval))+'.p'

    # Units
    if not isinstance(bval,u.quantity.Quantity):
        bval = bval * u.km/u.s # km/s

    # NHI
    nspl = 100
    log_NHI = 11.0 + 11*np.arange(nspl)/(nspl-1.)

    # Lines
    wrest= tau_eff_llist()

    # Output
    outp = {'wrest': wrest, 'tck': []}

    # Setup
    nvel = 60001
    velo = (-30000. + np.arange(nvel,dtype='float64'))*u.km/u.s # km/s
    dvel = 1. * u.km/u.s # km/s
    uval = velo / bval 

    # Loop
    for cnt,line in enumerate(wrest): 

        # Get atomic data
        abl_data = xsab.abs_line_data(line.value)

        # Wave array
        dwv = dvel.to(u.cm/u.s) * line / const.c.cgs  # Ang

        # Voigt
        vd = (bval/line).to(u.Hz)  # Frequency
        a = abl_data['gamma'] / (12.56637 * vd.value)
        vgt = xsv.voigtking(uval,a)

        # tau
        tau = 0.014971475*abl_data['fval']*vgt/vd  # Normalized to N_HI = 1 cm^-2

        # Flux
        tau_array = np.outer(tau, 10.**log_NHI)
        fx = np.exp(-1.*tau_array)

        # EW
        EW = np.sum(1.-fx, 0) * dwv
        #EW_SPLINE[qq].EW = EW
     
        # Spline
        #EW_SPLINE[qq].splint = spl_init(NHI, EW, /double)
        tck = interpolate.splrep(log_NHI, EW)#, s=0)

        # Check?
        chk=False
        #if line == (1215.6701*u.AA): 
        #    xdb.set_trace()
        #    chk=True
        if chk:
            from matplotlib import pyplot as plt
            plt.clf()
            plt.plot(log_NHI,EW,'o')
            # Spline
            xnew = np.linspace(np.amin(log_NHI),np.amax(log_NHI), nspl*10)
            ynew = interpolate.splev(xnew, tck, der=0)
            plt.plot(xnew, ynew, '-')
            plt.show()

        # Output
        print('line = %g' % line.value)
        outp['tck'].append(tck)

    # Write
    print('Writing %s' % ew_fil)
    pickle.dump( outp, open( ew_fil, "wb" ) )


# Line list for tau_eff HI Lyman calculatons
def tau_eff_llist():

    # Imports

    # Assumed Line list
    wrest = (np.array([1215.6701,  1025.7223,       972.53680,       949.74310,       937.80350, 
            930.74830,  926.22570,       923.15040,       920.96310,       919.35140, 
            918.12940,  917.18060,       916.42900,       915.82400,       915.32900,
            914.91900,  914.57600,       914.28600,       914.03900,       913.82600, 
            913.64100,  913.48000,       913.33900,       913.21500,       913.10400, 
            913.00600,  912.91800,       912.83900,       912.76800,       912.70300, 
            912.64500],dtype='float64' )) * u.AA
    wrest.sort()

    return wrest



## #################################    
## #################################    
## TESTING
## #################################    
if __name__ == '__main__':

    '''
    # Make EW spline file
    mk_ew_lyman_spline(24.)
    '''

    from xastropy.igm.fN import model as xifm
    import multiprocessing

    #xdb.set_trace()
    # read f(N)
    fN_model = xifm.default_model()
    print(fN_model)

    # tau_eff
    #tst_wv = tau_eff_llist()
    tst_wv = np.arange(915.,1255,1.)
    #lamb = 1215.6701*(1+2.4)
    adict = []
    for wrest in tst_wv:
        tdict = dict(ilambda=wrest*(1+2.4), zem=2.5, fN_model=fN_model)
        adict.append(tdict)

    pool = multiprocessing.Pool(4) # initialize thread pool N threads
    ateff = pool.map(map_etl, adict)
    # Plot
    xdb.xplot(tst_wv,np.exp(-np.array(ateff)))
    #xdb.set_trace()
    #teff = ew_teff_lyman(lamb, 2.5, fN_model, NHI_MIN=12., NHI_MAX=17.)
    #print('teff at z=2.4 :: %g' % teff)
    #teff = ew_teff_lyman(3400., 2.4, fN_model)
    #print('teff at 3400A = %g' % teff)
