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


#    Calculate tau_effective for the Lyman series using the EW
#    approximation (e.g. Zuo 93)
def ew_teff_lyman(ilambda, zem, fN_model, NHI_MIN=11.5, NHI_MAX=22.0, N_eval=5000,
                  EW_spline=None, bval=24., fNz=False, cosmo=None, debug=False):
    """ tau effective (follows ew_teff_lyman.pro from XIDL)
       teff = ew_teff_lyman(3400., 2.4)

    Parameters:
      ilambda: float
        Observed wavelength 
      zem: float 
        Emission redshift of the source [sets which Lyman lines are included]
      bva: float
         -- Characteristics Doppler parameter for the Lya forest
         -- [Options: 24, 35 km/s]
      NHI_MIN: float
         -- Minimum log HI column for integration [default = 11.5]
      NHI_MAX: float
         -- Maximum log HI column for integration [default = 22.0]
      fNz: Boolean (False)
         -- Inputs f(N,z) instead of f(N,X)
      cosmo: astropy.cosmology (None)
         -- Cosmological model to adopt (as needed)

    Returns:
      teff: 
        Total effective opacity of all lines contributing

    ToDo:
      1. Parallelize the Lyman loop

    JXP 07 Nov 2014
    """
    # Lambda
    if not isinstance(ilambda,float):
        raise ValueError('igm.tau_eff: ilambda must be a flaat for now')
    Lambda = ilambda
    if not isinstance(Lambda,u.quantity.Quantity):
        Lambda = Lambda * u.AA # Ang

    # Read in EW spline (if needed)
    if EW_spline == None:
        if int(bval) == 24: EW_FIL = xa_path+'/igm/EW_SPLINE_b24.p'
        elif int(bval) == 35: EW_FIL = os.environ.get('XIDL_DIR')+'/IGM/EW_SPLINE_b35.fits'
        else: 
            raise ValueError('igm.tau_eff: Not ready for this bvalue %g' % bval)
        EW_spline = pickle.load(open(EW_FIL,"rb"))

    # Lines
    wrest = tau_eff_llist()

    # Find the lines
    gd_Lyman = wrest[(Lambda/(1+zem)) < wrest]
    nlyman = len(gd_Lyman) 
    if nlyman == 0:
        print('igm.tau_eff: No Lyman lines covered at this wavelength')
        return -1

    # N_HI grid
    lgNval = NHI_MIN + (NHI_MAX-NHI_MIN)*np.arange(N_eval)/(N_eval-1) # Base 10 
    dlgN = lgNval[1]-lgNval[0]
    Nval = 10.**lgNval
    teff_lyman = np.zeros(nlyman)

    # Loop on the lines
    for qq,line in enumerate(gd_Lyman): # Would be great to do this in parallel... 
                             # (Can pack together and should)
        # Redshift
        zeval = (Lambda / line) - 1
        if fNz is False:
            if cosmo not in locals():
                cosmo = FlatLambdaCDM(H0=70, Om0=0.3) # Vanilla
            #dxdz = (np.fabs(xigmu.cosm_xz(zeval-0.1, cosmo=cosmo)- 
            #            xigmu.cosm_xz(zeval+0.1,cosmo=cosmo)) / 0.2 )
            #xdb.set_trace()
            dxdz = xigmu.cosm_xz(zeval,cosmo=cosmo,flg=1)
        else: dxdz = 1. # Code is using f(N,z)
        #print('dxdz = %g' % dxdz)

        # Get EW values (could pack these all together)
        idx = np.where(EW_spline['wrest'] == line)[0]
        if len(idx) != 1:
            raise ValueError('tau_eff: Line %g not included or over included?!' % line)
        restEW = interpolate.splev(lgNval, EW_spline['tck'][idx], der=0)

        # dz
        dz = (restEW*u.AA) * (1+zeval) / line

        # Evaluate f(N,X)
        log_fnX = fN_model.eval(lgNval,zeval)

        # Sum
        intgrnd = 10.**(log_fnX) * dxdz * dz * Nval
        teff_lyman[qq] = np.sum(intgrnd) * dlgN * np.log(10.)
        #xdb.set_trace()

        # Debug
        if debug==True:
            xdb.xplot(lgNval, np.log10(10.**(log_fnX) * dxdz * dz * Nval))
            #x_splot, lgNval, total(10.d^(log_fnX) * dxdz * dz * Nval,/cumul) * dlgN * alog(10.) / teff_lyman[qq], /bloc
            #printcol, lgnval, log_fnx, dz,  alog10(10.d^(log_fnX) * dxdz * dz * Nval)
            #writecol, 'debug_file'+strtrim(qq,2)+'.dat',  lgNval, restEW, log_fnX
            xdb.set_trace()

    #xdb.set_trace()
    return np.sum(teff_lyman)
     
    

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
    from barak import voigt as bv

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
        vgt = bv.voigt(a, uval)

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

    # Make EW spline file
    #mk_ew_lyman_spline(24.)

    # read f(N)
    
    from xastropy.igm.fN import model as xifm

    #xdb.set_trace()
    fN_model = xifm.default_model()
    print(fN_model)

    # tau_eff
    lamb = 1215.6701*(1+2.4)
    teff = ew_teff_lyman(lamb, 2.5, fN_model, NHI_MIN=12., NHI_MAX=17.)
    print('teff at z=2.4 :: %g' % teff)
    #teff = ew_teff_lyman(3400., 2.4, fN_model)
    #print('teff at 3400A = %g' % teff)
