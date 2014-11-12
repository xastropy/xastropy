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
from xastropy.igm import igm_utils as igmu
from xastropy.atomic import ionization as xai

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
    # l(X)
    def calc_lox(self, z, NHI_min, NHI_max=None, neval=10000, cumul=False):
        """ Calculate l(X) over an N_HI interval

        Parameters:
        z: float
          Redshift for evaluation
        NHI_min: float
          minimum NHI value
        NHI_max: float (Infinity)
          maximum NHI value for evaluation
        neval: int (10000)
          Discretization parameter
        cumul: boolean (False)
          Return a cumulative array?

        Returns:
        lX: float
          l(X) value

        JXP 10 Nov 2014
        """
        # Initial
        if NHI_max==None:
            NHI_max = 23.
            infinity=True
        else: infinity=False

        try:
            nz = len(z)
        except:
            nz=1
            z = [z]

        # Brute force (should be good to ~0.5%)
        lgNHI = NHI_min + (NHI_max-NHI_min)*np.arange(neval)/(neval-1.)
        dlgN = lgNHI[1]-lgNHI[0]

        # Evaluate f(N,X)
        lgfNX = np.zeros((neval,nz))
        lX = np.zeros(nz)
        for ii in range(nz): 
            lgfNX[:,ii] = self.eval(z[ii], lgNHI)

        # Sum
        for ii in range(nz): 
            lX[ii] = np.sum(10.**(lgfNX[:,ii]+lgNHI)) * dlgN * np.log(10.)
        if cumul==True: 
            if nz > 1: #; Have not modified this yet
                raise ValueError('fN.model: Not ready for this model type %s' % self.fN_mtype)
            cum_sum = np.cumsum(10.**(lgfNX+lgNHI)) * dlgN * np.log(10.)
        

        # Infinity?
        if infinity is True:
            neval2 = 1000L
            lgNHI2 = NHI_max + (99.-NHI_max)*np.arange(neval2)/(neval2-1.)
            dlgN = lgNHI2[1]-lgNHI2[0]
            lgfNX = np.zeros((neval2,nz))
            lX2 = np.zeros(nz)
            for ii in range(nz):
                lgfNX[:,ii] = self.eval(z[ii], lgNHI2)
                lX2[ii] = np.sum(10.**(lgfNX[:,ii]+lgNHI2)) * dlgN * np.log(10.)
            # 
            lX = lX + lX2

        # Return
        if nz==1:
            lX = lX[0]
        if cumul==True:
            return lX, cum_sum
        else:
            return lX
    ##
    # Evaluate
    def eval(self, z, NHI_values,vel_array=None):
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

        # Imports
        from astropy import constants as const

        # Evaluate
        if self.fN_mtype == 'Hspline': 
            log_fNX = self.model.__call__(NHI_values)
        else: 
            raise ValueError('fN.model: Not ready for this model type %s' % self.fN_mtype)

        # Redshift
        if vel_array==None:
            log_fNX += self.gamma * np.log10((1+z)/(1+self.zpivot))
        else:  # Here comes a grid..
            z_val = z + (1+z) * vel_array/(const.c.cgs.value/1e5)
            # 
            lgNHI_grid = np.outer(log_fNX, np.ones(len(z_val)))
            lenfX = len(log_fNX)
            #xdb.set_trace()
            # 
            z_grid1 = 10**( np.outer(np.ones(lenfX)*self.gamma,
                                     np.log10(1+z_val)) )  #; (1+z)^gamma
            z_grid2 = np.outer( np.ones(lenfX)*((1./(1+self.zpivot))**self.gamma), 
                        np.ones(len(z_val))  )
            log_fNX = lgNHI_grid + np.log10(z_grid1*z_grid2) 

        # Return
        return log_fNX
    ##
    # Mean Free Path
    def mfp(self, zem, neval=5000, cosmo=None, zmin=0.5):
        """ Calculate teff_LL 
        Effective opacity from LL absorption at z912 from zem

        Parameters:
        zem: float
          Redshift of source
        cosmo: astropy.cosmology (None)
          Cosmological model to adopt (as needed)
        neval: int (5000)
          Discretization parameter
        zmin: float (0.5)
          Minimum redshift in the calculation

        Returns:
        mfp : float
          Mean free path from zem (physical Mpc)

        JXP 11 Nov 2014
        """
        # Imports
        from astropy import constants as const
        from astropy import units as u
        from astropy import cosmology 

        # Cosmology
        if cosmo == None:
            cosmo = igmu.X_Cosmo(H0=70, Om0=0.3) # Vanilla

        # Calculate teff
        zval, teff_LL = self.teff_ll(zmin, zem, N_eval=neval, cosmo=cosmo)

        # Find tau=1
        imn = np.argmin( np.fabs(teff_LL-1.) )
        if np.fabs(teff_LL[imn]-1.) > 0.02:
            raise ValueError('fN.model.mfp: teff_LL too far from unity')

        # MFP
        mfp = np.fabs( cosmo.physical_distance(zval[imn]) -
                        cosmo.physical_distance(zem) ) # Mpc
        #xdb.set_trace()
        # Return
        return mfp
        

    ##
    # teff_LL
    def teff_ll(self, z912, zem, N_eval=5000, cosmo=None):
        """ Calculate teff_LL 
        Effective opacity from LL absorption at z912 from zem

        Parameters:
        z912: float
          Redshift for evaluation
        zem: float
          Redshift of source
        cosmo: astropy.cosmology (None)
          Cosmological model to adopt (as needed)
        N_eval: int (5000)
          Discretization parameter

        Returns:
        zval, teff_LL: array
          z values and Effective opacity from LL absorption from z912 to zem

        JXP 10 Nov 2014
        """
        # Imports
        from astropy import constants as const
        from astropy import units as u

        # NHI array
        lgNval = 11.5 + 10.5*np.arange(N_eval)/(N_eval-1.) #; This is base 10 [Max at 22]
        dlgN = lgNval[1]-lgNval[0]
        Nval = 10.**lgNval

        #; z array
        zval = z912 + (zem-z912)*np.arange(N_eval)/(N_eval-1.)
        dz = np.fabs(zval[1]-zval[0])

        teff_LL = np.zeros(N_eval)

        # dXdz
        dXdz = igmu.cosm_xz(zval, cosmo=cosmo, flg=1) 
        #if keyword_set(FNZ) then dXdz = replicate(1.,N_eval)

        # Evaluate f(N,X)
        velo = (zval-zem)/(1+zem) * (const.c.cgs.value/1e5) # Kludge for eval [km/s]

        log_fnX = self.eval(zem, lgNval, vel_array=velo)  
        log_fnz = log_fnX + np.outer(np.ones(N_eval), np.log10(dXdz))

        # Evaluate tau(z,N)
        teff_engy = (const.Ryd.to(u.eV,equivalencies=u.spectral()) /
                     ((1+zval)/(1+zem)) )
        sigma_z = xai.photo_cross(1,1,teff_engy)
        #xdb.set_trace()
        #sigma_z = teff_cross * ((1+zval)/(1+zem))**(2.75)  # Not exact but close
        tau_zN = np.outer(Nval, sigma_z)

        # Integrand
        intg = 10.**(log_fnz) * (1. - np.exp(-1.*tau_zN))

        # Sum
        sumz_first = False
        if sumz_first == False:
            #; Sum in N first
            N_summed = np.sum(intg * np.outer(Nval, np.ones(N_eval)),  0) * dlgN * np.log(10.)
            #xdb.set_trace()
            # Sum in z
            teff_LL = (np.cumsum(N_summed[::-1]))[::-1] * dz 
        #xdb.set_trace()

        # Debug
        debug=False
        if debug == True:
            #        x_splot, lgNval, alog10(10.d^(log_fnX) * dxdz * dz * Nval), /bloc
            #        x_splot, lgNval, total(10.d^(log_fnX) * dxdz * dz * Nval,/cumul) * dlgN * alog(10.) / teff_lyman[qq], /bloc
            #     printcol, lgnval, log_fnx, dz,  alog10(10.d^(log_fnX) * dxdz * dz * Nval)
            #     writecol, 'debug_file'+strtrim(qq,2)+'.dat', $
            #               lgNval, restEW, log_fnX
            xdb.set_trace()
        # Return
        return zval, teff_LL

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

    flg_test = 4+32
    
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
        fN_model = xifm.default_model()
        fN_data.tst_fn_data(fN_model=fN_model)

    # Check l(X)
    if (flg_test % 16) >= 8:
        fN_model = xifm.default_model()
        lX = fN_model.calc_lox(2.4, 17.19+np.log10(2.), 23.) 
        print('l(X) = %g' % lX)

    # Check teff_LL
    if (flg_test % 32) >= 16:
        fN_model = xifm.default_model()
        zval,teff_LL = fN_model.teff_ll(0.5, 2.45)
        xdb.xplot(zval,teff_LL)#,xlabel='z', ylabel=r'$\tau_{\rm LL}$')

    # Check MFP
    if (flg_test % 64) >= 32:
        fN_model = xifm.default_model()
        z = 2.44
        mfp = fN_model.mfp(z)
        print('MFP at z=%g is %g Mpc' % (z,mfp.value))
