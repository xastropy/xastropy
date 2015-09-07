"""
#;+
#; NAME:
#; continuum 
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for continuum code 
#;   20-Aug-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp
import astropy as apy

from astropy import units as u
from astropy import constants as const
from astropy.io import fits, ascii

from linetools.spectra.xspectrum1d import XSpectrum1D

from xastropy.xutils import xdebug as xdb

xa_path = imp.find_module('xastropy')[1]

def init_conti_dict(Norm=0., tilt=0., piv_wv=0., igm='None',
    fN_gamma=-1., LL_flatten='True'):
    '''Initialize a continuum conti_dict
    Parameters:
    ----------
    Norm: float, optional
      Normaliztion
    tilt: float, optional
      Power-law tilt to continuum
    piv_wv: float, optional
      Pivot wave for tilt.  Best kept *without* units
    igm: str, optional
      Adopt average IGM model? ['None']
    LL_flatten: bool, optional
      Set Telfer to a constant below the LL?

    Returns:
    ---------
    conti_dict: dict 
      Useful for simple modeling.  Keep as a dict for JSON writing
    '''
    conti_dict = dict(Norm=Norm, tilt=tilt, piv_wv=piv_wv, igm=igm,
        fN_gamma=fN_gamma, LL_flatten=LL_flatten)
    #
    return conti_dict

def get_telfer_spec(zqso=0., igm=False, fN_gamma=None, LL_flatten=True):
    '''Generate a Telfer QSO composite spectrum

    Paraemters:
    ----------
    zqso: float, optional
      Redshift of the QSO
    igm: bool, optional
      Include IGM opacity? [False]
    fN_gamma: float, optional
      Power-law evolution in f(N,X)
    LL_flatten: bool, optional
      Set Telfer to a constant below the LL?

    Returns:
    --------
    telfer_spec: XSpectrum1D
      Spectrum
    '''
    # Read
    telfer = ascii.read(
        xa_path+'/data/quasar/telfer_hst_comp01_rq.ascii', comment='#')
    scale = telfer['flux'][(telfer['wrest'] == 1450.)]
    telfer_spec = XSpectrum1D.from_tuple((telfer['wrest']*(1+zqso),
        telfer['flux']/scale[0])) # Observer frame

    # IGM?
    if igm is True:
        '''The following is quite experimental.
        Use at your own risk.
        '''
        import multiprocessing
        from xastropy.igm.fN import model as xifm
        from xastropy.igm import tau_eff as xit
        fN_model = xifm.default_model()
        # Expanding range of zmnx (risky)
        fN_model.zmnx = (0.,5.)
        if fN_gamma is not None:
            fN_model.gamma = fN_gamma
        # Parallel
        igm_wv = np.where(telfer['wrest']<1220.)[0]
        adict = []
        for wrest in telfer_spec.dispersion[igm_wv].value:
            tdict = dict(ilambda=wrest, zem=zqso, fN_model=fN_model)
            adict.append(tdict)
        # Run
        #xdb.set_trace()
        pool = multiprocessing.Pool(4) # initialize thread pool N threads
        ateff = pool.map(xit.map_etl, adict)
        # Apply
        telfer_spec.flux[igm_wv] *= np.exp(-1.*np.array(ateff))
        # Flatten?
        if LL_flatten:
            wv_LL = np.where(np.abs(telfer_spec.dispersion/(1+zqso)-914.*u.AA)<3.*u.AA)[0]
            f_LL = np.median(telfer_spec.flux[wv_LL])
            wv_low = np.where(telfer_spec.dispersion/(1+zqso)<911.7*u.AA)[0]
            telfer_spec.flux[wv_low] = f_LL

    # Return
    return telfer_spec
 
## #################################    
## #################################    
## TESTING
## #################################    
if __name__ == '__main__':

    flg_tst = 0 
    flg_tst += 2**0  # Simple Telfer

    #if (flg_fig % 2**4) >= 2**3:
