"""
#;+ 
#; NAME:
#; absline
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for absorption line kinematics
#;   29-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp
from astropy.io import fits, ascii
from astropy import units as u 
from astropy.convolution import convolve, Box1DKernel
#from astropy import constants as const

from xastropy.xutils import xdebug as xdb
from xastropy import spec as xspec

#def ion_name(ion):
#def photo_cross(Z, ion, E, datfil=None, silent=False):

########################## ##########################
########################## ##########################
class Kin_Abs(object):
    """ Class for kinematics on an absorption line

    Attributes
    ----------
    wrest: float
      Rest wavelength of line analyzed
    vmnx: tuple (vmin,vmax)
      Velocity range for analysis

    JXP on 11 Dec 2014
    """

    # Initialize 
    def __init__(self, wrest, vmnx):

        # Absorption system
        self.wrest = wrest
        self.vmnx = vmnx

        # Data
        self.kin_data = {}
        self.keys = ['flg', 'Dv', 'fedg', 'fmm']
        self.key_dtype = ['i4', 'f4', 'f4', 'f4']

        # Init
        for key in self.keys:
            self.kin_data[key] = 0
        #xdb.set_trace()

    # Access the data and return the value
    def __getitem__(self, item):
        try:
            return self.kin_data[item] 
        except KeyError:
           raise KeyError 

    # Output
    def __repr__(self):
        return ('[{:s}: {:g}]'.format(
                self.__class__.__name__, self.wrest) )

    ########################## ##########################
    def orig_kin(self, spec, kbin=22., per=0.05, fill=True, debug=False):
        """ Measure a standard suite of absorption line kinematics
    
        Parameters
        ----------
        spec: Spectrum1D class
          Input spectrum
        velo is expected to have been filled already
        fill: bool (True)
          Fill the dictionary with some items that other kin programs may need

        Returns
        -------
        out_kin : dict
        Dictionary of kinematic measurements
    
        JXP on 21 Nov 2014
        """

        # Calcualte dv
        imn = np.argmin( np.fabs(spec.velo) )
        dv = np.abs( spec.velo[imn] - spec.velo[imn+1] )

        # Test for bad pixels
        pixmin = np.argmin( np.fabs( spec.velo-self.vmnx[0] ) )
        pixmax = np.argmin( np.fabs( spec.velo-self.vmnx[1] ) )
        pix = np.arange(pixmin, pixmax+1)
        npix = len(pix)
        badzero=np.where((spec.flux[pix] == 0) & (spec.sig[pix] <= 0))[0]
        if len(badzero) > 0:
            if np.max(badzero)-np.min(badzero) >= 5: 
                raise ValueError('orig_kin: too many or too large sections of bad data')
            
            spec.flux[pix[badzero]] = np.mean(np.array([spec.flux[pix[np.min(badzero)-1]],
                                                        spec.flux[pix[np.max(badzero)+1]]]))
            xdb.set_trace() # Should add sig too

        # Generate the tau array
        tau = np.zeros(npix)
        gd = np.where((spec.flux[pix] > spec.sig[pix]/2.) &
                    (spec.sig[pix] > 0.) )
        if len(gd) == 0:
            raise ValueError('orig_kin: Profile too saturated.')

        tau[gd] = np.log(1./spec.flux[pix[gd]])
        sat = (pix == pix)
        sat[gd] = False
        tau[sat] = np.log(2./spec.sig[pix[sat]])

        # Smooth
        nbin = np.round(kbin/dv)
        kernel = Box1DKernel(nbin, mode='center')
        stau = convolve(tau, kernel, boundary='fill', fill_value=0.)
        if debug is True:
            xdb.xplot(spec.velo[pix], tau, stau)

        # Dv (usually dv90)
        tottau = np.sum( stau )
        cumtau = np.cumsum(stau) / tottau
        lft = (np.where(cumtau > per)[0])[0]
        rgt = (np.where(cumtau > (1.-per))[0])[0] - 1
        self.kin_data['Dv'] = np.round(np.abs(spec.velo[pix[rgt]]-spec.velo[pix[lft]]))
        #xdb.set_trace()

        # Mean/Median
        vcen = (spec.velo[pix[rgt]]+spec.velo[pix[lft]])/2.
        mean = self.kin_data['Dv']/2.
        imn = np.argmin( np.fabs(cumtau-0.5) )
        self.kin_data['fmm'] = np.abs( (spec.velo[pix[imn]]-vcen)/mean )
    
        # fedg
        imx = np.argmax(stau)
        self.kin_data['fedg'] = np.abs( (spec.velo[pix[imx]]-vcen) / mean )
    
        # Two-peak :: Not ported..  Not even to XIDL!

        # Set flag
        self.kin_data['flg'] = 1


        if fill is True:
            self.stau = stau
            self.pix = pix




    
#### ###############################
#### ###############################
# Testing
if __name__ == '__main__':

    flg_test = 0
    flg_test = 1  # First test of orig_kin

    # First test
    if (flg_test % 2**1) >= 2**0:
        # Grab spectrum
        spec_fil = '/u/xavier/Keck/HIRES/RedData/PH957/PH957_f.fits'
        spec = xspec.readwrite.readspec(spec_fil)
        vmnx = (-25., 80.)
        wrest = 1741.5490
        zabs = 2.309
        # Generate velo
        spec.velo = spec.relative_vel( (1+zabs)*wrest )
        kin = Kin_Abs(wrest, vmnx)
        # Call kin
        kin.orig_kin(spec)
        print('Kin results = {:g}, {:g}, {:g}'.format(kin['Dv'],
                                                      kin['fmm'], 
                                                      kin['fedg'] ))
