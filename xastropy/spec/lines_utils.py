"""
#;+ 
#; NAME:
#; lines_utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for SpectralLine class
#;   03-Mar-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
from abc import ABCMeta, abstractmethod

from astropy import constants as const
from astropy import units as u
from astropy.units import Quantity

import xastropy.atomic as xatom
from xastropy.spec import readwrite as xsr
from xastropy.stats import basic as xsb
from xastropy.xutils import xdebug as xdb

# class SpectralLine(object):
# class AbsLine(SpectralLine):

# Class for Spectral line
class SpectralLine(object):
    """Class for a spectral line.  Emission or absorption 

    Attributes:
        ltype: string
          type of line, e.g.  Abs, Emiss
        wrest: float
          Rest wavelength of the spectral feature
    """
    __metaclass__ = ABCMeta

    # Initialize with wavelength
    def __init__(self, ltype, wrest):
        """  Initiator

        Parameters
        ----------
        ltype : string
          Type of Spectral line, (Abs
        wrest : float
          Rest wavelength
        """

        # Required
        self.ltype = ltype
        if ltype not in ['Abs']:
            raise ValueError('spec/lines: Not ready for type {:s}'.format(ltype))

        #if type(wrest) is not float:
        #    raise ValueError('spec/lines: Rest wavelength must be a float!')
        self.wrest = wrest

        # Other
        self.atomic = {} # Atomic Data
        self.analy = {} # Analysis inputs (e.g. spectrum; from .clm file or AbsID)
        self.attrib = {} # Properties (e.g. column, EW, centroid)

        # Fill atomic data
        self.fill_data()

    # Output
    def __repr__(self):
        txt = '[{:s}:'.format(self.__class__.__name__)
        try:
            txt = txt+' {:s},'.format(self.atomic['name'])
        except KeyError:
            pass
        txt = txt + ' wrest={:g}'.format(self.wrest)
        txt = txt + ']'
        return (txt)

# Class for Generic Absorption Line System
class AbsLine(SpectralLine):
    """Spectral absorption line
    """
    # Initialize with a .dat file
    def __init__(self, wrest): 
        # Generate with type
        SpectralLine.__init__(self,'Abs', wrest)

    def print_specline_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'AbsLine'

    # Fill atomic data and setup analy
    def fill_data(self):
        import xastropy.spec.abs_line as xspa

        # Data
        self.atomic = xspa.abs_line_data(self.wrest)

        #
        self.analy['WVMNX'] = [0., 0.] # Wavelength interval about the line (observed)
        self.analy['VLIM'] = [0., 0.]*u.km/u.s # Velocity limit of line
        self.analy['FLG_ANLY'] = 1 # Analyze
        self.analy['FLG_EYE'] = 0
        self.analy['FLG_LIMIT'] = 0 # No limit
        self.analy['DATFIL'] = '' 
        self.analy['IONNM'] = self.atomic['name']
        self.analy['z'] = 0.  # Redshift
        # Characteristics
        self.attrib = {'N': 0., 'Nsig': 0., 'flgN': 0, # Column
                       'b': 0., 'bsig': 0.,  # Doppler
                       'EW': 0., 'EWsig': 0., 'flgEW': 0 # EW
                       }

    # Perform AODM on the line
    def aodm(self, **kwargs): 
        """  AODM calculation

        Parameters
        ----------
        spec : Spectrum1D (None)
          1D spectrum.  Required but often read in through the Class (self.spec)
        conti : np.array (None)
          Continuum array 

        Returns:
          N, sigN : Column and error in linear space
        """

        # Grab Spectrum
        spec = self.set_spec(**kwargs)

        # Velocity array
        spec.velo = spec.relative_vel(self.wrest*(1+self.analy['z']))

        # Pixels for evaluation
        pix = spec.pix_minmax(self.analy['z'], self.wrest,
                        self.analy['VLIM'].to('km/s'))[0]
                        #self.analy['VLIM'].to('km/s').value)[0]

        # For convenience + normalize
        velo = spec.velo[pix]
        fx, sig = parse_spec(spec, **kwargs)

        # dv
        delv = velo - np.roll(velo,1)
        delv[0] = delv[1]

        # Atomic data
        cst = (10.**14.5761)/(self.atomic['fval']*self.wrest) / (u.km/u.s) / u.cm * (u.AA/u.cm)

        # Mask
        mask = (pix == pix) # True = good
        nndt = Quantity(np.zeros(len(pix)), unit='s/(km cm cm)')

        # Saturated?
        satp = np.where( (fx <= sig/5.) | (fx < 0.05) )[0]
        if len(satp) > 0:
            mask[satp] = False
            lim = np.where(sig[sat] > 0.)[0]
            if len(lim) > 0:
                sub = np.maximum(0.05, sig[sat[lim]]/5.)
                nndt[sat[lim]] = np.log(1./sub)*cst
                flg_sat = len(lim)
        # AODM
        nndt[mask] = np.log(1./fx[mask])*cst

        # Sum it
        ntot = np.sum( nndt*delv )
        tvar = np.sum( (delv*cst*sig/fx)**2 )

        # Fill
        self.attrib['N'] = ntot
        self.attrib['sigN'] = np.sqrt(tvar)
        logN, sig_logN = xsb.lin_to_log(self.attrib['N'].value, self.attrib['sigN'].value)
        self.attrib['logN'] = logN
        self.attrib['sig_logN'] = sig_logN

        # Return
        return ntot, np.sqrt(tvar)

    # EW 
    def ew(self, **kwargs):
        """  EW calculation
        Default is simple boxcar integration
        Observer frame, not rest-frame
          WVMNX must be set!

        Parameters
        ----------
        spec : Spectrum1D (None)
          1D spectrum.  Required but often read in through the Class (self.spec)
        conti : np.array (None)
          Continuum array 

        Returns:
          EW, sigEW : EW and error in observer frame
        """

        # Check on WVMNX
        if np.sum(self.analy['WVMNX']) == 0.:
            raise ValueError('lines_utils.ew: Need to set WVMNX!')

        # Grab spectrum
        spec = self.set_spec(**kwargs)

        # Pixels for evaluation
        pix = spec.pix_minmax(self.analy['WVMNX'])[0]

        # Normalized + convenience
        fx, sig = parse_spec(spec, **kwargs)
        wv = spec.dispersion[pix]

        # dwv
        dwv = wv - np.roll(wv,1)
        dwv[0] = dwv[1]

        # Units
        if spec.wcs.unit == 1.:
            raise ValueError('Expecting a unit!')

        # Simple boxcar
        EW = np.sum( dwv * (1. - fx) ) 
        varEW = np.sum( dwv**2 * sig**2 )
        sigEW = np.sqrt(varEW) 


        # Fill
        self.attrib['EW'] = EW 
        self.attrib['sigEW'] = sigEW 

        # Return
        return EW, sigEW
            
    # EW 
    def restew(self, **kwargs):
        """  Rest EW calculation
        Return rest-frame.  See "ew" above for details
        """
        # Standard call
        EW,sigEW = self.ew(**kwargs)
        # Push to rest-frame
        self.attrib['EW'] = EW / (self.analy['z']+1)
        self.attrib['sigEW'] = sigEW / (self.analy['z']+1)

        # Return
        return self.attrib['EW'], self.attrib['sigEW'] 

    # Check for a spectrum
    def set_spec(self, **kwargs):
        ''' Try to grab a spectrum for analysis
        '''
        
        try: # Internal?
            spec = self.spec
        except AttributeError:
            # Look for it as a keyword
            try:
                spec = kwargs['spec']
            except KeyError:
                raise IOError('lines_utils: Need a spectrum!')
        return spec

    # Output
    def __repr__(self):
        txt = '[{:s}:'.format(self.__class__.__name__)
        # Name
        try:
            txt = txt+' {:s},'.format(self.atomic['name'])
        except KeyError:
            pass
        # wrest
        txt = txt + ' wrest={:.4f}'.format(self.wrest)
        # fval
        try:
            txt = txt+', f={:g}'.format(self.atomic['fval'])
        except KeyError:
            pass
        txt = txt + ']'
        return (txt)


# ######
# Check for a spectrum
def parse_spec(spec, **kwargs):
    ''' Splice the spectrum.
    Normalize too
    '''
    fx = spec.flux[spec.sub_pix]
    sig = spec.sig[spec.sub_pix] 

    # Normalize?
    try:
        conti = kwargs['conti']
    except KeyError:
        pass
    else:
        if len(conti) != len(spec.flux): # Check length
            raise ValueError('lines_utils.aodm: Continuum length must match input spectrum')
        fx = fx / conti[spec.sub_pix]
        sig = sig / conti[spec.sub_pix]

    return fx, sig

## #################################    
## #################################    
## TESTING
## #################################    
if __name__ == '__main__':

    flg_test = 0
    #flg_test += 2**0  # AbsLine
    flg_test += 2**1  # AODM
    #flg_test += 2**2  # EW

    # Test Absorption Line creation
    if (flg_test % 2**1) >= 2**0:
        print('-------------------------')
        aline = AbsLine(1215.6701*u.AA)
        print(aline)

    # Test AODM
    if (flg_test % 2**2) >= 2**1:
        print('------------ AODM -------------')
        # Spectrum
        fil = '~/PROGETTI/LLSZ3/data/normalize/UM669_nF.fits'
        aline = AbsLine(1302.1685*u.AA)
        aline.spec = xsr.readspec(fil)
        # Line info
        aline.analy['z'] = 2.92652
        aline.analy['VLIM'] = const.c.to('km/s') * (
                    ( np.array([5110.668, 5116.305])*u.AA/(
                        1+aline.analy['z']) - aline.wrest) / aline.wrest )
        # Evaluate
        N,sigN = aline.aodm(conti=np.ones(len(aline.spec.flux)))
        logN, sig_logN = xsb.lin_to_log(N,sigN)
        print('logN = {:g}, sig_logN = {:g}'.format(logN, sig_logN))

    # Test EW
    if (flg_test % 2**3) >= 2**2:
        print('------------ EW -------------')
        # Spectrum
        fil = '~/PROGETTI/LLSZ3/data/normalize/UM669_nF.fits'
        aline = AbsLine(1302.1685)
        aline.spec = xsr.readspec(fil)
        # Line info
        aline.analy['z'] = 2.92652
        aline.analy['WVMNX'] = [5110.668, 5116.305]

        # Evaluate
        EW,sigEW = aline.restew()
        print('Rest EW = {:g}, sig = {:g}'.format(EW, sigEW))
