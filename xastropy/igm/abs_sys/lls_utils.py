"""
#;+ 
#; NAME:
#; lls_utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for Lyman Limit Systems
#;   27-Oct-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import os, copy, sys, imp, glob
import numpy as np

from astropy import units as u
from astropy.io import ascii 

from linetools.analysis import voigt as lav
from linetools.lists.linelist import LineList

from xastropy.igm.abs_sys.abs_survey import AbslineSurvey
from xastropy.igm.abs_sys.abssys_utils import AbslineSystem, Abs_Sub_System
from xastropy.igm.abs_sys.ionclms import Ionic_Clm_File, IonClms
from xastropy.atomic import ionization as xatomi
from xastropy.xutils import xdebug as xdb

xa_path = imp.find_module('xastropy')[1]

#class LLSSystem(AbslineSystem):
#class LLS_Survey(Absline_Survey):

# Class for LLS Absorption Lines 
class LLSSystem(AbslineSystem):
    """An LLS absorption system

    Attributes:
        tau_ll: Opacity at the Lyman limit
    """
    # Create instance with a AbsID file
    @classmethod
    def from_absid_fil(cls, abs_fil, linelist=None):
        lls = cls() # Empty
        # Linelist (speeds things up)
        lls.linelist=linelist
        # Parse abs_fil
        lls.parse_absid_file(abs_fil)
        # Return
        return lls

    def __init__(self, dat_file=None, tree=None, **kwargs):
        # Generate with type
        AbslineSystem.__init__(self,'LLS', **kwargs)
        # Over-ride tree?
        if tree != None:
            self.tree = tree
        else:
            self.tree = ''
        # Parse .dat file
        if dat_file != None:
            self.parse_dat_file(self.tree+dat_file)
            self.dat_file = self.tree+dat_file

        # Set tau_LL
        self.tau_LL = (10.**self.NHI)*6.3391597e-18 # Should replace with photocross
        self.ions = None
        self.zpeak = None


    # Modify standard dat parsing
    def mk_subsys(self,nsub):
        '''Generate subsystems from Parent
        Parameters:
        ----------
        nsub: int
          Number to generate
        '''
        lbls= map(chr, range(65, 91))
        self.nsub = nsub
        self.subsys = {}
        for i in range(self.nsub):
            self.subsys[lbls[i]] = Abs_Sub_System('LLS')
            self.subsys[lbls[i]].name = self.name+lbls[i]
            self.subsys[lbls[i]].coord = self.coord
            self.subsys[lbls[i]].zem = self.zem
            self.subsys[lbls[i]].linelist = self.linelist

    # Modify standard dat parsing
    def parse_dat_file(self,dat_file):
        # Standard Call
        out_list = AbslineSystem.parse_dat_file(self,dat_file,flg_out=1)

        # LLS keys
        self.bgsrc = self.datdict['QSO name']
        self.zem = float(self.datdict['QSO zem'])  # Was zqso
        self.MH = float(self.datdict['[M/H] ave'])
        self.nsub = int(self.datdict['N subsys'])
        self.cldyfil = self.datdict['Cloudy Grid File']

        # LLS Subsystems
        if self.nsub > 0:
            self.subsys = {}
            lbls= map(chr, range(65, 91))
            # Dict
            keys = (['zabs','NHI','NHIsig','NH','NHsig','log x','sigx','b','bsig','Abund file',
                     'U','Usig','flg_low','flg_alpha','[alpha/H]','sig[a/H]',
                     'flg_Fe','[Fe/H]','sig[Fe/H]','VPFIT file'])
            att = (['zabs','NHI','NHIsig','NH','NHsig','logx','sigx','bval','bsig','clm_file',
                     'U','Usig','flg_low','flg_alpha','alpha_H','sig_a_H',
                     'flg_Fe','Fe_H','sig_Fe_H','VPFIT_file'])
            values = ([0., 0., np.zeros(2), 0., np.zeros(2), 0., np.zeros(2), 0., 0.,
                    '', 0., np.zeros(2), 0, 0, 0., 0., 0, 0., 0., ''])
            null_dict = dict(zip(keys,values))
            # Loop on subsystems
            for i in range(self.nsub):
                # Generate
                self.subsys[lbls[i]] = Abs_Sub_System('LLS')
                self.subsys[lbls[i]].name = self.name+lbls[i]
                self.subsys[lbls[i]].coord = self.coord
                self.subsys[lbls[i]].tree = self.tree
                self.subsys[lbls[i]].linelist = self.linelist
                # Fill in
                for ii,key in enumerate(keys):
                    try:
                        tmpc = self.datdict[lbls[i]+' '+key]
                    except:
                        raise ValueError('lls_utils: Key "{:s}" not found in {:s}'
                                         .format(lbls[i]+key,dat_file))
                    else:  # Convert
                        val = null_dict[key]
                        #pdb.set_trace()
                        if val.__class__ == np.ndarray:  
                            setattr(self.subsys[lbls[i]], att[ii],
                                    np.array(map(float,tmpc.split())) )
                        else: # Single value
                            setattr(self.subsys[lbls[i]], att[ii], (map(type(val),[tmpc]))[0] )

    # Fill up the ions
    def get_ions(self, idict=None, closest=False):
        """Parse the ions for each Subsystem
        And put them together for the full system
        Fills .ions with a IonsClm Class

        Parameters:
        -----------
        closest : bool, optional
          Take the closest line to input wavelength? [False]
        idict : dict, optional
          dict containing the IonClms info
        """
        if idict is not None:
            # Not setup for SubSystems
            self._ionclms = IonClms(idict=idict)
        else:
            # Subsystems
            if self.nsub > 0:  # This speeds things up (but is rarely used)
                if self.linelist is None:
                    self.linelist = LineList('ISM')
            lbls= map(chr, range(65, 91))
            for ii in range(self.nsub):
                clm_fil = self.tree+self.subsys[lbls[ii]].clm_file
                # Parse .clm and .all files
                self.subsys[lbls[ii]].clm_analy = Ionic_Clm_File(clm_fil, self.linelist)
                ion_fil = self.tree+self.subsys[lbls[ii]].clm_analy.ion_fil 
                all_fil = ion_fil.split('.ion')[0]+'.all'
                self.all_fil=all_fil #MF: useful to have
                self.subsys[lbls[ii]]._ionclms = IonClms(all_fil)
                # Linelist (for speed)
                if self.subsys[lbls[ii]].linelist is None:
                    self.subsys[lbls[ii]].linelist = self.linelist
                self.subsys[lbls[ii]].linelist.closest = closest
                # Parse .ion file
                self.subsys[lbls[ii]].read_ion_file(ion_fil)

            # Combine
            if self.nsub == 1:
                self._ionclms = self.subsys['A']._ionclms
                self.clm_analy = self.subsys['A'].clm_analy
                self.lines = self.subsys['A'].lines
                #xdb.set_trace()
            elif self.nsub == 0:
                raise ValueError('lls_utils.get_ions: Cannot have 0 subsystems..')
            else:
                self._ionclms = self.subsys['A']._ionclms
                self.clm_analy = self.subsys['A'].clm_analy
                self.lines = self.subsys['A'].lines
                print('lls_utils.get_ions: Need to update multiple subsystems!! Taking A.')
                

    # #############
    def fill_lls_lines(self, bval=20.*u.km/u.s):
        """
        Generate an HI line list for an LLS.
        Goes into self.lls_lines 

        Parameters
        ----------
        bval : float (20.)  Doppler parameter in km/s
        """
        from linetools.lists import linelist as lll
        from linetools.spectralline import AbsLine

        # May be replaced by component class (as NT desires)
        HIlines = lll.LineList('HI')

        self.lls_lines = []
        for lline in HIlines._data:
            aline = AbsLine(lline['wrest'],linelist=HIlines)
            # Attributes
            aline.attrib['N'] = self.NHI
            aline.attrib['b'] = bval
            aline.attrib['z'] = self.zabs
            # Could set RA and DEC too
            aline.attrib['RA'] = self.coord.ra
            aline.attrib['DEC'] = self.coord.dec
            self.lls_lines.append(aline)

    # Absorption model of the LLS (HI only)
    def flux_model(self,spec,smooth=0):
        """
        Generate a LLS model given an input spectrum

        Parameters:
          spec:  Spectrum1D
          smooth : (0) Number of pixels to smooth by

        Returns:
          Output model is passed back as a Spectrum 
        """
        # ########
        # LLS first

        # Energies in LLS rest-frame
        wv_rest = spec.dispersion / (self.zabs+1)
        energy = wv_rest.to(u.eV, equivalencies=u.spectral())

        # Get photo_cross and calculate tau
        tau_LL = (10.**self.NHI / u.cm**2) * xatomi.photo_cross(1,1,energy)

        # ########
        # Now the Lyman series

        # Check for lines
        if 'lls_lines' not in self.__dict__.keys():
            self.fill_lls_lines()

        #xdb.set_trace()
        tau_Lyman = lav.voigt_from_abslines(spec.dispersion, self.lls_lines, ret='tau')

        # Combine
        tau_model = tau_LL + tau_Lyman

        # Kludge around the limit
        pix_LL = np.argmin( np.fabs( wv_rest- 911.3*u.AA ) )
        pix_kludge = np.where( (wv_rest > 911.5*u.AA) & (wv_rest < 912.8*u.AA) )[0]
        tau_model[pix_kludge] = tau_model[pix_LL]
        
        # Fill in flux
        model = copy.deepcopy(spec)
        model.flux = np.exp(-1. * tau_model).value

        # Smooth?
        if smooth > 0:
            model.gauss_smooth(smooth)

        # Return
        return model

    # 
    def get_zpeak(self):
        ''' Measure zpeak from an ionic transition
        '''
        if self.ions is None:
            print('get_zpeak: Need to fill ions with get_ions first.')
            return

        # Ions for analysis
        low_ions = [ (14,2), (6,2), (13,2), (26,2), (13,3)]
        high_ions= [(14,4), (6,4)]

        for tt in range(4):
            if tt == 0:
                ions = low_ions
                iflg = 1 # Standard
            elif tt == 1:
                ions = low_ions
                iflg = 2 # Saturated
            elif tt == 2:
                ions = high_ions
                iflg = 1 # Standard
            elif tt == 3:
                ions = high_ions
                iflg = 2 # Standard
            else:
                raise ValueError('Bad value')

            # Search 
            for ion in ions:
                try:
                    t = self.ions[ion]
                except KeyError:
                    continue
                # Measurement?
                if t['flg_clm'] == iflg:
                # Identify the transition
                    gdi = np.where( (self.ions.trans['Z'] == ion[0]) &
                                (self.ions.trans['ion'] == ion[1]) &
                                (self.ions.trans['flg_clm'] <= iflg) )[0]
                    # Take the first one
                    gdt = self.ions.trans[gdi[0]]
                    wrest = gdt['wrest']
                    flgs = self.clm_analy.clm_lines[wrest].analy['FLAGS']
                    spec_file = self.clm_analy.fits_files[flgs[1] % 64]
                    # Generate an Abs_Line with spectrum
                    line = abs_line.Abs_Line(wrest, z=self.clm_analy.zsys, spec_file=spec_file)
                    # vpeak
                    from astropy.relativity import velocities as arv
                    vpeak = line.vpeak()
                    self.zpeak = arv.z_from_v(self.clm_analy.zsys, vpeak)
                    if tt == 3:
                        print('zpeak WARNING: Using saturated high-ions!!')
                    break
            else:
                continue
            # get out
            break

        # Error catching
        if self.zpeak is None:
            # Skip primordial LLS
            print('lls_utils.zpeak: No transition in {:s}'.format(self.clm_analy.clm_fil))
            xdb.set_trace()
            return (0,0), 0.
        # Return
        return ion, vpeak

    # Output
    def __repr__(self):
        return ('[{:s}: {:s} {:s}, z={:g}, NHI={:g}, tau={:g}, [M/H]={:g} dex]'.format(
                self.__class__.__name__,
                 self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True,alwayssign=True),
                 self.zabs, self.NHI, self.tau_LL, self.MH))

    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'LLS'

# #######################################################################
# #######################################################################
# #######################################################################
# Class for LLS Survey
class LLSSurvey(AbslineSurvey):
    """An LLS Survey class

    Attributes:
        
    """
    # Initialize with a .dat file
    def __init__(self, **kwargs): 
        # Generate with type
        AbslineSurvey.__init__(self, 'LLS', **kwargs)

    # Cut on NHI
    def cut_nhi_quality(self, sig_cut=0.4):
        """
        Cut the LLS on NHI quality.
        Could put this in Absline_Survey

        Parameters:
          sig_cut: float (0.4) 
            Limit to include as quality

        Returns:
          gdNHI, bdNHI
            Indices for those LLS that are good/bad
            gdNHI is a numpy array of indices
            bdNHI is a boolean array
        """
        # Cut
        gdNHI = np.where( (self.sigNHI[:,0] < sig_cut)
                        & (self.sigNHI[:,1] < sig_cut))[0] 
        # Mask
        bdNHI = (self.NHI == self.NHI)
        bdNHI[gdNHI] = False

        # Return
        return gdNHI, bdNHI

    # Default sample of LLS (HD-LLS, DR1)
    @classmethod
    def default_sample(cls):
        '''
        # Local
        sys.path.append(os.path.abspath(os.environ.get('LLSPAP')+
                                        "/Optical/Data/Analysis/py"))
        import lls_sample as lls_s

        lls = cls('Lists/lls_metals.lst', tree=os.environ.get('LLSTREE'))
        # Mask
        msk = lls_s.hdlls(lls)
        lls.upd_mask(msk)

        # Return
        return lls
        '''
        import urllib2
        # Pull from Internet (as necessary)
        summ_fil = glob.glob(xa_path+"/data/LLS/HD-LLS_DR1.fits")
        if len(summ_fil) > 0:
            summ_fil = summ_fil[0]
        else:
            url = 'http://www.ucolick.org/~xavier/HD-LLS/DR1/HD-LLS_DR1.fits'
            print('LLSSurvey: Grabbing summary file from {:s}'.format(url))
            f = urllib2.urlopen(url)
            summ_fil = xa_path+"/data/LLS/HD-LLS_DR1.fits"
            with open(summ_fil, "wb") as code:
                code.write(f.read())

        # Ions
        ions_fil = glob.glob(xa_path+"/data/LLS/HD-LLS_ions.json")
        if len(ions_fil) > 0:
            ions_fil = ions_fil[0]
        else:
            url = 'http://www.ucolick.org/~xavier/HD-LLS/DR1/HD-LLS_ions.json'
            print('LLSSurvey: Grabbing JSON ion file from {:s}'.format(url))
            f = urllib2.urlopen(url)
            ions_fil = xa_path+"/data/LLS/HD-LLS_ions.json"
            with open(ions_fil, "wb") as code:
                code.write(f.read())

        # Read
        lls = cls(summ_fits=summ_fil)
        lls.fill_ions(jfile=ions_fil)

        return lls

def tau_multi_lls(wave, all_lls):
    '''Calculate opacities on an input observed wavelength grid
    Parameters:
    -----------
    wave: Quantity array
      Wavelengths
    all_lls: List of LLS Class
    '''
    from xastropy.atomic import ionization as xai
    #
    all_tau_model = np.zeros(len(wave))
    # Loop on LLS
    for lls in all_lls:
        # LL
        wv_rest = wave / (lls.zabs+1)
        energy = wv_rest.to(u.eV, equivalencies=u.spectral())
        # Get photo_cross and calculate tau
        tau_LL = (10.**lls.NHI / u.cm**2) * xai.photo_cross(1,1,energy)

        # Lyman
        tau_Lyman = lav.voigt_from_abslines(wave, lls.lls_lines, ret='tau')
        tau_model = tau_LL + tau_Lyman

        # Kludge around the limit
        pix_LL = np.argmin( np.fabs( wv_rest- 911.3*u.AA ) )
        pix_kludge = np.where( (wv_rest > 911.5*u.AA) & (wv_rest < 912.8*u.AA) )[0]
        tau_model[pix_kludge] = tau_model[pix_LL]
        # Add
        all_tau_model += tau_model
    # Return
    return all_tau_model

#
def profile_read():
    '''
    Bit for profiling the main read commands
    '''
    # Tree
    lls = LLS_Survey('Lists/lls_metals.lst', tree=os.getenv('LLSTREE'))
    # Mask
    # Ions
    lls.fill_ions()
    return lls


## #################################    
## #################################    
## TESTING
## #################################    
if __name__ == '__main__':

    flg_test = 0
    #flg_test = 1  # ions
    #flg_test += 2 # LLS plot
    #flg_test += 2**2 # zpeak
    #flg_test += 2**3 # output .dat file
    flg_test += 2**4 # read/write AbsID
    #
    #flg_test += 2**9 # LLS Survey NHI
    #flg_test += 2**10 # LLS Survey ions
    #flg_test += 2**11 # Profiling the main read commands

    # Test Absorption System
    print('-------------------------')
    tmp1 = LLSSystem(dat_file='Data/HE0940-1050.z2916.dat',
                      tree=os.environ.get('LLSTREE'))
    print(tmp1)

    # Test ions
    if (flg_test % 2**1) >= 2**0:
        print('-------------------------')
        tmp1.get_ions()
        print('C IV: ')
        print(tmp1.ions[(6,4)])
        print(tmp1.ions.trans[5]) # CIV 1550
    
    # Plot the LLS
    if (flg_test % 2**2) >= 2**1:
        print('-------------------------')
        tmp1.fill_lls_lines()

        from specutils.spectrum1d import Spectrum1D
        wa=np.linspace(3400.0,5000.0,10000)
        spec = Spectrum1D.from_array(wa,np.ones(len(wa)))

        model = tmp1.flux_model(spec, smooth=4)
        model.qck_plot()
        #xdb.set_trace()

    # Test zpeak
    if (flg_test % 2**3) >= 2**2:
        print('-------------------------')
        tmp1.fill_ions()
        ion,vpeak = tmp1.get_zpeak()
        print('zpeak = {:g}'.format(tmp1.zpeak))

    # Write .dat
    if (flg_test % 2**4) >= 2**3:
        tmp1.write_dat_file()

    # Read and Write AbsID
    if (flg_test % 2**5) >= 2**4:
        abs_fil = '/Users/xavier/paper/LLS/Optical/Data/Analysis/MAGE/SDSSJ1004+0018_z2.746_id.fits'
        lls = LLSSystem.from_absid_fil(abs_fil)
        tmpfil= '/Users/xavier/Desktop/tmp.fits'
        lls.write_absid_file(tmpfil)
        lls = LLSSystem.from_absid_fil(tmpfil)
        xdb.set_trace()

    # #############################
    # LLS Survey
    if (flg_test % 2**10) >= 2**9:
        print('-------------------------')
        lls = LLSSurvey('Lists/lls_metals.lst', tree=os.environ.get('LLSTREE'))
        xdb.xhist(lls.NHI, binsz=0.30)

    # LLS Survey ions
    if (flg_test % 2**11) >= 2**10:
        lls = LLSSurvey('Lists/lls_metals.lst', tree=os.environ.get('LLSTREE'))
        lls.fill_ions()
        xdb.xhist(lls.ions((6,4),skip_null=True)['clm'], binsz=0.3,
                  xlabel=r'$\log_{10} N({\rm C}^{+3})$')
        xdb.xhist(lls.ions((14,2),skip_null=True)['clm'], binsz=0.3,
                  xlabel=r'$\log_{10} N({\rm Si}^{+})$')

    # Profiling the Read
    if (flg_test % 2**12) >= 2**11:
        # python -m cProfile -o lls_profile.dat lls_utils.py
        import time
        t0 = time.clock()
        lls = profile_read()
        t1 = time.clock()
        print('Total time = {:g}'.format(t1-t0))

    # All done
    print('-------------------------')
    print('lls_utils: All done testing..')
