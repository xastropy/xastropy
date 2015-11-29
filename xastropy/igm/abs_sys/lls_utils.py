""" Subclasses for LLS AbsSystem and AbsSurvey
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import os, copy, imp, glob
import numpy as np
import warnings

from astropy import units as u
from astropy.table import Table, Column

from linetools.lists.linelist import LineList
from linetools.isgm.abssystem import AbsSystem, AbsSubSystem
from linetools.isgm import utils as ltiu

from xastropy.igm.abs_sys.abssurvey import AbslineSurvey
from xastropy.igm.abs_sys import abssys_utils as xabsu
from xastropy.igm.abs_sys import ionclms as xiai
from xastropy.atomic import ionization as xatomi
from xastropy.xutils import xdebug as xdb

xa_path = imp.find_module('xastropy')[1]

#class LLSSystem(AbslineSystem):
#class LLS_Survey(Absline_Survey):

# Class for LLS Absorption Lines 
class LLSSystem(AbsSystem):
    """
    Class for an LLS absorption system

    Attributes
    ----------
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

    @classmethod
    def from_datfile(cls, dat_file, tree=None, **kwargs):
        """ Read from dat_file (historical JXP format)

        Parameters
        ----------
        dat_file : str
          dat file
        tree : str, optional
          Path to data files
        kwargs :
          Passed to __init__
        """
        if tree is None:
            tree = ''
        # Read datfile
        datdict = xabsu.read_dat_file(tree+dat_file)
        # Parse
        coord, zabs, name, NHI, sigNHI, clm_fil = xabsu.parse_datdict(datdict)
        kwargs['NHI'] = NHI
        kwargs['sig_NHI'] = sigNHI
        # Generate with type
        vlim = None
        slf = cls(coord, zabs, vlim, **kwargs)

        # Fill files
        slf.tree = tree
        slf.dat_file = slf.tree+dat_file

        # Parse datdict
        #   Includes Sub systems
        slf._datdict = datdict
        slf.parse_dat_file()

        return slf

    def __init__(self, radec, zabs, vlim, **kwargs):
        """Standard init

        NHI keyword is required

        Parameters
        ----------
        radec : tuple or coordinate
            RA/Dec of the sightline or astropy.coordinate
        zabs : float
          Absorption redshift
        vlim : Quantity array (2)
          Velocity limits of the system
          Defaulted to +/- 500 km/s if None (see Prochaska et al. 2016 HDLLS)
        NHI= : float, required despite being a keyword
          log10 of HI column density
        **kwargs : keywords
          passed to AbsSystem.__init__

        """
        # NHI
        try:
            NHI = kwargs['NHI']
        except KeyError:
            raise ValueError("NHI must be specified for LLSSystem")
        else:
            kwargs.pop('NHI')
        # vlim
        if vlim is None:
            vlim = [-500.,500.]*u.km/u.s
        # Generate with type
        AbsSystem.__init__(self, 'LLS', radec, zabs, vlim, NHI=NHI, **kwargs)

        # Set tau_LL
        self.tau_LL = (10.**self.NHI)*6.3391597e-18 # Should replace with photocross

        # Other
        self.zpeak = None  # Optical depth weighted redshift
        self.MH = 0.

        # Subsystems
        self.nsub = 0
        self.subsys = {}

    """
    # Modify standard dat parsing
    def mk_subsys(self,nsub):
        '''Generate subsystems from Parent
        DEPRECATED
        Parameters:
        ----------
        nsub: int
          Number to generate
        '''
        assert False # DEPRECATED
        lbls= map(chr, range(65, 91))
        for i in range(self.nsub):
            self.subsys[lbls[i]] = Abs_Sub_System('LLS')
            self.subsys[lbls[i]].name = self.name+lbls[i]
            self.subsys[lbls[i]].coord = self.coord
            self.subsys[lbls[i]].zem = self.zem
            self.subsys[lbls[i]].linelist = self.linelist
    """

    # Modify standard dat parsing
    def parse_dat_file(self,vlim=[-300.,300]*u.km/u.s):
        """ Parse the datdict read from the .dat file

        Parameters
        ----------
        datdict : OrderedDict
          info from the .dat file
        vlim : Quantity array (2), optional
          Velocity limits of the subsystems
          Should be pulled from the .clm files
        """

        # LLS keys
        self.bgsrc = self._datdict['QSO name']
        self.zem = float(self._datdict['QSO zem'])  # Was zqso
        self.MH = float(self._datdict['[M/H] ave'])
        self.nsub = int(self._datdict['N subsys'])
        self.cldyfil = self._datdict['Cloudy Grid File']

        # LLS Subsystems
        if self.nsub > 0:
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
                zabs = float(self._datdict[lbls[i]+ ' zabs'])
                self.subsys[lbls[i]] = AbsSubSystem(self,zabs,vlim,lbls[i])
                self.subsys[lbls[i]]._datdict = {}
                # Fill in dict
                for ii,key in enumerate(keys):
                    try:
                        tmpc = self._datdict[lbls[i]+' '+key]
                    except:
                        raise ValueError('lls_utils: Key "{:s}" not found in {:s}'
                                         .format(lbls[i]+key,self.dat_file))
                    else:  # Convert
                        val = null_dict[key]
                        #pdb.set_trace()
                        if val.__class__ == np.ndarray:
                            self.subsys[lbls[i]]._datdict[att[ii]] = np.array(map(float,tmpc.split()))
                        else: # Single value
                            self.subsys[lbls[i]]._datdict[att[ii]] = (map(type(val),[tmpc]))[0]
                # Set a few special ones as attributes
                self.subsys[lbls[i]].NHI = self.subsys[lbls[i]]._datdict['NHI']
                self.subsys[lbls[i]].sig_NHI = self.subsys[lbls[i]]._datdict['NHIsig']

    # Fill up the ions
    def get_ions(self, use_clmfile=False, idict=None, closest=False, update_zvlim=True):
        """Parse the ions for each Subsystem

        And put them together for the full system
        Fills ._ionclms with a QTable

        Parameters
        ----------
        closest : bool, optional
          Take the closest line to input wavelength? [False]
        idict : dict, optional
          dict containing the IonClms info
        use_clmfile : bool, optional
          Parse ions from a .clm file (JXP historical)
        update_zvlim : bool, optional
          Update zvlim from lines in .clm (as applicable)
        """
        reload(ltiu)
        from linetools.abund import ions as ltai
        if idict is not None:
            # Manipulate for astropy Table
            #  Could probably use add_row or dict instantiation
            table = None
            for ion in idict.keys():
                Zion = ltai.name_ion(ion)
                if table is None:
                    tkeys = idict[ion].keys()
                    lst = [[idict[ion][tkey]] for tkey in tkeys]
                    table = Table(lst, names=tkeys)
                    # Extra columns
                    if 'Z' not in tkeys:
                        table.add_column(Column([Zion[0]],name='Z'))
                        table.add_column(Column([Zion[1]],name='ion'))
                else:
                    tdict = idict[ion]
                    if 'Z' not in tkeys:
                        tdict['Z'] = Zion[0]
                        tdict['ion'] = Zion[1]
                    # Add
                    table.add_row(tdict)
            # Finish
            try:  # Historical keys
                table.rename_column('clm','logN')
            except:
                pass
            else:
                table.rename_column('sig_clm','sig_logN')
                table.rename_column('flg_clm','flag_N')
            self._ionclms = table
        elif use_clmfile:
            # Subsystems
            if self.nsub > 0:  # This speeds things up (but is rarely used)
                linelist = LineList('ISM')
            for lbl in self.subsys.keys():
                clm_fil = self.tree+self.subsys[lbl]._datdict['clm_file']
                # Parse .clm file
                self.subsys[lbl]._clmdict = xiai.read_clmfile(clm_fil,linelist=linelist)
                # Build components from lines
                abslines = []
                vmin,vmax = 9999., -9999.
                for wrest in self.subsys[lbl]._clmdict['lines']:
                    vmin = min(vmin, self.subsys[lbl]._clmdict['lines'][wrest].analy['vlim'][0].value)
                    vmax = max(vmax, self.subsys[lbl]._clmdict['lines'][wrest].analy['vlim'][1].value)
                    self.subsys[lbl]._clmdict['lines'][wrest].attrib['coord'] = self.coord
                    abslines.append(self.subsys[lbl]._clmdict['lines'][wrest])
                components = ltiu.build_components_from_abslines(abslines)
                # Update z, vlim
                if update_zvlim:
                    self.subsys[lbl].zabs = self.subsys[lbl]._clmdict['zsys']
                    self.subsys[lbl].vlim = [vmin,vmax]*u.km/u.s
                # Read .ion file and fill in components
                ion_fil = self.tree+self.subsys[lbl]._clmdict['ion_fil']
                self.subsys[lbl]._indiv_ionclms = xiai.read_ion_file(ion_fil,components=components)
                # Parse .all file
                all_file = ion_fil.split('.ion')[0]+'.all'
                self.subsys[lbl].all_file=all_file #MF: useful to have
                _ = xiai.read_all_file(all_file,components=components)
                # Build table
                self.subsys[lbl]._ionclms = ltiu.iontable_from_components(components,ztbl=self.subsys[lbl].zabs)
                # Add to AbsSystem
                for comp in components:
                    self.add_component(comp)

            # Combine
            if self.nsub == 1:
                self._ionclms = self.subsys['A']._ionclms
                self._clmdict = self.subsys['A']._clmdict
                #xdb.set_trace()
            elif self.nsub == 0:
                raise ValueError('lls_utils.get_ions: Cannot have 0 subsystems..')
            else:
                self._ionclms = self.subsys['A']._ionclms
                self._clmdict = self.subsys['A']._clmdict
                warnings.warn('lls_utils.get_ions: Need to update multiple subsystems!! Taking A.')
        else:
            raise ValueError("Need an option in get_ions")


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
            aline.attrib['N'] = 10**self.NHI / u.cm**2
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
        from linetools.analysis import voigt as lav
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
        return ('[{:s}: {:s} {:s}, zabs={:g}, NHI={:g}, tau_LL={:g}, [M/H]={:g} dex]'.format(
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
    """
    An LLS Survey class

    Attributes:
        
    """
    @classmethod
    def load_HDLLS(cls):
        """ Default sample of LLS (HD-LLS, DR1)
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
        """
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
        lls_survey = cls.from_sfits(summ_fil)
        # Load ions
        lls_survey.fill_ions(jfile=ions_fil)
        # Set data path (may be None)
        for lls in lls_survey._abs_sys:
            lls.spec_path = os.getenv('HDLLS_DATA')

        return lls_survey

    """
    @classmethod
    def from_flist(cls, flist, tree=None, **kwargs):
        xdb.set_trace()
        slf = AbslineSurvey.from_flist(self,'LLS', flist, tree=tree, **kwargs)
        return slf
    """

    def __init__(self, **kwargs):
        # Generate with type
        AbslineSurvey.__init__(self, 'LLS', **kwargs)

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

def tau_multi_lls(wave, all_lls, **kwargs):
    '''Calculate opacities on an input observed wavelength grid
    Parameters:
    -----------
    wave : Quantity array
      Wavelengths
    all_lls : List of LLS Class
    **kwargs : extra keywords go to lav.voigt_from_abslines

    Returns:
    --------
    tau : ndarray
      Optical depth values at input wavelengths
    '''
    from xastropy.atomic import ionization as xai
    from linetools.analysis import voigt as lav
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
        tau_Lyman = lav.voigt_from_abslines(wave, lls.lls_lines, ret='tau', **kwargs)
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
