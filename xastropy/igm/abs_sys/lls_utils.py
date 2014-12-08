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

import os, copy, sys
import numpy as np
import yaml

from astropy import units as u
from astropy.io import ascii 

from xastropy.igm.abs_sys.abssys_utils import Absline_System, Abs_Sub_System
from xastropy.igm.abs_sys.abs_survey import Absline_Survey
from xastropy.igm.abs_sys.ionic_clm import Ionic_Clm_File, Ions_Clm
from xastropy.spec import abs_line, voigt
from xastropy.atomic import ionization as xatomi
from xastropy.xutils import xdebug as xdb

#class LLS_System(Absline_System):
#class LLS_Survey(Absline_Survey):

# Class for LLS Absorption Lines 
class LLS_System(Absline_System):
    """An LLS absorption system

    Attributes:
        tau_ll: Opacity at the Lyman limit
    """
    # Initialize with a .dat file
    def __init__(self, dat_file=None, tree=None):
        # Generate with type
        Absline_System.__init__(self,'LLS')
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

        # Name

    # Modify standard dat parsing
    def parse_dat_file(self,dat_file):
        # Standard Call
        out_list = Absline_System.parse_dat_file(self,dat_file,flg_out=1)

        # LLS keys
        self.qso = self.datdict['QSO name']
        self.zqso = float(self.datdict['QSO zem'])
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
                self.subsys[lbls[i]].name = self.name
                self.subsys[lbls[i]].coord = self.coord
                self.subsys[lbls[i]].tree = self.tree
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
                            setattr(self.subsys[lbls[i]], att[ii], np.array(map(float,tmpc.split())) )
                        else: # Single value
                            setattr(self.subsys[lbls[i]], att[ii], (map(type(val),[tmpc]))[0] )

    # Fill up the ions
    def get_ions(self):
        """
        Parse the ions for each Subsystem
        And put them together for the full system

        Fills .ions with a Ions_Clm Class
        """
        # Subsystems
        lbls= map(chr, range(65, 91))
        for ii in range(self.nsub):
            clm_fil = self.tree+self.subsys[lbls[ii]].clm_file
            # Parse .clm and .all files
            self.subsys[lbls[ii]].clm_analy = Ionic_Clm_File(clm_fil)
            ion_fil = self.tree+self.subsys[lbls[ii]].clm_analy.ion_fil 
            all_fil = ion_fil.split('.ion')[0]+'.all'
            self.subsys[lbls[ii]].ions = Ions_Clm(all_fil, trans_file=ion_fil)

        # Combine
        if self.nsub == 1:
            self.ions = self.subsys['A'].ions
            self.clm_analy = self.subsys['A'].clm_analy
            #xdb.set_trace()
        elif self.nsub == 0:
            raise ValueError('lls_utils.fill_ions: Cannot have 0 subsystems..')
        else:
            self.ions = self.subsys['A'].ions
            self.clm_analy = self.subsys['A'].clm_analy
            print('lls_utils.fill_ions: Need to update multiple subsystems!! Taking A.')


    # #############
    def fill_lls_lines(self, bval=20.):
        """
        Generate an HI line list for an LLS.
        Goes into self.lls_lines 

        Parameters
        ----------
        bval : float (20.)  Doppler parameter in km/s
        """
        from barak import absorb as ba

        atom = ba.readatom()
        self.lls_lines = []
        for line in atom['HI']:
            tmp = abs_line.Abs_Line(line['wa'],fill=False)
            # Atomic data
            tmp.atomic = {'fval': line['osc'], 'gamma': line['gam'],
                          'name': 'HI %s' % line['wa'], 'wrest': line['wa']}
            tmp.name = tmp.atomic['name']
            # Attributes
            tmp.attrib['N'] = self.NHI
            tmp.attrib['b'] = bval
            tmp.z = self.zabs

            self.lls_lines.append(tmp)
            #xdb.set_trace()
        

    # Absorption model of the LLS (HI only)
    def flux_model(self,spec,smooth=0):
        """
        Generate a LLS model given an input spectrum

        Parameters:
          spec:  Barak Spectrum (will migrate to specutils.Spectrum1D)
          smooth : (0) Number of pixels to smooth by

        Returns:
          Output model is passed back as a Spectrum 
        """
        
        # ########
        # LLS first

        # Energies in LLS rest-frame
        wv_rest = spec.dispersion * u.AA / (self.zabs+1)
        energy = wv_rest.to(u.eV, equivalencies=u.spectral())

        # Get photo_cross and calcualte tau
        tau_LL = (10.**self.NHI / u.cm**2) * xatomi.photo_cross(1,1,energy)

        # ########
        # Now the Lyman series

        # Check for lines
        if 'lls_lines' not in self.__dict__.keys():
            self.fill_lls_lines()

        #xdb.set_trace()
        tau_Lyman = voigt.voigt_model(spec.dispersion, self.lls_lines, flg_ret=2)

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
            model.gauss_smooth(npix=smooth)

        # Return
        return model

    # 
    def get_zpeak(self):
        ''' Measure zpeak from an ionic transition
        '''
        if self.ions is None:
            print('get_zpeak: Need to fill ions with fill_ions first.')
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
        return ('[{:s}: {:s} {:s}, {:g}, NHI={:g}, tau={:g}, M/H={:g}]'.format(
                self.__class__.__name__,
                 self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True),
                 self.zabs, self.NHI, self.tau_LL, self.MH))

    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'LLS'

# #######################################################################
# #######################################################################
# #######################################################################
# Class for LLS Survey
class LLS_Survey(Absline_Survey):
    """An LLS Survey class

    Attributes:
        
    """
    # Initialize with a .dat file
    def __init__(self, dat_file, tree=None):
        # Generate with type
        Absline_Survey.__init__(self,dat_file,abs_type='LLS', tree=tree)

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

    # Default sample of LLS (high-z)
    @classmethod
    def default_sample(cls):
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
    #
    #flg_test += 2**9 # LLS Survey NHI
    flg_test += 2**10 # LLS Survey ions

    # Test Absorption System
    print('-------------------------')
    tmp1 = LLS_System(dat_file='Data/HE0940-1050.z2916.dat',
                      tree=os.environ.get('LLSTREE'))
    print(tmp1)

    # Test ions
    if (flg_test % 2**1) >= 2**0:
        print('-------------------------')
        tmp1.fill_ions()
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
        '''
        with open('tmp.yml', 'w') as outfile:
            outfile.write( yaml.safe_dump(tmp1.datdic, default_flow_style=False))
        dd = {'A':'a', 'B':{'C':'c', 'D':'d', 'E':'e'}}
        with open('result.yml', 'w') as yaml_file:
            yaml_file.write( yaml.safe_dump(dd, default_flow_style=False))
        xdb.set_trace()
        '''

    # #############################
    # LLS Survey
    if (flg_test % 2**10) >= 2**9:
        print('-------------------------')
        lls = LLS_Survey('Lists/lls_metals.lst', tree=os.environ.get('LLSTREE'))
        xdb.xhist(lls.NHI, binsz=0.30)

    # LLS Survey ions
    if (flg_test % 2**11) >= 2**10:
        lls = LLS_Survey('Lists/lls_metals.lst', tree=os.environ.get('LLSTREE'))
        lls.fill_ions()
        xdb.xhist(lls.ions((6,4),skip_null=True)['clm'], binsz=0.3,
                  xlabel=r'$\log_{10} N({\rm C}^{+3})$')
        xdb.xhist(lls.ions((14,2),skip_null=True)['clm'], binsz=0.3,
                  xlabel=r'$\log_{10} N({\rm Si}^{+})$')

    # All done
    print('-------------------------')
    print('lls_utils: All done testing..')
