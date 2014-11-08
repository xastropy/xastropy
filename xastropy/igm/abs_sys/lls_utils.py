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

from __future__ import print_function

import os, copy

import numpy as np

from astropy import units as u
from astropy.io import ascii 

from abssys_utils import Absline_System
from ionic_clm import Ionic_Clm, Ionic_Clm_File
from xastropy.spec import abs_line, voigt
from xastropy.atomic import ionization as xatomi
from xastropy.xutils import xdebug as xdb

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
        if tree != None: self.tree = tree
        # Parse .dat file
        if dat_file != None:
            self.parse_dat_file(tree+dat_file)

        # Set tau_LL
        self.tau_LL = (10.**self.NHI)*6.3391597e-18 # Should replace with photocross
        self.ionic = {}

        # Name

    # Modify standard dat parsing
    def parse_dat_file(self,dat_file):
        # Standard Call
        out_list = Absline_System.parse_dat_file(self,dat_file,flg_out=1)
        datdic = out_list[0]

        # LLS keys
        self.MH = float(datdic['[M/H]ave'])
        self.nsub = int(datdic['Nsubsys'])
        self.cldyfil = datdic['CloudyGridFile']

        # LLS Subsystems
        if self.nsub > 0:
            subsys = {}
            lbls= map(chr, range(65, 91))
            keys = (['zabs','NHI','NHIsig','NH','NHsig','logx','sigx','b','bsig','Abundfile',
                     'U','Usig','flg_low','flg_alpha','[alpha/H]','sig[a/H]',
                     'flg_Fe','[Fe/H]','sig[Fe/H]','VPFITfile'])
            for i in range(self.nsub):
                # Generate
                subsys[lbls[i]] = self.subsys()
                # Fill in
                for key in list(subsys[lbls[i]].keys()):
                    try:
                        tmpc = datdic[lbls[i]+key]
                    except:
                        print('lls_utils: Key %s not found', lbls[i]+key)
                    else:  # Convert
                        val = subsys[lbls[i]][key]
                        #pdb.set_trace()
                        if val.__class__ == np.ndarray:  
                            subsys[lbls[i]][key] = np.array(map(float,tmpc.split()))
                        else: # Single value
                            subsys[lbls[i]][key] = (map(type(val),[tmpc]))[0]
            # Encode
            self.subsys = subsys

    # #############
    def fill_lls_lines(self, bval=20.):
        """
        Generate a line list for an LLS.
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
        
    # Generate the Ionic dictionary
    def get_ions(self):
        lbls= map(chr, range(65, 91))
        # Loop on Sub-Systems
        for kk in range(self.nsub):
            # Read .clm file
            tmp = Ionic_Clm_File(self.tree+self.subsys[lbls[kk]]['Abundfile'])
            # Fill it up
            print('lls_utils.get_ions: The next line needs to be changed!')
            Dumb_Class = type('Dummy_Object', (object,), {})
            self.subsys[lbls[kk]]['Ionic'] = Dumb_Class()
            self.subsys[lbls[kk]]['Ionic'].analy = tmp # THIS NEEDS TO BE CHANGED
            
            #pdb.set_trace()
            #self.ionic[lbls[kk],
        #pdb.set_trace()

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
        wv_rest = spec.wa * u.AA / (self.zabs+1)
        energy = wv_rest.to(u.eV, equivalencies=u.spectral())

        # Get photo_cross and calcualte tau
        tau_LL = (10.**self.NHI / u.cm**2) * xatomi.photo_cross(1,1,energy)

        # ########
        # Now the Lyman series

        # Check for lines
        if 'lls_lines' not in self.__dict__.keys():
            self.fill_lls_lines()

        #xdb.set_trace()
        tau_Lyman = voigt.voigt_model(spec.wa, self.lls_lines, flg_ret=2)

        # Combine
        tau_model = tau_LL + tau_Lyman

        # Kludge around the limit
        pix_LL = np.argmin( np.fabs( wv_rest- 911.3*u.AA ) )
        pix_kludge = np.where( (wv_rest > 911.5*u.AA) & (wv_rest < 912.8*u.AA) )[0]
        tau_model[pix_kludge] = tau_model[pix_LL]
        
        # Fill in flux
        model = copy.deepcopy(spec)
        model.fl = np.exp(-1. * tau_model)

        # Smooth?
        if smooth > 0:
            model.gauss_smooth(npix=smooth)

        # Return
        return model

        #spec.qck_plot()
        #xdb.set_trace()


    # Subsystem Dict
    def subsys(self):
        keys = (['zabs','NHI','NHIsig','NH','NHsig','logx','sigx','b','bsig','Abundfile',
                'U','Usig','flg_low','flg_alpha','[alpha/H]','sig[a/H]',
                'flg_Fe','[Fe/H]','sig[Fe/H]','VPFITfile'])
        values = ([0., 0., np.zeros(2), 0., np.zeros(2), 0., np.zeros(2), 0., 0.,
                   '', 0., np.zeros(2), 0, 0, 0., 0., 0, 0., 0., ''])
        return dict(zip(keys,values))

    # Output
    def __repr__(self):
        return ('[%s: %s %s, %g, NHI=%g, tau=%g, M/H=%g]' %
                (self.__class__.__name__,
                 self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True),
                 self.zabs, self.NHI, self.tau_LL, self.MH))

















## #################################    
## #################################    
## TESTING
## #################################    
if __name__ == '__main__':

    # Test Absorption System
    tmp1 = LLS_System(dat_file='Data/HE0940-1050.z2916.dat',
                      tree=os.environ.get('LLSTREE'))
    print(tmp1)
    print(tmp1.subsys)
    tmp1.get_ions()
    #tmp1.ionic['1215.6701'] = Ionic_Clm(1215.6701)
    #print(tmp1.ionic)
    #print(tmp1.ionic['1215.6701'].wave)
    
    tmp1.fill_lls_lines()
    #print(tmp1.lls_lines)

    from barak import spec as bs
    spec = bs.Spectrum(wa=np.linspace(3400.0,5000.0,10000))

    model = tmp1.flux_model(spec, smooth=4)
    model.qck_plot()

    print('lls_utils: All done testing..')
