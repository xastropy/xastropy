"""
#;+ 
#; NAME:
#; cos_halos
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for COS-Halos analysis
#;   29-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp, pickle, sys, glob
from astropy.io import fits, ascii
from astropy import units as u 
from astropy.table import QTable, Table, Column
#from astropy import constants as const

from linetools.spectra import io as lsio

from xastropy.galaxy.core import Galaxy
#from xastropy.cgm.core import CGM_Abs, CGM_Abs_Survey
from xastropy.cgm.core import CGMAbsSurvey, CGMSys
from xastropy.xutils import xdebug as xdb
from xastropy.igm.abs_sys.ionclms import IonClms
from xastropy.kinematics.absline import Kin_Abs

from astropy.utils.misc import isiterable

# Path for xastropy
#xa_path = imp.find_module('xastropy')[1]

#def ion_name(ion):
#def photo_cross(Z, ion, E, datfil=None, silent=False):

# Class for COS_Halos Survey
class COSHalos(CGMAbsSurvey):
    """Inherits CGM Abs Survey

    Attributes:
    """
    # Initialize with a .dat file
    def __init__(self, tree=None):

        # Generate with type
        CGMAbsSurvey.__init__(self)
        self.survey = 'COS-Halos'
        self.ref = 'Tumlinson+11; Werk+12; Tumlinson+13; Werk+13'


    # Load from mega structure
    def load_mega(self,flg=1, data_file=None,cosh_dct=None, pckl_fil=None,
                  skip_ions=False, test=False):
        """ Load the data for COS-Halos

        Paramaeters
        ----------
        flg: integer (1)
          Flag indicating how to load the data
          0 = IDL mega structure
          1 = FITS files from Dropbox
        data_file: string
          Name of data file
        pckl_fil: string
          Name of file for pickling

        JXP on 30 Nov 2014
        """
        #from xastropy.cgm import core as xcc
        #reload(xcc)

        # IDL save file
        if flg == 0:
            if data_file is None:
                data_file = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Halos/lowions/'+
                                            'coshalos_lowmetals_mega.sav')
            '''
            from scipy.io import readsav
            print('cos_halos.load:  Be patient...')
            if cosh_dct is None:
                cosh_dct = readsav(data_file)
    
            # Generate the CGM Survey
            ncos = len(cosh_dct['megastruct'])
            self.nsys = ncos
            for kk in range(ncos):
            #  
                self.cgm_abs.append(CGM_Abs(
                    ras=cosh_dct['megastruct'][kk]['galaxy']['qsora'][0],
                    decs=cosh_dct['megastruct'][kk]['galaxy']['qsodec'][0],
                    g_ras=cosh_dct['megastruct'][kk]['galaxy']['ra'][0],
                    g_decs=cosh_dct['megastruct'][kk]['galaxy']['dec'][0],
                    zgal=cosh_dct['megastruct'][kk]['galaxy']['zspec'][0]
                    ))
            '''
        elif flg == 1: # FITS files
            fits_path = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Halos/lowions/FITS')
            # Loop
            if test is True:
                cos_files = glob.glob(fits_path+'/J091*.fits') # For testing
            else:
                cos_files = glob.glob(fits_path+'/J*.fits')
            # Setup
            self.nsys = len(cos_files)
            # Read
            for fil in cos_files:
                print('cos_halos: Reading {:s}'.format(fil))
                mm = cos_files.index(fil)
                hdu = fits.open(fil)
                summ = hdu[1].data
                galx = hdu[2].data
                self.cgm_abs.append( CGMSys(ras=galx['qsora'][0],
                    decs=galx['qsodec'][0],
                    g_ras=galx['ra'][0],
                    g_decs=galx['dec'][0],
                    zgal=summ['zfinal'][0]
                    ))
                # COS-Halos naming
                self.cgm_abs[mm].field = galx['field'][0]
                self.cgm_abs[mm].gal_id = galx['galid'][0]
                # Galxy properties
                self.cgm_abs[mm].galaxy.halo_mass = summ['LOGMHALO'][0] 
                self.cgm_abs[mm].galaxy.stellar_mass = summ['LOGMFINAL'][0] 
                # Ions
                if skip_ions is True:
                    continue
                self.cgm_abs[mm].abs_sys.ions = IonClms()
                all_Z = []
                all_ion = []
                for jj in range(summ['nion'][0]):
                    iont = hdu[3+jj].data
                    if jj == 0: # Generate new Table
                        dat_tab = Table(iont)
                    else:
                        try:
                            dat_tab.add_row(Table(iont)[0])
                        except:
                            xdb.set_trace()
                    all_Z.append(iont['zion'][0][0])
                    all_ion.append(iont['zion'][0][1])
                    '''
                    for key in self.cgm_abs[mm].abs_sys.ions.keys:
                        try:
                            self.cgm_abs[mm].abs_sys.ions.ion_data[zion][key] = iont[key][0]
                        except KeyError:
                            if key == 'flg_inst':
                                self.cgm_abs[mm].abs_sys.ions.ion_data[zion][key] = 0
                            else:
                                xdb.set_trace()
                    '''
                # Add Z,ion
                dat_tab.add_column(Column(all_Z,name='Z'))
                dat_tab.add_column(Column(all_ion,name='ion'))
                # Set
                self.cgm_abs[mm].abs_sys.ions._data = dat_tab
                # NHI
                self.cgm_abs[mm].abs_sys.NHI = self.cgm_abs[mm].abs_sys.ions[(1,1)]['CLM']
            # Mask
            self.mask = np.ones(self.nsys, dtype=bool)
        else:
            raise ValueError('cos_halos.load: Not ready for this flag {:d}'.format(flg))

        '''
        # Pickle?
        if pckl_fil is not None:
            xdb.set_trace()  # NOT GOING TO WORK
            pfil = open(pckl_fil, "wb")
            sys.setrecursionlimit(20000)
            pickle.dump(cos_halos,pfil,-1)
            pfil.close()
            print('cos_halos.load: Wrote pickle file {:s}'.format(pckl_fil))
        '''
    
    
    ########################## ##########################
    def load_abskin(self,flg=1,kin_init_file=None):
        """ Load the absorption-line kinematic data for COS-Halos
        Calculate from scratch if needed

        Paramaeters
        ----------
        flg: integer (1)
          Flag indicating how to load the data
        0 = Load from file
        1 = Generate
        kin_init_file: string
          Name of kinematics driver file
    
        JXP on 10 Dec 2014
        """
    
        if flg == 1: # Generate
            # Read init file
            if kin_init_file is None:
                kin_init_file = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Halos/Kin/'+
                                                  'coshalo_kin_driver.dat')
            kin_init = ascii.read(kin_init_file,guess=False)
    
            # Loop to my loop
            fgal = zip(self.field, self.gal_id)
            for cgm_abs in self.cgm_abs:
                # Match to kin_init
                mt = np.where( (cgm_abs.field == kin_init['QSO']) &
                               (cgm_abs.gal_id == kin_init['Galaxy']) )[0]
                if len(mt) == 0:
                    print('load_kin: No kinemtaics for {:s}, {:s}'.format(cgm_abs.field,
                                                                          cgm_abs.gal_id))
                    continue
                mt = mt[0]

                # Metals
                if kin_init['flgL'][mt] > 0:
                    wrest = kin_init['mtl_wr'][mt]*u.AA 
                    if wrest.value <= 1:
                        xdb.set_trace()
                    spec = get_coshalo_spec( cgm_abs, wrest )
                    vmnx = (kin_init['L_vmn'][mt]*u.km/u.s, kin_init['L_vmx'][mt]*u.km/u.s)
                    # Process
                    cgm_abs.abs_sys.kin['Metal'] = Kin_Abs(wrest, vmnx)
                    cgm_abs.abs_sys.kin['Metal'].fill_kin(spec, per=0.07)
                    # Save spec
                    cgm_abs.abs_sys.kin['Metal'].spec = spec
                else:
                    # Fill with zeros (for the keys)
                    cgm_abs.abs_sys.kin['Metal'] = Kin_Abs(0.*u.AA, (0., 0.))

                # HI
                if kin_init['flgH'][mt] > 0:
                    wrest = kin_init['HI_wrest'][mt]*u.AA 
                    if wrest.value <= 1:
                        xdb.set_trace()
                    spec = get_coshalo_spec( cgm_abs, wrest )
                    vmnx = (kin_init['HIvmn'][mt]*u.km/u.s, kin_init['HIvmx'][mt]*u.km/u.s) 
                    # Process
                    cgm_abs.abs_sys.kin['HI'] = Kin_Abs(wrest, vmnx)
                    cgm_abs.abs_sys.kin['HI'].fill_kin(spec, per=0.07)
                    cgm_abs.abs_sys.kin['HI'].spec = spec
                else:
                    # Fill with zeros (for the keys)
                    cgm_abs.abs_sys.kin['HI'] = Kin_Abs(0.*u.AA, (0., 0.))


            #tmp = cos_halos.abs_kin('Metal')['Dv']
            #xdb.set_trace()


# 
def get_coshalo_spec(cgm_abs, wrest):
        """ Load the absorption-line kinematic data for COS-Halos
        Calculate from scratch if needed

        Paramaeters
        ----------
        cgm_abs: CGM_Abs Class
        wrest: float
          Rest wavelength for spectrum of interest
    
        JXP on 11 Dec 2014
        """
        # Directories
        cdir = os.environ.get('DROPBOX_DIR')+'/COS-Halos/'
        fielddir = 'Targets/'+cgm_abs.field+'/'
        sysdir = cgm_abs.gal_id+'_z{:5.3f}'.format(cgm_abs.galaxy.z)
        sysname = cgm_abs.field+'_'+sysdir

        # Transition
        templ_fil = os.environ.get('DROPBOX_DIR')+'/COS-Halos/Targets/system_template.lst'
        tab = ascii.read(templ_fil)
        mt = np.where( np.fabs(tab['col1']-wrest.value) < 1e-3)[0]
        if len(mt) == 0:
            raise ValueError('get_coshalo_spec: wrest={:g} not found!'.format(wrest))
        mt = mt[0]
        trans = tab['col2'][mt]+tab['col3'][mt]

        # Read
        slicedir = cdir+fielddir+sysdir+'/fitting/'
        slicename = sysname+'_'+trans+'_slice.fits'
        spec = lsio.readspec(slicedir+slicename, flux_tags=['FNORM'], sig_tags=['ENORM'])
        # Fill velocity
        spec.velo = spec.relative_vel((cgm_abs.galaxy.z+1)*wrest)
    
        #spec.qck_plot()
        return spec

########################## ##########################
# Testing
if __name__ == '__main__':

    flg_fig = 0 
    #flg_fig += 1  # Load FITS
    #flg_fig += 2  # NHI plot
    flg_fig += 2**2  # Simple Kinematics

    # Load FITS
    if (flg_fig % 2) == 1:
        cos_halos = COSHalos()
        cos_halos.load_mega()
        print(cos_halos)
    
    # Simple rho vs NHI plot
    if (flg_fig % 2**2) >= 2**1:
        cos_halos = COSHalos()
        cos_halos.load_mega()
        x= cos_halos.rho
        y= cos_halos.NHI
        xdb.xplot(x, y, scatter=True)
    #
    # Simple kinematics
    if (flg_fig % 2**3) >= 2**2:
        cos_halos = COSHalos()
        cos_halos.load_mega()#test=True)
        cos_halos.load_abskin()
        # Plot
        mtl_kin = cos_halos.abs_kin('Metal')
        gd = np.where(mtl_kin['flg'] > 0)
        xdb.xplot(cos_halos.NHI[gd], mtl_kin['Dv'][gd], scatter=True)

        HI_kin = cos_halos.abs_kin('HI')
        gd = np.where(HI_kin['flg'] > 0)
        xdb.xplot(cos_halos.NHI[gd], HI_kin['Dv'][gd], scatter=True)
    print('All done')
