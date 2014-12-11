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
#from astropy import constants as const

from xastropy.igm.abs_sys.abssys_utils import Absline_System
from xastropy.galaxy.core import Galaxy
#from xastropy.cgm.core import CGM_Abs, CGM_Abs_Survey
from xastropy.igm.abs_sys.abs_survey import Absline_Survey
from xastropy.igm.abs_sys.ionic_clm import Ions_Clm
from xastropy.cgm.core import CGM_Abs_Survey, CGM_Sys
from xastropy.xutils import xdebug as xdb
from xastropy import spec as xspec
from xastropy.kinematics.absline import Kin_Abs

from astropy.utils.misc import isiterable

# Path for xastropy
#xa_path = imp.find_module('xastropy')[1]

#def ion_name(ion):
#def photo_cross(Z, ion, E, datfil=None, silent=False):

# Class for COS_Halos Survey
class COS_Halos(CGM_Abs_Survey):
    """Inherits CGM Abs Survey

    Attributes:
    """
    # Initialize with a .dat file
    def __init__(self, tree=None):

        # Generate with type
        CGM_Abs_Survey.__init__(self)
        self.survey = 'COS-Halos'
        self.ref = 'Tumlinson+11; Werk+12; Tumlinson+13; Werk+13'


    # Load from mega structure
    def load_mega(self,flg=1, data_file=None,cosh_dct=None, pckl_fil=None,
                  skip_ions=False):
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
            cos_files = glob.glob(fits_path+'/J*.fits')
            #cos_files = glob.glob(fits_path+'/J091*.fits') # For testing
            # Setup
            self.nsys = len(cos_files)
            # Read
            for fil in cos_files:
                print('cos_halos: Reading {:s}'.format(fil))
                mm = cos_files.index(fil)
                hdu = fits.open(fil)
                summ = hdu[1].data
                galx = hdu[2].data
                self.cgm_abs.append(CGM_Sys(
                    ras=galx['qsora'][0],
                    decs=galx['qsodec'][0],
                    g_ras=galx['ra'][0],
                    g_decs=galx['dec'][0],
                    zgal=summ['zfinal'][0]
                    ))
                # COS-Halos naming
                self.cgm_abs[mm].field = galx['field'][0]
                self.cgm_abs[mm].gal_id = galx['galid'][0]
                # Ions
                if skip_ions is True:
                    continue
                self.cgm_abs[mm].abs_sys.ions = Ions_Clm()
                self.cgm_abs[mm].abs_sys.ions.ion_data = {}
                for jj in range(summ['nion'][0]):
                    iont = hdu[3+jj].data
                    zion = (iont['zion'][0][0], iont['zion'][0][1])
                    self.cgm_abs[mm].abs_sys.ions.ion_data[zion] = {}
                    for key in self.cgm_abs[mm].abs_sys.ions.keys:
                        try:
                            self.cgm_abs[mm].abs_sys.ions.ion_data[zion][key] = iont[key][0]
                        except KeyError:
                            if key == 'flg_inst':
                                self.cgm_abs[mm].abs_sys.ions.ion_data[zion][key] = 0
                            else:
                                xdb.set_trace()
                # NHI
                self.cgm_abs[mm].abs_sys.NHI = self.cgm_abs[mm].abs_sys.ions.ion_data[(1,1)]['clm']
            # Mask
            self.mask = np.ones(cos_halos.nsys, dtype=bool)
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
                    wrest = kin_init['mtl_wr'][mt] 
                    if wrest <= 1:
                        xdb.set_trace()
                    spec = get_coshalo_spec( cgm_abs, wrest )
                    vmnx = (kin_init['L_vmn'][mt], kin_init['L_vmx'][mt]) 
                    # Process
                    cgm_abs.abs_sys.kin['Metal'] = Kin_Abs(wrest, vmnx)
                    cgm_abs.abs_sys.kin['Metal'].orig_kin(spec)
                else:
                    # Fill with zeros (for the keys)
                    cgm_abs.abs_sys.kin['Metal'] = Kin_Abs(0., (0., 0.))

                # HI
                if kin_init['flgH'][mt] > 0:
                    wrest = kin_init['HI_wrest'][mt] 
                    if wrest <= 1:
                        xdb.set_trace()
                    spec = get_coshalo_spec( cgm_abs, wrest )
                    vmnx = (kin_init['HIvmn'][mt], kin_init['HIvmx'][mt]) 
                    # Process
                    cgm_abs.abs_sys.kin['HI'] = Kin_Abs(wrest, vmnx)
                    cgm_abs.abs_sys.kin['HI'].orig_kin(spec)
                else:
                    # Fill with zeros (for the keys)
                    cgm_abs.abs_sys.kin['HI'] = Kin_Abs(0., (0., 0.))


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
        mt = np.where( np.fabs(tab['col1']-wrest) < 1e-3)[0]
        if len(mt) == 0:
            raise ValueError('get_coshalo_spec: wrest={:g} not found!'.format(wrest))
        mt = mt[0]
        trans = tab['col2'][mt]+tab['col3'][mt]

        # Read
        slicedir = cdir+fielddir+sysdir+'/fitting/'
        slicename = sysname+'_'+trans+'_slice.fits'
        spec = xspec.readwrite.readspec(slicedir+slicename,
                                        flux_tags=['FNORM'], sig_tags=['ENORM'])
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
        cos_halos = COS_Halos()
        cos_halos.load_mega()
        print(cos_halos)
    
    # Simple rho vs NHI plot
    if (flg_fig % 2**2) >= 2**1:
        cos_halos = COS_Halos()
        cos_halos.load_mega()
        x= cos_halos.rho
        y= cos_halos.NHI
        xdb.xplot(x, y, scatter=True)
    #
    # Simple kinematics
    if (flg_fig % 2**3) >= 2**2:
        cos_halos = COS_Halos()
        cos_halos.load_mega()
        cos_halos.load_abskin()
        # Plot
        mtl_kin = cos_halos.abs_kin('Metal')
        gd = np.where(mtl_kin['flg'] == 1)
        xdb.xplot(cos_halos.NHI[gd], mtl_kin['Dv'][gd], scatter=True)

        HI_kin = cos_halos.abs_kin('HI')
        gd = np.where(HI_kin['flg'] == 1)
        xdb.xplot(cos_halos.NHI[gd], HI_kin['Dv'][gd], scatter=True)
    print('All done')
