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

from xastropy.xutils import xdebug as xdb

from astropy.utils.misc import isiterable

# Path for xastropy
#xa_path = imp.find_module('xastropy')[1]

#def ion_name(ion):
#def photo_cross(Z, ion, E, datfil=None, silent=False):

########################## ##########################
########################## ##########################
def load(flg=0, data_file=None,cosh_dct=None, pckl_fil=None):
    """ Load the data for COS-Halos

    Paramaeters
    ----------
    flg: integer
      Flag indicating how to load the data
      0 = IDL mega structure
      1 = FITS files
    data_file: string
      Name of data file
    pckl_fil: string
      Name of file for pickling

    JXP on 30 Nov 2014
    """
    from xastropy.cgm import core as xcc
    reload(xcc)

    # IDL save file
    if flg == 0:
        if data_file is None:
            data_file = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Halos/lowions/'+
                                        'coshalos_lowmetals_mega.sav')
        from scipy.io import readsav
        print('cos_halos.load:  Be patient...')
        if cosh_dct is None:
            cosh_dct = readsav(data_file)

        # Generate the CGM Survey
        cos_halos = xcc.CGM_Abs_Survey()
        ncos = len(cosh_dct['megastruct'])
        cos_halos.nsys = ncos
        for kk in range(ncos):
           #  
            cos_halos.abs_sys.append(xcc.CGM_Abs(
                ras=cosh_dct['megastruct'][kk]['galaxy']['qsora'][0],
                decs=cosh_dct['megastruct'][kk]['galaxy']['qsodec'][0],
                g_ras=cosh_dct['megastruct'][kk]['galaxy']['ra'][0],
                g_decs=cosh_dct['megastruct'][kk]['galaxy']['dec'][0],
                zgal=cosh_dct['megastruct'][kk]['galaxy']['zspec'][0]
                ))
    elif flg == 1: # FITS files
        fits_path = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Halos/lowions/FITS/')
        # Loop
        cos_files = glob.glob(fits_path+'J*.fits')
        # Setup
        cos_halos = xcc.CGM_Abs_Survey()
        cos_halos.nsys = len(cos_files)
        # Read
        for fil in cos_files:
            hdu = fits.open(fil)
            summ = hdu[0].data
            galx = hdu[1].data
            xdb.set_trace()
        
    else:
        raise ValueError('cos_halos.load: Not read for this flag {:d}'.format(flg))

    # Pickle?
    if pckl_fil is not None:
        xdb.set_trace()  # NOT GOING TO WORK
        pfil = open(pckl_fil, "wb")
        sys.setrecursionlimit(20000)
        pickle.dump(cos_halos,pfil,-1)
        pfil.close()
        print('cos_halos.load: Wrote pickle file {:s}'.format(pckl_fil))
    
    #
    return cos_halos
    


# Testing
if __name__ == '__main__':

    # Load pickle
    cos_halos = load(flg=1)
    print(cos_halos)
    
    #
    print('All done')
