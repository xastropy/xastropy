"""
#;+ 
#; NAME:
#; casbah.igm_spec
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for analyzing, plotting, etc IGM Spectra for CASBAH
#;      Not much done yet (nothing really)
#;   13-Jan-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp

import matplotlib as mpl
from matplotlib import pyplot as plt

from astropy.io import fits, ascii
from astropy import units as u 
from astropy.table.table import Table
#from astropy import constants as const

from linetools.spectra import io as lsio

#from xastropy.spec import readwrite as xsr

from xastropy.xutils import xdebug as xdb

def spec_plot(spec_fil, ID_path=None, ambig_fil=None, dwv=25., outfil='tmp.pdf'):
    mpl.rcParams['font.family'] = 'stixgeneral'
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.gridspec as gridspec

    '''
    spec_fil: string
      Filename of spectrum
    outfil: string ('tmp.pdf')
      Filename for output
    dwv: float (25.)
      Ang per window
    '''

    npx = 1
    npy = 4

    # Read spectrum
    spec = lsio.readspec(spec_file)

    # Start plot
    pp = PdfPages(outfil)

#def spec_plot(spec_fil, ID_path=None, ambig_fil=None, dwv=25., outfil='tmp.pdf'):
    
# ################
def main(args):

    if args[0] == 0:
        # Generate PDF of IGM spectrum with IDs
        spec_plot(args[1:])

##
if __name__ == "__main__":

    if len(sys.argv) == 1: # Testing
        flg_fig = 0 
        flg_fig += 2**0  # SLLS 
    else:
        id_job = int(sys.argv[1])  # 1 = spec_plot
        args = sys.argv[1:]

        if id_job == 1:
            run_specplot()

