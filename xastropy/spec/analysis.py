"""
#;+ 
#; NAME:
#; analysis
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for Analysis of Spectra
#;   07-Sep-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""

import barak
import xastropy
import numpy as np
import matplotlib.pyplot as plt

#name = "SDSSJ114436.66+095904.9"
#redshift = 3.1483297
#sp = read(name + ".fits")
##sp = read(name + ".fits")
#contfit(name, sp, redshift,forest_divmult=3)



#### ###############################
#  Calls Barak routines to fit the continuum
#    Stolen from N. Tejos by JXP
#
def x_contifit(specfil, outfil=None, savfil=None, redshift=0., divmult=1, forest_divmult=1):

    import os
    import barak.fitcont as bf
    from barak.spec import read
    from barak.io import saveobj, loadobj
    import xastropy.spec.readwrite as xsr
    reload(xsr)
    reload(bf)

    # Initialize
    if savfil == None:
        savfil = 'conti.sav'
    if outfil == None:
        outfil = 'conti.fits'
        
    # Read spectrum + convert to Barak format
    sp = xsr.readspec(specfil)
    
    confirm = 'n'

    # Fit spline continuum:
    while(confirm == 'n'):
        if os.path.lexists(savfil): #'contfit_' + name + '.sav'):
            option = raw_input('Adjust old continuum? (y)/n: ')
            if option.lower() != 'n':
                co_old, knots_old = loadobj(savfil) #'contfit_' + name + '.sav')
                co, knots = bf.fitqsocont(sp.wa, sp.fl, sp.er, redshift,
                                       oldco=co_old, knots=knots_old,
                                       divmult=divmult,
                                       forest_divmult=forest_divmult)
            else:
                co, knots = bf.fitqsocont(sp.wa, sp.fl, sp.er, redshift,
                                       divmult=divmult,
                                       forest_divmult=forest_divmult)
        else:
            co, knots = bf.fitqsocont(sp.wa, sp.fl, sp.er, redshift,
                                   divmult=divmult,
                                   forest_divmult=forest_divmult)

        os.remove('_knots.sav')

        # Save continuum:
        saveobj(savfil, (co, knots), overwrite=1)

        # Check continuum:
        print('Plotting new continuum')
        #sp_nyquist = read(specfil)
        #co_new = np.interp(sp_nyquist.wa, sp.wa, co)
        #sp_new = np.rec.fromarrays([sp_nyquist.wa, sp_nyquist.fl,
        #                            sp_nyquist.er, co_new],
        #                           names='wa, fl, er, co')
        plt.clf()
        plt.plot(sp.wa, sp.fl, drawstyle='steps-mid')
        plt.plot(sp.wa, sp.co, color='r')
        plt.show()

        # Repeat?
        confirm = 'y'
        confirm = raw_input('Keep continuum? (y)/n: ')
    #print name

    ## Output
    # Data file with continuum
    fits.writeto(outfil, sp, clobber=True)

#### ###############################
#  Calls plotvel (Crighton)
#    Adapted from N. Tejos scripts
#
def velplt(specfil):

    # Imports
    from plotspec import plotvel_util as pspv
    reload(pspv)
    import xastropy as xa
    from subprocess import Popen

    # Initialize
    if 'f26_fil' not in locals():
        f26_fil = 'tmp.f26'
        command = ['touch',f26_fil]
        print Popen(command)
        print 'xa.spec.analysis.velplt: Generated a dummy f26 file -- ', f26_fil
    if 'transfil' not in locals():
        path = xa.__path__
        transfil = path[0]+'/spec/Data/initial_search.lines'
    
    # Call
    pspv.main([specfil, 'f26='+f26_fil, 'transitions='+transfil])
