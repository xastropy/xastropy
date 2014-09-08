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
    from barak.fitcont import fitqsocont
    from barak.spec import read
    from barak.io import saveobj, loadobj

    # Initialize
    if savfil == None:
        savfil = 'conti.sav'
    if outfil == None:
        outfil = 'conti.fits'
        
    # Read spectrum + convert to Barak format
    sp = xastropy.spec.readwrite.readspec(specfil)
    
    confirm = 'n'

    # Fit spline continuum:
    while(confirm == 'n'):
        if os.path.lexists(savfil): #'contfit_' + name + '.sav'):
            option = raw_input('Adjust old continuum? (y)/n: ')
            if option.lower() != 'n':
                co_old, knots_old = loadobj(savfil) #'contfit_' + name + '.sav')
                co, knots = fitqsocont(sp.wa, sp.fl, sp.er, redshift,
                                       oldco=co_old, knots=knots_old,
                                       divmult=divmult,
                                       forest_divmult=forest_divmult)
            else:
                co, knots = fitqsocont(sp.wa, sp.fl, sp.er, redshift,
                                       divmult=divmult,
                                       forest_divmult=forest_divmult)
        else:
            co, knots = fitqsocont(sp.wa, sp.fl, sp.er, redshift,
                                   divmult=divmult,
                                   forest_divmult=forest_divmult)

        os.remove('_knots.sav')

        # Save continuum:
        saveobj(savfil, (co, knots), overwrite=1)

        # Check continuum:
        print('Plotting new continuum')
        sp_nyquist = read(specfil)
        co_new = np.interp(sp_nyquist.wa, sp.wa, co)
        sp_new = np.rec.fromarrays([sp_nyquist.wa, sp_nyquist.fl,
                                    sp_nyquist.er, co_new],
                                   names='wa, fl, er, co')
        pl.clf()
        pl.plot(sp_new.wa, sp_new.fl, drawstyle='steps-mid')
        pl.plot(sp_new.wa, sp_new.co, color='r')
        pl.show()

        # Repeat?
        confirm = 'y'
        confirm = raw_input('Keep continuum? (y)/n: ')
    print name

    # Output

    # Data file with continuum
    fits.writeto(outfil, sp_new, clobber=True)

    #ascii.write(sp_new, name + '.txt', Writer=ascii.CommentedHeader)
    
    try:
        sp_unbinned = read(name + '_unbinned.fits')
        co_new = np.interp(sp_unbinned.wa, sp.wa, co)
        sp_new = np.rec.fromarrays([sp_unbinned.wa, sp_unbinned.fl,
                                    sp_unbinned.er, co_new],
                                   names='wa, fl, er, co')
        writetable(name + '_unbinned.fits', sp_new,
                   units=['angstrom', 'erg /s /cm**2 /angstrom',
                          'erg /s /cm**2 /angstrom',
                          'erg /s /cm**2 /angstrom'], overwrite=True)
    except:
        pass
