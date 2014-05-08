"""Module to perform the in-class arc lamp exercise for PH136.
"""

# Import libraries
import numpy as np
from astropy.io import fits
from matplotlib import pyplot
import pdb

#####################################
# Show a 1D spectrum 
#   import xastropy.PH136.exercises.arclin_exerc as ale 
#   reload(ale)
#   ale.plot_spec('tst.00000011.FIT.gz')
def plot_spec(fits_fil,prow=None,give_spec=False, noplot=False):
    from astropy.io.fits import getdata

    # Read
    arr,head = getdata(fits_fil,0,header=True)
    siz = arr.shape
    
    if prow == None:
        prow = int(siz[0]/2.)

    # Define the spectrum
    spec = arr[prow,:]
    npix = len(spec)

    #pdb.set_trace()
    # Plot
    if noplot:
        pyplot.plot(np.arange(npix), spec)
        pyplot.show()

    if give_spec:
        return spec
    else: 
        return

    #pdb.set_trace()

################################################
# Fit a wavelength solution to hard-coded values, and plot
#   import xastropy.PH136.exercises.arclin_exerc as ale 
#   reload(ale)
#   ale.fit_lines()
#   ale.fit_lines(fits_fil='tst.00000011.FIT.gz')
def fit_lines(xquery=None,show_plot=True, plot_spec=True, fits_fil=None):

    import xastropy.PH136.exercises.arclin_exerc as ale 

    # Generate the arrays
    pix_val = np.array( [70.71, 75.5, 147.2, 403.26] )
    wav_val = np.array( [5790.663, 5769.698, 5460.735, 4358.327] )

    # Fit
    fit = np.polyfit(pix_val, wav_val, 1)
    print 'Fit (dlam, w0): ', fit

    # Setup for plot
    pv = np.poly1d(fit)
    xval = np.linspace(1., 500, 100)
    yval = pv(xval)

    # Plot?
    if show_plot:
        pyplot.clf()
        pyplot.plot(pix_val, wav_val, 'o')
        pyplot.plot(xval, yval)
        pyplot.xlabel('pixel')
        pyplot.ylabel('Wave (Ang)')
        #pyplot.show()
        pyplot.savefig('arclin_fit.pdf')

    # Plot the spectrum
    if plot_spec and (fits_fil != None):
        spec = ale.plot_spec(fits_fil, give_spec=True, noplot=True)
        npix = len(spec)
        xval = np.arange(npix)
        wave = pv(xval)
        pyplot.clf()
        pyplot.plot(wave, spec,drawstyle="steps-mid", ls='-')
        pyplot.xlim([4000., 6000])
        pyplot.xlabel('Wavelength (Ang)')
        pyplot.ylabel('Counts')
        pyplot.savefig('arclin_spec.pdf')
        

    # Print a value
    if xquery != None:
        wquery = pv(xquery)
        print 'Wavelength for pixel = ', xquery, ' is wave = ', wquery
