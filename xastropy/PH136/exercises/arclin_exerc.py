"""Module to perform the in-class arc lamp exercise for PH136.
"""

# Import libraries
import numpy as np
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot
import pdb

#####################################
# Show a 1D spectrum 
#   import xastropy.PH136.exercises.arclin_exerc as ale 
#   reload(ale)
#   ale.plot_spec('tst.00000011.FIT.gz')
def plot_spec(fits_fil,prow=None):
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
    pyplot.plot(np.arange(npix), spec)
    pyplot.show()

    #pdb.set_trace()

################################################
# Fit a wavelength solution to hard-coded values, and plot
#   import xastropy.PH136.exercises.arclin_exerc as ale 
#   reload(ale)
#   ale.fit_lines()
def fit_lines():

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

    # Plot
    pyplot.plot(pix_val, wav_val, 'o')
    pyplot.plot(xval, yval)
    pyplot.show()
