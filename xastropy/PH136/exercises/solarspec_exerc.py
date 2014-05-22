"""Module to perform the Solar Spectrum exercise exercise for PH136.
      Currently tuned to Lick spectra from Jan 2013
"""

# Import libraries
import numpy as np
from astropy.io import fits
from matplotlib import pyplot
import pdb

#####################################
# Show a 1D spectrum 
#   Useful to eyeball the pixel values of a few key lines
#   import xastropy.PH136.exercises.solarspec_exerc as ssp
#   reload(ssp)
#   ssp.plot_spec('b1014.fits')
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
        pyplot.clf()
        pyplot.plot(np.arange(npix), spec)
        pyplot.show()

    if give_spec:
        return spec
    else: 
        return

    #pdb.set_trace()

################################################
# Define a Gaussian plus a floor offset
def gauss_off(x, Z, A, x0, sigma):
    return Z + A*np.exp(- (x-x0)**2 / (2.*sigma**2) )

################################################
# Fit a wavelength solution to hard-coded values, and plot
def gaussfit_line(xval, yval, xguess, xrng=15.,debug=False):
    from scipy.optimize import curve_fit
    import xastropy.PH136.exercises.solarspec_exerc as ssp

    # Grab the region to fit
    gdx = np.where( np.fabs( xval-xguess ) < xrng )
    newx = xval[gdx]
    newy = yval[gdx]

    # Guess
    Aguess = np.max(newy)
    sguess = 2.
    Z = np.median(newy)
    pguess = [Z, Aguess, xguess, sguess]
    
    # Fit
    popt, pcov = curve_fit(ssp.gauss_off, newx, newy, p0=pguess)

    # Debug
    if debug:
        pyplot.clf()
        #data
        pyplot.plot(newx, newy, 'o')
        #curve
        xdum = np.linspace(np.min(newx), np.max(newx), num=300)
        Z,A,x0,sigma = popt
        ydum = Z + A*np.exp(- (xdum-x0)**2 / (2.*sigma**2) )
        # plot
        pyplot.plot(xdum, ydum)
        pyplot.show()

    return popt[2]

################################################
# Fit a wavelength solution using hard-coded guesses, and plot
#   import xastropy.PH136.exercises.solarspec_exerc as ssp
#   reload(ssp)
#   fit = ssp.fit_lines(fits_fil='b1014.fits')
def fit_lines(fits_fil=None, xquery=None,show_plot=False, plot_spec=False):

    import xastropy.PH136.exercises.solarspec_exerc as ssp
    from astropy.io.fits import getdata

    if fits_fil == None:
        fits_fil = 'b1014.fits'

    spec = ssp.plot_spec(fits_fil, give_spec=True)
    npix = len(spec)
    xpix = np.arange(npix)

    # Generate the arrays
    wav_val = np.array( [5085.82, 4799.92, 4358.33, 3466.55])
    guess_pix_val = np.array( [1930.5, 1664.72, 1241.46, 316.8])

    # Centroid those lines!
    nlin = len(guess_pix_val)
    pix_val = np.zeros(nlin)
    ii=0
    for gpix in guess_pix_val:
        pix_val[ii] = ssp.gaussfit_line(xpix,spec,gpix)#,debug=True)
        ii += 1
        #pdb.set_trace()

    # Fit
    fit = np.polyfit(pix_val, wav_val, 2)
    print 'Fit (dlam, w0): ', fit

    # Setup for plot
    pv = np.poly1d(fit)
    xval = np.linspace(1., 2000, 100)
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
        spec = ssp.plot_spec(fits_fil, give_spec=True, noplot=True)
        npix = len(spec)
        xval = np.arange(npix)
        wave = pv(xval)
        pyplot.clf()
        pyplot.plot(wave, spec,drawstyle="steps-mid", ls='-')
        pyplot.xlim([3000., 5500])
        pyplot.xlabel('Wavelength (Ang)')
        pyplot.ylabel('Counts')
        pyplot.savefig('arclin_spec.pdf')
        

    # Print a value
    if xquery != None:
        wquery = pv(xquery)
        print 'Wavelength for pixel = ', xquery, ' is wave = ', wquery

    return fit

################################################
# Extract and show the solar spectrum
#   import xastropy.PH136.exercises.solarspec_exerc as ssp
#   reload(ssp)
#   ssp.sol_spec()
def sol_spec(fits_fil=None, xquery=None,show_plot=False, plot_spec=False, arc_fil=None):

    import xastropy.PH136.exercises.solarspec_exerc as ssp

    # Get wavelength solution
    if arc_fil == None:
        arc_fil = 'b1014.fits'
    fit = ssp.fit_lines(fits_fil=arc_fil)
    pv = np.poly1d(fit)

    # Read Solar spectrum
    if fits_fil == None:
        fits_fil = 'b1029.fits'

    spec = ssp.plot_spec(fits_fil, give_spec=True)
    npix = len(spec)

    xpix = np.arange(npix)
    wave = pv(xpix)

    # Plot
    pyplot.clf()
    pyplot.plot(wave, spec,drawstyle="steps-mid", ls='-')
    pyplot.xlim([3000., 5500])
    pyplot.xlabel('Wavelength (Ang)')
    pyplot.ylabel('Counts')
    # Ca lines
    pyplot.axvline(3933.68, color='r')
    pyplot.axvline(3968.47, color='r')
    pyplot.show()

    # Ca H+K
    # 3955.5, 3991.
