"""Module to perform the CCD exercise for PH136.
"""

# Import libraries
from numpy import *
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages

# ###############
# Read FITS file, plot image
def plot_img(pp,fil):

    # Imports
    from matplotlib import pyplot

    # Read the image
    from astropy.io.fits import getdata
    arr,head = getdata(fil,0,header=True)

    # Quick stats
    rms = std(arr)
    medi = median(arr)

    # Plot away
    pyplot.clf()
    pyplot.imshow(arr, aspect='equal',vmin=medi-2*rms,vmax=medi+2*rms)
    pyplot.colorbar()

    # Title (strip out path)
    ipos = fil.rfind('/')
    if ipos >= 0:
        name=fil[ipos+1:]
    else:
        name=fil
    pyplot.title(name)

    # Save to PDF
    pp.savefig()

    
# ###############
# Read FITS file, plot histogram
def hist_plot(pp,fil):

    # Imports
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib import pyplot

    # Read FITS file
    hdulist=fits.open(fil)
    head = hdulist[0].header
    arr = hdulist[0].data  # this is 1D, I think
    hdulist.close

    binsz = 1.  # insure it is a float

    # Find the range
    minv = amin(arr)
    maxv = amax(arr)

    # Set the boundaries sensibly given binsz
    i0 = int( minv / binsz) - 1
    i1 = int( maxv / binsz) + 1
    rng = tuple( binsz*array([i0,i1]) )
    nbin = i1-i0

    # Histogram
    hist, edges = histogram(arr, range=rng, bins=nbin)

    # plot the histogram
    pyplot.clf()
    pyplot.bar(edges[:-1], hist, width=binsz)

    # Time for stats
    rms = std(arr)
    avg = mean(arr)
    mx = amax(arr)
    mn = amin(arr)
    print 'Stats: ', avg, rms

    # Reset the plot range
    pyplot.xlim([avg-3.*rms, avg+3*rms])
    # Axis values
    axisv = pyplot.axis()

    # Label
    pyplot.xlabel('counts')
    pyplot.ylabel('N')
    s = "AVG= {:6.1f}".format(avg)
    s2 = "RMS= {:6.2f}".format(rms)
    s3 = "MIN= {}".format(mn)
    s4 = "MAX= {}".format(mx)

    pyplot.text(avg-2.5*rms, 0.9*axisv[3], s)
    pyplot.text(avg-2.5*rms, 0.8*axisv[3], s2)
    pyplot.text(avg-2.5*rms, 0.7*axisv[3], s3)
    pyplot.text(avg-2.5*rms, 0.6*axisv[3], s4)

    # Generate the Gaussian
    xval = arange(float(rng[0])-2.*rms, float(rng[1])+2*rms, binsz/5.)
    yval = amax(hist) * exp(-1. * (xval-avg)**2 / (2 * rms**2) )

    # plot
    pyplot.plot(xval, yval, 'r', lw=2)

    # Save to PDF
    pp.savefig()

#### ###############################
#  Loop through all .FIT files in the path
#  cep.make_all_plot('/Users/xavier/Class/PH136/Data/Tests/CCD/')
def make_all_plot(path,OUTFIL=0):
    import glob
    from matplotlib.backends.backend_pdf import PdfPages
    import ccd_exerc_plots as cep

    # Set output
    #if outfil == 0: 
    #    outfil = 'tmp.pdf'

    # Find all the files
    files = glob.glob(path+'*.FIT')
    #files = ['/Users/xavier/Class/PH136/Data/Tests/CCD/dark_test.00001011.DARK.FIT', '/Users/xavier/Class/PH136/Data/Tests/CCD/dark_test.00001012.DARK.FIT'] 

    pp = PdfPages('tmp.pdf')
    #if len(files) == 20:
    # Loop
    for ff in files:
        print 'File = ', ff
        cep.hist_plot(pp,ff)

    # Close plot
    pp.close()
            
    
# Read header, print exptime (one file)
def show_exptime(fil):
    # Error checking?

    # Read FITS file
    hdulist=fits.open(fil)
    head = hdulist[0].header
    hdulist.close

    # Read header
    expt = hdulist[0].header['EXPTIME']
    print fil
    print 'EXPTIME=', expt, ' sec'

#  Loop through all .FIT files in the path
#  cee.show_all_fit('/Users/xavier/Class/PH136/Data/Tests/CCD/')
def show_all_fit(path):
    import glob
    import ccd_exerc_exptime as cee
    # Find all the files
    files = glob.glob(path+'/*.FIT')
    if len(files) == 20:
        # Loop
        for ff in files:
            cee.show_exptime(ff)
            
    
#### ###############################
#  Loop through all .FIT files in the path
def stack_img(files, outfil):

    from astropy.io.fits import getdata

    # Create the image array
    flg = 0
    for ff in files:
        print 'File = ', ff
        arr,head = getdata(ff,0,header=True)

        # Save
        if flg == 0:
            flg=1
            all_arr = [arr] 
        else:
            all_arr.append(arr)
            #print 'shape', array(all_bias).shape
    
    # Stack away!
    all_arr = array(all_arr)
    final_arr = median(all_arr,0)
    #print 'shape', final_bias.shape

    # Write
    print 'Writing stacked frame to ', outfil
    fits.writeto(outfil, final_arr, clobber=True)

    # Return
    return final_arr


# #############################
#  Main driver
#  import ccd_exerc as cec
#  cec.make_stacks('../Lists/bias_files.lst', '../Lists/dark_files.lst')
def make_stacks(biaslist,darklist,path=0,biasout=0,darkout=0):

    import ccd_exerc as cec
    import pdb

    # Set path
    if path == 0: 
        path=''

    # Set bias output
    if biasout == 0: 
        biasout='final_bias.fits'

    # Set dark output
    if darkout == 0: 
        darkout='final_dark.fits'

    # ######
    # Bias first

    # Read the list of bias frames (and strip /n)
    with open(biaslist) as f:
        biasf = [line.rstrip() for line in f]
        #  This would work too: biasf = open(biaslist).read().splitlines()

    # Stack
    final_bias = cec.stack_img(biasf, biasout)

    # Stats
    bias_lvl = median(final_bias)
    print 'Bias Level = ', bias_lvl
    read_noise = std(final_bias)
    print 'Read Noise = ', read_noise
    
    # Open PDF for plots
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages('cec_plots.pdf')

    # Plot histogram
    import ccd_exerc_plots as cep
    cep.hist_plot(pp2,biasout)

    # Plot image
    cep.plot_img(pp,biasout)


    # ######
    # Now the darks
    # Read the list of dark frames (and strip /n)
    darkf = open(darklist).read().splitlines()
    #with open(darklist) as f:
    #    darkf = [line.rstrip() for line in f]

    # Stack
    final_dark = cec.stack_img(darkf, darkout)

    # Stats
    dark_lvl = median(final_dark) - bias_lvl
    print 'Dark Level = ', dark_lvl
    from astropy.io.fits import getheader
    head = getheader(darkf[0])
    expt = head['EXPTIME']
    dark_curr = dark_lvl / expt
    print 'Dark Current (counts/s) = ', dark_curr
    
    # Plot histogram
    cep.hist_plot(pp,darkout)

    # Plot image
    cep.plot_img(pp,darkout)


    # Close plot
    pp.close()
    print 'All done with make_stacks'
