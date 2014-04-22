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
            
    

