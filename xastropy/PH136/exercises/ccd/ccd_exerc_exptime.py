# Import libraries
from numpy import *
from astropy.io import fits

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
            
    
