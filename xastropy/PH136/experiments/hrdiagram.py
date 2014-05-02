# Module for the HR Diagram experiment
#  Run on 2014 Apr 30 data from the Nickel

import numpy as np
import pdb
import glob
from astropy.io import fits

# Generate a simple ASCII log from the data
def mklog(file_path=None,outfil=None):

    # Output file
    if outfil == None: 
        outfil = 'simple.log'
    fpt = open(outfil,'w')

    # Data path
    if file_path == None: 
        file_path = 'Raw/'

    files = glob.glob(file_path+'/*.fits.gz')
    for ff in files:
        print ff
        # Open the file and grab the header + data
        hdulist=fits.open(ff)
        head = hdulist[0].header
        dat = hdulist[0].data
        hdulist.close
        # Generate the line
        lin = str(ff).ljust(20,' ')
        lin += str(head['OBSTYPE']).ljust(7, ' ')
        lin += ' {:6.1f}'.format(head['EXPTIME'])
        lin += '\n'
        # Write
        fpt.write(lin)

    # Close the file
    fpt.close()
        

def mkbias(file_list=None,file_path=None):
    # There is no overscan region
    #  So, we need a bias image

    # Defaults
    if file_path == None: 
        file_path = 'Raw/'

    # Files
    if file_list == None:
        # Generate them ourself
        biasfrm = 1 + np.arange(10)
        bias_fil = file_path+'d'+str(biasfrm)+'.fits'
