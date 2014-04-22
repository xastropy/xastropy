# Import libraries
from numpy import *
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages

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
#  import ccd_exerc_stack as ces
#  ces.make_stacks('../Lists/bias_files.lst', '../Lists/dark_files.lst')
def make_stacks(biaslist,darklist,path=0,biasout=0,darkout=0):

    import ccd_exerc_stack as ces
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
    final_bias = ces.stack_img(biasf, biasout)

    # Stats
    bias_lvl = median(final_bias)
    print 'Bias Level = ', bias_lvl
    read_noise = std(final_bias)
    print 'Read Noise = ', read_noise
    
    # Open PDF for plots
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages('ces_plots.pdf')

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
    final_dark = ces.stack_img(darkf, darkout)

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
