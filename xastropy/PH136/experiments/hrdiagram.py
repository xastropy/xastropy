# Module for the HR Diagram experiment
#  Run on 2014 Apr 30 data from the Nickel

import numpy as np
import pdb
import glob
from astropy.io import fits
import ds9

####################################
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
    files.sort(key=len) # Sort by file length
    lin = 'Index| File      | Descript | Type | RA | DEC | Exp | Filter\n'
    fpt.write(lin)
    for ff in files:
        print ff
        # Open the file and grab the header + data
        hdulist=fits.open(ff)
        head = hdulist[0].header
        dat = hdulist[0].data
        hdulist.close
        # Generate the line
        lin = str(files.index(ff)).ljust(3, ' ')+'|'
        lin += str(ff).ljust(20,' ')+'|'
        lin += str(head['OBJECT']).ljust(15, ' ')+'|'
        lin += str(head['OBSTYPE']).ljust(7, ' ')+'|'
        lin += str(head['RA']).strip()+'|'
        lin += (str(head['DEC']).strip()).rjust(11,'+')+'|'
        lin += ' {:6.1f}'.format(head['EXPTIME'])+'|'
        lin += ' '+str(head['FILTNAM']).ljust(2, ' ')
        lin += '\n'
        # Write
        fpt.write(lin)

    # Close the file
    fpt.close()
        

####################################
# Generate the Bias frame and also Estimate+Record Read Noise
def mkbias(file_list=None,file_path=None,outfil=None):
    # There is no overscan region (yes there is!)
    #  So, we need a bias image
    import xastropy.PH136.experiments.hrdiagram as hrd
    from astropy.io.fits import getdata
    from astropy.stats import sigma_clip 

    # Defaults
    if file_path == None: 
        file_path = 'Raw/'
    
    if outfil == None: 
        outfil = 'Bias.fits'

    # Files
    bias_fil = None
    if file_list == None:
        # Generate them ourself
        biasfrm = 2 + np.arange(10)
        bias_fil = hrd.mk_file_list(biasfrm, file_path=file_path)

    # Read Noise
    arr,head = getdata(str(bias_fil[0]),0,header=True)
    clip_arr = sigma_clip(arr, 2.5, None)
    rnoise = np.std(clip_arr,dtype=np.float64)
    print 'Read Noise = ', rnoise, ' counts'

    #pdb.set_trace()
    #dpt = ds9.ds9()
    #dpt.set_np2arr(arr)

    # Stack the frames
    img = hrd.stack_img(bias_fil)

    # Write
    head.update('RNOISE', rnoise, 'READ NOISE')
    fits.writeto(outfil, img, head, clobber=True)
    return

####################################
# Generate the Sky flats for each Filter
def mk_skyflats(file_list=None,file_path=None, bias_fil=None):
    import xastropy.PH136.experiments.hrdiagram as hrd
    from astropy.io.fits import getdata

    # Frames
    B_sky = 32 + np.arange(5)
    V_sky = 37 + np.arange(5)
    R_sky = 43 + np.arange(5)
    all_sky = [B_sky, V_sky, R_sky]

    # Outfiles
    outfil = ['Sky_B.fits', 'Sky_V.fits', 'Sky_R.fits']
    filters = ['B','V','R']

    # Bias
    if bias_fil == None:
        bias_fil = 'Bias.fits'
    bias_img,bias_head = getdata(bias_fil,0,header=True)
    

    # Loop on Filters
    for ff in filters:
        # Index
        idx = filters.index(ff)

        # Generate file names
        files= hrd.mk_file_list(all_sky[idx])

        # Stack with scaling
        img = hrd.stack_img(files, bias_img=bias_img, norm=True)

        # Trim
        trim_img = hrd.trim_img(img)

        # Write
        print 'Sky Flats: Writing ', outfil[idx]
        fits.writeto(outfil[idx], trim_img, clobber=True)

    print 'Sky Flats: All done'
    return
        
    

######
# Stack a given set of images
######
def stack_img(file_list, outfil=None, norm=False, bias_img=None):
    from astropy.io.fits import getdata

    # file_list is a List of filenames
    flg = 0
    for ff in file_list:
        print 'Reading ', ff
        arr,head = getdata(ff,0,header=True)

        # Bias subtract?
        if bias_img != None:
            if arr.shape != bias_img.shape:
                raise NameError('stack_img: Bad shapes!')
            arr = arr - bias_img

        # Normalize (for flats)
        if norm:
            norm_val = np.median(arr)
        else: 
            norm_val = 1.
        arr = arr/norm_val

        # Save
        if flg == 0:
            flg=1
            all_arr = [arr] 
        else:
            all_arr.append(arr)
            #print 'shape', array(all_bias).shape
    
    # Median Stack 
    all_arr = np.array(all_arr)
    final_arr = np.median(all_arr,0)

    # Write
    if outfil != None:
        print 'Writing stacked frame to ', outfil
        fits.writeto(outfil, final_arr, clobber=True)

    return final_arr

# Generate a list of filenames from the frame list
def mk_file_list(frames, file_path=None):

    # Path
    if file_path == None: 
        file_path = 'Raw/'

    all_fil = None
    for ii in frames: 
        fil = file_path+'d'+str(ii)+'.fits.gz'
        if all_fil == None: 
            all_fil = [fil]
        else: 
            all_fil.append(fil)

    return all_fil

# Trim the images down (remove overscan, etc.)
def trim_img(img):

    # Bottom row and overscan
    trim = img[1:,:1021]
    #pdb.set_trace()

    return trim

####################################
# Process M67 images
def proc_m67(file_path=None,outdir=None):

    from astropy.coordinates import ICRS, Distance
    from astropy.io import ascii
    from astropy import units as u
    from astropy.coordinates import angles

    # Defaults
    if file_path == None: 
        file_path = 'Raw/'
    if outdir == None: 
        outdir = 'Science/'

    # Read Log
    data = ascii.read('simple.log',delimiter='|')
    nfil = len(data)
    all_coord = ICRS(ra=data['RA'], dec=data['DEC'], unit=(u.hour,u.degree))

    # M67 coords
    m67_ra = '08:54:24'
    m67_dec = '+11:49:00'
    m67 = ICRS(m67_ra, m67_dec, unit=(u.hour,u.degree))

    # Find all separations
    sep = (m67.separation(all_coord)).degree
    im67, = np.where( sep < 1. )

    # Defaults
    return
