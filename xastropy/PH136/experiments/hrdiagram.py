# Module for the HR Diagram experiment
#  Run on 2014 Apr 30 data from the Nickel

import numpy as np
import pdb
import glob
from astropy.io import fits
import ds9

####################################
# Generate a simple ASCII log from the data
def mk_log(file_path=None,outfil=None):

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
def mk_bias(file_list=None,file_path=None,outfil=None):
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
        trim_img = hrd.trimflip_img(img)

        # Deal with zeros
        zro = np.where( trim_img == 0.)
        trim_img[zro] = 1.

        # Write
        print 'Sky Flats: Writing ', outfil[idx]
        fits.writeto(outfil[idx], trim_img, clobber=True)

    print 'Sky Flats: All done'
    return
        
    

####################################
# Stack a given set of images
######
def stack_img(file_list, outfil=None, norm=False, bias_img=None):
    from astropy.io.fits import getdata

    # file_list is a List of filenames
    all_arr = []    
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
        all_arr.append(arr)
    
    # Median Stack 
    all_arr = np.array(all_arr)
    final_arr = np.median(all_arr,0)


    # Write
    if outfil != None:
        print 'Writing stacked frame to ', outfil
        fits.writeto(outfil, final_arr, clobber=True)

    return final_arr

####################################
# Generate a list of filenames from the frame list
def mk_file_list(frames, file_path=None):

    # Path
    if file_path == None: 
        file_path = 'Raw/'

    all_fil = []
    for ii in frames: 
        fil = file_path+'d'+str(ii)+'.fits.gz'
        all_fil.append(fil)

    return all_fil

####################################
# Trim and flip the images down (remove overscan, etc.)
def trimflip_img(img):

    # Bottom row and overscan
    trim = img[1:,:1000]
    #pdb.set_trace()

    # flip
    #pdb.set_trace()
    newimg = np.flipud(trim)

    return newimg

####################################
# Process M67 images
def proc_m67(file_path=None,outdir=None, bias_fil=None):

    from astropy.coordinates import ICRS 
    from astropy.io import ascii
    from astropy import units as u
    from astropy.io.fits import getdata
    import xastropy.PH136.experiments.hrdiagram as hrd

    # Defaults
    if file_path == None: 
        file_path = 'Raw/'
    if outdir == None: 
        outdir = 'Science/'

    # Bias frame
    if bias_fil == None:
        bias_fil = 'Bias.fits'
    bias_img,bias_head = getdata(bias_fil,0,header=True)
    

    # Read Log
    data = ascii.read('simple.log',delimiter='|')
    nfil = len(data)
    all_coord = ICRS(ra=data['RA'], dec=data['DEC'], unit=(u.hour,u.degree))

    # M67 coords
    m67_rac = '08:54:24'
    m67_decc = '+11:49:00'
    m67_c = ICRS(m67_rac, m67_decc, unit=(u.hour,u.degree))

    # Find all M67
    sep = (m67_c.separation(all_coord)).degree
    im67, = np.where( sep < 1. ) # 1 degree
    m67 = data[im67]

    # 5 positions
    m67_ra = ['08:52:02.2', '08:52:15.3', '08:51:49.9', '08:51:50.0', '08:52:16.2']
    m67_dec = ['+11:52:41.0', '+11:55:51.0', '+11:55:53.0', '+11:49:38.0', '+11:49:40.0']
    m67_pointings = ICRS(ra=m67_ra, dec=m67_dec, unit=(u.hour, u.degree))

    # Filters
    all_filt=np.array(m67['Filter'])
    filters,ifilt = np.unique(all_filt,return_index=True)

    # Loop on Filterse
    all_fil = []
    for ff in filters:

        # Load Sky frame
        skyfil = 'Sky_'+ff+'.fits'
        sky_img,head = getdata(skyfil,0,header=True)

        # Images
        idx = np.where(m67['Filter'] == ff) 

        # Loop on images
        for kk in np.concatenate(idx,axis=0):
            # Read
            img,head = getdata(m67[kk]['File'],0,header=True)

            # Bias subtract
            img = img - bias_img

            # Trim
            timg = hrd.trimflip_img(img)

            # Flat field
            timg = timg / sky_img

            # Normalize by exposure
            timg = timg / m67[kk]['Exp']

            # Filename
            coord = ICRS(head['RA'], head['DEC'], unit=(u.hour,u.degree))
            sep = (coord.separation(m67_pointings)).degree
            ipos = np.argmin(sep)
            outfil = outdir+'M67_C'+str(ipos)+'_t'+str(int(m67[kk]['Exp']))+'_'+ff+'.fits'

            # Check for duplicate
            flg_skip = 0
            mt = [i for i in range(len(all_fil)) if all_fil[i] == outfil] 
            if len(mt) > 0:
                print 'Duplicate image', outfil
                print 'Skipping...'
                continue
            all_fil.append(outfil)
            

            # Write
            print 'Writing ', outfil
            fits.writeto(outfil, timg, clobber=True)
            
    return

####################################
# Process SA 104 images
def proc_sa104(file_path=None,outdir=None, bias_fil=None):

    from astropy.coordinates import ICRS 
    from astropy.io import ascii
    from astropy import units as u
    from astropy.io.fits import getdata
    import xastropy.PH136.experiments.hrdiagram as hrd

    # Defaults
    if file_path == None: 
        file_path = 'Raw/'
    if outdir == None: 
        outdir = 'Std/'

    # Bias frame
    if bias_fil == None:
        bias_fil = 'Bias.fits'
    bias_img,bias_head = getdata(bias_fil,0,header=True)
    

    # Read Log
    data = ascii.read('simple.log',delimiter='|')
    nfil = len(data)
    all_coord = ICRS(ra=data['RA'], dec=data['DEC'], unit=(u.hour,u.degree))

    # M67 coords
    sa104_rac = '12:43:44.3'
    sa104_decc = '-00:29:40.0'
    sa104_c = ICRS(sa104_rac, sa104_decc, unit=(u.hour,u.degree))

    # Find all SA 104
    sep = (sa104_c.separation(all_coord)).degree
    isa104, = np.where( sep < 1. ) # 1 degree
    sa104 = data[isa104]

    # Filters
    all_filt=np.array(sa104['Filter'])
    filters,ifilt = np.unique(all_filt,return_index=True)

    # Loop on Filterse
    all_fil = []
    for ff in filters:

        # Load Sky frame
        skyfil = 'Sky_'+ff+'.fits'
        sky_img,head = getdata(skyfil,0,header=True)

        # Images
        idx = np.where(sa104['Filter'] == ff) 

        # Loop on images
        for kk in np.concatenate(idx,axis=0):
            # Read
            img,head = getdata(sa104[kk]['File'],0,header=True)

            # Bias subtract
            img = img - bias_img

            # Trim
            timg = hrd.trimflip_img(img)

            # Flat field
            timg = timg / sky_img

            # Normalize by exposure
            timg = timg / sa104[kk]['Exp']

            # Filename
            outfil = outdir+'SA104_t'+str(int(sa104[kk]['Exp']))+'_'+ff+'.fits'

            # Check for duplicate
            flg_skip = 0
            mt = [i for i in range(len(all_fil)) if all_fil[i] == outfil] 
            if len(mt) > 0:
                print 'Duplicate image', outfil
                print 'Skipping...'
                continue
            all_fil.append(outfil)
            

            # Write
            print 'Writing ', outfil
            fits.writeto(outfil, timg, clobber=True)
            
    return

####################################
# Process SA 104 images
#def phot_sa104(std_path=None,bias_fil=None):

    # SA 104 stars
#    sa104 = { 

####################################
# Process M67 images
#def sex_all(file_path=None,outdir=None, bias_fil=None):


class star_data:
    """A simple class to track SED data for an object"""

    __slots__ = ['RA', 'DEC', 'zqso', 'FUV', 'NUV', 'umag', 'gmag']

    def __init__(self, RA, DEC,  zqso, FUV, NUV, umag, gmag):
        self.RA = RA
        self.DEC = DEC
        self.zqso = zqso
        self.FUV = FUV
        self.NUV = NUV
        self.umag = umag
        self.gmag = gmag

class Landolt_data:
    """A simple class for Landolt data"""

    __slots__ = ['RA', 'DEC', 'V', 'B-V', 'U-B', 'V-R', 'R-I', 'V-I', 'n', 'm', 'xpix', 'ypix']
                 

    def __init__(self, RA, DEC, V, BV, UB, VR, RI, VI, n=0, m=0, xpix=0., ypix=0.):
        self.RA = RA
        self.DEC = DEC
        self.V = V
        self.BV = BV
        self.UB = UB
        self.VR = VR
        self.RI = RI
        self.VI = VI
        self.n = n
        self.m = m
        self.xpix = xpix
        self.ypix = ypix

