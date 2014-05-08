"""Python module for the HR diagram experiment of PH136
Run on 2014 Apr 30 data from the Nickel
"""
# Module for the HR Diagram experiment

import numpy as np
import pdb
import glob
from astropy.io import fits
import ds9

####################################
# CLASSES
###########
class Landolt_data:
    """A simple class for Landolt data"""

    __slots__ = ['Name', 'RA', 'DEC', 'V', 'B-V', 'U-B', 'V-R', 'R-I', 'V-I', 'n', 'm']
                 

    def __init__(self, Name, RA, DEC, V, BV, UB, VR, RI, VI, n=0, m=0):
        self.Name = Name
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

    # Push back the apparent magnitude
    def mAB(self, Filter):
        # Data
        allf = ['U','B','V','R','I']
        allmAB = [ self.UB+self.BV+self.V, 
                   self.BV+self.V, 
                   self.V,
                   self.V-self.VR,
                   self.V-self.VI]
        # Index
        #pdb.set_trace()
        try:
            idx = allf.index(Filter)
        except: 
            pdb.set_trace()
        return allmAB[idx]

#
class standard_star:
    """A simple class for a standard star and its analysis"""

    __slots__ = ['Name', 'RA', 'DEC', 'xpix', 'ypix', 'Filter', 'mI', 'sigmI', 'mAB']
                 

    def __init__(self, Name, xpix=0., ypix=0., Filter='', mI=0., sigmI=0., mAB=0.):
        self.Name = Name # Landolt name (Field_star)
        self.xpix = xpix # Detector coordinates
        self.ypix = ypix # Detector coordinates
        self.mI = mI
        self.sigmI = sigmI
        self.Filter = Filter
        self.mAB = mAB
        #self.sigmAB = sigmAB
        #self.RA = RA
        #self.DEC = DEC

    def __str__(self): # Return info
        lin= 'Standard Star: {} at pixel ({:.1f},{:.1f})\n'.format(self.Name,
                                                                   self.xpix, 
                                                                   self.ypix)
        lin+='Instrument magnitude: {:.2f}, {:.2f}'.format(self.mI,self.sigmI)
        lin+=' for filter {}'.format(self.Filter)
        if self.mAB > 0.:
            lin+=' Landolt magnitude = {:.2f}'.format(self.mAB)
        return lin

    # Center the star on an image
    def centroid(self,img,win_xy=None, No_Recenter=None, mask=None, 
                 weight=None, Silent=None): 

        # Region to center about
        if win_xy != None: 
            win_xy = (20, 20)

        # Cut image
        shpe = img.shape
        i0 = np.max( [int(self.xpix-win_xy[0]), 0] )
        i1 = np.min( [int(self.xpix+win_xy[0]), shpe[0]] )

        i2 = np.max( [int(self.ypix-win_xy[1]), 0] )
        i3 = np.min( [int(self.ypix+win_xy[1]), shpe[1]] )
        tmp_img = img[i0:i1, i2:i3]

        # Centroid
        import xastropy.phot.ian_phot as iph
        reload(iph)
        #pdb.set_trace()
        x0,y0 = iph.centroid(tmp_img,mask=mask, w=weight)

        # Offset based on windowing
        x0 = x0 + i0
        y0 = y0 + i2

        if Silent == None:
            print 'Original position: ({:.1f},{:.1f})'.format(self.xpix,self.ypix)
            print 'New position: ({:.1f},{:.1f})'.format(x0,y0)

        # Save?
        if No_Recenter == None: 
            self.xpix = x0
            self.ypix = y0

    # Turn into a FITS table

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
    from astropy.io.fits import Column
    from astropy.io import fits 
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
            fits.writeto(outfil, timg, head, clobber=True)
            
    return

####################################
# Perform photometry on SA104 and calculate ZP
def phot_sa104(outfil=None):

    import xastropy.PH136.experiments.hrdiagram as hrd
    from astropy.io.fits import getdata
    import xastropy.phot.ian_phot as iph
    from astropy.stats import sigma_clip 
    reload(iph)

    # Outfil
    if outfil == None:
        outfil = 'Std/ZP_SA104.fits'
    # SA 104 stars (470, 350, 461)
    sa104 = [ hrd.Landolt_data('104_470', '12:43:22.314', '-00:29:52.83', 14.310, 0.732, 
                               0.101, 0.295, 0.356, 0.649), 
              hrd.Landolt_data('104_350', '12:43:14.204', '-00:33:20.54', 13.634, 0.673, 
                               0.165, 0.383, 0.353, 0.736), 
              hrd.Landolt_data('104_461', '12:43:06.031', '-00:32:18.01', 9.705, 0.476, 
                               -0.035, 0.288, 0.289, 0.579)]

    # Set (approximate) pixel values in image 
    sa104_xypix = np.array( ((58.1, 976.2), 
                    (308.0, 427.), 
                    (645.9, 583.)))
    gdstar = [1,2] # First one is too close to the edge (I think)

    # Aperture:  Nickel binned 2x2 = 0.
    arcpix = 0.368 # arcsec/pix
    aper = (np.array( (7., 15., 25.) )/arcpix).tolist()

    # Grab files
    std_files = glob.glob('Std/SA104*fits')
    
    std_stars = []
    afilt = []
    # Loop on images
    for ff in std_files:
        # Read
        img,head = getdata(ff,0,header=True)
        filt = str(head['FILTNAM']).strip()
        afilt.append(filt)
        # Loop on stars
        for ii in gdstar: 
            # Construct
            std_star = hrd.standard_star(sa104[ii].Name, 
                                         sa104_xypix[ii,0], sa104_xypix[ii,1],
                                         Filter=filt, mAB=sa104[ii].mAB(filt))
            #pdb.set_trace()
            # Centroid
            std_star.centroid(img,win_xy=(20,20))

            # Photometry
            iphot =  iph.aperphot(ff, pos=[std_star.xpix,std_star.ypix], dap=aper)

            # Push into our data
            std_star.mI = -2.5*np.log10(iphot.phot)
            std_star.sigmI = 2.5 * iphot.ephot / iphot.phot / np.log(10.)

            # Add to lists
            std_stars.append(std_star)
        
    # Calculate the Zero point
    all_filt=np.array(afilt)
    filters,ifilt = np.unique(all_filt,return_index=True)
    nfilt = len(filters)
    mZP = np.zeros(nfilt)
    sig_mZP = np.zeros(nfilt)

    # Loop on filter
    idx = -1
    for ff in filters:
        zp = []
        idx = idx + 1
        # Loop on standar star obs
        for std in std_stars:
            if std.Filter == ff:
                zp.append( std.mAB-std.mI )
        if len(zp) == 0:
            pdb.set_trace()
        # Combine
        clipzp = sigma_clip(zp,2.5,None)
        mZP[idx] = np.mean(clipzp)
        sig_mZP[idx] = np.std(clipzp)

    # Save as FITS table
    c1 = fits.Column(name='Filter',format='1A',array=filters)
    c2 = fits.Column(name='ZP',format='E',array=mZP)
    c3 = fits.Column(name='sig_ZP',format='E',array=sig_mZP)
    tbhdu = fits.new_table([c1, c2, c3])
    tbhdu.writeto(outfil, clobber=True)
    print 'phot_sa104: Wrote ', outfil

    #pdb.set_trace()

    # 'Raw' Data too?
    import pickle
    rawfil = 'Std/Raw_SA104.pkl'
    f = open(rawfil,'wb')
    pickle.dump(std_stars,f)
    f.close()

    #
    print 'phot_sa104: All done!'

    # Loop
    #for obj in sa104:
        # Refine centroid
    return

####################################
# Process M67 images
def sex_m67():

    from subprocess import Popen, PIPE

    # Get the list
    m67_files = glob.glob('Science/M67*.fits')

    # Loop away on M67 images
    for ff in m67_files:
        # Run SExtractor
        p = Popen(['sex', ff, '-c', 'Sex/m67_config.sex']).wait() #, stdout=PIPE,stderr=PIPE)

        # New files into
        newdat = 'Sex/'+'sex_'+ff[8:-5]+'.dat'
        newseg = 'Sex/'+'seg_'+ff[8:]
        print 'Writing: ', newdat, newseg

        # Push them
        p2 = Popen(['mv', 'm67.dat', newdat]).wait()#, stdout=PIPE,stderr=PIPE)
        p3 = Popen(['mv', 'check.fits', newseg]).wait()#, stdout=PIPE,stderr=PIPE)
        #pdb.set_trace()
    
    print 'sex_m67: All Done'
    return

####################################
# Generate the M67 catalog
#def cat_m67():

    # Add filter column
    # Add field column
