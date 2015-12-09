# Module survey figures

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import glob, copy, os, sys
import json


import matplotlib as mpl
mpl.rcParams['font.family'] = 'stixgeneral'
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

import aplpy

from astropy import units as u
from astropy.units import Unit
from astropy.io import ascii, fits
from astropy.table import QTable, Table, Column
from astropy.coordinates import SkyCoord
from astropy import constants as const
from astropy import wcs

from linetools.spectra.xspectrum1d import XSpectrum1D

from xastropy.obs import x_getsdssimg as xosdss
from xastropy.obs import radec as xra
from xastropy.plotting import utils as xputils
from xastropy.xutils import xdebug as xdb
from xastropy.casbah import utils as xcasu
from xastropy.casbah import load_casbah as xcasl

# Local 
#sys.path.append(os.path.abspath("./py"))
#import qpq_spec as qpqs

#### ########################## #########################
def hectospec_targets(field, outfil=None):

    # Init
    fcoord = SkyCoord(ra=field[1],dec=field[2])
    if outfil is None:
        outfil = xcasu.get_filename(field,'HECTO_TARG_FIG')
    # Load field
    lfield = xcasl.load_field(field)
    targ_coord = SkyCoord(ra=lfield.targets['TARG_RA']*u.deg,
        dec=lfield.targets['TARG_DEC']*u.deg)
    all_pa = fcoord.position_angle(targ_coord)
    all_sep = fcoord.separation(targ_coord).to('arcmin')

    # Start the plot
    if outfil is not None: 
        pp = PdfPages(outfil)

    # Targets only
    plt.figure(figsize=(8, 4.5))
    plt.clf()
    gs = gridspec.GridSpec(1,2)

    plt.suptitle('{:s}: Hectospec Targets (from SDSS imaging)'.format(field[0])
        ,fontsize=19.)
    ##
    # Hectospec first
    for tt in range(2):
        ax_hecto = plt.subplot(gs[tt])

        # Read SDSS Image
        if tt == 0:
            imsize=60. # arcmin
        else:
            imsize=10. # arcmin

        #Configs
        if tt == 0:
            hecto_obs = lfield.observing[np.where(
                lfield.observing['INSTR']=='HECTOSPEC')[0]]
            unimsk = np.unique(np.array(hecto_obs['MASK_NAME']))
            for msk in unimsk:
                mt = np.where(hecto_obs['MASK_NAME'] == msk)[0]
                if not hecto_obs['DATE_OBS'].mask[mt[0]]:
                    # RA,DEC
                    rd_off, PA = xra.offsets(fcoord, 
                        (hecto_obs['MASK_RA'][mt[0]],
                        hecto_obs['MASK_DEC'][mt[0]]),verbose=False)
                    # Plot
                    circ = plt.Circle((rd_off[0].to('arcmin').value,
                        rd_off[1].to('arcmin').value), 30., 
                        color='y', fill=False, alpha=0.5)
                    ax_hecto.add_artist(circ)

        # Plot SDSS image
        sdss_img, _ = xosdss.getimg(field[1],field[2], imsize)
        ax_hecto.imshow(sdss_img,aspect='equal',
            extent=(imsize/2., -imsize/2, -imsize/2.,imsize/2))

        # Targets
        hecto_targ = np.where(lfield.targets['INSTR'] == 'HECTOSPEC')[0]
        ddec = all_sep[hecto_targ]*np.cos(all_pa[hecto_targ])
        dra = all_sep[hecto_targ]*np.sin(all_pa[hecto_targ])
        #xdb.set_trace()
        ax_hecto.scatter(dra,ddec, marker='o',color='gray',s=10.,
            facecolor='none',linewidth=0.3, alpha=0.5)
        # Observed
        obs_targ_tab, obs_dict, obs_idx = lfield.get_observed(imsize*u.arcmin,
            subtab=lfield.targets[hecto_targ])
        ax_hecto.scatter(dra[obs_idx],ddec[obs_idx], marker='o',color='gray',s=10.,
            facecolor='none',linewidth=0.3)

        # Labels
        ax_hecto.set_xlabel('Relative ArcMin', fontsize=13)
        ax_hecto.set_ylabel('Relative ArcMin', fontsize=13)
        ax_hecto.set_xlim(imsize/2., -imsize/2.)
        ax_hecto.set_ylim(-imsize/2., imsize/2.)

    plt.tight_layout(pad=0.2,h_pad=0.3,w_pad=0.0)
    if outfil is not None:
        pp.savefig()
        print('Wrote figure: {:s}'.format(outfil))
        pp.close()
    else: 
        plt.show()

#### ########################## #########################
def deimos_targets(field, gc=None, outfil=None, replot=True,
    reset_layers=False):
    '''
    Call with gc to speed things up
    '''

    import matplotlib.cm as cm
    import matplotlib.patches as patches
    reload(xra)

    # Init
    # Mask dimensions
    smsk_dim = (7., 16.6) # arcmin (approximate)
    smsk_ang = np.arctan(smsk_dim[0]/smsk_dim[1])
    smsk_diag = np.sqrt(smsk_dim[0]**2 + smsk_dim[1]**2)/2.
    all_smsk_ang = [-1.*smsk_ang, smsk_ang, np.pi-smsk_ang, np.pi+smsk_ang]

    # Coordinates
    fcoord = SkyCoord(ra=field[1],dec=field[2])
    if outfil is None:
        outfil = xcasu.get_filename(field,'DEIMOS_TARG_FIG')
    # Load field
    lfield = xcasl.load_field(field)
    targ_coord = SkyCoord(ra=lfield.targets['TARG_RA']*u.deg,
        dec=lfield.targets['TARG_DEC']*u.deg)
    all_pa = fcoord.position_angle(targ_coord)
    all_sep = fcoord.separation(targ_coord).to('arcmin')

    # Targets
    deimos_targ = np.where(lfield.targets['INSTR'] == 'DEIMOS')[0]
    #dra =  all_sep[deimos_targ]*np.sin(all_pa[deimos_targ])
    #ddec = all_sep[deimos_targ]*np.cos(all_pa[deimos_targ])

    # Image
    path = os.getenv('CASBAH_GALAXIES')
    img_fil = path+'/'+field[0]+'/'+lfield.targets[deimos_targ]['TARG_IMG'][0]
    gfil = glob.glob(img_fil+'*')
    if len(gfil) != 1:
        raise ValueError('Image not found! {:s}'.format(img_fil))
    else:
        img_fil = gfil[0]
    #hdu = fits.open(img_fil)
    #w = wcs.WCS(hdu[0].header,fix=True)
    #img = hdu[0].data
    '''
    # Size
    pix_qso = w.wcs_world2pix(np.array([[field[1].value,field[2].value]]), 1)
    corners = [ [0.,0.], [0.,img.shape[0]], [img.shape[1],0.], [img.shape[1],img.shape[0]]]
    rad_corners =  w.wcs_pix2world(np.array(corners), 1) # deg but no units
    mx_rad = [0.,0.]
    for rad in rad_corners:
        c_off, _ = xra.offsets(fcoord, tuple(rad))
        for ii in range(2):
            mx_rad[ii] = np.maximum(
                mx_rad[ii],np.abs(c_off[ii].to('arcmin').value))
    imsize = 2.*np.maximum(mx_rad[0],mx_rad[1])
    '''

    # Start the plot
    #if outfil is not None: 
    #    pp = PdfPages(outfil)

    # Targets only
    plt.figure(figsize=(8, 4.9))
    plt.clf()

    plt.suptitle('{:s}: DEIMOS Targets (from LBT imaging)'.format(field[0])
        ,fontsize=19.)

    # Plot the Image
    if gc is None:
        gc = aplpy.FITSFigure(img_fil)

    radius = 25./60/2. #0.01 # deg
    if replot:
        gc.recenter(field[1].value, field[2].value,radius=radius)
        # This next one is a bit expensive
        gc.show_grayscale(vmin=300., vmax=3000.,stretch='arcsinh',invert='True')

    # Masks
    msk_clrs = ['blue', 'red', 'green', 'orange']
    deimos_obs = lfield.observing[np.where(
        lfield.observing['INSTR']=='DEIMOS')[0]]
    unimsk = np.unique(np.array(deimos_obs['MASK_NAME']))
    for ii,msk in enumerate(unimsk):
        mt = np.where(deimos_obs['MASK_NAME'] == msk)[0]
        if not deimos_obs['DATE_OBS'].mask[mt[0]]:
            # RA,DEC offset of center
            mask_rad = xra.stod1((deimos_obs['MASK_RA'][mt[0]],
                deimos_obs['MASK_DEC'][mt[0]]))
            mask_pa = deimos_obs['MASK_PA'][mt[0]]*np.pi/180.
            # RA/DEC offsets of corners
            rect_x = [smsk_diag*np.sin(smsk_ang+mask_pa) for smsk_ang in all_smsk_ang]
            rect_y = [smsk_diag*np.cos(smsk_ang+mask_pa) for smsk_ang in all_smsk_ang]
            #xdb.set_trace()
            prect_y = mask_rad[1].value + np.array(rect_y+[rect_y[0]])/60.
            prect_x = mask_rad[0].value + np.array(rect_x+[rect_x[0]])/60./np.cos(prect_y*np.pi/180.)
            '''
            coord0 = SkyCoord(ra=prect_x[0]*u.deg,dec=prect_y[0]*u.deg)            
            coord2 = SkyCoord(ra=prect_x[2]*u.deg,dec=prect_y[2]*u.deg)
            rcoord0 = SkyCoord(ra=rect_x[0]*u.arcmin,dec=rect_y[0]*u.arcmin)            
            rcoord2 = SkyCoord(ra=rect_x[2]*u.arcmin,dec=rect_y[2]*u.arcmin)
            sep02 = coord0.separation(coord2).to('arcmin')
            rsep02 = rcoord0.separation(rcoord2).to('arcmin')
            print(sep02)
            print(rsep02)
            print(smsk_diag*2.)
            xdb.set_trace()
            '''
            # Reshape
            xy = np.zeros((2,5))
            xy[0,:] = prect_x
            xy[1,:] = prect_y
            #xdb.set_trace()
            gc.show_lines([xy],color=msk_clrs[ii])#,alpha=0.5)

            # Markers
            # Observed targets
            sub_targ = np.where((lfield.targets['MASK_NAME']==msk)
                & (lfield.targets['INSTR']=='DEIMOS'))[0]
            obs_targ_tab, obs_dict, obs_idx = lfield.get_observed(radius*u.deg,
                subtab=lfield.targets[sub_targ])
            ra_obs_targ = lfield.targets[sub_targ[obs_idx]]['TARG_RA']
            dec_obs_targ = lfield.targets[sub_targ[obs_idx]]['TARG_DEC']
            gc.show_markers(ra_obs_targ, dec_obs_targ, edgecolor=msk_clrs[ii],
                facecolor='none', marker='o', s=50)#, alpha=0.5)
    # QSO
    gc.show_markers(field[1].value, field[2].value, edgecolor='lightgreen', facecolor='none',
        marker='s', s=10, alpha=0.5)
    # Label
    #gc.add_label(field[1].value, field[2].value, 'Q', color='green')

    if outfil != None:
        gc.save(outfil) 

    '''
    ##
    # Loop
    for tt in range(2):
        ax_deimos = plt.subplot(gs[tt])

        if tt == 0:
            imsize=25. # arcmin
        else:
            imsize=3. # arcmin

        #Configs
        if tt == 0:
            deimos_obs = lfield.observing[np.where(
                lfield.observing['INSTR']=='DEIMOS')[0]]
            unimsk = np.unique(np.array(deimos_obs['MASK_NAME']))
            for msk in unimsk:
                mt = np.where(deimos_obs['MASK_NAME'] == msk)[0]
                if not deimos_obs['DATE_OBS'].mask[mt[0]]:
                    # RA,DEC offset
                    rd_off, PA_off = xra.offsets(fcoord, 
                        (deimos_obs['MASK_RA'][mt[0]],
                        deimos_obs['MASK_DEC'][mt[0]]),verbose=False)
                    mask_pa = deimos_obs['MASK_PA'][mt[0]]*np.pi/180.
                    # Plot
                    rect_x = [smsk_diag*np.cos(smsk_ang+mask_pa) for smsk_ang in all_smsk_ang]
                    rect_y = [smsk_diag*np.sin(smsk_ang+mask_pa) for smsk_ang in all_smsk_ang]
                    prect_x = rd_off[0].to('arcmin').value + np.array(rect_x+[rect_x[0]]) # Extra for plotting
                    prect_y = rd_off[1].to('arcmin').value + np.array(rect_y+[rect_y[0]])
                    #xdb.set_trace()
                    ax_deimos.plot(prect_x,prect_y,'y-',alpha=0.5)

        # Plot Target image
        ax_deimos.imshow(img,aspect='equal',
            extent=(-imsize/2., imsize/2, -imsize/2.,imsize/2))

        #xdb.set_trace()
        ax_deimos.scatter(dra,ddec, marker='o',color='gray',s=10.,
            facecolor='none',linewidth=0.3, alpha=0.1)
        # Observed
        obs_targ_tab, obs_dict, obs_idx = lfield.get_observed(imsize*u.arcmin,
            subtab=lfield.targets[deimos_targ])
        ax_deimos.scatter(dra[obs_idx],ddec[obs_idx], marker='o',color='gray',s=10.,
            facecolor='none',linewidth=0.3)

        # Labels
        ax_deimos.set_xlabel('Relative ArcMin', fontsize=13)
        ax_deimos.set_ylabel('Relative ArcMin', fontsize=13)
        ax_deimos.set_xlim(-1*imsize/2., imsize/2.)
        ax_deimos.set_ylim(-1*imsize/2., imsize/2.)

    plt.tight_layout(pad=0.2,h_pad=0.3,w_pad=0.0)
    if outfil is not None:
        pp.savefig()
        print('Wrote figure: {:s}'.format(outfil))
        pp.close()
    else: 
        plt.show()
    '''

    
# ##################################################
# ##################################################
# ##################################################
# Command line execution for testing
# ##################################################
if __name__ == '__main__':


    if len(sys.argv) == 1: # TESTING

        flg_fig = 0 
        flg_fig += 2**0  # XSpec
        #flg_fig += 2**1  # XAbsID
        #flg_fig += 2**2  # XVelPlt Gui
        #flg_fig += 2**3  # XVelPlt Gui without ID list; Also tests select wave
        #flg_fig += 2**4  # XAODM Gui
        #flg_fig += 2**5  # Fit LLS GUI

        if (flg_fig % 2**1) >= 2**0:
            compare_models()
    else: # Run something
        idv = int(sys.argv[1])  # 1 = XSpec, 2=XAbsId

        if idv == 1:
            # python py/xq100_lls 1 jxp
            fit_lls_loop()
        else:
            pass