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

from astropy import units as u
from astropy.units import Unit
from astropy.io import ascii, fits
from astropy.table import QTable, Table, Column
from astropy.coordinates import SkyCoord
from astropy import constants as const

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
def field_targets(field, outfil=None):

    # Init
    fcoord = SkyCoord(ra=field[1],dec=field[2])
    if outfil is None:
        outfil = xcasu.get_filename(field,'TARG_FIG')
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
    plt.figure(figsize=(7, 5))
    plt.clf()
    gs = gridspec.GridSpec(1,2)

    plt.suptitle('Hectospec Targets (from SDSS imaging)',fontsize=19.)
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
            extent=(-imsize/2., imsize/2, -imsize/2.,imsize/2))

        # Targets
        hecto_targ = np.where(lfield.targets['INSTR'] == 'HECTOSPEC')[0]
        ddec = all_sep[hecto_targ]*np.cos(all_pa[hecto_targ])
        dra = -1.*all_sep[hecto_targ]*np.sin(all_pa[hecto_targ])
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
        ax_hecto.set_xlim(-1*imsize/2., imsize/2.)
        ax_hecto.set_ylim(-1*imsize/2., imsize/2.)




    #ax_full.set_xlim(wvmnx)
    #ax_full.set_ylim(-0.05*perc[1], 1.1*perc[1])

    # Label
    #ax_full.set_xlabel('Wavelength (AA)')
    #ax_full.set_ylabel('Relative Flux')
    #ax_full.text(0.05, 0.8, ifile[i0:i1], transform=ax_full.transAxes, 
    #    color='black',
    #    size='smaller', ha='left', va='center', bbox={'facecolor':'white'})
    # Font size
    #xputils.set_fontsize(ax_full,fsz)

    plt.tight_layout(pad=0.2,h_pad=0.3,w_pad=0.0)
    if outfil is not None:
        pp.savefig()
        print('Wrote figure: {:s}'.format(outfil))
        pp.close()
    else: 
        plt.show()


    
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