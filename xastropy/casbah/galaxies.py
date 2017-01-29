"""
#;+ 
#; NAME:
#; build_casbah_galaxies
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for buildling CASBAH galaxy database
#;   02-Jan-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, glob, copy
import warnings
from os import makedirs
from os.path import exists
from astropy.io import fits, ascii
from astropy import units as u 
from astropy.table import Table, Column, MaskedColumn, vstack

from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck15 as cosmo

from linetools import utils as ltu

#from astropy import constants as const
from xastropy.casbah import utils as xcasbahu
from xastropy.xutils import lists as xxul
from xastropy.obs import radec as xra

from xastropy.xutils import xdebug as xdb

# SDSS
def build_sdss(field, radius=2.0*u.deg):
    """ Grab SDSS photometry and spectra for those fields in the footprint

    Includes BOSS data.

    Parameters
    ----------
    field : tuple
      (name, ra_deg, dec_deg)
    radius : Angle or Quantity, optional
    """

    # Directory
    if not os.path.exists(field[0]):
        os.makedirs(field[0])
    # Grab SDSS data + write to folder
    sdss_fil = xcasbahu.get_filename(field,'SDSS')
    print('CASBAH_SDSS: Building {:s}'.format(sdss_fil))
    print('CASBAH_SDSS: Be patient..')
    grab_sdss_spectra((field[1],field[2]),
        radius=radius, outfil=sdss_fil, maxsep=20., zmin=500./3e5)
        #outfig = os.environ.get('DROPBOX_DIR')+'/CASBAH/Galaxies/SDSS/PG1407+265_SDSS.pdf'

def build_spectra(field, obs_path=None, path='./'):
    """Top-level program to build spectra files

    Parameters
    ----------
    field : tuple
      (Name, ra, dec)
    """
    if obs_path is None:
        obs_path = os.getenv('DROPBOX_DIR')+'CASBAH_Observing/'
    """
    Hectospec
    """
    ## DEAL WITH DUPLICATES (TAKE HIGHER ZQ)
    #HDU0: wavelengths (Angstroms)
    #HDU1: sky-subtracted, variance-weighted coadded spectra (total counts)
    #HDU2: inverse variance (counts)
    #HDU3: AND bad pixel mask
    #HDU4: OR bad pixel mask
    #HDU5: Plugmap structure (fiber info)
    #HDU6: Combined sky spectra
    #HDU7: Summed (unweighted) spectra

    # Load up the data
    hecto_path = '/Galx_Spectra/Hectospec/'
    spfiles = glob.glob(obs_path+field[0]+hecto_path+'spHect-*')
    spfiles.sort()
    for spfile in spfiles:
        if 'zcat' not in spfile:  # Spectra
            hdu = fits.open(spfile)
            print('Reading {:s}'.format(spfile))
            wave = hdu[0].data
            flux = hdu[1].data
            ivar = hdu[2].data
            sig = np.zeros_like(flux)
            gd = ivar > 0.
            sig[gd] = np.sqrt(ivar[gd])
            tbl = Table(hdu[5].data)
            if 'hecto_wave' not in locals():
                hecto_wave, hecto_flux, hecto_sig, hecto_stbl = wave, flux, sig, tbl
            else:
                hecto_wave = np.concatenate((hecto_wave, wave))
                hecto_flux = np.concatenate((hecto_flux, flux))
                hecto_sig = np.concatenate((hecto_sig, sig))
                hecto_stbl = vstack([hecto_stbl,tbl])
        else:
            tmp = Table.read(spfile)  # z values
            if 'hecto_ztbl' not in locals():
                hecto_ztbl = tmp
            else:
                hecto_ztbl = vstack([hecto_ztbl,tmp])
    # Check
    if len(hecto_stbl) != len(hecto_ztbl):
        raise ValueError("Bad Hecto tables..")
    # Objects only
    gdobj = np.where(hecto_ztbl['MAG'] > 1.)[0]
    nobj = len(gdobj)
    # Check for duplicates
    idval = np.array(hecto_ztbl[gdobj]['ID']).astype(int)
    uni, counts = np.unique(idval, return_counts=True)
    if len(uni) != nobj:
        warnings.warn("Resetting duplicate ID values using the targs table")
        # Load targets file
        targ_file = xcasbahu.get_filename(field, 'TARGETS')
        targs = Table.read(targ_file,delimiter='|', format='ascii.fixed_width',
                                fill_values=[('--','0','MASK_NAME')])
        tcoord = SkyCoord(ra=targs['TARG_RA']*u.deg, dec=targs['TARG_DEC']*u.deg)
        # Loop on duplicates
        dup = np.where(counts>1)[0]
        for idup in dup:
            dobj = np.where(hecto_ztbl['ID'] == str(uni[idup]))[0]
            # Loop on objects
            for idobj in dobj:
                dcoord = SkyCoord(ra=hecto_stbl['RA'][idobj]*u.deg,
                                  dec=hecto_stbl['DEC'][idobj]*u.deg)
                # Match by RA/DEC
                mt = np.argmin(dcoord.separation(tcoord))
                # Reset ID
                #xdb.set_trace()
                print('Setting ID to {:s} from {:s}'.format(
                        str(targs['TARG_ID'][mt]), hecto_ztbl['ID'][idobj]))
                hecto_ztbl['ID'][idobj] = str(targs['TARG_ID'][mt])
    # Double check
    idval = np.array(hecto_ztbl[gdobj]['ID']).astype(int)
    uni, counts = np.unique(idval, return_counts=True)
    if len(uni) != nobj:
        xdb.set_trace()
        raise ValueError("Should not get here")

    # Generate the final Table
    hecto_spec = Table()
    hecto_spec.add_column(hecto_stbl['RA'][gdobj])
    hecto_spec.add_column(hecto_stbl['DEC'][gdobj])
    hecto_spec['RA'].unit = u.deg
    hecto_spec['DEC'].unit = u.deg
    hecto_spec.add_column(hecto_ztbl['Z'][gdobj])
    hecto_spec.add_column(hecto_ztbl['Z_ERR'][gdobj])
    hecto_spec.add_column(hecto_ztbl['ZQ'][gdobj])
    hecto_spec.add_column(hecto_ztbl['APERTURE'][gdobj])
    hecto_spec.add_column(hecto_ztbl['ID'][gdobj])  # May wish to recast as int
    hecto_spec.add_column(hecto_ztbl['MAG'][gdobj])
    hecto_spec.add_column(Column(['MMT']*nobj, name='TELESCOPE'))
    hecto_spec.add_column(Column(['Hectospec']*nobj, name='INSTRUMENT'))

    hecto_spec.add_column(Column(hecto_wave[gdobj,:], name='WAVE'))
    hecto_spec.add_column(Column(hecto_flux[gdobj,:], name='FLUX'))
    hecto_spec.add_column(Column(hecto_sig[gdobj,:], name='SIG'))
    # Write
    hectospec_file = xcasbahu.get_filename(field,'HECTOSPEC')
    hecto_spec.write(hectospec_file, overwrite=True)


def build_imaging(field, obs_path=None, path='./'):
    """Top-level program to build images

    Parameters
    ----------
    field : tuple
      (Name, ra, dec)
    obs_path : str, optional
    """
    if obs_path is None:
        obs_path = os.getenv('DROPBOX_DIR')+'CASBAH_Observing/'
    import shutil
    # DEIMOS mask image
    targ_file = xcasbahu.get_filename(field,'TARGETS')
    targets = Table.read(targ_file,delimiter='|',
        format='ascii.fixed_width', 
        fill_values=[('--','0','MASK_NAME')])
    deimos_targ = np.where(targets['INSTR'] == 'DEIMOS')[0]
    if len(deimos_targ) > 0:
        # Search for LBT image
        msk_img = targets[deimos_targ]['TARG_IMG'][0]
        img_fil = glob.glob(obs_path+field[0]+'/IMG/LBT/'+msk_img+'*')
        if len(img_fil) == 1:
            # Copy
            path = os.getenv('CASBAH_GALAXIES')
            new_img_file = xcasbahu.get_filename(field, 'IMAGING', orig_file=img_fil[0], rename_deimos=True)
            # Rename target imaging
            targets['TARG_IMG'][np.where(deimos_targ)[0]] = new_img_file[new_img_file.rfind('/')+1:]
            xdb.set_trace()
            if 'gz' in img_fil[0]:
                new_img_file += '.gz'
            shutil.copy2(img_fil[0], new_img_file)
            print('Copied {:s} to {:s}'.format(img_fil[0], new_img_file))
        else:
            raise ValueError('Need to provide the image! {:s}'.format(
                field[0]+'/IMG/LBT/'+msk_img))

def build_targets(field, obs_path=None, path='./'):
    """Top-level program to build target info

    Parameters
    ----------
    field : tuple
      (Name, ra, dec)
    """
    if obs_path is None:
        obs_path = os.getenv('DROPBOX_DIR')+'CASBAH_Observing/'

    targ_list = []
    # MMT
    hecto_masks, hecto_obs, hecto_targs = hecto_targets(field, obs_path)
    if hecto_targs is not None:
        targ_list = targ_list + [hecto_targs]

    # DEIMOS
    deimos_sex, deimos_masks, deimos_obs, deimos_targs = deimos_targets(field, obs_path)
    if deimos_masks is None:
        deimos_masks = []
        deimos_obs = []
    else:
        targ_list = [deimos_sex] + targ_list

    # COLLATE
    all_masks = deimos_masks + hecto_masks
    all_obs = deimos_obs + hecto_obs
    all_sex = vstack(targ_list, join_type='inner')  # Use vstack when needed

    # Out path
    outpath = xcasbahu.get_filename(field, 'FIELD_PATH')
    if not exists(outpath):
        makedirs(outpath)

    # Generate Target table
    targ_file = xcasbahu.get_filename(field, 'TARGETS')
    cut_sex = all_sex[['TARG_RA','TARG_DEC','EPOCH','TARG_ID',
        'TARG_MAG','TARG_IMG','INSTR','MASK_NAME']]
    #cut_sex.write(targ_file,overwrite=True)
    # Rename DEIMOS mask image
    deimos_targ = np.where(targets['INSTR'] == 'DEIMOS')[0]
    cut_sex.write(targ_file,format='ascii.fixed_width',delimiter='|')
    print('Wrote file {:s}'.format(targ_file))

    # Generate MULTI_OBJ file
    multi_file = xcasbahu.get_filename(field,'MULTI_OBJ')
    tab_names=('INSTR', 'MASK_NAME', 'MASK_RA', 'MASK_DEC', 'MASK_EPOCH', 
        'MASK_PA', 'DATE_OBS', 'DISPERSER', 'TEXP', 'CONDITIONS')
    dtypes=('S12', 'S25', 'S12', 'S13', 'f4',
        'f4', 'S11', 'S10', 'f4', 'S25')
    multi_tab = Table(names=tab_names, dtype=dtypes, masked=True)
    for kk,mask in enumerate(all_masks):
        # Loop on observations
        if all_obs[kk] is not None:
            for obs in all_obs[kk]:
                maskobs = copy.deepcopy(mask)
                maskobs['DATE_OBS'] = obs['DATE']
                maskobs['TEXP'] = obs['TEXP']
                maskobs['DISPERSER'] = obs['DISPERSER']
                maskobs['CONDITIONS'] = obs['CONDITIONS']
                # Add row
                multi_tab.add_row(maskobs)
        else:
            multi_tab.add_row(mask)

    multi_tab.write(multi_file,format='ascii.fixed_width',delimiter='|')
    print('Wrote file {:s}'.format(multi_file))
    # Reading the MULTI-OBJ file
    #mtab = Table.read(ifile,delimiter='|',format='ascii.fixed_width',
             #fill_values=[('--','0','DATE_OBS','TEXP')])


# SExtrator
def parse_sex_file(field,targ_yaml_file):
    '''Parse SExtractor file for targets
    Parameters:
    ----------
    field: tuple
      (Name, ra, dec)
    targ_file: str
      Name of targetting info file.  Yaml
    '''
    qso_coord = SkyCoord(ra=field[1],dec=field[2])
    # Read yaml file
    import yaml
    stram = open(targ_yaml_file,"r")
    targ_dict = yaml.load(stram)
    # Read SExtractor output file
    i0 = targ_yaml_file.rfind('/')
    sex_file = targ_yaml_file[0:i0+1]+targ_dict['SEX_FILE']
    sex_hdu = fits.open(sex_file)
    head = sex_hdu[0].header
    sex_tab = Table(fits.open(sex_file)[1].data)

    # Coordinates
    targ_coord = SkyCoord(ra=sex_tab['ALPHA_J2000']*u.deg,
        dec=sex_tab['DELTA_J2000']*u.deg)
    sep = qso_coord.separation(targ_coord)

    # Build up list
    targ_mask = np.array([False]*len(sex_tab))

    # Star/galaxy
    try:
        star_gal_obj = np.where(sex_tab['CLASS_STAR'] < 
            targ_dict['STAR_GAL'])[0]
    except KeyError:
        pass
    else:
        print('Using STAR_GAL<{:g}'.format(targ_dict['STAR_GAL']))
        targ_mask[star_gal_obj] = True

    # Close separation?
    try:
        close_sep = np.where(sep < targ_dict['ALL_WITHIN_SEP']*u.arcsec)[0]
    except KeyError:
        pass
    else:
        print('Taking all sources within {:g}arcsec'.format(
            targ_dict['ALL_WITHIN_SEP']))
        targ_mask[close_sep] = True

    if len(np.where(targ_mask == True)) == 0:
        raise ValueError('No targets found!  Probably a problem.')

    # ####
    # Now perform cuts

    # QSO
    qsoid = np.where(sep < 1.0*u.arcsec)[0]
    if len(qsoid) > 0:
        targ_mask[qsoid] = False

    # Magnitude
    try:
        mag_cut = np.where(sex_tab[targ_dict['PHOTOM']] > 
            targ_dict['MAX_MAG'])[0]
    except KeyError:
        pass
    else:
        print('Cutting on MAG_MAX={:g}'.format(targ_dict['MAX_MAG']))
        targ_mask[mag_cut] = False
    try:
        mag_cut2 = np.where(sex_tab[targ_dict['PHOTOM']] < 
            targ_dict['MIN_MAG'])[0]
    except KeyError:
        pass
    else:
        print('Cutting on MIN_MAX={:g}'.format(targ_dict['MIN_MAG']))
        targ_mask[mag_cut2] = False

    # Separation
    try:
        sep_cut = np.where(sep > targ_dict['MAX_SEP']*u.arcmin)[0]
    except KeyError:
        pass
    else:
        print('Cutting on MAX_SEP={:g}'.format(targ_dict['MAX_SEP']))
        targ_mask[sep_cut] = False

    # Generate Table and polish
    targ_tab = sex_tab[targ_mask]
    ntarg = len(targ_tab)
    targ_tab.rename_column('ALPHA_J2000','TARG_RA')
    targ_tab.rename_column('DELTA_J2000','TARG_DEC')
    targ_tab.add_column(Column([2000.]*ntarg,name='EPOCH'))
    targ_tab.add_column(Column([head['FITSFILE']]*ntarg,name='TARG_IMG'))
    targ_tab.rename_column(targ_dict['PHOTOM'], 'TARG_MAG')
    targ_tab.rename_column('NUMBER', 'TARG_ID')

    # Return
    return targ_tab


def deimos_targets(field, obs_path, deimos_path=None):
    """Generate files related to DEIMOS deimos_targets

    Parameters:
    -----------
    field : tuple
      (Name, ra, dec)

    Returns:
    ----------
    sex_targ, all_masks, all_obs, all_masktarg
    """
    if deimos_path is None:
        deimos_path = '/Galx_Spectra/DEIMOS/'

    # Loop on Fields
    mask_path = obs_path+field[0]+deimos_path+'/Masks/'
    # SExtractor targeting
    targetting_file = glob.glob(mask_path+'*targ.yaml')
    if len(targetting_file) == 1:
        sex_targ = parse_sex_file(field, targetting_file[0])
        sex_targ.add_column(Column(['DEIMOS']*len(sex_targ),name='INSTR'))
        # Setup for mask matching
        sex_coord = SkyCoord(ra=sex_targ['TARG_RA']*u.deg, 
            dec=sex_targ['TARG_DEC']*u.deg)
        sex_msk_clms = {}
        cnames = ['MASK_NAME', 'MASK_ID']
        smsk = '--'
        msk_val = [smsk]*len(cnames)
        for kk,cname in enumerate(cnames):
            sex_msk_clms[cname] = [msk_val[kk]]*len(sex_targ)
    elif len(targetting_file) == 0:
        print('WARNING: No SExtractor info for mask path {:s}'.format(mask_path))
        print('Assuming no DEIMOS files')
        return [None]*4
        #sex_targ = None
    else:
        raise ValueError('Found multiple targ.yaml files!!')

    # Mask info
    all_masks = []
    all_masktarg = []
    all_obs = []
    files = glob.glob(mask_path+'*.out')
    for msk_file in files:
        print('Reading DEIMOS mask file: {:s}'.format(msk_file))
        # Parse 
        mask_dict, targ_tab, obs_tab = parse_deimos_mask_file(msk_file)
        # Fill up SEx file
        if sex_targ is not None:
            for targ in targ_tab:
                targ_coord = xra.to_coord((targ['RAS'], targ['DECS']))
                sep = targ_coord.separation(sex_coord)
                isep = np.argmin(sep)
                if sep[isep] > 0.5*u.arcsec:
                    raise ValueError('No match in SExtractor?!')
                else: # Fill
                    if sex_msk_clms['MASK_NAME'][isep] == smsk:
                        sex_msk_clms['MASK_NAME'][isep] = mask_dict['MASK_NAME']
                    else: # Already full 
                        sex_targ.add_row(sex_targ[isep])
                        sex_msk_clms['MASK_NAME'].append(mask_dict['MASK_NAME'])
        # Append
        all_masks.append(mask_dict)
        all_masktarg.append(targ_tab)
        all_obs.append(obs_tab)

    # Add columns to sex_targ
    for tt, cname in enumerate(cnames):
        # Mask
        mask = np.array([False]*len(sex_targ))
        bad = np.where(np.array(sex_msk_clms[cname])==msk_val[tt])[0]
        if len(bad)>0:
            mask[bad]=True
        #
        clm = MaskedColumn(sex_msk_clms[cname],name=cname,mask=mask)
        sex_targ.add_column(clm)
    # Return
    return sex_targ, all_masks, all_obs, all_masktarg


# DEIMOS
def parse_deimos_mask_file(msk_file):
    '''Parse the standard DEIMOS mask file
    '''
    # Read
    f = open(msk_file, 'r')
    lines = f.readlines()
    f.close()
    #
    flg_mask = 0
    flg_select = 1
    mask_dict = dict(INSTR='DEIMOS')
    all_targ = []
    all_obs = []
    for kk,line in enumerate(lines):
        # Searching for mask name, center
        if flg_mask==0:
            if 'Mask' in line:
                flg_mask=1
                continue
        # Mask info
        if flg_mask==1:
            prs = line.strip().split(' ')
            # 
            cnt = 0
            for jj,iprs in enumerate(prs):
                if len(iprs) == 0:
                    continue
                if cnt==0:
                    mask_dict['MASK_NAME'] = iprs
                elif cnt==1:
                    mask_dict['MASK_RA'] = iprs
                elif cnt==2:
                    mask_dict['MASK_DEC'] = iprs
                elif cnt==3:
                    mask_dict['MASK_EPOCH'] = float(iprs)
                elif cnt==4:
                    mask_dict['MASK_PA'] = float(iprs[3:])
                cnt += 1
            flg_mask = 2
            continue
        # Observation info
        if 'OBS' in line:
            prs = line.strip().split(' ')
            gdprs = [iprs for iprs in prs if len(iprs)>0]
            obs_dict = {}
            obs_dict['DATE'] = gdprs[2]
            obs_dict['TEXP'] = float(gdprs[3])
            obs_dict['DISPERSER'] = gdprs[4]
            obs_dict['CONDITIONS'] = gdprs[5]
            # 
            all_obs.append(obs_dict)


        # Done
        if 'Non' in line:
            break
        # Dummy line?
        if (len(line.strip())==0) or (line[0]=='#'):
            continue
        # Selected object
        prs = line.strip().split(' ')
        gdprs = [iprs for iprs in prs if len(iprs)>0]
        if int(gdprs[6]) < 0: # QSO or Star
            continue
        targ_dict = {}
        targ_dict['ID'] = int(gdprs[0])
        targ_dict['RAS'] = gdprs[1]
        targ_dict['DECS'] = gdprs[2]
        targ_dict['EPOCH'] = float(gdprs[3])
        # Add to list
        all_targ.append(targ_dict)

    # Generate Tables
    # Obs
    if len(all_obs) > 0:
        obs_tab = xxul.dict_list_to_table(all_obs)
        obs_tab['TEXP'].unit = u.s
    else:
        obs_tab = None
    # Targ
    targ_tab = xxul.dict_list_to_table(all_targ)

    # Return
    return mask_dict, targ_tab, obs_tab 

def hecto_targets(field, obs_path, hecto_path=None):
    '''Read files related to Hectospec targets

    Parameters:
    -----------
    field : tuple
      (Name, ra, dec)
    obs_path : str, optional
      Path to the observing tree
    hecto_path : str, optional
      Path within the file tree to Hectospec data

    Returns:
    ----------
    Target and observing info 
    '''
    if hecto_path is None:
        hecto_path = '/Galx_Spectra/Hectospec/'

    # Targets
    targ_path = obs_path+field[0]+hecto_path

    # Target file
    targ_file = glob.glob(targ_path+'*.targ')
    if len(targ_file) != 1:
        raise ValueError('Wrong number of Hectospec target files')
    else:
        targ_file = targ_file[0]

    # Read PI, program info [NOT IMPLEMENTED]
    #f = open(msk_file, 'r')
    #lines = f.readlines()
    #f.close()

    # Read target table
    tab = ascii.read(targ_file,comment='#')
    # Restrict to targets
    itarg = np.where(tab['type']=='TARGET')
    targs = tab[itarg]
    # Polish
    nrow = len(targs)
    targs.rename_column('ra','RAS')
    targs.rename_column('dec','DECS')
    targs.add_column(Column([0.]*nrow,name='TARG_RA'))
    targs.add_column(Column([0.]*nrow,name='TARG_DEC'))
    # Get RA/DEC in degrees
    for k,row in enumerate(targs):
        coord = ltu.radec_to_coord((row['RAS'], row['DECS']))
        targs[k]['TARG_RA'] = coord.ra.value
        targs[k]['TARG_DEC'] = coord.dec.value
    # ID/Mag (not always present)
    targ_coord = SkyCoord(ra=targs['TARG_RA']*u.deg, dec=targs['TARG_DEC']*u.deg)
    try:
        targs.rename_column('objid','TARG_ID')
    except KeyError:
        targs.add_column(Column([0]*nrow,name='TARG_ID'))
        targs.add_column(Column([0.]*nrow,name='TARG_MAG'))
        flg_id = 0
    else:
        flg_id = 1
        targs.rename_column('mag','TARG_MAG')
    targs.add_column(Column([0.]*nrow,name='EPOCH'))
    targs.add_column(Column(['SDSS']*nrow,name='TARG_IMG'))
    targs.add_column(Column(['HECTOSPEC']*nrow,name='INSTR'))

    targ_mask = {}
    cnames = ['MASK_NAME', 'MASK_ID']
    smsk = '--'
    msk_val = [smsk]*len(cnames)
    for kk,cname in enumerate(cnames):
        targ_mask[cname] = [msk_val[kk]]*nrow

    # Now the 'mask' files
    mask_files = glob.glob(targ_path+'*.cat')
    all_obs = []
    all_masks = []
    for mask_file in mask_files:
        print('Reading Hectospec mask file: {:s}'.format(mask_file))
        i0 = mask_file.rfind('/')
        mask_nm = mask_file[i0+1:mask_file.find('.cat')]
        # Grab info from spectrum file
        #xdb.set_trace()
        spec_fil = glob.glob(mask_file[:i0+1]+'spHect-'+mask_nm+'.*.fits.gz')
        if len(spec_fil) == 0:
            raise ValueError('Mask not found! {:s}'.format(spec_fil))
            #ras, decs = xra.dtos1((field[1],field[2]))
            #pa=0.
        else:
            header = fits.open(spec_fil[0])[0].header
            if header['APERTURE'] != mask_nm:
                raise ValueError('Mask doesnt match!')
            pa = header['POSANGLE']
            ras = header['CAT-RA']
            decs = header['CAT-DEC']
        # Continuing
        mask_dict = dict(INSTR='HECTOSPEC',MASK_NAME=mask_nm,
            MASK_RA=ras, MASK_DEC=decs, MASK_EPOCH=2000.,
            MASK_PA=pa) # SHOULD GRAB PA, RA, DEC FROM SPECTRA FITS HEADER
        all_masks.append(mask_dict)
        # Read obs
        f = open(mask_file, 'r')
        lines = f.readlines()
        f.close()
        iall_obs = []
        for line in lines:
            if 'OBS' in line:
                prs = line.strip().split(' ')
                gdprs = [iprs for iprs in prs if len(iprs)>0]
                obs_dict = {}
                obs_dict['DATE'] = gdprs[2]
                obs_dict['TEXP'] = float(gdprs[3])
                obs_dict['DISPERSER'] = gdprs[4]
                obs_dict['CONDITIONS'] = gdprs[5]
                # 
                iall_obs.append(obs_dict) 
        obs_tab = xxul.dict_list_to_table(iall_obs)
        obs_tab['TEXP'].unit = u.s
        # Read observed targets
        obs_targ = ascii.read(mask_file,comment='#')
        gdt = np.where(obs_targ['flag'] == 1)[0]
        # Match to target list
        obs_coord = SkyCoord(ra=obs_targ['ra'][gdt]*u.hour, dec=obs_targ['dec'][gdt]*u.deg)
        idx, d2d, d3d = coords.match_coordinates_sky(obs_coord, targ_coord, nthneighbor=1)
        gdm = np.where(d2d < 1.*u.arcsec)[0]
        if len(gdm) != len(gdt):
            raise ValueError('No match')
        else:
            for ii in range(len(gdm)):
                targ_mask['MASK_NAME'][idx[ii]] = mask_nm
                if flg_id == 0:
                    targs['TARG_ID'][idx[ii]] = int(obs_targ['objid'][gdt[ii]])
        """
        for gdi in gdt:
            mtt = np.where(targs['TARG_ID']==
                int(obs_targ['objid'][gdi]))[0]
            if len(mtt) != 1:
                raise ValueError('Multiple matches?!')
            targ_mask['MASK_NAME'][mtt[0]] = mask_nm
        """
        all_obs.append(obs_tab)
    # Add columns to targs
    for tt,cname in enumerate(cnames):
        mask = np.array([False]*len(targs))
        bad = np.where(np.array(targ_mask[cname])==msk_val[tt])[0]
        if len(bad)>0:
            mask[bad]=True
        #
        clm = MaskedColumn(targ_mask[cname],name=cname, mask=mask)
        targs.add_column(clm)

    # Look for ID duplicates (rare)
    gdobj = targs['TARG_ID'] > 0
    idval = np.array(targs[gdobj]['TARG_ID']).astype(int)
    uni, counts = np.unique(idval, return_counts=True)
    if len(uni) != np.sum(gdobj):
        warnings.warn("Found duplicated ID values in Hectospect cat files")
        warnings.warn("Modifying these by hand!")
        dup = np.where(counts>1)[0]
        # Fix by-hand
        for idup in dup:
            dobj = np.where(targs['TARG_ID'] == uni[idup])[0]
            if len(dobj) == 1:
                xdb.set_trace()
            # Confirm RA/DEC are different
            dcoord = SkyCoord(ra=targs['TARG_RA'][dobj]*u.deg,
                              dec=targs['TARG_DEC'][dobj]*u.deg)
            idx, d2d, d3d = coords.match_coordinates_sky(dcoord, dcoord, nthneighbor=2)
            if np.sum(d2d < 1*u.arcsec) > 0:
                raise ValueError("Two with the same RA/DEC.  Deal")
            else:
                for ii in range(1,len(dobj)):
                    # Increment
                    print('Setting TARG_ID to {:d} from {:d}'.format(
                            (ii+1)*targs['TARG_ID'][dobj[ii]],targs['TARG_ID'][dobj[ii]]))
                    targs['TARG_ID'][dobj[ii]] = (ii+1)*targs['TARG_ID'][dobj[ii]]
    # Double check
    idval = np.array(targs[gdobj]['TARG_ID']).astype(int)
    uni, counts = np.unique(idval, return_counts=True)
    if len(uni) != np.sum(gdobj):
        raise ValueError("Cannot happen")

    # Finish
    return all_masks, all_obs, targs


def galaxy_attrib():
    """ List of properties expected to be stored for CASBAH galaxies

    """
    attrib = [ (str('RA'), float),                 # RA (J2000)
               (str('DEC'), float),                # DEC (J2000)
               (str('Z'), float),                  # Redshift
               (str('Z_ERR'), float),              # Redshift uncertainty
               (str('SDSS_MAG'), float, (5,)),     # ugriz photometry from SDSS
               (str('SDSS_MAGERR'), float, (5,)),    # ugriz photometry uncertainties
               (str('TELESCOPE'), '|S80'),            # Telescope(s) used
               (str('INSTRUMENT'), '|S80')            # Instrument(s) used
               ]
    # Return
    return attrib


def grab_sdss_spectra(radec, radius=0.1*u.deg, outfil=None,
                      debug=False, maxsep=None, timeout=600., zmin=None):
    """ Grab SDSS spectra

    Parameters
    ----------
    radec : tuple
      RA, DEC in deg
    radius : float, optional (0.1*u.deg)
      Search radius -- Astroquery actually makes a box, not a circle
    timeout : float, optional
      Timeout limit for connection with SDSS
    outfil : str ('tmp.fits')
      Name of output file for FITS table
    maxsep : float (None) :: Mpc
      Maximum separation to include
    zmin : float (None)
      Minimum redshift to include

    Returns
    -------
    tbl : Table

    """

    cC = coords.SkyCoord(ra=radec[0], dec=radec[1])

    # Query
    photoobj_fs = ['ra', 'dec', 'objid', 'run', 'rerun', 'camcol', 'field']
    mags = ['petroMag_u', 'petroMag_g', 'petroMag_r', 'petroMag_i', 'petroMag_z']
    magsErr = ['petroMagErr_u', 'petroMagErr_g', 'petroMagErr_r', 'petroMagErr_i', 'petroMagErr_z']

    phot_catalog = SDSS.query_region(cC,spectro=True,radius=radius, timeout=timeout,
                                     photoobj_fields=photoobj_fs+mags+magsErr) # Unique
    spec_catalog = SDSS.query_region(cC,spectro=True, radius=radius, timeout=timeout) # Duplicates exist
    nobj = len(phot_catalog)

    #
    print('grab_sdss_spectra: Found {:d} sources in the search box.'.format(nobj))

    # Coordinates
    cgal = SkyCoord(ra=phot_catalog['ra']*u.degree, dec=phot_catalog['dec']*u.degree)
    sgal = SkyCoord(ra=spec_catalog['ra']*u.degree, dec=spec_catalog['dec']*u.degree)
    sepgal = cgal.separation(cC) #in degrees

    # Check for problems and parse z
    zobj = np.zeros(nobj)
    idx, d2d, d3d = coords.match_coordinates_sky(cgal, sgal, nthneighbor=1)
    if np.max(d2d) > 1.*u.arcsec:
        print('No spectral match!')
        xdb.set_trace()
    else:
        zobj = spec_catalog['z'][idx]

    idx, d2d, d3d = coords.match_coordinates_sky(cgal, cgal, nthneighbor=2)
    if np.min(d2d.to('arcsec')) < 1.*u.arcsec:
        print('Two photometric sources with same RA/DEC')
        xdb.set_trace()

    #xdb.set_trace()


    # Cut on Separation
    if not maxsep is None:
        print('grab_sdss_spectra: Restricting to {:g} Mpc separation.'.format(maxsep))
        sepgal_kpc = cosmo.kpc_comoving_per_arcmin(zobj) * sepgal.to('arcmin')
        sepgal_mpc = sepgal_kpc.to('Mpc')
        gdg = np.where( sepgal_mpc < (maxsep * u.Unit('Mpc')))[0]
        phot_catalog = phot_catalog[gdg]
        #xdb.set_trace()

    nobj = len(phot_catalog)
    print('grab_sdss_spectra: Grabbing data for {:d} sources.'.format(nobj))

    # Grab Spectra from SDSS

    # Generate output table
    attribs = galaxy_attrib()
    npix = 5000 #len( spec_hdus[0][1].data.flux )
    spec_attrib = [(str('FLUX'), np.float32, (npix,)),
                   (str('SIG'), np.float32, (npix,)),
                   (str('WAVE'), np.float64, (npix,))]
    tbl = np.recarray( (nobj,), dtype=attribs+spec_attrib)

    tbl['RA'] = phot_catalog['ra']
    tbl['DEC'] = phot_catalog['dec']
    tbl['TELESCOPE'] = str('SDSS 2.5-M')

    # Deal with spectra separately (for now)
    npix = 5000 #len( spec_hdus[0][1].data.flux )

    for idx,obj in enumerate(phot_catalog):
        #print('idx = {:d}'.format(idx))

        # Grab spectra (there may be duplicates)
        mt = np.where( sgal.separation(cgal[idx]).to('arcsec') < 1.*u.Unit('arcsec'))[0]
        if len(mt) > 1:
            # Use BOSS if you have it
            mmt = np.where( spec_catalog[mt]['instrument'] == 'BOSS')[0]
            if len(mmt) > 0:
                mt = mt[mmt[0]]
            else:
                mt = mt[0]
        elif len(mt) == 0:
            xdb.set_trace()
        else:
            mt = mt[0]

        # Grab spectra
        spec_hdus = SDSS.get_spectra(matches=Table(spec_catalog[mt]))

        tbl[idx]['INSTRUMENT'] = spec_catalog[mt]['instrument']
        spec = spec_hdus[0][1].data
        npp = len(spec.flux)
        tbl[idx]['FLUX'][0:npp] = spec.flux
        sig = np.zeros(npp)
        gdi = np.where(spec.ivar > 0.)[0]
        if len(gdi) > 0:
            sig[gdi] = np.sqrt( 1./spec.ivar[gdi] )
        tbl[idx]['SIG'][0:npp] = sig
        tbl[idx]['WAVE'][0:npp] = 10.**spec.loglam

        # Redshifts
        meta = spec_hdus[0][2].data
        for attrib in ['Z','Z_ERR']:
            tbl[idx][attrib] = meta[attrib]

        if debug:
            sep_to_qso = cgal[idx].separation(cC).to('arcmin')
            print('z = {:g}, Separation = {:g}'.format(tbl[idx].Z, sep_to_qso))
            xdb.set_trace()

        # Fill in rest
        tbl[idx].SDSS_MAG = np.array( [obj[phot] for phot in mags])
        tbl[idx].SDSS_MAGERR = np.array( [obj[phot] for phot in magsErr])

    # Clip on redshift to excise stars/quasars
    if zmin is not None:
        gd = np.where(tbl['Z'] > zmin)[0]
        tbl = tbl[gd]

    # Write to FITS file
    if outfil is not None:
        prihdr = fits.Header()
        prihdr['COMMENT'] = 'SDSS Spectra'
        prihdu = fits.PrimaryHDU(header=prihdr)

        tbhdu = fits.BinTableHDU(tbl)

        thdulist = fits.HDUList([prihdu, tbhdu])
        thdulist.writeto(outfil,clobber=True)

    print('Wrote SDSS table to {:s}'.format(outfil))
    return tbl


# ################
if __name__ == "__main__":

    flg_fig = 0
    flg_fig += 2**0  # SDSS search

    # XSpec
    if (flg_fig % 2**1) >= 2**0:
        radec = (212.34957*u.deg,26.30585*u.deg)
        grab_sdss_spectra(radec, radius=1.*u.degree/12.)
