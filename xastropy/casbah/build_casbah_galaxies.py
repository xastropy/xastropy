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
import os, imp, glob, copy
from astropy.io import fits, ascii
from astropy import units as u 
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column, MaskedColumn, vstack

#from astropy import constants as const
from xastropy.casbah import galaxy_data as xcgd
from xastropy.casbah import utils as xcasbahu
from xastropy.xutils import lists as xxul
from xastropy.obs import radec as xra

from xastropy.xutils import xdebug as xdb

# SDSS
def build_sdss(radius=2.0*u.deg):    
    ''' Grab SDSS photometry and spectra for those fields in the footprint
    Includes BOSS data.
    '''
    fields = [('PG1407+265',212.34957*u.deg,26.30585*u.deg)]

    for field in fields:
        # Directory
        if not os.path.exists(field[0]):
            os.makedirs(field[0])
        # Grab SDSS data + write to folder
        sdss_fil = xcasbahu.get_filename(field,'SDSS')
        print('CASBAH_SDSS: Building {:s}'.format(sdss_fil))
        print('CASBAH_SDSS: Be patient..')
        xcgd.grab_sdss_spectra( (field[1],field[2]), 
            radius=radius, outfil=sdss_fil, maxsep=20., zmin=500./3e5)
        #outfig = os.environ.get('DROPBOX_DIR')+'/CASBAH/Galaxies/SDSS/PG1407+265_SDSS.pdf'

def build_spectra(field,path='./'):
    '''Top-level program to build spectra files

    Parameters:
    -----------
    field: tuple
      (Name, ra, dec)
    '''
    # MMT

def build_imaging(field,path='./'):
    '''Top-level program to build images

    Parameters:
    -----------
    field: tuple
      (Name, ra, dec)
    '''
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
        img_fil = glob.glob(field[0]+'/IMG/LBT/'+msk_img+'*')
        if len(img_fil) == 1:
            # Copy
            path = os.getenv('CASBAH_GALAXIES')
            shutil.copy2(img_fil[0], path+'/'+field[0]+'/')
            print('Copied {:s}'.format(img_fil[0]))
        else:
            raise ValueError('Need to provide the image! {:s}'.format(
                field[0]+'/IMG/LBT/'+msk_img))

def build_targets(field,path='./'):
    '''Top-level program to build target info 

    Parameters:
    -----------
    field: tuple
      (Name, ra, dec)
    '''
    # MMT
    mmt_masks, mmt_obs, mmt_targs = mmt_targets(field)

    # DEIMOS
    deimos_sex, deimos_masks, deimos_obs, deimos_targs = deimos_targets(field)

    # COLLATE
    all_masks = deimos_masks + mmt_masks
    all_obs = deimos_obs + mmt_obs
    all_sex = vstack([deimos_sex,mmt_targs],join_type='inner')  # Use vstack when needed

    # Generate Target table
    targ_file = xcasbahu.get_filename(field,'TARGETS')
    cut_sex = all_sex[['TARG_RA','TARG_DEC','EPOCH','TARG_ID',
        'TARG_MAG','TARG_IMG','INSTR','MASK_NAME']]
    #cut_sex.write(targ_file,overwrite=True)
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


def deimos_targets(field,path=None):
    '''Generate files related to DEIMOS deimos_targets

    Parameters:
    -----------
    field: tuple
      (Name, ra, dec)

    Returns:
    ----------
    '''
    if path is None:
        path = '/Galx_Spectra/DEIMOS/'

    # Loop on Fields
    mask_path = field[0]+path+'/Masks/'
    # SExtractor targeting
    targetting_file = glob.glob(mask_path+'*targ.yaml')
    if len(targetting_file) == 1:
        sex_targ = parse_sex_file(field,targetting_file[0])
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
        xdb.set_trace()
        sex_targ = None
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
    for tt,cname in enumerate(cnames):
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

def mmt_targets(field,path=None):
    '''Read files related to MMT targets

    Parameters:
    -----------
    field: tuple
      (Name, ra, dec)

    Returns:
    ----------
    Target and observing info 
    '''
    if path is None:
        path = '/Galx_Spectra/MMT/'

    # Targets
    targ_path = field[0]+path

    # Target file
    targ_file = glob.glob(targ_path+'*.targ')
    if len(targ_file) != 1:
        raise ValueError('Wrong number of MMT target files')
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
        rad, decd = xra.stod1( (row['RAS'], row['DECS']) )
        targs[k]['TARG_RA'] = rad.value
        targs[k]['TARG_DEC'] = decd.value
    targs.rename_column('objid','TARG_ID')
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
        print('Reading MMT mask file: {:s}'.format(mask_file))
        i0 = mask_file.rfind('/')
        mask_nm = mask_file[i0+1:mask_file.find('.cat')]
        # Grab info from spectrum file
        #xdb.set_trace()
        spec_fil = glob.glob(mask_file[:i0+1]+'spHect-'+mask_nm+'*.fits.gz')
        if len(spec_fil) != 1:
            print('spec_fil -- Not found!'.format(spec_fil))
            ras, decs = xra.dtos1((field[1],field[2]))
            pa=0.
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
        gdt = np.where(obs_targ['flag']==1)[0]
        for gdi in gdt:
            mtt = np.where(targs['TARG_ID']==
                int(obs_targ['objid'][gdi]))[0]
            if len(mtt) != 1:
                raise ValueError('Multiple matches?!')
            targ_mask['MASK_NAME'][mtt[0]] = mask_nm
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
    # Finish
    return all_masks, all_obs, targs
