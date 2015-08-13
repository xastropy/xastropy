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
from astropy.table import Table, Column, MaskedColumn

#from astropy import constants as const
from xastropy.casbah import galaxy_data as xcgd
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
		outfil = field[0]+'/'+field[0]+'_SDSS.fits'
		print('CASBAH_SDSS: Building {:s}'.format(outfil))
		print('CASBAH_SDSS: Be patient..')
		xcgd.grab_sdss_spectra( (field[1],field[2]), 
			radius=radius, outfil=outfil, maxsep=20.)
	    #outfig = os.environ.get('DROPBOX_DIR')+'/CASBAH/Galaxies/SDSS/PG1407+265_SDSS.pdf'

def build_targets(field,path='./'):
	'''Top-level program to build target info 

	Parameters:
	-----------
	field: tuple
	  (Name, ra, dec)
	'''
	# DEIMOS
	deimos_sex, deimos_masks, deimos_obs, deimos_targs = deimos_targets(field)

	# MMT
	mmt_masks = []
	mmt_obs = []

	# COLLATE
	all_masks = deimos_masks + mmt_masks
	all_obs = deimos_obs + mmt_obs
	all_sex = deimos_sex  # Use vstack when needed

	# Generate Target table
	targ_file = path+field[0]+'/'+field[0]+'_targets.fits'
	xdb.set_trace()
	all_sex.write(targ_file,overwrite=True)
	print('Wrote file {:s}'.format(targ_file))

	# Generate MULTI_OBJ file
	multi_file = path+field[0]+'/'+field[0]+'_MULTI_OBJ.ascii'
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
		# Setup for mask matching
		sex_coord = SkyCoord(ra=sex_targ['TARG_RA']*u.deg, 
			dec=sex_targ['TARG_DEC']*u.deg)
		sex_msk_clms = {}
		cnames = ['INSTR', 'MASK_ID']
		smsk = '--'
		msk_val = [smsk, smsk]
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
					if sex_msk_clms['INSTR'][isep] == smsk:
						sex_msk_clms['INSTR'][isep] = mask_dict['INSTR']
						sex_msk_clms['MASK_ID'][isep] = targ['ID']
					else: # Already full 
						sex_targ.add_row(sex_tarb[isep])
						sex_msk_clms['INSTR'].append(mask_dict['INSTR'])
						sex_msk_clms['MASK_ID'].append(targ['ID'])
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

