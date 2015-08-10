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
import os, imp, glob
from astropy.io import fits, ascii
from astropy import units as u 
from astropy.coordinates import SkyCoord
#from astropy import constants as const
from xastropy.casbah import galaxy_data as xcgd

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
	mask_dict = {}
	targ_dict = {}
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
					mask_dict['NAME'] = iprs
				elif cnt==1:
					mask_dict['RAS'] = iprs
				elif cnt==2:
					mask_dict['DECS'] = iprs
				elif cnt==3:
					mask_dict['EPOCH'] = float(iprs)
				elif cnt==4:
					mask_dict['PA'] = float(iprs[3:])
				cnt += 1
			flg_mask = 2
			continue
		# Done?
		if 'Non' in line:
			break
		# Dummy line?
		if (len(line.strip())==0) or (line[0]=='#'):
			continue
		# Selected object
		prs = line.strip().split(' ')
		gdprs = [iprs for iprs in prs if len(iprs)>0]
		targ_dict['ID'] = int(gdprs[0])
		targ_dict['RAS'] = gdprs[1]
		targ_dict['DECS'] = gdprs[2]
		targ_dict['EPOCH'] = float(gdprs[3])
		# Add to Table
		xdb.set_trace()

def deimos_targets(path=None):
	'''Generate files related to DEIMOS deimos_targets
	'''
	if path is None:
		path = './Observing/DEIMOS/'
	fields = [('PG1407+265',212.34957*u.deg,26.30585*u.deg)]

	# Mask info
	for field in fields:
		mask_path = path+field[0]+'/Masks/'
		files = glob.glob(mask_path+'*.out')
		for msk_file in files:
			# Parse
			obj_tab = parse_deimos_mask_file(msk_file)
