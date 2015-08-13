"""
#;+ 
#; NAME:
#; utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for CASBAH utilities
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

from xastropy.obs import radec as xra

from xastropy.xutils import xdebug as xdb

def get_filename(field, ftype):
	'''Generate a CASBAH file given field and type

	Parameters:
	-----------
	field: tuple
	  (Name, ra, dec)
	ftype: str
	  Type of filename desired
	'''
	if ftype == 'MULTI_OBJ':
		path = os.getenv('CASBAH_GALAXIES')
		filename = path+'/'+field[0]+'/'+field[0]+'_MULTI_OBJ.ascii'
	elif ftype == 'TARGETS':
		path = os.getenv('CASBAH_GALAXIES')
		filename = path+'/'+field[0]+'/'+field[0]+'_targets.ascii'
	else:
		raise ValueError('Not ready for this ftype: {:s}'.format(ftype))
	# Return
	return filename

