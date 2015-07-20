# Tests the codes in xastropy.obs 
#   Requires "requests", "astropy", "aplpy", "PIL"
import numpy as np
import os, pdb
import pytest

from astropy import units as u

from xastropy.obs import finder as xf
from xastropy.obs import radec as x_r

'''
def test_finder():
    # SDSS
    xf.main(['TST', '10:31:38.87', '+25:59:02.3'],radec=1, imsize=8.)#,BW=1)
    # DSS
    xf.main(['TST2', '10:31:38.87', '-25:59:02.3'],radec=1, DSS=1, imsize=8.,BW=1)
'''
def test_stod():
	radec = x_r.stod1('J103138.87+255902.3')
	np.testing.assert_allclose(radec[0].value, 157.91195833333336)
	assert radec[0].unit == u.deg

def test_tocoord():
	from astropy.coordinates import SkyCoord
	radec = x_r.stod1('J103138.87+255902.3')
	coord = x_r.to_coord(radec)
	# 
	assert isinstance(coord,SkyCoord)
