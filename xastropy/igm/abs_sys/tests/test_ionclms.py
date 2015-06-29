# Module to run tests on initializing AbsLine

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
import pytest
from astropy import units as u

from xastropy.igm.abs_sys.ionclms import IonClms

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_init_ionclms_all():
	# Init from all file
    ioncs = IonClms(all_file=data_path('UM184.z2929_MAGE.all'))

    assert len(ioncs._data) == 14

def test_ionclms_items():
	# Init from all file
	ioncs = IonClms(all_file=data_path('UM184.z2929_MAGE.all'))

	np.testing.assert_allclose(ioncs['SiII']['clm'], 13.7)
	np.testing.assert_allclose(ioncs[(14,2)]['clm'], 13.7)

def test_ionclms_attr():
	# Init from all file
	ioncs = IonClms(all_file=data_path('UM184.z2929_MAGE.all'))

	assert ioncs.Z[3] == 6
	