# Tests the codes in SDSS Quasars
import numpy as np
import os, pdb
import pytest

from astropy import units as u

from xastropy.sdss import quasars as sdssq

def test_init_class():
	sdss_dr7 = sdssq.SdssQuasars()
	assert sdss_dr7._version == 'DR7'

def test_parse_by_plate_fiber():
    if os.getenv('SDSSPATH') is None:
        assert True
        return
	sdss_dr7 = sdssq.SdssQuasars()
	row = sdss_dr7.get_qso((287,264))
	np.testing.assert_allclose(row['Z'], 0.331188)

def test_parse_by_name():
    if os.getenv('SDSSPATH') is None:
        assert True
        return
	sdss_dr7 = sdssq.SdssQuasars()
	row = sdss_dr7.get_qso('J000009.42-102751.9')
	np.testing.assert_allclose(row['Z'], 1.84493)
