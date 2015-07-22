# Tests the codes in SDSS Quasars
import numpy as np
import os, pdb
import pytest

from astropy import units as u

from xastropy.sdss import quasars as sdssq

def test_parse_by_plate_fiber():
	#
	row = sdssq.get_qso((287,264))
	np.testing.assert_allclose(row['Z'], 0.331188)

def test_parse_by_name():
	row = sdssq.get_qso('J000009.42-102751.9')
	np.testing.assert_allclose(row['Z'], 1.84493)
