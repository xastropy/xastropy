# Tests the codes in SDSS Quasars
import numpy as np
import os, pdb
import pytest

from astropy import units as u

from xastropy.sdss import quasars as sdssq
from xastropy.sdss import qso as sdssqso

def test_init_class():
	qso = sdssqso.SdssQso()
	np.testing.assert_allclose(qso.z,0.)

def test_from_database():
    if os.getenv('SDSSPATH') is None:
        assert True
        return
    #
	sdss_dr7 = sdssq.SdssQuasars()
	qso = sdss_dr7[(287,264)]
	np.testing.assert_allclose(qso.z, 0.331188)

