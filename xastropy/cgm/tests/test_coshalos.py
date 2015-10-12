# Module to run tests on initializing AbslineSystem

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
import pytest
from astropy import units as u

from xastropy.cgm.cos_halos import COSHalos

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_load_sngl():
    # Class
    cos_halos = COSHalos(fits_path='files/')
    # Load
    cos_halos.load_single( ('J0950+4831','177_27'))

def test_load_survey():
    # Class
    cos_halos = COSHalos(fits_path='files/')
    # Load
    cos_halos.load_mega()  # Only reads one file, actually
    cos_halos.load_mega(skip_ions=True)

def test_getitem():
    # Class
    cos_halos = COSHalos(fits_path='files/')
    # Load
    cos_halos.load_single( ('J0950+4831','177_27'))
    # Get item
    cgm_abs = cos_halos[('J0950+4831','177_27')]
    assert cgm_abs.field == 'J0950+4831'

