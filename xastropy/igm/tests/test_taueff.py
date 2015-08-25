# Module to run tests on initializing AbslineSystem

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
import pytest
from astropy import units as u

from xastropy.igm import tau_eff as xit
from xastropy.igm.fN import model as xifm

from linetools.spectralline import AbsLine

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_lya():
    # f(N)
    fN_model = xifm.default_model()

    # tau_eff
    lamb = 1215.6701*(1+2.4)
    teff = xit.ew_teff_lyman(lamb, 2.5, fN_model, NHI_MIN=12., NHI_MAX=17.)
    # Test
    np.testing.assert_allclose(teff, 0.19821448846)

'''
def test_lyx():
    # f(N)
    fN_model = xifm.default_model()

    # tau_eff
    lamb = 917.*(1+2.4)
    teff = xit.ew_teff_lyman(lamb, 2.5, fN_model, NHI_MIN=12., NHI_MAX=17.)
    # Test
    pdb.set_trace()
    np.testing.assert_allclose(teff, 0.19821448846)
'''

def test_parallel():
    import multiprocessing
    # f(N)
    fN_model = xifm.default_model()
    # Generate dicts
    tst_wv = xit.tau_eff_llist()
    adict = []
    for wrest in tst_wv:
        tdict = dict(ilambda=wrest.value*(1+2.4), zem=2.5, fN_model=fN_model)
        adict.append(tdict)

    pool = multiprocessing.Pool(2) # initialize thread pool N threads
    ateff = pool.map(xit.map_etl, adict)
    #pdb.set_trace()
    np.testing.assert_allclose(ateff[-3:], [0.23437970917295792, 0.20262106481068218, 0.21927058621246398])
