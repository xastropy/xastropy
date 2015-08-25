# Module to run tests on continuum codes

## # TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
import pytest
from astropy import units as u

from xastropy.spec import continuum as xconti

from linetools.spectralline import AbsLine

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''


def test_dict():
    # Init 
    cdict = xconti.init_conti_dict(Norm=1.)

    assert isinstance(cdict,dict)
    np.testing.assert_allclose(cdict['Norm'], 1.)

def test_telfer():
    # 
    telfer = xconti.get_telfer_spec(3.)
    # Grab tau
    np.testing.assert_allclose(telfer.flux[100].value, 2.297435281318983)

def test_igm_telfer():
    telfer = xconti.get_telfer_spec(3.,igm=True)
    #telfer2 = xconti.get_telfer_spec(3.)
    np.testing.assert_allclose(telfer.flux[500].value, 1.1857421030955235)
    #from xastropy.xutils import xdebug as xdb
    #xdb.set_trace()
