# Module to run tests on initializing AbsSurvey

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from xastropy.igm.abs_sys.lls_utils import LLSSurvey


'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_dat_list():
    '''JXP format :: Likely to be Deprecated
    '''
    # LLS Survey
    if os.getenv('LLSTREE') is None:
        assert True
    # Load
    lls = LLSSurvey(flist='Lists/lls_metals.lst', tree=os.getenv('LLSTREE'))
    # tests
    np.testing.assert_allclose(lls.NHI[0], 19.25)
    assert lls.nsys == 165