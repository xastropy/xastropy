# Module to run tests on initializing AbsLine

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
import pytest
from astropy import units as u

from xastropy.igm.abs_sys import abssys_utils as xabsys

from linetools.spectralline import AbsLine

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_generic():
	# Init 
    gensys = xabsys.Generic_System(NHI=16., zabs=1.244)
    #
    assert gensys.abs_type == 'Generic'

def test_lines_generic():
    # Init 
    gensys = xabsys.Generic_System()
    #
    few_lines = [1215.6700, 1334.5323]*u.AA
    for ilin in few_lines:
        gensys.lines.append(AbsLine(ilin))

    # Test
    Lya = gensys[1215.670*u.AA]
    assert Lya[0].trans == 'HI 1215'
