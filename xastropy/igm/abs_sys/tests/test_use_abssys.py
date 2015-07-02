# Module to run tests on initializing AbslineSystem

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

def test_grab_line():
    # Init 
    gensys = xabsys.GenericAbsSystem()
    #
    few_lines = [1215.6700, 1334.5323]*u.AA
    z=1.34
    for ilin in few_lines:
        gensys.lines.append(AbsLine(ilin,z=z))

    # Test grab line
    Lya = gensys.grab_line((z,1215.670*u.AA))
    assert Lya.trans == 'HI 1215'
