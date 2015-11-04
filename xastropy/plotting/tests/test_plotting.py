# Module to run tests on simple plotting code

# TEST_UNICODE_LITERALS

import matplotlib
matplotlib.use('Agg')  # For Travis

import numpy as np
import os, pdb
import pytest

from xastropy.plotting import simple 

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)

'''
def test_load_kin():
    # Class
    cos_halos = COSHalos()
    cos_halos.load_mega(skip_ions=True)
    # Load kin
    cos_halos.load_abskin()
'''

'''
def test_simple_array():
    # test
    x = np.arange(100)
    simple.plot_1d_arrays(x,outfil='tst.pdf')
    assert True

def test_simple_xtwo():
    # Class
    x = np.linspace(0,2*np.pi,100)
    y = np.sin(x)
    xtwo = np.linspace(0,2*np.pi,50)
    ytwo = np.cos(xtwo)
    simple.plot_1d_arrays(x,y,xtwo=xtwo,ytwo=ytwo,outfil='tst.pdf')
    # 
    xtwo = [np.argmax(y), np.argmin(y)]
    ytwo = y[xtwo]
    simple.plot_1d_arrays(x,y,xtwo=x[xtwo],ytwo=ytwo,mtwo='o',outfil='tst.pdf')
    assert True
'''


