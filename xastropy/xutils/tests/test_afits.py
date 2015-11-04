# Module to run tests on simple fitting routines for arrays

### TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
import pytest

from xastropy.xutils import afits as xafits
from xastropy.xutils import xdebug as xdb

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)

def test_poly_fit():
    # Generate data
    x = np.linspace(0,np.pi,50)
    y = np.sin(x)
    # Fit
    dfit = xafits.func_fit(x, y, 'polynomial', 3)
    x2 = np.linspace(0,np.pi,100)
    y2 = xafits.func_val(x2,dfit)
    np.testing.assert_allclose(y2[50], 0.97854984428713754)
    #xdb.xplot(x2,y2, xtwo=x,ytwo=y,mtwo='o') 

def test_legend_fit():
    # Generate data
    x = np.linspace(0,np.pi,50)
    y = np.sin(x)
    # Fit
    dfit = xafits.func_fit(x, y, 'legendre', 4)
    x2 = np.linspace(0,np.pi,100)
    y2 = xafits.func_val(x2,dfit)
    np.testing.assert_allclose(y2[50], 0.99940823486206976)
    #xdb.xplot(x2,y2, xtwo=x,ytwo=y,mtwo='o') 

def test_iter_fit():
    # Generate data
    x = np.linspace(0,np.pi,100)
    y = np.sin(x)
    # 
    y[50] = 3.
    # Fit
    dfit, mask = xafits.iter_fit(x, y, 'legendre', 4)
    assert np.sum(mask) == 1
    x2 = np.linspace(0,np.pi,100)
    y2 = xafits.func_val(x2,dfit)
    np.testing.assert_allclose(y2[50], 0.99941444872371643)
