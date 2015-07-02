"""
#;+ 
#; NAME:
#; stats.basic
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for likelihood stat calculations
#;   01-Jul-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
from scipy.interpolate import interp1d

from xastropy.xutils import xdebug as xdb

# def perc
# def lin_to_log

def cl_image(lnL, sigma=False):
    """ Calculate a confidence level image from a lnL image
    Simple area under the curve with cubic spline interpolation

    Parameters:
      lnL: np.array
        log-Likelihood image
      sigma: bool, optional
        Return as sigma values [not implemented]

    Returns:
      cl_img: np.array
        Image with the same dimensions with confidence levels
    """
    # Max
    mxL = np.max(lnL)

    # Noramlize and flatten
    norm_img = lnL-mxL
    flat_img = norm_img.flatten()

    # Sort
    srt = np.argsort(flat_img)
    norm_lnL = flat_img[srt]

    # Sum
    cumul_area = np.cumsum(np.exp(np.maximum(norm_lnL,-15.)))
    tot_area = np.max(cumul_area)
    cumul_area = cumul_area/tot_area

    # Interpolation (smoothing a bit)
    f_area = interp1d(norm_lnL, cumul_area)

    # Finish
    area_img = f_area(norm_img)
    cl_img = 1.-area_img

    # Return
    return cl_img
