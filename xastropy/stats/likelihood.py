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
import os, copy
from scipy.interpolate import interp1d

from xastropy.xutils import xdebug as xdb

# def cl_image
# def cl_interval

def cl_image(lnL, sigma=False):
    """ Calculate a confidence level image from a lnL image
    Simple area under the curve with cubic spline interpolation

    Parameters:
      lnL: np.array
        log-Likelihood image
        Should probably be 2D
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

def cl_interval(lnL, sigma=None, CL=0.68):
    """ Calculate a confidence level interval from a log-likelihood image
    Simple area under the curve with the image collapsed along each
    dimension

    Parameters:
      lnL: np.array
        log-Likelihood image
      CL: float, optional
      sigma: float, optional
        Use to calculate confindence interval

    Returns:
      best_idx, all_error: Lists
        [best] [-, +] indices for each dimension
    """
    # Confidence limits
    if sigma is None:
        c0 = (1. - CL)/2.
        c1 = 1.-c0
    # Image dimensions
    shape = lnL.shape
    ndim = len(shape)
    slc = [slice(None)]*ndim 
    # Find max
    norm_L = np.exp(np.maximum(lnL - np.max(lnL),-15.))
    # Find best indices 
    indices = np.where(lnL == np.max(lnL))
    best_idx = [bi[0] for bi in indices]

    # Error intervals
    all_error = []
    for kk in range(ndim):
        # Collapse on this dimension
        slc = copy.deepcopy(best_idx)
        slc[kk] = slice(None)
        Lslice = norm_L[slc].flatten()
        # Interpolate and go
        cumul_area = np.cumsum(Lslice)
        f_area = interp1d(cumul_area/cumul_area[-1], np.arange(len(Lslice)))
        # Here we go
        idx0 = int(np.round(f_area(c0)))
        idx1 = int(np.round(f_area(c1)))
        all_error.append([idx0,idx1])

    # Return
    return best_idx, all_error


def cl_indices(lnL, cl, sigma=False):
    """ Find the indices of a log-Likelihood grid encompassing a 
    given confidence interval

    Parameters:
      lnL: np.array
        log-Likelihood image
      sigma: bool, optional
        Return as sigma values [not implemented]

    Returns:
      indices: Tuple of np.where output
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
    cumulsum = np.cumsum(np.exp(np.maximum(norm_lnL,-15.)))
    cumul = cumulsum/cumulsum[-1]

    # Interpolation (smoothing a bit)
    fsum = interp1d(norm_lnL, cumul)

    # Finish
    indices = np.where(fsum(norm_img) > (1-cl))

    # Return
    return indices
