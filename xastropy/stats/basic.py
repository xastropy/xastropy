"""
#;+ 
#; NAME:
#; stats.basic
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for basic stat calculations
#;   04-Dec-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
from scipy.interpolate import interp1d

#from xastropy.xutils import xdebug as xdb

# def perc
# def lin_to_log

def lin_to_log(x, sig):
    """ Convert linear value+error to log 

    Parameters:
      x: float
      sig: float 

    Returns:
      logx, sig_logx
        Log value and error in log

    JXP 26 Mar 2015
    """
    logx = np.log10( x ) 
    lgvar = ((1. / (np.log(10.0)*x))**2) * sig**2
    sig_logx = np.sqrt(lgvar)

    return logx, sig_logx

def perc(x, per=0.68):
    """ Calculate the percentile bounds of a distribution, 
    i.e. for per=0.68, the code returns the upper and lower bounds
    that encompass 68percent of the distribution.

    Parameters:
      x: float
        numpy array of values
      per: float (0.68)
          Percentile for the calulation

    Returns:
      xper: array
        Value at lower, value at upper

    JXP 04 Dec 2014
    """
    #
    npt = len(x)

    # Sort
    xsort = np.sort(x)
    perx = (np.arange(npt)+1) / npt

    f = interp1d(perx,xsort)

    frac = (1.-per) / 2.

    # Fill
    xper = np.zeros(2)
    try:
        xper[0] = f( frac )
    except ValueError:
        xper[0] = np.min(x)

    try:
        xper[1] = f( 1.-frac )
    except ValueError:
        xper[1] = np.max(x)

    #xdb.set_trace()

    # Return
    return xper

def poisson_interval(k, cl=0.95, sigma=None): 
    """Uses chisquared info to get the poisson interval. Uses scipy.stats
    (imports in function). 
    Taken from http://stackoverflow.com/questions/14813530/poisson-confidence-interval-with-numpy
    Checked against my own x_poisscl.pro code in XIDL

    Parameters:
    -----------
    cl: float
      Confidence limit
    """
    from scipy.stats import norm, chi2
    if sigma is not None:
        icl = norm.cdf(sigma)
        cl = 1. - 2*(1.-icl)
    #
    alpha = 1. - cl
    a = alpha
    low, high = (chi2.ppf(a/2, 2*k) / 2, chi2.ppf(1-a/2, 2*k + 2) / 2)
    if k == 0: 
        low = 0.0
    return low, high


def binomial_ci(mle, N, alpha=0.05):
    """ One sided confidence interval for a binomial test.
    To find the two sided interval, call with (1-alpha/2) and alpha/2 as arguments

    Parameters
    ----------
    mle : float
      Fraction of successes
    N : int
      Number of trials

    If after N trials we obtain mle as the proportion of those
    trials that resulted in success, find c such that

    P(k/N < mle; theta = c) = alpha

    where k/N is the proportion of successes in the set of trials,
    and theta is the success probability for each trial.
    """
    from scipy.stats import binom
    from scipy.optimize import bisect


    to_minimise = lambda c: binom.cdf(mle*N,N,c)-alpha
    return bisect(to_minimise,0,1)

