"""
#;+ 
#; NAME:
#; afits
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for simple fitting routines
#;   27-Oct-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import copy
import warnings

from scipy.interpolate import splev, splrep # bspline

from xastropy.xutils import xdebug as xdb

# def bintab_to_table(fits_fil,exten=1, silent=False):
# def table_to_fits(table, outfil, compress=False, comment=None):

def bspline_inner_knots(all_knots):
    '''Trim to the inner knots.  Used in bspline_magfit
    Might be useful elsewhere
    Assumes the outer knots are identical

    Parameters:
    ---------
    all_knots: ndarray
      Array of knots returned by scipy.  Includes outer knots

    Reults:
    ---------
    inner_knots: ndarray
      Trimmed down array.  
    '''
    diff = all_knots - np.roll(all_knots,1)
    pos = np.where(diff>0.)[0]
    i0=pos[0]
    i1=pos[-1]
    return all_knots[i0:i1]

def bspline_fit(x,y,order=3,w=None, knots=None,everyn=None,bkspace=None,
    xmin=None,xmax=None):
    ''' bspline fit to x,y
    Should only be called from func_fit

    Parameters:
    ---------
    x: ndarray
    y: ndarray
    func: str
      Name of the fitting function:  polynomial, legendre, chebyshev, bspline
    deg: int 
      deg of the spline.  Default=3 (cubic)
    xmin: float, optional
      Minimum value in the array  [both must be set to normalize]
    xmax: float, optional
      Maximum value in the array  [both must be set to normalize]
    w: ndarray, optional
      weights to be used in the fitting (weights = 1/sigma)
    everyn: int 
      Knot everyn good pixels, if used
    bkspace: float 
      Spacing of breakpoints in units of x

    Returns:
    ---------
    fit_dict: dict  
      dict describing the bspline fit 
    ''' 
    # Save args for later
    args = locals().copy()
    # Generate dict
    fit_dict = dict(func='bspline')
    for key in args.keys():
        if key in ['x','y','w']:
            continue
        else:
            fit_dict[key] = args[key]
    # Normalize?
    if xmin is not None and xmin is not None:
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
    else:
        xv = x
    #
    if w is None:
        ngd = xv.size
        gd = np.arange(ngd)
        weights = None
    else:
        gd = np.where(w > 0.)[0]
        weights = w[gd]
    # Make the knots
    if knots is None:
        if bkspace is not None: 
            xrnge = (np.max(x) - np.min(x))
            startx = np.min(x)
            nbkpts = max(int(xrnge/bkspace) + 1,2)
            tempbkspace = xrnge/(nbkpts-1)
            knots = np.arange(1,nbkpts-1)*tempbkspace + startx
        elif everyn is not None:
            idx_knots = np.arange(10, ngd-10, everyn) # A knot every good N pixels
            knots = xv[gd[idx_knots]]
        else:
            raise IOError("No method specified to generate knots")
    # Generate spline
    tck = splrep( xv[gd], y[gd], w=weights, k=order, t=knots)
    # Update dict
    fit_dict['tck'] = tck
    fit_dict['knots'] = knots

    return fit_dict

#
def func_fit(x,y,func,deg,xmin=None,xmax=None,w=None, **kwargs):
    ''' Simple function fit to 2 arrays
    Modified code originally from Ryan Cooke (PYPIT)

    Parameters:
    ---------
    x: ndarray
    y: ndarray
    func: str
      Name of the fitting function:  polynomial, legendre, chebyshev, bspline
    deg: int or dict
      Order of the fit or a dict for bspline fits
    xmin: float, optional
      Minimum value in the array (or the left limit for a legendre/chebyshev polynomial)
    xmax: float, optional
      Maximum value in the array (or the left limit for a legendre/chebyshev polynomial)
    w: ndarray, optional
      weights to be used in the fitting (weights = 1/sigma)

    Returns:
    ---------
    fit_dict: dict  
      dict describing the Fit including the coefficients
    ''' 
    # Normalize
    if xmin is None or xmax is None:
        if np.size(x) == 1:
            xmin, xmax = -1.0, 1.0
        else:
            xmin, xmax = np.min(x), np.max(x)    
    xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
    # Fit
    if func == "polynomial":
        fit = np.polynomial.polynomial.polyfit(xv,y,deg,w=w)
    elif func == "legendre":
        fit = np.polynomial.legendre.legfit(xv,y,deg,w=w)
    elif func == "chebyshev":
        fit = np.polynomial.chebyshev.chebfit(xv,y,deg,w=w)
    elif func == "bspline":
        fit_dict = bspline_fit(xv,y,order=deg,w=w,**kwargs)
        fit_dict['xmin'] = xmin
        fit_dict['xmax'] = xmax
        return fit_dict
    else:
        raise ValueError("Fitting function '{0:s}' is not implemented yet".format(func))
    # Finish    
    fit_dict = dict(coeff=fit, order=deg, func=func, xmin=xmin, xmax=xmax, **kwargs)
    return fit_dict

def func_val(x,fit_dict):
    ''' Get values from a fit_dict
    Modified code originally from Ryan Cooke (PYPIT)

    Parameters:
    ---------
    x: ndarray

    Returns:
    ---------
    fit_dict: dict  
      dict describing the Fit including the coefficients
    '''
    xv = 2.0 * (x-fit_dict['xmin'])/(fit_dict['xmax']-fit_dict['xmin']) - 1.0
    if fit_dict['func'] == "polynomial":
        return np.polynomial.polynomial.polyval(xv,fit_dict['coeff'])
    elif fit_dict['func'] == "legendre":
        return np.polynomial.legendre.legval(xv,fit_dict['coeff'])
    elif fit_dict['func'] == "chebyshev":
        return np.polynomial.chebyshev.chebval(xv,fit_dict['coeff'])
    elif fit_dict['func'] == "bspline":
        return splev(xv, fit_dict['tck'],ext=1)
    else:
        raise ValueError("Fitting function '{0:s}' is not implemented yet".format(fit_dict['func']))

def iter_fit(xarray, yarray, func, order, weights=None, sigma=None, max_rej=None,  
    maxone=True, sig_rej=3.0, initialmask=None, forceimask=False, 
    xmin=None, xmax=None, niter=999, debug=False, **kwargs):
    """A "robust" fit with iterative rejection is performed to the xarray, yarray pairs
    Modified code originally from Ryan Cooke (PYPIT)

    Parameters:
    ----------
    xarray: ndarray
      independent variable values
    yarray: ndarray
      dependent variable values
    func: str
      Name of the fitting function:  polynomial, legendre, chebyshev
    order: int
      the order of the function to be used in the fitting
    sigma: ndarray, optional
      Error in the yvalues.  Used only for rejection
    weights: ndarray, optional
      weights to be used in the fitting (weights = 1/sigma)
    maxone: bool, optional [True]
      If True, only the most deviant point in a given iteration will be removed
    sig_rej: float, optional [3.0]
      confidence interval for rejection 
    max_rej: int, optional [None]
      Maximum number of points to reject
    initialmask: ndarray
       a mask can be supplied as input, these values will be masked for the first iteration. 1 = value masked
    forceimask: bool, optional [False]
      If True, the initialmask will be forced for all iterations
    niter: int, optional [999]
      Maximum number of iterations 
    xmin: float
      minimum value in the array (or the left limit for a legendre/chebyshev polynomial)
    xmax: float
      maximum value in the array (or the right limit for a legendre/chebyshev polynomial)
    debug: bool, optional

    Returns:
    -------
    return: mask, fit_dict
       mask is an array of the masked values, 
       fit_dict is a dict containing the fit 
    """
    # Setup the initial mask
    if initialmask is None:
        mask = np.zeros(xarray.size,dtype=np.int)
        if forceimask:
            warnings.warn("Initial mask cannot be enforced -- no initital mask supplied")
            forceimask = False
    else:
        mask = initialmask.copy()
    # Avoid zero or negative weights
    if weights is not None:
        mask[weights <= 0.] = 1
    mskcnt=np.sum(mask)
    imskcnt=copy.copy(mskcnt)
    # Iterate, and mask out new values on each iteration
    iiter = 0
    while True:
        iiter += 1
        if iiter > niter:
            warnings.warn("Reached maximum number of iterations")
            break
        # Mask
        w = np.where(mask==0)
        xfit = xarray[w]
        yfit = yarray[w]
        # Fit
        dfit = func_fit(xfit,yfit,func,order,xmin=xmin,xmax=xmax, **kwargs)
        yrng = func_val(xarray, dfit) 
        # Reject
        sigmed = 1.4826*np.median(np.abs(yfit-yrng[w]))
        if debug:
            import xpdb
            xpdb.set_trace()
        # Check number of parameters
        if xarray.size-np.sum(mask) <= order+2:
            warnings.warn("More parameters than data points - fit might be undesirable")
            break # More data was masked than allowed by order
        if maxone: # Only remove the most deviant point
            if sigma is not None:
                tst = np.abs(yarray[w]-yrng[w])/sigma[w]
                m = np.argmax(tst)
                if tst[m] > sig_rej:
                    mask[w[0][m]] = 1
            else:
                tst = np.abs(yarray[w]-yrng[w])
                m = np.argmax(tst)
                if tst[m] > sig_rej*sigmed:
                    mask[w[0][m]] = 1
        else:
            if forceimask:
                if sigma is not None:
                    w = np.where((np.abs(yarray-yrng) > sig_rej*sigma) | (initialmask==1))
                else:
                    w = np.where((np.abs(yarray-yrng) > sig_rej*sigmed) | (initialmask==1))
            else:
                if sigma is not None:
                    w = np.where(np.abs(yarray-yrng) > sig_rej*sigma)
                else:
                    w = np.where(np.abs(yarray-yrng) > sig_rej*sigmed)
            mask[w] = 1
        if mskcnt == np.sum(mask): break # No new values have been included in the mask
        if max_rej is not None:
            if mskcnt-imskcnt > max_rej:
                break
        mskcnt = np.sum(mask)
        w = np.where(mask==0)
    # Final fit
    xfit = xarray[w]
    yfit = yarray[w]
    fdict = func_fit(xfit,yfit,func,order,xmin=xmin,xmax=xmax,**kwargs)
    return fdict, mask

def normalize(x,xmin,xmax):
    '''Normalize an array

    Parameters:
    -----------
    x: ndarray
      Array to normalize
    xmin: float
      minimum value in the array (or the left limit for a legendre/chebyshev polynomial)
    xmax: float
      maximum value in the array (or the right limit for a legendre/chebyshev polynomial)

    Returns:
    --------
    xnorm: ndarray
    '''
