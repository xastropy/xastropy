"""
#;+ 
#; NAME:
#; utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Plotting utilities
#;   24-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

try:
    import bokeh
except ImportError:
    print('-----------------------------------------------------------')
    print('-----------------------------------------------------------')
    print('WARNING: Not loading bokeh module in xastropy.plotting   \n Install bokeh if you want it')
    print('-----------------------------------------------------------')


from bokeh.io import output_notebook, show
from bokeh.plotting import figure
from bokeh.models import Range1d


def plot_spec_notebook(spec, title=None, xmnx=None, ymnx=None):
    """ Simple spectrum plot in a Notebook

    Parameters
    ----------
    spec : XSpectrum1D
    title : str, optional
    xmnx : list or tuple or ndarray
      xmin, xmax values
    ymnx : list or tuple or ndarray
      ymin, ymax values

    """
    p = figure(plot_width=900, plot_height=500, title=title)
    # Data
    p.line(spec.dispersion.value, spec.flux.value, color='black', line_width=2)
    p.line(spec.dispersion.value, spec.sig, color='red', line_width=0.5)
    # Labels
    p.xaxis.axis_label = "Wavelength"
    p.yaxis.axis_label = "Flux"
    # Axes
    if xmnx is not None:
        p.set(x_range=Range1d(xmnx[0], xmnx[1]))
    if ymnx is not None:
        p.set(y_range=Range1d(ymnx[0], ymnx[1]))
    # Show
    show(p)


def histogram(arr, xlbl, xrng=None, nbins=20, alpha=1.):
    """ Generate a bokeh histogram

    Parameters
    ----------
    arr : ndarray
    xlbl : str
      Label for the x-axis
    xrng : tuple, optional
      xmin, xmax
    nbins : int, optional
      number of bins
    alpha : float, optional
      alpha for
    """
    if xrng is None:
        xrng = (np.min(arr),np.max(arr))
    p = figure(plot_width=600, plot_height=400)
    # Histogram
    hist, edges = np.histogram(arr, range=xrng, density=True, bins=nbins)
    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], fill_color='blue', alpha=alpha)
    # Label
    p.xaxis.axis_label = xlbl
    # Show
    show(p)


def scatter(xarr, yarr, xlbl=None, ylbl=None, pw=600, ph=400):
    """ Generate a bokeh scatter plot

    Parameters
    ----------
    xarr : ndarray
    yarr : ndarray
    xlbl : str, optional
    ylbl : str, optional
    pw : int, optional
      Width in pixels
    ph : int, optional
      Height in pixels
    """
    p = figure(plot_width=pw, plot_height=ph)
    # Model
    p.circle(xarr, yarr, color='black')#, legend='data')
    # Label
    if xlbl is not None:
        p.xaxis.axis_label = xlbl
    if ylbl is not None:
        p.yaxis.axis_label = ylbl
    # Show
    show(p)