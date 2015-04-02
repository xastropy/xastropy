# import code_exercise as ce
#  x=random.normal(scale=1.0, size=1000) + 1.
#  ce.myplot(x, 0.3)
##
#  Import libraries
from __future__ import print_function, absolute_import, division, unicode_literals

assert False # Code is now back in PH136 SVN as code_exercise.py

import numpy as np
import glob, os, sys

from xastropy.xutils import xdebug as xdb

import matplotlib as mpl
mpl.rcParams['font.family'] = 'stixgeneral'
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

# Here's the code to do the plotting
def myplot(arr,ibinsz=0):
    if ibinsz == 0: 
        ibinsz = std(arr)/5.
    binsz = float(ibinsz)  # insure it is a float
    # Find the range
    minv = amin(arr)
    maxv = amax(arr)
    # Set the boundaries sensibly given binsz
    i0 = int( minv / binsz) - 1
    i1 = int( maxv / binsz) + 1
    rng = tuple( binsz*array([i0,i1]) )
    #print 'range', rng
    nbin = i1-i0
    # Histogram
    hist, edges = histogram(arr, range=rng, bins=nbin)

    # plot the histogram
    clf()
    bar(edges[:-1], hist, width=binsz)

    # Time for stats
    rms = std(arr)
    avg = mean(arr)

    # Generate the Gaussian
    xval = arange(float(rng[0])-2.*rms, float(rng[1])+2*rms, binsz/5.)
    yval = amax(hist) * exp(-1. * (xval-avg)**2 / (2 * rms**2) )
    # plot
    plot(xval, yval, 'r', lw=3)

# Here's the code to generate the random numbers
def myran(n,mean,rms):
    ran = random.normal(scale=rms, size=n) + mean
    return(ran)
