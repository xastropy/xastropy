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
import pdb

#def set_fontsize

def set_fontsize(ax,fsz):
    '''
    Generate a Table of columns and so on
    Restrict to those systems where flg_clm > 0

    Parameters
    ----------
    ax : Matplotlib ax class
    fsz : float
      Font size
    '''
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fsz)

#
def whisker_box(ax,xval,yval,per=0.5,color='gray',alpha=0.2):
    '''
    Plot a simple 2D 'whisker' box for a set of points

    Parameters
    ----------
    ax : Matplotlib ax class
    xval: array
      x values
    yval: array
      y values
    '''
    from xastropy import stats as xstat

    xper = xstat.basic.perc(xval, per=per)
    yper = xstat.basic.perc(yval, per=per)
    #pdb.set_trace()
    ax.fill_between(xper, yper[0],
                    np.array((yper[1],yper[1])), color=color, alpha=alpha)

#
def plt_arrows(ax,xval,yval, color='black', up=False, csz=1.5, alen=1.):
    '''
    Plot arrows (upper or lower) on a given axis

    Parameters
    ----------
    ax : Matplotlib ax class
    xval : float or array
    yval : float or array
    up: Bool (False)
      Arrow up or down?
    alen: float (1.5)
      Length of arrow
    csz: float (1.5)
      Cap size
    '''

    ax.errorbar(xval, yval, yerr=alen, ecolor=color, lolims=(up is True),
                uplims=(up is False), fmt='none', capsize=csz)
