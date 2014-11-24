"""
#;+ 
#; NAME:
#; x_guis
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for Plotting GUIs
#;   07-Sep-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""

# Import libraries
import numpy as np
import pdb
import matplotlib.pyplot as plt

#from xastropy.xutils import xdebug as xdb

#def plot_1d_arrays -- Plot a series of arrays (as many as you want!!)
#def plot_hist -- Plot a simple histogram 

#### ###############################
#  Simplest quick plot
#   Plot a series of arrays (as many as you want!!)
def plot_1d_arrays(*args):
    # Error checking
    if len(args) == 0:
        print 'x_guis.simple_splot: No arguments!'
        return
    
    if not isinstance(args[0],np.ndarray):
        print 'x_guis: Input array is not a numpy.ndarray!'
        return

    # Clear
    plt.clf()
    # Plot it right up
    if len(args) == 1:
        plt.plot(args[0].flatten())
    else: 
        for kk in range(1,len(args)):
            plt.plot(args[0].flatten(),args[kk].flatten())

    # Finish
    #mng = plt.get_current_fig_manager()
    #mng.full_screen_toggle() # Mac specific?! Not working anyhow
    plt.show()

    return
    
#### ###############################
#Plot a simple histogram 
def plot_hist(*args,**kwargs):
    """
    Quick histogram plot

    Parameters
    ----------
    binsz : float 
      Width of the bin
    ax : Matplot pyplot
      Useful for sub plots
    noshow : boolean (False) 
      Set keyword to True to not show to screen
    noclear : boolean (False) 
      Set keyword to True to not clear the figure
    """
    # Error checking
    if len(args) == 0:
        print 'x_guis.simple_splot: No arguments!'
        return
    
    if not isinstance(args[0],np.ndarray):
        print 'x_guis: Input array is not a numpy.ndarray!'
        return

    # Bin size
    if not 'binsz' in kwargs:
        kwargs['binsz'] = float(np.std(args[0].flatten())/5.)

    # Clear
    if not 'noclear' in kwargs:
        plt.clf()

    # Ax
    if not 'ax' in kwargs:
        ax = plt
    else:
        ax = kwargs['ax']
    
    # Plot 
    #pdb.set_trace()
    if len(args) == 1:
        arr = args[0].flatten()

        # Find the range
        minv = np.amin(arr)
        maxv = np.amax(arr)
        # Set the boundaries sensibly given binsz
        i0 = int( minv / kwargs['binsz']) - 1
        i1 = int( maxv / kwargs['binsz']) + 1
        rng = tuple( kwargs['binsz']*np.array([i0,i1]) )
        nbin = i1-i0
        # Histogram
        hist, edges = np.histogram(arr, range=rng, bins=nbin)
        ax.bar(edges[:-1], hist, width=kwargs['binsz'])
        # Labels
        if 'xlabel' in kwargs:
            try:
                ax.set_xlabel(kwargs['xlabel'])
            except: 
                ax.xlabel(kwargs['xlabel'])
    else: 
        pdb.set_trace() # Not ready for this yet
        for kk in range(1,len(args)):
            fig.plot(args[0].flatten(),args[kk].flatten())

    # Finish
    if not 'noshow' in kwargs:
        plt.show()

    return
    
