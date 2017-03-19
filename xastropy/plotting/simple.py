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
# def add_subplot_axes -- Add a subplot to a given axis
#### ###############################
#  Simplest quick plot
#   Plot a series of arrays (as many as you want!!)
def plot_1d_arrays(*args,**kwargs):
    """
    Plot arrays

    Parameters
    ----------
    outfil= : string
      Outfil
    xlbl,ylbl= : string
      Labels for x,y axes
    xrng= list
      Range of x limits
    yrng= list
      Range of y limits
    xtwo= : ndarray
      x-values for a second array
    ytwo= : ndarray 
      y-values for a second array
    mtwo= : str
      marker for xtwo
    scatter= : Bool
      True for a scatter plot
    NOTE: Any extra parameters are fed as kwargs to plt.plot()
    """
    # Error checking
    if len(args) == 0:
        print('x_guis.simple_splot: No arguments!')
        return
    
    if not isinstance(args[0],np.ndarray):
        print('x_guis: Input array is not a numpy.ndarray!')
        return

    plt_dict = {}

    # Outfil
    if ('outfil' in kwargs):
        plt_dict['outfil'] = kwargs['outfil']
        kwargs.pop('outfil')
    else:
        plt_dict['outfil'] = None

    # Scatter plot?
    if ('scatter' in kwargs):
        #kwargs['marker'] = 'o'
        kwargs.pop('scatter')
        plt_dict['flg_scatt'] = 1
    else:
        plt_dict['flg_scatt'] = 0

    # Second array?
    if ('xtwo' in kwargs) & ('ytwo' in kwargs):
        plt_dict['xtwo'] = kwargs['xtwo']
        kwargs.pop('xtwo')
        plt_dict['ytwo'] = kwargs['ytwo']
        kwargs.pop('ytwo')
        plt_dict['flg_two'] = 1
        # mtwo
        if 'mtwo' in kwargs:
            plt_dict['mtwo']=kwargs['mtwo']
            kwargs.pop('mtwo')
        else:
            plt_dict['mtwo']=''
    else:
        plt_dict['flg_two'] = 0

    # Limits
    for irng in ['xrng','yrng']:
        try:
            plt_dict[irng] = kwargs[irng]
        except KeyError:
            plt_dict[irng] = None
        else:
            kwargs.pop(irng)

    # Labels
    for ilbl in ['xlbl','ylbl']:
        try:
            plt_dict[ilbl] = kwargs[ilbl]
        except KeyError:
            plt_dict[ilbl] = None
        else:
            kwargs.pop(ilbl)

    # Clear
    plt.clf()
    # Plot it right up
    if len(args) == 1:
        plt.plot(args[0].flatten(), **kwargs)
    else: 
        for kk in range(1,len(args)):
            if plt_dict['flg_scatt'] == 0:
                plt.plot(args[0].flatten(),args[kk].flatten(), **kwargs)
            else:
                plt.scatter(args[0].flatten(),args[kk].flatten(), marker='o', **kwargs)

    if plt_dict['flg_two'] == 1:
        plt.plot(plt_dict['xtwo'], plt_dict['ytwo'], plt_dict['mtwo'], color='red', **kwargs)

    # Limits
    if plt_dict['xrng'] is not None:
        plt.xlim(plt_dict['xrng'])
    if plt_dict['yrng'] is not None:
        plt.ylim(plt_dict['yrng'])

    # Label
    if plt_dict['xlbl'] is not None:
        plt.xlabel(plt_dict['xlbl'])
    if plt_dict['ylbl'] is not None:
        plt.ylabel(plt_dict['ylbl'])

    # Output?
    if plt_dict['outfil'] is not None:
        plt.savefig(plt_dict['outfil']) 
        print('Wrote figure to {:s}'.format(plt_dict['outfil']))
    else: # Show
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
    xmnx : tuple, optional
      (xmin, xmax) for plotting
    """
    # Error checking
    if len(args) == 0:
        print('x_guis.simple_splot: No arguments!')
        return
    
    if not isinstance(args[0],np.ndarray):
        print('x_guis: Input array is not a numpy.ndarray!')
        return

    plt_dict = {}

    # Bin size
    if not 'alpha' in kwargs:
        kwargs['alpha'] = 1.

    # Bin size
    if not 'binsz' in kwargs:
        if 'xrng' in kwargs:
            tmp = args[0].flatten()
            idx = np.where((tmp > kwargs['xrng'][0]) & (tmp < kwargs['xrng'][1]))[0]
            kwargs['binsz'] = float(np.std(tmp[idx])/5.)
        else:
            kwargs['binsz'] = float(np.std(args[0].flatten())/5.)
    #pdb.set_trace()

    # Clear
    if (not 'noclear' in kwargs) and (not 'ax' in kwargs):
        plt.clf()

    # Ax
    if not 'ax' in kwargs:
        ax = plt.gca()
    else:
        ax = kwargs['ax']
    
    # Plot 
    #pdb.set_trace()
    if len(args) == 1:
        arr = args[0].flatten()

        # Find the range
        if 'xrng' not in kwargs:
            minv = np.amin(arr)
            maxv = np.amax(arr)
        else:
            minv = kwargs['xrng'][0]
            maxv = kwargs['xrng'][1]
        # Set the boundaries sensibly given binsz
        i0 = int( minv / kwargs['binsz']) - 1
        i1 = int( maxv / kwargs['binsz']) + 1
        rng = tuple( kwargs['binsz']*np.array([i0,i1]) )
        nbin = i1-i0
        # Histogram
        hist, edges = np.histogram(arr, range=rng, bins=nbin)
        ax.bar(edges[:-1], hist, width=kwargs['binsz'], alpha=kwargs['alpha'])
        # Labels
        if 'xlabel' in kwargs:
            ax.set_xlabel(kwargs['xlabel'])
        if 'xmnx' in kwargs:
            ax.set_xlim(kwargs['xmnx'])
    else: 
        pdb.set_trace() # Not ready for this yet
        for kk in range(1,len(args)):
            fig.plot(args[0].flatten(),args[kk].flatten())

    # Finish
    if (not 'noshow' in kwargs) and (not 'ax' in kwargs):
        plt.show()

    return
    
#### ###############################
#Inset (stolen from http://stackoverflow.com/questions/
    #    17458580/embedding-small-plots-inside-subplots-in-matplotlib )
def add_subplot_axes(ax, rect, axisbg='w'):
    """
    Add an inset to a subplot.
    This can run slowly..

    Parameters
    ----------
    ax: Axis
      matplotlib subplot
    rect: list of 4 floats
      x, y, Dx, Dy of the inset
    """
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    #
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    #
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    # Fonts
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    #
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize) 
    # Return
    return subax

#### ###############################
#Inset (stolen from http://stackoverflow.com/questions/
    #    17458580/embedding-small-plots-inside-subplots-in-matplotlib )
def dark_bkgd(matplt):
    '''Set up matplotlib for a dark background'''
    matplt.rcParams['lines.color']= 'white'
    matplt.rcParams['patch.edgecolor']= 'white'

    matplt.rcParams['text.color']= 'white'

    matplt.rcParams['axes.facecolor']= 'black'
    matplt.rcParams['axes.edgecolor']= 'white'
    matplt.rcParams['axes.labelcolor']= 'white'

    matplt.rcParams['xtick.color']= 'white'
    matplt.rcParams['ytick.color']= 'white'

    matplt.rcParams['grid.color']= 'white'

    matplt.rcParams['figure.facecolor']= 'black'
    matplt.rcParams['figure.edgecolor']= 'black'

    matplt.rcParams['savefig.facecolor']= 'black'
    matplt.rcParams['savefig.edgecolor']= 'black'
    return

## #################################    
## #################################    
## TESTING
## #################################    
if __name__ == '__main__':

    flg_test = 0
    flg_test = 1  # simple plots
    #flg_test += 2 # LLS plot
    #flg_test += 2**2 # zpeak
    #
    #flg_test += 2**9 # LLS Survey NHI
    #flg_test += 2**10 # LLS Survey ions

    # simple plots
    if (flg_test % 2**1) >= 2**0:
        xval = np.arange(10)
        yval = xval**2
        #plot_1d_arrays(xval, yval) # line
        plot_1d_arrays(xval, yval, scatter='True') #marker='o', linestyle='None') # dots
