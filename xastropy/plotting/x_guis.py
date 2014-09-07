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
from matplotlib import pyplot

#### ###############################
#  Simplest quick plot
#   Plot a series of arrays (as many as you want!!)
#
def plot_1d_arrays(*args):
    # Error checking
    if len(args) == 0:
        print 'x_guis.simple_splot: No arguments!'
        return
    
    if not isinstance(args[0],np.ndarray):
        print 'x_guis: Input array is not a numpy.ndarray!'
        return

    # Clear
    pyplot.clf()
    # Plot it right up
    if len(args) == 1:
        pyplot.plot(args[0].flatten())
    else: 
        for kk in range(1,len(args)):
            pyplot.plot(args[0].flatten(),args[kk].flatten())

    # Finish
    #mng = pyplot.get_current_fig_manager()
    #mng.full_screen_toggle() # Mac specific?! Not working anyhow
    pyplot.show()

    return
    
