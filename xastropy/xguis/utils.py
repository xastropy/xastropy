"""
#;+ 
#; NAME:
#; utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for Simple Guis with QT
#;   27-Dec-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Import libraries
import numpy as np
import os, sys, copy
import matplotlib.pyplot as plt
import glob

from PyQt4 import QtGui
from PyQt4 import QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
# Matplotlib Figure object
from matplotlib.figure import Figure

from astropy import units as u

from xastropy.xutils import xdebug as xdb





# ################
# TESTING
if __name__ == "__main__":

    flg_fig = 0 
    #flg_fig += 2**0  # Warning
    flg_fig += 2**1  # AnsBox

    if (flg_fig % 2**1) == 2**0:
        app = QtGui.QApplication(sys.argv)
        app.setApplicationName('Warning')
        main = WarningWidg('Will remove all lines. \n  Continue??')
        main.show()
        app.exec_()
        if main.ans:
            print('You answered yes!')
        else:
            print('You answered no!')

    if (flg_fig % 2**2) >= 2**1:
        app = QtGui.QApplication(sys.argv)
        app.setApplicationName('Ans')
        main = AnsBox('Enter redshift',float)
        main.show()
        app.exec_()
        print(main.value)
