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
import os, sys
import matplotlib.pyplot as plt
import glob

from PyQt4 import QtGui
from PyQt4 import QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
# Matplotlib Figure object
from matplotlib.figure import Figure

from xastropy.xutils import xdebug as xdb

# ##################################
# GUI for velocity plot
class WarningWidg(QtGui.QDialog):
    ''' GUI to warn user about coming action and solicit response
        24-Dec-2014 by JXP
    '''
    def __init__(self, message, parent=None):
        '''
        message = str
          Message to display
        '''
        super(WarningWidg, self).__init__(parent)

        # Initialize

        # Grab the pieces and tie together
        z_label = QtGui.QLabel('Warning: {:s}'.format(message))

        # Quit
        nbtn = QtGui.QPushButton('No', self)
        nbtn.clicked.connect(self.touch_no)
        ybtn = QtGui.QPushButton('Yes', self) 
        ybtn.clicked.connect(self.touch_yes)

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(z_label)
        vbox.addWidget(nbtn)
        vbox.addWidget(ybtn)
        self.setLayout(vbox)

    def touch_yes(self):
        self.ans = True
        self.done(0)

    def touch_no(self):
        self.ans = False
        self.done(0)

# ################
# TESTING
if __name__ == "__main__":

    flg_fig = 0 
    flg_fig += 2**0  # Warning

    if (flg_fig % 2) == 1:
        app = QtGui.QApplication(sys.argv)
        app.setApplicationName('Warning')
        main = WarningWidg('Will remove all lines. \n  Continue??')
        main.show()
        app.exec_()
        if main.ans:
            print('You answered yes!')
        else:
            print('You answered no!')
