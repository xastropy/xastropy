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

# def EditBox
# def WriteQuitWidget
# def WarningWidg


class AnsBox(QtGui.QDialog):
    '''Solicit an input answer from the User
    lbl: str
    format: str
      Format for value
    '''
    def __init__(self, lbl, format=str, parent=None):
        '''
        '''
        super(AnsBox, self).__init__(parent)

        self.format=format
        # 
        label = QtGui.QLabel(lbl) 
        self.box = QtGui.QLineEdit()
        self.box.setMinimumWidth(90)
        # Connect
        self.connect(self.box, 
            QtCore.SIGNAL('editingFinished ()'), self.setv)
        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(label)
        vbox.addWidget(self.box)
        self.setLayout(vbox)

    def setv(self):
        try:
            self.value = self.format(unicode(self.box.text()))
        except ValueError:
            print('Bad input value! Try again with right type')
        else:
            self.done(0)

class EditBox(QtGui.QWidget):
    '''
    initv: Initial value
    lbl: str
    format: str
      Format for value
    '''
    def __init__(self, initv, lbl, format, parent=None):
        '''
        '''
        super(EditBox, self).__init__(parent)

        self.value = initv
        # 
        label = QtGui.QLabel(lbl) 
        self.box = QtGui.QLineEdit()
        # Format
        self.box.frmt = format
        self.box.setText(self.box.frmt.format(self.value))
        self.box.setMinimumWidth(90)
        # Connect
        self.connect(self.box, 
            QtCore.SIGNAL('editingFinished ()'), self.setv)
        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(label)
        vbox.addWidget(self.box)
        self.setLayout(vbox)

    def setv(self):
        self.value = unicode(self.box.text())

    def set_text(self,value):
        self.value = value
        self.box.setText(self.box.frmt.format(self.value))

# ##################################
class WriteQuitWidget(QtGui.QWidget):
    def __init__(self, parent=None):
        '''
        '''
        super(WriteQuitWidget, self).__init__(parent)
        self.parent = parent

        # Generate Buttons
        wbtn = QtGui.QPushButton('Write', self)
        wbtn.setAutoDefault(False)
        wbtn.clicked.connect(self.parent.write_out)

        wqbtn = QtGui.QPushButton('Write\n Quit', self)
        wqbtn.setAutoDefault(False)
        wqbtn.clicked.connect(self.parent.write_quit)

        qbtn = QtGui.QPushButton('Quit', self)
        qbtn.setAutoDefault(False)
        qbtn.clicked.connect(self.parent.quit)

        # Layout
        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(wbtn)
        hbox.addWidget(wqbtn)
        hbox.addWidget(qbtn)
        self.setLayout(hbox)


# ##################################
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
