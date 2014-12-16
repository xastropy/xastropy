"""
#;+ 
#; NAME:
#; spec_guis
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for Spectroscopy Guis with QT
#;      These call pieces from spec_widgets
#;   12-Dec-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Import libraries
import numpy as np
import os, sys
import matplotlib.pyplot as plt

from PyQt4 import QtGui
from PyQt4 import QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
# Matplotlib Figure object
from matplotlib.figure import Figure

from xastropy.xutils import xdebug as xdb

from xastropy.xguis import spec_widgets as xspw

# x_specplot replacement
class XSpecGui(QtGui.QMainWindow):
    ''' GUI to replace XIDL x_specplot

        12-Dec-2014 by JXP
    '''
    def __init__(self, spec, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        '''
        spec = Spectrum1D
        '''
        # Build a widget combining several others
        self.main_widget = QtGui.QWidget()

        # Status bar
        self.create_status_bar()

        # Grab the pieces and tie together
        self.pltline_widg = xspw.PlotLinesWidget(status=self.statusBar)
        self.spec_widg = xspw.ExamineSpecWidget(spec,status=self.statusBar,
                                                llist=self.pltline_widg.llist)
        self.pltline_widg.spec_widg = self.spec_widg

        self.spec_widg.canvas.mpl_connect('button_press_event', self.on_click)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(self.spec_widg)
        hbox.addWidget(self.pltline_widg)

        self.main_widget.setLayout(hbox)

        # Point MainWindow
        self.setCentralWidget(self.main_widget)

    def create_status_bar(self):
        self.status_text = QtGui.QLabel("XSpec")
        self.statusBar().addWidget(self.status_text, 1)

    def on_click(self,event):
        if event.button == 3: # Set redshift
            if self.pltline_widg.llist['List'] is None:
                return
            self.select_line_widg = xspw.SelectLineWidget(
                self.pltline_widg.llist[self.pltline_widg.llist['List']])
            self.select_line_widg.exec_()
            line = self.select_line_widg.line
            if line.strip() == 'None':
                return
            #
            wrest = float(line.split('::')[1].lstrip())
            z = event.xdata/wrest - 1.
            self.pltline_widg.llist['z'] = z
            self.statusBar().showMessage('z = {:f}'.format(z))

            self.pltline_widg.zbox.setText('{:.5f}'.format(self.pltline_widg.llist['z']))
    
            # Draw
            self.spec_widg.on_draw()



# Run
def run_xspec(spec):
    ''' Expecting a Spectrum1D class
    '''
    reload(xspw)
    app = QtGui.QApplication(sys.argv)
    gui = XSpecGui(spec)
    gui.show()
    app.exec_()

# ################
# TESTING
if __name__ == "__main__":
    import sys
    from xastropy import spec as xspec

    # Read spectrum
    spec_fil = '/u/xavier/Keck/HIRES/RedData/PH957/PH957_f.fits'
    spec = xspec.readwrite.readspec(spec_fil)

    run_xspec(spec)
