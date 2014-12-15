"""
#;+ 
#; NAME:
#; spec_widgets
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for Spectroscopy widgets with QT
#;   12-Dec-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Import libraries
import numpy as np
import os, sys, imp
import matplotlib.pyplot as plt

from PyQt4 import QtGui
from PyQt4 import QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
# Matplotlib Figure object
from matplotlib.figure import Figure

from xastropy import spec as xspec 
from xastropy.xutils import xdebug as xdb

xa_path = imp.find_module('xastropy')[1]

class ExamineSpecWidget(QtGui.QWidget):
    ''' Widget to plot a spectrum and interactively
        fiddle about.  Akin to XIDL/x_specplot.pro

        12-Dec-2014 by JXP
    '''
    def __init__(self, spec, parent=None, status=None, llist=None):
        '''
        spec = Spectrum1D
        '''
        super(ExamineSpecWidget, self).__init__(parent)

        self.spec = spec
        self.psdict = {} # Dict for spectra plotting
        self.init_spec() 

        # Status Bar?
        if not status is None:
            self.statusBar = status

        # Line List?
        if llist is None:
            self.llist = {'Plot': False}
        else:
            self.llist = llist
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 150
        self.fig = Figure((8.0, 4.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)

        #self.canvas.mpl_connect('button_press_event', self.onclick)
        self.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
        self.canvas.setFocus()
        self.canvas.mpl_connect('key_press_event', self.on_key)
        self.canvas.mpl_connect('button_press_event', self.on_click)

        # Make two plots

        self.ax = self.fig.add_subplot(1,1,1)

        self.fig.subplots_adjust(hspace=0.1, wspace=0.1)
        
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.canvas)
        
        self.setLayout(vbox)

        # Draw on init
        self.on_draw()

    # Setup the spectrum plotting info
    def init_spec(self):
        #xy min/max
        xmin = np.min(self.spec.dispersion)
        xmax = np.max(self.spec.dispersion)
        ymed = np.median(self.spec.flux).value
        ymin = 0. - 0.1*ymed
        ymax = ymed * 1.5
        self.psdict['xmnx'] = [xmin,xmax]
        self.psdict['ymnx'] = [ymin,ymax]
        self.psdict['sv_xy'] = [ [xmin,xmax], [ymin,ymax] ]
        
    # Main Driver
    def on_key(self,event):

        flg = 0
        ## NAVIGATING
        if event.key in ['l','r','b','t','i','o','[',']','W','Z']:  # Set left
            flg = navigate(self.psdict,event)
        ## DOUBLET
        if event.key in ['C','M','V','A']:  # Set left
            wave = set_doublet(self, event)
            #print('wave = {:g},{:g}'.format(wave[0], wave[1]))
            self.ax.plot( [wave[0],wave[0]], self.psdict['ymnx'], '--', color='red')
            self.ax.plot( [wave[1],wave[1]], self.psdict['ymnx'], '--', color='red')
            flg = 2 # Layer
        if event.key == 'R': # Clear lines
            flg = 1 

        if flg==1: # Default is not to redraw
            self.on_draw()
        elif flg==2: # Layer (no clear)
            self.on_draw(replot=False) 

    # Click of main mouse button
    def on_click(self,event):
        # Draw line
        self.ax.plot( [event.xdata,event.xdata], self.psdict['ymnx'], ':', color='green')
        self.on_draw(replot=False) 

        # Print values
        print('button={:d}, x={:f}, y={:f}, xdata={:f}, ydata={:f}'.format(
            event.button, event.x, event.y, event.xdata, event.ydata))
        try:
            self.statusBar().showMessage('x,y = {:f}, {:f}'.format(event.xdata,event.ydata))
        except AttributeError:
            return

    def on_draw(self, replot=True):
        """ Redraws the figure
        """
        #
        if replot is True:
            self.ax.clear()        
            self.ax.plot(self.spec.dispersion, self.spec.flux, 'k-',drawstyle='steps-mid')
            self.ax.set_xlabel('Wavelength')
            self.ax.set_ylabel('Flux')

            # Spectral lines?
            if self.llist['Plot'] is True:
                z = self.llist['z']
                wvobs = np.array((1+z) * self.llist[self.llist['List']]['wrest'])
                gdwv = np.where( (wvobs > self.psdict['xmnx'][0]) &
                                 (wvobs < self.psdict['xmnx'][1]))[0]
                ylbl = self.psdict['ymnx'][1]-0.2*(self.psdict['ymnx'][1]-self.psdict['ymnx'][0])
                for kk in range(len(gdwv)): 
                    jj = gdwv[kk]
                    wrest = self.llist[self.llist['List']]['wrest'][jj]
                    lbl = self.llist[self.llist['List']]['name'][jj]
                    #QtCore.pyqtRemoveInputHook()
                    #xdb.set_trace()
                    #QtCore.pyqtRestoreInputHook()
                    # Plot
                    self.ax.plot(wrest*np.array([z+1,z+1]), self.psdict['ymnx'], 'b:')
                    # Label
                    self.ax.text(wrest*(z+1), ylbl, lbl, color='blue', rotation=90., size='small')
        
        self.ax.set_xlim(self.psdict['xmnx'])
        self.ax.set_ylim(self.psdict['ymnx'])



        # Draw
        self.canvas.draw()
    

        
# #####
class PlotLinesWidget(QtGui.QWidget):
    ''' Widget to set up spectral lines for plotting 

        13-Dec-2014 by JXP
    '''
    def __init__(self, parent=None, status=None):
        '''
        '''
        super(PlotLinesWidget, self).__init__(parent)

        if not status is None:
            self.statusBar = status
        
        # Create a dialog window for redshift
        z_label = QtGui.QLabel('z=')
        self.zbox = QtGui.QLineEdit()
        self.zbox.setText('0.')
        self.zbox.setMinimumWidth(50)
        self.connect(self.zbox, QtCore.SIGNAL('editingFinished ()'), self.setz)

        # Create the line list 
        list_label = QtGui.QLabel('Line Lists:')
        self.llist_widget = QtGui.QListWidget(self) 
        self.llist_widget.addItem('None')
        self.llist_widget.addItem('grb.lst')
        self.llist_widget.setCurrentRow(0)
        self.llist_widget.currentItemChanged.connect(self.on_list_change)

        self.llist = {} # Dict for the line lists
        self.llist['Plot'] = False
        self.llist['z'] = 0.

        # Create the selectable line
        line_label = QtGui.QLabel('Set Line:')
        self.scrollArea = QtGui.QScrollArea()

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(z_label)
        vbox.addWidget(self.zbox)
        vbox.addWidget(list_label)
        vbox.addWidget(self.llist_widget)
        
        self.setLayout(vbox)

    def on_list_change(self,curr,prev):
        llist = str(curr.text())
        # Print
        try:
            self.statusBar().showMessage('You chose: {:s}'.format(llist))
        except AttributeError:
            print('You chose: {:s}'.format(curr.text()))

        # Set line list
        self.llist['List'] = llist
        if llist == 'None':
            self.llist['Plot'] = False
        else:
            self.llist['Plot'] = True
            # Load?
            if not (llist in self.llist):
                line_file = xa_path+'/data/spec_lines/'+llist
                llist_cls = xspec.abs_line.Abs_Line_List(line_file)
                self.llist[llist] = llist_cls.data
            # Try to draw
            try:
                self.spec_widg.on_draw()
            except AttributeError:
                return

    def setz(self):
        sstr = unicode(self.zbox.text())
        try:
            self.llist['z'] = float(sstr)
        except ValueError:
            try:
                self.statusBar().showMessage('ERROR: z Input must be a float! Try again..')
            except AttributeError:
                print('ERROR: z Input must be a float! Try again..')
            self.zbox.setText(str(self.llist['z']))
            return
            
        # Report
        try:
            self.statusBar().showMessage('z = {:g}'.format(self.llist['z']))
        except AttributeError:
            print('z = {:g}'.format(self.llist['z']))

        # Try to draw
        try:
            self.spec_widg.on_draw()
        except AttributeError:
            return

# Plot Doublet
def set_doublet(iself,event):
    ''' Set z and plot doublet
    '''
    wv_dict = {'C': (1548.195, 1550.770, 'CIV'), 'M': (2796.352, 2803.531, 'MgII')}
    wrest = wv_dict[event.key]
    #QtCore.pyqtRemoveInputHook()
    #xdb.set_trace()
    #QtCore.pyqtRestoreInputHook()

    # Set z
    iself.zabs = event.xdata/wrest[0] - 1.
    try:
        iself.statusBar().showMessage('z = {:g} for {:s}'.format(iself.zabs, wrest[2]))
    except AttributeError:
        print('z = {:g} for {:s}'.format(iself.zabs, wrest[2]))

    return np.array(wrest[0:2])*(1.+iself.zabs)

# Navigate
def navigate(psdict,event):
    ''' Method to Navigate spectrum
    '''
    if event.key == 'l':  # Set left
        psdict['xmnx'][0] = event.xdata
    elif event.key == 'r':  # Set Right
        psdict['xmnx'][1] = event.xdata
    elif event.key == 'b':  # Set Bottom
        psdict['ymnx'][0] = event.ydata
    elif event.key == 't':  # Set Top
        psdict['ymnx'][1] = event.ydata
    elif event.key == 'i':  # Zoom in (and center)
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        deltx = (psdict['xmnx'][1]-psdict['xmnx'][0])/4.
        psdict['xmnx'] = [event.xdata-deltx, event.xdata+deltx]
    elif event.key == 'o':  # Zoom in (and center)
        deltx = psdict['xmnx'][1]-psdict['xmnx'][0]
        psdict['xmnx'] = [event.xdata-deltx, event.xdata+deltx]
    elif event.key in ['[',']']:  # Pan 
        center = (psdict['xmnx'][1]+psdict['xmnx'][0])/2.
        deltx = (psdict['xmnx'][1]-psdict['xmnx'][0])/2.
        if event.key == '[':
            new_center = center - deltx
        else:
            new_center = center + deltx
        psdict['xmnx'] = [new_center-deltx, new_center+deltx]
    elif event.key == 'W': # Reset the Window
        psdict['xmnx'] = psdict['sv_xy'][0]
        psdict['ymnx'] = psdict['sv_xy'][1]
    elif event.key == 'Z': # Zero
        psdict['ymnx'][0] = 0.
    else:
        if not (event.key in ['shift']):
            rstr = 'Key {:s} not supported.'.format(event.key)
            print(rstr)
        return 0
    return 1





# ################
# TESTING
if __name__ == "__main__":
    from xastropy import spec as xspec

    flg_fig = 0 
    flg_fig += 2**0  # ExamineSpecWidget
    #flg_fig += 2**1  # PlotLinesWidget

    # ExamineSpec
    if (flg_fig % 2) == 1:
        app = QtGui.QApplication(sys.argv)
        spec_fil = '/u/xavier/Keck/HIRES/RedData/PH957/PH957_f.fits'
        spec = xspec.readwrite.readspec(spec_fil)
        app.setApplicationName('XSpec')
        main = ExamineSpecWidget(spec)
        main.show()
        sys.exit(app.exec_())

    # PltLineWidget
    if (flg_fig % 2**2) >= 2**1:
        app = QtGui.QApplication(sys.argv)
        app.setApplicationName('PltLine')
        main = PlotLinesWidget()
        main.show()
        sys.exit(app.exec_())
