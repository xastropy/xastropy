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

from astropy.table.table import Table

from xastropy import spec as xspec 
from xastropy.xutils import xdebug as xdb

xa_path = imp.find_module('xastropy')[1]

class ExamineSpecWidget(QtGui.QWidget):
    ''' Widget to plot a spectrum and interactively
        fiddle about.  Akin to XIDL/x_specplot.pro

        12-Dec-2014 by JXP
    '''
    def __init__(self, spec, parent=None, status=None, llist=None,
                 abs_sys=None):
        '''
        spec = Spectrum1D
        '''
        super(ExamineSpecWidget, self).__init__(parent)

        self.spec = spec
        self.abs_sys = abs_sys
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

        #

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
        print('button={:d}, x={:f}, y={:f}, xdata={:f}, ydata={:f}'.format(
            event.button, event.x, event.y, event.xdata, event.ydata))
        if event.button == 1: # Draw line
            self.ax.plot( [event.xdata,event.xdata], self.psdict['ymnx'], ':', color='green')
            self.on_draw(replot=False) 
    
            # Print values
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
                ylbl = self.psdict['ymnx'][1]-0.2*(self.psdict['ymnx'][1]-self.psdict['ymnx'][0])
                z = self.llist['z']
                wvobs = np.array((1+z) * self.llist[self.llist['List']]['wrest'])
                gdwv = np.where( (wvobs > self.psdict['xmnx'][0]) &
                                 (wvobs < self.psdict['xmnx'][1]))[0]
                for kk in range(len(gdwv)): 
                    jj = gdwv[kk]
                    wrest = self.llist[self.llist['List']]['wrest'][jj]
                    lbl = self.llist[self.llist['List']]['name'][jj]
                    # Plot
                    self.ax.plot(wrest*np.array([z+1,z+1]), self.psdict['ymnx'], 'b:')
                    # Label
                    self.ax.text(wrest*(z+1), ylbl, lbl, color='blue', rotation=90., size='small')

            # Abs Sys?
            if not self.abs_sys is None:
                ylbl = self.psdict['ymnx'][0]+0.2*(self.psdict['ymnx'][1]-self.psdict['ymnx'][0])
                for abs_sys in self.abs_sys:
                    wrest = np.array(abs_sys.lines.keys()) 
                    wvobs = wrest * (abs_sys.zabs+1)
                    gdwv = np.where( ((wvobs+5) > self.psdict['xmnx'][0]) &  # Buffer for region
                                    ((wvobs-5) < self.psdict['xmnx'][1]))[0]
                    for kk in range(len(gdwv)): 
                        jj = gdwv[kk]
                        #QtCore.pyqtRemoveInputHook()
                        #xdb.set_trace()
                        #QtCore.pyqtRestoreInputHook()
                        # Paint spectrum red
                        wvlim = wvobs[jj]*(1 + abs_sys.lines[wrest[jj]].analy['VLIM']/3e5)
                        pix = np.where( (self.spec.dispersion > wvlim[0]) & (self.spec.dispersion < wvlim[1]))[0]
                        self.ax.plot(self.spec.dispersion[pix], self.spec.flux[pix], 'r-',drawstyle='steps-mid')
                        # Label
                        lbl = abs_sys.lines[wrest[jj]].analy['IONNM']
                        self.ax.text(wvobs[jj], ylbl, lbl, color='red', rotation=90., size='small')
        
        # Reset window limits
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
            self.zbox.setText('{:.5f}'.format(self.llist['z']))
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

# #####
class SelectLineWidget(QtGui.QDialog):
    ''' Widget to select a spectral line
    inp: string or dict
      Input line list

    15-Dec-2014 by JXP
    '''
    def __init__(self, inp, parent=None):
        '''
        '''
        super(SelectLineWidget, self).__init__(parent)

        # Line list dict
        if isinstance(inp,Table):
            lines = inp

        self.resize(250, 800)

        # Create the line list 
        line_label = QtGui.QLabel('Lines:')
        self.lines_widget = QtGui.QListWidget(self) 
        self.lines_widget.addItem('None')
        self.lines_widget.setCurrentRow(0)

        #xdb.set_trace()
        # Loop on lines (could put a preferred list first)
        nlin = len(lines['wrest'])
        for ii in range(nlin):
            self.lines_widget.addItem('{:s} :: {:.4f}'.format(lines['name'][ii],
                                                         lines['wrest'][ii]))
        self.lines_widget.currentItemChanged.connect(self.on_list_change)
        #self.scrollArea = QtGui.QScrollArea()

        # Quit
        qbtn = QtGui.QPushButton('Quit', self)
        qbtn.clicked.connect(self.close)

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(line_label)
        vbox.addWidget(self.lines_widget)
        vbox.addWidget(qbtn)
        
        self.setLayout(vbox)

    def on_list_change(self,curr,prev):
        self.line = str(curr.text())
        # Print
        print('You chose: {:s}'.format(curr.text()))

# #####
class AbsSysWidget(QtGui.QWidget):
    ''' Widget to organize AbsSys along a given sightline

    Parameters:
    -----------
    abssys_list: List
      String list of abssys files

    16-Dec-2014 by JXP
    '''
    def __init__(self, abssys_list, parent=None):
        '''
        '''
        super(AbsSysWidget, self).__init__(parent)

        #if not status is None:
        #    self.statusBar = status
        self.abssys_list = abssys_list
        
        # Create the line list 
        list_label = QtGui.QLabel('Abs Systems:')
        self.abslist_widget = QtGui.QListWidget(self) 
        self.abslist_widget.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.abslist_widget.addItem('None')
        #self.abslist_widget.addItem('Test')

        self.all_items = []
        for abssys_fil in self.abssys_list:
            ipos0 = abssys_fil.rfind('/') + 1
            ipos1 = abssys_fil.rfind('.fits')
            self.all_items.append( abssys_fil[ipos0:ipos1] )
            self.abslist_widget.addItem(abssys_fil[ipos0:ipos1] )

        self.abslist_widget.setCurrentRow(0)
        self.abslist_widget.itemSelectionChanged.connect(self.on_list_change)

        # List for Abs Sys
        self.abs_sys = []
        self.items = []


        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(list_label)
        vbox.addWidget(self.abslist_widget)
        
        self.setLayout(vbox)

    # ##
    def on_list_change(self):
        
        from xastropy.igm.abs_sys.lls_utils import LLS_System

        items = self.abslist_widget.selectedItems()
        # Empty the list
        if len(self.abs_sys) > 0:
            for ii in range(len(self.abs_sys)-1,-1,-1):
                self.abs_sys.pop(ii)
        # Load up abs_sys (as need be)
        new_items = []
        for item in items:
            txt = item.text()
            # Dummy
            if txt == 'None':
                continue
            print('Including {:s} in the list'.format(txt))
            # Using LLS for now.  Might change to generic
            new_items.append(txt)
            ii = self.all_items.index(txt)
            self.abs_sys.append(LLS_System.from_absid_fil(self.abssys_list[ii]))

        # Pass back
        self.items = new_items
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()


# ######
# Plot Doublet
def set_doublet(iself,event):
    ''' Set z and plot doublet
    '''
    wv_dict = {'C': (1548.195, 1550.770, 'CIV'), 'M': (2796.352, 2803.531, 'MgII')}
    wrest = wv_dict[event.key]

    # Set z
    iself.zabs = event.xdata/wrest[0] - 1.
    try:
        iself.statusBar().showMessage('z = {:g} for {:s}'.format(iself.zabs, wrest[2]))
    except AttributeError:
        print('z = {:g} for {:s}'.format(iself.zabs, wrest[2]))

    return np.array(wrest[0:2])*(1.+iself.zabs)

# ######
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
    #flg_fig += 2**0  # ExamineSpecWidget
    #flg_fig += 2**1  # PlotLinesWidget
    #flg_fig += 2**2  # SelectLineWidget
    flg_fig += 2**3  # AbsSysWidget

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

    # SelectLineWidget
    if (flg_fig % 2**3) >= 2**2:
        line_file = xa_path+'/data/spec_lines/grb.lst'
        llist_cls = xspec.abs_line.Abs_Line_List(line_file)

        app = QtGui.QApplication(sys.argv)
        app.setApplicationName('SelectLine')
        main = SelectLineWidget(llist_cls.data)
        main.show()
        app.exec_()
        print(main.line)
        sys.exit()
    
    # AbsSys Widget
    if (flg_fig % 2**4) >= 2**3:
        abs_fil = '/Users/xavier/paper/LLS/Optical/Data/Analysis/MAGE/SDSSJ1004+0018_z2.746_id.fits'
        abs_fil2 = '/Users/xavier/paper/LLS/Optical/Data/Analysis/MAGE/SDSSJ2319-1040_z2.675_id.fits'
        app = QtGui.QApplication(sys.argv)
        app.setApplicationName('AbsSys')
        main = AbsSysWidget([abs_fil,abs_fil2])
        main.show()
        sys.exit(app.exec_())

