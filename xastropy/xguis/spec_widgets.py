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
from astropy import constants as const
from astropy import units as u
from astropy.units import Quantity
u.def_unit(['mAA', 'milliAngstrom'], 0.001 * u.AA, namespace=globals()) # mA
from astropy.nddata import StdDevUncertainty

from specutils.spectrum1d import Spectrum1D

from linetools.spectra import io as lsi
from linetools.lists.linelist import LineList

from xastropy import spec as xspec 
from xastropy import stats as xstats
from xastropy import xutils 
from xastropy.xutils import xdebug as xdb
from xastropy.plotting import utils as xputils
from xastropy.igm.abs_sys import abssys_utils as xiaa
from xastropy.igm.abs_sys.lls_utils import LLS_System
from xastropy.xguis import utils as xguiu

xa_path = imp.find_module('xastropy')[1]

# class ExamineSpecWidget
# class PlotLinesWidget

class ExamineSpecWidget(QtGui.QWidget):
    ''' Widget to plot a spectrum and interactively
        fiddle about.  Akin to XIDL/x_specplot.pro

        12-Dec-2014 by JXP
    '''
    def __init__(self, ispec, parent=None, status=None, llist=None,
                 abs_sys=None, norm=True, second_file=None, zsys=None):
        '''
        spec = Spectrum1D
        '''
        super(ExamineSpecWidget, self).__init__(parent)

        # Spectrum
        spec, spec_fil = read_spec(ispec, second_file=second_file)
        self.orig_spec = spec # For smoothing
        self.spec = self.orig_spec 

        # Abs Systems
        if abs_sys is None:
            self.abs_sys = []
        else:
            self.abs_sys = abs_sys
        self.norm = norm
        self.psdict = {} # Dict for spectra plotting
        self.adict = {}  # Dict for analysis
        self.init_spec() 
        self.xval = None # Used with velplt

        # Status Bar?
        if not status is None:
            self.statusBar = status

        # Line List?
        if llist is None:
            self.llist = {'Plot': False, 'List': 'None', 'z': 0.}
        else:
            self.llist = llist

        # zsys
        if not zsys is None:
            self.llist['z'] = zsys
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 150 # 150
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
        xmin = np.min(self.spec.dispersion).value
        xmax = np.max(self.spec.dispersion).value
        ymed = np.median(self.spec.flux).value
        ymin = 0. - 0.1*ymed
        ymax = ymed * 1.5
        #
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        self.psdict['xmnx'] = np.array([xmin,xmax])
        self.psdict['ymnx'] = [ymin,ymax]
        self.psdict['sv_xy'] = [ [xmin,xmax], [ymin,ymax] ]
        self.psdict['nav'] = navigate(0,0,init=True)
        # Analysis dict
        self.adict['flg'] = 0 # Column density flag
        
        
    # Main Driver
    def on_key(self,event):

        flg = -1

        ## NAVIGATING
        if event.key in self.psdict['nav']: 
            flg = navigate(self.psdict,event)

        ## DOUBLETS
        if event.key in ['C','M','X','4','8','B']:  # Set left
            wave = set_doublet(self, event)
            #print('wave = {:g},{:g}'.format(wave[0], wave[1]))
            self.ax.plot( [wave[0],wave[0]], self.psdict['ymnx'], '--', color='red')
            self.ax.plot( [wave[1],wave[1]], self.psdict['ymnx'], '--', color='red')
            flg = 2 # Layer

        ## SMOOTH
        if event.key == 'S':
            self.spec = self.spec.box_smooth(2)
            flg = 1 
        if event.key == 'U':
            self.spec = self.orig_spec
            flg = 1 

        ## Lya Profiles 
        if event.key in ['D', 'R']:
            # Set NHI
            if event.key == 'D':
                NHI = 20.3
            elif event.key == 'R':
                NHI = 19.0
            zlya = event.xdata/1215.6701 - 1.
            self.llist['z'] = zlya
            # Generate Lya profile
            lya_line = xspec.lines_utils.AbsLine(1215.6701*u.AA)
            lya_line.z = zlya
            lya_line.attrib['N'] = NHI
            lya_line.attrib['b'] = 30.
            self.lya_line = xspec.voigt.voigt_model(self.spec.dispersion, lya_line, Npix=3.)
            self.adict['flg'] = 4
            flg = 1 

        ## ANALYSIS:  EW, AODM column density
        if event.key in ['N', 'E', '$']: 
            # If column check for line list
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            if (event.key in ['N','E']) & (self.llist['List'] == 'None'):
                print('xspec: Choose a Line list first!')
                try:
                    self.statusBar().showMessage('Choose a Line list first!')
                except AttributeError:
                    pass
                self.adict['flg'] = 0
                return
            flg = 1 
                
            if self.adict['flg'] == 0:
                self.adict['wv_1'] = event.xdata # wavelength
                self.adict['C_1'] = event.ydata # continuum
                self.adict['flg'] = 1 # Plot dot
            else:
                self.adict['wv_2'] = event.xdata # wavelength
                self.adict['C_2'] = event.ydata # continuum
                self.adict['flg'] = 2 # Ready to plot + print

                # Sort em + make arrays
                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()
                iwv = np.array(sorted([self.adict['wv_1'], self.adict['wv_2']])) * self.spec.wcs.unit
                ic = np.array(sorted([self.adict['C_1'], self.adict['C_2']]))

                # Calculate the continuum (linear fit)
                param = np.polyfit(iwv, ic, 1)
                cfunc = np.poly1d(param)
                conti = cfunc(self.spec.dispersion)

                if event.key == '$': # Simple stats
                    pix = self.spec.pix_minmax(iwv)[0]
                    mean = np.mean(self.spec.flux[pix])
                    median = np.median(self.spec.flux[pix])
                    stdv = np.std(self.spec.flux[pix]-conti[pix])
                    S2N = median / stdv
                    mssg = 'Mean={:g}, Median={:g}, S/N={:g}'.format(mean,median,S2N)
                else:
                    # Find the spectral line (or request it!)
                    rng_wrest = iwv / (self.llist['z']+1)
                    #QtCore.pyqtRemoveInputHook()
                    #xdb.set_trace()
                    #QtCore.pyqtRestoreInputHook()
                    gdl = np.where( (self.llist[self.llist['List']].wrest-rng_wrest[0]) *
                                    (self.llist[self.llist['List']].wrest-rng_wrest[1]) < 0.)[0]
                    if len(gdl) == 1:
                        wrest = self.llist[self.llist['List']].wrest[gdl[0]]
                    else:
                        if len(gdl) == 0: # Search through them all
                            gdl = np.arange(len(self.llist[self.llist['List']]))
                        sel_widg = SelectLineWidget(self.llist[self.llist['List']]._data[gdl])
                        sel_widg.exec_()
                        line = sel_widg.line
                        #wrest = float(line.split('::')[1].lstrip())
                        quant = line.split('::')[1].lstrip()
                        spltw = quant.split(' ')
                        wrest = Quantity(float(spltw[0]), unit=spltw[1])
                    # Units
                    if not hasattr(wrest,'unit'):
                        # Assume Ang
                        wrest = wrest * u.AA
    
                    # Generate the Spectral Line
                    from xastropy.spec.lines_utils import AbsLine
                    raise ValueError('DEPRECATE THE NEXT LINE! OR ADD LineList')
                    aline = AbsLine(wrest)
                    aline.analy['z'] = self.llist['z']
                    aline.spec = self.spec
    
                    # AODM
                    if event.key == 'N': 
                        # Calculate the velocity limits and load-up
                        aline.analy['VLIM'] = const.c.to('km/s') * (
                            ( iwv/(1+self.llist['z']) - wrest) / wrest )
        
                        # AODM
                        #QtCore.pyqtRemoveInputHook()
                        #xdb.set_trace()
                        #QtCore.pyqtRestoreInputHook()
                        aline.aodm(conti=conti)
                        mssg = 'Using '+ aline.__repr__()
                        mssg = mssg + ' ::  logN = {:g} +/- {:g}'.format(aline.attrib['logN'],
                                                                        aline.attrib['sig_logN'])
                    elif event.key == 'E':  #EW
                        aline.analy['WVMNX'] = iwv
                        aline.restew(conti=conti)
                        mssg = 'Using '+ aline.__repr__()
                        mssg = mssg + ' ::  EW = {:g} +/- {:g}'.format(aline.attrib['EW'].to(mAA),
                                                                        aline.attrib['sigEW'].to(mAA))
                # Display values
                try:
                    self.statusBar().showMessage(mssg)
                except AttributeError:
                    pass
                print(mssg)

                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()


        ## Velocity plot
        if event.key == 'v': 
            flg = 0
            from xastropy.xguis import spec_guis as xsgui
            z=self.llist['z']
            # Check for a match in existing list and use it if so
            if len(self.abs_sys) > 0:
                zabs = np.array([abs_sys.zabs for abs_sys in self.abs_sys])
                mt = np.where( np.abs(zabs-z) < 1e-4)[0]
            else:
                mt = []
            if len(mt) == 1:
                ini_abs_sys = self.abs_sys[mt[0]]
                outfil = ini_abs_sys.absid_file
                self.vplt_flg = 0 # Old one
                print('Using existing ID file {:s}'.format(outfil))
            else: 
                ini_abs_sys = None
                outfil = None
                self.vplt_flg = 1 # New one

            # Outfil
            if outfil is None:
                i0 = self.spec.filename.rfind('/')
                i1 = self.spec.filename.rfind('.')
                if i0 < 0:
                    path = './ID_LINES/'
                else:
                    path = self.spec.filename[0:i0]+'/ID_LINES/'
                outfil = path + self.spec.filename[i0+1:i1]+'_z'+'{:.4f}'.format(z)+'_id.fits'
                xutils.files.ensure_dir(outfil)
                self.outfil = outfil
                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()

            # Launch
            gui = xsgui.XVelPltGui(self.spec, z=z, outfil=outfil, llist=self.llist,
                                   abs_sys=ini_abs_sys, norm=self.norm, sel_wv=self.xval)
            gui.exec_()
            if gui.flg_quit == 0: # Quit without saving (i.e. discarded)
                self.vplt_flg = 0 
            else:
                # Push to Abs_Sys
                if len(mt) == 1:
                    self.abs_sys[mt[0]] = gui.abs_sys
                else:
                    self.abs_sys.append(gui.abs_sys)
                    print('Adding new abs system')
            # Redraw
            flg=1

        # Dummy keys
        if event.key in ['shift', 'control', 'shift+super', 'super+shift']:
            flg = 0

        # Draw
        if flg==1: # Default is not to redraw
            self.on_draw()
        elif flg==2: # Layer (no clear)
            self.on_draw(replot=False) 
        elif flg==-1: # Layer (no clear)
            try:
                self.statusBar().showMessage('Not a valid key!  {:s}'.format(event.key))
            except AttributeError:
                pass

    # Click of main mouse button
    def on_click(self,event):
        try:
            print('button={:d}, x={:f}, y={:f}, xdata={:f}, ydata={:f}'.format(
                event.button, event.x, event.y, event.xdata, event.ydata))
        except ValueError:
            print('Out of bounds')
            return
        if event.button == 1: # Draw line
            self.xval = event.xdata
            self.ax.plot( [event.xdata,event.xdata], self.psdict['ymnx'], ':', color='green')
            self.on_draw(replot=False) 
    
            # Print values
            try:
                self.statusBar().showMessage('x,y = {:f}, {:f}'.format(event.xdata,event.ydata))
            except AttributeError:
                return

    # ######
    def on_draw(self, replot=True):
        """ Redraws the spectrum
        """
        #

        if replot is True:
            self.ax.clear()        
            self.ax.plot(self.spec.dispersion, self.spec.flux, 'k-',drawstyle='steps-mid')
            try: 
                self.ax.plot(self.spec.dispersion, self.spec.sig, 'r:')
            except ValueError:
                pass
            self.ax.set_xlabel('Wavelength')
            self.ax.set_ylabel('Flux')

            # Spectral lines?
            if self.llist['Plot'] is True:
                ylbl = self.psdict['ymnx'][1]-0.2*(self.psdict['ymnx'][1]-self.psdict['ymnx'][0])
                z = self.llist['z']
                wvobs = np.array((1+z) * self.llist[self.llist['List']].wrest)
                gdwv = np.where( (wvobs > self.psdict['xmnx'][0]) &
                                 (wvobs < self.psdict['xmnx'][1]))[0]
                for kk in range(len(gdwv)): 
                    jj = gdwv[kk]
                    wrest = self.llist[self.llist['List']].wrest[jj].value
                    lbl = self.llist[self.llist['List']].name[jj]
                    # Plot
                    self.ax.plot(wrest*np.array([z+1,z+1]), self.psdict['ymnx'], 'b--')
                    # Label
                    self.ax.text(wrest*(z+1), ylbl, lbl, color='blue', rotation=90., size='small')

            # Abs Sys?
            if not self.abs_sys is None:
                ylbl = self.psdict['ymnx'][0]+0.2*(self.psdict['ymnx'][1]-self.psdict['ymnx'][0])
                clrs = ['red', 'green', 'cyan', 'orange', 'gray', 'purple']*10
                for abs_sys in self.abs_sys:
                    ii = self.abs_sys.index(abs_sys)
                    kwrest = np.array(abs_sys.lines.keys()) 
                    wvobs = kwrest * (abs_sys.zabs+1) * u.AA
                    gdwv = np.where( ((wvobs.value+5) > self.psdict['xmnx'][0]) &  # Buffer for region
                                    ((wvobs.value-5) < self.psdict['xmnx'][1]))[0]
                    for kk in range(len(gdwv)): 
                        jj = gdwv[kk]
                        if abs_sys.lines[kwrest[jj]].analy['FLG_ANLY'] == 0:
                            continue
                        # Paint spectrum red
                        wvlim = wvobs[jj]*(1 + abs_sys.lines[kwrest[jj]].analy['VLIM']/3e5)
                        pix = np.where( (self.spec.dispersion > wvlim[0]) & (self.spec.dispersion < wvlim[1]))[0]
                        self.ax.plot(self.spec.dispersion[pix], self.spec.flux[pix], '-',drawstyle='steps-mid',
                                     color=clrs[ii])
                        # Label
                        lbl = abs_sys.lines[kwrest[jj]].analy['IONNM']+' z={:g}'.format(abs_sys.zabs)
                        self.ax.text(wvobs[jj].value, ylbl, lbl, color=clrs[ii], rotation=90., size='x-small')
            # Analysis? EW, Column
            if self.adict['flg'] == 1:
                self.ax.plot(self.adict['wv_1'], self.adict['C_1'], 'go')
            elif self.adict['flg'] == 2:
                self.ax.plot([self.adict['wv_1'], self.adict['wv_2']],
                             [self.adict['C_1'], self.adict['C_2']], 'g--', marker='o')
                self.adict['flg'] = 0
            # Lya line?
            if self.adict['flg'] == 4:
                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()
                self.ax.plot(self.spec.dispersion, self.lya_line.flux, color='green')
        
        # Reset window limits
        self.ax.set_xlim(self.psdict['xmnx'])
        self.ax.set_ylim(self.psdict['ymnx'])



        # Draw
        self.canvas.draw()

    # Notes on usage
    def help_notes():
        doublets = [ 'Doublets --------',
                     'C: CIV',
                     'M: MgII', 
                     'O: OVI',
                     '8: NeVIII',
                     'B: Lyb/Lya'
                     ]
        analysis = [ 'Analysis --------',
                     'N/N: Column density (AODM)',
                     'E/E: EW (boxcar)',
                     '$/$: stats on spectrum'
                     ]
    

        
# #####
class PlotLinesWidget(QtGui.QWidget):
    ''' Widget to set up spectral lines for plotting 

        13-Dec-2014 by JXP
    '''
    def __init__(self, parent=None, status=None, init_llist=None, init_z=None):
        '''
        '''
        super(PlotLinesWidget, self).__init__(parent)

        # Initialize
        if not status is None:
            self.statusBar = status
        if init_z is None:
            init_z = 0.
            
        
        # Create a dialog window for redshift
        z_label = QtGui.QLabel('z=')
        self.zbox = QtGui.QLineEdit()
        self.zbox.z_frmt = '{:.7f}'
        self.zbox.setText(self.zbox.z_frmt.format(init_z))
        self.zbox.setMinimumWidth(50)
        self.connect(self.zbox, QtCore.SIGNAL('editingFinished ()'), self.setz)

        # Create the line list 
        self.lists = ['None', 'ISM', 'Strong', 'H2'] 
        #'grb.lst', 'dla.lst', 'lls.lst', 'subLLS.lst', 
#                      'lyman.lst', 'Dlyman.lst', 'gal_vac.lst', 'ne8.lst',
#                      'lowz_ovi.lst', 'casbah.lst', 'H2.lst']
        list_label = QtGui.QLabel('Line Lists:')
        self.llist_widget = QtGui.QListWidget(self) 
        for ilist in self.lists:
            self.llist_widget.addItem(ilist)
        self.llist_widget.setCurrentRow(0)
        self.llist_widget.currentItemChanged.connect(self.on_list_change)
        self.llist_widget.setMaximumHeight(100)

        # Input line list?
        if init_llist is None:
            self.llist = {} # Dict for the line lists
            self.llist['Plot'] = False
            self.llist['z'] = 0.
            self.llist['List'] = 'None'
        else: # Fill it all up and select
            self.llist = init_llist
            if not init_llist['List'] in self.lists:
                self.lists.append(init_llist['List'])
                self.llist_widget.addItem(init_llist['List'])
                self.llist_widget.setCurrentRow(len(self.lists)-1)
            else:
                idx = self.lists.index(init_llist['List'])
                self.llist_widget.setCurrentRow(idx)
            try:
                self.zbox.setText(self.zbox.z_frmt.format(init_llist['z']))
            except KeyError:
                pass

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(z_label)
        vbox.addWidget(self.zbox)
        vbox.addWidget(list_label)
        vbox.addWidget(self.llist_widget)
        
        self.setLayout(vbox)
        self.setMaximumHeight(200)

    def on_list_change(self,curr,prev):
        llist = str(curr.text())
        # Print
        try:
            self.statusBar().showMessage('You chose: {:s}'.format(llist))
        except AttributeError:
            print('You chose: {:s}'.format(curr.text()))

        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        self.llist = set_llist(llist,in_dict=self.llist)

        # Try to draw
        if self.llist['Plot'] is True:
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
            self.zbox.setText(self.zbox.z_frmt.format(self.llist['z']))
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
    inp: string or dict or Table
      Input line list

    15-Dec-2014 by JXP
    '''
    def __init__(self, inp, parent=None):
        '''
        '''
        super(SelectLineWidget, self).__init__(parent)

        # Line list Table
        if isinstance(inp,Table):
            lines = inp
        else:
            raise ValueError('SelectLineWidget: Wrong type of input')

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
class SelectedLinesWidget(QtGui.QWidget):
    ''' Widget to show and enable lines to be selected
    inp: LineList
      Input LineList

    24-Dec-2014 by JXP
    '''
    def __init__(self, inp, parent=None, init_select=None, plot_widget=None):
        '''
        '''
        super(SelectedLinesWidget, self).__init__(parent)

        # Line list Table
        if isinstance(inp,LineList):
            self.lines = inp._data
            self.llst = inp
        elif isinstance(inp,Table):
            raise ValueError('SelectLineWidget: DEPRECATED')
        else:
            raise ValueError('SelectLineWidget: Wrong type of input')

        self.plot_widget = plot_widget

        # Create the line list 
        line_label = QtGui.QLabel('Lines:')
        self.lines_widget = QtGui.QListWidget(self) 
        self.lines_widget.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)

        # Initialize list
        self.item_flg = 0 
        self.init_list()

        # Initial selection
        if init_select is None:
            self.selected = [0]
        else:
            self.selected = init_select

        for iselect in self.selected:
            self.lines_widget.item(iselect).setSelected(True)

        self.lines_widget.scrollToItem( self.lines_widget.item( self.selected[0] ) )

        # Events
        #self.lines_widget.itemClicked.connect(self.on_list_change)
        self.lines_widget.itemSelectionChanged.connect(self.on_item_change)

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(line_label)
        vbox.addWidget(self.lines_widget)
        
        self.setLayout(vbox)

    def init_list(self):
        nlin = len(self.lines['wrest'])
        for ii in range(nlin):
            self.lines_widget.addItem('{:s} :: {:.3f}'.format(self.lines['name'][ii],
                                                         self.lines['wrest'][ii].value))

    def on_item_change(self): #,item):
        # For big changes
        if self.item_flg == 1:
            return
        all_items = [self.lines_widget.item(ii) for ii in range(self.lines_widget.count())]
        sel_items = self.lines_widget.selectedItems()
        self.selected = [all_items.index(isel) for isel in sel_items]
        self.selected.sort()

        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

        # Update llist
        try:
            self.plot_widget.llist['show_line'] = self.selected
        except AttributeError:
            return
        else:
            self.plot_widget.on_draw()

    def on_list_change(self,lines): 
        # Clear
        self.item_flg = 1
        self.lines_widget.clear()
        # Initialize
        self.lines = lines
        self.init_list()
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        # Set selected
        for iselect in self.selected:
            self.lines_widget.item(iselect).setSelected(True)
        self.lines_widget.scrollToItem( self.lines_widget.item( self.selected[0] ) )
        self.item_flg = 0

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

        # Lists
        self.abs_sys = []
        self.items = []
        self.all_items = []
        self.all_abssys = []
        for abssys_fil in self.abssys_list:
            self.all_abssys.append(LLS_System.from_absid_fil(abssys_fil))
            self.add_item(abssys_fil)

        self.abslist_widget.setCurrentRow(0)
        self.abslist_widget.itemSelectionChanged.connect(self.on_list_change)

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(list_label)

        # Buttons
        buttons = QtGui.QWidget()
        self.refine_button = QtGui.QPushButton('Refine', self)
        #self.refine_button.clicked.connect(self.refine) # CONNECTS TO A PARENT
        reload_btn = QtGui.QPushButton('Reload', self)
        reload_btn.clicked.connect(self.reload)
        hbox1 = QtGui.QHBoxLayout()
        hbox1.addWidget(self.refine_button)
        hbox1.addWidget(reload_btn)
        buttons.setLayout(hbox1)
        vbox.addWidget(buttons)

        vbox.addWidget(self.abslist_widget)
        self.setLayout(vbox)

    # ##
    def on_list_change(self):
        
        items = self.abslist_widget.selectedItems()
        # Empty the list
        #self.abs_sys = []
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
            self.abs_sys.append(self.all_abssys[ii])

        # Pass back
        self.items = new_items
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

    def add_fil(self,abssys_fil):
        self.abssys_list.append( abssys_fil )
        self.add_item(abssys_fil)

    def add_item(self,abssys_fil):
        ipos0 = abssys_fil.rfind('/') + 1
        ipos1 = abssys_fil.rfind('.fits')
        self.all_items.append( abssys_fil[ipos0:ipos1] )
        self.abslist_widget.addItem(abssys_fil[ipos0:ipos1] )

    def reload(self):
        print('AbsSysWidget: Reloading systems..')
        self.all_abssys = []
        for abssys_fil in self.abssys_list:
            self.all_abssys.append(LLS_System.from_absid_fil(abssys_fil))
            #self.add_item(abssys_fil)
        self.on_list_change()

# ######################
class VelPlotWidget(QtGui.QWidget):
    ''' Widget for a velocity plot with interaction.

        19-Dec-2014 by JXP
    '''
    def __init__(self, ispec, z=None, parent=None, llist=None, norm=True,
                 vmnx=[-300., 300.], abs_sys=None):
        '''
        spec = Spectrum1D
        Norm: Bool (False)
          Normalized spectrum?
        abs_sys: AbsSystem
          Absorption system class
        '''
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        super(VelPlotWidget, self).__init__(parent)

        # Initialize
        spec, spec_fil = read_spec(ispec)
        
        self.spec = spec
        self.spec_fil = spec_fil
        self.z = z
        self.vmnx = vmnx
        self.norm = norm

        # Abs_System 
        self.abs_sys = abs_sys
        if self.abs_sys is None:
            self.abs_sys = xiaa.Generic_System(None)
            self.abs_sys.zabs = self.z
        else:
            self.z = self.abs_sys.zabs
            # Line list
            if llist is None:
                try:
                    lwrest = self.abs_sys.lines.keys()
                except AttributeError:
                    lwrest = None
                if not lwrest is None:
                    llist = set_llist(lwrest)


        self.psdict = {} # Dict for spectra plotting
        self.psdict['xmnx'] = self.vmnx
        self.psdict['ymnx'] = [-0.1, 1.1]
        self.psdict['nav'] = navigate(0,0,init=True)

        # Status Bar?
        #if not status is None:
        #    self.statusBar = status

        # Line List
        if llist is None:
            self.llist = set_llist('lls.lst')
        else:
            self.llist = llist
        self.llist['z'] = self.z

        # Indexing for line plotting
        self.idx_line = 0

        self.init_lines()
        
        # Create the mpl Figure and FigCanvas objects. 
        #
        self.dpi = 150
        self.fig = Figure((8.0, 4.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)

        self.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
        self.canvas.setFocus()
        self.canvas.mpl_connect('key_press_event', self.on_key)
        self.canvas.mpl_connect('button_press_event', self.on_click)

        # Sub_plots
        self.sub_xy = [3,4]

        self.fig.subplots_adjust(hspace=0.0, wspace=0.1)
        
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.canvas)
        
        self.setLayout(vbox)

        # Draw on init
        self.on_draw()

    # Load them up for display
    def init_lines(self):
        wvmin = np.min(self.spec.dispersion)
        wvmax = np.max(self.spec.dispersion)
        #
        wrest = self.llist[self.llist['List']].wrest
        wvobs = (1+self.z) * wrest
        gdlin = np.where( (wvobs > wvmin) & (wvobs < wvmax) )[0]
        self.llist['show_line'] = gdlin

        existing_lines = self.abs_sys.lines.keys()

        # Update/generate lines
        for idx in gdlin:
            # Generate?
            kwrest = wrest[idx].value
            if not kwrest in existing_lines:
                self.abs_sys.lines[kwrest] = xspec.analysis.Spectral_Line(kwrest)
                print('VelPlot: Generating line {:g}'.format(kwrest))
                self.abs_sys.lines[kwrest].analy['VLIM'] = np.array([self.vmnx[0]/2.,
                                                                self.vmnx[1]/2.])
                self.abs_sys.lines[kwrest].analy['FLG_ANLY'] = 2 # Init to ok
            # Spec file
            if not self.spec_fil is None:
                self.abs_sys.lines[kwrest].analy['DATFIL'] = self.spec_fil
            
        
    # Key stroke 
    def on_key(self,event):

        # Init
        rescale = True
        fig_clear = False
        wrest = None
        flg = 0
        sv_idx = self.idx_line

        ## Change rows/columns
        if event.key == 'k':
            self.sub_xy[0] = max(0, self.sub_xy[0]-1)
        if event.key == 'K':
            self.sub_xy[0] = self.sub_xy[0]+1
        if event.key == 'c':
            self.sub_xy[1] = max(0, self.sub_xy[1]-1)
        if event.key == 'C':
            self.sub_xy[1] = max(0, self.sub_xy[1]+1)

        ## NAVIGATING
        if event.key in self.psdict['nav']: 
            flg = navigate(self.psdict,event)
        if event.key == '-':
            self.idx_line = max(0, self.idx_line-self.sub_xy[0]*self.sub_xy[1]) # Min=0
            if self.idx_line == sv_idx:
                print('Edge of list')
        if event.key == '=':
            self.idx_line = min(len(self.llist['show_line'])-self.sub_xy[0]*self.sub_xy[1],
                                self.idx_line + self.sub_xy[0]*self.sub_xy[1]) 
            if self.idx_line == sv_idx:
                print('Edge of list')

        ## Reset z
        if event.key == 'z': 
            from astropy.relativity import velocities
            newz = velocities.z_from_v(self.z, event.xdata)
            self.z = newz
            self.abs_sys.zabs = newz
            # Drawing
            self.psdict['xmnx'] = self.vmnx

        # Single line command
        if event.key in ['1','2','B','U','L','N','V','A', 'x', 'X']:
            try:
                wrest = event.inaxes.get_gid()
            except AttributeError:
                return
            else:
                kwrest = wrest.value

        ## Velocity limits
        if event.key == '1': 
            self.abs_sys.lines[kwrest].analy['VLIM'][0] = event.xdata
        if event.key == '2': 
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            self.abs_sys.lines[kwrest].analy['VLIM'][1] = event.xdata
        if event.key == '!': 
            for key in self.abs_sys.lines.keys():
                try:
                    self.abs_sys.lines[key].analy['VLIM'][0] = event.xdata
                except KeyError:
                    print('Not setting VLIM for {:g}'.format(key))
        if event.key == '@': 
            for key in self.abs_sys.lines.keys():
                try:
                    self.abs_sys.lines[key].analy['VLIM'][1] = event.xdata
                except KeyError:
                    print('Not setting VLIM for {:g}'.format(key))
        ## Line type
        if event.key == 'A': # Add to lines
            if not kwrest in self.abs_sys.lines.keys():
                self.abs_sys.lines[kwrest] = xspec.analysis.Spectral_Line(wrest)
                print('VelPlot: Generating line {:g}'.format(kwrest))
                self.abs_sys.lines[kwrest].analy['VLIM'] = np.array([self.vmnx[0]/2.,
                                                                self.vmnx[1]/2.])
                self.abs_sys.lines[kwrest].analy['FLG_ANLY'] = 2 # Init to ok
                self.abs_sys.lines[kwrest].analy['DATFIL'] = self.spec_fil
        if event.key == 'x': # Remove line
            if kwrest in self.abs_sys.lines.keys():
                self.abs_sys.lines.pop(kwrest)
                print('VelPlot: Removed line {:g}'.format(wrest))
        if event.key == 'X': # Remove all lines (might add warning widget)
            # Double check
            gui = xguiu.WarningWidg('About to remove all lines. \n  Continue??')
            gui.exec_()
            if gui.ans is False:
                return
            #
            for kwrest in self.abs_sys.lines.keys():
                self.abs_sys.lines.pop(wrest)
                print('VelPlot: Removed line {:g}'.format(wrest))
        if event.key == 'B':  # Toggle blend
            try:
                feye = self.abs_sys.lines[kwrest].analy['FLG_EYE'] 
            except KeyError:
                feye = 0
            feye = (feye + 1) % 2
            self.abs_sys.lines[kwrest].analy['FLG_EYE']  = feye
        if event.key == 'N':  # Toggle NG
            try:
                fanly = self.abs_sys.lines[kwrest].analy['FLG_ANLY'] 
            except KeyError:
                fanly = 2
            if fanly == 0:
                fanly = 2 # Not using 1 anymore..
            else:
                fanly = 0
            self.abs_sys.lines[kwrest].analy['FLG_ANLY']  = fanly
        if event.key == 'V':  # Normal
            self.abs_sys.lines[kwrest].analy['FLG_LIMIT'] = 1
        if event.key == 'L':  # Lower limit
            self.abs_sys.lines[kwrest].analy['FLG_LIMIT'] = 2
        if event.key == 'U':  # Upper limit
            self.abs_sys.lines[kwrest].analy['FLG_LIMIT'] = 3
            
        # AODM plot
        if event.key == ':':  # 
            # Grab good lines
            from xastropy.xguis import spec_guis as xsgui
            gdl = []
            for iwr in self.abs_sys.lines.keys():
                if self.abs_sys.lines[iwr].analy['FLG_ANLY'] > 0:
                    gdl.append(iwr*u.AA)
            # Launch AODM
            gui = xsgui.XAODMGui(self.spec, self.z, gdl, vmnx=self.vmnx, norm=self.norm)
            gui.exec_()

            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()

        if not wrest is None: # Single window
            flg = 3
        if event.key in ['c','C','k','K','W','!', '@', '=', '-', 'X', 'z','R']: # Redraw all
            flg = 1 
        if event.key in ['Y']:
            rescale = False
        if event.key in ['k','c','C','K', 'R']:
            fig_clear = True

        if flg==1: # Default is not to redraw
            self.on_draw(rescale=rescale, fig_clear=fig_clear)
        elif flg==2: # Layer (no clear)
            self.on_draw(replot=False, rescale=rescale) 
        elif flg==3: # Layer (no clear)
            self.on_draw(in_wrest=wrest, rescale=rescale)

    # Click of main mouse button
    def on_click(self,event):
        try:
            print('button={:d}, x={:f}, y={:f}, xdata={:f}, ydata={:f}'.format(
                event.button, event.x, event.y, event.xdata, event.ydata))
        except ValueError:
            return
        if event.button == 1: # Draw line
            self.ax.plot( [event.xdata,event.xdata], self.psdict['ymnx'], ':', color='green')
            self.on_draw(replot=False) 
    
            # Print values
            try:
                self.statusBar().showMessage('x,y = {:f}, {:f}'.format(event.xdata,event.ydata))
            except AttributeError:
                return

    def on_draw(self, replot=True, in_wrest=None, rescale=True, fig_clear=False):
        """ Redraws the figure
        """
        #
        if replot is True:
            if fig_clear:
                self.fig.clf()
            # Loop on windows
            all_idx = self.llist['show_line']
            nplt = self.sub_xy[0]*self.sub_xy[1]
            if len(all_idx) <= nplt:
                self.idx_line = 0
            subp = np.arange(nplt) + 1
            subp_idx = np.hstack(subp.reshape(self.sub_xy[0],self.sub_xy[1]).T)
            for jj in range(min(nplt, len(all_idx))):
                try:
                    idx = all_idx[jj+self.idx_line]
                except IndexError:
                    continue # Likely too few lines
                # Grab line
                #wvobs = np.array((1+self.z) * self.llist[self.llist['List']]['wrest'][idx])
                wrest = self.llist[self.llist['List']].wrest[idx] # *
                        # self.llist[self.llist['List']].wrest.unit)
                kwrest = wrest.value # For the Dict
                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()
                # Single window?
                if not in_wrest is None:
                    if np.abs(wrest-in_wrest) > (1e-3*u.AA):
                        continue
                # Generate plot
                self.ax = self.fig.add_subplot(self.sub_xy[0],self.sub_xy[1], subp_idx[jj])
                self.ax.clear()        

                #print('Plotting {:g}, {:d}'.format(wrest,subp_idx[jj]))

                # Zero line
                self.ax.plot( [0., 0.], [-1e9, 1e9], ':', color='gray')
                # Velocity
                wvobs = (1+self.z) * wrest
                velo = (self.spec.dispersion/wvobs - 1.)*const.c.to('km/s').value
                
                # Plot
                self.ax.plot(velo, self.spec.flux, 'k-',drawstyle='steps-mid')

                # GID for referencing
                self.ax.set_gid(wrest)

                # Labels
                #if jj >= (self.sub_xy[0]-1)*(self.sub_xy[1]):
                if ((jj+1) % self.sub_xy[0]) == 0:
                    self.ax.set_xlabel('Relative Velocity (km/s)')
                else:
                    self.ax.get_xaxis().set_ticks([])
                #if ((jj+1) // 2 == 0) & (jj < self.sub_xy[0]):
                #    self.ax.set_ylabel('Relative Flux')
                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()
                lbl = self.llist[self.llist['List']].name[idx]
                self.ax.text(0.1, 0.05, lbl, color='blue', transform=self.ax.transAxes,
                             size='x-small', ha='left')

                # Reset window limits
                self.ax.set_xlim(self.psdict['xmnx'])

                # Rescale?
                if (rescale is True) & (self.norm is False):
                    gdp = np.where( (velo > self.psdict['xmnx'][0]) &
                                    (velo < self.psdict['xmnx'][1]))[0]
                    if len(gdp) > 5:
                        per = xstats.basic.perc(self.spec.flux[gdp])
                        self.ax.set_ylim((0., 1.1*per[1]))
                    else:
                        self.ax.set_ylim(self.psdict['ymnx'])
                else:
                    self.ax.set_ylim(self.psdict['ymnx'])

                # Fonts
                xputils.set_fontsize(self.ax,6.)

                # Abs_Sys: Color the lines
                if not self.abs_sys is None:
                    #QtCore.pyqtRemoveInputHook()
                    #xdb.set_trace()
                    #QtCore.pyqtRestoreInputHook()
                    try:
                        vlim = self.abs_sys.lines[kwrest].analy['VLIM']
                    except KeyError:
                        continue
                    # Color coding
                    clr = 'black'
                    try:  # .clm style
                        flag = self.abs_sys.lines[kwrest].analy['FLAGS'][0]
                    except KeyError:
                        flag = None
                    else:
                        if flag <= 1: # Standard detection
                            clr = 'green'
                        elif flag in [2,3]:
                            clr = 'blue'
                        elif flag in [4,5]:
                            clr = 'purple'
                    # ABS ID
                    try: # NG?
                        flagA = self.abs_sys.lines[kwrest].analy['FLG_ANLY']
                    except KeyError:
                        flagA = None
                    else:
                        if (flagA>0) & (clr == 'black'):
                            clr = 'green'
                    try: # Limit?
                        flagL = self.abs_sys.lines[kwrest].analy['FLG_LIMIT']
                    except KeyError:
                        flagL = None
                    else:
                        if flagL == 2:
                            clr = 'blue'
                        if flagL == 3:
                            clr = 'purple'
                    try: # Blends?
                        flagE = self.abs_sys.lines[kwrest].analy['FLG_EYE']
                    except KeyError:
                        flagE = None
                    else:
                        if flagE == 1:
                            clr = 'orange'
                    if flagA == 0:
                        clr = 'red'

                    pix = np.where( (velo > vlim[0]) & (velo < vlim[1]))[0]
                    self.ax.plot(velo[pix], self.spec.flux[pix], '-',
                                 drawstyle='steps-mid', color=clr)


        # Draw
        self.canvas.draw()
    
# ######################
class AODMWidget(QtGui.QWidget):
    ''' Widget for comparing tau_AODM profiles

        19-Dec-2014 by JXP
    '''
    def __init__(self, spec, z, wrest, parent=None, vmnx=[-300., 300.],
                 norm=True):
        '''
        spec = Spectrum1D
        '''
        super(AODMWidget, self).__init__(parent)

        # Initialize
        self.spec = spec
        self.norm = norm
        self.z = z
        self.vmnx = vmnx
        self.wrest = wrest  # Expecting (requires) units
        self.lines = []
        for iwrest in self.wrest:
            self.lines.append(xspec.analysis.Spectral_Line(iwrest))


        self.psdict = {} # Dict for spectra plotting
        self.psdict['xmnx'] = self.vmnx
        self.psdict['ymnx'] = [-0.1, 1.1]
        self.psdict['nav'] = navigate(0,0,init=True)

        # Create the mpl Figure and FigCanvas objects. 
        #
        self.dpi = 150
        self.fig = Figure((8.0, 4.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)

        self.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
        self.canvas.setFocus()
        self.canvas.mpl_connect('key_press_event', self.on_key)
        self.canvas.mpl_connect('button_press_event', self.on_click)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.canvas)
        
        self.setLayout(vbox)

        # Draw on init
        self.on_draw()

    # Key stroke 
    def on_key(self,event):

        # Init
        rescale = True
        flg = 0

        ## NAVIGATING
        if event.key in self.psdict['nav']: 
            flg = navigate(self.psdict,event)
        if event.key in ['b','t','W','Z','Y','l','r']:  
            rescale = False

        self.on_draw(rescale=rescale)

    # Click of main mouse button
    def on_click(self,event):
        return # DO NOTHING FOR NOW
        try:
            print('button={:d}, x={:f}, y={:f}, xdata={:f}, ydata={:f}'.format(
                event.button, event.x, event.y, event.xdata, event.ydata))
        except ValueError:
            return
        if event.button == 1: # Draw line
            self.ax.plot( [event.xdata,event.xdata], self.psdict['ymnx'], ':', color='green')
            self.on_draw()
    
            # Print values
            try:
                self.statusBar().showMessage('x,y = {:f}, {:f}'.format(event.xdata,event.ydata))
            except AttributeError:
                return

    def on_draw(self, rescale=True):
        """ Redraws the figure
        """
        #
        self.ax = self.fig.add_subplot(1,1,1)
        self.ax.clear()

        ymx = 0.
        for ii,iwrest in enumerate(self.wrest):

            # Velocity
            wvobs = (1+self.z) * iwrest
            velo = (self.spec.dispersion/wvobs - 1.)*const.c.to('km/s').value
            gdp = np.where((velo > self.psdict['xmnx'][0]) &
                           (velo < self.psdict['xmnx'][1]))[0]

            # Normalize?
            if self.norm is False:
                per = xstats.basic.perc(self.spec.flux[gdp])
                fsplice = per[1] / self.spec.flux[gdp] 
            else:
                fsplice = 1./ self.spec.flux[gdp]

            # AODM
            cst = (10.**14.5761)/(self.lines[ii].atomic['fval']*iwrest.value)
            Naodm = np.log(fsplice)*cst
            ymx = max(ymx,np.max(Naodm))
                
            # Plot
            line, = self.ax.plot(velo[gdp], Naodm, '-', drawstyle='steps-mid')

            # Labels
            lbl = '{:g}'.format(iwrest)
            clr = plt.getp(line, 'color') 
            self.ax.text(0.1, 1.-(0.05+0.05*ii), lbl, color=clr,
                         transform=self.ax.transAxes, size='small', ha='left')

        self.ax.set_xlabel('Relative Velocity (km/s)')
        self.ax.set_ylabel('N(AODM)')
        # Zero line
        self.ax.plot( [0., 0.], [-1e29, 1e29], ':', color='gray')

        # Reset window limits
        self.ax.set_xlim(self.psdict['xmnx'])
        if rescale:
            self.psdict['ymnx'] = [0.05*ymx, ymx*1.1]
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        self.ax.set_ylim(self.psdict['ymnx'])

        # Fonts
        #xputils.set_fontsize(self.ax,6.)

        # Draw
        self.canvas.draw()
    


# ######
# Plot Doublet
def set_doublet(iself,event):
    ''' Set z and plot doublet
    '''
    wv_dict = {'C': (1548.195, 1550.770, 'CIV'), 'M': (2796.352, 2803.531, 'MgII'),
               '4': (1393.755, 1402.770, 'SiIV'),
               'X': (1031.9261, 1037.6167, 'OVI'), '8': (770.409, 780.324, 'NeVIII'),
               'B': (1025.4433, 1215.6701, 'Lyba')}
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
def navigate(psdict,event,init=False):
    ''' Method to Navigate spectrum
    init:  (False) Initialize
      Just pass back valid key strokes
    '''
    # Initalize
    if init is True:
        return ['l','r','b','t','T','i','I', 'o','O', '[',']','W','Z', 'Y', '{', '}']

    #
    if (not isinstance(event.xdata,float)) or (not isinstance(event.ydata,float)):
        print('Navigate: You entered the {:s} key out of bounds'.format(event.key))
        return 0

    if event.key == 'l':  # Set left
        psdict['xmnx'][0] = event.xdata
    elif event.key == 'r':  # Set Right
        psdict['xmnx'][1] = event.xdata
    elif event.key == 'b':  # Set Bottom
        psdict['ymnx'][0] = event.ydata
    elif event.key == 't':  # Set Top
        psdict['ymnx'][1] = event.ydata
    elif event.key == 'T':  # Set Top to 1.1
        psdict['ymnx'][1] = 1.1
    elif event.key == 'i':  # Zoom in (and center)
        deltx = (psdict['xmnx'][1]-psdict['xmnx'][0])/4.
        psdict['xmnx'] = [event.xdata-deltx, event.xdata+deltx]
    elif event.key == 'I':  # Zoom in (and center)
        deltx = (psdict['xmnx'][1]-psdict['xmnx'][0])/16.
        psdict['xmnx'] = [event.xdata-deltx, event.xdata+deltx]
    elif event.key == 'o':  # Zoom in (and center)
        deltx = psdict['xmnx'][1]-psdict['xmnx'][0]
        psdict['xmnx'] = [event.xdata-deltx, event.xdata+deltx]
    elif event.key == 'O':  # Zoom in (and center)
        deltx = psdict['xmnx'][1]-psdict['xmnx'][0]
        psdict['xmnx'] = [event.xdata-2*deltx, event.xdata+2*deltx]
    elif event.key == 'Y':  # Zoom in (and center)
        delty = psdict['ymnx'][1]-psdict['ymnx'][0]
        psdict['ymnx'] = [event.ydata-delty, event.ydata+delty]
    elif event.key in ['[',']','{','}']:  # Pan 
        center = (psdict['xmnx'][1]+psdict['xmnx'][0])/2.
        deltx = (psdict['xmnx'][1]-psdict['xmnx'][0])/2.
        if event.key == '[':
            new_center = center - deltx
        elif event.key == ']':
            new_center = center + deltx
        elif event.key == '{':
            new_center = center - 4*deltx
        elif event.key == '}':
            new_center = center + 4*deltx
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


# ######
# 
def set_llist(llist,in_dict=None):
    ''' Method to set a line list dict for the Widgets
    '''
    from linetools.lists.linelist import LineList
    if in_dict is None:
        in_dict = {}

    if type(llist) in [str,unicode]: # Set line list from a file
        in_dict['List'] = llist
        if llist == 'None':
            in_dict['Plot'] = False
        else:
            in_dict['Plot'] = True
            # Load?
            if not (llist in in_dict):
                #line_file = xa_path+'/data/spec_lines/'+llist
                #llist_cls = xspec.abs_line.Abs_Line_List(llist)
                llist_cls = LineList(llist)
                in_dict[llist] = llist_cls
    elif isinstance(llist,list): # Set from a list of wrest

        from astropy.table import Column

        in_dict['List'] = 'input.lst'
        in_dict['Plot'] = True
        # Fill
        llist.sort()
        tmp_dict = {}
        # Parse from grb.lst
        llist_cls = LineList('ISM', gd_lines=llist)
        in_dict['input.lst'] = llist_cls
        '''
        line_file = xa_path+'/data/spec_lines/grb.lst'
        llist_cls = xspec.abs_line.Abs_Line_List(line_file)
        adict = llist_cls.data
        # Fill 
        names = []
        fval = []
        for wrest in llist:
            mt = np.where(np.abs(wrest-adict['wrest']) < 1e-3)[0]
            if len(mt) != 1:
                raise ValueError('Problem!')
            names.append(adict['name'][mt][0])
            fval.append(adict['fval'][mt][0])
        # Set
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        # Generate a Table
        col0 = Column(np.array(llist), name='wrest', unit=u.AA) # Assumed Angstroms
        col1 = Column(np.array(names), name='name')
        col2 = Column(np.array(fval), name='fval')
        in_dict['input.lst'] = Table( (col0,col1,col2) )
        '''

    # Return
    return in_dict

# Read spectrum, pass back it and spec_file name
def read_spec(ispec, second_file=None):
    #
    if isinstance(ispec,str) or isinstance(ispec,unicode):
        spec_fil = ispec
        spec = lsi.readspec(spec_fil)
        # Second file?
        if not second_file is None:
            spec2 = lsi.readspec(second_file)
            if spec2.sig is None:
                spec2.sig = np.zeros(spec.flux.size)
            # Scale for convenience of plotting
            xper1 = xstats.basic.perc(spec.flux, per=0.9)
            xper2 = xstats.basic.perc(spec2.flux, per=0.9)
            scl = xper1[1]/xper2[1]
            # Stitch together
            wave3 = np.append(spec.dispersion, spec2.dispersion)
            flux3 = np.append(spec.flux, spec2.flux*scl)
            sig3 = np.append(spec.sig, spec2.sig*scl)
            spec3 = Spectrum1D.from_array(wave3, flux3, uncertainty=StdDevUncertainty(sig3))
            # Overwrite
            spec = spec3
            spec.filename = spec_fil
    else:
        spec = ispec # Assuming Spectrum1D
        spec_fil = spec.filename # Grab from Spectrum1D 
    # Return
    return spec, spec_fil

# ################
# TESTING
if __name__ == "__main__":
    from xastropy import spec as xspec

    if len(sys.argv) == 1: #
        flg_tst = 0
        flg_tst += 2**0  # ExamineSpecWidget
        #flg_tst += 2**1  # PlotLinesWidget
        #flg_tst += 2**2  # SelectLineWidget
        #flg_tst += 2**3  # AbsSysWidget
        #flg_tst += 2**4  # VelPltWidget
        #flg_tst += 2**5  # SelectedLinesWidget
        #flg_tst += 2**6  # AODMWidget
    else:
        flg_tst = int(sys.argv[1])

    # ExamineSpec
    if (flg_tst % 2) == 1:
        app = QtGui.QApplication(sys.argv)
        spec_fil = '/u/xavier/Keck/HIRES/RedData/PH957/PH957_f.fits'
        spec = lsi.readspec(spec_fil)
        app.setApplicationName('XSpec')
        main = ExamineSpecWidget(spec)
        main.show()
        sys.exit(app.exec_())

    # PltLineWidget
    if (flg_tst % 2**2) >= 2**1:
        app = QtGui.QApplication(sys.argv)
        app.setApplicationName('PltLine')
        main = PlotLinesWidget()
        main.show()
        sys.exit(app.exec_())

    # SelectLineWidget
    if (flg_tst % 2**3) >= 2**2:
        orig = False
        llist_cls = LineList('ISM')

        app = QtGui.QApplication(sys.argv)
        app.setApplicationName('SelectLine')
        main = SelectLineWidget(llist_cls._data)
        main.show()
        app.exec_()
        print(main.line)
        # Another test
        quant = main.line.split('::')[1].lstrip()
        spltw = quant.split(' ')
        wrest = Quantity(float(spltw[0]), unit=spltw[1])
        print(wrest)
        sys.exit()
    
    # AbsSys Widget
    if (flg_tst % 2**4) >= 2**3:
        abs_fil = '/Users/xavier/paper/LLS/Optical/Data/Analysis/MAGE/SDSSJ1004+0018_z2.746_id.fits'
        abs_fil2 = '/Users/xavier/paper/LLS/Optical/Data/Analysis/MAGE/SDSSJ2319-1040_z2.675_id.fits'
        app = QtGui.QApplication(sys.argv)
        app.setApplicationName('AbsSys')
        main = AbsSysWidget([abs_fil,abs_fil2])
        main.show()
        sys.exit(app.exec_())

    # VelPlt Widget
    if (flg_tst % 2**5) >= 2**4:
        specf = 1
        if specf == 0: # PH957 DLA
            # Spectrum
            spec_fil = '/u/xavier/Keck/HIRES/RedData/PH957/PH957_f.fits'
            spec = lsi.readspec(spec_fil)
            # Abs_sys
            abs_sys = xiaa.Generic_System(None)
            abs_sys.clm_fil = '/Users/xavier/DLA/Abund/PH957.z2309.clm'
            abs_sys.get_ions(skip_ions=True, fill_lines=True)
            abs_sys.zabs = abs_sys.clm_analy.zsys
        elif specf == 1: # UM184 LLS
            # Spectrum
            spec_fil = '/Users/xavier/PROGETTI/LLSZ3/data/normalize/UM184_nF.fits'
            spec = lsi.readspec(spec_fil)
            # Abs_sys
            abs_fil = '/Users/xavier/paper/LLS/Optical/Data/Analysis/MAGE/UM184_z2.930_id.fits'
            abs_sys = xiaa.Generic_System(None)
            abs_sys.parse_absid_file(abs_fil)
        # Launch
        app = QtGui.QApplication(sys.argv)
        app.setApplicationName('VelPlot')
        main = VelPlotWidget(spec, abs_sys=abs_sys)
        main.show()
        sys.exit(app.exec_())

    # SelectedLines Widget
    if (flg_tst % 2**6) >= 2**5:
        print('Test: SelectedLines Widget')
        llist = set_llist('ISM')
        # Launch
        app = QtGui.QApplication(sys.argv)
        app.setApplicationName('SelectedLines')
        main = SelectedLinesWidget(llist['ISM'])#._data)
        main.show()
        sys.exit(app.exec_())

    # AODM Widget
    if (flg_tst % 2**7) >= 2**6:
        spec_fil = '/Users/xavier/PROGETTI/LLSZ3/data/normalize/UM184_nF.fits'
        spec = lsi.readspec(spec_fil)
        z=2.96916
        lines = np.array([1548.195, 1550.770]) * u.AA
        # Launch
        app = QtGui.QApplication(sys.argv)
        app.setApplicationName('AODM')
        main = AODMWidget(spec, z, lines)
        main.show()
        sys.exit(app.exec_())
