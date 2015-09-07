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
import os, sys, warnings, imp
import matplotlib.pyplot as plt
import glob

from PyQt4 import QtGui
from PyQt4 import QtCore

from matplotlib import mpl
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
# Matplotlib Figure object
from matplotlib.figure import Figure

from astropy.units import Quantity
from astropy.nddata import StdDevUncertainty
from astropy import units as u
from astropy.io import fits, ascii

from linetools.lists.linelist import LineList
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectra import convolve as lsc
from linetools.spectralline import AbsLine

from xastropy.xutils import xdebug as xdb
from xastropy.xguis import spec_widgets as xspw
from xastropy.xguis import utils as xxgu
from xastropy.igm.abs_sys.lls_utils import LLSSystem
from xastropy.igm.abs_sys import lls_utils as xialu
from xastropy.atomic import ionization as xatomi
from xastropy.spec import voigt as xsv
from xastropy.spec import continuum as xspc

xa_path = imp.find_module('xastropy')[1]

#class XSpecGui(QtGui.QMainWindow):
#class XAbsIDGui(QtGui.QMainWindow):
#class XVelPltGui(QtGui.QDialog):
# class XFitLLSGUI(QtGui.QMainWindow):

# x_specplot replacement
class XSpecGui(QtGui.QMainWindow):
    ''' GUI to replace XIDL x_specplot

        12-Dec-2014 by JXP v1.0
        27-Mar-2015 by JXP v2.0 :: EW, column, better zooming + panning
    '''
    def __init__(self, ispec, parent=None, zsys=None, norm=None):
        QtGui.QMainWindow.__init__(self, parent)
        '''
        Parameters
        --------------
        spec = Spectrum1D or tuple of arrays
          Input spectrum.  If tuple then (wave,fx) or (wave,fx,sig)
        '''
        from linetools.lists.linelist import LineList
        #reload(xxgu)
        #reload(xspw)
        # INIT
        spec,_ = xxgu.read_spec(ispec)
        #xdb.set_trace()

        # 
        mpl.rcParams['agg.path.chunksize'] = 20000 # Needed to avoid carsh in large spectral files
        
        # Build a widget combining several others
        self.main_widget = QtGui.QWidget()

        # Status bar
        self.create_status_bar()

        # Grab the pieces and tie together
        self.pltline_widg = xspw.PlotLinesWidget(status=self.statusBar, init_z=zsys)
        self.pltline_widg.setMaximumWidth(300)

        # Hook the spec widget to Plot Line
        self.spec_widg = xspw.ExamineSpecWidget(spec,status=self.statusBar,
                                                llist=self.pltline_widg.llist,
                                                zsys=zsys, norm=norm)
        self.pltline_widg.spec_widg = self.spec_widg

        self.spec_widg.canvas.mpl_connect('button_press_event', self.on_click)

        extras = QtGui.QWidget()
        extras.setMaximumWidth(130)
        vbox = QtGui.QVBoxLayout()
        qbtn = QtGui.QPushButton('Quit', self)
        qbtn.clicked.connect(self.quit)
        vbox.addWidget(self.pltline_widg)
        vbox.addWidget(qbtn)
        extras.setLayout(vbox)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(self.spec_widg)
        hbox.addWidget(extras)

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
                self.pltline_widg.llist[self.pltline_widg.llist['List']]._data)
            self.select_line_widg.exec_()
            line = self.select_line_widg.line
            if line.strip() == 'None':
                return
            #
            quant = line.split('::')[1].lstrip()
            spltw = quant.split(' ')
            wrest = Quantity(float(spltw[0]), unit=spltw[1])
            z = event.xdata/wrest.value - 1.
            self.pltline_widg.llist['z'] = z
            self.statusBar().showMessage('z = {:f}'.format(z))

            self.pltline_widg.zbox.setText('{:.5f}'.format(self.pltline_widg.llist['z']))
    
            # Draw
            self.spec_widg.on_draw()
    # Quit
    def quit(self):
        self.close()

# GUI for Identifying many (all) Abs Systems in a Spectrum
class XAbsIDGui(QtGui.QMainWindow):
    ''' GUI to analyze absorption systems in a spectrum

        16-Dec-2014 by JXP
    '''
    def __init__(self, inspec, parent=None, abssys_dir=None, absid_list=None, norm=True,
                 srch_id=True, id_dir='ID_LINES/'):
        QtGui.QMainWindow.__init__(self, parent)
        '''
        spec = Spectrum1D
        second_file = Second spectrum file  (e.g. COS + STIS)
        '''
        # Build a widget combining several others
        self.main_widget = QtGui.QWidget()

        # Status bar
        self.create_status_bar()

        # Initialize
        if absid_list is None:
            # Automatically search for ID files
            if srch_id:
                absid_list = glob.glob(id_dir+'*id.fits')
            else:
                absid_list = []

        # Grab the pieces and tie together
        self.abssys_widg = xspw.AbsSysWidget(absid_list)
        self.pltline_widg = xspw.PlotLinesWidget(status=self.statusBar)
        self.spec_widg = xspw.ExamineSpecWidget(inspec,status=self.statusBar,
                                                llist=self.pltline_widg.llist, norm=norm,
                                                abs_sys=self.abssys_widg.abs_sys)
        self.pltline_widg.spec_widg = self.spec_widg

        # Connections
        self.spec_widg.canvas.mpl_connect('button_press_event', self.on_click)
        self.spec_widg.canvas.mpl_connect('key_press_event', self.on_key)
        self.abssys_widg.refine_button.clicked.connect(self.refine_abssys) 

        # Layout
        anly_widg = QtGui.QWidget()
        anly_widg.setMaximumWidth(300)
        anly_widg.setMinimumWidth(150)
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.pltline_widg)
        vbox.addWidget(self.abssys_widg)
        anly_widg.setLayout(vbox)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(self.spec_widg)
        hbox.addWidget(anly_widg)

        self.main_widget.setLayout(hbox)

        # Point MainWindow
        self.setCentralWidget(self.main_widget)

    def create_status_bar(self):
        self.status_text = QtGui.QLabel("XAbsID")
        self.statusBar().addWidget(self.status_text, 1)

    def on_key(self,event):
        if event.key == 'v': # Stack plot
            if self.spec_widg.vplt_flg == 1:
                self.abssys_widg.add_fil(self.spec_widg.outfil)
                self.abssys_widg.reload()
            elif self.spec_widg.vplt_flg == -1:
                return
            # Update line list
            idx = self.pltline_widg.lists.index(self.spec_widg.llist['List'])
            self.pltline_widg.llist_widget.setCurrentRow(idx)
        elif event.key == '?': # Check for a match with known systems
            wv_chosen = event.xdata
            # Load grb 
            llist = xspw.set_llist('grb.lst')

            # Loop through systems
            for iabs_sys in self.abssys_widg.all_abssys:
                z = iabs_sys.zabs
                wvobs = np.array((1+z) * llist['grb.lst']['wrest'])
                mtwv = np.where( np.abs( wvobs-wv_chosen ) < 0.2 )[0]
                for imt in mtwv:
                    #QtCore.pyqtRemoveInputHook()
                    #xdb.set_trace()
                    #QtCore.pyqtRestoreInputHook()
                    print('z={:g},  {:s},  f={:g}'.format(z, 
                        llist['grb.lst']['name'][imt],
                        llist['grb.lst']['fval'][imt]))
                if len(mtwv) == 0:
                    print('No match. wrest={:g} for z={:g}'.format(wv_chosen/(1+z), z))

    def on_click(self,event):
        if event.button == 3: # Set redshift
            # Line list?
            try:
                self.pltline_widg.llist['List']
            except KeyError:
                print('Set a line list first!!')
                return
            # 
            if self.pltline_widg.llist[self.pltline_widg.llist['List']] == 'None':
                return
            self.select_line_widg = xspw.SelectLineWidget(
                self.pltline_widg.llist[self.pltline_widg.llist['List']]._data)
            self.select_line_widg.exec_()
            line = self.select_line_widg.line
            if line.strip() == 'None':
                return
            #
            quant = line.split('::')[1].lstrip()
            spltw = quant.split(' ')
            wrest = Quantity(float(spltw[0]), unit=spltw[1])
            z = event.xdata/wrest.value - 1.
            self.pltline_widg.llist['z'] = z
            self.statusBar().showMessage('z = {:f}'.format(z))

            self.pltline_widg.zbox.setText(self.pltline_widg.zbox.z_frmt.format(
                self.pltline_widg.llist['z']))
    
            # Draw
            self.spec_widg.on_draw()

    def refine_abssys(self):
        item = self.abssys_widg.abslist_widget.selectedItems()
        if len(item) != 1:
            self.statusBar().showMessage('AbsSys: Must select only 1 system!')
            print('AbsSys: Must select only 1 system!')
        txt = item[0].text()
        ii = self.abssys_widg.all_items.index(txt)
        iabs_sys = self.abssys_widg.all_abssys[ii]
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        # Launch
        gui = XVelPltGui(self.spec_widg.spec, outfil=iabs_sys.absid_file,
                               abs_sys=iabs_sys, norm=self.spec_widg.norm)
        gui.exec_()

# ##################################
# GUI for velocity plot
class XVelPltGui(QtGui.QDialog):
    ''' GUI to analyze absorption systems in a spectrum
        24-Dec-2014 by JXP
    '''
    def __init__(self, ispec, z=None, parent=None, llist=None, norm=True,
                 vmnx=[-300., 300.]*u.km/u.s, abs_sys=None, outfil='dum_ID.fits',
                 sel_wv=None):
        '''
        spec = Filename or Spectrum1D
        Norm: Bool (False)
          Normalized spectrum?
        abs_sys: AbsSystem
          Absorption system class
        sel_wv: Selected wavelength.  Used to inspect a single, unknown line
        '''
        super(XVelPltGui, self).__init__(parent)

        # Initialize
        self.abs_sys = abs_sys
        if self.abs_sys is not None:
            self.z = self.abs_sys.zabs
        else:
            if z is None:
                raise ValueError('XVelPlt: Need to set abs_sys or z!')
            self.z = z
        self.vmnx = vmnx
        self.outfil = outfil
        self.norm = norm
        self.sel_wv = sel_wv

        # Grab the pieces and tie together
        self.vplt_widg = xspw.VelPlotWidget(ispec, abs_sys=self.abs_sys, llist=llist, 
                                            vmnx=self.vmnx, z=self.z, norm=self.norm)
        self.pltline_widg = xspw.PlotLinesWidget(init_llist=self.vplt_widg.llist,
                                                 init_z=self.z)
        #self.pltline_widg.spec_widg = self.vplt_widg

        self.slines = xspw.SelectedLinesWidget(self.vplt_widg.llist[self.vplt_widg.llist['List']],
                                               init_select=self.vplt_widg.llist['show_line'],
                                               plot_widget=self.vplt_widg)

        # Connections
        self.pltline_widg.llist_widget.currentItemChanged.connect(self.on_llist_change)
        self.connect(self.pltline_widg.zbox, QtCore.SIGNAL('editingFinished ()'), self.setz)
        self.vplt_widg.canvas.mpl_connect('key_press_event', self.on_key)

        # Outfil
        wbtn = QtGui.QPushButton('Write', self)
        wbtn.setAutoDefault(False)
        wbtn.clicked.connect(self.write_out)
        self.out_box = QtGui.QLineEdit()
        self.out_box.setText(self.outfil)
        self.connect(self.out_box, QtCore.SIGNAL('editingFinished ()'), self.set_outfil)

        # Quit
        buttons = QtGui.QWidget()
        wqbtn = QtGui.QPushButton('Write+Quit', self)
        wqbtn.setAutoDefault(False)
        wqbtn.clicked.connect(self.write_quit)
        qbtn = QtGui.QPushButton('Quit', self)
        qbtn.setAutoDefault(False)
        qbtn.clicked.connect(self.quit)

        # Sizes
        lines_widg = QtGui.QWidget()
        lines_widg.setMaximumWidth(300)
        lines_widg.setMinimumWidth(200)

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.pltline_widg)
        vbox.addWidget(self.slines)
        vbox.addWidget(wbtn)
        vbox.addWidget(self.out_box)
        # Write/Quit buttons
        hbox1 = QtGui.QHBoxLayout()
        hbox1.addWidget(wqbtn)
        hbox1.addWidget(qbtn)
        buttons.setLayout(hbox1)
        #
        vbox.addWidget(buttons)
        lines_widg.setLayout(vbox)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(self.vplt_widg)
        hbox.addWidget(lines_widg)

        self.setLayout(hbox)
        # Initial draw
        self.vplt_widg.on_draw()

    # Change z
    def on_key(self,event):
        if event.key == 'z': 
            self.z = self.vplt_widg.z
            self.pltline_widg.llist['z'] = self.z
            self.pltline_widg.zbox.setText(self.pltline_widg.zbox.z_frmt.format(self.z))
        if event.key == 'T':  # Try another rest wavelength for input line
            # Get line from User
            self.select_line_widg = xspw.SelectLineWidget(
                self.pltline_widg.llist[self.pltline_widg.llist['List']]._data)
            self.select_line_widg.exec_()
            line = self.select_line_widg.line
            quant = line.split('::')[1].lstrip()
            spltw = quant.split(' ')
            wrest = Quantity(float(spltw[0]), unit=spltw[1])
            # Set redshift
            self.z = self.sel_wv / wrest - 1.
            print('Setting z = {:g}'.format(self.z))
            self.pltline_widg.llist['z'] = self.z
            self.pltline_widg.zbox.setText(self.pltline_widg.zbox.z_frmt.format(self.z))
            self.vplt_widg.z = self.pltline_widg.llist['z']
            # Reset
            self.vplt_widg.init_lines()
            self.vplt_widg.on_draw()

    # Set z from pltline_widg
    def setz(self):
        self.vplt_widg.abs_sys.zabs = self.pltline_widg.llist['z']
        self.vplt_widg.z = self.pltline_widg.llist['z']
        self.z = self.pltline_widg.llist['z']
        self.vplt_widg.on_draw()

    # Change list of lines to choose from
    def on_llist_change(self):
        llist = self.pltline_widg.llist
        all_lines = list( llist[llist['List']]._data['wrest'] )
        # Set selected
        abs_sys = self.vplt_widg.abs_sys
        wrest = [line.wrest for line in abs_sys.lines]
        select = []
        for iwrest in wrest:
            try:
                select.append(all_lines.index(iwrest))
            except ValueError:
                pass
        select.sort()
        # GUIs
        self.vplt_widg.llist['List'] = llist['List']
        self.vplt_widg.llist['show_line'] = select
        self.vplt_widg.idx_line = 0
        self.slines.selected = select
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        self.slines.on_list_change(llist[llist['List']])

    # Write
    def set_outfil(self):
        self.outfil = str(self.out_box.text())
        print('XVelPlot: Will write to {:s}'.format(self.outfil))

    # Write
    def write_out(self):
        self.vplt_widg.abs_sys.absid_file = self.outfil
        self.vplt_widg.abs_sys.write_absid_file()

    # Write + Quit
    def write_quit(self):
        self.write_out()
        self.flg_quit = 1
        self.abs_sys = self.vplt_widg.abs_sys
        self.done(1)

    # Write + Quit
    def quit(self):
        #self.abs_sys = self.vplt_widg.abs_sys # Have to write to pass back
        self.flg_quit = 0
        self.done(1)

# x_specplot replacement
class XAODMGui(QtGui.QDialog):
    ''' GUI to show AODM plots
        28-Dec-2014 by JXP
    '''
    def __init__(self, spec, z, wrest, vmnx=[-300., 300.]*u.km/u.s, parent=None, norm=True):
        super(XAODMGui, self).__init__(parent)
        '''
        spec = Spectrum1D
        '''
        # Grab the pieces and tie together
        self.aodm_widg = xspw.AODMWidget(spec,z,wrest,vmnx=vmnx,norm=norm)
        self.aodm_widg.canvas.mpl_connect('key_press_event', self.on_key)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.aodm_widg)

        self.setLayout(vbox)
        self.aodm_widg.on_draw()

    def on_key(self,event):
        if event.key == 'q': # Quit
            self.done(1)

# GUI for fitting LLS in a spectrum
class XFitLLSGUI(QtGui.QMainWindow):
    ''' GUI to fit LLS in a given spectrum
        v0.6
        30-Jul-2015 by JXP
    '''
    def __init__(self, ispec, parent=None, lls_fit_file=None, 
        outfil=None, smooth=3., zqso=None, fN_gamma=None): 
        QtGui.QMainWindow.__init__(self, parent)
        '''
        ispec = Spectrum1D or specfil
        lls_fit_file: str, optional
          Name of the LLS fit file to input
        smooth: float, optional
          Number of pixels to smooth on (FWHM)
        zqso: float, optional
          Redshift of the quasar.  If input, a Telfer continuum is used
        fN_gamma: float, optional
          Redshift evolution of f(N) or IGM fiddled continuum
        '''

        # Build a widget combining several others
        self.main_widget = QtGui.QWidget()

        # Status bar
        self.create_status_bar()

        # Initialize
        if outfil is None:
            self.outfil = 'LLS_fit.json'
        else:
            self.outfil = outfil
        self.count_lls = 0
        self.lls_model = None
        self.all_forest = []
        self.flag_write = False

        # Spectrum
        spec, spec_fil = xxgu.read_spec(ispec)

        # Continuum
        self.conti_dict = xspc.init_conti_dict(
            Norm=float(np.median(spec.flux.value)),
            piv_wv=np.median(spec.dispersion.value),
            igm='True')
        if zqso is not None:
            self.zqso = zqso
            # Read Telfer and apply IGM
            tspec = xspc.get_telfer_spec(zqso=zqso,
                igm=(self.conti_dict['igm']=='True'))
            # Rebin
            self.continuum = tspec.rebin(spec.dispersion)
            # Reset pivot wave
            self.conti_dict['piv_wv'] = 915.*(1+zqso)
        else:
            self.zqso = None
            self.continuum = XSpectrum1D.from_tuple((
                spec.dispersion,np.ones(len(spec.dispersion))))
        self.base_continuum = self.continuum.flux
        self.update_conti()

        # Full Model (LLS+continuum)
        self.full_model = XSpectrum1D.from_tuple((
            spec.dispersion,np.ones(len(spec.dispersion))))
        self.smooth = smooth

        # LineList
        self.llist = xxgu.set_llist('Strong') 
        self.llist['z'] = 0.
        self.plt_wv = zip(np.array([911.7,972.5367,1025.7222,1215.6700])*u.AA,
            ['LL','Lyg','Lyb','Lya'])

        # z and N boxes
        self.zwidget = xxgu.EditBox(-1., 'z_LLS=', '{:0.5f}')
        self.Nwidget = xxgu.EditBox(-1., 'NHI=', '{:0.2f}')
        self.bwidget = xxgu.EditBox(-1., 'b=', '{:0.1f}')
        self.Cwidget = xxgu.EditBox('None', 'Comment=', '{:s}')

        # Grab the pieces and tie together
        self.abssys_widg = xspw.AbsSysWidget([],only_one=True,
            no_buttons=True, linelist=self.llist[self.llist['List']])
        self.spec_widg = xspw.ExamineSpecWidget(spec,status=self.statusBar,
                                                llist=self.llist, key_events=False,
                                                abs_sys=self.abssys_widg.abs_sys)
        self.spec_widg.continuum = self.continuum

        #if other_spec is not None:
        #    ospec, ospec_fil = xxgu.read_spec(other_spec)
        #    self.spec_widg.other_spec = ospec

        # Initial file
        if lls_fit_file is not None:
            self.init_LLS(lls_fit_file)

        # Outfil
        wbtn = QtGui.QPushButton('Write', self)
        wbtn.setAutoDefault(False)
        wbtn.clicked.connect(self.write_out)
        #self.out_box = QtGui.QLineEdit()
        #self.out_box.setText(self.outfil)
        #self.connect(self.out_box, QtCore.SIGNAL('editingFinished ()'), self.set_outfil)

        # Quit
        buttons = QtGui.QWidget()
        wqbtn = QtGui.QPushButton('Write\n Quit', self)
        wqbtn.setAutoDefault(False)
        wqbtn.clicked.connect(self.write_quit)
        qbtn = QtGui.QPushButton('Quit', self)
        qbtn.setAutoDefault(False)
        qbtn.clicked.connect(self.quit)

        # Connections (buttons are above)
        self.spec_widg.canvas.mpl_connect('key_press_event', self.on_key)
        self.abssys_widg.abslist_widget.itemSelectionChanged.connect(
            self.on_list_change)
        self.connect(self.Nwidget.box, 
            QtCore.SIGNAL('editingFinished ()'), self.setbzN)
        self.connect(self.zwidget.box, 
            QtCore.SIGNAL('editingFinished ()'), self.setbzN)
        self.connect(self.bwidget.box, 
            QtCore.SIGNAL('editingFinished ()'), self.setbzN)
        self.connect(self.Cwidget.box, 
            QtCore.SIGNAL('editingFinished ()'), self.setbzN)

        # Layout
        anly_widg = QtGui.QWidget()
        anly_widg.setMaximumWidth(400)
        anly_widg.setMinimumWidth(250)

        # Write/Quit buttons
        hbox1 = QtGui.QHBoxLayout()
        hbox1.addWidget(wbtn)
        hbox1.addWidget(wqbtn)
        hbox1.addWidget(qbtn)
        buttons.setLayout(hbox1)
 
        # z,N
        zNwidg = QtGui.QWidget()
        hbox2 = QtGui.QHBoxLayout()
        hbox2.addWidget(self.zwidget)
        hbox2.addWidget(self.Nwidget)
        hbox2.addWidget(self.bwidget)
        zNwidg.setLayout(hbox2)
        #vbox.addWidget(self.pltline_widg)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(zNwidg)
        vbox.addWidget(self.Cwidget)
        vbox.addWidget(self.abssys_widg)
        vbox.addWidget(buttons)
        anly_widg.setLayout(vbox)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(self.spec_widg)
        hbox.addWidget(anly_widg)

        self.main_widget.setLayout(hbox)

        # Point MainWindow
        self.setCentralWidget(self.main_widget)

    def on_list_change(self):
        self.update_boxes()

    def create_status_bar(self):
        self.status_text = QtGui.QLabel("XFitLLS v0.4.2")
        self.statusBar().addWidget(self.status_text, 1)

    def setbzN(self):
        '''Set the column density or redshift from the box
        '''
        idx = self.get_sngl_sel_sys()
        if idx is None:
            return
        self.abssys_widg.all_abssys[idx].NHI = (
            float(self.Nwidget.box.text()))
        self.abssys_widg.all_abssys[idx].zabs = (
            float(self.zwidget.box.text()))
        self.abssys_widg.all_abssys[idx].bval = (
            float(self.bwidget.box.text()))*u.km/u.s
        self.abssys_widg.all_abssys[idx].comment = (
            self.Cwidget.box.text())
        # Update the lines
        for iline in self.abssys_widg.all_abssys[idx].lls_lines:
            iline.attrib['z'] = self.abssys_widg.all_abssys[idx].zabs 
            iline.attrib['N'] = self.abssys_widg.all_abssys[idx].NHI
            iline.attrib['b'] = self.abssys_widg.all_abssys[idx].bval
        # Update the rest
        self.update_model()
        self.draw()

    def update_boxes(self):
        '''Update Nbz boxes'''
        idx = self.get_sngl_sel_sys()
        if idx is None:
            return
        # z
        self.zwidget.box.setText(
            self.zwidget.box.frmt.format(
                self.abssys_widg.all_abssys[idx].zabs)) 
        # N
        self.Nwidget.box.setText(
            self.Nwidget.box.frmt.format(
                self.abssys_widg.all_abssys[idx].NHI)) 
        # b
        self.bwidget.box.setText(
            self.bwidget.box.frmt.format(
                self.abssys_widg.all_abssys[idx].bval.value)) 
        # Comment
        self.Cwidget.box.setText(
            self.Cwidget.box.frmt.format(
                self.abssys_widg.all_abssys[idx].comment))

    def update_conti(self):
        '''Update continuum '''
        self.continuum.flux = (self.base_continuum * self.conti_dict['Norm'] * 
            (self.continuum.dispersion.value/
                self.conti_dict['piv_wv'])**self.conti_dict['tilt'])
        if self.lls_model is not None:
            self.full_model.flux = self.lls_model * self.continuum.flux

    def update_model(self):
        '''Update absorption model '''
        if len(self.abssys_widg.all_abssys) == 0:
            self.lls_model = None
            self.spec_widg.model = None
            return
        #
        all_tau_model = xialu.tau_multi_lls(self.full_model.dispersion,
            self.abssys_widg.all_abssys)

        # Loop on forest lines
        for forest in self.all_forest:
            tau_Lyman = xsv.voigt_model(self.full_model.dispersion, 
                forest.lines, flg_ret=2)
            all_tau_model += tau_Lyman

        # Flux and smooth
        flux = np.exp(-1. * all_tau_model)
        if self.smooth > 0:
            self.lls_model = lsc.convolve_psf(flux, self.smooth)
        else:
            self.lls_model = flux

        # Finish
        self.full_model.flux = self.lls_model * self.continuum.flux
        # Over-absorbed
        self.spec_widg.bad_model = np.where( (self.lls_model < 0.7) &
            (self.full_model.flux < (self.spec_widg.spec.flux-
                self.spec_widg.spec.sig*1.5)))[0] 
        # Model
        self.spec_widg.model = self.full_model

    def get_sngl_sel_sys(self):
        '''Grab selected system
        '''
        items = self.abssys_widg.abslist_widget.selectedItems()
        if len(items) > 1:
            print('Need to select only 1 system!')
            return None
        # 
        item = items[0]
        txt = item.text()
        if txt == 'None':
            return None
        else:
            idx = self.abssys_widg.all_items.index(item.text())
            return idx

    def on_key(self,event):
        if event.key in ['C','1','2']: # Set continuum level
            if event.key == 'C':
                imin = np.argmin(np.abs(
                    self.continuum.dispersion.value-event.xdata))
                self.conti_dict['Norm'] = float(event.ydata / 
                    (self.base_continuum[imin].value*(event.xdata/
                        self.conti_dict['piv_wv'])**self.conti_dict['tilt']))
            elif event.key == '1':
                self.conti_dict['tilt'] += 0.1
            elif event.key == '2':
                self.conti_dict['tilt'] -= 0.1
            self.update_conti()
        elif event.key == 'A': # New LLS
            # Generate
            z = event.xdata/911.7633 - 1.
            self.add_LLS(z, bval=20.*u.km/u.s, NHI=17.3)
            # Update
            self.llist['Plot'] = False # Turn off metal-lines
            self.update_model()
        elif event.key in ['L','a','N','n','v','V','D','@','g']: # LLS-centric
            idx = self.get_sngl_sel_sys()
            if idx is None:
                return
            if event.key == 'L': #LLS
                self.abssys_widg.all_abssys[idx].zabs = event.xdata/911.7633 - 1.
            elif event.key == 'a': #Lya
                self.abssys_widg.all_abssys[idx].zabs = event.xdata/1215.6700-1.
            elif event.key == 'g': # Move nearest line to cursor
                wrest = event.xdata/(1+self.abssys_widg.all_abssys[idx].zabs)
                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()
                awrest = np.array([iline.wrest.value for iline in self.abssys_widg.all_abssys[idx].lls_lines])
                imn = np.argmin(np.abs(wrest-awrest))
                self.abssys_widg.all_abssys[idx].zabs = event.xdata/awrest[imn]-1.
            elif event.key == 'N': #Add to NHI
                self.abssys_widg.all_abssys[idx].NHI += 0.05
            elif event.key == 'n': #Subtract from NHI
                self.abssys_widg.all_abssys[idx].NHI -= 0.05
            elif event.key == 'v': #Subtract from bval
                self.abssys_widg.all_abssys[idx].bval -= 2*u.km/u.s
            elif event.key == 'V': #Add to bval
                self.abssys_widg.all_abssys[idx].bval += 2*u.km/u.s
            elif event.key == 'D': # Delete system
                self.abssys_widg.remove_item(idx)
                idx = None
            elif event.key == '@': # Toggle metal-lines
                self.llist['Plot'] = not self.llist['Plot']
                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()
            else:
                raise ValueError('Not ready for this keystroke')
            # Update the lines 
            if idx is not None:
                self.llist['z'] = self.abssys_widg.all_abssys[idx].zabs
                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()
                for iline in self.abssys_widg.all_abssys[idx].lls_lines:
                    iline.attrib['z'] = self.abssys_widg.all_abssys[idx].zabs 
                    iline.attrib['N'] = self.abssys_widg.all_abssys[idx].NHI
                    iline.attrib['b'] = self.abssys_widg.all_abssys[idx].bval
            # Update the model
            self.update_model()
        elif event.key in ['6','7','8','9']: # Add forest line
            self.add_forest(event.key,event.xdata/1215.6701 - 1.)
            self.update_model()
        else:
            self.spec_widg.on_key(event)

        # Draw by default
        self.update_boxes()
        self.draw()

    def draw(self):
        self.spec_widg.on_draw(no_draw=True)
        # Add text?
        for kk,lls in enumerate(self.abssys_widg.all_abssys):
            # Label
            ipos = self.abssys_widg.all_items[kk].rfind('_')
            ilbl = self.abssys_widg.all_items[kk][ipos+1:]
            # Add text
            for wv,lbl in self.plt_wv:
                idx = np.argmin(np.abs(self.continuum.dispersion-wv*(1+lls.zabs)))
                self.spec_widg.ax.text(wv.value*(1+lls.zabs), 
                    self.continuum.flux[idx],
                    '{:s}_{:s}'.format(ilbl,lbl), ha='center', 
                    color='blue', size='small', rotation=90.)
        # Draw
        self.spec_widg.canvas.draw()

    def add_forest(self,inp,z):
        '''Add a Lya/Lyb forest line
        '''
        from xastropy.igm.abs_sys.abssys_utils import GenericAbsSystem
        forest = GenericAbsSystem(zabs=z)
        # NHI
        NHI_dict = {'6':12.,'7':13.,'8':14.,'9':15.}
        forest.NHI=NHI_dict[inp]
        # Lines
        for name in ['HI 1215','HI 1025', 'HI 972']:
            aline = AbsLine(name,
                linelist=self.llist[self.llist['List']])
            # Attributes
            aline.attrib['N'] = forest.NHI
            aline.attrib['b'] = 20.*u.km/u.s
            aline.attrib['z'] = forest.zabs
            # Append
            forest.lines.append(aline)
        # Append to forest lines
        self.all_forest.append(forest)

    def add_LLS(self,z,NHI=17.3,bval=20.*u.km/u.s,comment='None'):
        '''Generate a new LLS
        '''
        #
        new_sys = LLSSystem(NHI=NHI)
        new_sys.zabs = z
        new_sys.bval = bval # This is not standard, but for convenience
        new_sys.comment = comment
        new_sys.fill_lls_lines(bval=bval)
        # Name
        self.count_lls += 1
        new_sys.label = 'LLS_Sys_{:d}'.format(self.count_lls)
        # Add
        self.abssys_widg.add_fil(new_sys.label)
        self.abssys_widg.all_abssys.append(new_sys)
        self.abssys_widg.abslist_widget.item(
            len(self.abssys_widg.all_abssys)).setSelected(True)

    def refine_abssys(self):
        item = self.abssys_widg.abslist_widget.selectedItems()
        if len(item) != 1:
            self.statusBar().showMessage('AbsSys: Must select only 1 system!')
            print('AbsSys: Must select only 1 system!')
        txt = item[0].text()
        ii = self.abssys_widg.all_items.index(txt)
        iabs_sys = self.abssys_widg.all_abssys[ii]
        # Launch
        gui = XVelPltGui(self.spec_widg.spec, outfil=iabs_sys.absid_file,
                               abs_sys=iabs_sys, norm=self.spec_widg.norm)
        gui.exec_()

    # Read from a JSON file
    def init_LLS(self,fit_file):
        import json
        # Read the JSON file
        with open(fit_file) as data_file:    
            lls_dict = json.load(data_file)
        # Init
        self.conti_dict = lls_dict['conti']
        self.update_conti()
        # Check spectra names
        if self.spec_widg.spec.filename != lls_dict['spec_file']:
            warnings.warn('Spec file names do not match!')
        # LLS
        for key in lls_dict['LLS'].keys():
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            self.add_LLS(lls_dict['LLS'][key]['z'],
                NHI=lls_dict['LLS'][key]['NHI'],
                bval=lls_dict['LLS'][key]['bval']*u.km/u.s,
                comment=lls_dict['LLS'][key]['comment'])
        self.smooth = lls_dict['smooth']
        # Updates
        self.update_boxes()
        self.update_model()
        QtCore.pyqtRemoveInputHook()
        xdb.set_trace()
        QtCore.pyqtRestoreInputHook()

    # Write
    def write_out(self):
        import json, io
        # Create dict
        out_dict = dict(LLS={},conti=self.conti_dict,
            spec_file=self.spec_widg.spec.filename,smooth=self.smooth)
        if self.zqso is not None:
            out_dict['zqso'] = self.zqso
        # Load
        for kk,lls in enumerate(self.abssys_widg.all_abssys):
            key = '{:d}'.format(kk+1)
            out_dict['LLS'][key] = {}
            out_dict['LLS'][key]['z'] = lls.zabs
            out_dict['LLS'][key]['NHI'] = lls.NHI
            out_dict['LLS'][key]['bval'] = lls.lls_lines[0].attrib['b'].value
            out_dict['LLS'][key]['comment'] = str(lls.comment).strip()
        # Write
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        with io.open(self.outfil, 'w', encoding='utf-8') as f:
            f.write(unicode(json.dumps(out_dict, sort_keys=True, indent=4, 
                separators=(',', ': '))))
        self.flag_write = True

    # Write + Quit
    def write_quit(self):
        self.write_out()
        self.quit()

    # Quit
    def quit(self):
        self.close()


# Script to run XSpec from the command line or ipython
def run_xspec(*args, **kwargs):
    '''
    Runs the XSpecGui

    Command line
    or from Python
    Examples:
      1.  python ~/xastropy/xastropy/xguis/spec_guis.py 1
      2.  spec_guis.run_xspec(filename)
      3.  spec_guis.run_xspec(spec1d)
    '''

    import argparse
    from specutils import Spectrum1D

    parser = argparse.ArgumentParser(description='Parse for XSpec')
    parser.add_argument("flag", type=int, help="GUI flag (ignored)")
    parser.add_argument("file", type=str, help="Spectral file")
    parser.add_argument("-zsys", type=float, help="System Redshift")
    parser.add_argument("--un_norm", help="Spectrum is NOT normalized",
                        action="store_true")

    if len(args) == 0:
        pargs = parser.parse_args()
    else: # better know what you are doing!
        if isinstance(args[0],(Spectrum1D,tuple)):
            app = QtGui.QApplication(sys.argv)
            gui = XSpecGui(args[0], **kwargs)
            gui.show()
            app.exec_()
            return
        else: # String parsing 
            largs = ['1'] + [iargs for iargs in args]
            pargs = parser.parse_args(largs)
    
    # Normalized?
    norm=True
    if pargs.un_norm:
        norm=False

    # Second spectral file?
    try:
        zsys = pargs.zsys
    except AttributeError:
        zsys=None

    app = QtGui.QApplication(sys.argv)
    gui = XSpecGui(pargs.file, zsys=zsys, norm=norm)
    gui.show()
    app.exec_()

# Script to run XAbsID from the command line
def run_xabsid():

    import argparse

    parser = argparse.ArgumentParser(description='Script for XSpec')
    parser.add_argument("flag", type=int, help="GUI flag (ignored)")
    parser.add_argument("file", type=str, help="Spectral file")
    parser.add_argument("--un_norm", help="Spectrum is NOT normalized",
                        action="store_true")
    parser.add_argument("-id_dir", type=str,
                        help="Directory for ID files (ID_LINES is default)")
    parser.add_argument("-secondfile", type=str, help="Second spectral file")
    parser.add_argument("-thirdfile", type=str, help="Third spectral file")
    
    args = parser.parse_args()

    # Normalized?
    norm=True
    if args.un_norm:
        norm=False

    infile = args.file
    # Second spectral file?
    if args.secondfile:
        infile = [infile, args.secondfile]
    if args.thirdfile:
        infile.append(args.thirdfile)

    # Launch
    app = QtGui.QApplication(sys.argv)
    gui = XAbsIDGui(infile, norm=norm)
    gui.show()
    app.exec_()

# Script to run XVelPltGui from the command line or ipython
def run_xvelp(*args, **kwargs):
    '''
    Runs the XVelPltGui

    Command line
    or from Python
    Examples:
      1.  python ~/xastropy/xastropy/xguis/spec_guis.py 3
      2.  spec_guis.run_xvelp(filename)
      3.  spec_guis.run_xvelp(spec1d)
    '''

    import argparse
    from specutils import Spectrum1D

    xdb.set_trace() # DEPRECATED FOR NOW

    parser = argparse.ArgumentParser(description='Parse for XVelPlt')
    parser.add_argument("flag", type=int, help="GUI flag (ignored)")
    parser.add_argument("file", type=str, help="Spectral file")
    parser.add_argument("-zsys", type=float, help="System Redshift")
    parser.add_argument("--un_norm", help="Spectrum is NOT normalized",
                        action="store_true")

    if len(args) == 0:
        pargs = parser.parse_args()
    else: # better know what you are doing!
        if isinstance(args[0],(Spectrum1D,tuple)):
            if not kwargs['rerun']:
                app = QtGui.QApplication(sys.argv)
            xdb.set_trace()
            gui = XVelPltGui(args[0], **kwargs)
            gui.exec_()
            #gui.show()
            #app.exec_()
            return gui, app
        else: # String parsing 
            largs = ['1'] + [iargs for iargs in args]
            pargs = parser.parse_args(largs)
    
    # Normalized?
    norm=True
    if pargs.un_norm:
        norm=False

    # Second spectral file?
    try:
        zsys = pargs.zsys
    except AttributeError:
        zsys=None

    xdb.set_trace() # Not setup for command line yet
    app = QtGui.QApplication(sys.argv)
    gui = XSpecGui(pargs.file, zsys=zsys, norm=norm)
    gui.show()
    app.exec_()

    return gui

# Script to run XSpec from the command line or ipython
def run_fitlls(*args, **kwargs):
    '''
    Runs the XFitLLSGUI

    Command line or from Python
    Examples:
      1.  python ~/xastropy/xastropy/xguis/spec_guis.py 1
      2.  spec_guis.run_fitlls(filename)
      3.  spec_guis.run_fitlls(spec1d)
    '''

    import argparse
    from specutils import Spectrum1D

    parser = argparse.ArgumentParser(description='Parser for XFitLLSGUI')
    parser.add_argument("flag", type=int, help="GUI flag (ignored)")
    parser.add_argument("in_file", type=str, help="Spectral file")
    parser.add_argument("-out_file", type=str, help="Output LLS Fit file")
    parser.add_argument("-smooth", type=float, help="Smoothing (pixels)")
    parser.add_argument("-lls_fit_file", type=str, help="Input LLS Fit file")
    parser.add_argument("-zqso", type=float, help="Use Telfer template with zqso")

    if len(args) == 0:
        pargs = parser.parse_args()
    else: # better know what you are doing!
        if isinstance(args[0],(Spectrum1D,tuple)):
            app = QtGui.QApplication(sys.argv)
            gui = XFitLLSGUI(args[0], **kwargs)
            gui.show()
            app.exec_()
            return
        else: # String parsing 
            largs = ['1'] + [iargs for iargs in args]
            pargs = parser.parse_args(largs)

    # Output file
    try:
        outfil = pargs.out_file
    except AttributeError:
        outfil=None

    # Input LLS file
    try:
        lls_fit_file = pargs.lls_fit_file
    except AttributeError:
        lls_fit_file=None

    # Smoothing parameter
    try:
        smooth = pargs.smooth
    except AttributeError:
        smooth=3.

    # Smoothing parameter
    try:
        zqso = pargs.zqso
    except AttributeError:
        zqso=None
    
    app = QtGui.QApplication(sys.argv)
    gui = XFitLLSGUI(pargs.in_file,outfil=outfil,smooth=smooth,
        lls_fit_file=lls_fit_file, zqso=zqso)
    gui.show()
    app.exec_()

# ################
if __name__ == "__main__":
    import sys
    from linetools.spectra import io as lsi
    from xastropy.igm import abs_sys as xiabs

    if len(sys.argv) == 1: # TESTING

        flg_fig = 0 
        #flg_fig += 2**0  # XSpec
        #flg_fig += 2**1  # XAbsID
        #flg_fig += 2**2  # XVelPlt Gui
        #flg_fig += 2**3  # XVelPlt Gui without ID list; Also tests select wave
        #flg_fig += 2**4  # XAODM Gui
        flg_fig += 2**5  # Fit LLS GUI
    
        # Read spectrum
        spec_fil = '/u/xavier/Keck/HIRES/RedData/PH957/PH957_f.fits'
        spec = lsi.readspec(spec_fil)
    
        # XSpec
        if (flg_fig % 2) == 1:
            app = QtGui.QApplication(sys.argv)
            gui = XSpecGui(spec)
            gui.show()
            app.exec_()
    
        # XAbsID
        if (flg_fig % 2**2) >= 2**1:
            #spec_fil = '/u/xavier/PROGETTI/LLSZ3/data/normalize/SDSSJ1004+0018_nF.fits'
            #spec = xspec.readwrite.readspec(spec_fil)
            #norm = True
            spec_fil = '/Users/xavier/Dropbox/CASBAH/jxp_analysis/FBQS0751+2919/fbqs0751_nov2014bin.fits'
            norm = False
            absid_fil = '/Users/xavier/paper/LLS/Optical/Data/Analysis/MAGE/SDSSJ1004+0018_z2.746_id.fits'
            absid_fil2 = '/Users/xavier/paper/LLS/Optical/Data/Analysis/MAGE/SDSSJ2348-1041_z2.997_id.fits'

            app = QtGui.QApplication(sys.argv)
            gui = XAbsIDGui(spec_fil,norm=norm) #,absid_list=[absid_fil, absid_fil2])
            gui.show()
            app.exec_()

        # XVelPlt with existing AbsID file
        if (flg_fig % 2**3) >= 2**2:
            spec_fil = '/u/xavier/PROGETTI/LLSZ3/data/normalize/SDSSJ1004+0018_nF.fits'
            #spec = xspec.readwrite.readspec(spec_fil)
            absid_fil = '/Users/xavier/paper/LLS/Optical/Data/Analysis/MAGE/SDSSJ1004+0018_z2.746_id.fits'
            abs_sys = xiabs.abssys_utils.Generic_System(None)
            abs_sys.parse_absid_file(absid_fil)
            #
            app = QtGui.QApplication(sys.argv)
            app.setApplicationName('XVelPlt')
            gui = XVelPltGui(spec_fil,abs_sys=abs_sys,
                             outfil='/Users/xavier/Desktop/tmp.fits')
            gui.show()
            sys.exit(app.exec_())

        # XVelPlt without existing AbsID file
        if (flg_fig % 2**4) >= 2**3:
            #spec_fil = '/u/xavier/PROGETTI/LLSZ3/data/normalize/SDSSJ1004+0018_nF.fits'
            #z=2.746
            #outfil='/Users/xavier/Desktop/J1004+0018_z2.746_id.fits'
            spec_fil = '/Users/xavier/Dropbox/CASBAH/jxp_analysis/FBQS0751+2919/fbqs0751_nov2014bin.fits'
            z=0.
            outfil='/Users/xavier/Desktop/tmp.fits'
            #
            app = QtGui.QApplication(sys.argv)
            app.setApplicationName('XVelPlt')
            gui = XVelPltGui(spec_fil, z=z, outfil=outfil,norm=False, sel_wv=1526.80)
            gui.show()
            sys.exit(app.exec_())

        # AODM GUI
        if (flg_fig % 2**5) >= 2**4:
            #spec_fil = '/Users/xavier/PROGETTI/LLSZ3/data/normalize/UM184_nF.fits'
            #z=2.96916
            #lines = [1548.195, 1550.770]
            norm = True
            spec_fil = '/Users/xavier/Dropbox/CASBAH/jxp_analysis/FBQS0751+2919/fbqs0751_nov2014bin.fits'
            z=0.4391
            lines = [1215.6701, 1025.7223] * u.AA
            norm = False
            # Launch
            spec = lsi.readspec(spec_fil)
            app = QtGui.QApplication(sys.argv)
            app.setApplicationName('AODM')
            main = XAODMGui(spec, z, lines, norm=norm)
            main.show()
            sys.exit(app.exec_())

        # LLS
        if (flg_fig % 2**6) >= 2**5:
            #spec_fil = '/Users/xavier/PROGETTI/LLSZ3/data/normalize/UM184_nF.fits'
            #z=2.96916
            #lines = [1548.195, 1550.770]
            spec_fil = '/Users/xavier/VLT/XShooter/LP/idl_reduced_frames/0952-0115_uvb_coadd_vbin_flx.fits'
            # Launch
            spec = lsi.readspec(spec_fil)
            app = QtGui.QApplication(sys.argv)
            app.setApplicationName('FitLLS')
            main = XFitLLSGUI(spec) 
            main.show()
            sys.exit(app.exec_())

    else: # RUN A GUI
        id_gui = int(sys.argv[1])  # 1 = XSpec, 2=XAbsId

        if id_gui == 1:
            run_xspec()
        elif id_gui == 2:
            run_xabsid()
        elif id_gui == 3:
            run_xvelp()
        elif id_gui == 4:
            run_fitlls()
        else:
            raise ValueError('Unsupported flag for spec_guis')
            