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
import glob

from PyQt4 import QtGui
from PyQt4 import QtCore

from matplotlib import mpl
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
# Matplotlib Figure object
from matplotlib.figure import Figure

from astropy.units import Quantity
from astropy import units as u

from linetools.lists.linelist import LineList

from xastropy.xutils import xdebug as xdb
from xastropy.xguis import spec_widgets as xspw

#class XSpecGui(QtGui.QMainWindow):
#class XAbsIDGui(QtGui.QMainWindow):
#class XVelPltGui(QtGui.QDialog):

# x_specplot replacement
class XSpecGui(QtGui.QMainWindow):
    ''' GUI to replace XIDL x_specplot

        12-Dec-2014 by JXP v1.0
        27-Mar-2015 by JXP v2.0 :: EW, column, better zooming + panning
    '''
    def __init__(self, spec, parent=None, zsys=None, norm=None):
        QtGui.QMainWindow.__init__(self, parent)
        '''
        spec = Spectrum1D
        '''
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
    def __init__(self, spec, parent=None, abssys_dir=None, absid_list=None, norm=True,
                 srch_id=True, id_dir='ID_LINES/', second_file=None):
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
        self.spec_widg = xspw.ExamineSpecWidget(spec,status=self.statusBar,
                                                llist=self.pltline_widg.llist, norm=norm,
                                                second_file=second_file, 
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
                 vmnx=[-300., 300.], abs_sys=None, outfil='dum_ID.fits',
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
        if not self.abs_sys is None:
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
        # Quit buttons
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
        all_lines = list( llist[llist['List']]['wrest'] )
        # Set selected
        abs_sys = self.vplt_widg.abs_sys
        wrest = abs_sys.lines.keys()
        wrest.sort()
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
    def __init__(self, spec, z, wrest, vmnx=[-300., 300.], parent=None, norm=True):
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



# Script to run XSpec from the command line
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
    from xastropy.spec.utils import XSpectrum1D

    parser = argparse.ArgumentParser(description='Parse for XSpec')
    parser.add_argument("flag", type=int, help="GUI flag (ignored)")
    parser.add_argument("file", type=str, help="Spectral file")
    parser.add_argument("-zsys", type=float, help="System Redshift")
    parser.add_argument("--un_norm", help="Spectrum is NOT normalized",
                        action="store_true")

    if len(args) == 0:
        pargs = parser.parse_args()
    else: # better know what you are doing!
        #xdb.set_trace()
        if type(args[0]) in [XSpectrum1D, Spectrum1D]: 
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
    
    args = parser.parse_args()

    # Normalized?
    norm=True
    if args.un_norm:
        norm=False

    # Second spectral file?
    second_file=None
    if args.secondfile:
        second_file=args.secondfile

    # Launch
    app = QtGui.QApplication(sys.argv)
    gui = XAbsIDGui(args.file, norm=norm, second_file=second_file)
    gui.show()
    app.exec_()


# ################
if __name__ == "__main__":
    import sys
    from linetools.spectra import io as lsi
    from xastropy.igm import abs_sys as xiabs

    if len(sys.argv) == 1: # TESTING

        flg_fig = 0 
        flg_fig += 2**0  # XSpec
        #flg_fig += 2**1  # XAbsID
        #flg_fig += 2**2  # XVelPlt Gui
        #flg_fig += 2**3  # XVelPlt Gui without ID list; Also tests select wave
        #flg_fig += 2**4  # XAODM Gui
    
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

    else: # RUN A GUI
        id_gui = int(sys.argv[1])  # 1 = XSpec, 2=XAbsId

        if id_gui == 1:
            run_xspec()
        elif id_gui == 2:
            run_xabsid()
