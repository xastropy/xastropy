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
import linetools.spectra.io as lsi
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

# class XFitLLSGUI(QtGui.QMainWindow):

'''
=======
Analyzing spectra with auto_plls

Here is now my preferred approach to searching for
LLS with auto_plls:

1.  Load up the spectrum.  Fiddle with the continuum
normalization (and tilt, if necessary).

2.  Eyeball search for a putative break that one might
associate with a PLLS.

3.  Put the cursor a bit to the left (blueward) of the break
and at the approximate flux level of the data, post-break.
The former sets the starting redshift to search and the latter
sets a guess for NHI.

4. Hit "F" and wait for magic to happen.

5. If an LLS satisfying a rather simple criterion is found, it 
should appear.  Otherwise, a message prints to the terminal
stating none found.  You should then inspect the model,
including at Lyb and Lya and modify it (or even delete it).

6. Once you are happy with what you got (hopefully you are),
move on to the next putative break and try again, i.e.
repeat steps 2-5.
'''

# GUI for fitting LLS in a spectrum
class XFitLLSGUI(QtGui.QMainWindow):
    ''' GUI to fit LLS in a given spectrum
        v1.1
        30-Jul-2015 by JXP
    '''
    def __init__(self, ispec, parent=None, lls_fit_file=None, 
        outfil=None, smooth=3., zqso=None, fN_gamma=None, template=None): 
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
        template: str, optional
          Filename of a QSO template to use instead of the Telfer
          continuum. Only used if zqso is also given.
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
            if template is not None:
                tspec = lsi.readspec(template)
                # assume wavelengths
                tspec = XSpectrum1D.from_tuple(
                    (tspec.dispersion.value * (1 + zqso),
                    tspec.flux.value))
            else:
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
        self.plt_wv = zip(np.array([911.7, 949.743, 972.5367,1025.7222,1215.6700])*u.AA,
            ['LL','Lyd', 'Lyg','Lyb','Lya'])

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
        if len(items) == 0:
            return None
        elif len(items) > 1:
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
        elif event.key == 'F': # New LLS
            self.auto_plls(event.xdata, event.ydata)
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
        # Ticks for selected LLS
        idxl = self.get_sngl_sel_sys()
        if idxl is not None:
            lls = self.abssys_widg.all_abssys[idxl]
            # Label
            ipos = self.abssys_widg.all_items[idxl].rfind('_')
            ilbl = self.abssys_widg.all_items[idxl][ipos+1:]
            for line in lls.lls_lines:
                if line.wrest < 915.*u.AA:
                    continue
                idx = np.argmin(np.abs(self.continuum.dispersion-
                    line.wrest*(1+lls.zabs)))
                self.spec_widg.ax.text(line.wrest.value*(1+lls.zabs), 
                    self.continuum.flux[idx],
                    '-{:s}'.format(ilbl), ha='center', 
                    color='red', size='small', rotation=90.)
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
        
        # Update
        self.llist['Plot'] = False # Turn off metal-lines
        self.update_model()

    def auto_plls(self,x,y):
        '''Automatically fit a pLLS
        Parameters:
        ----------
        x,y: floats
          x,y values in the GUI
        '''
        spec = self.spec_widg.spec # For convenience
        if len(self.abssys_widg.all_abssys) > 0:
            conti= self.full_model
        else:
            conti= self.continuum
        # Generate toy LLS from click
        ximn = np.argmin(np.abs(spec.dispersion.value-x))
        NHI = 16.29 - np.log(y/self.continuum.flux.value[ximn])
        #print('NHI={:g}'.format(NHI))
        plls = LLSSystem(NHI=NHI)
        plls.zabs = x/(911.7)-1
        plls.bval = 20*u.km/u.s
        plls.fill_lls_lines(bval=20*u.km/u.s)

        # wrest, Tau model, flux
        wrest = spec.dispersion/(1+plls.zabs)
        tau = xialu.tau_multi_lls(spec.dispersion,[plls])
        emtau = np.exp(-1. * tau)
        lls_flux = lsc.convolve_psf(emtau, 3.)
#xdb.xplot(wrest, lls_flux)

        # zmin
        if len(self.abssys_widg.all_abssys) != 0:
            zlls = [lls.zabs for lls in self.abssys_widg.all_abssys]
            zmin = np.min(np.array(zlls)) - 0.01
        else:
            zmin = self.zqso+0.01

        # Pixels for analysis and rolling
        # NEED TO CUT ON X-Shooter ARM
        apix = np.where( (wrest > 914*u.AA) & #(spec.dispersion<5600*u.AA) &
                        (spec.dispersion<(1+zmin)*1026.*u.AA))[0] # Might go to Lyb
        nroll = (np.argmin(np.abs(spec.dispersion-(911.7*u.AA*(1+zmin))))- # Extra 0.01 for bad z
                   np.argmin(np.abs(spec.dispersion-(911.7*u.AA*(1+plls.zabs)))))
        gdpix = np.arange(np.min(apix)-nroll,np.max(apix)+nroll+1)
        #print(len(apix), nroll)
        roll_flux = np.concatenate([np.ones(nroll),lls_flux[apix], np.ones(nroll)])
        roll_msk = roll_flux < 0.7

        # Generate data arrays
        wave_pad = spec.dispersion[gdpix]
        flux_pad = spec.flux[gdpix]
        sig_pad = spec.sig[gdpix]
        if len(self.abssys_widg.all_abssys) > 0:
            conti_pad = conti.flux[gdpix]
        else:
            conti_pad = conti.flux[gdpix]

        # Generate matricies
        flux_matrix = np.zeros((len(roll_flux),nroll))
        sig_matrix = np.zeros((len(roll_flux),nroll))
        conti_matrix = np.zeros((len(roll_flux),nroll))

        roll_matrix = np.zeros((len(roll_flux),nroll))
        mask_matrix = np.zeros((len(roll_flux),nroll))
        for kk in range(nroll):
            roll_matrix[:,kk] = np.roll(roll_flux,kk)
            mask_matrix[:,kk] = np.roll(roll_msk,kk)
            flux_matrix[:,kk] = flux_pad
            conti_matrix[:,kk] = conti_pad
            sig_matrix[:,kk] = sig_pad

        # Model -- Multiply by continuum
        model = roll_matrix * conti_matrix

        # Condition
        idx = np.where( (model < (flux_matrix-sig_matrix*1.5)) & (mask_matrix==True))
        bad_matrix = np.zeros((len(roll_flux),nroll))
        bad_matrix[idx] = 1

        # Sum on offsets and get redshift
        bad = np.sum(bad_matrix,0)
        ibest = np.argmin(bad)
        zbest = spec.dispersion[ibest+ximn]/(911.7*u.AA)-1 # Quantity

        # Add pLLS?
        if bad[ibest] < 10:
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            self.add_LLS(zbest.value, bval=20.*u.km/u.s, NHI=NHI)
        else:
            print('No viable pLLS found with our criteria!')


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
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

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
    from xastropy.igm import abs_sys as xiabs

    if len(sys.argv) == 1: # TESTING

        flg_tst = 0 
        flg_tst += 2**0  # Fit LLS GUI
    
        # LLS
        if (flg_tst % 2**1) >= 2**0:
            spec_fil = '/Users/xavier/VLT/XShooter/LP/idl_reduced_frames/0952-0115_uvb_coadd_vbin_flx.fits'
            # Launch
            spec = lsi.readspec(spec_fil)
            app = QtGui.QApplication(sys.argv)
            app.setApplicationName('FitLLS')
            main = XFitLLSGUI(spec) 
            main.show()
            sys.exit(app.exec_())

    else: # RUN A GUI
        run_fitlls()
            
