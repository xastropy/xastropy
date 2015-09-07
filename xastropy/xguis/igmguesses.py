"""
#;+ 
#; NAME:
#; igmguesses
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for LineIDs and Initial guesses in IGM spectra with QT
#;      Likely only useful for lowz-IGM
#;   14-Aug-2015 by JXP
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
from astropy import units as u
from astropy.io import fits, ascii
from astropy import constants as const

from linetools.lists.linelist import LineList
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectra import convolve as lsc
from linetools.spectralline import AbsLine

#from xastropy.atomic import ionization as xatomi
from xastropy.plotting import utils as xputils
from xastropy.xguis import spec_widgets as xspw
from xastropy.xguis import utils as xxgu
from xastropy.spec import voigt as xsv

from xastropy.xutils import xdebug as xdb

xa_path = imp.find_module('xastropy')[1]



#class IGMGuessesGui(QtGui.QMainWindow):

# GUI for fitting LLS in a spectrum
class IGMGuessesGui(QtGui.QMainWindow):
    ''' GUI to fit LLS in a given spectrum
        v0.5
        30-Jul-2015 by JXP
    '''
    def __init__(self, ispec, parent=None, previous_file=None, 
        srch_id=True, outfil=None, fwhm=None, zqso=None,
        plot_residuals=True):
        QtGui.QMainWindow.__init__(self, parent)
        '''
        spec = Spectrum1D
        previous_file: str, optional
          Name of the previous guesses file
        smooth: float, optional
          Number of pixels to smooth on
        zqso: float, optional
          Redshift of the quasar.  If input, a Telfer continuum is used
        plot_residuals : bool, optional
          Whether to plot residuals

        '''
        # TODO
        # 1. Fix convolve window size
            # 2. Avoid sorting of wavelengths
            # 3. Remove astropy.relativity 
            # 4. Display self.z 
            # 5. Recenter on 'z'
        # 6. Add COS LSF
            # 7. Refit component key 'F'
            # 8. Plot only selected components [DEPRECATED]
            # 9. Avoid shifting z of component outside its velocity range
            # 10. Write Component vlim to JSON
            # 11. Key to set line as some other transition (e.g. RMB in XSpec)
            # 12. Mask array with points
            # 13. Toggle line ID names
            # 14. Add error + residual arrays [NT]
        # 15. Adjust component redshift by keystroke
            # 16. Input redshift value via Widget
            # 17. Use Component list to jump between components (like 'S')

        # Build a widget combining several others
        self.main_widget = QtGui.QWidget()

        # Status bar
        self.create_status_bar()

        # Initialize
        self.previous_file = previous_file
        if outfil is None:
            self.outfil = 'IGM_model.json'
        else:
            self.outfil = outfil
        if fwhm is None:
            self.fwhm = 3.
        else:
            self.fwhm = fwhm
        self.plot_residuals = plot_residuals

        # Spectrum
        spec, spec_fil = xxgu.read_spec(ispec)
        spec.normalize() # Need to change this..
        spec.mask = np.zeros(len(spec.dispersion),dtype=int)

        # Normalize

        # Full Model 
        self.model = XSpectrum1D.from_tuple((
            spec.dispersion,np.ones(len(spec.dispersion))))

        # LineList (Grab ISM and HI as defaults)
        self.llist = xxgu.set_llist('ISM')
        # Load others
        self.llist['HI'] = LineList('HI')
        self.llist['Strong'] = LineList('Strong')
        self.llist['Lists'].append('HI')
        self.llist['HI']._data = self.llist['HI']._data[::-1] #invert order of Lyman series
        #self.llist['show_line'] = np.arange(10) #maximum 10 to show for Lyman series
        
        #define initial redshift
        z=0.0
        self.llist['z'] = z
        
        # Grab the pieces and tie together
        self.slines_widg = xspw.SelectedLinesWidget(
            self.llist[self.llist['List']], parent=self, init_select='All')
        self.fiddle_widg = FiddleComponentWidget(parent=self)
        self.comps_widg = ComponentListWidget([], parent=self)
        self.velplot_widg = IGGVelPlotWidget(spec, z, 
            parent=self, llist=self.llist, fwhm=self.fwhm,plot_residuals=self.plot_residuals)
        self.wq_widg = xxgu.WriteQuitWidget(parent=self)
        

        # Setup strongest LineList
        self.llist['strongest'] = LineList('ISM')
        self.llist['Lists'].append('strongest')
        self.update_strongest_lines()
        self.slines_widg.selected = self.llist['show_line']
        self.slines_widg.on_list_change(
            self.llist[self.llist['List']])

        # Load prevoius file
        if self.previous_file is not None:
            self.read_previous()

        # Connections (buttons are above)
        #self.spec_widg.canvas.mpl_connect('key_press_event', self.on_key)
        #self.abssys_widg.abslist_widget.itemSelectionChanged.connect(
        #    self.on_list_change)

        # Layout
        anly_widg = QtGui.QWidget()
        anly_widg.setMaximumWidth(400)
        anly_widg.setMinimumWidth(250)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.fiddle_widg)
        vbox.addWidget(self.comps_widg)
        vbox.addWidget(self.slines_widg)
        vbox.addWidget(self.wq_widg)
        anly_widg.setLayout(vbox)
        
        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(self.velplot_widg)
        hbox.addWidget(anly_widg)

        self.main_widget.setLayout(hbox)

        # Point MainWindow
        self.setCentralWidget(self.main_widget)

    def update_strongest_lines(self):
        '''Grab the strongest lines in the spectrum at the current
        redshift.
        '''
        z = self.velplot_widg.z
        wvmin = np.min(self.velplot_widg.spec.dispersion) 
        wvmax = np.max(self.velplot_widg.spec.dispersion) 
        wvlims = (wvmin/(1+z),wvmax/(1+z))
        transitions = self.llist['ISM'].available_transitions(
            wvlims,n_max=100, n_max_tuple=4,min_strength=5.)

        if transitions is not None:
            names = list(np.array(transitions['name']))
        else:
            names = ['HI 1215']
        self.llist['strongest'].subset_lines(reset_data=True,subset=names)
        self.llist['show_line'] = np.arange(len(self.llist['strongest']._data))
        self.llist['List'] = 'strongest'
        # self.llist['strongest'] = self.llist['ISM'].subset(names)


    def on_list_change(self):
        self.update_boxes()

    def create_status_bar(self):
        self.status_text = QtGui.QLabel("IGMGuessesGui")
        self.statusBar().addWidget(self.status_text, 1)

    def delete_component(self, component):
        '''Remove component'''
        # Component list
        self.comps_widg.remove_item(component.name)
        # Fiddle query (need to come back to it)
        if component is self.fiddle_widg.component:
            self.fiddle_widg.reset()
        # Mask
        for line in component.lines:
            wvmnx = line.wrest*(1+component.zcomp)*(1 + component.vlim.value/3e5)
            gdp = np.where((self.velplot_widg.spec.dispersion>wvmnx[0])&
                (self.velplot_widg.spec.dispersion<wvmnx[1]))[0]
            self.velplot_widg.spec.mask[gdp] = 0
        # Delete
        del component
        # Update
        self.velplot_widg.update_model()
        self.velplot_widg.on_draw(fig_clear=True)

    def updated_slines(self, selected):
        self.llist['show_line'] = selected
        self.velplot_widg.on_draw(fig_clear=True)

    def updated_component(self):
        '''Component attrib was updated. Deal with it'''
        self.fiddle_widg.component.sync_lines()
        self.velplot_widg.update_model()
        self.velplot_widg.on_draw(fig_clear=True)

    def updated_compslist(self,component):
        '''Component list was updated'''
        if component is None:
            self.fiddle_widg.reset()
        else:
            self.fiddle_widg.init_component(component)
        #self.velplot_widg.update_model()
        #self.velplot_widg.on_draw(fig_clear=True)

    def read_previous(self):
        ''' Read from a previous guesses file'''
        import json
        # Read the JSON file
        with open(self.previous_file) as data_file:    
            igmg_dict = json.load(data_file)
        # Check FWHM
        if igmg_dict['fwhm'] != self.fwhm:
            raise ValueError('Input FWHMs do not match.  Fix')
        # Mask
        msk = igmg_dict['mask']
        if len(msk) > 0:
            self.velplot_widg.spec.mask[np.array(msk)] = 1
        # Check spectra names
        if self.velplot_widg.spec.filename != igmg_dict['spec_file']:
            warnings.warn('Spec file names do not match! Could just be path..')
        # Components
        for key in igmg_dict['cmps'].keys():
            self.velplot_widg.add_component(
                igmg_dict['cmps'][key]['wrest']*u.AA, 
                zcomp=igmg_dict['cmps'][key]['zcomp'],
                vlim=igmg_dict['cmps'][key]['vlim']*u.km/u.s,
                no_fit_mask=True)
            # Name
            self.velplot_widg.current_comp.name = key
            # Set N,b,z
            self.velplot_widg.current_comp.attrib['z']= igmg_dict['cmps'][key]['zfit']
            self.velplot_widg.current_comp.attrib['b']= igmg_dict['cmps'][key]['bval']*u.km/u.s
            self.velplot_widg.current_comp.attrib['N']= igmg_dict['cmps'][key]['N']
            self.velplot_widg.current_comp.attrib['Quality']= igmg_dict['cmps'][key]['Quality']
            self.velplot_widg.current_comp.comment = igmg_dict['cmps'][key]['comment']
            # Sync
            self.velplot_widg.current_comp.sync_lines()
        # Updates
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        self.velplot_widg.update_model()
        self.fiddle_widg.init_component(self.velplot_widg.current_comp)


    def write_out(self):
        ''' Write to a JSON file'''
        import json, io
        # Create dict of the components
        out_dict = dict(cmps={},
            spec_file=self.velplot_widg.spec.filename,
            fwhm=self.fwhm)
        mskp = np.where(self.velplot_widg.spec.mask == 1)[0]
        if len(mskp) > 0:
            out_dict['mask'] = list(mskp)
        # Load
        for kk,comp in enumerate(self.comps_widg.all_comp):
            key = comp.name
            out_dict['cmps'][key] = {}
            out_dict['cmps'][key]['zcomp'] = comp.zcomp
            out_dict['cmps'][key]['zfit'] = comp.attrib['z']
            out_dict['cmps'][key]['N'] = comp.attrib['N']
            out_dict['cmps'][key]['bval'] = comp.attrib['b'].value
            out_dict['cmps'][key]['wrest'] = comp.init_wrest.value
            out_dict['cmps'][key]['vlim'] = list(comp.vlim.value)
            out_dict['cmps'][key]['Quality'] = str(comp.attrib['Quality'])
            out_dict['cmps'][key]['comment'] = str(comp.comment)
        # Write
        print('Wrote: {:s}'.format(self.outfil))
        with io.open(self.outfil, 'w', encoding='utf-8') as f:
            f.write(unicode(json.dumps(out_dict, sort_keys=True, indent=4, 
                separators=(',', ': '))))

    # Write + Quit
    def write_quit(self):
        self.write_out()
        self.quit()

    # Quit
    def quit(self):
        self.close()

 ######################
class IGGVelPlotWidget(QtGui.QWidget):
    ''' Widget for a velocity plot with interaction.
          Adapted from VelPlotWidget in spec_guis
        14-Aug-2015 by JXP
    '''
    def __init__(self, ispec, z, parent=None, llist=None, norm=True,
                 vmnx=[-500., 500.]*u.km/u.s, fwhm=0.,plot_residuals=True):
        '''
        spec = Spectrum1D
        Norm: Bool (False)
          Normalized spectrum?
        abs_sys: AbsSystem
          Absorption system class
        '''
        super(IGGVelPlotWidget, self).__init__(parent)

        # Initialize
        self.parent = parent
        spec, spec_fil = xxgu.read_spec(ispec)
        
        self.spec = spec
        self.spec_fil = spec_fil
        self.fwhm = fwhm
        self.z = z
        self.vmnx = vmnx
        self.norm = norm
        # Init
        self.flag_add = False
        self.flag_idlbl = False
        self.flag_mask = False
        self.wrest = 0.
        self.avmnx = np.array([0.,0.])*u.km/u.s
        self.model = XSpectrum1D.from_tuple((
            spec.dispersion,np.ones(len(spec.dispersion))))

        self.plot_residuals = plot_residuals
        #Define arrays for plotting residuals
        if self.plot_residuals:
            self.residual_normalization_factor = 0.02/np.median(self.spec.sig)
            self.residual_limit = self.spec.sig * self.residual_normalization_factor
            self.residual = (self.spec.flux - self.model.flux) * self.residual_normalization_factor


        self.psdict = {} # Dict for spectra plotting
        self.psdict['xmnx'] = self.vmnx.value # Too much pain to use units with this
        self.psdict['ymnx'] = [-0.1, 1.1]
        self.psdict['nav'] = xxgu.navigate(0,0,init=True)

        

        # Status Bar?
        #if not status is None:
        #    self.statusBar = status

        # Line List
        if llist is None:
            self.llist = xxgu.set_llist(['HI 1215', 'HI 1025'])
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
        self.sub_xy = [5,3]
        self.subxy_state = 'Out'

        self.fig.subplots_adjust(hspace=0.0, wspace=0.1,left=0.04,right=0.975)
        
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
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        wrest = self.llist[self.llist['List']].wrest
        wvobs = (1+self.z) * wrest
        gdlin = np.where( (wvobs > wvmin) & (wvobs < wvmax) )[0]
        self.llist['show_line'] = gdlin
        # Update GUI
        self.parent.slines_widg.selected = self.llist['show_line']
        self.parent.slines_widg.on_list_change(
            self.llist[self.llist['List']])


    # Add a component
    def update_model(self):
        if self.parent is None:
            return
        all_comp = self.parent.comps_widg.all_comp #selected_components()
        if len(all_comp) == 0:
            self.model.flux[:] = 1.
            return
        # Setup lines
        wvmin, wvmax = np.min(self.spec.dispersion), np.max(self.spec.dispersion)
        gdlin = []
        for comp in all_comp:
            for line in comp.lines:
                wvobs = (1+line.attrib['z'])*line.wrest 
                if (wvobs>wvmin) & (wvobs<wvmax):
                    gdlin.append(line)
        # Voigt
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        self.model = xsv.voigt_model(self.spec.dispersion,gdlin, fwhm=self.fwhm)#,debug=True)
        
        #Define arrays for plotting residuals
        if self.plot_residuals:
            self.residual_limit = self.spec.sig * self.residual_normalization_factor
            self.residual = (self.spec.flux - self.model.flux) * self.residual_normalization_factor

    # Add a component
    def add_component(self, wrest, vlim=None, zcomp=None, no_fit_mask=False):
        '''Generate a component and fit with Gaussian
        Parameters:
        ------------
        no_fit_mask: bool, optional
          Skip fit + masking (mainly for reading in a previous file)
        '''
        # Center z and reset vmin/vmax
        if zcomp is None:
            zmin,zmax = self.z + (1+self.z)*(self.avmnx.value/3e5)
            zcomp = (zmin+zmax)/2.
        if vlim is None:
            vlim = self.avmnx - (self.avmnx[1]+self.avmnx[0])/2.
        new_comp = Component(zcomp, wrest,vlim=vlim,
            linelist=self.llist['ISM']) 
        # Fit
        #print('doing fit for {:g}'.format(wrest))
        if not no_fit_mask:
            self.fit_component(new_comp)

            # Mask for analysis
            for line in new_comp.lines:
                #print('masking {:g}'.format(line.wrest))
                wvmnx = line.wrest*(1+new_comp.zcomp)*(1 + vlim.value/3e5)
                gdp = np.where((self.spec.dispersion>wvmnx[0])&
                    (self.spec.dispersion<wvmnx[1]))[0]
                if len(gdp) > 0:
                    self.spec.mask[gdp] = 1

        # Add to component list and Fiddle
        if self.parent is not None:
            self.parent.comps_widg.add_component(new_comp)
            self.parent.fiddle_widg.init_component(new_comp)

        # Update model
        self.current_comp = new_comp
        if not no_fit_mask:
            self.update_model()

    def fit_component(self,component):
        '''Fit the component and save values'''
        from astropy.modeling import fitting
        # Generate Fit line
        fit_line = AbsLine(component.init_wrest,
            linelist=self.llist[self.llist['List']])
        fit_line.analy['vlim'] = component.vlim
        fit_line.analy['spec'] = self.spec
        fit_line.attrib['z'] = component.zcomp
        fit_line.measure_aodm()
        # Guesses
        fmin = np.argmin(self.spec.flux[fit_line.analy['pix']])
        zguess = self.spec.dispersion[fit_line.analy['pix'][fmin]]/component.init_wrest - 1.
        bguess = (component.vlim[1]-component.vlim[0])/2.
        Nguess = fit_line.attrib['logN']
        # Voigt model
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        fitvoigt = xsv.single_voigt_model(logN=Nguess,b=bguess.value,
                                z=zguess, wrest=component.init_wrest.value,
                                gamma=fit_line.data['gamma'].value, 
                                f=fit_line.data['f'], fwhm=self.fwhm)
        # Restrict z range
        fitvoigt.logN.min = 10.
        fitvoigt.b.min = 1.
        fitvoigt.z.min = component.zcomp+component.vlim[0].value/3e5/(1+component.zcomp)
        fitvoigt.z.max = component.zcomp+component.vlim[1].value/3e5/(1+component.zcomp)
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        
        # Fit
        fitter = fitting.LevMarLSQFitter()
        parm = fitter(fitvoigt,self.spec.dispersion[fit_line.analy['pix']],
            self.spec.flux[fit_line.analy['pix']].value)

        # Save and sync
        component.attrib['N'] = parm.logN.value
        component.attrib['z'] = parm.z.value
        component.attrib['b'] = parm.b.value * u.km/u.s
        component.sync_lines()

    def out_of_bounds(self,coord):
        '''Check for out of bounds
        '''
        # Default is x
        if ((coord < np.min(self.spec.dispersion)) 
            or (coord > np.max(self.spec.dispersion))):
            print('Out of bounds!')
            return True
        else:
            return False

    # Key stroke 
    def on_key(self,event):
        # Init
        rescale = True
        fig_clear = False
        wrest = None
        flg = 1
        sv_idx = self.idx_line

        ## Change rows/columns
        if event.key == 'k':
            self.sub_xy[0] = max(0, self.sub_xy[0]-1)
        if event.key == 'K':
            self.sub_xy[0] = self.sub_xy[0]+1
        if event.key == 'c':
            self.sub_xy[1] = max(1, self.sub_xy[1]-1)
        if event.key == 'C':
            self.sub_xy[1] = max(1, self.sub_xy[1]+1)
        if event.key == '(':
            if self.subxy_state == 'Out':
                self.sub_xy = [3,2]
                self.subxy_state = 'In'
            else:
                self.sub_xy = [5,3]
                self.subxy_state = 'Out'

        ## NAVIGATING
        if event.key in self.psdict['nav']: 
            flg = xxgu.navigate(self.psdict,event)
        if event.key == '-':
            self.idx_line = max(0, self.idx_line-self.sub_xy[0]*self.sub_xy[1]) # Min=0
            if self.idx_line == sv_idx:
                print('Edge of list')
        if event.key == '=':
            self.idx_line = min(len(self.llist['show_line'])-self.sub_xy[0]*self.sub_xy[1],
                                self.idx_line + self.sub_xy[0]*self.sub_xy[1]) 
            if self.idx_line == sv_idx:
                print('Edge of list')

        # Find line
        try:
            wrest = event.inaxes.get_gid()
        except AttributeError:
            return
        else:
            wvobs = wrest*(1+self.z)
            pass

        ## Fiddle with a Component
        if event.key in ['N','n','v','V','<','>','R']:
            if self.parent.fiddle_widg.component is None:
                print('Need to generate a component first!')
                return
            if event.key == 'N': 
                self.parent.fiddle_widg.component.attrib['N'] += 0.05
            elif event.key == 'n': 
                self.parent.fiddle_widg.component.attrib['N'] -= 0.05
            elif event.key == 'v': 
                self.parent.fiddle_widg.component.attrib['b'] -= 2*u.km/u.s
            elif event.key == 'V': 
                self.parent.fiddle_widg.component.attrib['b'] += 2*u.km/u.s
            elif event.key == '<': 
                self.parent.fiddle_widg.component.attrib['z'] -= 2e-5 #should be a fraction of pixel size
            elif event.key == '>': 
                self.parent.fiddle_widg.component.attrib['z'] += 2e-5
    
            elif event.key == 'R': # Refit
                self.fit_component(self.parent.fiddle_widg.component)
            # Updates (this captures them all and redraws)
            self.parent.fiddle_widg.update_component()
        ## Grab/Delete a component
        if event.key in ['D','S','d']:
            # Delete selected component
            if event.key == 'd': 
                self.parent.delete_component(self.parent.fiddle_widg.component)
                return
            #
            components = self.parent.comps_widg.all_comp
            iwrest = np.array([comp.init_wrest.value for comp in components])*u.AA
            mtc = np.where(wrest == iwrest)[0]
            if len(mtc) == 0:
                return
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            dvz = np.array([3e5*(self.z- components[mt].zcomp)/(1+self.z) for mt in mtc])
            # Find minimum
            mindvz = np.argmin(np.abs(dvz+event.xdata))
            if event.key == 'S':
                self.parent.fiddle_widg.init_component(components[mtc[mindvz]])
            elif event.key == 'D': # Delete nearest component to cursor
                self.parent.delete_component(components[mtc[mindvz]])

        ## Reset z
        if event.key == ' ': #space to move redshift
            #from xastropy.relativity import velocities
            #newz = velocities.z_from_v(self.z, event.xdata)
            self.z = self.z + event.xdata*(1+self.z)/3e5
            #self.abs_sys.zabs = newz
            # Drawing
            self.psdict['xmnx'] = self.vmnx.value
        if event.key == '^': 
            zgui = xxgu.AnsBox('Enter redshift:',float)
            zgui.exec_()
            self.z = zgui.value
            self.psdict['xmnx'] = self.vmnx.value

        # Choose line
        if event.key == "%":
            # GUI
            self.select_line_widg = xspw.SelectLineWidget(
                self.llist[self.llist['List']]._data)
            self.select_line_widg.exec_()
            line = self.select_line_widg.line
            if line.strip() == 'None':
                return
            #
            quant = line.split('::')[1].lstrip()
            spltw = quant.split(' ')
            wrest = Quantity(float(spltw[0]), unit=spltw[1])
            #
            self.z = (wvobs/wrest - 1.).value
            #self.statusBar().showMessage('z = {:f}'.format(z))
            self.init_lines()

        # Toggle line lists
        if event.key == 'H':
            self.llist['List'] = 'HI'
            self.init_lines()
        if event.key == 'U':
            self.parent.update_strongest_lines()
            self.init_lines()

        ## Velocity limits
        unit = u.km/u.s
        if event.key in ['1','2']:
            if event.key == '1': 
                self.vmin = event.xdata*unit
            if event.key == '2': 
                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()
                absline.analy['vlim'][1] = event.xdata*unit
            self.update_component()

        ## Add component
        if event.key == 'A': # Add to lines
            if self.out_of_bounds(wvobs*(1+event.xdata/3e5)):
                return
            if self.flag_add is False:
                self.vtmp = event.xdata
                self.flag_add = True
                self.wrest = wrest
            else:
                self.avmnx = np.array([np.minimum(self.vtmp,event.xdata), 
                    np.maximum(self.vtmp,event.xdata)])*unit
                self.add_component(wrest)
                # Reset
                self.flag_add = False
                self.wrest = 0.

        # Fiddle with analysis mask        
        if event.key in ['x','X']:  
            # x = Delete mask
            # X = Add to mask
            if self.flag_mask is False:
                self.wrest = wrest
                self.wtmp = wvobs*(1+event.xdata/3e5)
                self.vtmp = event.xdata
                self.flag_mask = True
            else:
                wtmp2 = wvobs*(1+event.xdata/3e5)
                twvmnx = [np.minimum(self.wtmp,wtmp2), np.maximum(self.wtmp,wtmp2)]
                # Modify mask
                mskp = np.where((self.spec.dispersion>twvmnx[0])&
                    (self.spec.dispersion<twvmnx[1]))[0]
                #print(twvmnx,len(mskp))
                if event.key == 'x':
                    self.spec.mask[mskp] = 0
                elif event.key == 'X':
                    self.spec.mask[mskp] = 1
                # Reset
                self.flag_mask = False
                self.wrest = 0.

        # Labels
        if event.key == 'L': # Toggle ID lines
            self.flag_idlbl = ~self.flag_idlbl
            
        # AODM plot
        if event.key == ':':  # 
            # Grab good lines
            from xastropy.xguis import spec_guis as xsgui
            gdl = [iline.wrest for iline in self.abs_sys.lines 
                if iline.analy['do_analysis'] > 0]
            # Launch AODM
            if len(gdl) > 0:
                gui = xsgui.XAODMGui(self.spec, self.z, gdl, vmnx=self.vmnx, norm=self.norm)
                gui.exec_()
            else:
                print('VelPlot.AODM: No good lines to plot')

            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()

        #if wrest is not None: # Single window
        #    flg = 3
        if event.key in ['c','C','k','K','W','!', '@', '=', '-', 'X', ' ','R']: # Redraw all
            flg = 1 
        if event.key in ['Y']:
            rescale = False
        if event.key in ['c','C','k','K', 'R', '(']:
            fig_clear = True

        if flg==1: # Default is to redraw
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
            # Title
            self.fig.suptitle('z={:.5f}'.format(self.z),fontsize='large')
            # Components
            components = self.parent.comps_widg.all_comp 
            iwrest = np.array([comp.init_wrest.value for comp in components])*u.AA
            # Loop on windows
            all_idx = self.llist['show_line']
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            # Labels
            if self.flag_idlbl:
                line_wvobs = []
                line_lbl = []
                for comp in components:
                    if comp.attrib['Quality'] == 'None':
                        la = ''
                    else: 
                        la = comp.attrib['Quality']
                    for line in comp.lines:
                        line_wvobs.append(line.wrest.value*(line.attrib['z']+1))
                        line_lbl.append(line.trans+',{:.3f}{:s}'.format(line.attrib['z'],la))
                line_wvobs = np.array(line_wvobs)*u.AA
                line_lbl = np.array(line_lbl)
            # Subplots
            nplt = self.sub_xy[0]*self.sub_xy[1]
            if len(all_idx) <= nplt:
                self.idx_line = 0
            subp = np.arange(nplt) + 1
            subp_idx = np.hstack(subp.reshape(self.sub_xy[0],self.sub_xy[1]).T)
            #print('idx_l={:d}, nplt={:d}, lall={:d}'.format(self.idx_line,nplt,len(all_idx)))
            
            #try different color per ion species
            color_model = '#999966'
            colors = ['#0066FF','#339933','#CC3300','#660066','#FF9900','#B20047']
            color_ind = 0

            #loop over individual velplot axes
            for jj in range(min(nplt, len(all_idx))):
                try:
                    idx = all_idx[jj+self.idx_line]
                except IndexError:
                    continue # Likely too few lines
                #print('jj={:d}, idx={:d}'.format(jj,idx))
                # Grab line
                wrest = self.llist[self.llist['List']].wrest[idx]
                kwrest = wrest.value # For the Dict
                
                #define colors for visually grouping same species
                if jj > 0:
                    name_aux = self.llist[self.llist['List']].name[idx].split(' ')[0]
                    name_aux2 = self.llist[self.llist['List']].name[idx-1].split(' ')[0]
                    if name_aux != name_aux2:
                        color_ind += 1
                color = colors[color_ind % len(colors)]

                # Single window?
                #if in_wrest is not None:
                #    if np.abs(wrest-in_wrest) > (1e-3*u.AA):
                #        continue
                # Generate plot
                self.ax = self.fig.add_subplot(self.sub_xy[0],self.sub_xy[1], subp_idx[jj])
                self.ax.clear()        

                # GID for referencing
                self.ax.set_gid(wrest)

                # Zero velocity line
                self.ax.plot( [0., 0.], [-1e9, 1e9], ':', color='gray')
                # Velocity
                wvobs = (1+self.z) * wrest
                wvmnx = wvobs*(1 + np.array(self.psdict['xmnx'])/3e5)
                velo = (self.spec.dispersion/wvobs - 1.)*const.c.to('km/s')
                
                # Plot
                self.ax.plot(velo, self.spec.flux, '-',color=color,drawstyle='steps-mid',lw=0.5)
                # Model
                self.ax.plot(velo, self.model.flux, '-',color=color_model,lw=0.5)

                #Error & residuals
                if self.plot_residuals:
                    self.ax.plot(velo, self.residual_limit, 'k-',drawstyle='steps-mid',lw=0.5)
                    self.ax.plot(velo, -self.residual_limit, 'k-',drawstyle='steps-mid',lw=0.5)
                    self.ax.plot(velo, self.residual, '.',color='grey',ms=2)

                #import pdb
                #pdb.set_trace()

                # Labels
                if (((jj+1) % self.sub_xy[0]) == 0) or ((jj+1) == len(all_idx)):
                    self.ax.set_xlabel('Relative Velocity (km/s)')
                else:
                    self.ax.get_xaxis().set_ticks([])
                lbl = self.llist[self.llist['List']].name[idx]
                self.ax.text(0.01, 0.15, lbl, color=color, transform=self.ax.transAxes,
                             size='x-small', ha='left',va='center',backgroundcolor='w',bbox={'pad':0,'edgecolor':'none'})
                if self.flag_idlbl:
                    # Any lines inside?
                    mtw = np.where((line_wvobs > wvmnx[0]) & (line_wvobs<wvmnx[1]))[0]
                    for imt in mtw:
                        v = 3e5*(line_wvobs[imt]/wvobs - 1)
                        self.ax.text(v, 0.5, line_lbl[imt], color=color_model,backgroundcolor='w',
                            bbox={'pad':0,'edgecolor':'none'}, size='xx-small', rotation=90.,ha='center',va='center')

                # Analysis regions
                if np.sum(self.spec.mask) > 0.:
                    gdp = self.spec.mask==1
                    if len(gdp) > 0:
                        self.ax.scatter(velo[gdp],self.spec.flux[gdp],
                            marker='o',color=color,s=3.,alpha=0.5)

                # Reset window limits
                self.ax.set_ylim(self.psdict['ymnx'])
                self.ax.set_xlim(self.psdict['xmnx'])

                # Add line?
                if self.wrest == wrest:
                    self.ax.plot([self.vtmp]*2,self.psdict['ymnx'], '--',
                        color='red')

                # Components
                mtc = np.where(wrest == iwrest)[0]
                if len(mtc) > 0:
                    #QtCore.pyqtRemoveInputHook()
                    #xdb.set_trace()
                    #QtCore.pyqtRestoreInputHook()
                    for mt in mtc:
                        comp = components[mt]
                        #QtCore.pyqtRemoveInputHook()
                        #xdb.set_trace()
                        #QtCore.pyqtRestoreInputHook()
                        dvz = const.c.to('km/s')*(self.z-comp.zcomp)/(1+self.z)
                        if dvz.value < np.max(np.abs(self.psdict['xmnx'])):
                            if comp is self.parent.fiddle_widg.component:
                                lw = 1.5
                            else:
                                lw = 1.
                            # Plot
                            for vlim in comp.vlim:
                                self.ax.plot([vlim.value-dvz.value]*2,self.psdict['ymnx'], 
                                    '--', color='r',linewidth=lw)
                            self.ax.plot([-1.*dvz.value]*2,[1.0,1.05],
                                '-', color='grey',linewidth=lw)

                # Fonts
                xputils.set_fontsize(self.ax,6.)
        # Draw
        self.canvas.draw()
############        
class FiddleComponentWidget(QtGui.QWidget):
    ''' Widget to fiddle with a given component
    '''
    def __init__(self, component=None, parent=None):
        '''
        '''
        super(FiddleComponentWidget, self).__init__(parent)

        self.parent = parent
        #if not status is None:
        #    self.statusBar = status
        self.label = QtGui.QLabel('Component:',self)
        self.zwidget = xxgu.EditBox(-1., 'zc=', '{:0.5f}')
        self.Nwidget = xxgu.EditBox(-1., 'Nc=', '{:0.2f}')
        self.bwidget = xxgu.EditBox(-1., 'bc=', '{:0.1f}')

        self.ddlbl = QtGui.QLabel('Quality')
        self.ddlist = QtGui.QComboBox(self)
        self.ddlist.addItem('None')
        self.ddlist.addItem('a')
        self.ddlist.addItem('b')
        self.ddlist.addItem('c')
        self.Cwidget = xxgu.EditBox('None', 'Comment=', '{:s}')

        # Init further
        if component is not None:
            self.init_component(component)
        else:
            self.component = component

        # Connect
        self.ddlist.activated[str].connect(self.setQuality)
        self.connect(self.Nwidget.box, 
            QtCore.SIGNAL('editingFinished ()'), self.setbzN)
        self.connect(self.zwidget.box, 
            QtCore.SIGNAL('editingFinished ()'), self.setbzN)
        self.connect(self.bwidget.box, 
            QtCore.SIGNAL('editingFinished ()'), self.setbzN)
        self.connect(self.Cwidget.box, 
            QtCore.SIGNAL('editingFinished ()'), self.setbzN)

        # Layout
        zNbwidg = QtGui.QWidget()
        hbox2 = QtGui.QHBoxLayout()
        hbox2.addWidget(self.zwidget)
        hbox2.addWidget(self.Nwidget)
        hbox2.addWidget(self.bwidget)
        zNbwidg.setLayout(hbox2)

        ddwidg = QtGui.QWidget()
        vbox1 = QtGui.QVBoxLayout()
        vbox1.addWidget(self.ddlbl)
        vbox1.addWidget(self.ddlist)
        ddwidg.setLayout(vbox1)

        commwidg = QtGui.QWidget()
        hbox3 = QtGui.QHBoxLayout()
        hbox3.addWidget(ddwidg)
        hbox3.addWidget(self.Cwidget)
        commwidg.setLayout(hbox3)

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.label)
        vbox.addWidget(zNbwidg)
        vbox.addWidget(commwidg)
        self.setLayout(vbox)

    def init_component(self,component):
        '''Setup Widget for the input component'''
        self.component = component
        # Values
        self.Nwidget.set_text(self.component.attrib['N'])
        self.zwidget.set_text(self.component.attrib['z'])
        self.bwidget.set_text(self.component.attrib['b'].value)
        self.Cwidget.set_text(self.component.comment)
        # Quality
        idx = self.ddlist.findText(self.component.attrib['Quality'])
        self.ddlist.setCurrentIndex(idx)
        # Label
        self.set_label()

    def setQuality(self,text):
        if self.component is not None:
            self.component.attrib['Quality'] = text

    def reset(self):
        #
        self.component = None
        #  Values
        self.Nwidget.set_text(-1.)
        self.zwidget.set_text(-1.)
        self.bwidget.set_text(-1.)
        self.Cwidget.set_text('None')
        idx = self.ddlist.findText('None')
        self.ddlist.setCurrentIndex(idx)
        # Label
        self.set_label()

    def update_component(self):
        '''Values have changed'''
        self.Nwidget.set_text(self.component.attrib['N'])
        self.zwidget.set_text(self.component.attrib['z'])
        self.bwidget.set_text(self.component.attrib['b'].value)
        self.Cwidget.set_text(self.component.comment)
        if self.parent is not None:
            self.parent.updated_component()

    def set_label(self):
        '''Sets the label for the Widget'''
        if self.component is not None:
            self.label.setText('Component: {:s}'.format(self.component.name))            
        else:
            self.label.setText('Component:')

    def setbzN(self):
        '''Set the component column density or redshift from the boxes'''
        if self.component is None:
            print('Need to generate a component first!')
        else:
            # Grab values
            self.component.attrib['N'] = (float(self.Nwidget.box.text()))
            self.component.attrib['z'] = (float(self.zwidget.box.text()))
            self.component.attrib['b'] = (float(self.bwidget.box.text()))*u.km/u.s
            self.component.comment = str(self.Cwidget.box.text())
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            # Update beyond
            if self.parent is not None:
                self.parent.updated_component()

# #####
class ComponentListWidget(QtGui.QWidget):
    ''' Widget to organize components on a sightline

    Parameters:
    -----------
    components: List
      List of components

    16-Dec-2014 by JXP
    '''
    def __init__(self, components, parent=None, no_buttons=False):
        '''
        only_one: bool, optional
          Restrict to one selection at a time? [False]
        no_buttons: bool, optional
          Eliminate Refine/Reload buttons?
        '''
        super(ComponentListWidget, self).__init__(parent)

        self.parent = parent

        #if not status is None:
        #    self.statusBar = status
        self.all_comp = components  # Actual components

        list_label = QtGui.QLabel('Components:')
        self.complist_widget = QtGui.QListWidget(self) 
        #self.complist_widget.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.complist_widget.addItem('None')
        #self.abslist_widget.addItem('Test')

        # Lists
        self.items = []     # Selected
        self.all_items = [] # Names

        self.complist_widget.setCurrentRow(0)
        self.complist_widget.itemSelectionChanged.connect(self.on_list_change)

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(list_label)
        vbox.addWidget(self.complist_widget)
        self.setLayout(vbox)

    # ##
    def on_list_change(self):
        '''
        Changed an item in the list
        '''
        item = self.complist_widget.selectedItems()
        txt = item[0].text()
        if txt == 'None':
            if self.parent is not None:
                self.parent.updated_compslist(None)
        else:
            ii = self.all_items.index(txt)
            if self.parent is not None:
                self.parent.updated_compslist(self.all_comp[ii])

        '''
        items = self.complist_widget.selectedItems()
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
        '''
    '''
    def selected_components(self):
        items = self.complist_widget.selectedItems()
        selc = []
        for item in items:
            txt = item.text()
            if txt == 'None':
                continue
            ii = self.all_items.index(txt)
            selc.append(self.all_comp[ii])
        # Return
        return selc
    '''

    def add_component(self,component):
        self.all_comp.append( component )
        self.add_item(component.name)

    def add_item(self,comp_name):
        #
        self.all_items.append(comp_name) 
        self.complist_widget.addItem(comp_name)
        self.complist_widget.item(len(self.all_items)).setSelected(True)

    def remove_item(self,comp_name):
        # Delete
        idx = self.all_items.index(comp_name)
        del self.all_items[idx]
        self.all_comp.pop(idx)
        tmp = self.complist_widget.takeItem(idx+1) # 1 for None
        #self.on_list_change()


class Component(object):
    def __init__(self, z, wrest, vlim=[-300.,300]*u.km/u.s,
        linelist=None):
        # Init
        self.init_wrest = wrest
        self.zcomp = z
        self.vlim = vlim
        self.attrib = {'N': 0., 'Nsig': 0., 'flagN': 0, # Column
                       'b': 0.*u.km/u.s, 'bsig': 0.*u.km/u.s,  # Doppler
                       'z': self.zcomp, 'zsig': 0.,
                       'Quality': 'None'}
        self.comment = 'None'
        #
        self.linelist = linelist
        self.lines = []
        self.init_lines()
        #
        self.name = 'z{:.5f}_{:s}'.format(
            self.zcomp,self.lines[0].data['name'].split(' ')[0])
        #
    def init_lines(self):
        '''Fill up the component lines
        '''
        if self.linelist is None:
            self.linelist = LineList('Strong')
        # Get the lines

        all_trans = self.linelist.all_transitions(self.init_wrest)
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        if isinstance(all_trans,dict):
            all_trans = [all_trans]
        for trans in all_trans:
            self.lines.append(AbsLine(trans['wrest'],
                linelist=self.linelist))

        # Sync
        self.sync_lines()

    def sync_lines(self):
        '''Synchronize attributes of the lines
        '''
        for line in self.lines:
            line.attrib['N'] = self.attrib['N']
            line.attrib['b'] = self.attrib['b']
            line.attrib['z'] = self.attrib['z']

# Script to run XSpec from the command line or ipython
def run_gui(*args, **kwargs):
    '''
    Runs the IGMGuessesGui

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
    parser.add_argument("-out_file", type=str, help="Output Guesses file")
    parser.add_argument("-fwhm", type=float, help="FWHM smoothing (pixels)")
    parser.add_argument("-previous_file", type=str, help="Input Guesses file")
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
        previous_file = pargs.previous_file
    except AttributeError:
        previous_file=None

    # Smoothing parameter
    try:
        fwhm = pargs.smooth
    except AttributeError:
        fwhm=None

    # zqso
    try:
        zqso = pargs.zqso
    except AttributeError:
        zqso=None
    
    app = QtGui.QApplication(sys.argv)
    gui = IGMGuessesGui(pargs.in_file, outfil=outfil, fwhm=fwhm,
        previous_file=previous_file, zqso=zqso)
    gui.show()
    app.exec_()

# ################
if __name__ == "__main__":
    import sys, os
    from linetools.spectra import io as lsi
    from xastropy.igm import abs_sys as xiabs

    if len(sys.argv) == 1: # TESTING

        flg_fig = 0 
        flg_fig += 2**0  # Fit LLS GUI
    
        # LLS
        if (flg_fig % 2**1) >= 2**0:
            #spec_fil = '/Users/xavier/Keck/ESI/RedData/PSS0133+0400/PSS0133+0400_f.fits'
            # spec_fil = os.getenv('DROPBOX_DIR')+'/Tejos_X/COS-Clusters/J1018+0546.txt'
            # spec_fil = os.getenv('DROPBOX_DIR')+'/Tejos_X/COS-Filaments/q1410.fits'
            spec_fil = os.getenv('DROPBOX_DIR')+'/Tejos_X/COS-Filaments/J1619+2543.fits'
            
            spec = lsi.readspec(spec_fil)
            spec.normalize()
            #spec.plot()
            #xdb.set_trace()
            # Launch
            app = QtGui.QApplication(sys.argv)
            app.setApplicationName('IGMGuesses')
            main = IGMGuessesGui(spec) 
            main.show()
            sys.exit(app.exec_())

    else: # RUN A GUI
        run_gui()
            