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
#;- NT: New version using linetools' AbsComponent
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Import libraries
import numpy as np
import warnings, imp

from PyQt4 import QtGui
from PyQt4 import QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
# Matplotlib Figure object
from matplotlib.figure import Figure

from astropy.units import Quantity
from astropy import units as u
from astropy import constants as const

from linetools.analysis import voigt as lav
from linetools.lists.linelist import LineList
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectralline import AbsLine
from linetools.isgm.abscomponent import AbsComponent
from linetools.guis import utils as ltgu
from linetools.guis import line_widgets as ltgl
from linetools.guis import simple_widgets as ltgsm

#from xastropy.atomic import ionization as xatomi
from xastropy.plotting import utils as xputils
from xastropy.xguis import spec_widgets as xspw
#from xastropy.xguis import utils as xxgu

from xastropy.xutils import xdebug as xdb

xa_path = imp.find_module('xastropy')[1]

c_mks = const.c.to('km/s')

#class IGMGuessesGui(QtGui.QMainWindow):

# GUI for fitting LLS in a spectrum
class IGMGuessesGui(QtGui.QMainWindow):
    ''' GUI to identify absorption features and provide reasonable
        first guesses of (z, logN, b) for subsequent Voigt profile
        fitting.

        v0.5
        30-Jul-2015 by JXP
    '''
    def __init__(self, ispec, parent=None, previous_file=None, 
        srch_id=True, outfil=None, fwhm=None, zqso=None,
        plot_residuals=True,n_max_tuple=None, min_strength=0.):
        QtGui.QMainWindow.__init__(self, parent)
        """
        ispec : str
            Name of the spectrum file to load
        previous_file: str, optional
            Name of the previous IGMguesses json file
        smooth: float, optional
            Number of pixels to smooth on
        zqso: float, optional
            Redshift of the quasar.  If input, a Telfer continuum is used
        plot_residuals : bool, optional
            Whether to plot residuals
        n_max_tuple : int, optional
            Maximum number of transitions per ion species to consider for plotting display.
        min_strength : float, optional
            Minimum strength for a transition to be considered in the analysis.
            The value should lie between (0,14.7), where 0. means 
            include everything, and 14.7 corresponds to the strength of 
            HI Lya transition assuming solar abundance.


        """
        # TODO
        # 1. Fix convolve window size
        # 2. Add COS LSF (?)

        self.help_message = """
i,o       : zoom in/out x limits
y         : zoom out y limits
Y         : guess y limits
t,b       : set y top/bottom limit
l,r       : set left/right x limit
[,]       : pan left/right
C,c       : add/remove column
K,k       : add/remove row
(         : toggle between many/few (15 or 6) panels per page
=,-       : move to next/previous page
f         : move to the first page
Space bar : set redshift from cursor position
^         : set redshift by hand
U         : update the main LineList at current redshift
H         : update to Lyman series LineList at current redshift
            (type `U` to get metals back)
A         : set limits for fitting an absorption component
            from cursor position (need to be pressed twice:
            once for left and once for right limit, respectively)
S         : select an absorption component from cursor position
D         : delete currently selected absorption component
d         : delete absorption component selected from component widget
N,n       : slightly increase/decrease column density in initial guess
V,v       : slightly increase/decrease b-value in initial guess
<,>       : slightly increase/decrease redshift in initial guess
R         : refit
X,x       : add/remove `good pixels` to keep for subsequent VP fitting
            (works as `A` command, i.e. need to define two limits)
L         : toggle between displaying/hiding labels of currently
            identified lines
%         : guess a transition and redshift for a given feature at
            the cursor's position
?         : print help message
"""

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
        self.n_max_tuple = n_max_tuple
        self.min_strength = min_strength

        # Load spectrum
        spec, spec_fil = ltgu.read_spec(ispec)
        # Normalize
        spec.normalize(co=spec.data[0]['co'])
        # make sure there are no nans in uncertainty, which affects the display of residuals
        spec.data[0]['sig'] = np.where(np.isnan(spec.data[0]['sig']), 0, spec.data[0]['sig'])

        # This attribute will store `good pixels` for subsequent Voigt Profile fitting
        spec.mask = np.zeros(len(spec.wavelength),dtype=int)

        # Full spectrum model
        self.model = XSpectrum1D.from_tuple(
            (spec.wavelength, np.ones(len(spec.wavelength))))

        # LineList (Grab ISM and HI as defaults)
        self.llist = ltgu.set_llist('ISM')
        self.llist['HI'] = LineList('HI')
        # self.llist['Strong'] = LineList('Strong')
        self.llist['Lists'].append('HI')
        self.llist['HI']._data = self.llist['HI']._data[::-1] #invert order of Lyman series
        #self.llist['show_line'] = np.arange(10) #maximum 10 to show for Lyman series
        
        # Define initial redshift
        z=0.0
        self.llist['z'] = z
        
        # Grab the pieces and tie together
        self.slines_widg = ltgl.SelectedLinesWidget(
            self.llist[self.llist['List']], parent=self, init_select='All')
        self.fiddle_widg = FiddleComponentWidget(parent=self)
        self.comps_widg = ComponentListWidget([], parent=self)
        self.velplot_widg = IGGVelPlotWidget(spec, z, 
            parent=self, llist=self.llist, fwhm=self.fwhm,plot_residuals=self.plot_residuals)
        self.wq_widg = ltgsm.WriteQuitWidget(parent=self)
        
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

        # Attempt to initialize
        self.update_strongest_lines()
        self.velplot_widg.init_lines()
        self.velplot_widg.on_draw(rescale=True, fig_clear=True)

        # Point MainWindow
        self.setCentralWidget(self.main_widget)

    def update_strongest_lines(self):
        '''Grab the strongest lines in the spectrum at the current
        redshift.
        '''
        z = self.velplot_widg.z
        wvmin = np.min(self.velplot_widg.spec.wavelength)
        wvmax = np.max(self.velplot_widg.spec.wavelength)
        wvlims = (wvmin/(1+z),wvmax/(1+z))
        transitions = self.llist['ISM'].available_transitions(
            wvlims,n_max=None, n_max_tuple=self.n_max_tuple,min_strength=self.min_strength)

        if transitions is not None:
            names = list(np.array(transitions['name']))
        else:
            names = ['HI 1215']
        self.llist['strongest'] = self.llist['strongest'].subset_lines(reset_data=True,subset=names)
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
            gdp = np.where((self.velplot_widg.spec.wavelength>wvmnx[0])&
                (self.velplot_widg.spec.wavelength<wvmnx[1]))[0]
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
            raise ValueError('Input FWHMs do not match. Please fix it!')
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
            self.velplot_widg.current_comp.attrib['b']= igmg_dict['cmps'][key]['bfit']*u.km/u.s
            self.velplot_widg.current_comp.attrib['logN']= igmg_dict['cmps'][key]['Nfit']
            self.velplot_widg.current_comp.attrib['Quality']= igmg_dict['cmps'][key]['Quality']
            self.velplot_widg.current_comp.comment = igmg_dict['cmps'][key]['Comment']
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
            out_dict['cmps'][key]['Nfit'] = comp.attrib['logN']
            out_dict['cmps'][key]['bfit'] = comp.attrib['b'].value
            out_dict['cmps'][key]['wrest'] = comp.init_wrest.value
            out_dict['cmps'][key]['vlim'] = list(comp.vlim.value)
            out_dict['cmps'][key]['Quality'] = str(comp.attrib['Quality'])
            out_dict['cmps'][key]['Comment'] = str(comp.comment)
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
    """ Widget for a velocity plot with interaction.
          Adapted from VelPlotWidget in spec_guis
        14-Aug-2015 by JXP
    """
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

        # init help message
        self.help_message = parent.help_message

        # Initialize
        self.parent = parent
        spec, spec_fil = ltgu.read_spec(ispec)
        
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
        self.model = XSpectrum1D.from_tuple(
            (spec.wavelength, np.ones(len(spec.wavelength))))

        self.plot_residuals = plot_residuals
        #Define arrays for plotting residuals
        if self.plot_residuals:
            self.residual_normalization_factor = 0.02/np.median(self.spec.sig)
            self.residual_limit = self.spec.sig * self.residual_normalization_factor
            self.residual = (self.spec.flux - self.model.flux) * self.residual_normalization_factor


        self.psdict = {} # Dict for spectra plotting
        self.psdict['x_minmax'] = self.vmnx.value # Too much pain to use units with this
        self.psdict['y_minmax'] = [-0.1, 1.1]
        self.psdict['nav'] = ltgu.navigate(0,0,init=True)

        

        # Status Bar?
        #if not status is None:
        #    self.statusBar = status

        # Line List
        if llist is None:
            self.llist = ltgu.set_llist(['HI 1215', 'HI 1025'])
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
        wvmin = np.min(self.spec.wavelength)
        wvmax = np.max(self.spec.wavelength)
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


    # Update model
    def update_model(self):
        if self.parent is None:
            return
        all_comp = self.parent.comps_widg.all_comp #selected_components()
        if len(all_comp) == 0:
            self.model.flux[:] = 1.
            return
        # Setup lines
        wvmin, wvmax = np.min(self.spec.wavelength), np.max(self.spec.wavelength)
        gdlin = []
        for comp in all_comp:
            for line in comp.lines:
                wvobs = (1+line.attrib['z'])*line.wrest 
                if (wvobs>wvmin) & (wvobs<wvmax):
                    line.attrib['N'] = 10.**line.attrib['logN'] / u.cm**2
                    gdlin.append(line)
        # Voigt
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        self.model = lav.voigt_from_abslines(self.spec.wavelength, gdlin, fwhm=self.fwhm)#,debug=True)
        
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

            # For Lyman series only mask pixels for fitting 
            # up to Ly-gamma; the rest should be done manually 
            # if wanted
            if new_comp.lines[0].name.startswith('HI '):
                aux_comp_list = new_comp.lines[::-1][:3] #invert order from ISM LineList and truncate
            else:
                aux_comp_list = new_comp.lines

            # Mask for analysis
            for line in aux_comp_list:
                #print('masking {:g}'.format(line.wrest))
                wvmnx = line.wrest*(1+new_comp.zcomp)*(1 + vlim.value/3e5)
                gdp = np.where((self.spec.wavelength>wvmnx[0])&
                    (self.spec.wavelength<wvmnx[1]))[0]
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

    def fit_component(self, component):
        '''Fit the component and save values'''
        from astropy.modeling import fitting
        # Generate Fit line
        fit_line = AbsLine(component.init_wrest,
            linelist=self.llist[self.llist['List']])
        fit_line.analy['vlim'] = component.vlim
        fit_line.analy['spec'] = self.spec
        fit_line.attrib['z'] = component.zcomp
        fit_line.measure_aodm(normalize=False)  # Already normalized
        # Guesses
        fmin = np.argmin(self.spec.flux[fit_line.analy['pix']])
        zguess = self.spec.wavelength[fit_line.analy['pix'][fmin]]/component.init_wrest - 1.
        bguess = (component.vlim[1]-component.vlim[0])/2.
        Nguess = np.log10(fit_line.attrib['N'].to('cm**-2').value)
        # Voigt model
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        fitvoigt = lav.single_voigt_model(logN=Nguess,b=bguess.value,
                                z=zguess, wrest=component.init_wrest.value,
                                gamma=fit_line.data['gamma'].value, 
                                f=fit_line.data['f'], fwhm=self.fwhm)
        # Restrict z range
        fitvoigt.logN.min = 10.
        fitvoigt.b.min = 1.
        fitvoigt.z.min = component.zcomp+component.vlim[0].value/3e5/(1+component.zcomp)
        fitvoigt.z.max = component.zcomp+component.vlim[1].value/3e5/(1+component.zcomp)

        # Fit
        fitter = fitting.LevMarLSQFitter()
        parm = fitter(fitvoigt,self.spec.wavelength[fit_line.analy['pix']],
            self.spec.flux[fit_line.analy['pix']].value)
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

        # Save and sync
        component.attrib['logN'] = parm.logN.value
        component.attrib['N'] = 10**parm.logN.value / u.cm**2
        component.attrib['z'] = parm.z.value
        component.attrib['b'] = parm.b.value * u.km/u.s
        component.sync_lines()

    def out_of_bounds(self,coord):
        '''Check for out of bounds
        '''
        # Default is x
        if ((coord < np.min(self.spec.wavelength))
            or (coord > np.max(self.spec.wavelength))):
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
            flg = ltgu.navigate(self.psdict,event)
        if event.key == '-':
            self.idx_line = max(0, self.idx_line-self.sub_xy[0]*self.sub_xy[1]) # Min=0
            if self.idx_line == sv_idx:
                print('Edge of list')
        if event.key == '=':
            self.idx_line = min(len(self.llist['show_line'])-self.sub_xy[0]*self.sub_xy[1],
                                self.idx_line + self.sub_xy[0]*self.sub_xy[1])
            if self.idx_line == sv_idx:
                print('Edge of list')
        if event.key == 'f':
            self.idx_line = 0
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
                self.parent.fiddle_widg.component.attrib['logN'] += 0.05
            elif event.key == 'n':
                self.parent.fiddle_widg.component.attrib['logN'] -= 0.05
            elif event.key == 'v':
                self.parent.fiddle_widg.component.attrib['b'] -= 5*u.km/u.s
            elif event.key == 'V':
                self.parent.fiddle_widg.component.attrib['b'] += 5*u.km/u.s
            elif event.key == '<':
                self.parent.fiddle_widg.component.attrib['z'] -= 4e-5 # should be a fraction of pixel size
            elif event.key == '>':
                self.parent.fiddle_widg.component.attrib['z'] += 4e-5

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
            self.psdict['x_minmax'] = self.vmnx.value
        if event.key == '^':
            zgui = ltgsm.AnsBox('Enter redshift:',float)
            zgui.exec_()
            self.z = zgui.value
            self.psdict['x_minmax'] = self.vmnx.value

        # Choose line
        if event.key == "%":
            # GUI
            self.select_line_widg = ltgl.SelectLineWidget(
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
                self.avmnx[0] = event.xdata*unit
            elif event.key == '2':
                self.avmnx[1] = event.xdata*unit
            # todo: we need to update the fit with new edges here

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

        # Fiddle with analysis pixel mask
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
                mskp = np.where((self.spec.wavelength>twvmnx[0])&
                    (self.spec.wavelength<twvmnx[1]))[0]
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

        if event.key == '?':
            print(self.help_message)

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
            self.ax.plot( [event.xdata,event.xdata], self.psdict['y_minmax'], ':', color='green')
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
                        line_lbl.append(line.name+',{:.3f}{:s}'.format(line.attrib['z'],la))
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
                wvmnx = wvobs*(1 + np.array(self.psdict['x_minmax'])/3e5)
                velo = (self.spec.wavelength/wvobs - 1.) * c_mks
                
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
                             size='x-small', ha='left',va='center',backgroundcolor='w',bbox={'pad':0,'edgecolor':'none',
                                                                                             'facecolor':'w'})
                if self.flag_idlbl:
                    # Any lines inside?
                    mtw = np.where((line_wvobs > wvmnx[0]) & (line_wvobs<wvmnx[1]))[0]
                    for imt in mtw:
                        v = 3e5*(line_wvobs[imt]/wvobs - 1)
                        self.ax.text(v, 0.5, line_lbl[imt], color=color_model,backgroundcolor='w',
                            bbox={'pad':0,'edgecolor':'none', 'facecolor':'w'}, size='xx-small', rotation=90.,ha='center',va='center')

                # Analysis regions
                if np.sum(self.spec.mask) > 0.:
                    gdp = self.spec.mask==1
                    if len(gdp) > 0:
                        self.ax.scatter(velo[gdp],self.spec.flux[gdp],
                            marker='o',color=color,s=3.,alpha=0.5)

                # Reset window limits
                self.ax.set_ylim(self.psdict['y_minmax'])
                self.ax.set_xlim(self.psdict['x_minmax'])

                # Add line?
                if self.wrest == wrest:
                    self.ax.plot([self.vtmp]*2,self.psdict['y_minmax'], '--',
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
                        dvz = c_mks * (self.z - comp.zcomp) / (1 + self.z)
                        if dvz.value < np.max(np.abs(self.psdict['x_minmax'])):
                            if comp is self.parent.fiddle_widg.component:
                                lw = 1.5
                            else:
                                lw = 1.
                            # Plot
                            for vlim in comp.vlim:
                                self.ax.plot([vlim.value-dvz.value]*2,self.psdict['y_minmax'],
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
        self.zwidget = ltgsm.EditBox(-1., 'zc=', '{:0.5f}')
        self.Nwidget = ltgsm.EditBox(-1., 'Nc=', '{:0.2f}')
        self.bwidget = ltgsm.EditBox(-1., 'bc=', '{:0.1f}')

        self.ddlbl = QtGui.QLabel('Quality')
        self.ddlist = QtGui.QComboBox(self)
        self.ddlist.addItem('None')
        self.ddlist.addItem('a')
        self.ddlist.addItem('b')
        self.ddlist.addItem('c')
        self.Cwidget = ltgsm.EditBox('None', 'Comment=', '{:s}')

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
        self.Nwidget.set_text(self.component.attrib['logN'])
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
        self.Nwidget.set_text(self.component.attrib['logN'])
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
            self.component.attrib['logN'] = (float(self.Nwidget.box.text()))
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


class Component(AbsComponent):
    def __init__(self, z, wrest, vlim=[-300.,300]*u.km/u.s,
        linelist=None):

        # Init
        self.init_wrest = wrest
        self.linelist = linelist
        self.lines = []
        self.init_lines()

        # Generate with type
        radec = (0*u.deg,0*u.deg)
        Zion = (self.lines[0].data['Z'],self.lines[0].data['ion'])
        Ej = self.lines[0].data['Ej']

        # Name for fine-structure
        if Ej.value > 0.:
            stars = '*'*(len(self.lines[0].name.split('*'))-1)
        else:
            stars = None
        AbsComponent.__init__(self,radec, Zion, z, vlim, Ej, comment='None', stars=stars)

        # Init cont.
        self.attrib = {'N': 0./u.cm**2, 'Nsig': 0./u.cm**2, 'flagN': 0,  # Column
                       'logN': 0., 'sig_logN': 0.,
                       'b': 0.*u.km/u.s, 'bsig': 0.*u.km/u.s,  # Doppler
                       'z': self.zcomp, 'zsig': 0.,
                       'Quality': 'None'}

        # Sync
        self.sync_lines()

        # Use different naming convention here
        self.name = 'z{:.5f}_{:s}'.format(
            self.zcomp,self.lines[0].data['name'].split(' ')[0])



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


    def sync_lines(self):
        '''Synchronize attributes of the lines
        '''
        for line in self.lines:
            line.attrib['logN'] = self.attrib['logN']
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

    parser = argparse.ArgumentParser(description='Parser for IGMGuesses')
    parser.add_argument("in_file", type=str, help="Spectral file")
    parser.add_argument("-out_file", type=str, help="Output Guesses file")
    parser.add_argument("-fwhm", type=float, help="FWHM smoothing (pixels)")
    parser.add_argument("-previous_file", type=str, help="Input Guesses file")
    parser.add_argument("-n_max_tuple", type=int, help="Maximum number of transitions per ion species to display")
    parser.add_argument("-min_strength", type=float, help="Minimum strength for transitions to be considered; choose values (0,14.7)")


    if len(args) == 0:
        pargs = parser.parse_args()
    else: # better know what you are doing!
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

    # n_max_tuple
    try:
        n_max_tuple = pargs.n_max_tuple
    except AttributeError:
        n_max_tuple=None

    # min_strength
    try:
        min_strength = pargs.min_strength
    except AttributeError:
        min_strength = 0.
    if min_strength is None:
        min_strength = 0.
    
    app = QtGui.QApplication(sys.argv)
    gui = IGMGuessesGui(pargs.in_file, outfil=outfil, fwhm=fwhm,
        previous_file=previous_file, zqso=zqso,n_max_tuple=n_max_tuple,min_strength=min_strength)
    gui.show()
    app.exec_()

# ################
if __name__ == "__main__":
    import sys, os
    from linetools.spectra import io as lsi

    if len(sys.argv) == 1: # TESTING

        xdb.set_trace()  # Do the line below
        # python igmguesses.py ~/Dropbox/BHB_abs/COSdata/visit01/tmp_nrm.fits

    else: # RUN A GUI
        run_gui()
            