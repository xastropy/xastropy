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
import warnings
import imp

from PyQt4 import QtGui
from PyQt4 import QtCore

# Matplotlib Figure object
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

from astropy.units import Quantity
from astropy import units as u

from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectra import convolve as lsc
import linetools.spectra.io as lsi
from linetools.spectralline import AbsLine

from pyigm.abssys.lls import LLSSystem
from pyigm.abssys import lls as igmlls

from xastropy.xutils import xdebug as xdb
from xastropy.xguis import spec_widgets as xspw
from xastropy.xguis import spec_guis as xspg
from xastropy.xguis import utils as xxgu
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

def set_llist(llist, in_dict=None, sort=True):
    ''' Method to set a line list dict for the Widgets
    sort: bool, optional [DEPRECATED CURRENTLY]
      Sort lines by rest wavelength [True]
    '''
    from linetools.lists.linelist import LineList
    from astropy.units.quantity import Quantity

    if in_dict is None:
        in_dict = dict(Lists=[])

    if isinstance(llist,basestring): # Set line list from a file
        in_dict['List'] = llist
        in_dict['Lists'].append(llist)
        if llist == 'None':
            in_dict['Plot'] = False
        else:
            in_dict['Plot'] = True
            # Load?
            if not (llist in in_dict):
                # Homebrew
                if llist == 'OVI':
                    gdlines = u.AA*[629.730, 702.332, 770.409, 780.324, 787.711, 832.927, 972.5367, 977.0201,
                        1025.7222, 1031.9261, 1037.6167, 1206.5, 1215.6700, 1260.4221]
                    llist_cls = LineList('Strong', subset=gdlines)
                    in_dict[llist] = llist_cls
                else:
                    llist_cls = LineList(llist)
                    # Sort
                    llist_cls._data.sort('wrest')
                    # Load
                    in_dict[llist] = llist_cls
    elif isinstance(llist,(Quantity,list)): # Set from a list of wrest
        in_dict['List'] = 'input.lst'
        in_dict['Lists'].append('input.lst')
        in_dict['Plot'] = True
        # Fill
        if sort:
            llist.sort()
        llist_cls = LineList('ISM', subset=llist)
        in_dict['input.lst'] = llist_cls
    else:
        raise IOError('Not ready for this type of input')

    # Return
    return in_dict

def navigate(psdict,event,init=False):
    ''' Method to Navigate spectrum
    init:  (False) Initialize
      Just pass back valid key strokes
    '''
    # Initalize
    if init is True:
        return ['l','r','b','t','T','i','I', 'o','O',
        '[',']','W','Z', 'Y', '{', '}', 's']

    #
    if (not isinstance(event.xdata,float)) or (not isinstance(event.ydata,float)):
        print('Navigate: You entered the {:s} key out of bounds'.format(
            event.key))
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
    elif event.key == 's':  # Select window (i.e. zoom-in)
        if psdict['tmp_xy'] is None:
            psdict['tmp_xy'] = [event.xdata,event.ydata]
            print('Press another s to set the zoom-in window')
        else:
            psdict['xmnx'][0] = np.minimum(event.xdata,psdict['tmp_xy'][0])
            psdict['xmnx'][1] = np.maximum(event.xdata,psdict['tmp_xy'][0])
            psdict['ymnx'][0] = np.minimum(event.ydata,psdict['tmp_xy'][1])
            psdict['ymnx'][1] = np.maximum(event.ydata,psdict['tmp_xy'][1])
            psdict['tmp_xy'] = None
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
        psdict['xmnx'] = copy.deepcopy(psdict['sv_xy'][0])
        psdict['ymnx'] = copy.deepcopy(psdict['sv_xy'][1])
    elif event.key == 'Z': # Zero
        psdict['ymnx'][0] = 0.
    else:
        if not (event.key in ['shift']):
            rstr = 'Key {:s} not supported.'.format(event.key)
            print(rstr)
        return 0
    return 1


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



class ExamineSpecWidget(QtGui.QWidget):
    ''' Widget to plot a spectrum and interactively
        fiddle about.  Akin to XIDL/x_specplot.pro

        12-Dec-2014 by JXP
        Parameters:
        ------------
        key_events: bool, optional
          Use key events? [True]
          Useful when coupling to other widgets
    '''
    def __init__(self, ispec, parent=None, status=None, llist=None,
                 abs_sys=None, norm=True, second_file=None, zsys=None,
                 key_events=True, vlines=None, plotzero=False, exten=None):
        '''
        spec = Spectrum1D
        '''
        super(ExamineSpecWidget, self).__init__(parent)

        # Spectrum
        if isinstance(ispec, XSpectrum1D):
            spec = ispec
            spec_fil = spec.filename
        else:
            # this is broken
            spec, spec_fil = xxgu.read_spec(ispec)

        self.orig_spec = spec # For smoothing
        self.spec = self.orig_spec

        self.vlines = []
        if vlines is not None:
            self.vlines.extend(vlines)

        self.plotzero = plotzero

        # Other bits (modified by other widgets)
        self.continuum = None
        self.model = None
        self.bad_model = None  # Discrepant pixels in model
        self.use_event = 1

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
            self.llist = {'Plot': False, 'List': 'None', 'z': 0., 'Lists': []}
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
        if key_events:
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
        xmin = np.min(self.spec.dispersion.value)
        xmax = np.max(self.spec.dispersion.value)
        from linetools.spectra.plotting import get_flux_plotrange
        ymin, ymax = get_flux_plotrange(self.spec.flux.value)#, mult_pos=1)
        # ymed = np.median(self.spec.flux.value)
        # ymin = 0. - 0.1*ymed
        # ymax = ymed * 1.5
        #
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        self.psdict['xmnx'] = np.array([xmin,xmax])
        self.psdict['ymnx'] = [ymin,ymax]
        self.psdict['sv_xy'] = [ [xmin,xmax], [ymin,ymax] ]
        self.psdict['tmp_xy'] = None
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
            lya_line = AbsLine(1215.6701*u.AA)
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
                iwv = np.array(sorted([self.adict['wv_1'], self.adict['wv_2']])) * self.spec.wcs.unit
                ic = np.array(sorted([self.adict['C_1'], self.adict['C_2']]))

                # Calculate the continuum (linear fit)
                param = np.polyfit(iwv, ic, 1)
                cfunc = np.poly1d(param)
                self.spec.conti = cfunc(self.spec.dispersion)

                if event.key == '$': # Simple stats
                    pix = self.spec.pix_minmax(iwv)[0]
                    mean = np.mean(self.spec.flux[pix])
                    median = np.median(self.spec.flux[pix])
                    stdv = np.std(self.spec.flux[pix]-self.spec.conti[pix])
                    S2N = median / stdv
                    mssg = 'Mean={:g}, Median={:g}, S/N={:g}'.format(mean,median,S2N)
                else:
                    # Find the spectral line (or request it!)
                    rng_wrest = iwv / (self.llist['z']+1)
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
                    aline = AbsLine(wrest,linelist=self.llist[self.llist['List']])
                    aline.attrib['z'] = self.llist['z']
                    aline.analy['spec'] = self.spec

                    # AODM
                    if event.key == 'N':
                        # Calculate the velocity limits and load-up
                        aline.analy['vlim'] = const.c.to('km/s') * (
                            ( iwv/(1+self.llist['z']) - wrest) / wrest )

                        # AODM
                        #QtCore.pyqtRemoveInputHook()
                        #xdb.set_trace()
                        #QtCore.pyqtRestoreInputHook()
                        aline.measure_aodm()
                        mssg = 'Using '+ aline.__repr__()
                        mssg = mssg + ' ::  logN = {:g} +/- {:g}'.format(
                            aline.attrib['logN'], aline.attrib['sig_logN'])
                    elif event.key == 'E':  #EW
                        aline.analy['wvlim'] = iwv
                        aline.measure_restew()
                        mssg = 'Using '+ aline.__repr__()
                        mssg = mssg + ' ::  Rest EW = {:g} +/- {:g}'.format(
                            aline.attrib['EW'].to(mAA), aline.attrib['sig_EW'].to(mAA))
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
                if self.llist['List'] == 'None':
                    print('Need to set a line list first!!')
                    self.vplt_flg = -1 # Nothing to do here
                    return
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
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            gui = xsgui.XVelPltGui(self.spec, z=z, outfil=outfil, llist=self.llist,
                                   abs_sys=ini_abs_sys, norm=self.norm,
                                   sel_wv=self.xval*self.spec.wcs.unit)
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
    def on_draw(self, replot=True, no_draw=False):
        """ Redraws the spectrum
        no_draw: bool, optional
          Draw the screen on the canvas?
        """
        #
        if replot is True:
            self.ax.clear()
            self.ax.plot(self.spec.dispersion, self.spec.flux,
                'k-',drawstyle='steps-mid')
            try:
                self.ax.plot(self.spec.dispersion, self.spec.sig, 'r:')
            except ValueError:
                pass
            self.ax.set_xlabel('Wavelength')
            self.ax.set_ylabel('Flux')

            # Continuum?
            if self.continuum is not None:
                self.ax.plot(self.continuum.dispersion, self.continuum.flux,
                    color='purple')

            # Model?
            if self.model is not None:
                self.ax.plot(self.model.dispersion, self.model.flux,
                    color='cyan')
                if self.bad_model is not None:
                    self.ax.scatter(self.model.dispersion[self.bad_model],
                        self.model.flux[self.bad_model],  marker='o',
                        color='red', s=3.)


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
                ii=-1
                for abs_sys in self.abs_sys:
                    ii+=1
                    lines = abs_sys.list_of_abslines()
                    #QtCore.pyqtRemoveInputHook()
                    #xdb.set_trace()
                    #QtCore.pyqtRestoreInputHook()
                    wrest = Quantity([line.wrest for line in lines])
                    wvobs = wrest * (abs_sys.zabs+1)
                    gdwv = np.where( ((wvobs.value+5) > self.psdict['xmnx'][0]) &  # Buffer for region
                                    ((wvobs.value-5) < self.psdict['xmnx'][1]))[0]
                    #for kk in range(len(gdwv)):
                    for jj in gdwv:
                        if lines[jj].analy['do_analysis'] == 0:
                            continue
                        # Paint spectrum red
                        wvlim = wvobs[jj]*(1 + lines[jj].analy['vlim']/const.c.to('km/s'))
                        pix = np.where( (self.spec.dispersion > wvlim[0]) & (self.spec.dispersion < wvlim[1]))[0]
                        self.ax.plot(self.spec.dispersion[pix], self.spec.flux[pix], '-',drawstyle='steps-mid',
                                     color=clrs[ii])
                        # Label
                        lbl = lines[jj].analy['name']+' z={:g}'.format(abs_sys.zabs)
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
                self.ax.plot(self.spec.dispersion,
                    self.lya_line.flux, color='green')

        # Reset window limits
        self.ax.set_xlim(self.psdict['xmnx'])
        self.ax.set_ylim(self.psdict['ymnx'])

        if self.plotzero:
            self.ax.axhline(0, lw=0.3, color='k')

        for line in self.vlines:
            self.ax.axvline(line, color='k', ls=':')

        # Draw
        if not no_draw:
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

# GUI for fitting LLS in a spectrum
class XFitLLSGUI(QtGui.QMainWindow):
    ''' GUI to fit LLS in a given spectrum
        v1.2
        30-Jul-2015 by JXP
    '''
    def __init__(self, ispec, parent=None, lls_fit_file=None,
        outfil=None, smooth=3., zqso=None, fN_gamma=None, template=None,
        dw=0.1, skip_wveval=False):
        QtGui.QMainWindow.__init__(self, parent)
        '''
        ispec : Spectrum1D or specfil
        lls_fit_file: str, optional
          Name of the LLS fit file to input
        smooth : float, optional
          Number of pixels to smooth on (FWHM)
        zqso : float, optional
          Redshift of the quasar.  If input, a Telfer continuum is used
        fN_gamma : float, optional
          Redshift evolution of f(N) or IGM fiddled continuum
        template : str, optional
          Filename of a QSO template to use instead of the Telfer
          continuum. Only used if zqso is also given.
        dw : float, optional
          Pixel width in Angstroms for the wavelength array used to
          generate optical depths. Default is 0.1.
        skip_wveval : bool, optional
          Skip rebinning of wavelengths in the Voigt profile generation.
          This can speed up the code considerably, but use it wisely.
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
        self.smooth = None
        self.base_continuum = None
        self.all_forest = []
        self.flag_write = False
        self.dw = float(dw)
        self.skip_wveval = skip_wveval
        if skip_wveval:
            warnings.warn("Skipping wavelength rebinning in Voigt.")
            warnings.warn("Make sure you know what you are doing!")

        # Spectrum
        if isinstance(ispec, XSpectrum1D):
            spec = ispec
            spec_fil = spec.filename
        else:
            # this is broken
            spec, spec_fil = xxgu.read_spec(ispec)


        # LineList
        self.llist = set_llist('Strong')
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

        vlines = [(912 * (1 + zqso) if zqso is not None else None)]
        self.spec_widg = ExamineSpecWidget(spec,status=self.statusBar,
                                           llist=self.llist, key_events=False,
                                           abs_sys=self.abssys_widg.abs_sys,
                                           vlines=vlines, plotzero=1)
        # Initialize continuum (and LLS if from a file)
        if lls_fit_file is not None:
            self.init_LLS(lls_fit_file,spec)
        else:
            self.conti_dict = xspc.init_conti_dict(
                Norm=float(np.median(spec.flux.value)),
                piv_wv=1215.*(1+zqso),
                #piv_wv2=915.*(1+zqso),
                igm='True')
        if self.base_continuum is None:
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
                #self.conti_dict['piv_wv'] = 1215.*(1+zqso)
                #self.conti_dict['piv_wv2'] = 915.*(1+zqso)
            else:
                self.zqso = None
                self.continuum = XSpectrum1D.from_tuple((
                    spec.dispersion,np.ones(len(spec.dispersion))))
            self.base_continuum = self.continuum.flux
        self.update_conti()

        self.spec_widg.continuum = self.continuum

        # Full Model (LLS+continuum)
        self.full_model = XSpectrum1D.from_tuple((
            spec.dispersion,np.ones(len(spec.dispersion))))
        if self.smooth is None:
            self.smooth = smooth

        # Initialize as needed
        if lls_fit_file is not None:
            self.update_boxes()
            self.update_model()

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

        #self.spec_widg.setFixedWidth(900)
        self.spec_widg.setMinimumWidth(900)

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
            iline.attrib['N'] = 10**self.abssys_widg.all_abssys[idx].NHI * u.cm**-2
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
        """Update continuum
        """
        cflux = self.base_continuum * self.conti_dict['Norm']
        # Double tilt
        if 'piv_wv2' in self.conti_dict.keys():
            lowwv = self.continuum.dispersion.value < self.conti_dict['piv_wv2']
            self.continuum.flux[lowwv] = (cflux[lowwv] * (self.continuum.dispersion.value[lowwv]/
                                            self.conti_dict['piv_wv2'])**self.conti_dict['tilt2'])
            self.continuum.flux[~lowwv] = (cflux[~lowwv] * (self.continuum.dispersion.value[~lowwv]/
                                            self.conti_dict['piv_wv'])**self.conti_dict['tilt'])
        else:
            self.continuum.flux = (cflux * (self.continuum.dispersion.value/
                    self.conti_dict['piv_wv'])**self.conti_dict['tilt'])
        if self.lls_model is not None:
            self.full_model.flux = self.lls_model * self.continuum.flux

    def update_model(self):
        '''Update absorption model '''
        from linetools.analysis import voigt as lav

        if len(self.abssys_widg.all_abssys) == 0:
            self.lls_model = None
            self.spec_widg.model = None
            return
        # use finer wavelength array to resolve absorption features.
        wa = self.full_model.dispersion
        # Angstroms
        # should really make this a constant velocity width array instead.
        if not self.skip_wveval:
            wa1 = np.arange(wa[0].value, wa[-1].value, self.dw) * wa.unit
        else:
            wa1 = wa
        all_tau_model = igmlls.tau_multi_lls(wa1,
           self.abssys_widg.all_abssys, skip_wveval=self.skip_wveval)
        #QtCore.pyqtRemoveInputHook()
        #import pdb; pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

        # Loop on forest lines
        for forest in self.all_forest:
            tau_Lyman = lav.voigt_from_abslines(wa1, forest.lines,
                ret='tau', skip_wveval=self.skip_wveval)
            all_tau_model += tau_Lyman

        # Flux and smooth
        flux = np.exp(-1. * all_tau_model)
        if self.smooth > 0:
            if not self.skip_wveval:
                mult = np.median(np.diff(wa.value)) / self.dw
                flux = lsc.convolve_psf(flux, self.smooth * mult)
            else:
                flux = lsc.convolve_psf(flux, self.smooth)
        if not self.skip_wveval:
            self.lls_model = np.interp(wa.value, wa1.value, flux)
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
        if event.key in ['C','1','2','!','@']: # Set continuum level
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
            elif event.key == '!':
                self.conti_dict['tilt2'] += 0.1
            elif event.key == '@':
                self.conti_dict['tilt2'] -= 0.1
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
            elif event.key == '$': # Toggle metal-lines
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
                    iline.attrib['N'] = 10**self.abssys_widg.all_abssys[idx].NHI * u.cm**-2
                    iline.attrib['b'] = self.abssys_widg.all_abssys[idx].bval
            # Update the model
            self.update_model()
        elif event.key in ['6','7','8','9']: # Add forest line
            self.add_forest(event.key,event.xdata/1215.6701 - 1.)
            self.update_model()
        elif event.key == 'Q':
            # set an attribute to tell calling script to abort
            print("Setting the quit attribute, the calling script should "
                  "abort after you close the GUI")
            self.script_quit = True
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
        forest = GenericAbsSystem((0.*u.deg,0.*u.deg), z, [-300.,300.]*u.km/u.s)
        # NHI
        NHI_dict = {'6':12.,'7':13.,'8':14.,'9':15.}
        forest.NHI=NHI_dict[inp]
        # Lines
        for name in ['HI 1215','HI 1025', 'HI 972']:
            aline = AbsLine(name,
                linelist=self.llist[self.llist['List']])
            # Attributes
            aline.attrib['N'] = 10**forest.NHI * u.cm**-2
            aline.attrib['b'] = 20.*u.km/u.s
            aline.attrib['z'] = forest.zabs
            # Append
            forest.lines.append(aline)
        # Append to forest lines
        self.all_forest.append(forest)

    def add_LLS(self,z, NHI=17.3,bval=20.*u.km/u.s, comment='None', model=True):
        """Generate a new LLS
        """
        #
        new_sys = LLSSystem((0*u.deg,0*u.deg),z,[-300.,300]*u.km/u.s,NHI=NHI)
        new_sys.bval = bval # This is not standard, but for convenience
        new_sys.comment = comment
        new_sys.fill_lls_lines(bval=bval, do_analysis=0)
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
        if model:  # For dealing with initialization
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
        NHI = 17.29 + np.log10(-1.*np.log(y/conti.flux.value[ximn]))
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

        #print('NHI={:g}'.format(NHI))
        z = x/(911.7)-1
        plls = LLSSystem((0*u.deg,0*u.deg),z,[-300.,300]*u.km/u.s,NHI=NHI)
        plls.bval = 20*u.km/u.s
        plls.fill_lls_lines(bval=20*u.km/u.s, do_analysis=0)

        # wrest, Tau model, flux
        wrest = spec.dispersion/(1+plls.zabs)
        tau = igmlls.tau_multi_lls(spec.dispersion,[plls])
        emtau = np.exp(-1. * tau)
        lls_flux = lsc.convolve_psf(emtau, 3.)
#xdb.xplot(wrest, lls_flux)

        # zmin (next highest LLS or zem)
        if len(self.abssys_widg.all_abssys) != 0:
            zlls = [lls.zabs for lls in self.abssys_widg.all_abssys if lls.zabs > plls.zabs]
            if len(zlls) == 0:
                zmin = self.zqso+0.01
            else:
                zmin = np.min(np.array(zlls)) - 0.01
        else:
            zmin = self.zqso+0.01

        # Pixels for analysis and rolling
        # NEED TO CUT ON X-Shooter ARM
        apix = np.where( (wrest > 914*u.AA) & #(spec.dispersion<5600*u.AA) &
                        (spec.dispersion<(1+zmin)*1026.*u.AA))[0] # Might go to Lyb
        nroll = (np.argmin(np.abs(spec.dispersion-(911.7*u.AA*(1+zmin))))- # Extra 0.01 for bad z
                   np.argmin(np.abs(spec.dispersion-(911.7*u.AA*(1+plls.zabs)))))
        # Require nroll does not exceed length of spectrum
        if np.max(apix)+nroll > len(spec.dispersion):
            nroll = len(spec.dispersion) - np.max(apix) - 1
        gdpix = np.arange(np.min(apix)-nroll,np.max(apix)+nroll+1)
        roll_flux = np.concatenate([np.ones(nroll),lls_flux[apix], np.ones(nroll)])
        roll_msk = roll_flux < 0.7

        # Generate data arrays
        wave_pad = spec.dispersion[gdpix]
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
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
        gui = xspg.XVelPltGui(self.spec_widg.spec, outfil=iabs_sys.absid_file,
                               abs_sys=iabs_sys, norm=self.spec_widg.norm)
        gui.exec_()

    # Read from a JSON file
    def init_LLS(self,fit_file,spec):
        import json
        # Read the JSON file
        with open(fit_file) as data_file:
            lls_dict = json.load(data_file)
        # Init continuum
        try:
            self.conti_dict = lls_dict['conti_model']
        except KeyError: # Historic
            self.conti_dict = lls_dict['conti']
        else:
            try:
                self.base_continuum = Quantity(lls_dict['conti'])
            except:
                print('Will generate a new base continuum')
                self.base_continuum = None
            else:
                self.continuum = XSpectrum1D.from_tuple((
                    spec.dispersion,np.ones(len(spec.dispersion))))
        #self.update_conti()
        # Check spectra names
        if spec.filename != lls_dict['spec_file']:
            warnings.warn('Spec file names do not match!')
        # LLS
        for key in lls_dict['LLS'].keys():
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            self.add_LLS(lls_dict['LLS'][key]['z'],
                NHI=lls_dict['LLS'][key]['NHI'],
                bval=lls_dict['LLS'][key]['bval']*u.km/u.s,
                comment=lls_dict['LLS'][key]['comment'], model=False)
        self.smooth = lls_dict['smooth']
        try:
            self.zqso = lls_dict['zqso']
        except KeyError:
            self.zqso = None
        # Updates
        #self.update_boxes()
        #self.update_model()
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

    # Write
    def write_out(self):
        import json, io
        # Create dict
        out_dict = dict(LLS={},conti_model=self.conti_dict, conti=list(self.base_continuum.value),
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
