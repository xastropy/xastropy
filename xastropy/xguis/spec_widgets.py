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
# Matplotlib Figure object
from matplotlib.figure import Figure

from astropy import constants as const
from astropy import units as u
from astropy.units import Quantity
u.def_unit(['mAA', 'milliAngstrom'], 0.001 * u.AA, namespace=globals()) # mA

from linetools.spectra import io as lsi
from linetools.spectralline import AbsLine
from linetools.lists.linelist import LineList
from linetools.guis import utils as ltgu
from linetools import utils as ltu
from linetools.isgm.abssystem import GenericAbsSystem

from xastropy import stats as xstats
from xastropy.xutils import xdebug as xdb
from xastropy.plotting import utils as xputils
from xastropy.igm.abs_sys import abssys_utils as xiaa
from pyigm.abssys.lls import LLSSystem
from xastropy.xguis import utils as xguiu

xa_path = imp.find_module('xastropy')[1]

# class ExamineSpecWidget
# class PlotLinesWidget
# class SelectLineWidget
# class SelectedLinesWidget
# class AbsSysWidget
# class VelPlotWidget
# class AODMWidget




# #####


# #####
class AbsSysWidget(QtGui.QWidget):
    ''' Widget to organize AbsSys along a given sightline

    Parameters:
    -----------
    abssys_list: List
      String list of abssys files

    16-Dec-2014 by JXP
    '''
    def __init__(self, abssys_list, parent=None, 
        only_one=False, linelist=None, no_buttons=False):
        '''
        only_one: bool, optional
          Restrict to one selection at a time? [False]
        no_buttons: bool, optional
          Eliminate Refine/Reload buttons?
        '''
        super(AbsSysWidget, self).__init__(parent)

        #if not status is None:
        #    self.statusBar = status
        self.abssys_list = abssys_list

        # Speeds things up
        if linelist is None:
            self.linelist = LineList('ISM')
        else:
            self.linelist = linelist

        # Create the line list 
        list_label = QtGui.QLabel('Abs Systems:')
        self.abslist_widget = QtGui.QListWidget(self) 
        if not only_one:
            self.abslist_widget.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.abslist_widget.addItem('None')
        #self.abslist_widget.addItem('Test')

        # Lists
        self.abs_sys = []
        self.items = []
        self.all_items = []
        self.all_abssys = []
        for abssys_fil in self.abssys_list:
            self.all_abssys.append(LLSSystem.from_absid_fil(abssys_fil,
                linelist=self.linelist))
            self.add_item(abssys_fil)

        self.abslist_widget.setCurrentRow(0)
        self.abslist_widget.itemSelectionChanged.connect(self.on_list_change)

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(list_label)

        # Buttons
        if not no_buttons:
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
        if ipos1 == -1:
            ipos1 = len(abssys_fil)
        #
        self.all_items.append( abssys_fil[ipos0:ipos1] )
        self.abslist_widget.addItem(abssys_fil[ipos0:ipos1] )

    def remove_item(self,idx):
        # Delete
        del self.all_items[idx]
        del self.all_abssys[idx]
        tmp = self.abslist_widget.takeItem(idx+1) # 1 for None
        self.on_list_change()

    def reload(self):
        print('AbsSysWidget: Reloading systems..')
        self.all_abssys = []
        for abssys_fil in self.abssys_list:
            self.all_abssys.append(LLSSystem.from_absid_fil(abssys_fil,
                linelist=self.linelist))
            #self.add_item(abssys_fil)
        self.on_list_change()

# ######################
class VelPlotWidget(QtGui.QWidget):
    ''' Widget for a velocity plot with interaction.

        19-Dec-2014 by JXP
    '''
    def __init__(self, ispec, z=None, parent=None, llist=None, norm=True,
                 vmnx=[-300., 300.]*u.km/u.s, abs_sys=None):
        '''
        spec = Spectrum1D
        Norm: Bool (False)
          Normalized spectrum?
        abs_sys: AbsSystem
          Absorption system class
        '''
        super(VelPlotWidget, self).__init__(parent)

        # Initialize
        spec, spec_fil = ltgu.read_spec(ispec)
        
        self.spec = spec
        self.spec_fil = spec_fil
        self.z = z
        self.vmnx = vmnx
        self.norm = norm

        # Abs_System 
        self.abs_sys = abs_sys
        if self.abs_sys is None:
            self.abs_sys = GenericAbsSystem((0.*u.deg,0.*u.deg), self.z, self.vmnx)
            self.abs_lines = []
        else:
            self.z = self.abs_sys.zabs
            # Line list
            if llist is None:
                self.abs_lines = self.abs_sys.list_of_abslines()
                if len(self.abs_lines)>0:
                    lwrest = [iline.wrest for iline in self.abs_lines]
                else:
                    lwrest = None
                if lwrest is not None:
                    llist = ltgu.set_llist(lwrest) # Not sure this is working..

        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

        self.psdict = {} # Dict for spectra plotting
        self.psdict['xmnx'] = self.vmnx.value # Too much pain to use units with this
        self.psdict['ymnx'] = [-0.1, 1.1]
        self.psdict['nav'] = ltgu.navigate(0,0,init=True)

        # Status Bar?
        #if not status is None:
        #    self.statusBar = status

        # Line List
        if llist is None:
            self.llist = ltgu.set_llist('Strong')
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

        # Update/generate lines [will not update] 
        for idx in gdlin:
            self.generate_line((self.z,wrest[idx])) 

    def grab_line(self, wrest):
        """ Grab a line from the list
        Parameters
        ----------
        wrest

        Returns
        -------
        iline : AbsLine object
        """
        awrest = [iline.wrest for iline in self.abs_lines]
        try:
            idx = awrest.index(wrest)
        except ValueError:
            return None
        else:
            return self.abs_lines[idx]

    def generate_line(self, inp):
        ''' Generate a new line, if it doesn't exist
        Parameters:
        ----------
        inp: tuple
          (z,wrest)
        '''
        # Generate?
        if self.grab_line(inp[1]) is None:
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            newline = AbsLine(inp[1],linelist=self.llist[self.llist['List']])
            print('VelPlot: Generating line {:g}'.format(inp[1]))
            newline.analy['vlim'] = self.vmnx/2.
            newline.attrib['z'] = self.abs_sys.zabs
            newline.analy['do_analysis'] = 1 # Init to ok
            # Spec file
            if self.spec_fil is not None:
                newline.analy['datafile'] = self.spec_fil
            # Append
            self.abs_lines.append(newline)

    def remove_line(self, wrest):
        """ Remove a line, if it exists
        Parameters
        ----------
        wrest : Quantity
        """
        awrest = [iline.wrest for iline in self.abs_lines]
        try:
            idx = awrest.index(wrest)
        except ValueError:
            return None
        else:
            _ = self.abs_lines.pop(idx)

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
            flg = ltgu.navigate(self.psdict,event)
        if event.key == '-':
            self.idx_line = max(0, self.idx_line-self.sub_xy[0]*self.sub_xy[1]) # Min=0
            if self.idx_line == sv_idx:
                print('Edge of list')
        if event.key == '=':
            self.idx_line = min(len(self.llist['show_line'])-self.sub_xy[0]*self.sub_xy[1],
                                self.idx_line + self.sub_xy[0]*self.sub_xy[1]) 
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            if self.idx_line == sv_idx:
                print('Edge of list')

        ## Reset z
        if event.key == 'z': 
            newz = ltu.z_from_v(self.z, event.xdata)
            self.z = newz
            self.abs_sys.zabs = newz
            # Drawing
            self.psdict['xmnx'] = self.vmnx.value

        # Single line command
        if event.key in ['1','2','B','U','L','N','V','A', 'x', 'X',
                         '^', '&']:
            try:
                wrest = event.inaxes.get_gid()
            except AttributeError:
                return
            else:
                absline = self.grab_line(wrest)
                kwrest = wrest.value

        ## Velocity limits
        unit = u.km/u.s
        if event.key == '1': 
            absline.analy['vlim'][0] = event.xdata*unit
        if event.key == '2': 
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            absline.analy['vlim'][1] = event.xdata*unit
        if event.key == '!': 
            for iline in self.abs_sys.lines:
                iline.analy['vlim'][0] = event.xdata*unit
        if event.key == '@': 
            for iline in self.abs_sys.lines:
                iline.analy['vlim'][1] = event.xdata*unit
        ## Line type
        if event.key == 'A': # Add to lines
            self.generate_line((self.z,wrest)) 
        if event.key == 'x': # Remove line
            if self.remove_line(wrest):
                print('VelPlot: Removed line {:g}'.format(wrest))
        if event.key == 'X': # Remove all lines 
            # Double check
            gui = xguiu.WarningWidg('About to remove all lines. \n  Continue??')
            gui.exec_()
            if gui.ans is False:
                return
            #
            self.abs_lines = []  # Flush??
        # Kinematics
        if event.key == '^':  # Low-Ion
            try:
                fkin = absline.analy['flag_kin']
            except KeyError:
                fkin = 0
            fkin += (-1)**(fkin % 2**1 >= 2**0) * 2**0
            absline.analy['flag_kin'] = fkin
        if event.key == '&':  # High-Ion
            try:
                fkin = absline.analy['flag_kin']
            except KeyError:
                fkin = 0
            fkin += (-1)**(fkin % 2**2 >= 2**1) * 2**1
            absline.analy['flag_kin'] = fkin
        # Toggle blend
        if event.key == 'B':
            try:
                feye = absline.analy['flg_eye'] 
            except KeyError:
                feye = 0
            feye = (feye + 1) % 2
            absline.analy['flg_eye']  = feye
        # Toggle NG
        if event.key == 'N':
            try:
                fanly = absline.analy['do_analysis'] 
            except KeyError:
                fanly = 1
            if fanly == 0:
                fanly = 1
            else:
                fanly = 0
            absline.analy['do_analysis']  = fanly
        if event.key == 'V':  # Normal
            absline.analy['flg_limit'] = 1
        if event.key == 'L':  # Lower limit
            absline.analy['flg_limit'] = 2
        if event.key == 'U':  # Upper limit
            absline.analy['flg_limit'] = 3
            
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
            #print('idx_l={:d}, nplt={:d}, lall={:d}'.format(self.idx_line,nplt,len(all_idx)))
            for jj in range(min(nplt, len(all_idx))):
                try:
                    idx = all_idx[jj+self.idx_line]
                except IndexError:
                    continue # Likely too few lines
                #print('jj={:d}, idx={:d}'.format(jj,idx))
                # Grab line
                wrest = self.llist[self.llist['List']].wrest[idx] 
                kwrest = wrest.value # For the Dict
                # Single window?
                if in_wrest is not None:
                    if np.abs(wrest-in_wrest) > (1e-3*u.AA):
                        continue

                # Abs_Sys: Color the lines
                if self.abs_sys is not None:
                    absline = self.grab_line(wrest)

                # Generate plot
                self.ax = self.fig.add_subplot(self.sub_xy[0],self.sub_xy[1], subp_idx[jj])
                self.ax.clear()        
                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()

                # Zero line
                self.ax.plot( [0., 0.], [-1e9, 1e9], ':', color='gray')
                # Velocity
                wvobs = (1+self.z) * wrest
                velo = (self.spec.dispersion/wvobs - 1.)*const.c.to('km/s')
                
                # Plot
                self.ax.plot(velo, self.spec.flux, 'k-',drawstyle='steps-mid')

                # GID for referencing
                self.ax.set_gid(wrest)

                # Labels
                #if jj >= (self.sub_xy[0]-1)*(self.sub_xy[1]):
                if (((jj+1) % self.sub_xy[0]) == 0) or ((jj+1) == len(all_idx)):
                    self.ax.set_xlabel('Relative Velocity (km/s)')
                else:
                    self.ax.get_xaxis().set_ticks([])
                lbl = self.llist[self.llist['List']].name[idx]
                # Kinematics
                kinl = ''
                if absline is not None:
                    if (absline.analy['flag_kin'] % 2) >= 1:
                        kinl = kinl + 'L'
                    if (absline.analy['flag_kin'] % 4) >= 2:
                        kinl = kinl + 'H'
                self.ax.text(0.1, 0.05, lbl+kinl, color='blue', transform=self.ax.transAxes,
                             size='x-small', ha='left')

                # Reset window limits
                #QtCore.pyqtRemoveInputHook()
                #xdb.set_trace()
                #QtCore.pyqtRestoreInputHook()
                self.ax.set_xlim(self.psdict['xmnx'])

                # Rescale?
                if (rescale is True) & (self.norm is False):
                    gdp = np.where( (velo.value > self.psdict['xmnx'][0]) &
                                    (velo.value < self.psdict['xmnx'][1]))[0]
                    if len(gdp) > 5:
                        per = xstats.basic.perc(self.spec.flux[gdp])
                        self.ax.set_ylim((0., 1.1*per[1]))
                    else:
                        self.ax.set_ylim(self.psdict['ymnx'])
                else:
                    self.ax.set_ylim(self.psdict['ymnx'])

                # Fonts
                xputils.set_fontsize(self.ax,6.)


                clr='black'
                if absline is not None:
                    try:
                        vlim = absline.analy['vlim']
                    except KeyError:
                        pass
                    # Color coding
                    try:  # .clm style
                        flag = absline.analy['FLAGS'][0]
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
                        flagA = absline.analy['do_analysis']
                    except KeyError:
                        flagA = None
                    else:
                        if (flagA>0) & (clr == 'black'):
                            clr = 'green'
                    try: # Limit?
                        flagL = absline.analy['flg_limit']
                    except KeyError:
                        flagL = None
                    else:
                        if flagL == 2:
                            clr = 'blue'
                        if flagL == 3:
                            clr = 'purple'
                    try: # Blends?
                        flagE = absline.analy['flg_eye']
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
    def __init__(self, spec, z, wrest, parent=None, vmnx=[-300., 300.]*u.km/u.s,
                 norm=True, linelist=None):
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
        if linelist is None:
            self.linelist = LineList('ISM')
        for iwrest in self.wrest:
            self.lines.append(AbsLine(iwrest,linelist=self.linelist))


        self.psdict = {} # Dict for spectra plotting
        self.psdict['xmnx'] = self.vmnx.value # Too painful to use units here
        self.psdict['ymnx'] = [-0.1, 1.1]
        self.psdict['nav'] = ltgu.navigate(0,0,init=True)

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
            flg = ltgu.navigate(self.psdict,event)
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
            velo = (self.spec.dispersion/wvobs - 1.)*const.c.to('km/s')
            gdp = np.where((velo.value > self.psdict['xmnx'][0]) &
                           (velo.value < self.psdict['xmnx'][1]))[0]

            # Normalize?
            if self.norm is False:
                per = xstats.basic.perc(self.spec.flux[gdp])
                fsplice = per[1] / self.spec.flux[gdp] 
            else:
                fsplice = 1./ self.spec.flux[gdp]

            # AODM
            cst = (10.**14.5761)/(self.lines[ii].data['f']*iwrest.value)
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

        # Draw
        self.canvas.draw()
    
# ##################################
# GUI for simple text entering
class EnterTextGUI(QtGui.QDialog):
    ''' GUI to grab text from the user
        29-Julc-2014 by JXP
    '''
    def __init__(self, directions='Enter:', parent=None):
        '''
        message = str
          Message to display
        '''
        super(EnterTextGUI, self).__init__(parent)

        # Initialize
        self.text = ''

        #age textbox
        textW = QtGui.QWidget()
        textlabel = QtGui.QLabel(directions)
        self.textbox = QtGui.QLineEdit()
        self.connect(self.textbox,QtCore.SIGNAL('editingFinished ()'), 
            self.set_text)
       # self.ageerror = QtGui.QLabel('')

        textvbox = QtGui.QVBoxLayout()
        textvbox.addWidget(textlabel)
        textvbox.addWidget(self.textbox)
        textW.setLayout(textvbox)

        # Buttons
        donebtn = QtGui.QPushButton('Done', self)
        donebtn.clicked.connect(self.touch_done)
        donebtn.setAutoDefault(False)

        # Main Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(textW)
        vbox.addWidget(donebtn)
        #vbox.addWidget(self.cntdwn)
        self.setLayout(vbox)

    def set_text(self):
        self.text = str(self.textbox.text())
            
    def touch_done(self):
        self.done(0)



# ################
# TESTING
if __name__ == "__main__":
    from xastropy import spec as xspec

    if len(sys.argv) == 1: #
        flg_tst = 0
        #flg_tst += 2**0  # ExamineSpecWidget
        #flg_tst += 2**1  # PlotLinesWidget
        #flg_tst += 2**2  # SelectLineWidget
        #flg_tst += 2**3  # AbsSysWidget
        #flg_tst += 2**4  # VelPltWidget
        #flg_tst += 2**5  # SelectedLinesWidget
        #flg_tst += 2**6  # AODMWidget
        flg_tst += 2**7  # Simple Text Widget
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
        specf = 0
        if specf == 0: # PH957 DLA
            # Spectrum
            spec_fil = '/u/xavier/Keck/HIRES/RedData/PH957/PH957_f.fits'
            spec = lsi.readspec(spec_fil)
            # Abs_sys
            abs_sys = xiaa.GenericAbsSystem()
            ion_fil = '/Users/xavier/DLA/Abund/Tables/PH957.z2309.ion'
            abs_sys.zabs = 2.309
            abs_sys.read_ion_file(ion_fil)
        elif specf == 1: # UM184 LLS
            # Spectrum
            spec_fil = '/Users/xavier/PROGETTI/LLSZ3/data/normalize/UM184_nF.fits'
            spec = lsi.readspec(spec_fil)
            # Abs_sys
            abs_fil = '/Users/xavier/paper/LLS/Optical/Data/Analysis/MAGE/UM184_z2.930_id.fits'
            abs_sys = xiaa.GenericAbsSystem()
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
        llist = ltgu.set_llist('ISM')
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

    # Simple text input widget
    if (flg_tst % 2**8) >= 2**7:
        app = QtGui.QApplication(sys.argv)
        app.setApplicationName('TEXT')
        main = EnterTextGUI('Enter some text:')
        main.exec_()
        print('You entered: {:s}'.format(main.text))
        sys.exit()
