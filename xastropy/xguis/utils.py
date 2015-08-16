"""
#;+ 
#; NAME:
#; utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for Simple Guis with QT
#;   27-Dec-2014 by JXP
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

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
# Matplotlib Figure object
from matplotlib.figure import Figure

from astropy import units as u

from xastropy.xutils import xdebug as xdb

# def EditBox
# def WarningWidg

class EditBox(QtGui.QWidget):
    '''
    initv: Initial value
    lbl: str
    format: str
      Format for value
    '''
    def __init__(self, initv, lbl, format, parent=None):
        '''
        '''
        super(EditBox, self).__init__(parent)

        self.value = initv
        # 
        label = QtGui.QLabel(lbl) 
        self.box = QtGui.QLineEdit()
        # Format
        self.box.frmt = format
        self.box.setText(self.box.frmt.format(self.value))
        self.box.setMinimumWidth(90)
        # Connect
        self.connect(self.box, 
            QtCore.SIGNAL('editingFinished ()'), self.setv)
        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(label)
        vbox.addWidget(self.box)
        self.setLayout(vbox)

    def setv(self):
        self.value = unicode(self.box.text())

    def set_text(self,value):
        self.value = value
        self.box.setText(self.box.frmt.format(self.value))


# ##################################
# GUI for velocity plot
class WarningWidg(QtGui.QDialog):
    ''' GUI to warn user about coming action and solicit response
        24-Dec-2014 by JXP
    '''
    def __init__(self, message, parent=None):
        '''
        message = str
          Message to display
        '''
        super(WarningWidg, self).__init__(parent)

        # Initialize

        # Grab the pieces and tie together
        z_label = QtGui.QLabel('Warning: {:s}'.format(message))

        # Quit
        nbtn = QtGui.QPushButton('No', self)
        nbtn.clicked.connect(self.touch_no)
        ybtn = QtGui.QPushButton('Yes', self) 
        ybtn.clicked.connect(self.touch_yes)

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(z_label)
        vbox.addWidget(nbtn)
        vbox.addWidget(ybtn)
        self.setLayout(vbox)

    def touch_yes(self):
        self.ans = True
        self.done(0)

    def touch_no(self):
        self.ans = False
        self.done(0)


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
def set_llist(llist,in_dict=None,sort=True):
    ''' Method to set a line list dict for the Widgets
    sort: bool, optional
      Sort lines by rest wavelength [True]
    '''
    from linetools.lists.linelist import LineList
    from astropy.units.quantity import Quantity

    if in_dict is None:
        in_dict = {}

    if isinstance(llist,basestring): # Set line list from a file
        in_dict['List'] = llist
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
                    llist_cls = LineList('Strong', gd_lines=gdlines) 
                    in_dict[llist] = llist_cls
                else:
                    llist_cls = LineList(llist)
                    # Sort
                    llist_cls._data.sort('wrest')
                    # Load
                    in_dict[llist] = llist_cls
    elif isinstance(llist,(Quantity,list)): # Set from a list of wrest
        in_dict['List'] = 'input.lst'
        in_dict['Plot'] = True
        # Fill
        if sort:
            llist.sort()
        llist_cls = LineList('ISM', gd_lines=llist) # May need to let ISM be a choice
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
    else:
        raise IOError('Not ready for this type of input')

    # Return
    return in_dict

# Read spectrum, pass back it and spec_file name
def read_spec(ispec):
    '''Parse spectrum out of the input
    Parameters:
    -----------
    ispec: Spectrum1D, str, list of files (ordered blue to red), 
       or tuple of arrays

    Returns:
    -----------
    spec: XSpectrum1D 
    spec_file: str
    '''
    from specutils.spectrum1d import Spectrum1D
    from astropy.nddata import StdDevUncertainty
    from linetools.spectra import xspectrum1d as lsx 
    from linetools.spectra import io as lsi 
    from xastropy.stats import basic as xsb
    #
    if isinstance(ispec,basestring):
        spec_fil = ispec
        spec = lsx.XSpectrum1D.from_file(spec_fil)
    elif isinstance(ispec,Spectrum1D):
        spec = ispec 
        spec_fil = spec.filename # Grab from Spectrum1D 
    elif isinstance(ispec,tuple):
        spec = lsx.XSpectrum1D.from_tuple(ispec)
        spec_fil = 'none'
    elif isinstance(ispec,list): # Multiple file names
        # Loop on the files
        for kk,ispecf in enumerate(ispec):
            jspec = lsx.XSpectrum1D.from_file(ispecf)
            if kk == 0:
                spec = jspec
                xper1 = xsb.perc(spec.flux, per=0.9)
            else: 
                # Scale flux for convenience of plotting (sig is not scaled)
                xper2 = xsb.perc(jspec.flux, per=0.9)
                scl = xper1[1]/xper2[1]
                # Splice
                spec = spec.splice(jspec,scale=scl)
            # Filename
            spec_fil = ispec[0]
            spec.filename=spec_fil
    else:
        raise ValueError('Bad input to read_spec')

    # Return
    return spec, spec_fil


# ################
# TESTING
if __name__ == "__main__":

    flg_fig = 0 
    flg_fig += 2**0  # Warning

    if (flg_fig % 2) == 1:
        app = QtGui.QApplication(sys.argv)
        app.setApplicationName('Warning')
        main = WarningWidg('Will remove all lines. \n  Continue??')
        main.show()
        app.exec_()
        if main.ans:
            print('You answered yes!')
        else:
            print('You answered no!')
