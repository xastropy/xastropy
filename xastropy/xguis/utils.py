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
                    gdlines = u.AA*[702.332, 787.711, 832.927, 972.5367, 977.0201, 
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
def read_spec(ispec, second_file=None):
    '''Parse spectrum out of the input
    Parameters:
    -----------
    ispec: Spectrum1D, str, or tuple

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
        # Second file?
        if not second_file is None:
            spec2 = lsi.readspec(second_file)
            if spec2.sig is None:
                spec2.sig = np.zeros(spec.flux.size)
            # Scale for convenience of plotting
            xper1 = xsb.perc(spec.flux, per=0.9)
            xper2 = xsb.perc(spec2.flux, per=0.9)
            scl = xper1[1]/xper2[1]
            # Stitch together
            wave3 = np.append(spec.dispersion.value, spec2.dispersion.value)
            flux3 = np.append(spec.flux, spec2.flux*scl)
            sig3 = np.append(spec.sig, spec2.sig*scl)
            #QtCore.pyqtRemoveInputHook()
            #xdb.set_trace()
            #QtCore.pyqtRestoreInputHook()
            spec3 = lsx.XSpectrum1D.from_array(wave3*u.AA, flux3, 
                uncertainty=StdDevUncertainty(sig3))
            # Overwrite
            spec = spec3
            spec.filename = spec_fil
    elif isinstance(ispec,Spectrum1D):
        spec = ispec # Assuming Spectrum1D
        spec_fil = spec.filename # Grab from Spectrum1D 
    elif isinstance(ispec,tuple):
        spec = lsx.XSpectrum1D.from_tuple(tuple)
        spec_fil = 'none'
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
