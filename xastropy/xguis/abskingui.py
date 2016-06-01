""" AbsKinGui for setting lines for Kinematic analysis
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Import libraries
import numpy as np
import warnings
import io
import json

from PyQt4 import QtGui
from PyQt4 import QtCore

# Matplotlib Figure object

from astropy import units as u

from linetools.guis import line_widgets as ltgl
from linetools.isgm import utils as ltiu
#from linetools.guis import spec_widgets as lspw

from xastropy.xutils import xdebug as xdb
from xastropy.xguis import spec_widgets as xspw

'''
=======
Analyzing system for future kinematic analysis

Here is now my preferred approach to perform the
analysis:

1.  Inspect the velocity plots.
2.  Identify the best low-ion transition for the analysis.
  a. High S/N
  b. Strong, but not saturated (or just barely)
  c. Preferably SiII, SII, ZnII, MgII (i.e. highly not refractory)
3.  Hit "^" on the line for a low-ion kinematic tracer
  a.  Adjust velocity limits if need be (1, 2)
4.  Hit "&" on the line for a high-ion kinematic tracer
'''


class AbsKinGui(QtGui.QDialog):
    """ GUI to analyze absorption lines for future kinematic analysis
    """
    def __init__(self, ispec, z=None, parent=None, llist=None, norm=True,
                 vmnx=[-300., 300.]*u.km/u.s, abs_sys=None, outfil=None,
                 sel_wv=None, name=''):
        """
        spec : Filename or Spectrum1D
        Norm : Bool (False)
          Normalized spectrum?
        abs_sys : AbsSystem
          Absorption system class
        sel_wv : Selected wavelength.  Used to inspect a single, unknown line
        """
        super(AbsKinGui, self).__init__(parent)

        # Initialize
        self.abs_sys = abs_sys
        if self.abs_sys is not None:
            self.z = self.abs_sys.zabs
        else:
            if z is None:
                raise ValueError('AbsKin: Need to set abs_sys or z!')
            self.z = z
        self.vmnx = vmnx
        self.outfil = outfil
        self.norm = norm
        self.sel_wv = sel_wv

        # Grab the pieces and tie together
        newfont = QtGui.QFont("Times", 10, QtGui.QFont.Bold)
        sys_label = QtGui.QLabel('Name: \n {:s}'.format(name))
        sys_label.setFont(newfont)
        self.vplt_widg = xspw.VelPlotWidget(ispec, abs_sys=self.abs_sys, llist=llist,
                                            vmnx=self.vmnx, z=self.z, norm=self.norm)
        self.pltline_widg = ltgl.PlotLinesWidget(init_llist=self.vplt_widg.llist,
                                                 init_z=self.z)
        #self.pltline_widg.spec_widg = self.vplt_widg

        self.slines = ltgl.SelectedLinesWidget(self.vplt_widg.llist[self.vplt_widg.llist['List']],
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

        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

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
        vbox.addWidget(sys_label)
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

    # Overload, as needed
    def on_key(self, event):
        pass

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
        self.slines.on_list_change(llist[llist['List']])

    # Write
    def set_outfil(self):
        self.outfil = str(self.out_box.text())
        print('AbsKin: Will write to {:s}'.format(self.outfil))

    # Set z from pltline_widg
    def setz(self):
        self.vplt_widg.abs_sys.zabs = self.pltline_widg.llist['z']
        self.vplt_widg.z = self.pltline_widg.llist['z']
        self.z = self.pltline_widg.llist['z']
        self.vplt_widg.on_draw()

    # Write
    def write_out(self):
        # Add components
        #comps = ltiu.build_components_from_abslines(self.vplt_widg.abs_lines)
        #self.vplt_widg.abs_sys._components = comps
        # Dict
        adict = self.vplt_widg.abs_sys.to_dict()
        if self.outfil is None:
            outfil = 'tmp_abskin.json'
            warnings.warn("Outfil not specified.  Using {:s}".format(outfil))
        else:
            outfil = self.outfil

        QtCore.pyqtRemoveInputHook()
        xdb.set_trace()
        QtCore.pyqtRestoreInputHook()
        print("Wrote abs_sys to {:s}".format(outfil))
        with io.open(outfil, 'w', encoding='utf-8') as f:
            f.write(unicode(json.dumps(adict, sort_keys=True, indent=4,
                                       separators=(',', ': '))))

    # Write + Quit
    def write_quit(self):
        self.write_out()
        self.flg_quit = 1
        self.abs_sys = self.vplt_widg.abs_sys
        self.done(1)

    # Write + Quit
    def quit(self):
        self.abs_sys = self.vplt_widg.abs_sys # Have to write to pass back
        self.flg_quit = 0
        self.done(1)


# Script to run XVelPltGui from the command line or ipython
def main(*args, **kwargs):
    """ Runs the AbsKinGui

    Command line
    or from Python
    Examples:
      1.  python ~/xastropy/xastropy/xguis/abskingui.py
      2.  abskingui.main(filename)
      3.  abskingui.main(spec1d)
    """
    import sys
    import argparse

    parser = argparse.ArgumentParser(description='Parse for AbsKingGui')
    parser.add_argument("file", type=str, help="Spectral file")
    parser.add_argument("-sysfile", type=str, help="System JSON file")
    parser.add_argument("-zsys", type=float, help="System Redshift")
    parser.add_argument("-outfil", type=str, help="Output filename")
    parser.add_argument("--un_norm", help="Spectrum is NOT normalized",
                        action="store_true")

    if len(args) == 0:
        pargs = parser.parse_args()
    else: # better know what you are doing!
        if isinstance(args[0],(Spectrum1D, tuple)):
            if not kwargs['rerun']:
                app = QtGui.QApplication(sys.argv)
            xdb.set_trace()
            gui = AbsKinGui(args[0], **kwargs)
            gui.exec_()
            #gui.show()
            #app.exec_()
            return gui, app
        else: # String parsing
            largs = [iargs for iargs in args]
            pargs = parser.parse_args(largs)
            xdb.set_trace() # Not setup for command line yet

    # Normalized?
    norm = True
    if pargs.un_norm:
        norm = False

    # Read AbsSystem
    from linetools.isgm.abssystem import GenericAbsSystem
    if pargs.sysfile is not None:
        abs_sys = GenericAbsSystem.from_json(pargs.sysfile, chk_vel=False)
    else:
        abs_sys = None

    app = QtGui.QApplication(sys.argv)
    gui = AbsKinGui(pargs.file, z=pargs.zsys, norm=norm, abs_sys=abs_sys, outfil=pargs.outfil)
    gui.show()
    app.exec_()

    return gui, app

if __name__ == "__main__":
    main()

    # python abskingui.py /Users/xavier/Dropbox/CASBAH/jxp_analysis/FBQS0751+2919/fbqs0751_nov2014bin.fits -zsys 0. -outfil /Users/xavier/Desktop/tmp.fits -unnorm
