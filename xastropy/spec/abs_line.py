"""
#;+ 
#; NAME:
#; abs_line
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for individual absorption lines
#;   01-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function

import numpy as np
import os, pdb
from astropy.io import fits, ascii

# Class for Absorption Line 
class Abs_Line(object):
    """An absorption line

    Attributes:
        name: Name of transition
        wrest: Rest wavelength (Ang)
    """
    # Init
    def __init__(self, wrest):
        self.wrest = wrest

        # Fill in atomic data from Table (default)
        self.atomic = abs_line_data(self.wrest)

    # Printing
    def __repr__(self):
        return '[Abs_Line: %s %g]' % (self.name, self.wrest)

# Class for Absorption Line List 
class Abs_Line_List(object):
    """An absorption line list

    Attributes:
        source: Data source(s)
        lines: List of absorption lines
    """
    # Init
    def __init__(self,linelist):
        self.lines = []

        # Fill her up
        self.source = [linelist]
        self.read_llist(linelist)

    # Read standard line list
    def read_llist(self,llist):

        # Read with Fixed Format (astropy Table)
        self.data = ascii.read(llist_file(llist), format='fixed_width_no_header',data_start=1,
                        names=('wrest', 'name', 'fval'),
                        col_starts=(0,10,22), col_ends=(8,20,32))
        return

    # Add additional atomic data
    def fill_atomic(self):
        raise Exception('fill_atomic: Not ready for this yet')



## ##############
# Get line list path
def llist_file(llist):

    # Get the right file
    if os.path.isfile(llist): fil = llist
    elif os.path.isfile(xastropy.__path__[0]+'/spec/Data/'+llist):
        fil = os.path.isfile(xastropy.__path__[0]+'/spec/Data/'+llist)
    else:
        # XIDL?
        fil = os.getenv('XIDL_DIR')+'/Spec/Lines/Lists/'+llist
        if not os.path.isfile(fil): 
            raise Exception('llist_file: File does not exist %s' % llist)
    return fil

## ##############
# Create line list 
def mk_line_list_fits_table(outfil=None,XIDL=True):
    from barak import absorb as ba

    if XIDL is True:
        lindat = 'grb.lst' 
        finedat = os.getenv('XIDL_DIR')+'/Spec/Lines/Lists/fine_strct.lst'
  
    # Read XIDL line list
    data = read_llist(lindat) # name, fval, wrest
    ndata = len(data)

    # Add columns
    gamma = Column(np.zeros(ndata),name='gamma')
    A = Column(np.zeros(ndata),name='A') # Einstein coefficient
    j = Column(np.zeros(ndata),name='j') # Tot ang mom (z projection)
    Ex = Column(np.zeros(ndata),name='Ex') # Excitation energy (cm^-1)
    Elow = Column(np.zeros(ndata),name='Elow') # Energy of lower level
    Eup = Column(np.zeros(ndata),name='Eup') # Energy of upper level

    data.add_columns([gamma,A,j,Ex,Elow,Eup])

    # #######
    # Fill in matches
    
    # Read atom.dat using barak
    atom, atomflat = ba.readatom(flat=True)

    # Fine structure
    fdata = ascii.read(finedat)

    # Loop
    for ii in range(ndata):
        # Atom.dat
        mt = np.where(np.fabs(data['wrest'][ii]-atomflat['wa']) < 1e-3)[0]
        if len(mt) > 0: data['gamma'][ii] = atomflat['gam'][mt[0]] # Takes the first match

        # Fine structure
        mt = np.where(np.fabs(data['wrest'][ii]-fdata['wrest']) < 1e-3)[0]
        if len(mt) > 0:
            data['A'][ii] = fdata['A'][mt[0]] # Takes the first match
    
    # Output
    if outfil == None:
        outfil = xastropy.__path__[0]+'/spec/Data/spec_lines.fits'
        data.write(outfil, format='fits')

## ##############
# Test.
#  Also generates spec_lines.fits
if __name__ == '__main__':
    # Generate spec_lines.fits
    mk_line_list_fits_table()
