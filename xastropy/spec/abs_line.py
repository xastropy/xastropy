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
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp
from astropy import units as u
from astropy.io import fits, ascii
from astropy.utils.misc import isiterable

from xastropy.outils import roman
from xastropy.atomic.elements import ELEMENTS
from xastropy.xutils import xdebug as xdb

# Path for xastropy
xa_path = imp.find_module('xastropy')[1]

#class Abs_Line(object):
#def mk_line_list_fits_table(outfil=None,XIDL=True):

# Class for Absorption Line  --- SHOULD COMBINE WITH Spectral_Line Class
class Abs_Line(object):
    """An absorption line

    Attributes:
        name: string
          Name of transition
        wrest: float
          Rest wavelength (Ang)
        z: float
          Redshift
    """
    # Init
    def __init__(self, wrest, z=0., N=0., fill=True, spec_file=None):
        self.wrest = wrest
        self.z = z

        # Fill in atomic data from Table (default)
        # Dict : fval, gamma, A, Elow, Eup, Ex, j, name, wrest 
        if fill is True:
            self.atomic = abs_line_data(self.wrest) # Dict
            self.name = self.atomic['name']
        else: self.name = ''

        # Attribute dict
        self.attrib = {'N': 0., 'Nsig': 0., 'flgN': 0,
                       'b': 0., 'bsig': 0., 'vmin': 0., 'vmax': 0., 
                       'EW': 0., 'EWsig': 0., 'flgEW': 0}

        # Spectrum info (optional)
        if spec_file is None:
            self.spec = None
        else:
            import xastropy.spec.readwrite as xspec_rw
            self.spec = xspec_rw.readspec(spec_file)

    # Method to find the flux-weighted optical-depth velocity (requires spectrum)
    def vpeak(self, smooth=0):
        '''
        Parameters:
        smooth: int (0)
          Number of pixels to smooth over

        Returns:
        vpeak: float [km/s]
          Flux-weighted velocity of optical depth (aka, minimum flux)

        JXP 06 Dec 2014
        '''
        # Get pixels covering the line
        import xastropy.spec.analysis as xspec_anly
        pix = xspec_anly.pixminmax(self.spec, self.z, self.wrest,
                                    (self.attrib['vmin'],self.attrib['vmax']))
        # Smooth?
        if smooth > 0:
            raise ValueError('abs_line.vpeak: Not ready for smoothing')

        # Good data
        gdp = np.where(self.spec.sig[pix] > 0)[0]
        gdpix = np.array(pix)[gdp]

        # Set tau with a floor of 0.05 in the flux
        tau = -1.* np.log( np.maximum( self.spec.flux[gdpix], 0.05*np.ones(len(gdp)) ) )

        # Weight
        vpeak = np.sum(self.spec.velo[gdpix] * tau) / np.sum( tau ) 

        # Return with units
        return vpeak * u.Unit('km/s')
            
    # Printing
    def __repr__(self):
        return '[%s: %s, %.4f]' % (self.__class__.__name__,
                                   self.name, self.wrest)

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
        self.sources = [linelist]
        self.read_llist(linelist)

    # Read standard line list
    def read_llist(self,llist, fmt=0):
        '''
        fmt: int (0)
           Format of line list.  Follows XIDL formatting..
           0: Standard absorption lines
           1: Galaxy/Quasar lines
        '''

        # Path + Format
        gdfil,fmt = llist_file(llist)
        print('gdfil = {:s}, fmt={:d}'.format(gdfil,fmt))

        if fmt == 0:
            # Read Absorption Lines with Fixed Format (astropy Table)
            self.data = ascii.read(gdfil, format='fixed_width_no_header',data_start=1,
                            names=('wrest', 'name', 'fval'),
                            col_starts=(0,9,22), col_ends=(8,20,32))
        elif fmt == 1:
            self.data = ascii.read(gdfil, format='fixed_width_no_header',data_start=1,
                            names=('wrest', 'flg', 'name'),
                            col_starts=(0,10,13), col_ends=(7,11,23))
        return

    # Add additional atomic data
    def fill_atomic(self):
        raise Exception('fill_atomic: Not ready for this yet')



## ##############
# Grab atomic data
def abs_line_data(wrest, datfil=None, ret_flg=0, tol=2e-3):
    """
    wrest : float or array
      -- Input wavelength (Ang)
    tol : float (2e-3)
      Tolerance for finding a match in wrest
    ret_flg : int (0)
      0: Return a dictionary
      1: Return an astropy Table
    """
    # Data file
    if datfil == None:
        datfil = xa_path+'/data/atomic/spec_atomic_lines.fits'
    # Read
    hdu = fits.open(datfil)
    data = hdu[1].data

    if not isiterable(wrest):
        wrest = [wrest]

    # Loop
    all_row = []
    for iwrest in wrest:
        mt = np.where(np.fabs(data['wrest']-iwrest) < tol)[0]
        nm = len(mt)
        # Found?
        if nm == 0:
            raise ValueError('abs_line_data: {:g} not in our table {:s}'.format(iwrest,datfil))
        elif nm == 1:
            # Grab
            all_row.append(mt[0])
        else:
            raise ValueError('abs_line_data: {:g} appears {:d} times in our table {:s}'.format(
                iwrest,nm,datfil))

    tab = data[all_row]

    # Return
    if ret_flg == 0: # Dictionary(ies)
        adict = []
        for row in all_row:
            adict.append(dict(zip(data.dtype.names,data[row])))
        if len(wrest) == 1:
            return adict[0]
        else:
            return adict
    elif ret_flg == 1:
        return tab
    else:
        raise Exception('abs_line_data: Not ready for this..')
    


## ##############
# Get line list path
def llist_file(llist):


    fmt = 0
    # Get the right file
    if os.path.isfile(llist): fil = llist
    elif os.path.isfile(xa_path+'/data/spec_lines/'+llist):
        fil = xa_path+'/data/spec_lines/'+llist
    else:
        # XIDL?
        fil = os.getenv('XIDL_DIR')+'/Spec/Lines/Lists/'+llist
        if not os.path.isfile(fil): 
            raise Exception('llist_file: File does not exist %s' % llist)

    # Format
    tfil = llist[max(0,llist.rfind('/')):]
    if tfil in ['gal_vac.lst']:
        fmt = 1

    return fil, fmt

## ##############
# Create line list 
def mk_line_list_fits_table(outfil=None,XIDL=False):
    from barak import absorb as ba

    if XIDL is True:
        lindat =  os.getenv('XIDL_DIR')+'/Spec/Lines/Lists/grb.lst'
        finedat = os.getenv('XIDL_DIR')+'/Spec/Lines/Lists/fine_strct.lst'
    else:
        lindat = 'grb.lst'  # This pulls from xastropy/data/spec_lines first
        finedat = os.getenv('XIDL_DIR')+'/Spec/Lines/Lists/fine_strct.lst'
  
    # Read XIDL line list
    llist = Abs_Line_List(lindat)
    ndata = len(llist.data)

    # Add columns
    from astropy.table import Column
    gamma = Column(np.zeros(ndata),name='gamma')
    A = Column(np.zeros(ndata),name='A') # Einstein coefficient
    j = Column(np.zeros(ndata),name='j') # Tot ang mom (z projection)
    Ex = Column(np.zeros(ndata),name='Ex') # Excitation energy (cm^-1)
    Elow = Column(np.zeros(ndata),name='Elow') # Energy of lower level
    Eup = Column(np.zeros(ndata),name='Eup') # Energy of upper level
    Z = Column(np.zeros(ndata,dtype='int'),name='Z') # Atomic number
    ion = Column(np.zeros(ndata,dtype='int'),name='ion') # Ionic state

    llist.data.add_columns([gamma,A,j,Ex,Elow,Eup,Z,ion])

    # Z,ion
    for ii in range(ndata):
        nm = llist.data['name'][ii]
        # Z
        if nm[1] == 'I' or nm[1] == 'V': 
            ielm = 1
        else:
            ielm = 2
        elm = nm[:ielm]
        try:
            Zv = ELEMENTS[elm].number
        except KeyError:
            if elm in ['CO','CC','HH']: # Molecules
                Zv = 999
            elif elm in ['D']: # Deuterium
                Zv = 1
            else:
                xdb.set_trace()
        llist.data['Z'][ii] = Zv
        # ion
        ispc = nm.find(' ')
        cion = nm[ielm:ispc].strip('*')
        if len(cion) == 0:
            ionv =0
        else:
            ionv = roman.fromRoman(cion)
        llist.data['ion'][ii] = ionv

    # #######
    # Fill in matches
    
    # Read atom.dat using barak
    atom, atomflat = ba.readatom(flat=True)
    #pdb.set_trace()
    llist.sources.append('atom.dat') # I wish I could pull this from ba.readatom

    # Fine structure
    fdata = ascii.read(finedat)
    llist.sources.append(finedat)

    # Loop
    for ii in range(ndata):
        # Atom.dat
        mt = np.where(np.fabs(llist.data['wrest'][ii]-atomflat['wa']) < 1e-3)[0]
        if len(mt) > 0: llist.data['gamma'][ii] = atomflat['gam'][mt[0]] # Takes the first match

        # Fine structure
        mt = np.where(np.fabs(llist.data['wrest'][ii]-fdata['wrest']) < 1e-3)[0]
        if len(mt) > 0:
            llist.data['A'][ii] = fdata['A'][mt[0]] # Takes the first match
    
    # Output file
    if outfil is None:
        outfil = xa_path+'/data/atomic/spec_atomic_lines.fits'

    # Header
    prihdr = fits.Header()
    prihdr['COMMENT'] = "Above are the data sources"
    for ii in range(len(llist.sources)):
        card = 'SOURCE'+str(ii+1)
        prihdr[card] = llist.sources[ii]
    prihdu = fits.PrimaryHDU(header=prihdr)

    # Table
    table_hdu = fits.BinTableHDU.from_columns(np.array(llist.data.filled()))

    # Write
    thdulist = fits.HDUList([prihdu, table_hdu])
    thdulist.writeto(outfil,clobber=True)
    print('mk_line_list: Wrote {:s}'.format(outfil))


## ##############
# Test.
#  Also generates spec_lines.fits
if __name__ == '__main__':

    flg_test = 0
    flg_test = 1  # Generate Line list
    #flg_test += 2 # Abs_Line Class
    #flg_test += 2**2 # Abs_Line with spectra
    #flg_test += 2**3 # abs_line_data
    #
    #flg_test += 2**9 # LLS Survey NHI
    #flg_test += 2**10 # LLS Survey ions

    
    # Generate spec_lines.fits
    if (flg_test % 2**1) >= 2**0:
        mk_line_list_fits_table()

    # Line
    if (flg_test % 2**2) >= 2**1:
        line = Abs_Line(1215.6701)
        print(line)

    # Spectrum tests for Abs_Line
    if (flg_test % 2**3) >= 2**2:
        spec_fil = '/Users/xavier/Keck/HIRES/RedData/PH957/PH957_f.fits'
        line = Abs_Line(1608.4511, z=2.309, spec_file=spec_fil)
        line.attrib['vmin'] = -20.
        line.attrib['vmax'] = 90.
        print(line.vpeak())

    # abs_line_data
    if (flg_test % 2**4) >= 2**3:
        lines = abs_line_data([1215.6701,1206.500])
        print(lines)
        lines = abs_line_data([1215.6701,1206.500], ret_flg=1)
        print(lines)

