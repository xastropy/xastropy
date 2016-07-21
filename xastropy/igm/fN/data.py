"""
#;+ 
#; NAME:
#; fN.data
#;    Version 2.0
#;
#; PURPOSE:
#;    Module for fN data constraints
#;   12-Mar-2015 by JXP edited by Alix Feinsod
#;-
#;------------------------------------------------------------------------------
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp
from xastropy.xutils import xdebug as xdb
from xastropy.igm import tau_eff

from pyigm.fN.fnmodel import FNModel
from pyigm.fN import tau_eff as pyiteff

from astropy.io import fits

# Path for xastropy
xa_path = imp.find_module('xastropy')[1]

#class fN_Constraint(object):

class fN_Constraint(object):
    """A Class for fN constraints

    Attributes:
       fN_dtype: string
          Constraint type for the fN
           'fN' -- Standard f(N) evaluation
           'MFP' -- MFP
           'LLS' -- LLS incidence 
           'teff' -- tau effective
           'beta' -- slope constraint
       flavor: string
          Specific type of constraint
       comment: string
       ref: string
          Reference
       cosm: string
          Cosmology used (e.g. WMAP05)
       zeval: float
          Redshift where the constraint is evaluated
       data: dict
          Dictionary containing the constraints
    """

    # Initialize with type
    def __init__(self, fN_dtype, zeval=0., ref='', flavor=''):
        self.fN_dtype = fN_dtype  # Should probably check the choice
        self.zeval = zeval
        self.ref = ref
        self.flavor = flavor

    # Read from binary FITS table
    def from_fits_table(self, row):
        # Common items
        common = ['REF','COSM','TYPE','COMMENT']
        self.ref = row['REF']
        self.cosm = row['COSM']
        self.flavor = row['TYPE']
        self.comment = row['COMMENT']

        # zeval
        if 'ZEVAL' in row.array.names: self.zeval = row['ZEVAL']
        elif 'Z_LLS' in row.array.names: self.zeval = row['Z_LLS']
        elif 'Z_MFP' in row.array.names: self.zeval = row['Z_MFP']
        elif 'Z_TEFF' in row.array.names: self.zeval = row['Z_TEFF']
        else:
            raise ValueError('fN.data: No redshift info!')

        # zip the rest
        self.data = dict(zip(row.array.names,row))
        for item in common: self.data.pop(item) # No need to duplicate

    # Output
    def __repr__(self):
        return ('[%s: %s_%s z=%g, ref=%s]' %
                (self.__class__.__name__,
                 self.fN_dtype, self.flavor,
                 self.zeval, self.ref) )


# ###################### ###############
# ###################### ###############
# Read from ASCII file
def fN_data_from_ascii_file(infile):

    #makes new fN constraint with data type fN
    fNc = fN_Constraint('fN')
    ftype = fNc.fN_dtype.encode('ascii')
    fNc.fN_dtype = ftype
    fNc.ref=infile.encode('ascii')
    
    # Open file
    f = open(infile, 'r')

    # Read and ignore header lines
    firstline = f.readline()
    # get rid of newline /n symbol
    firstline =firstline.strip()
    #get zeval and DX from first line
    values = firstline.split()
    fNc.zeval = float(values[0])
    ZEVAL = float(values[0])
    DX = float(values[1])

    #declaration of variables
    BINS1 =[]
    BINS2 = []
    fn = []
    SIG_FN1 = []
    SIG_FN2 = []
    count = 0
    numlines=0

    # Loop over lines and extract info
    for line in f:
    	line = line.strip()
    	columns = line.split()
    	BINS1.append(float(columns[0]))
    	BINS2.append(float(columns[1]))
	fn.append(float(columns[2]))
	SIG_FN1.append(float(columns[3]))
	SIG_FN2.append(float(columns[3]))
	numlines +=1
	if (float(columns[0])!=0) or (float(columns[1])!=0) or (float(columns[2])!=0) or (float(columns[3])!=0):
	    count +=1
    f.close()

    NPT = int(count)
    bins = []
    bins.append(BINS1)
    bins.append(BINS2)
    sig_fn = []
    sig_fn.append(SIG_FN1)
    sig_fn.append(SIG_FN2)
    
    BINS = np.ndarray(shape=(2, numlines), dtype=float, buffer=np.array(bins))
    SIG_FN = np.ndarray(shape=(2, numlines), dtype=float, buffer=np.array(sig_fn))
    FN = np.ndarray(shape=(numlines,), dtype=float, buffer=np.array(fn))
    
    #makes array with names in ASCII not unicode
    arrayofnames = ['BINS','FN','SIG_FN','DX','NPT','ZEVAL']
    names = []
    for name in arrayofnames:
    	newname = name.encode('ascii')
    	names.append(newname)
    	
    values = [BINS,FN,SIG_FN,DX,NPT,ZEVAL]
    	
    fNc.data = dict(zip(names, values))
    
    return fNc

def fn_data_from_fits(fits_file):
    """ Build up a list of fN constraints from a multi-extension FITS file

    Parameters:
       fits_file: string
          Name of FITS file

    Returns:
       fN_list: list
          List of fN_Constraint objects

    JXP 07 Nov 2014
    """

    # List of constraints
    fN_cs = []

    # Read
    if isinstance(fits_file,list):
        for ifile in fits_file:
            tmp_cs = fn_data_from_fits(ifile)
            for cs in tmp_cs: fN_cs.append(cs)
    else:

        hdus = fits.open(fits_file)
        if len(hdus) == 1:
            raise ValueError('fN.data: Expecting a multi-extension fits file -- %s' % fits_file)
    
    
        # Loop through hdu
        for hdu in hdus[1:]:
            data = hdu.data
            # Get ftype
            if 'FN' in data.dtype.names: ftype = 'fN' # Standard f(N) data
            elif 'TAU_LIM' in data.dtype.names: ftype = 'LLS' # LLS survey
            elif 'MFP' in data.dtype.names: ftype = 'MFP' # MFP measurement
            elif 'TEFF' in data.dtype.names: ftype = 'teff' # tau effective (Lya)
            else: 
                raise ValueError('fN.data: Cannot figure out ftype')
    
            # Loop on the Table
            for row in data:
                fNc = fN_Constraint(ftype)
                fNc.from_fits_table(row)
                fN_cs.append(fNc)

    # Return
    return fN_cs


## #################################    
## #################################    
## TESTING
## #################################    
if __name__ == '__main__':

    # Read a dataset
    fn_file = xa_path+'/igm/fN/fn_constraints_z2.5_vanilla.fits'
    k13r13_file = xa_path+'/igm/fN/fn_constraints_K13R13_vanilla.fits'
    n12_file = xa_path+'/igm/fN/fn_constraints_N12_vanilla.fits'
    all_fN_cs = fn_data_from_fits([fn_file, k13r13_file, n12_file])
    #ascii_file = xa_path+'/igm/fN/asciidatan12'
    #ascii_data = fN_data_from_ascii_file(ascii_file)
    #all_fN_cs.append(ascii_data)
    
    print(all_fN_cs)
    for fN_c in all_fN_cs: print(fN_c)

    # Plot with model
    fnmodel = FNModel.default_model()
    #tst_fn_data(fN_model=fnmodel)
    xdb.set_trace()
    print('fN.data: All done testing..')

