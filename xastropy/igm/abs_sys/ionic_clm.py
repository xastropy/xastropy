"""
#;+ 
#; NAME:
#; ionic_clm
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for ionic column densities in Abs Systems
#;   28-Oct-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

from astropy.io import fits, ascii

from xastropy.atomic import ionization as xai
import xastropy as xa
from xastropy.xutils import xdebug as xdb

#class Ion_Clm(object):
#class Ions_Clm(object):
#class Ionic_Clm_File(object):
#def fits_flag(idx):

# Class for Ionic columns -- one ion at at time
class Ion_Clm(object):
    """Ionic column densities for an absorption system

    Attributes:
       ion: tuple (Z,ion) 
       name: string
         e.g. Si II
       flg_clm: int
         Flag describing the measurement
       clm: float
         log10 column density
       sigclm: float
         error in log10 column density
    """

    # Initialize with wavelength
    def __init__(self, ion):
        self.iZion = ion  
        self.name = xaa.ion_name(ion)
        self.lines = [] # List of transitions contributing
        #self.analy = {} # Analysis inputs (from .clm file)
        # Data
        self.flg_clm = 0
        self.clm = 0.
        self.sigclm = 0.

# ###################
# Class for Ionic columns 
class Ions_Clm(object):
    """Set of Ionic column densities for a given system

    Attributes:
    ion_data -- Dict containing the Ion info
    
    """

    # Initialize with wavelength
    def __init__(self, all_file=None, trans_file=None):
        '''
        all_file -- .all file
           File for ionic column values 
           Generally a .all file for parsing
        trans_file -- string
           File for transition-by-transition measurements
           Usually has extension .ion
        '''
        # Generate -- Other options will appear

        # Dictionary stuff
        self.keys= ('clm', 'sig_clm', 'flg_clm', 'flg_inst') 
        self.key_dtype=('f4','f4','i4','i4')
        if all_file is not None:
            self.read_all_file(all_file)
            self.all_file = all_file

        # Transitions?
        if trans_file is not None:
            self.read_ion_file(trans_file)

    # Access the data and return a dict
    def __getitem__(self, iZion):
        try:
            return self.ion_data[iZion] 
        except KeyError:
           raise KeyError 

    # Read a .all file
    def read_all_file(self,all_fil):
        """ 
        Read in the .all file in an appropriate manner
        NOTEIf program breaks in this function, check the all file 
        to see if it is properly formatted.
        """
        # Read
        print('Reading {:s}'.format(all_fil))
        names=('Z', 'ion', 'clm', 'sig_clm', 'flg_clm', 'flg_inst') 
        table = ascii.read(all_fil, format='no_header', names=names) 
        # Convert to dict
        tmp = {}
        for row in table:
            tmp[(row['Z'],row['ion'])] = {}
            for key in self.keys:
                tmp[(row['Z'],row['ion'])][key] = row[key]
        # Write
        self.ion_data = tmp

    # Read a .ion file (transitions)
    def read_ion_file(self,ion_fil):
        """ 
        Read in the .ion file in an appropriate manner
        """
        # Read
        names=('wrest', 'clm', 'sig_clm', 'flg_clm', 'flg_inst') 
        table = ascii.read(ion_fil, format='no_header', names=names) 

        # Get ion info
        adata = xa.spec.abs_line.abs_line_data( table['wrest'], ret_flg=1)

        # Add
        from astropy.table import Column
        Z = Column(adata['Z'], name='Z') # Atomic number
        ion = Column(adata['ion'], name='ion') # Atomic number
        table.add_columns([Z,ion])

        # Save
        self.trans = table


    # Printing
    def __repr__(self):
        tmp = '[Ions_Clm]\n'
        tmp += 'Z  ion  logN  sigN flgN flgI\n'
        tmp += '----------------------------\n'
        for keys in self.ion_data:
            tmp += '{:2d} {:2d} {:.3f} {:.3f}  {:d} {:3d}'.format(keys[0], keys[1],
                                          self.ion_data[keys][self.keys[0]],
                                          self.ion_data[keys][self.keys[1]],
                                          self.ion_data[keys][self.keys[2]],
                                          self.ion_data[keys][self.keys[3]] )
            tmp += '\n'
        return tmp

## ###################
##
# Class generated when parsing (Mainly useful for AbsSys)
class Ionic_Clm_File(object):
    """Ionic column densities for an absorption system

    Attributes:
        clm_fil: Systemic redshift
    """
    # Initialize with a .clm file
    def __init__(self, clm_fil):
        #
        self.clm_fil = clm_fil
        # Parse
        self.read_clmfil()

    # Read a .CLM file and return a Class
    def read_clmfil(self,linedic=None):
        """ 
        Read in the .CLM file in an appropriate manner
        NOTEIf program breaks in this function, check the clm to see if it is properly formatted.
    
        RETURNS two dictionaries CLM and LINEDIC. CLM contains the contents of CLM
        for the given DLA. THe LINEDIC that is passed (when not None) is updated appropriately.

        Keys in the CLM dictionary are:
		  INST - Instrument used
		  FITS - a list of fits files
		  ZABS - absorption redshift
		  ION - .ION file location
		  HI - THe HI column and error; [HI, HIerr]
		  FIX - Any abundances that need fixing from the ION file
		  VELS - Dictioanry of velocity limits, which is keyed by
			FLAGS - Any measurment flags assosicated with VLIM
			VLIM - velocity limits in km/s [vmin,vmax]
			ELEM - ELement (from get_elem)

        See get_elem for properties of LINEDIC
        """

        # Read file
        f=open(self.clm_fil, 'r')
        arr=f.readlines()
        f.close()
        nline = len(arr)
        #
        source=arr[0][:-1]
        # Data files
        self.flg_data = int(arr[1][:-1])
        self.fits_files={}
        ii=2
        for jj in range(0,6):
            if (self.flg_data % (2**(jj+1))) > (2**jj - 1):
                self.fits_files[2**jj] = arr[ii].strip()
                ii += 1

        # Redshift
        self.zsys=float(arr[ii][:-1]) ; ii+=1
        self.ion_fil=arr[ii].strip() ; ii+=1
        # NHI
        tmp = arr[ii].split(',') ; ii+=1
        if len(tmp) != 2:
            raise ValueError('ionic_clm: Bad formatting {:s} in {:s}'
                                           .format(arr[ii-1],self.clm_fil))
        self.NHI=float(tmp[0])
        self.sigNHI=float(tmp[1])
        # Abundances by hand
        numhand=int(arr[ii][:-1]) ; ii+=1
        self.fixabund={}
        if numhand>0:
            for jj in range(numhand):
                # Atomic number
                atom=int(arr[ii][:-1]) ; ii+=1
                # Values
                tmp = arr[ii].strip().split(',') ; ii+=1
                self.fixabund[atom]= float(tmp[0]), float(tmp[1]), int(tmp[2])
        # Loop on lines
        self.clm_lines = {}
        while ii < (nline-1):
            # No empty lines allowed
            if len(arr[ii].strip()) == 0:
               break
            # Read flag
            ionflg = int(arr[ii].strip()); ii+=1
            # Read the rest
            tmp = arr[ii].split(',') ; ii+=1
            if len(tmp) != 4: raise ValueError('ionic_clm: Bad formatting {:s} in {:s}'
                                            .format(arr[ii-1],self.clm_fil))
            vmin = float(tmp[1].strip())
            vmax = float(tmp[2].strip())
            key = float(tmp[0].strip()) # Using a float not string!
            # Generate
            self.clm_lines[key] = xa.spec.analysis.Spectral_Line(key)
            self.clm_lines[key].analy['FLAGS'] = ionflg, int(tmp[3].strip())
            # By-hand
            if ionflg >= 8:
                self.clm_lines[key].measure['N'] = 10.**vmin
                self.clm_lines[key].measure['SIGN'] = (10.**(vmin+vmax) - 10.**(vmin-vmax))/2	      
            else:
                self.clm_lines[key].analy['VLIM']= [vmin,vmax]


# Converts Flag to Instrument
def fits_flag(idx):
        # Standard dict
        fits_list = dict(zip(list(map((lambda x: 2**x),range(6))),
                             ['HIRES','ESI','UVES','XX','MIKEb','MIKEr']))
        try:
            return fits_list[idx]
        except:
            return 'Unknown'

        
