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
import pdb

from xastropy.atomic import ionization as xai

#class Ionic_Clm(object):
#class Line_Clm(object):
#class Ionic_Clm_File(object):

# Class for Ionic columns -- ion at at time
class Ionic_Clm(object):
    """Ionic column densities for an absorption system

    Attributes:
       ion: tuple (Z,ion) 
       name: string   ,  e.g. Si II
       flg_clm: int
         Flag describing the measurement
       clm: float
         log10 column density
       sigclm: float
         error in log10 column density
    """

    # Initialize with wavelength
    def __init__(self, ion):
        self.ion = ion  
        self.name = xaa.ion_name(ion)
        self.lines = [] # List of transitions contributing
        self.analy = {} # Analysis inputs (from .clm file)
        # Data
        self.flg_clm = 0
        self.clm = 0.
        self.sigclm = 0.

# Class for Ionic columns of a given line
class Line_Clm(object):
    """Ionic column densities for an absorption system

    Attributes:
        zsys: Systemic redshift
    """

    # Initialize with wavelength
    def __init__(self, wave, clm_file=None):
        self.wave = wave
        self.atomic = {} # Atomic Data
        self.analy = {} # Analysis inputs (from .clm file)
        self.measure = {} # Measured quantities

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
        self.NHI=float(tmp[0])
        self.sigNHI=float(tmp[1])
        # Abundances by hand
        numhand=int(arr[ii][:-1]) ; ii+=1
        self.fixabund={}
        if numhand>0:
            for ii in range(ii,ii+numhand):
                # Atomic number
                atom=int(arr[ii][:-1]) ; ii+=1
                # Values
                tmp = arr[ii].split(',')
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
            vmin = float(tmp[1].strip())
            vmax = float(tmp[2].strip())
            key = float(tmp[0].strip()) # Using a float not string!
            # Generate
            self.clm_lines[key] = Line_Clm(float(tmp[0]))
            self.clm_lines[key].analy['FLAGS'] = ionflg, int(tmp[3].strip())
            # By-hand
            if ionflg >= 8:
                self.clm_lines[key].measure['N'] = 10.**vmin
                self.clm_lines[key].measure['SIGN'] = (10.**(vmin+vmax) - 10.**(vmin-vmax))/2	      
            else:
                self.clm_lines[key].analy['VLIM']= [vmin,vmax]

    # Read a .all file and return a Class
    def read_allfil(self,linedic=None):
        """ 
        Read in the .all file in an appropriate manner
        NOTEIf program breaks in this function, check the clm to see if it is properly formatted.
    
        RETURNS one dictionary for the given system. 

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

# Converts Flag to Instrument
def fits_flag(idx):
        # Standard dict
        fits_list = dict(zip(list(map((lambda x: 2**x),range(6))),
                             ['HIRES','ESI','UVES','XX','MIKEb','MIKEr']))
        try:
            return fits_list[idx]
        except:
            return 'Unknown'

        
