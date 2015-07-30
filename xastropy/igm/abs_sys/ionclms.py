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
import copy

from astropy.io import fits, ascii
from astropy.units.quantity import Quantity
from astropy.table import QTable, Table, Column

from linetools.spectralline import AbsLine

from xastropy.atomic import ionization as xai
import xastropy as xa
from xastropy.xutils import xdebug as xdb

#class Ion_Clm(object):
#class Ions_Clm(object):
#class Ionic_Clm_File(object):
#def fits_flag(idx):


# ###################
# Class for Ionic columns 
class IonClms(object):
    """Set of Ionic column densities for an absorption system

    Attributes:
    _data -- Dict containing the Ion info
    
    """

    # Initialize with wavelength
    def __init__(self, all_file=None, trans_file=None, idict=None):
        '''
        all_file -- .all file
           File for ionic column values 
           Generally a .all file for parsing
        idict -- Input dict
           Often from a JSON file, e.g. "HD-LLS_ions.json"
        trans_file -- string  [DEPRECATED]
           File for transition-by-transition measurements
           Usually has extension .ion
        '''
        # Generate -- Other options will appear

        # Dictionary stuff
        #self.keys= ('clm', 'sig_clm', 'flg_clm', 'flg_inst') 
        #self.key_dtype=('f4','f4','i4','i4')

        # .all file?
        if all_file is not None:
            self.from_all_file(all_file)
            self.all_file = all_file

        # dict?
        if idict is not None:
            self.from_dict(idict)

        # Transitions? [DEPRECATED]
        if trans_file is not None:
            raise ValueError('Read these in as AbsLine(s)')
            #self.read_ion_file(trans_file)

    # Read a .all file
    def from_dict(self,idict,verbose=False):
        # Manipulate for astropy Table
        table = None
        for ion in idict.keys():
            Zion = xai.name_ion(ion)
            if table is None:
                tkeys = idict[ion].keys()
                lst = [[idict[ion][tkey]] for tkey in tkeys]
                table = Table(lst, names=tkeys)
                # Extra columns
                if 'Z' not in tkeys:
                    table.add_column(Column([Zion[0]],name='Z'))
                    table.add_column(Column([Zion[1]],name='ion'))
            else:
                tdict = idict[ion]
                if 'Z' not in tkeys:
                    tdict['Z'] = Zion[0]
                    tdict['ion'] = Zion[1]
                # Add
                table.add_row(tdict)
        # Finish
        self._data = table

    # Read a .all file
    def from_all_file(self,all_fil,verbose=False):
        """Read in JXP-style .all file in an appropriate manner
        
        NOTE: If program breaks in this function, check the all file 
        to see if it is properly formatted.
        """
        # Read
        if verbose:
            print('Reading {:s}'.format(all_fil))
        names=('Z', 'ion', 'clm', 'sig_clm', 'flg_clm', 'flg_inst') 
        table = ascii.read(all_fil, format='no_header', names=names) 

        # Write
        #self._data = tmp
        self._data = table


    ##
    def sum(self,other):
        '''Sum two IonClms classes
        Parameters:
        ----------
        other: IonClms Class  
          Another one

        Returns:
        --------
        A new instance of IonClms with the column densities summed
        '''
        # Instantiate and use data form original as starting point
        newIC = IonClms()
        newIC._data = copy.deepcopy(self._data)
        # Loop through other
        for row in other._data:
            # New?
            Zion = (row['Z'], row['ion'])
            try:
                sdict = newIC[Zion]
            except KeyError:
                # Add in the new row
                newIC._data.add_row(row)
            else:
                idx = np.where((newIC.Z==Zion[0]) & (newIC.ion==Zion[1]))[0][0]
                # Clm
                logN, siglogN = sum_logN(sdict,row)
                newIC._data['clm'][idx] = logN
                # Error
                newIC._data['sig_clm'][idx] = siglogN
                '''
                np.sqrt(
                    np.sum([(sdict['sig_clm']*(10.**sdict['clm']))**2,
                    (row['sig_clm']*(10.**row['clm']))**2]))/(10.**newIC._data['clm'][idx])
                '''
                # Flag
                flags = [sdict['flg_clm'], row['flg_clm']]
                if 2 in flags:   # At least one saturated
                    flag = 2
                elif 1 in flags: # None saturated; at least one detection
                    flag = 1
                else:            # Both upper limits
                    flag = 3  
                newIC._data['flg_clm'][idx] = flag
                # Instrument (assuming binary flag)
                if 'flg_inst' in sdict.keys():
                    binflg = [0]*10
                    for jj in range(10):
                        if (row['flg_inst'] % 2**(jj+1)) >= 2**jj:
                            binflg[jj] = 1
                        if (sdict['flg_inst'] % 2**(jj+1)) >= 2**jj:
                            binflg[jj] = 1
                    newIC._data['flg_inst'][idx] = int(np.sum(
                        [2**kk for kk,ibinf in enumerate(binflg) if ibinf==1]))
            # Return
        return newIC

    #####
    def __getattr__(self,k):
        ''' Passback an array or Column of the data 
        k: str
          Must be a Column name in the data Table
        '''
        if not isinstance(k, basestring): 
            raise ValueError('Entry must be a basestring')

        # Deal with QTable
        colm = self._data[k]
        if isinstance(colm[0], Quantity):
            return self._data[k]
        else:
            return np.array(self._data[k])



    def __getitem__(self, ion):
        '''Passback a dict of measured data on a given ion

        Parameters:
        -----------
        ion: tuple or str
          tuple:  (Z,ion_state) e.g. (14,2) 
          str:  Name, e.g. 'SiII'

        Returns:
        ----------
        Dict (from row in the data table)
        '''
        if isinstance(ion,tuple):
            mt = np.where((self._data['Z'] == ion[0]) & (self._data['ion'] == ion[1]))[0]
        elif isinstance(ion,basestring):
            # Convert to tuple
            return self.__getitem__(xai.name_ion(ion))
        else:
            raise ValueError('Not prepared for this type')

        if len(mt) == 0:
            raise KeyError
        else:
            return dict(zip(self._data.dtype.names,self._data[mt][0]))

    # Printing
    def __repr__(self):
        tmp = '[IonClms]\n'
        '''
        tmp += 'Z  ion  logN  sigN flgN flgI\n'
        tmp += '----------------------------\n'
        for keys in self._data:
            tmp += '{:2d} {:2d} {:.3f} {:.3f}  {:d} {:3d}'.format(keys[0], keys[1],
                                          self._data[keys][self.keys[0]],
                                          self._data[keys][self.keys[1]],
                                          self._data[keys][self.keys[2]],
                                          self._data[keys][self.keys[3]] )
            tmp += '\n'
        '''
        tmp += self._data.__repr__()
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
            self.clm_lines[key] = AbsLine(key*u.AA)
            #self.clm_lines[key] = xa.spec.analysis.Spectral_Line(key)
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

# Converts Flag to Instrument
def sum_logN(obj1,obj2):
    '''Add log columns and return value and errors
    Parameters:
    -----------
    obj1: object
      An object with tags appropriate for the analysis
      Assumes 'clm' for column and 'sig_clm' for error for now
    obj2: object
      Another object with tags appropriate for the analysis

    Returns:
    --------
    logN, siglogN
    '''
    # Calculate
    logN = np.log10(np.sum(10.**np.array([obj1['clm'],obj2['clm']])))
    siglogN = np.sqrt(
        np.sum([(obj1['sig_clm']*(10.**obj1['clm']))**2,
        (obj2['sig_clm']*(10.**obj2['clm']))**2]))/(10.**logN)
    # Return
    return logN, siglogN



# Class for Ionic columns -- one ion at at time
class Ion_Clm(object):
    """Ionic column densities for an absorption system
    DEPRECATED

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
        
