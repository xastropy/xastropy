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

from astropy.io import ascii
from astropy.units.quantity import Quantity
from astropy import units as u
from astropy.table import Table, Column

from linetools.spectralline import AbsLine
from linetools.analysis import absline as ltaa

from xastropy.atomic import ionization as xai
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
    _data -- Table containing the Ion info
    
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

    # Read a dict of ions info
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
        names=('Z', 'ion', 'logN', 'sig_logN', 'flg_clm', 'flg_inst')
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
                logN, siglogN = ltaa.sum_logN(sdict,row)
                newIC._data['logN'][idx] = logN
                # Error
                newIC._data['sig_logN'][idx] = siglogN
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
        tmp += self._data.__repr__()
        return tmp

## ###################
##
# Class generated when parsing (Mainly useful for AbsSys)
class Ionic_Clm_File(object):
    """Ionic column densities for an absorption system
    DEPRECATED

    Attributes:
        clm_fil: Systemic redshift
    """
    # Initialize with a .clm file
    def __init__(self, clm_fil, linelist):
        #
        self.clm_fil = clm_fil
        self.linelist = linelist
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
            self.clm_lines[key] = AbsLine(key*u.AA,closest=True,linelist=self.linelist)
            #self.clm_lines[key] = xa.spec.analysis.Spectral_Line(key)
            self.clm_lines[key].analy['FLAGS'] = ionflg, int(tmp[3].strip())
            # By-hand
            if ionflg >= 8:
                self.clm_lines[key].attrib['N'] = 10.**vmin
                self.clm_lines[key].attrib['SIGN'] = (10.**(vmin+vmax) - 10.**(vmin-vmax))/2	      
            else:
                self.clm_lines[key].analy['VLIM']= [vmin,vmax]

# Read a .all file
def read_all_file(all_file,components=None,verbose=False):
    """Read in JXP-style .all file in an appropriate manner

    NOTE: If program breaks in this function, check the all file
    to see if it is properly formatted.

    Fills components if inputted

    Parameters
    ----------
    all_file : str
      Full path to the .all file
    components : list, optional
      List of AbsComponent objects
    """
    # Read
    if verbose:
        print('Reading {:s}'.format(all_file))
    names=('Z', 'ion', 'logN', 'sig_logN', 'flag_N', 'flg_inst') # was using flg_clm
    table = ascii.read(all_file, format='no_header', names=names)

    # Fill components
    if components is not None:
        allZ = np.array([comp.Zion[0] for comp in components])
        allion = np.array([comp.Zion[1] for comp in components])
        # Loop
        for row in table:
            mt = np.where((allZ==row['Z'])&(allion==row['ion']))[0]
            if len(mt) == 0:
                pass
            elif len(mt) == 1:
                # Fill
                components[mt[0]].flag_N = row['flag_N']
                components[mt[0]].logN = row['logN']
                components[mt[0]].sig_logN = row['sig_logN']
            else:
                raise ValueError("Found multiple component matches in read_all_file")
    # Write
    return table

def read_clmfile(clm_file,linelist=None):
    """ Read in a .CLM file in an appropriate manner

    NOTE: If program breaks in this function, check the clm to see if it is properly formatted.


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

    Parameters
    ----------
    clm_file : str
      Full path to the .clm file
    linelist : LineList
      can speed up performance
    """
    clm_dict = {}
    # Read file
    f=open(clm_file, 'r')
    arr=f.readlines()
    f.close()
    nline = len(arr)
    # Data files
    clm_dict['flg_data'] = int(arr[1][:-1])
    clm_dict['fits_files']={}
    ii=2
    for jj in range(0,6):
        if (clm_dict['flg_data'] % (2**(jj+1))) > (2**jj - 1):
            clm_dict['fits_files'][2**jj] = arr[ii].strip()
            ii += 1

    # Redshift
    clm_dict['zsys']=float(arr[ii][:-1]) ; ii+=1
    clm_dict['ion_fil']=arr[ii].strip() ; ii+=1
    # NHI
    tmp = arr[ii].split(',') ; ii+=1
    if len(tmp) != 2:
        raise ValueError('ionic_clm: Bad formatting {:s} in {:s}'
                                       .format(arr[ii-1],clm_file))
    clm_dict['NHI']=float(tmp[0])
    clm_dict['sigNHI']=float(tmp[1])
    # Abundances by hand
    numhand=int(arr[ii][:-1]) ; ii+=1
    clm_dict['fixabund']={}
    if numhand>0:
        for jj in range(numhand):
            # Atomic number
            atom=int(arr[ii][:-1]) ; ii+=1
            # Values
            tmp = arr[ii].strip().split(',') ; ii+=1
            clm_dict['fixabund'][atom]= float(tmp[0]), float(tmp[1]), int(tmp[2])
    # Loop on lines
    clm_dict['lines'] = {}
    while ii < (nline-1):
        # No empty lines allowed
        if len(arr[ii].strip()) == 0:
           break
        # Read flag
        ionflg = int(arr[ii].strip()); ii+=1
        # Read the rest
        tmp = arr[ii].split(',') ; ii+=1
        if len(tmp) != 4: raise ValueError('ionic_clm: Bad formatting {:s} in {:s}'
                                        .format(arr[ii-1],clm_file))
        vmin = float(tmp[1].strip())
        vmax = float(tmp[2].strip())
        key = float(tmp[0].strip()) # Using a float not string!
        # Generate
        clm_dict['lines'][key] = AbsLine(key*u.AA,closest=True,linelist=linelist)
        clm_dict['lines'][key].attrib['z'] = clm_dict['zsys']
        clm_dict['lines'][key].analy['FLAGS'] = ionflg, int(tmp[3].strip())
        # By-hand
        if ionflg >= 8:
            clm_dict['lines'][key].attrib['N'] = 10.**vmin / u.cm**2
            clm_dict['lines'][key].attrib['sig_N'] = (10.**(vmin+vmax) - 10.**(vmin-vmax))/2/u.cm**2
        else:
            clm_dict['lines'][key].analy['vlim']= [vmin,vmax]*u.km/u.s
    # Return
    return clm_dict

def read_ion_file(ion_fil,lines=None,components=None,linelist=None,toler=0.05*u.AA):
    """ Read in JXP-style .ion file in an appropriate manner

    NOTE: If program breaks in this function, check the .ion file
    to see if it is properly formatted.

    If components is passed in, these are filled as applicable.

    Parameters
    ----------
    ion_fil : str
      Full path to .ion file
    lines : list, optional
      List of AbsLine objects [used for historical reasons, mainly]
    components : list, optional
      List of AbsComponent objects
    linelist : LineList
      May speed up performance
    toler : Quantity, optional
      Tolerance for matching wrest
    """
    # Read
    names=('wrest', 'logN', 'sig_logN', 'flag_N', 'flg_inst')
    table = ascii.read(ion_fil, format='no_header', names=names)

    if components is None:
        if lines is None:
            lines = []
        # Generate AbsLine's
        for row in table:
            # Generate the line
            aline = AbsLine(row['wrest']*u.AA, linelist=linelist, closest=True)
            # Set z, RA, DEC, etc.
            aline.attrib['z'] = self.zabs
            aline.attrib['RA'] = self.coord.ra
            aline.attrib['Dec'] = self.coord.dec
            aline.attrib['coord'] = self.coord
            aline.attrib['logN'] = row['logN']
            aline.attrib['sig_logN'] = row['sig_logN']
            aline.attrib['flag_N'] = row['flag_N']
            aline.analy['flg_inst'] = row['flg_inst']
            # Check against existing lines
            mt = [kk for kk,oline in enumerate(lines) if oline.ismatch(aline)]
            if len(mt) > 0:
                mt.reverse()
                for imt in mt:
                    print('read_ion_file: Removing line {:g}'.format(lines[imt].wrest))
                    lines.pop(imt)
            # Append
            lines.append(aline)
            return lines
    else: # Fill entries in components
        # Generate look-up table for quick searching
        all_wv = []
        all_idx = []
        for jj,comp in enumerate(components):
            for kk,iline in enumerate(comp._abslines):
                all_wv.append(iline.wrest)
                all_idx.append((jj,kk))
        all_wv = u.Quantity(all_wv)
        # Loop now
        for row in table:
            mt = np.where(np.abs(all_wv-row['wrest']*u.AA)<toler)[0]
            if len(mt) == 0:
                pass
            elif len(mt) == 1:
                # Fill
                jj = all_idx[mt[0]][0]
                kk = all_idx[mt[0]][1]
                components[jj]._abslines[kk].attrib['flag_N'] = row['flag_N']
                components[jj]._abslines[kk].attrib['logN'] = row['logN']
                components[jj]._abslines[kk].attrib['sig_logN'] = row['sig_logN']
                components[jj]._abslines[kk].analy['flg_inst'] = row['flg_inst']
            else:
                raise ValueError("Matched multiple lines in read_ion_file")
        # Return
        return table


# Converts Flag to Instrument
def fits_flag(idx):
        # Standard dict
        fits_list = dict(zip(list(map((lambda x: 2**x),range(6))),
                             ['HIRES','ESI','UVES','XX','MIKEb','MIKEr']))
        try:
            return fits_list[idx]
        except:
            return 'Unknown'


