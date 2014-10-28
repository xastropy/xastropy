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
from __future__ import print_function

import numpy as np
import pdb

# Top-Level Class for Ionic columns
class Ionic_Clm(object):
    """Ionic column densities for an absorption system

    Attributes:
        zsys: Systemic redshift
    """

    # Initialize with a .clm file
    def __init__(self, wave, clm_file=None):
        self.wave = wave
        self.data = {}

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

        clm={}
        # Read file
        f=open(self.clm_fil, 'r')
        arr=f.readlines()
        f.close()
        #
        source=arr[0][:-1]
        # Data files
        self.flg_data = int(arr[1][:-1])
        self.fits_files=[]

        anotherspec=True
        ii=0
        while anotherspec:
            if len(arr[2+ii+1][:-1].split('.fit'))>1:
                ii+=1
                fits.append(arr[3+ii][:-1])
                arr.pop(3)
            else: anotherspec=False
        clm['FITS']=fits
        zabs=float(arr[3][:-1])
        clm['ZABS']=zabs
        fion=arr[4][:-1]
        clm['ION']=fion
        if len(arr[5][:-1].split(','))>1: #When there is a comma
            NHI=float(arr[5][:-1].split(',')[0])
            HIerr=float(arr[5][:-1].split(',')[1])
            HI=[NHI,HIerr]
            clm['HI']=HI
        else:#When there is no comma (space)
            NHI=float(arr[5][:-1].split()[0])
            HIerr=float(arr[5][:-1].split()[1])
            HI=[NHI,HIerr]
            clm['HI']=HI
        numlines=int(arr[6][:-1])#Number of lines that follow with 'by hand' abund
        fixabund={}
        if numlines>0:
            for ii in range(numlines):
                atom=int(arr[6+2*ii+1][:-1])
                N=float(arr[6+2*ii+2][:-1].split(',')[0])
                Nerr=float(arr[6+2*ii+2][:-1].split(',')[1])
                inst=find_dflag(int(arr[6+2*ii+2][:-1].split(',')[2]))
                fixabund[atom]=N, Nerr, inst
        clm['FIX']=fixabund
        vels={}
        velline=6+2*numlines+1#line number where velocity limits start
        if len(arr[velline:])>1:
            for ii in range(len(arr[velline:])/2):
                tmpdic={}
                mflag=int(arr[velline+2*ii][:-1])
                wl,vmin,vmax,dflag=arr[velline+2*ii+1][:-1].split(',')

                """ 
                wl,vmin,vmax,dflag=arr[velline+2*ii+1][:-1].split()
                wl=wl.split(',')[0]
                vmin=wl.split(',')[0]
                vmax=wl.split(',')[0]
                """
                if vmin!='  -UU' and vmax!=' LL':
                    restwvl=wl.split('.')[0]+'.'+str(wl.split('.')[1][0:3])
                    elem,linedic,wl=get_elem(linelist='pybund_linelist.dat', linedic=linedic, wvl=restwvl)
                    tmpdic['FLAGS']=find_mflag(mflag), find_dflag(dflag)
                    tmpdic['VLIM']=[float(vmin),float(vmax)]
                    tmpdic['ELEM']=elem
                vels[wl]=tmpdic
        clm['VELS']=vels
        return clm, linedic

