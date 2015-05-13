"""
#;+ 
#; NAME:
#; abssys_utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for Absorption Systems
#;   23-Oct-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
from abc import ABCMeta, abstractmethod
from collections import OrderedDict

from astropy.io import ascii, fits 
from astropy import units as u

from xastropy.igm.abs_sys.ionic_clm import Ions_Clm, Ionic_Clm_File
from xastropy.xutils import xdebug as xdb
from xastropy import spec as xspec 
from xastropy import kinematics as xkin
from astropy.coordinates import SkyCoord

###################### ######################
###################### ######################
###################### ######################
# Class for Absorption Line System
class Absline_System(object):
    """An absorption line system

    Attributes:
        name: Coordinates
        coord: Coordinates
        zabs : float
          Absorption redshift
        NHI:  float
          Log10 of the HI column density
        sigNHI:  np.array(2)
          Log10 error of the HI column density (-/+)
        ions:  Ions_Clm Class
    """

    __metaclass__ = ABCMeta

    # Init
    def __init__(self, abs_type, zabs=0., NHI=0., MH=0., dat_file=None, tree=None, verbose=False):
        """  Initiator

        Parameters
        ----------
        abs_type : string
          Type of Abs Line System, e.g.  MgII, DLA, LLS, CGM
        dat_file : string
          ASCII .dat file summarizing the system
        """

        self.zabs = zabs
        self.NHI = NHI
        self.MH = MH
        self.coord = None

        # Abs type
        if abs_type == None:
            self.abs_type = 'NONE'
        else:
            self.abs_type = abs_type
        # Tree
        if tree == None: tree = ''
        self.tree = tree

        # Lines
        self.lines = {}  # Dict of Spectra_Line classes
        self.absid_file = None

        # Kinematics
        self.kin = {}

        # Fill in
        if dat_file != None:
            if verbose:
                print('absys_utils: Reading {:s} file'.format(dat_file))
            self.parse_dat_file(dat_file)
            self.dat_file = dat_file

        # Initialize coord
        if self.coord is None:
            ras, decs = ('00 00 00', '+00 00 00')
            self.coord = SkyCoord(ras, decs, 'icrs', unit=(u.hour, u.deg))

    # Read a .dat file
    def parse_dat_file(self,dat_file,verbose=False,flg_out=None):
        '''
        Parameters
        flg_out: int
          1: Return the dictionary
        '''
        # Define
        datdict = OrderedDict()
        # Open
        f=open(dat_file,'r')
        for line in f:
            tmp=line.split('! ')
            #tmp=line.split(' ! ')
            tkey=tmp[1].strip()
            key=tkey
            #key=tkey.replace(' ','')
            val=tmp[0].strip()
            datdict[key]=val
        f.close()
        #pdb.set_trace()

        self.datdict = datdict

        #  #########
        # Pull attributes

        # RA/DEC
        try:
            ras,decs = (datdict['RA (2000)'], datdict['DEC (2000)'])
            #print(datdict['RA(2000)'], datdict['DEC(2000)'])
            #pdb.set_trace()
        except:
            ras, decs = ('00 00 00', '+00 00 00')
        self.coord = SkyCoord(ras, decs, 'icrs', unit=(u.hour, u.deg))

        # Name
        self.name = ('J'+
                    self.coord.ra.to_string(unit=u.hour,sep='',pad=True)+
                    self.coord.dec.to_string(sep='',pad=True,alwayssign=True))

        # zabs
        try: 
            self.zabs = float(datdict['zabs'])
        except: self.zabs=0.

        # NHI
        try: 
            self.NHI = float(datdict['NHI']) # DLA format
        except:
            try:
                self.NHI = float(datdict['NHI tot']) # LLS format
            except: self.NHI=0.

        # NHIsig
        try: 
            key_sigNHI = datdict['sig(NHI)'] # DLA format
        except:
            try:
                key_sigNHI = datdict['NHI sig'] # LLS format
            except:
                key_sigNHI='0.0 0.0'
        self.sigNHI = np.array(map(float,key_sigNHI.split()))

        # Abund file
        try: 
            key_clmfil = datdict['Abund file'] # DLA format
        except:
            key_clmfil=''
        self.clm_fil = key_clmfil.strip()
        #xdb.set_trace()

        # Finish
        if verbose: print(datdict)
        if flg_out != None:
            if (flg_out % 2) == 1: ret_val = [datdict]
            else: ret_val = [0]
            return ret_val

    # Write a .dat file
    def write_dat_file(self):
        # Assuming an OrderedDict
        f=open(self.dat_file,'w')
        for key in self.datdict:
            sv = '{:60s}! {:s}\n'.format(self.datdict[key],key)
            f.write(str(sv)) # Avoids unicode
        f.close()
        print('abssys_utils.write_dat_file: Wrote {:s}'.format(self.dat_file))

    # ##
    # Parse AbsID file
    def parse_absid_file(self, abs_fil):

        # FITS binary table
        hdu = fits.open(abs_fil)
        table = hdu[1].data
        newz = table[0]['ZABS']
        if (self.zabs > 0.) & (np.abs(self.zabs-newz) > 1e-4):
            print('WARNING: Updating zabs from {:s}'.format(abs_fil))
        self.zabs = newz
        self.absid_file = abs_fil

        # Load up lines
        for row in table:
            self.lines[row['WREST']] = xspec.analysis.Spectral_Line(row['WREST'])
            # Velocity limits and flags
            try:
                self.lines[row['WREST']].analy['VLIM'] = row['VLIM']
            except KeyError:
                self.lines[row['WREST']].analy['VLIM'] = row['DV']
            self.lines[row['WREST']].analy['FLG_ANLY'] = row['FLG_ANLY']
            self.lines[row['WREST']].analy['FLG_EYE'] = row['FLG_EYE']
            self.lines[row['WREST']].analy['FLG_LIMIT'] = row['FLG_LIMIT']
            self.lines[row['WREST']].analy['DATFIL'] = row['DATFIL']
            self.lines[row['WREST']].analy['IONNM'] = row['IONNM']

    # ##
    # Write AbsID file
    def write_absid_file(self, outfil=None):

        from astropy.table import Column
        from astropy.table.table import Table

        wrest = self.lines.keys()
        wrest.sort()

        if outfil is None:
            outfil = self.absid_file

        # Columns
        cols = [Column(np.array(wrest), name='WREST')] 
        clm_nms = self.lines[wrest[0]].analy.keys()
        for clm_nm in clm_nms:
            clist = [self.lines[iwrest].analy[clm_nm] for iwrest in wrest]
            cols.append( Column(np.array(clist), name=clm_nm) )
        cols.append( Column(np.ones(len(cols[0]))*self.zabs, name='ZABS') ) 
        
        table = Table(cols)

        prihdr = fits.Header()
        prihdr['COMMENT'] = "Above are the data sources"
        prihdu = fits.PrimaryHDU(header=prihdr)
        table_hdu = fits.BinTableHDU.from_columns(np.array(table.filled()))

        thdulist = fits.HDUList([prihdu, table_hdu])
        thdulist.writeto(outfil,clobber=True)
        print('Wrote AbsID file: {:s}'.format(outfil))

    # #################
    # Parse the ion files
    def get_ions(self, skip_ions=False, fill_lines=False):
        # Read .clm file
        clm_fil=self.tree+self.clm_fil
        self.clm_analy = Ionic_Clm_File(clm_fil)
        if fill_lines is True:
            self.lines = self.clm_analy.clm_lines
        # Read .all file
        ion_fil = self.tree+self.clm_analy.ion_fil # Should check for existence
        all_fil = ion_fil.split('.ion')[0]+'.all'
        if skip_ions is False:
            self.ions = Ions_Clm(all_fil, trans_file=ion_fil)

    # #################
    # Load low_ion kinematics
    def load_low_kin(self):
        # Grab spectrum from ions
        xdb.set_trace()
        out_kin = xkin.orig_kin(spec, vmnx)

    @abstractmethod
    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        pass

    # #############
    def __repr__(self):
        return ('[Absline_System: %s %s %s %s, %g, NHI=%g]' %
                (self.name, self.abs_type,
                 self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True),
                 self.zabs, self.NHI))


# Class for Generic Absorption Line System
class Generic_System(Absline_System):
    """A simple absorption system

    """
    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'Generic'

# Class for Generic Absorption Line System
class Abs_Sub_System(Absline_System):
    """A simple absorption system

    """
    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'SubSystem'






    
###################### ###################### ######################
###################### ###################### ######################
###################### ###################### ######################
# Testing
###################### ###################### ######################
if __name__ == '__main__':

    # Test Absorption System
    tmp1 = Absline_System('LLS')
    tmp1.parse_dat_file('/Users/xavier/LLS/Data/UM669.z2927.dat')
    print(tmp1)

    #pdb.set_trace()

    # Test the Survey
    tmp = Absline_Survey('Lists/lls_metals.lst',abs_type='LLS',
                         tree='/Users/xavier/LLS/')
    print(tmp)
    print('z  NHI')
    xdb.xpcol(tmp.zabs, tmp.NHI)
    
    #xdb.set_trace()
    
    print('abssys_utils: All done testing..')
        
