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
from astropy.units import Quantity
from astropy.coordinates import SkyCoord

from linetools.spectralline import AbsLine

from xastropy.igm.abs_sys.ionclms import IonClms, Ionic_Clm_File
from xastropy.xutils import xdebug as xdb
from xastropy.atomic import ionization as xai

###################### ######################
###################### ######################
###################### ######################
# Class for Absorption Line System
class AbslineSystem(object):
    """An absorption line system

    Attributes:
        abs_type: str
        name: str
            Name of the system
        coord: SkyCoord
            RA/Dec of the sightline
        zabs : float
          Absorption redshift
        zem : float
          Emission redshift of background source
        vlim : Quantity array (2) 
          Velocity limits of the system
        NHI:  float
          Log10 of the HI column density
        sigNHI:  np.array(2)
          Log10 error of the HI column density (-/+)
        MH:  float
          Metallicity (log10)
        RA:  Quantity
          RA of the System (deg)
        Dec:  Quantity
          Dec of the System (deg)
    """

    __metaclass__ = ABCMeta

    # Init
    def __init__(self, abs_type, zabs=0., vlim=np.zeros(2)*u.km/u.s, NHI=0., MH=0., 
        name='No_Name',
        dat_file=None, tree=None, verbose=False, linelist=None, zem=0., 
        RA=None, Dec=None, sigNHI=np.zeros(2)):
        """  Initiator

        Parameters
        ----------
        abs_type : str
          Type of Abs Line System, e.g.  MgII, DLA, LLS, CGM
        dat_file : str, optional
          ASCII .dat file summarizing the system
        """

        self.name = name
        self.zabs = zabs
        self.zem = zem
        self.vlim = vlim
        self.NHI = NHI
        self.sigNHI = sigNHI
        self.MH = MH
        if (RA is not None) & (Dec is not None):
            self.coord = SkyCoord(ra=RA, dec=Dec)
        else:
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
        self.linelist = linelist
        self.lines = []  # List of SpectraLine classes
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
            radec = (0.*u.deg, 0.*u.deg)
            self.coord = SkyCoord(ra=radec[0], dec=radec[0])

        # Refs (list of references)
        self.Refs = []

    def grab_line(self,inp):
        '''Search for line in the AbslineSystem
        Parameters:
        -----------
        inp: tuple or AbsLine
          tuple -- (z,wrest) or (z,wrest,RA,Dec)

        Returns:
        -----------
        First matching AbsLine or None
        '''
        # Loop on lines
        for iline in self.lines:
            # Match?
            if iline.ismatch(inp):
                return iline 

    def remove_line(self,inp):
        '''Search for line in the AbslineSystem and 'pop' it
        Parameters:
        -----------
        inp: tuple or AbsLine
          tuple -- (z,wrest) or (z,wrest,RA,Dec)

        Returns:
        -----------
        bool -- True if popped 
        '''
        iline = self.grab_line(inp) 
        if iline is not None:
            self.lines.remove(iline)
            return True
        else:
            return False
             
    # Read a .dat file
    def parse_dat_file(self,dat_file,verbose=False,flg_out=None):
        ''' Parse an ASCII ".dat" file from JXP format 'database'
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

        # zabs
        try: 
            self.zabs = float(datdict['zabs'])
        except: self.zabs=0.

        # Name
        self.name = ('J'+
                    self.coord.ra.to_string(unit=u.hour,sep='',pad=True)+
                    self.coord.dec.to_string(sep='',pad=True,alwayssign=True)+
                    '_z{:0.3f}'.format(self.zabs))

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

        from xastropy import spec as xxspec 
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
            self.lines[row['WREST']] = xxspec.analysis.Spectral_Line(row['WREST'])
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

    # Read a .ion file (transitions)
    def read_ion_file(self,ion_fil,zabs=0.,RA=0.*u.deg, Dec=0.*u.deg):
        """Read in JXP-style .ion file in an appropriate manner

        NOTE: If program breaks in this function, check the all file 
        to see if it is properly formatted.
        """
        # Read
        names=('wrest', 'clm', 'sig_clm', 'flg_clm', 'flg_inst') 
        table = ascii.read(ion_fil, format='no_header', names=names) 

        if self.linelist is None:
            self.linelist = LineList('ISM')

        # Generate AbsLine's
        for row in table:
            # Generate the line
            aline = AbsLine(row['wrest']*u.AA, linelist=self.linelist, closest=True)
            # Set z, RA, DEC, etc.
            aline.attrib['z'] = self.zabs
            aline.attrib['RA'] = self.coord.ra
            aline.attrib['Dec'] = self.coord.dec
            # Check against existing lines
            mt = [kk for kk,oline in enumerate(self.lines) if oline.ismatch(aline)]
            if len(mt) > 0:
                mt.reverse()
                for imt in mt:
                    print('read_ion_file: Removing line {:g}'.format(self.lines[imt].wrest))
                    self.lines.pop(imt)
            # Append
            self.lines.append(aline)

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
            self._ionclms = Ions_Clm(all_fil, trans_file=ion_fil)

    # #################
    # Load low_ion kinematics
    def load_low_kin(self):
        from xastropy import kinematics as xkin
        # Grab spectrum from ions
        xdb.set_trace()
        out_kin = xkin.orig_kin(spec, vmnx)

    @abstractmethod
    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        pass

    #
    def __getitem__(self, k):
        '''Passback list of lines with this input
        See also grab_line method for quries on individual lines

        Parameters:
        -----------
        k: Quantity, tuple, or str
          float -- Rest wavelength, e.g. 1215.6701*u.AA
          tuple   -- Zion, e.g. (14,2)
          str     -- Name of the ion, e.g. 'SiII'

        Returns:
        ----------
        float: List of all matching absorption lines 
        str:   Dict of info on that ion and wavelengths of matching lines
        '''
        if isinstance(k,Quantity):  # List of AbsLines
            return [ilin for ilin in self.lines if np.abs(ilin.wrest-k)<1e-4*u.AA]
        elif isinstance(k,(basestring,tuple)):  # 
            # Column densities
            idict = self._ionclms[k]
            # Matching Absorption lines matching
            if isinstance(k,basestring):
                Zion = xai.name_ion(k)
            else:
                Zion = k
            idict[str('lines')] = [ilin.wrest for ilin in self.lines if ( 
                (ilin.data['Z']==Zion[0]) & (ilin.data['ion']==Zion[1]))]
            return idict
            #lines = [ilin for ilin in self.lines if ilin.trans==k]
        else:
            raise ValueError('Not prepared for this type')

    # #############
    def __repr__(self):
        return ('[AbslineSystem: %s %s %s %s, z=%g, NHI=%g]' %
                (self.name, self.abs_type,
                 self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True),
                 self.zabs, self.NHI))


# Class for Generic Absorption Line System
class GenericAbsSystem(AbslineSystem):
    """A simple absorption system
    """
    def __init__(self, **kwargs):
        AbslineSystem.__init__(self,'Generic',**kwargs)
        self.name = 'Foo'

    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'Generic'

# Class for Generic Absorption Line System
class Abs_Sub_System(AbslineSystem):
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
        
