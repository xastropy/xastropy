'''
#;+ 
#; NAME:
#; radec
#;    Version 1.1
#;
#; PURPOSE:
#;   2014 Written by JXP
#;-
#;------------------------------------------------------------------------------
'''

# Import libraries
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

from xastropy.xutils import xdebug as xdb

# def stod1 :: Input one RA/DEC pair as strings and return RA/DEC in decimal degrees
# def to_coord :: Input RA/DEC in one of several formats and return SkyCoord

#### ###############################
#  Main driver (return decimal values for inputted string)
def stod1(rads):
    """
    Input RA/DEC as strings and return RA/DEC in decimal degrees
    """
    # RA
    #pdb.set_trace()
    ra = np.array(rads[0].split(':'),dtype='float')
    rad = (360./24.)*(ra[0] + ra[1]/60. + ra[2]/3600.)
    # DEC
    dec = np.array(rads[1].split(':'),dtype='float')
    decd = abs(dec[0]) + dec[1]/60. + dec[2]/3600.

    #pdb.set_trace()
    # Deal with negative sign
    flg_neg = rads[1][0].strip() == '-'
    if flg_neg:
        decd = -1. * decd

    return rad*u.degree, decd*u.degree

#### ###############################
#  Decimal degress to string
def dtos1(irad, fmt=0):
    '''
    Converts a tuple of RA/DEC into J## format

    Parameters
    ----------
    rad: tuple (RA, DEC in decimal degrees [with units!])
    fmt: int (0)
      0: colon deliminated, e.g. 'J11:23:21.23', '+23:11:45.0'
      1: J name, e.g. 'J112321.23+231145.0'
    '''
    rad = list(irad)
    for ii in range(2):
        try:
            rad[ii].to('degrees')
        except AttributeError:
            rad[ii] = rad[ii] * u.degree
        
    coord = SkyCoord(ra=rad[0], dec=rad[1])
    if fmt == 0:
        ras = coord.ra.to_string(unit=u.hour,sep=':',pad=True, precision=2)
        decs = coord.dec.to_string(sep=':',pad=True, alwayssign=True, precision=1)
        return ras, decs
    elif fmt == 1:
        ras = coord.ra.to_string(unit=u.hour,sep='',pad=True, precision=2)
        decs = coord.dec.to_string(sep='',pad=True, alwayssign=True, precision=1)
        return str('J'+ras+decs)

#### ###############################
#  Deal with arrays of RA, DEC
#    Not coded yet.  Not sure we'll need it
#def stod_array(in_ra, in_dec):

#    import x_radec as x_r
#    rad = np.zeros(nrow)
#    decd = np.zeros(nrow)
#    for k in range(nrow):
#        rad[k], decd[k] = x_r.stod1( in_arr[k] )

#    return rad, decd

#### ###############################
#  Loop on rows in the Table
def stod_table(table):

    import x_radec as x_r
    # Loop on rows
    for k in range(len(table)):
        rad, decd = x_r.stod1( (table['RA'][k], table['DEC'][k]) )
        table['RAD'][k] = rad
        table['DECD'][k] = decd

        return rad, decd

#### ###############################
#  String to decimal degress
def stod(in_rads, radec=None):

    import x_radec as x_r
    from astropy.table import Table

    x_r.stod_table(in_rads)


    
#### ###############################
#  
def to_coord(irad):
    """
    Input RA/DEC as a tuple and return a SkyCoord
    """
    if not irad.__class__ is tuple:
        raise TypeError('x_radec.to_coord: Requires tuple input!')
    if len(irad) != 2:
        raise TypeError('x_radec.to_coord: Requires length two (RA,DEC)')

    # String?
    if type(irad[0]) in [str,unicode]:
        rad = stod1(irad)
    else:
        rad = irad # Assuming decimal degrees

    # Return
    return SkyCoord(ra=rad[0], dec=rad[1])
