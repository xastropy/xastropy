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
from astropy.table import QTable, Column
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.units import Quantity

from xastropy.xutils import xdebug as xdb

# def stod1 :: Input one RA/DEC pair as strings and return RA/DEC in decimal degrees
# def to_coord :: Input RA/DEC in one of several formats and return SkyCoord

#### ###############################
#  Main driver (return decimal values for inputted string)
def stod1(rads):
    """
    Input RA/DEC as strings and return RA/DEC in decimal degrees

    Parameters:
    ----------
    rads:
        tuple (RA, DEC as a string with colon format)
        string (JXXXXXX.X+XXXXXX.X format)

    Returns:
    ----------
    rad: tuple (RA, DEC in decimal degrees with units)

    """
    # JXXXXXXX+XXXXXX format
    if len(rads) > 2:
        i0 = rads.find('J')+1
        isign = np.max( [rads.find('+'), rads.find('-')] )
        # Parse
        ra = np.array( [rads[i0:i0+2], rads[i0+2:i0+4], rads[i0+4:isign]], dtype='float')
        dec = np.array( [rads[isign:isign+3], rads[isign+3:isign+5], rads[isign+5:]], dtype='float')
        flg_neg = rads[isign] == '-'
    else:
        # Look for a colon
        if rads[0].find(':') == -1:
            # No colons
            ra = np.array(rads[0].split(' '),dtype='float')
            dec = np.array(rads[1].split(' '),dtype='float')
        else:
            ra = np.array(rads[0].split(':'),dtype='float')
            dec = np.array(rads[1].split(':'),dtype='float')
        # Goofy length
        if len(ra) == 2:
            ra = np.append(ra, [0.])
        # Sign
        #flg_neg = (rads[1].strip())[0] == '-'
        flg_neg = '-' in rads[1].strip() # Handles $-$ too

    #xdb.set_trace()
    # RA
    rad = (360./24.)*(ra[0] + ra[1]/60. + ra[2]/3600.)
    # DEC
    decd = abs(dec[0]) + dec[1]/60. + dec[2]/3600.

    # Deal with negative sign
    if flg_neg:
        decd = -1. * decd

    return rad*u.degree, decd*u.degree

#### ###############################
#  Decimal degress or SkyCoord to string
def dtos1(irad, fmt=0):
    '''
    Converts a tuple of RA/DEC into J## format

    Parameters
    ----------
    rad: tuple (RA, DEC in decimal degrees [with units!]) or SkyCoord
    fmt: int (0)
      0: colon deliminated, e.g. '11:23:21.23', '+23:11:45.0'
      1: J name, e.g. 'J112321.23+231145.0'
    '''
    # Get to SkyCoord
    if type(irad) is SkyCoord:
        coord = irad
    else:
        rad = list(irad)
        for ii in range(2):
            try:
                rad[ii].to('degree')
            except AttributeError:
                rad[ii] = rad[ii] * u.degree
            
        coord = to_coord(rad) 
    # String
    if fmt == 0:
        ras = coord.ra.to_string(unit=u.hour,sep=':',pad=True, precision=2)
        decs = coord.dec.to_string(sep=':',pad=True, alwayssign=True, precision=1)
        return str(ras), str(decs)
    elif fmt == 1:
        ras = coord.ra.to_string(unit=u.hour,sep='',pad=True, precision=2)
        decs = coord.dec.to_string(sep='',pad=True, alwayssign=True, precision=1)
        return str('J'+ras+decs)


#### ###############################
#  Add RA, DEC to a Table
def stod_table(table):
    ''' Converts RAS, DECS columns in a Table to RA, DEC

    Parameters
    ----------
    table:  QTable (needs to handle units)
    '''

    if type(table) is not QTable:
        raise TypeError('radec.stod_table: table needs to be a QTable')

    # Generate Columns (as needed)
    for card in ['RA','DEC']:
        if card not in table.dtype.names:
            table.add_column( Column( np.zeros(len(table)), name=card, unit=u.degree) )

    # Loop on rows
    for k,row in enumerate(table):
        rad, decd = stod1( (row['RAS'], row['DECS']) )
        table['RA'][k] = rad
        table['DEC'][k] = decd

#### ###############################
#  String to decimal degress
def stod(in_rads, radec=None):

    import radec as x_r
    from astropy.table import Table

    x_r.stod_table(in_rads)


    
#### ###############################
#  Main conversion
def to_coord(irad):
    """
    Input RA/DEC as a tuple or SkyCoord and return a SkyCoord
    """
    # SkyCoord
    if type(irad) is SkyCoord:
        return irad
    if not type(irad) in [tuple,list]: 
        if len(irad) < 10:
            raise TypeError('radec.to_coord: Expecting JXXXXXX.X+XXXXXX.X format')
    else:
        if len(irad) != 2:
            raise TypeError('radec.to_coord: Requires length two (RA,DEC)')

    # String?
    if type(irad[0]) in [str,unicode]:
        rad = stod1(irad)
    elif type(irad[0]) is Quantity:
        rad = irad 
    else: # Assuming two floats
        rad = [iirad*u.degree for iirad in irad] 

    # Return
    return SkyCoord(ra=rad[0], dec=rad[1])


#### ###############################
#  Offsets
def offsets(irad1, irad2):
    """
    Input a pair of RA/DEC and calculate the RA/DEC offsets between them

    Parameters:
    ----------
    irad1 : RA/DEC of source 1 (origin)
    irad2 : RA/DEC of source 2 (destination)

    Returns:
    -------
    offsets, PA : Tuple of offsets (itself a Tuple in arcsec) and Position Angle (degrees)
    """

    # Convert to SkyCoord
    coord1 = to_coord(irad1)
    coord2 = to_coord(irad2)

    # Angular separation
    sep = coord1.separation(coord2).to('arcsec')

    # PA
    PA = coord1.position_angle(coord2)

    # RA/DEC
    dec_off = np.cos(PA) * sep # arcsec
    ra_off = np.sin(PA) * sep # arcsec

    # Print
    print('RA Offset from 1 to 2 is {:g}'.format(ra_off))
    print('DEC Offset from 1 to 2 is {:g}'.format(dec_off))
    print('PA = {:g}'.format(PA.degree*u.degree))

    # Return
    return (ra_off, dec_off), PA.degree * u.degree
