#;+ 
#; NAME:
#; x_fndchrt   
#;    Version 1.1
#;
#; PURPOSE:
#;    Given an array of Obj name, RA, and DEC create a set of postscript
#;  finding charts. 
#;
#; CALLING SEQUENCE:
#;  x_fndchrt, targlist, OUTDIR=, IMSIZE=, SURVEY=, /ESO
#;
#; INPUTS:
#;  targlist  -- ASCII file containing  (QSO,  RA,  DEC)
#;
#; RETURNS:
#;
#; OUTPUTS:
#;
#; OPTIONAL KEYWORDS:
#;  imsize - Arcmin of image [default is 5']
#;  circ  -- Radius of green to draw about the target [default: 5"]
#;  /ESO -- Use the DSS at ESO
#;  /radec -- Input is ['Name', 'RA:RA', 'DEC:DEC']
#;  EPOCH=  -- Epoch of RA/DEC [default: 2000]
#;  SKIP=   -- Skip lines at the start of the file
#;  TWOCIRC = (x,y) offset in arcmin from field center
#;  ADDCIRC = [x,y] offset in arcmin from field center for multiple circles
#;
#; OPTIONAL OUTPUTS:
#;  OUTDIR=  -- Name of output directory
#;
#; COMMENTS:
#;
#; EXAMPLES:
#;   x_fndchrt, 'targets.list'
#;
#; PROCEDURES/FUNCTIONS CALLED:
#;  showfits
#;  querydss
#;  sdss_queryimage
#;
#; REVISION HISTORY:
#;   21-Nov-2003 Written by JXP
#;-
#;------------------------------------------------------------------------------

# Import libraries
import numpy as np
from astropy.table import Table
import pdb

#### ###############################
#  Main driver
def stod1(rads):
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

    return rad, decd

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
    #    options = {'astropy.table.table.Table': x_r.stod_table(in_rads),
        #           'list': x_r.stod1(in_rads),    # 2 elements only
        #}
    #               'numpy.ndarray': x_r.stod_array(in_rads),

    # Do the right operation
    #ty = type(in_rads)
    #print 'x_radec: ', ty, type(ty) is Table
    #if type(ty) is Table:   # isintance!!
    x_r.stod_table(in_rads)
    #elif type(ty) is 'numpy.ndarray':   # this won't work
    #    return -1
    #    # Not ready for this
    #    x_r.stod_array(in_rads)
    #else:
    #    return -1

    #print 'Success!!'


    
