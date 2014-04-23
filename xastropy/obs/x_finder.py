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
from astropy.io import fits
from astropy.io import ascii 
from astropy.table import Table
from astropy.table import Column


#### ###############################
#  Deal with the RA/DEC
def get_coord(targ_file, radec=None):

    from astropy.io import ascii 
    # Import Tables
    if radec == None:
        # Read 
        ra_tab = ascii.read(targ_file) #, names=('Name','RA','DEC','Epoch'))
        # Rename the columns
        ra_tab.rename_column('col1','Name')
        ra_tab.rename_column('col2','RA')
        ra_tab.rename_column('col3','DEC')
    else: 
        # Error check
        #if len(targ_file) != 3 then stop
        # Generate the Table
        ra_tab = Table( targ_file, names=('Name','RA','DEC') )

    # Add dummy columns for decimal degrees and EPOCH
    nrow = len(ra_tab)
    col_RAD = Column(name='RAD', data=np.zeros(nrow))
    col_DECD = Column(name='DECD', data=np.zeros(nrow))
    col_EPOCH = Column(name='EPOCH', data=np.zeros(nrow))
    ra_tab.add_columns( [col_RAD, col_DECD, col_EPOCH] )
    # Assume 2000 for now
    ra_tab['EPOCH'] = 2000.
        
    return ra_tab

#### ###############################
#  Main driver
def main(targ_file, survey='2r', radec=None, deci=None, EPOCH=0.):

    # Check for wget
    #if keyword_set(SDSS) then begin
        #  ;; Check for wget
        #  spawn, 'which wget', blah
        #  if strmid(blah, 0, 1) NE '/' then begin
        #      print, 'x_fndchrt:  No wget on your machine. SDSS will not work'
        #      print, 'x_fndchrt:  Continue only if you are sure it is in your path.'
        #  endif
    # endif

    # Read in the Target list
    import x_finder as x_f
    ra_tab = x_f.get_coord(targ_file, radec=radec)

    # Grab ra, dec in decimal degrees
    if deci != None: 
        return
    # Convert to decimal degress 
    x_radec.stod(ra_tab) #ra_tab['RA'][q], ra_tab['DEC'][q], TABL)

    # Precess (as necessary)
    if EPOCH > 1000.:
        from astropy import units as u
        from astropy.coordinates import FK5
        from astropy.time import Time
        # Precess to 2000.
        tEPOCH = Time(EPOCH, format='jyear', scale='utc')
        # Load into astropy
        fk5c = FK5(ra=ra_tab['RAD'], dec=ra_tab['DECD'], equinox=tEPOCH, unit=(u.degree,u.degree))
        # Precess
        newEPOCH = Time(2000., format='jyear', scale='utc')
        newfk5 = fk5c.precess_to(newEPOCH)
        # Save
        ra_tab['RAD'] = newfk5.ra.degree
        ra_tab['DECD'] = newfk5.dec.degree
        # Strings too?
        ra_tab['RA'] = str(newfk5.ra.to_string(unit=u.hour,sep=':'))
        ra_tab['DEC'] = str(newfk5.dec.to_string(unit=u.hour,sep=':'))
            
    ## 
    # Main Loop
    nobj = len(ra_tab) 

    for q in range(nobj):

        # Grab the Image
