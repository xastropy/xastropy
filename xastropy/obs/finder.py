"""
#;+ 
#; NAME:
#; finder   
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
"""

# Import libraries
import numpy as np

from astropy.io import fits, ascii
from astropy import units as u
from astropy.table import QTable, Column

import matplotlib

from xastropy.xutils import xdebug as xdb
from xastropy.obs import radec as x_radec

#### ###############################
#  Deal with the RA/DEC
def get_coord(targ_file, radec=None):
    '''
    radec: int (None)
      None: Read from ASCII file
      1: List of [Name, RA, DEC] with RA/DEC as : separated strings
      2: List of [Name, RA, DEC] with RA/DEC as decimal degrees
    '''

    from astropy.io import ascii 
    # Import Tables
    if radec == None:
        # Read 
        ra_tab = ascii.read(targ_file) #, names=('Name','RA','DEC','Epoch'))
        # Rename the columns
        ra_tab.rename_column('col1','Name')
        ra_tab.rename_column('col2','RA')
        ra_tab.rename_column('col3','DEC')
    elif radec == 1: 
        # Error check
        if len(targ_file) != 3:
            return -1
        # Manipulate
        arr = np.array(targ_file).reshape(1,3)
        # Generate the Table
        ra_tab = QTable( arr, names=('Name','RA','DEC') )
    elif radec == 2: 
        # Error check
        if len(targ_file) != 3:
            return -1
        # Manipulate
        ras, decs = x_radec.dtos1((targ_file[1], targ_file[2]))
        # Generate the Table
        ra_tab = QTable( [ [targ_file[0]], [ras], [decs] ], names=('Name','RA','DEC') )
    else:
        raise ValueError('get_coord: Bad flag')

    # Add dummy columns for decimal degrees and EPOCH
    nrow = len(ra_tab)
    col_RAD = Column(name='RAD', data=np.zeros(nrow), unit=u.degree)
    col_DECD = Column(name='DECD', data=np.zeros(nrow), unit=u.degree)
    col_EPOCH = Column(name='EPOCH', data=np.zeros(nrow))
    ra_tab.add_columns( [col_RAD, col_DECD, col_EPOCH] )
    # Assume 2000 for now
    ra_tab['EPOCH'] = 2000.
        
    return ra_tab

#### ###############################
#  Main driver
#  finder.main(['TST', '10:31:38.87', '+25:59:02.3'], radec=1)
#  imsize is in arcmin
def main(targ_file, survey='2r', radec=None, deci=None, fpath=None,
         EPOCH=0., DSS=None, BW=False, imsize=5., show_spec=False):
    '''
    Parameters:
    ---------
    targ_file: string or List of string
       ASCII file for targets (Name, RA, DEC)
       or a List 
       Colon separated RA, DEC expected
    radec: integer (0)
       Flag indicating type of input
       0 = ASCII file
       1 = List or ['Name', 'RA', 'DEC']  
    BW: bool (False)
       B&W image?
    show_spec: bool (False)
       Try to grab and show an SDSS spectrum 
    '''
    import radec as x_r
    reload(x_r)
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    # Init
    if fpath is None:
        fpath = './'
    cradius = imsize / 50. 

    # Read in the Target list
    ra_tab = get_coord(targ_file, radec=radec)

    # Grab ra, dec in decimal degrees
    if deci != None: 
        return
    # Convert to decimal degress 

    x_r.stod(ra_tab) #ra_tab['RA'][q], ra_tab['DEC'][q], TABL)

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
    for qq in range(nobj):

        # Outfil
        nm = "".join(ra_tab['Name'][qq].split()) 
        outfil = fpath+ nm + '.pdf'
        print(outfil)

        # Grab the Image
        from xastropy.obs import x_getsdssimg as xgs
        reload(xgs)
        img, oBW = xgs.getimg(ra_tab['RAD'][qq], ra_tab['DECD'][qq], imsize, BW=BW,DSS=DSS)

        # Generate the plot
        plt.clf()
        fig = plt.figure(dpi=1200)
        fig.set_size_inches(8.0,10.5)

        # Font
        plt.rcParams['font.family']= 'times new roman'
        ticks_font = matplotlib.font_manager.FontProperties(family='times new roman', 
           style='normal', size=16, weight='normal', stretch='normal')
        ax = plt.gca()
        for label in ax.get_yticklabels() :
            label.set_fontproperties(ticks_font)
        for label in ax.get_xticklabels() :
            label.set_fontproperties(ticks_font)

        # Image
        if oBW == 1: 
            cmm = cm.Greys_r
        else: 
            cmm = None 
        plt.imshow(img,cmap=cmm,aspect='equal',extent=(-imsize/2., imsize/2, -imsize/2.,imsize/2))

        # Axes
        plt.xlim(-imsize/2., imsize/2.)
        plt.ylim(-imsize/2., imsize/2.)

        # Label
        plt.xlabel('Relative ArcMin', fontsize=20)
        xpos = 0.12*imsize
        ypos = 0.02*imsize
        plt.text(-imsize/2.-xpos, 0., 'EAST', rotation=90.,fontsize=20)
        plt.text(0.,imsize/2.+ypos, 'NORTH', fontsize=20, horizontalalignment='center')

        # Title
        plt.text(0.5,1.24, str(nm), fontsize=32, 
        horizontalalignment='center',transform=ax.transAxes)
        plt.text(0.5,1.16, 'RA (J2000) = '+str(ra_tab['RA'][qq]), fontsize=28, 
        horizontalalignment='center',transform=ax.transAxes)
        plt.text(0.5,1.10, 'DEC (J2000) = '+str(ra_tab['DEC'][qq]), fontsize=28, 
        horizontalalignment='center',transform=ax.transAxes)
        #import pdb; pdb.set_trace()

        # Circle
        circle=plt.Circle((0,0),cradius,color='y', fill=False)
        plt.gca().add_artist(circle)

        # Spectrum??
        if show_spec:
            spec_img = xgs.get_spec_img(ra_tab['RAD'][qq], ra_tab['DECD'][qq]) 
            plt.imshow(spec_img,extent=(-imsize/2.1, imsize*(-0.1), -imsize/2.1, imsize*(-0.2)))

        # Write
        if show_spec:
            plt.savefig(outfil, dpi=300)
        else:
            plt.savefig(outfil)
        print 'finder: Wrote '+outfil
        #xdb.set_trace()

    print 'finder: All done.'
    return

# ################
if __name__ == "__main__":

    flg_fig = 0 
    flg_fig += 2**0  # Test standard chart
    flg_fig += 2**1  # Test with spectrum
    
    # 
    if (flg_fig % 2**1) >= 2**0:
        main(['TST', '10:31:38.87', '+25:59:02.3'], radec=1)

    # 
    if (flg_fig % 2**2) >= 2**1:
        main(['TST2', '16:11:51.946', '+49:45:32.0'], radec=1, show_spec=True)
