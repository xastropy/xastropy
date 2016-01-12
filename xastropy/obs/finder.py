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
from astropy import units as astrou
from astropy.table import QTable, Column
from astropy.coordinates import SkyCoord

import matplotlib

from xastropy.obs import x_getsdssimg as xgs
from xastropy.xutils import xdebug as xdb
from xastropy.obs import radec as x_r

#### ###############################
#  Deal with the RA/DEC
def get_coord(targ_file, radec=None):
    '''
    radec: int (None)
      None: Read from ASCII file
      1: List of [Name, RA, DEC] with RA/DEC as : separated strings
      2: List of [Name, RA, DEC] with RA/DEC as decimal degrees
    '''

    if not isinstance(targ_file,basestring):
        raise IOError('Bad input to finder.get_coord!')

    from astropy.io import ascii 
    # Import Tables
    if radec == None:
        # Read 
        ra_tab = ascii.read(targ_file) #, names=('Name','RA','DEC','Epoch'))
        # Rename the columns
        ra_tab.rename_column('col1','Name')
        if isinstance(ra_tab['col2'][0],basestring):
            ra_tab.rename_column('col2','RAS')
            ra_tab.rename_column('col3','DECS')
        else:
            ra_tab.rename_column('col2','RA')
            ra_tab.rename_column('col3','DEC')
    elif radec == 1: 
        # Error check
        if len(targ_file) != 3:
            return -1
        # Manipulate
        arr = np.array(targ_file).reshape(1,3)
        # Generate the Table
        ra_tab = QTable( arr, names=('Name','RAS','DECS') )
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
#  finder.main(['TST', '10:31:38.87', '+25:59:02.3'])
#  imsize is in arcmin
def main(inp, survey='2r', radec=None, deci=None, fpath=None, show_circ=True,
         EPOCH=0., DSS=None, BW=False, imsize=5.*astrou.arcmin, show_spec=False,
         OUT_TYPE='PDF', show_another=None):
    '''
    Parameters:
    ---------
    inp: Input
       string or List of strings or List of several items
       'ra_dec_list.txt' -- ASCII file with columns of Name,RA,DEC and RA,DEC are string or float (deg)
        ['NAME_OF_TARG', '10:31:38.87', '+25:59:02.3']
        ['NAME_OF_TARG', 124.24*u.deg, -23.244*u.deg]
        ['NAME_OF_TARG', SkyCoord]
    radec: integer (0) [DEPRECATED!]
       Flag indicating type of input
       0 = ASCII file with columns of Name,RA,DEC and RA,DEC are string or float (deg)
       1 = List of string ['Name', 'RA', 'DEC']  
       2 = ['Name', ra_deg, dec_deg]
    BW: bool (False)
       B&W image?
    show_circ: bool (True)
       Show a yellow circle on the target
    show_another : tuple, optional
       RA,DEC for another target to circle (e.g. offset star)
    show_spec: bool (False)
       Try to grab and show an SDSS spectrum 
    imsize: Quantity, optional
       Image size 
    OUT_TYPE: str, optional  
       File type -- 'PDF', 'PNG'
    '''
    reload(x_r)
    reload(xgs)
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    # Init
    if fpath is None:
        fpath = './'
    try:
        imsize=imsize.to('arcmin').value
    except AttributeError:
        raise AttributeError('finder: Input imsize needs to be an Angle')
    cradius = imsize / 50. 

    # Read in the Target list
    if isinstance(inp,basestring):
        ra_tab = get_coord(targ_file, radec=radec)
    else:
        ira_tab = {}
        ira_tab['Name'] = inp[0]
        if isinstance(inp[1],basestring):
            ra, dec = x_r.stod1((inp[1],inp[2]))
            ira_tab['RA'] = ra
            ira_tab['DEC'] = dec
        elif isinstance(inp[1],float):
            ira_tab['RA'] = inp[1] * astrou.deg
            ira_tab['DEC'] = inp[2]* astrou.deg
        elif isinstance(inp[1],SkyCoord):
            ira_tab['RA'] = inp[1].ra.deg
            ira_tab['DEC'] = inp[1].dec.deg
        else: # Should check it is a Quantity
            ira_tab['RA'] = inp[1]
            ira_tab['DEC'] = inp[2]
        # Strings
        ras,decs = x_r.dtos1((ira_tab['RA'],ira_tab['DEC']))
        ira_tab['RAS'] = ras
        ira_tab['DECS'] = decs
        # Make a list
        ra_tab = [ira_tab]

    # Grab ra, dec in decimal degrees
    if deci is not None: 
        return

    #xdb.set_trace()
    #x_r.stod(ra_tab) #ra_tab['RA'][q], ra_tab['DEC'][q], TABL)

    # Precess (as necessary)
    if EPOCH > 1000.:
        from astropy import units as u
        from astropy.coordinates import FK5
        from astropy.time import Time
        # Precess to 2000.
        tEPOCH = Time(EPOCH, format='jyear', scale='utc')
        # Load into astropy
        fk5c = FK5(ra=ra_tab['RA'], dec=ra_tab['DEC'], equinox=tEPOCH, unit=(u.degree,u.degree))
        # Precess
        newEPOCH = Time(2000., format='jyear', scale='utc')
        newfk5 = fk5c.precess_to(newEPOCH)
        # Save
        ra_tab['RA'] = newfk5.ra.degree
        ra_tab['DEC'] = newfk5.dec.degree
        # Strings too?
        ra_tab['RAS'] = str(newfk5.ra.to_string(unit=u.hour,sep=':'))
        ra_tab['DECS'] = str(newfk5.dec.to_string(unit=u.hour,sep=':'))
            
    
    ## 
    # Main Loop
    for obj in ra_tab: 

        # Outfil
        nm = "".join(obj['Name'].split()) 
        if OUT_TYPE=='PNG':
            outfil = fpath+ nm + '.png'
        else:
            outfil = fpath+ nm + '.pdf'
        print(outfil)

        # Grab the Image
        reload(xgs)
        img, oBW = xgs.getimg(obj['RA'], obj['DEC'], imsize, BW=BW,DSS=DSS)

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

        #import pdb; pdb.set_trace()

        # Circle
        if show_circ:
            circle=plt.Circle((0,0),cradius,color='y', fill=False)
            plt.gca().add_artist(circle)

        # Second Circle
        if show_another is not None:
            # Coordinates
            cobj = x_r.to_coord((obj['RA'],obj['DEC']))
            canother = x_r.to_coord(show_another)
            # Offsets
            off, PA = x_r.offsets(cobj, canother)
            xanother = -1*off[0].to('arcmin').value
            yanother = off[1].to('arcmin').value
            square=matplotlib.patches.Rectangle((xanother-cradius,
                                                 yanother-cradius),
                                                cradius*2,cradius*2,color='cyan', fill=False)
            plt.gca().add_artist(square)
            plt.text(0.5, 1.24, str(nm), fontsize=32,
                 horizontalalignment='center',transform=ax.transAxes)
            plt.text(0.5, 1.18, 'RA (J2000) = '+str(obj['RAS'])+
                     '  DEC (J2000) = '+str(obj['DECS']), fontsize=22,
                     horizontalalignment='center',transform=ax.transAxes)
            plt.text(0.5, 1.12, 'RA(offset) = {:s}  DEC(offset) = {:s}'.format(
                     canother.ra.to_string(unit=astrou.hour,pad=True,sep=':', precision=2),
                     canother.dec.to_string(pad=True, alwayssign=True, sep=':', precision=1)),
                     fontsize=22, horizontalalignment='center',transform=ax.transAxes,
                     color='blue')
            plt.text(0.5, 1.06, 'RA(offset to obj) = {:g}  DEC(offset to obj) = {:g}'.format(
                     -1*off[0].to('arcsec'), -1*off[1].to('arcsec')),
                      fontsize=18, horizontalalignment='center',transform=ax.transAxes)
        else:
            # Title
            plt.text(0.5,1.24, str(nm), fontsize=32,
                     horizontalalignment='center',transform=ax.transAxes)
            plt.text(0.5,1.16, 'RA (J2000) = '+str(obj['RAS']), fontsize=28,
                     horizontalalignment='center',transform=ax.transAxes)
            plt.text(0.5,1.10, 'DEC (J2000) = '+str(obj['DECS']), fontsize=28,
                     horizontalalignment='center',transform=ax.transAxes)

        # Spectrum??
        if show_spec:
            spec_img = xgs.get_spec_img(obj['RA'], obj['DEC']) 
            plt.imshow(spec_img,extent=(-imsize/2.1, imsize*(-0.1), -imsize/2.1, imsize*(-0.2)))

        # Write
        if show_spec:
            plt.savefig(outfil, dpi=300)
        else:
            plt.savefig(outfil)
        print 'finder: Wrote '+outfil
        plt.close()
        #xdb.set_trace()

    print 'finder: All done.'
    return oBW

# ################
if __name__ == "__main__":

    flg_fig = 0 
    flg_fig += 2**0  # Test standard chart
    flg_fig += 2**1  # Test with spectrum
    
    # 
    if (flg_fig % 2**1) >= 2**0:
        main(['TST', '10:31:38.87', '+25:59:02.3'])

    # 
    if (flg_fig % 2**2) >= 2**1:
        main(['TST2', '16:11:51.946', '+49:45:32.0'], show_spec=True)
