#;+ 
#; NAME:
#; x_getsdssimg
#;    Version 1.1
#;
#; PURPOSE:
#;    Returns an Image by querying the SDSS website
#;      Will use DSS2-red as a backup 
#;
#; CALLING SEQUENCE:
#;
#; INPUTS:
#;
#; RETURNS:
#;
#; OUTPUTS:
#;
#; OPTIONAL KEYWORDS:
#;
#; OPTIONAL OUTPUTS:
#;
#; COMMENTS:
#;
#; EXAMPLES:
#;
#; PROCEDURES/FUNCTIONS CALLED:
#;
#; REVISION HISTORY:
#;   23-Apr-2014 Written by JXP
#;-
#;------------------------------------------------------------------------------

# Import libraries
from __future__ import print_function, absolute_import, division#, unicode_literals

import requests
from cStringIO import StringIO

from astroquery.sdss import SDSS

from astropy.coordinates import SkyCoord
from astropy import units as u

from xastropy.xutils import xdebug as xdb


# Generate the SDSS URL (default is 202" on a side)
def sdsshttp(ra, dec, imsize, scale=0.39612, grid=None, label=None, invert=None):#, xs, ys):
    '''
    Parameters:
    ----------
    imsize: float
      Image size (rectangular) in arcmin and without units
    '''

    # Pixels
    npix = round(imsize*60./scale)
    xs = npix
    ys = npix
    #from StringIO import StringIO

    # Generate the http call
    name1='http://skyservice.pha.jhu.edu/DR12/ImgCutout/'
    name='getjpeg.aspx?ra='
    
    name+=str(ra) 	#setting the ra
    name+='&dec='
    name+=str(dec)	#setting the declination
    name+='&scale='
    name+=str(scale) #setting the scale
    name+='&width='
    name+=str(int(xs))	#setting the width
    name+='&height='
    name+=str(int(ys)) 	#setting the height
    
    #------ Options
    options = ''
    if grid != None:
        options+='G'
    if label != None: 
        options+='L'
    if invert != None: 
        options+='I'
    if len(options) > 0: 
        name+='&opt='+options
        
    name+='&query='

    url = name1+name
    return url

# Generate the SDSS URL (default is 202" on a side)
def dsshttp(ra, dec, imsize):

    #https://archive.stsci.edu/cgi-bin/dss_search?v=poss2ukstu_red&r=00:42:44.35&d=+41:16:08.6&e=J2000&h=15.0&w=15.0&f=gif&c=none&fov=NONE&v3=

    Equinox = 'J2000'
    dss = 'poss2ukstu_red'
    url = "http://archive.stsci.edu/cgi-bin/dss_search?"
    url += "v="+dss+'&r='+str(ra)+'&d='+str(dec)
    url += "&e="+Equinox
    url += '&h='+str(imsize)+"&w="+str(imsize)
    url += "&f=gif"
    url += "&c=none"
    url += "&fov=NONE"
    url += "&v3="

    return url
    

# ##########################################
def getimg(ira, idec, imsize, BW=False, DSS=None):
    ''' Grab an SDSS image from the given URL, if possible

    Parameters:
    ----------
    ira: (float or Quantity) RA in decimal degrees
    idec: (float or Quantity) DEC in decimal degrees
    imsize: Image size in arcmin (without units)
    '''
    import PIL
    from PIL import Image
    # Strip units as need be
    try:
        ra = ira.value
    except AttributeError:
        ra = ira
        dec = idec
    else:
        dec = idec.value
    
    # Get URL
    if DSS == None:  # Default
        url = sdsshttp(ra,dec,imsize)
    else:
        url = dsshttp(ra,dec,imsize) # DSS

    # Request
    rtv = requests.get(url)

    # Query for photometry :: Requires updated SDSS core
    coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    phot = SDSS.query_region(coord, radius=0.02*u.deg)
    if phot is None:
        print('getimg: Pulling from DSS instead of SDSS')
        BW = 1
        url = dsshttp(ra,dec,imsize) # DSS
        rtv = requests.get(url)

    img = Image.open(StringIO(rtv.content))

    # B&W ?
    if BW:
        import PIL.ImageOps
        img2 = img.convert("L")
        img2 = PIL.ImageOps.invert(img2)
        img = img2

    return img, BW

# ##########################################
def get_spec_img(ra, dec):

    from PIL import Image
    from cStringIO import StringIO

    # Coord
    if hasattr(ra,'unit'):
        coord = SkyCoord(ra=ra, dec=dec)
    else:
        coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)

    # Query database
    radius = 1*u.arcsec
    spec_catalog = SDSS.query_region(coord,spectro=True, radius=radius.to('degree'))

    # Request
    url = 'http://skyserver.sdss.org/dr12/en/get/SpecById.ashx?id='+str(int(spec_catalog['specobjid']))
    rtv = requests.get(url) 
    img = Image.open(StringIO(rtv.content))

    return img


# #############
# Call with RA/DEC (decimal degrees)
def radecd(ra, dec): 
    import x_getsdssimg as x_gsdss
    img = x_gsdss.getimg(ra,dec)
    return img
