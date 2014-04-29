#;+ 
#; NAME:
#; x_getsdssimg
#;    Version 1.1
#;
#; PURPOSE:
#;    Returns an Image by querying the SDSS website
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

# Grab the Image (default is 202" on a side)
def getimg(ra, dec, scale=0.3961, xs=512, ys=512, grid=None,label=None,invert=None, BW=None):

    import requests
    import PIL
    from PIL import Image
    from cStringIO import StringIO

    # Generate the http call
    name1='http://skyservice.pha.jhu.edu/DR10/ImgCutout/'
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
    print 'xs = ', int(xs)
    
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

    # Request
    rtv = requests.get(name1+name)
    img = Image.open(StringIO(rtv.content))

	
    # B&W ?
    if BW != None:
        import PIL.ImageOps
        img2 = img.convert("L")
        img2 = PIL.ImageOps.invert(img2)
        img = img2

    return img

# #############
# Call with RA/DEC (decimal degrees)
def radecd(ra, dec): 
    import x_getsdssimg as x_gsdss
    img = x_gsdss.getimg(ra,dec)
    return img
