"""
#;+ 
#; NAME:
#; general
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for monkeying with files and filenames 
#;   172Sep-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""

# Import libraries
import numpy as np
from astropy.io import fits
from astropy.io import ascii 
import os, pdb

#### ###############################
# Deal with .gz extensions, usually on FITS files
#   See if filenm exists, if so pass it back
#
def chk_for_gz(filenm,chk=None):

    import os, pdb

    # File exist?
    if os.path.lexists(filenm): 
        chk=1
        return filenm, chk

    # .gz already
    if filenm.find('.gz') > 0:
        chk=0
        return filenm, chk

    # Add .gz
    if os.path.lexists(filenm+'.gz'): 
        chk=1
        return filenm+'.gz', chk
    else:
        chk=0
        return filenm, chk
