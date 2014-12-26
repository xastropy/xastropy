"""
#;+ 
#; NAME:
#; files
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for file utilities
#;   10-Dec-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os

from astropy import units as u

from xastropy.xutils import xdebug as xdb

#
def ensure_dir(fil):  # Stolen from the Web
    ''' Make sure a directory exists
    -- includes quantities
    '''
    #
    d = os.path.dirname(fil)
    if not os.path.exists(d):
        os.mkdir(d)
