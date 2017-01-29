"""
#;+ 
#; NAME:
#; utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for CASBAH utilities
#;   02-Jan-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp, glob, copy
from astropy.io import fits, ascii
from astropy import units as u 
from astropy.coordinates import SkyCoord

#from xastropy.xutils import xdebug as xdb


