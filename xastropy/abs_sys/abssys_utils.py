"""
#;+ 
#; NAME:
#; abssys_utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for Absorption Systems
#;   23-Oct-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""

import numpy as np
import pdb
from astropy.io import ascii 

# Class for Absorption Line Survey
class Absline_Survey(object):
    """A survey of absorption line systems

    Attributes:
        nsys: An integer representing the number of absorption systems
        abs_type: Type of Absorption system (DLA, LLS)
        ref: Reference to the Survey
    """

    # Number of systems

    def __init__(self, flist, abs_type=None, ref=None):
        # Expecting a list of files describing the absorption systems
        data = ascii.read(flist, data_start=0, guess=False,format='no_header')
        self.files = list(data['col1'])
        self.nsys = len(self.files)
        self.abs_type = abs_type
        self.ref = ref
        #
        print('Read %s',self.nsys)
        
