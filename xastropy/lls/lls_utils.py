"""Python module for LLS Analysis
General Utilities
"""
import numpy as np
import pdb
import glob
from astropy.io import fits

####################################
# CLASSES
###########
class lls_sys:
    """A simple class for Landolt data"""

    __slots__ = ['Name', 'RA', 'DEC', 'V', 'B-V', 'U-B', 'V-R', 'R-I', 'V-I', 'n', 'm']
                 

    def __init__(self, Name, RA, DEC, V, BV, UB, VR, RI, VI, n=0, m=0):
        self.Name = Name
        self.RA = RA
        self.DEC = DEC
        self.V = V
        self.BV = BV
        self.UB = UB
        self.VR = VR
        self.RI = RI
        self.VI = VI
        self.n = n
        self.m = m
