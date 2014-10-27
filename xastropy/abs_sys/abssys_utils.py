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

# Class for Absorption Line Survey
class Absline_Survey(object):
    """A survey of absorption line systems

    Attributes:
        nsys: An integer representing the number of absorption systems
        miles:Theintegralnumberofmilesdrivenonthevehicle.
    """

    # Number of systems
    nsys = 0

    def __init__(self, filename):
        
