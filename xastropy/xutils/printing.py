"""
#;+ 
#; NAME:
#; printing
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for printing routines
#;   05-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""

from __future__ import print_function

import os, pdb, imp
import numpy as np

def printcol(*args):
    """ Print column data.  Akin to XIDL printcol
    JXP on 05 Nov 2014

    Parameters
    ----------
    *args : List of vectors to print
    """

    pad = ' '

    if len(args) < 2:
        print('xutils.printing: No arguments!')
        return

    # Column widths
    widths = [max(map(len, list(str(item) for item in col))) for col in args]

    for ii in range(len(args[0])):
        string = ''
        for jj in range(len(args)):
            if jj > 0: string += pad
            #print(val)
            string += str(args[jj][ii]).ljust(widths[jj])
        print(string)
            #pdb.set_trace()

# Testing
if __name__ == '__main__':

    # printcol
    a = np.arange(1,10)
    b = a**2
    printcol(a,b)
