### TESTING
import numpy as np
from astropy.io import fits
from astropy.io import ascii 
from astropy.table import Table
from astropy.table import Column
    
import xastropy
from xastropy.obs import x_radec as x_r
from xastropy.obs import x_finder as x_f

## Load the Table
#ra_tab = x_f.get_coord('Lick_2014A_z1qso.lst')

## 
#x_r.stod_table(ra_tab)
#zip(ra_tab['RAD'], ra_tab['DECD'])

