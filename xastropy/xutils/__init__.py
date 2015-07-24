import arrays
import files
import fits
import math
import printing
import xdebug
try:
    import ginga
except ImportError:
    print('-----------------------------------------------------------')
    print('-----------------------------------------------------------')
    print('WARNING: Not loading ginga model in xastropy.utils   \n Install ginga if you want it')
    print('-----------------------------------------------------------')
else:
    import xginga

