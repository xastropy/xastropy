from . import arrays
from . import files
from . import fits
from . import math
from . import printing
from . import xdebug
try:
    import ginga
except ImportError:
    print('-----------------------------------------------------------')
    print('-----------------------------------------------------------')
    print('WARNING: Not loading ginga model in xastropy.utils   \n Install ginga if you want it')
    print('-----------------------------------------------------------')
else:
    pass
    #import xginga

