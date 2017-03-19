# Dependent modules

try:
    import astroquery, PIL
except ImportError:
    print('-----------------------------------------------------------')
    print('-----------------------------------------------------------')
    print('WARNING: Not loading modules in xastropy.obs except radec.  \n Install astroquery and PIL if you want them')
    print('-----------------------------------------------------------')
else:
    from . import finder
    from . import keck
    from . import lick
    from . import x_getsdssimg

# Non-dependent modules
from . import radec
