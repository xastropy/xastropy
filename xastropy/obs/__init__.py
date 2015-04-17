# Dependent modules
try:
    import astroquery
except ImportError:
    print('-----------------------------------------------------------')
    print('-----------------------------------------------------------')
    print('WARNING: Not loading modules in xastropy.obs except radec.  \n Install astroquery if you want them')
    print('-----------------------------------------------------------')
else:
    import finder
    import keck
    import lick
    import x_getsdssimg

# Non-dependent modules
import radec
