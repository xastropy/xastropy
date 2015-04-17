# Non-dependent modules
import data

# Dependent modules
try:
    import pymc
except ImportError:
    print('-----------------------------------------------------------')
    print('-----------------------------------------------------------')
    print('WARNING: Not loading model in xastropy.fN   \n Install pymc if you want it')
    print('-----------------------------------------------------------')
else:
    import model
