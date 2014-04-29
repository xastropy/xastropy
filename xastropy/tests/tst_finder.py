# Tests the Finder chart program
#   Requires "requests", "astropy", "aplpy", "PIL"
from xastropy.obs import x_finder as xf
from xastropy.obs import x_radec as x_r
reload(x_r)
reload(xf)
xf.main(['TST', '10:31:38.87', '+25:59:02.3'],radec=1, SDSS=1, imsize=8.)#,BW=1)
exit()
