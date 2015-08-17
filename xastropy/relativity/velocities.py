# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Convenience functions for `astropy.relativity`.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import warnings
import numpy as np

# Bring these back when pushing
from astropy import constants as const
from astropy import units as u

# Remove these when pushing
#from astropy import constants as const
#from astropy import units as u


__all__ = ['v_from_z']

__doctest_requires__ = {'*': ['scipy.integrate']}


# #####
def v_from_z(z1, z2):
    """ Find the relativistic velocity between 2 redshifts.

    Parameters
    ----------
    z1 : float
       One redshift.  
    z2 : float or array
       Other redshift(s)

    Returns
    -------
    v : Quantity (km/s)
      Velocity

    Notes
    -----
    """
    R = (1+z1) / (1+z2)
    v = const.c * (R**2 - 1)/(1+R**2) 

    return v.to('km/s')

# #####
def z_from_v(z, v):
    """ Find the redshift given z and v

    Parameters
    ----------
    z : float
       Redshift
    v : float or array
       Velocities

    Returns
    -------
    z : float or array
      New redshifts

    Notes:
    -----
    """
    # Check for unit
    if not isinstance(v,u.quantity.Quantity):
        # Assume km/s
        v = v * u.Unit('km/s')
    
    # b
    bval = (v/const.c.to('km/s'))
    
    # R
    R = np.sqrt((1-bval)/(1+bval))
    # Finally
    znew = (1+z)/R - 1

    return znew.value


# ###################
# TESTING
if __name__ == '__main__':

    from astropy.relativity import velocities as arv
    # v_from_z
    z1=2.0
    z2=2.01
    v = arv.v_from_z(z1,z2)
    print('v_from_z: v = {:g}, for (z1,z2) = ({:g},{:g})'.format(v,z1,z2))

    # z_from_v
    z1 = arv.z_from_v(z2,v)
    print('z_from_v: v = {:g}, for (z1,z2) = ({:g},{:g})'.format(v,z1,z2))
    twov = np.array([-300., -3000.])
    z2 = arv.z_from_v(z1,twov)
    print(z2)
