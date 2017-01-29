"""
#;+ 
#; NAME:
#; build_casbah_galaxies
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for buildling CASBAH galaxy database
#;   02-Jan-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

from pyigm.field.igmfield import IgmGalaxyField

#from astropy import constants as const
from xastropy.casbah import utils as xcasbahu
from xastropy.xutils import lists as xxul

from xastropy.xutils import xdebug as xdb

# SDSS
def load_field(field):
    ''' Load up CASBAH data for a given field

    Parameters:
    -----------
    field: tuple
      (Name, ra, dec)

    Returns:
    --------
    lfield:     
      Loaded IgmGalaxyField class
    '''
    reload(xcasbahu)
    lfield = IgmGalaxyField((field[1],field[2]), name=field[0])

    # Load targets
    targ_file = xcasbahu.get_filename(field,'TARGETS')
    lfield.targets = Table.read(targ_file,delimiter='|',
        format='ascii.fixed_width', 
        fill_values=[('--','0','MASK_NAME')])

    # Load observing details for multi-object follow-up
    obs_file = xcasbahu.get_filename(field,'MULTI_OBJ')
    lfield.observing = Table.read(obs_file,delimiter='|',
        format='ascii.fixed_width', 
        fill_values=[('--','0','DATE_OBS','TEXP')])
    lfield.observing = Table(lfield.observing,masked=True) # Insist on Masked

    # Load galaxies
    sdss_file = xcasbahu.get_filename(field,'SDSS')
    sdss_tab = Table.read(sdss_file)
    # VSTACK
    hectospec_file = xcasbahu.get_filename(field,'HECTOSPEC')

    lfield.galaxies = sdss_tab

    # Return
    return lfield
