"""
#;+ 
#; NAME:
#; cgm.core 
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for core routines of CGM analysis
#;   29-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp
from astropy.io import fits, ascii
from astropy import units as u 
from astropy.coordinates import SkyCoord
from astropy import coordinates as coords
#from astropy import constants as const

#from xastropy.galaxy.core import Galaxy

from xastropy.atomic.elements import ELEMENTS
from xastropy.xutils import arrays as xu_array
from xastropy.obs import radec as xra

from astropy.utils.misc import isiterable

from xastropy.xutils import xdebug as xdb

# Path for xastropy
#xa_path = imp.find_module('xastropy')[1]

########################## ##########################
########################## ##########################
class IgmGalaxyField(object):
    """ Class for a field associating galaxies to the IGM/CGM
    """

    # Initialize 
    def __init__(self, name, radec, cosmo=None,
                 verbose=False):

        # Field
        self.name = name
        self.coord = xra.to_coord(radec)

        # Cosmology
        if cosmo is None:
            from astropy.cosmology import WMAP9 as cosmo
            if verbose is True:
                print('IgmGalxyField: Using WMAP9 cosmology')
        self.cosmo = cosmo

        # Init
        self.igm = None
        self.targets = None
        self.galaxies = None
        self.observing = None

    def get_observed(self,theta,subtab=None):
        '''Generate a Table of observed targets within an angular offset

        Parameters:
        ----------
        theta: Quantity
          Angular radius

        Returns:
        ------------
        obs_targ: Table
          Sub-table of targets that have been observed within this radius
        obs_dates: List
          List of observing dates [eventually might add to Table]
        '''
        if (self.targets is None) or (self.observing is None):
            raise ValueError('IgmGalaxyField: Need to fill the target and/or observing table first!')
        if subtab is None:
            # Trim on angular cut first
            targ_coord = SkyCoord(ra=self.targets['TARG_RA']*u.deg,
                dec=self.targets['TARG_DEC']*u.deg)
            sep = self.coord.separation(targ_coord)
            gdsep = np.where(sep < theta)[0]
            if len(gdsep) == 0:
                return None
            # Set all to False to start
            subtab = self.targets[gdsep]
        # Generate mask
        tmsk = np.array([False]*len(subtab))
        # Grab those with a MASK_NAME
        have_mask = np.where(subtab['MASK_NAME'].mask == False)[0]
        if len(have_mask) == 0:
            return None
        # Get unique mask values
        all_masks = subtab['MASK_NAME'][have_mask]
        uni_masks = np.unique(all_masks)
        obs_dict = {}
        # Loop on these
        for mask in uni_masks:
            obs_dates = self.get_mask_obsdate(mask)
            if len(obs_dates) > 0:
                mt2 = np.where(subtab['MASK_NAME'][have_mask]==mask)
                tmsk[have_mask[mt2]] = True
                obs_dict[mask] = obs_dates
        # Finish
        return subtab[tmsk], obs_dict

    def get_unobserved(self,theta):
        '''Generate a Table of observed targets within an angular offset

        Parameters:
        ----------
        theta: Quantity
          Angular radius

        Returns:
        ------------
        unobs_targ: Table
          Sub-table of targets that have been observed within this radius
        '''
        if (self.targets is None) or (self.observing is None):
            raise ValueError('IgmGalaxyField: Need to fill the target and/or observing table first!')
        # Trim on angular cut first
        targ_coord = SkyCoord(ra=self.targets['TARG_RA']*u.deg,
            dec=self.targets['TARG_DEC']*u.deg)
        sep = self.coord.separation(targ_coord)
        gdsep = np.where(sep < theta)[0]
        if len(gdsep) == 0:
            return None
        # Set all to False to start
        subtab = self.targets[gdsep]
        tmsk = np.array([True]*len(subtab))
        # Grab observed (short cut!)
        obs_tab, odict = self.get_observed(theta,subtab=subtab)
        # Remove those
        for kk,row in enumerate(subtab):
            if row['TARG_RA'] in obs_tab['TARG_RA']: # Could use DEC too
                tmsk[kk] = False
        # Return
        return subtab[tmsk]


    def get_mask_obsdate(self,mask_name):
        '''Given a mask name, find the observing dates
        Parameters:
        -----------
        mask_name: str
          Name of the mask

        Returns:
        -----------
        obs_dates: List
          List of the observing dates (can be empty)
        '''
        if self.observing is None:
            raise ValueError('Need to fill observing info!')
        #
        mt = np.where(self.observing['MASK_NAME'] == mask_name)[0]
        if self.observing['DATE_OBS'].mask[mt[0]]:
            return []
        obs_dates = [self.observing['DATE_OBS'][imt] for imt in mt]
        # Return
        return obs_dates

    #    
    def __repr__(self):
        return ('[{:s}: {:s} {:s} {:s}]'.format(
                self.__class__.__name__,
                 self.name,
                 self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True,alwayssign=True)))

