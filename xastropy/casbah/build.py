"""
#;+ 
#; NAME:
#; build
#;    Version 1.0
#;
#; PURPOSE:
#;    Simple script for running CASBAH builds
#;   02-Jan-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os
import warnings

from astropy import units as u

from xastropy.casbah import galaxies as xcasg
from xastropy.casbah import survey_figs as xcassf

from xastropy.xutils import xdebug as xdb

fields = [ ('PG1630+377', 248.00461542*u.deg, 37.63055111*u.deg),
           ('PG1522+101', 231.1023075*u.deg, 9.9749372*u.deg),
           ('PG1206+459', 182.24172042*u.deg, 45.67652417*u.deg),
           #('PG1148+549', 177.83526042*u.deg, 54.625855*u.deg),
           ('FBQS0751+2919', 117.80127583*u.deg, 29.3272975*u.deg),
           ('PG1407+265', 212.349634*u.deg, 26.3058650*u.deg),
          ]

warnings.warn('casbah.build: Deal with Cosmology')
# Galaxies first
path = os.getenv('CASBAH_GALAXIES')
skip_SDSS = True
mk_figs = False
for field in fields:

    # Targeting for other Telescopes
    xcasg.build_targets(field, path=path)

    # Target Imaging
    xcasg.build_imaging(field)

    # Target figures
    if mk_figs:
        xcassf.hectospec_targets(field)
        xcassf.deimos_targets(field)

    # Spectra+Redshifts
    xcasg.build_spectra(field)

    # SDSS/BOSS
    if not skip_SDSS:
        xcasg.build_sdss(radius=0.2*u.deg)



