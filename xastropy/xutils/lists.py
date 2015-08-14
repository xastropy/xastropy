"""
#;+ 
#; NAME:
#; lists
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for list utilities
#;   10-Dec-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os

from astropy import units as u
from astropy.table import QTable

from xastropy.xutils import xdebug as xdb

#
def dict_list_to_table(ldict):  
	''' Convert of list of dicts into a Table
	Parameters:
	-----------
	ldict: list  
	  List of dicts

	Returns:
	-----------
	table: Table  
	  Table constructed from the list
	'''
	# Init
	new_dict = {}
	keys = ldict[0].keys()
	for key in keys:
		new_dict[key] = []
	# Fill
	for idict in ldict:
		for key in idict.keys():
			new_dict[key].append(idict[key])
	# Table
	table = QTable(new_dict)
	# Return
	return table

