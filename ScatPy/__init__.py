# -*- coding: utf-8 -*-
"""
ScatPy is a set of tools for setting up DDSCAT jobs and analaysing the results.


It includes a number of submodules:
  -core: the central definiton for a ddscat run including the run settings
  -targets: target definition classes
  -results: classes for manipulating the ddscat output files
  -ranges: classes for defining DDScat wavelength, size, and other ranges
  -utils: a handful of common utilities

"""
import core
import utils

# Create profile folder in user home directory if one does not already exist.
try:
    utils.make_profile()
except IOError:
    pass

import fileio
import ranges
import results
import targets

from core import (DDscat, Settings, set_config)
from core import (pol_cL, pol_cR, pol_lV, pol_lH)

__all__=["DDscat", "Settings", "set_config", "ranges", "results", "targets", "fileio", "utils",
         'pol_cL', 'pol_cR', 'pol_lV', 'pol_lH']

try:
    from _version import __version__ as v
    __version__ = v
    del v
except ImportError:
    __version__ = "UNKNOWN"
    
