# -*- coding: utf-8 -*-
"""
A set of tools for setting up and working with the program DDScat.


It includes a number of submodules:
  -core: the central definiton for a ddscat run including the run settings
  -config: configuration settings
  -targets: target definition classes
  -results: classes for manipulating the ddscat output files
  -ranges: classes for defining DDScat wavelength, size, and other ranges
  -utils: a handful of common utilities

@author: andrewmark
"""

import config
import fileio
import ranges
import results
import targets
import utils
from core import DDscat as DDscat
from core import Settings as Settings
from core import (pol_cL, pol_cR, pol_lV, pol_lH)

__version__=0.1
__all__=["DDscat", "ranges", "results", "targets", "config", "fileio", "utils",
         'pol_cL', 'pol_cR', 'pol_lV', 'pol_lH']
