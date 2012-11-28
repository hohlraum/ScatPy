# -*- coding: utf-8 -*-
"""
A set of tools for setting up and working with the program DDScat.

It includes a number of submodules:
    core: the central definiton for a ddscat run including the run settings
    config: configuration settings
    targets: target definition classes
    results: classes for manipulating the ddscat output files
    ranges: classes for defining DDScat wavelength, size, and other ranges
    utils: a handful of common utilities

@author: andrewmark
"""

import config
import fileio
import ranges
import results
import targets
from core import DDscat as DDscat
from core import Settings as Settings

__version__=0.922
__all__=["DDscat", "ranges", "results", "targets", "config", "fileio", "utils"]



"""
Version Log.
Another test message. #2

0.924 
    -Started to use git for versioning controls

0.923 Nov 27th, 2012
    -Added utility to mix materials for alloy
    -Added save method to MInTable

0.922 Nov 19th, 2012
  -Added Target CYLINDER contributed by Sahand
  -Added luna and landau profiles for Sahand to main "trunk"
  -Changed io module name to fileio to avoid clashes with builtin module io

0.921 Oct 2nd, 2012
  -Added profile for luna
  -Added additional profile field for mpi_path
  -Modified sge_write to write to use mpi_paths





"""