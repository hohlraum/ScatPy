# -*- coding: utf-8 -*-
"""
The package configuration settings.

All paths should be specified using unix style forward slashes. These will
be written to file in the correct style.

"""
import ntpath
import posixpath

local_eslami={'mat_library':'h:/DDSCAT/material/',
        'ddscat_path':'c:/Program Files/ddscat/',
        'os':'windows',
        'path': ntpath}


#  landau eslami:
landau_eslami={
        'name':'landau',
        'mat_library':'/home/guests/eslami/mat_prop/',
        'ddscat_path':'/opt/local/lib/ddscat/src',
        'mpi_path': '/usr/mpi/gcc/openmpi-1.4.2/bin/',
        'os':'unix',
        'path': posixpath}

luna_eslami={
        'name': 'luna',
        'mat_library':'/home/guests/eslami/mat_prop',
        'ddscat_path':'/opt/local/lib/ddscat/src',
        'mpi_path': '/usr/mpi/gcc/openmpi-1.4.3/bin/',
        'os':'unix',
        'path': posixpath}


#  landau:
landau_mark={
        'name':'landau',
        'mat_library':'/home/guests/mark/mat_prop',
        'ddscat_path':'/opt/local/lib/ddscat/src',
        'mpi_path': '/usr/mpi/gcc/openmpi-1.4.2/bin/',
        'os':'unix',
        'path': posixpath}


luna_mark={
        'name': 'luna',
        'mat_library':'/home/guests/mark/mat_prop',
        'ddscat_path':'/opt/local/lib/ddscat/src',
        'mpi_path': '/usr/mpi/gcc/openmpi-1.4.3/bin/',
        'os':'unix',
        'path': posixpath}


#  local:
local_mark={'mat_library': '~/Documents/Analysis/ddscat/mp/',
            'ddscat_path':'/usr/local/bin/',
            'os':'unix',
            'path': posixpath}


#edit this to use the desired configuration
exec_settings=luna_mark


