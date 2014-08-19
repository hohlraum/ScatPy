************
User's Guide
************

Introduction
============
This user's guide is intended to help new users become proficient at using
ScatPy to manage their work with DDSCAT. It assumes you have a working copy of
DDSCAT and some understanding of how it works. 

DDSCAT
======
DDSCAT is an implementation of the discrete dipole approximation by Bruce T. Draine
and Piotr J. Flatau. The DDSCAT homepage is at:
* <http://www.astro.princeton.edu/~draine/DDSCAT.html>_

The latest version of the source code is available at:
* https://code.google.com/p/ddscat/

The authors have kindly provided a precompiled binary for Windows users. For 
other operating systems you will have to compile from source.

The indispensable DDSCAT user guide can be downloaded from:
* https://ddscat.googlecode.com/files/ddscat7.3.0_UserGuide_130529.pdf


Overview
========
The general workflow for ScatPy is to:

1. Create a target
#. Create a DDscat job
#. Select the desired run parameters for job (wavelength range, solver,)
#. Run the job, either locally or on a remote server
#. Load the result tables
#. Process and plot the results

Profiles
========
When ScatPy is imported for the first time it will attempt to create a folder ``.ScatPy``
in the current user's home directory. It will populate this folder with configuration files
that can be edited to suit your requirements. The two most important ones are ``default.par``
which is the base ``DDSCAT.par`` file for all DDSCAT jobs, and ``default.py`` which
contains the default configurations for file locations, and file system types. These files
can be edited to match your preferred configuration.

Multiple profile files mean that it is possible to setup and test a job on one computer 
(e.g. a local workstation) and then export the job settings for execution on a different
system (e.g. a remote cluster). 


Userguide Contents
==================
.. toctree::
    :maxdepth: 2
    
    targets_ug.rst
    settings_ug.rst
    jobs_ug.rst
    table_ug.rst
    collections_ug.rst
