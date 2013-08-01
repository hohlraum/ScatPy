*********************
Profiles and Defaults
*********************
.. currentmodule:: ScatPy

We have seen earlier that when you create new settings or targets using
``Settings()`` or ``Target()`` with no arguments ScatPy generates them
based on default files. This chapter will discuss the locations and contents
of those default. There are two distinct types of defaults used by ScatPy:

Profiles
    These are plain python scripts which tailor ScatPy's output for different
    platforms. For example, they allow you to generate jobs on a Windows computer
    that will be run on a remote cluster server.

Default Jobs
    These are ddscat.par files named ``default.par``. The default :class:`DDscat`,
    :class:`Settings` and :class:`target.Target` are all made to be like what is 
    defined in this file.

The Search Path
===============
The search scheme used to find a profile or default job is the same. The default
job is always named ``default.par``; profiles can have any name. The search
looks for the required file in three places:

1. In the current working directory
#. In the configuration folder in the user's home directory ``~/.ScatPy``
#. In the subfolder ``profiles/`` of the directory containing ``ScatPy.core``.
   Generally this will be something like ``site-packages/ScatPy/profiles


Switching Profiles
==================
The startup profile is ``default.py``, subject to the path resolution outline above.
Profiles can be switched on the fly with :func:`core.set_config(fname)`, where
``fname`` is the filename of the desired profile. This means that a job can be
setup and tested on a local machine with one configuration before changing profiles
to deploy and run it on a remote machine with a different setup.

Contents of Default Jobs
========================
The default job ``default.par`` is simply any valid ``ddscat.par`` file. If
the defined target requires it a ``shape.dat`` file should be included in
the same folder.

Contents of Profiles
====================
Here are the contents of a typical profile

.. literalinclude:: ../../../ScatPy/profiles/default.py

name
    A string describing the profile

os
    One of 'unix', 'mac', or 'windows' to describe the file system of the target
    computer. This determines things like the path separator to use.

mat_library
    The location of materials library on the target system

ddscat_path
    The path to the local copy of DDSCAT. Used by ``calculate()`` to find the
    DDSCAT executable

Here is a profile to generate jobs to be run on a SGE server.

.. literalinclude:: ../../../ScatPy/profiles/example_mpi.py


The main differences are that there is no ``ddscat_path`` since the jobs
are meant to be run remotely, and there is now a function ``write_script``.
This optional function will be run after calls to ``DDscat.write()``. Here it's
used to create a submission script that will queue the job to be run on 
a parallel SGE cluster.