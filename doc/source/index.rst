**********************************
Welcome to ScatPy's documentation!
**********************************

Introduction
============

ScatPy is a Python package for interfacing to the popular scattering simulator
`DDSCAT <http://www.astro.princeton.edu/~draine/DDSCAT.html>`_. ScatPy provides a rich toolset to:

* Create standard DDSCAT scattering targets based on physical (rather than dipole) dimensions
* Construct and visualize complex custom scattering targets
* Manage the job parameters found in the ddscat.par file
* Organize iterative jobs requiring multiple targets or input parameters
* Script job submission to cluster queue managers
* Maintain profiles and defaults for deployment on platforms other than the local machine
* Load, plot and manipulate DDSCAT output tables
* Manage the output from multiple jobs through results collections
* Work with and visualize nearfield results as multidimensional numpy arrays
* Suitable for interactive or scripted use

Gallery
=======

.. image:: /Gallery.png
    :target: examples.html


A Simple Example
================

An example

.. plot:: ../examples/example1.py
   :include-source:


Contents
========

.. toctree::
    :maxdepth: 2

    userguide/userguide.rst
    examples.rst
    versions.rst
    api/api.rst

Installation
============

1. Start by getting a working version of `Python <http://www.python.org/getit/>`_.
    For Windows users this is most easily done by installing
    `Python (x,y) <https://code.google.com/p/pythonxy/>`_\ . 

#. ScatPy has the following mandatory dependency (it's installed as standard by
    Python(x,y)):

    * `numpy <http://www.numpy.org/.>`_
    * `matplotlib <http://matplotlib.org/index.html>`_
    * `scipy <http://www.scipy.org/>`_

    Visualizing targets in 3D also requires: 

    * `mayavi <http://code.enthought.com/projects/mayavi/>`_

#. Download ScatPy from `PyPi <https://pypi.python.org/pypi/ScatPy>`_

#. Unzip the package and run ``python setup.py install`` from the command line.

#. Start Python and type ``import ScatPy``. If you don't recieve any import
   warnings then you're ready to go.

Getting Started
===============
Once you have ``ScatPy`` installed, start by reading the :doc:`userguide/userguide`
and reviewing the :doc:`examples`. More advanced users will benefit from referring
to the :doc:`api/api`. 

Acknowledgements and Citing
===========================
Work on ScatPy was carried out at the Max Planck Institute for Intelligent Systems, under funding provided by the DPG.

If this code contributes to a publication, please cite::

    A.G. Mark, S. Eslami, and P. Fischer; Python code for interacting with DDSCAT.
    Available at https://pypi.python.org/pypi/ScatPy.


License
=======
ScatPy is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

