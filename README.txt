
**************************************
ScatPy -- A Python interface to DDSCAT
**************************************


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

Documentation
=============

Complete documentation can be found at:
    http://pythonhosted.org/ScatPy


Download
========

The package can be downloaded for installation via easy_install at
    https://pypi.python.org/pypi/ScatPy


Example
=======

.. code:: python

    from ScatPy import *

    # Establish target geometry (in um)
    length = 0.100
    radius = 0.020
    target = targets.CYLNDRCAP(length, radius, d=0.005, material='Au_Palik.txt')

    # Create a job to be run in the subdirectory tmp/
    job = DDscat(folder = './tmp', target=target)

    # Change the range of calculated wavelengths and ambient index
    job.settings.wavelengths = ranges.How_Range(0.300, 0.600, 15)
    job.settings.NAMBIENT = 1.0

    # Run the job locally
    job.calculate()

    # Open the results qtable, plot Q_sca, and Q_abs, and add a legend
    ans = results.QTable(folder = './tmp')
    ax = ans.plot(['Q_sca', 'Q_abs'])
    ax.legend(loc=0)