***************
DDSCAT Settings
***************

Aside from the target definition, the various parameters which describe
a DDSCAT run found in a ddscat.par file are stored in ScatPy in a
:class:`core.Settings` object.

The Simplest Case
=================
Like targets, an empty call to :class:`core.Settings` yields a new ``Settings``
object based on the default:

>>> s = Settings()
>>> print s
<ScatPy.core.Settings object at 0x10a544550>


Settings Fields
===============
Most of the fields contained in the ``Settings`` object have names based on
the ones found in ddscat.par and the DDSCAT userguide. A complete list of the
fields, and the values they can take is found in the `ScatPy API`_ under 
:class:`core.Settings`. A few common ones, with examples of how to modify
them are::

    >>> s.NAMBIENT
    1.0
    >>> s.NAMBIENT = 1.33                   # Change medium to water
    >>>
    >>> s.NRFLD = True                      # Make nearfield calc
    >>> s.NRFLD_EXT = np.array((0.5,)*6)    # Set nearfield box to 0.5x on all sides
    >>>
    >>> s.IORTH = True                      # Do calculation of 2nd polarization
    >>> s.Epol = np.array((0, 1, 1))        # 1st polarization linearly polarized @ 45deg

Alternatively, fields can be modified on creation by passing new field values
by keyword:

>>> default = Settings()
>>> default.S_INDICES                          # Default value
[11, 12, 21, 22, 31, 41]
>>> modified = Settings(S_INDICES = [22,31])   # Modified during creation
>>> modified.S_INDICES
[22, 31]

Polarization
^^^^^^^^^^^^
ScatPy includes named polarizations that can be used to simplify polarization
specification

>>> s.Epol = pol_lH     # Linear horizontal
>>> s.Epol = pol_vH     # Linear vertical
>>> s.Epol = pol_cL     # Circular left
>>> s.Epol = pol_cR     # Circular right

The sense of circular polarization is defined in :mod:`core`. There are
hardcoded definitions for both the convention used by spectroscopists and
that documented in the DDSCAT user guide. The default is to use the spectroscopist's
convention.

Ranges
^^^^^^
ScatPy defines classes for managing ranges used to specify things like
wavelengths, scattering planes and target orientations. These are found in
the module :mod:`ranges`. The most basic type is a :class:`ranges.How_Range`,
used for defining the range of wavelengths over which to do calculations. The
final string parameter defines how to space the steps

>>> wave = ranges.How_Range(0.200, 1.100, 10, 'LIN')
>>> print wave
0.200000  1.100000  10  LIN
>>> wave.last = 1.5           # Change the final wavelength
>>> print wave
0.200000  1.500000  10  LIN

``How_Range`` can be used as an iterator

>>> [round(w,2) for w in wave]
[0.2, 0.34, 0.49, 0.63, 0.78, 0.92, 1.07, 1.21, 1.36, 1.5]

Other ranges are :class:`ranges.Lin_Range`, :class:`ranges.Scat_Range`,
:class:`ranges.Scat_Range_1dPBC`, and :class:`ranges.Scat_Range_2dPBC`.

Defining Ranges of Target Sizes
===============================
DDSCAT allows a many target sizes to be calculated in a single job. This
is done by specifying a range for ``aeff``. In ScatPy, targets have a defined
physical size, so multiple target size calculations are accomplised by scaling
that size with the :class:`Settings` field ``scale_range``. A single numeric value
for ``scale_range`` instructs DDSCAT to do a single size calculation with the
target geometry scaled by that value. Providing a range for ``scale_range`` will
do calculations with target geometries scaled by all of the values in the range.
So with this target (a 1um diamsphere)

>>> t = targets.ELLIPSOID((0.5, 0.5, 0.5))

This will do a calculation on a single 1um sphere:

>>> s.scale_range = 1.0

This will do a calculation on a single 2um sphere:

>>> s.scale_range = 2.0

Both calculations will be done with the same number of dipoles, but with sparser
dipole density for the second. The following will do calculations over spheres
with the diameters (1.0, 1.5, 2.0):

>>> s.scale_range = ranges.How_Range(1, 2, 3)

