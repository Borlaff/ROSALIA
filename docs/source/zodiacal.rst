Zodiacal Light
==============

Modeling the Zodiacal light background with ROSALIA
---------------------------------------------------

The Zodiacal light is the dominant source of background for most Roman observations, except at the reddest filters (greater than about 2 micron) where the thermal background can become comparable. It arises from sunlight scattered by interplanetary dust in the Solar System. Its intensity depends on the sky position, observation time, filter bandpass, and the observer geometry.

ROSALIA models the Zodiacal light background through the ``rosalia.core.exposure.zodiacal()`` method. This method uses the exposure geometry and computes the expected background for each detector, then writes out a FITS product and a drizzled summary image.

The implementation supports three available backends:

- ``zodi_mode="stsci"``
- ``zodi_mode="zodipy"``
- ``zodi_mode="irsa"``

The notebook tutorial ``notebooks/tutorials/R2_Rosalia_zody_example.ipynb`` demonstrates the standard usage for all three options.


Quickstart example
------------------

The simplest workflow is to build a mock Roman/WFI exposure and then call the ``zodiacal`` method on it.

.. code-block:: python

    import rosalia as rs
    from astropy.time import Time

    ra = 180
    dec = 0
    PA = 9
    date = Time("2026-12-01T00:00:00")
    bandpass = "F129"
    exptime = 600

    observer = {
        "TELESCOP": "Roman/WFI",
        "pointing": [ra, dec],
        "FILTER": bandpass,
        "PA_Y": PA,
        "EXPSTART": date.mjd,
        "EXPTIME": exptime,
    }

    prefix = "zodi_test"
    custom_roman_exposure = rs.core.exposure(observer=observer, prefix=prefix)
    custom_roman_exposure.plot_footprint()
    zodi_out = custom_roman_exposure.zodiacal(zodi_mode="stsci")

This computes the Zodiacal background model and writes FITS products in the working directory or next to the input product.


Zodiacal model options
----------------------

The tutorial covers all three supported options for the Zodiacal background model.

STScI model
~~~~~~~~~~

.. code-block:: python

    zodi_out = custom_roman_exposure.zodiacal(zodi_mode="stsci")

This uses the STScI background model for Roman/WFI and is the direct workflow shown in the tutorial. This is the default option if an uknown string is entered for zodi_mode.

Zodipy model
~~~~~~~~~~~~

.. code-block:: python

    zodipy_zodi_out = custom_roman_exposure.zodiacal(zodi_mode="zodipy")

This is the fastest option in the tutorial and is based on the COSMOGLOBE/Zodipy package. This uses the extrapolate option and the 'dirbe' model, the only configuration that works at Roman's wavelength coverage.

IRSA model
~~~~~~~~~

.. code-block:: python

    irsa_zodi_out = custom_roman_exposure.zodiacal(zodi_mode="irsa")

This uses the IRSA/IPAC Zodiacal light model and is typically the slowest option, but it can be useful for comparisons with external reference models.


Method behavior
---------------

The method is implemented in ``rosalia.core.exposure`` as:

.. code-block:: python

    def zodiacal(self, zodi_mode="zodipy", verbose=False, output_name=None,
                 output_units="e/s", resolution=1, grid_method="random"):

The main arguments are:

- ``zodi_mode``: model backend. Supported values are ``"stsci"``, ``"zodipy"``, and ``"irsa"``.
- ``verbose``: enable more detailed progress output.
- ``output_name``: optional name for the output FITS file.
- ``output_units``: output units for the model, default is ``"e/s"``.
- ``resolution``: drizzled output resolution.
- ``grid_method``: sampling strategy used to evaluate the background field.

Internally, the method:

1. checks the number of science extensions in the exposure,
2. computes the background for each detector using ``rs.sky.get_zodiacal_background()``,
3. stores the detector-level arrays in a list,
4. writes a FITS file named using the exposure root name plus ``_zody.fits``,
5. creates a drizzled product named with ``_drz.fits``,
6. generates summary plots in the requested units.

If ``output_units == "e/s"``, it creates:

- an electron-per-second image plot,
- a surface-brightness plot in mag arcsec\ :sup:`-2`.

If another unit is entered, it generates a plot in MJy/sr instead.

The method returns a dictionary of the form:

.. code-block:: python

    {
        "zodi_list": zodiacal_background_list,
        "output_name": output_name,
    }


Outputs produced
----------------

The method writes output products in the same directory as the exposure or the local working directory when using a synthetic observer-only setup. In the tutorial workflow, this includes:

- a Zodiacal light FITS file with a ``_zody.fits`` suffix
- a drizzled product with a ``_drz.fits`` suffix
- generated diagnostic plots

These products are designed to be used as background model components for science analysis, correction, or comparison with measured sky levels.


Scientific context
------------------

The Zodiacal light background is one of the dominant sky components in Roman/WFI observations. Accurately modeling it is critical for:

- background subtraction,
- low-surface-brightness analysis,
- detector-level calibration,
- comparison of exposure geometry across observing conditions.

ROSALIA makes this model available in a consistent exposure-based workflow so that the same pointing, date, and filter setup can be evaluated across different Zodiacal-light backends and compared directly.

