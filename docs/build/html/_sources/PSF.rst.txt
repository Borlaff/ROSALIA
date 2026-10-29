In-field Stars
==============

Modelling scattered emission (PSF) from stars inside the Roman field of view
-----------------------------------------------------------------------------

ROSALIA can model the contribution from stars inside the Roman/WFI field of view using the ``rosalia.core.exposure.psf()`` method. This is useful when you want to estimate the stellar PSF contribution across the detector footprint, either for science analysis or for subtracting in-field stellar contamination from an exposure.

The method constructs a PSF model for stars around the FOV using a source catalog and a selected PSF-generation method. By default, the method uses the `superback <https://science.nasa.gov/mission/roman-space-telescope/superback-building-an-empirically-motivated-background-model-to-optimize-roman-observations/>`_ model, and it returns the detector-by-detector PSF arrays, their WCS information, and the combined canvas model.


Quickstart example
------------------

The notebook tutorial ``notebooks/tutorials/R3_Rosalia_psf_example.ipynb`` demonstrates the simplest workflow.

.. code-block:: python

    import rosalia as rs
    from astropy.time import Time

    # Define the target
    target = 'alpha Centauri'
    ra = 239.92041
    dec = -60.835148
    PA = 30
    date = Time("2027-01-01T00:00:00")
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

    prefix = "alphacen"
    custom_roman_exposure = rs.core.exposure(observer=observer, prefix=prefix)
    psf_out = custom_roman_exposure.psf()

This produces a PSF model for the exposure. The returned object is a tuple of the form:

.. code-block:: python

    detectors, wcs_list, canvas = custom_roman_exposure.psf()

The first item contains the 18 detector-level PSF arrays, the second item contains their WCS metadata, and the third item is the global PSF canvas for the full WFI footprint.


Generating a mosaic from the model
----------------------------------

Once the PSF model is generated, it can be assembled into a mosaic for inspection or output.

.. code-block:: python

    mos_data, mos_wcs = rs.mosaics.generate_mosaic(psf_out[0], psf_out[1])
    rs.utils.save_fits(mos_data, "test_psf.fits", mos_wcs)

This is the pattern shown in the tutorial notebook, and it is useful for validating the PSF model or writing a combined FITS product for further analysis.


Method behaviour
----------------

The ``rosalia.core.exposure.psf()`` method is defined as:

.. code-block:: python

    def psf(self, g_mag_max=15, depth=31, catalog=None, method='superback', n_stars=None, n_workers=None, verbose=False):

The main parameters are:

- ``g_mag_max``: maximum Gaia magnitude used to select stars in the simulation.
- ``depth``: depth of the PSF model, expressed in a model-specific way by the underlying PSF implementation.
- ``catalog``: optional user-supplied catalog of stars. If not supplied, the method builds one automatically.
- ``method``: PSF modeling method. The default is ``'superback'``.
- ``n_stars``: optional limit on the number of stars included.
- ``n_workers``: number of worker processes used for parallel evaluation.
- ``verbose``: enable additional progress output.

The method selects stars that fall within the Roman/WFI footprint, maps them onto detector coordinates, and builds a PSF model across the full focal plane. It then reconstructs detector-level PSF arrays and their WCS objects.


Supported modeling approach
---------------------------

The code path used by the method currently supports a PSF model built around the ``superback`` implementation, with other options available in the lower-level PSF tooling. In the project’s current implementation, the higher-level exposure method defaults to ``method='superback'``.

This is appropriate when the goal is to model in-field stellar contamination across the Roman detector array in a compact and reproducible way.


Notes for use
-------------

- This workflow is intended for Roman/WFI exposures.
- The method uses source catalog information over the field of view.
- The output is not a single image product by itself; it is a model object built from the detector-level PSF arrays and WCS metadata.
- For downstream work, the reconstructed PSF canvas or detector mosaics can be saved to FITS for inspection.

In short, the feature is designed to turn a Roman exposure into a modeled stellar PSF field that can then be inspected, converted into a mosaic, or used as input for subsequent calibration and contamination analysis.