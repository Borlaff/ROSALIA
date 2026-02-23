How to use
==========

ROSALIA - From the command line
-------------------------------
Most users will utilize *ROSALIA* to generate sky-background models for images obtained with space telescopes. To make these tasks easier, we have developed a series of executables to analyze images directly from shell. Here we show the most common tools that will be used.

The primary tool is ``rosalia-sky``. ``rosalia-sky`` accepts a FITS or ASDF exposure image file and generates another image containing the expected background level per pixel. This program analyzes four different components of the sky background present for space telescopes, these are:

1. Zodiacal light
2. External stray light (i.e., stars, planets)
3. Diffracted light (PSF from bright stars)
4. Thermal self emission (also called internal stray light)

The expected brightness of these components vary with wavelength and observational conditions. Observations nearby a bright star might be dominated by stray-light, while in less crowded fields Zodiacal light will be the dominant contributor to the background flux. `Thermal self emission <https://roman.gsfc.nasa.gov/science/WFI_technical.html>`_ from the spacecraft starts to be comparable to the Zodiacal light the more we get into the infrared region of the electromagnetic spectrum. ``rosalia-sky`` generates one image per layer, allowing the user to combine them or subtract them from the real observed background.

Here is an example of how ``rosalia-sky`` works.

::

  rosalia-sky --all RST_WFI_ROSALIA_test_Orion_Belt_SCAWFI01.asdf


``rosalia-sky`` extracts all the necessary information from the exposure file metadata, finds the location of all-sky sources (stars, planets, Earth, Sun, the Moon), predicts the flux of stray-light that those sources will produce in the focal plane array, and finally, generates a series of ASDF and FITS files with the pixel-to-pixel flux (electrons per second) level expected for this particular exposure.

The ``rosalia-sky`` command line interface (CLI) has a number of arguments to support this functionality::

    rosalia-sky -h
    usage: rosalia-sky [-h] [--output OUTPUT] [--radius RADIUS]
                         [--g_mag_max G_MAG_MAX]
                         [--verbose VERBOSE]
                         input

    ROSALIA / rosalia-sky: Calculate stray-light level in
    Roman/WFI exposures

    positional arguments:
      input                 Input pattern that all asdf files
                            share. Example: input_pattern =
                            RST_WFI_SCA_*.asdf if your files are
                            RST_WFI_SCA_01.asdf,
                            RST_WFI_SCA_02.asdf,
                            RST_WFI_SCA_03.asdf [...]
                            RST_WFI_SCA_18.asdf.

    options:
      -h, --help            show this help message and exit
      --output OUTPUT       Name of the output FITS file
                            containing the stray-light level.
                            Default: "rosalia_stray_output.fits"
                            (default: default_rosalia-
                            stray_output.fits)
      --radius RADIUS       Maximum radius up to where all stars
                            are considered individually. The lower
                            the value, the faster the processing.
                            Default: 0.1 degree. (default: 0.1)
      --g_mag_max G_MAG_MAX
                            Maximum magnitude for stars considered
                            in the calculated. The lower the
                            value, the faster the processing.
                            Default: 15 mag. (default: 15)
      --verbose VERBOSE     Verbose option. Set True to see all
                            the information. (default: False)

    EXAMPLE: rosalia-stray input_image*.asdf


The mandatory input argument is the image filename. As of April 2025, Roman/WFI exposures are expected to be stored in individual ASDF files per detector (SCAs). To facilitate the analysis of complete exposures composed of multiple SCA images, ``rosalia-stray`` accepts `GNU shell patterns <https://www.gnu.org/software/findutils/manual/html_node/find_html/Shell-Pattern-Matching.html>`_,   (i.e., ``example_image_SCA*.asdf``). *ROSALIA* will identify and analyze all the files that correspond to the input pattern.

**Optional parameters**:

``--output``: [type: string] Use it to set a custom name for the sky background model.

``--radius``: [unit: degrees] *ROSALIA* automatically queries the Gaia, 2MASS, and WISE catalogs to construct an all-sky map of every single source in space, aiming to generate a stray-light model as precise as possible. However, for computational efficiency, *ROSALIA* assumes that stars at longer distances than **--radius** from the center of each detector do not need to be treated individually, and can be grouped together into *superstars* of 30 arcsec in size. Stars brighter than magnitude 6 (a total of 230) are not included in Gaia, so they are treated individually using `Hipparcos photometry <https://ui.adsabs.harvard.edu/abs/2016SPIE.9904E..2ES/abstract>`_. The higher the ``--radius`` value, the more stars are treated individually. Generally, 0.1 degrees is a minimum reasonable value. Higher quality models close to high stray-light contamination (i.e., an observation close to a bright group of stars) may require to set ``--radius`` to 1 degree or higher.

``--g_mag_max``: [unit: magnitudes] For efficiency, *ROSALIA* searches for stars up to a certain Gaia *g-band* magnitude. Stars with lower magnitudes than ``--g_mag_max`` are not considered in the stray-light analysis. The default value is *m=15* in the Gaia *g-band* AB system.

``--verbose``: [bool] Set True to return all the information about the internal steps made by *ROSALIA*. Useful for debugging.



ROSALIA - From Python
---------------------


To retrieve a list of random ingredients, you can use the ``rosalia.utils.sphere_dist()`` function:

.. autofunction:: rosalia.correct.
