Zodiacal Light
=========

Modeling the Zodiacal light background with ROSALIA
-----------------------------

*This is a work in progress section. We are currently (April 16, 2026) implementing a new Zodiacal light model from STScI team in place. However, the currently in place model is usable an it will provide reasonable solutions for most uses*.

The Zodiacal light is the dominant source of background for most of the observations with the Roman Space Telescope, with exception of the thermal background at the reddest filters. It is produced by the scattering of sunlight by interplanetary dust particles in the Solar System. The intensity of the Zodiacal light depends on the line of sight through the Solar System, and therefore on the position of the target in the sky, the time of the observation, and the wavelength. 

There are several models to estimate the Zodiacal light background, with different levels of complexity and accuracy. ROSALIA implements the model by  `COSMOGLOBE / Zodipy <https://cosmoglobe.github.io/zodipy/usage/>`_. If you use this model, please cite the following papers: `(San et al. 2022) <https://arxiv.org/abs/2205.12962>`_  and  `(San 2024) <https://joss.theoj.org/papers/10.21105/joss.06648#>`_.

The following code block shows a quick example to easily model the Zodiacal light with ROSALIA. 

.. code-block:: python

    # This notebook demonstrates how to use ROSALIA to measure the expected zodiacal background in a simple image. 
    import rosalia as rs
    from astropy.time import Time
    ra = 123 # Right ascension, in degrees. 
    dec = 45 # Declination, in degrees.
    PA = 32 # Position angle, in degrees.
    date = Time("2026-10-01T00:00:00") # Date of the observation, in Astropy Time YYYY-MM-DDTHH:MM:SS format.
    bandpass = "F129"
    exptime = 600 # Exposure time, in seconds.

    rosalia_zody = rs.correct.rosalia_zody(ra=ra, dec=dec, PA=PA, date=date, 
                                     bandpass=bandpass, exptime=exptime,
                                     verbose=False)


.. image:: ../../notebooks/tutorials/WFI_F129_RA_123.000_DEC_045.000_MJD_61314.00000_PA_067.00_zody_drz_scaled_fe2mu.png
  :alt: ROSALIA example of Zodiacal light model.

