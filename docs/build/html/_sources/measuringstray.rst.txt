Stray Light
=========

In-field vs. out-field stray-light
-----------------------------

Stray-light is the light that reaches the detector without following the nominal optical path. It can be produced by a variety of sources, like bright stars, planets, the Moon, or even the Earth. On first approximation, stray-light can be classified into two main categories: **in-field** and **out-field** stray-light.



Quickstart: Stray Light Analysis Example with ROSALIA
-----------------------------

The following Jupyter-notebook shows a quick example to estimate stray-light with ROSALIA. 

.. code-block:: python

    import rosalia as rs
    # Let's define the minimum parameters to generate a dummy Roman / WFI image

    # To get some challenging environment, let's target the Pleiades.
    ra = 56.6583333  # Right ascension, in degrees. 
    dec = +24.1780556 # Declination, in degrees.
    PA = 0 # Position angle, in degrees.
    
    from astropy.time import Time
    date = Time("2026-06-01T00:00:00") # Date of the observation, in Astropy Time YYYY-MM-DDTHH:MM:SS format.
    bandpass = "F129"
    exptime = 600 # Exposure time, in seconds.

    rosalia_stray = rs.correct.rosalia_stray(ra=ra, dec=dec, PA=PA, date=date, bandpass=bandpass, exptime=exptime, radius=1, g_mag_max=15, sun_block=False, verbose=False, catalog=None)

.. image:: ../../images/rosalia_loading.gif
  :width: 802
  :height: 512
  :alt: ROSALIA Stray-light loading bar

The notebook is available in the ``notebooks`` directory of the ROSALIA repository, and can be downloaded from the link below:
* :download:`R1_Rosalia_stray_example.ipynb <../../notebooks/R1_Rosalia_stray_example.ipynb>`. To run this notebook locally, make sure you have Jupyter installed and then:

.. code-block:: bash
    jupyter notebook notebooks/R1_Rosalia_stray_example.ipynb


Case study: Calibrating stray-light with a bright star
-----------------------------

A particular case of interest would be to produce an intense stray-light background signal on purpose, to calibrate the models against real observations. To gather sufficient signal-to-noise and minimize contamination from other background sources like the Zodiacal light, the selected star (*offending* source) would need to be as bright as possible. In addition, we need to place it in positions around Roman Space Telescope WFI that are known to produce the maximum stray light if hit by a source. With ROSALIA, it is possible to find out where should we point our telescope and which roll (position) angle we should have to maximize or minimize the stray-light from a star. 

First, let's visualize the Normalized Detector Irradiance (NDI) map for Roman Space Telescope Wide Field Instrument. 

.. image:: ../../notebooks/ndi_2deg_CARs_no_points.png
    :width: 49 %
.. image:: ../../notebooks/ndi_20deg_CARs_no_points.png
    :width: 49 %

The NDI is defined as the ratio of the stray-light irradiance (power per unit area) received at the detector to the irradiance of the source at the entrance of the telescope `(see Bely 2003) <https://ui.adsabs.harvard.edu/abs/2003dclo.book.....B/abstract>`_ . This function can be used to estimate the flux of photoelectrons that an off-axis source will generate in a certain region of the detector. Similarly to the Point Spread Function (PSF) defines the spread of point sources by diffraction effects, the sharpness of the NDI defines the effectiveness of the telescope's baffling. This function is strongly dependent on the optical setup and wavelength. For a given telescope, the NDI depends on the angular distance between the optical axis and the source, the position angle of the source in the focal plane reference frame, the observation wavelength, and the position on the FOV (x,y),

We can identify several points of interested were stray-light is expected to peak if illuminated: 

.. image:: ../../notebooks/ndi_2deg_CARs.png
    :width: 49 %
.. image:: ../../notebooks/ndi_20deg_CARs.png
    :width: 49 %


Now that we have found out where which orientations - defined as X and Y angular offsets from the center of the Roman WFI focal plane - should we use to produce stray-light, let's find out where should we point Roman/WFI to get a particular target oriented to those locations. This can be done using *rs.telescopes.Roman.find_wfi_center_for_offset_targets*:

.. code-block:: python

    # These are the coordinates of specially sensitive locations 
    # for stray-light, measured in offset degrees from the center of WFI focal plane
    dX = np.array([-3.7, 3.7, 6.05, 2.45, -2.45, -6.05, -2, 2, -1.3, 1.3, -0.0688, 0.0688, -0.138, 0.138])
    dY = np.array([5.27, 5.27, 1.18, -5.07, -5.07, 1.18, 6.15, 6.15, 1, 1, -0.0312, -0.0312, -0.0035, -0.0048])
    number_of_points_of_interest = len(dX)
    # To measure the stray-light from those sources, we need to place the center 
    # of Roman / WFI at a certain distance and position angle from the source 
    # that generates the stray-light. We will name that source the "offending" source. 

    # ROSALIA has a specific tool to compute those locations. 

    # The optimal position angle of Roman depends with time. 
    # Let's set an approximate time for the Commissioning Activities
    from astropy.time import Time
    date = Time('2026-11-21T00:00:00.0', format='isot', scale='utc')

    ra_wfi = np.zeros(number_of_points_of_interest)
    dec_wfi = np.zeros(number_of_points_of_interest)
    pa_wfi = np.zeros(number_of_points_of_interest)

    for i in range(number_of_points_of_interest):
        offset_pointing = rs.telescopes.Roman.find_wfi_center_for_offset_target(ra_target=target.ra.degree,
                                                          dec_target=target.dec.degree,
                                                          mjd=date.mjd,
                                                          dX=dX[i], dY=dY[i])
        ra_wfi[i] = offset_pointing["ra_wficen"]
        dec_wfi[i] = offset_pointing["dec_wficen"]
        pa_wfi[i] = offset_pointing["pa_wfi"]

Now that we have the exposure parameters for WFI, we only need to simulate the stray-light with *rosalia_stray*: 

.. code-block:: python

    for i in range(number_of_points_of_interest):
        ra = ra_wfi[i]
        dec = dec_wfi[i]
        pa = pa_wfi[i]
        mjd = mjd
        bandpass="F129"
        exptime=600
        rosalia_stray = rs.correct.rosalia_stray(ra=ra, dec=dec, 
                                                 PA=pa, date=date, bandpass=bandpass, 
                                                 exptime=exptime, radius=1,
                                                 g_mag_max=15)




The following Jupyter-notebook shows an example of how to do that with ROSALIA. The notebook is available in the ``notebooks`` directory of the ROSALIA repository, and can be downloaded from the link below:
