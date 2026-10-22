In-field Stars & PSF Modeling
==============================

Modelling scattered emission and diffraction from stars inside the Roman Field of View (WFI).

Overview
--------
The Nancy Grace Roman Space Telescope's Wide Field Instrument (WFI) features an 18-detector focal plane array with wide-angle diffraction spikes and extended point spread function (PSF) wings produced by the primary mirror aperture and support struts.

In crowded fields or observations containing bright stars ($m \lesssim 10$), in-field stellar emission contributes significant diffuse flux across large fractions of the detectors. *ROSALIA* provides analytical and numerical PSF modeling methods to calculate and subtract this contamination:

1. **SUPERBACK** (``method='superback'``, default):
   An analytical Fourier optics model for the Roman pupil geometry. Optimized with **PyTorch acceleration** for batched GPU and multi-core CPU execution.
2. **STPSF** (``method='stpsf'``):
   Numerical optical modeling utilizing `WebbPSF/STPSF <https://stpsf.readthedocs.io/>`_.
3. **GalSim** (``method='galsim'``):
   Parametric Moffat/Airy profile rendering.

---

PyTorch Acceleration for SUPERBACK
-----------------------------------

When ``method='superback'`` is selected, ROSALIA automatically leverages PyTorch for tensor-vectorized optical evaluations, star batching, and GPU acceleration.

Key Features:
~~~~~~~~~~~~~
- **Multi-Device Support**: Native execution on NVIDIA GPUs (CUDA), Apple Silicon GPUs (MPS), and multi-threaded CPUs.
- **Equal-Size Star Batching**: Stars are pre-sorted by their Moffat extension radius and partitioned into batches sharing identical stamp sizes ($S_{\text{batch}} \ge S_{\text{required}}$). This prevents wing truncation while fully utilizing GPU parallel cores.
- **Zero Gradient Tracking**: Executed under ``@torch.inference_mode()`` to eliminate computational graph overhead and conserve VRAM.
- **Out-of-Memory (OOM) Protection**: Extremely bright stars with gigantic stamps ($S > 2048$ pixels) stream anchor wavelength projections sequentially with eager garbage collection.

---

Quickstart: Using PyTorch with ``core.exposure.psf``
---------------------------------------------------

Default Execution (Automatic Hardware Detection)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
By default, ``exposure.psf()`` uses ``method='superback'`` and ``backend='pytorch'``, automatically detecting the fastest available hardware (CUDA -> MPS -> CPU):

.. code-block:: python

    import rosalia as rs

    # Create or load an exposure
    exposure = rs.core.exposure(filename="RST_WFI_example.fits")

    # Generate PSF models for in-field stars (auto-detects CUDA / MPS / CPU)
    detectors, wcs_list, canvas = exposure.psf(
        depth=27,          # Limiting surface brightness in mag arcsec^-2
        n_stars=50,        # Number of brightest stars to model (or None for all)
        verbose=True
    )

Targeting a Specific Device (``device`` parameter)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You can explicitly choose the hardware device using the ``device`` parameter:

**Apple Silicon GPU (MPS)**:

.. code-block:: python

    detectors, wcs_list, canvas = exposure.psf(
        device='mps',
        verbose=True
    )

**NVIDIA GPU (CUDA)**:

.. code-block:: python

    detectors, wcs_list, canvas = exposure.psf(
        device='cuda',
        verbose=True
    )

**Multi-core CPU**:

.. code-block:: python

    detectors, wcs_list, canvas = exposure.psf(
        device='cpu',
        verbose=True
    )

Selecting Precision (``dtype`` parameter)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- ``device='mps'``: Defaults to ``torch.float32`` (Apple Silicon MPS does not support float64).
- ``device='cuda'`` or ``device='cpu'``: Defaults to ``torch.float64`` for double-precision astronomy math.

You can explicitly set precision:

.. code-block:: python

    import torch

    detectors, wcs_list, canvas = exposure.psf(
        device='cpu',
        dtype=torch.float32  # Cuts CPU memory and runtime in half
    )

Switching Back to Legacy NumPy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To use the legacy NumPy/SciPy implementation without PyTorch:

.. code-block:: python

    detectors, wcs_list, canvas = exposure.psf(
        method='superback',
        backend='numpy',
        n_workers=4
    )

---

Parameters of ``core.exposure.psf``
-----------------------------------

.. list-table::
   :widths: 20 15 65
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``depth``
     - ``float``
     - Limiting surface brightness depth in $\text{mag}/\text{arcsec}^2$ (default: 31). Controls star extension radius.
   * - ``catalog``
     - ``DataFrame``
     - Optional custom star catalog. If ``None``, queries Gaia/2MASS/WISE automatically.
   * - ``method``
     - ``str``
     - PSF modeling engine: ``'superback'`` (default), ``'stpsf'``, or ``'galsim'``.
   * - ``backend``
     - ``str``
     - Backend for ``superback``: ``'pytorch'`` (default) or ``'numpy'``.
   * - ``device``
     - ``str``
     - PyTorch hardware device: ``'auto'`` (default), ``'cuda'``, ``'mps'``, or ``'cpu'``.
   * - ``dtype``
     - ``torch.dtype``
     - Tensor precision: ``torch.float32`` or ``torch.float64``.
   * - ``n_stars``
     - ``int``
     - Maximum number of brightest stars to model (default: ``None``, all in-footprint stars).
   * - ``verbose``
     - ``bool``
     - If ``True``, prints real-time progress bars and device information.

---

Stand-Alone Chromatic PSF Generation
------------------------------------
You can also generate an individual PSF FITS PrimaryHDU directly:

.. code-block:: python

    import rosalia as rs

    psf_hdu = rs.psf.superback_roman_psf(
        position=(1024, 1024),
        SCA=1,
        bandpass="F158",
        shape=(512, 512),
        device="mps"  # or 'cuda', 'cpu', 'auto'
    )

    # Save to FITS
    psf_hdu.writeto("roman_wfi_psf.fits", overwrite=True)
