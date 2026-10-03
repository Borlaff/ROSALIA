# Alejandro S. Borlaff. NASA Ames Research Center. a.s.borlaff@nasa.gov / asborlaff@gmail.com
# January 20, 2023.
#
# STRAYCOR/PSF module
# This module will hold all the programs related to the modelling and removal
# of the PSF.
#
# Version log:
# v.1.0 - 20 Enero 2023. First loading of programs inherited from former monolithic straycor.py
#
##########################################################

############################
import os
import pandas as pd
import numpy as np
import bottleneck as bn
from tqdm import tqdm
from astropy.io import fits
import astropy.wcs as astropy_wcs
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.interpolate import interp1d
import rosalia as rs

# Suppress warnings. Comment this out if you wish to see the warning messages
import warnings
warnings.filterwarnings('ignore')


#%%
class Star:
    """Represent a catalog star with broadband photometry and an interpolated SED.

    The input catalog row is expected to contain ``ra``, ``dec``, ``source_id``,
    and the AB-magnitude columns consumed by :meth:`set_magnitudes` for Gaia,
    2MASS, and WISE. If ``detector_id`` is present, it is also stored as
    ``WFI``. Photometry is associated with fixed representative wavelengths;
    :meth:`sed` linearly interpolates (and extrapolates beyond) those samples
    in AB magnitude as a function of wavelength.

    This interpolation is a convenient approximate stellar SED, not a
    spectrophotometric model. In particular, :meth:`get_magnitudes` returns a
    transmission-weighted average of interpolated magnitudes; it does not
    perform a flux-space synthetic-photometry integration.
    """
    def __init__(self, row):
        """Initialize a star from one catalog row.

        Parameters
        ----------
        row : pandas.Series
            Catalog record. Must provide ``ra``, ``dec``, ``source_id`` and
            all photometric columns accessed by :meth:`set_magnitudes`:
            ``phot_bp_mean_mag_AB``, ``phot_g_mean_mag_AB``,
            ``phot_rp_mean_mag_AB``, ``phot_w1_mean_mag_AB``,
            ``phot_w2_mean_mag_AB``, ``phot_w3_mean_mag_AB``,
            ``phot_w4_mean_mag_AB``, ``phot_j_mean_mag_AB``,
            ``phot_h_mean_mag_AB``, and ``phot_ks_mean_mag_AB``. An optional
            ``detector_id`` value is stored as ``self.WFI``.

        Notes
        -----
        Initialization stores the photometry and creates its wavelength
        interpolator. Missing required keys raise ``KeyError``.
        """
        self.ra = row["ra"]
        self.dec = row["dec"]
        self.source_id = row["source_id"]
        if 'detector_id' in row.keys():
            self.WFI = row["detector_id"]
        self.set_magnitudes(row)
        self.make_interpolator()

    def set_magnitudes(self, row):
        """Load the catalog's AB magnitudes and their representative wavelengths.

        The wavelengths are fixed reference values in Angstroms for the Gaia,
        2MASS, and WISE bands. This method populates the per-band attributes
        (for example, ``bp`` and ``bp_wave``), plus the ``mag`` and ``wave``
        arrays used to construct the interpolated SED.

        Parameters
        ----------
        row : pandas.Series
            Catalog record containing the photometry columns listed in
            :class:`Star`'s constructor documentation.
        """
        wavelengths = {
        'phot_bp_mean_mag_AB': 531.987,         # Gaia BP (blue photometer) ~505 nm
        'phot_g_mean_mag_AB': 671.955,          # Gaia G ~632 nm
        'phot_rp_mean_mag_AB': 793.910,         # Gaia RP (red photometer) ~797 nm
        'phot_w1_mean_mag_AB': 3352.60,        # WISE W1 ~3.4 μm
        'phot_w2_mean_mag_AB': 4602.80,        # WISE W2 ~4.6 μm
        'phot_w3_mean_mag_AB': 11560.8,       # WISE W3 ~12.0 μm
        'phot_w4_mean_mag_AB': 22088.3,       # WISE W4 ~22.0 μm
        'phot_j_mean_mag_AB': 1235.00,         # 2MASS J ~1.235 μm
        'phot_h_mean_mag_AB': 1662.00,         # 2MASS H ~1.662 μm
        'phot_ks_mean_mag_AB': 2159.00,        # 2MASS Ks ~2.159 μm
        }
        self.bp = row["phot_bp_mean_mag_AB"]
        self.bp_wave = wavelengths['phot_bp_mean_mag_AB']
        self.g = row["phot_g_mean_mag_AB"]
        self.g_wave = wavelengths['phot_g_mean_mag_AB']        
        self.rp = row["phot_rp_mean_mag_AB"]
        self.rp_wave = wavelengths['phot_rp_mean_mag_AB']
        self.w1 = row["phot_w1_mean_mag_AB"]
        self.w1_wave = wavelengths['phot_w1_mean_mag_AB']
        self.w2 = row["phot_w2_mean_mag_AB"]
        self.w2_wave = wavelengths['phot_w2_mean_mag_AB']
        self.w3 = row["phot_w3_mean_mag_AB"]
        self.w3_wave = wavelengths['phot_w3_mean_mag_AB']
        self.w4 = row["phot_w4_mean_mag_AB"]
        self.w4_wave = wavelengths['phot_w4_mean_mag_AB']
        self.j = row["phot_j_mean_mag_AB"]
        self.j_wave = wavelengths['phot_j_mean_mag_AB']
        self.h = row["phot_h_mean_mag_AB"]
        self.h_wave = wavelengths['phot_h_mean_mag_AB']
        self.ks = row["phot_ks_mean_mag_AB"]
        self.ks_wave = wavelengths['phot_ks_mean_mag_AB']

        self.mag = np.array([self.bp, self.g, self.rp, self.w1, self.w2, self.w3, self.w4, self.j, self.h, self.ks])
        self.wave = np.array([self.bp_wave, self.g_wave, self.rp_wave, self.w1_wave, self.w2_wave, self.w3_wave, self.w4_wave, self.j_wave, self.h_wave, self.ks_wave])
        self.wave_units = u.nm

        self.make_interpolator()

    def make_interpolator(self):
        """Create a linear AB-magnitude interpolator over the stored samples.

        The interpolator allows extrapolation outside the sampled wavelength
        range. It is stored as ``self._interpolator`` and used by :meth:`sed`.
        """
        self._interpolator = interp1d(self.wave, self.mag, kind='linear', bounds_error=False, fill_value="extrapolate")

    def sed(self, wave):
        """Evaluate the interpolated SED in AB magnitudes.
        
        Parameters
        ----------
            wave : astropy.units.Quantity or array-like
                Wavelength(s) at which to evaluate the SED. Any length unit
                supported by Astropy is accepted. Unitless inputs are assumed
                to be in ``self.wave_units``.
        Returns
        -------
        numpy.ndarray or numpy scalar
            Interpolated AB magnitudes. The returned values are plain numeric
            magnitudes, not an Astropy ``Quantity``. Values outside the
            photometric wavelength range are linearly extrapolated.
        """
        if not hasattr(self, '_interpolator'):
            self.make_interpolator()
        
        # Check units using astropy
        if not isinstance(wave, u.Quantity):
            wave = wave * self.wave_units
            warnings.warn(f"Input wavelength has no units; assuming {self.wave_units}.")
        
        wave = wave.to(self.wave_units).value
        return self._interpolator(wave)

    def get_magnitudes(self, bandpass, instrument='WFI', return_weights=False):
        """Estimate an AB magnitude for a Roman/WFI filter.
        
        Parameters
        ----------
            bandpass : str
                Roman/WFI filter name, such as ``"F158"``.
            instrument : str, optional
                Name of the instrument containing the specified filterband. Default is 'WFI'.
        
        Returns
        -------
        float
            Transmission-weighted mean of the interpolated AB magnitudes at
            the filter's sampled wavelengths.

        Notes
        -----
        This is an approximate average in magnitude space, calculated as
        ``sum(magnitude * transmission) / sum(transmission)``. It is not the
        standard flux-space synthetic AB magnitude. The filter curve is
        retrieved through :meth:`rosalia.telescopes.Roman.get_filter`.
        """

        band = rs.telescopes.Roman.get_filter(instrument=instrument, filter_name=bandpass)
        wave = band['wavelength_bins']
        transmission = band['transmission_bins']

        # Integrate the SED over the filter transmission to get the magnitude
        mags = self.sed(wave)
        mag_weights = mags * transmission/np.sum(transmission)
        mag_weighted = np.sum(mag_weights)

        output = (mag_weighted, (wave.to(u.nm).value, mag_weights)) if return_weights else mag_weighted
        return output

    def get_magnitude_weights(self, bandpass, instrument='WFI'):
        """Return each wavelength sample's contribution to the approximate magnitude.
        
        Parameters
        ----------
            bandpass : str
                Roman/WFI filter name, such as ``"F158"``.
            instrument : str, optional
                Name of the instrument containing the specified filterband. Default is 'WFI'.
        
        Returns
        ----------
            numpy.ndarray
                Values ``magnitude(wavelength) * transmission / sum(transmission)``
                at the filter's sampled wavelengths. Their sum equals the value
                returned by :meth:`get_magnitudes` (up to floating-point error).

        Notes
        -----
        Despite the method name, these are magnitude contributions, not
        dimensionless SED weights. They use the same approximate magnitude-
        space averaging as :meth:`get_magnitudes`.
        """
        
        band = rs.telescopes.Roman.get_filter(instrument=instrument, filter_name=bandpass)
        wave = band['wavelength_bins']
        transmission = band['transmission_bins']

        # Get the SED weights
        mags = self.sed(wave)
        mag_weights = mags * transmission/np.sum(transmission)

        return (wave.to(u.nm).value, mag_weights)

    def get_linear_weights(self, bandpass, instrument='WFI', zp=23.9, scale=0.11):
        """Return filter wavelengths and spectral flux contributions in nJy/pixel.

        Parameters
        ----------
        bandpass : str
            Roman/WFI filter name, such as ``"F158"``.
        instrument : str, optional
            Instrument containing the filter. Defaults to ``"WFI"``.
        zp : float, optional
            AB zero point in the nJy convention. Defaults to 23.9, for which
            ``f_nJy = 10**((zp - magnitude) / 2.5)``.
        scale : float, optional
            Pixel scale in arcsec/pixel. The returned flux contributions are
            scaled by ``scale**2`` to give nJy/pixel from nJy/arcsec**2.

        Returns
        -------
        tuple[numpy.ndarray, numpy.ndarray]
            ``(wavelength_nm, flux_contributions_nJy_per_pixel)``. The flux
            contributions sum to the flux corresponding to
            :meth:`get_magnitudes` multiplied by the pixel area.

        Notes
        -----
        The target total is kept consistent with the class's approximate
        magnitude-space filter average. The distribution across wavelength
        follows the interpolated SED, filter transmission, and wavelength-bin
        widths; it is intended for weighting monochromatic PSFs.
        """
        
        band = rs.telescopes.Roman.get_filter(instrument=instrument, filter_name=bandpass)
        wave = band['wavelength_bins']
        transmission = band['transmission_bins']

        # Preserve the same approximate band magnitude as get_magnitudes,
        # while distributing its flux over the band using the sampled SED.
        mags = self.sed(wave)
        wave_nm = wave.to(u.nm).value
        throughput = np.asarray(transmission, dtype=float)
        bin_width = np.gradient(wave_nm)
        response = np.clip(throughput, 0, None) * np.abs(bin_width)
        spectral_shape = 10**(-0.4 * (mags - np.nanmin(mags)))
        relative_flux = spectral_shape * response

        band_magnitude = self.get_magnitudes(bandpass, instrument=instrument)
        target_flux = 10**(0.4 * (zp - band_magnitude)) * scale**2
        if not np.isfinite(np.sum(relative_flux)) or np.sum(relative_flux) <= 0:
            raise ValueError(f"Filter {bandpass!r} has no usable transmission samples.")
        flux_weights = target_flux * relative_flux / np.sum(relative_flux)

        return (wave_nm, flux_weights)
    
    def get_galsim_sed(self):
        """Build and store a GalSim SED from the catalog photometry.

        Converts the stored AB magnitudes to flux densities and constructs a
        linearly interpolated GalSim spectrum. The resulting objects are
        stored on ``self`` as ``fnu``, ``spectrum_table``, ``galsim_sed``, and
        ``profile``; the latter is a delta-function point source multiplied
        by the SED.

        Returns
        -------
        None
            The GalSim objects are assigned to instance attributes.
        """
        import galsim

        # 0. Convert to flux density for galsim
        self.fnu = 3631e-23 * 10**(-0.4 * self.mag)

        # 1. Build a continuous spectrum via LookupTable interpolation
        self.spectrum_table = galsim.LookupTable(
            x=self.wave*1e9,  # Convert from meters to nanometers
            f=self.fnu, 
            interpolant='linear' # or 'cubic' if you have plenty of points and want smooth curves
        )

        # 2. Turn the table into a GalSim SED object
        self.galsim_sed = galsim.SED(
            spec=self.spectrum_table, 
            wave_type='nm', 
            flux_type='fnu',
        )

        # 3. Define a spatial star model (a point source) and apply the SED
        # This creates a ChromaticObject representing the star
        self.profile = galsim.DeltaFunction() * self.galsim_sed

    def plot_spectrum(self):
        """Plot the interpolated AB-magnitude SED and its catalog samples.

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
            Figure and axes containing the spectrum. The magnitude axis is
            inverted, as is conventional for plotting magnitudes.
        """
        import matplotlib.pyplot as plt

        dense_waves = np.linspace(self.wave.min(), self.wave.max(), 1000)*u.nm
        fig, ax = plt.subplots()
        ax.plot(dense_waves, self.sed(dense_waves),  label="Interpolator", color='royalblue', lw=2)
        ax.scatter(self.wave, self.mag, label='Photometry', color='darkorange', edgecolors='black', s=80, zorder=3)
        ax.invert_yaxis()
        ax.set_xlabel("Wavelength ($\mathrm{nm}$)", fontsize=12)
        ax.set_ylabel("AB Magnitude", fontsize=12)
        ax.legend()
        return fig, ax


#%%

############################

def recover_detector_frames(canvas, optimal_wcs, detector_wcs, detector_shapes=None):
    """
    Reproject a mosaic canvas into the pixel grid of each detector.

    The returned list follows the same ordering as ``detector_wcs`` and can be
    used like an exposure's ``DATA`` list.
    """
    from reproject import reproject_interp

    if np.ndim(canvas) != 2:
        raise ValueError("canvas must be a two-dimensional image")

    detector_wcs = list(detector_wcs)
    if detector_shapes is not None:
        detector_shapes = list(detector_shapes)
        if len(detector_shapes) != len(detector_wcs):
            raise ValueError("detector_shapes and detector_wcs must have the same length")

    detector_frames = []
    for i, wcs in enumerate(detector_wcs):
        shape = detector_shapes[i] if detector_shapes is not None else wcs.array_shape
        if shape is None:
            raise ValueError(
                "Detector image shape is required when its WCS has no array_shape"
            )

        frame, _ = reproject_interp(
            (canvas, optimal_wcs),
            output_projection=wcs,
            shape_out=tuple(shape),
        )
        detector_frames.append(frame)

    return detector_frames


ROMAN_N_DETECTORS = 18
ROMAN_DETECTOR_SIZE = 4088


def _roman_apertures():
    import pysiaf
    siaf = pysiaf.Siaf("Roman")
    detectors = [siaf["WFI%02d_FULL" % (i + 1)] for i in range(ROMAN_N_DETECTORS)]
    return siaf["WFI_CEN"], detectors


def _detector_to_canvas(center, detector, x, y, pixscale, origin):
    # x, y: 0-based detector pixels. Returns 0-based canvas pixels.
    v2, v3 = detector.sci_to_tel(np.asarray(x) + 1.0, np.asarray(y) + 1.0)
    xidl, yidl = center.tel_to_idl(v2, v3)
    return (xidl - origin[0]) / pixscale - 0.5, (yidl - origin[1]) / pixscale - 0.5


def _canvas_to_detector(center, detector, cx, cy, pixscale, origin):
    xidl = origin[0] + (np.asarray(cx) + 0.5) * pixscale
    yidl = origin[1] + (np.asarray(cy) + 0.5) * pixscale
    v2, v3 = center.idl_to_tel(xidl, yidl)
    x, y = detector.tel_to_sci(v2, v3)
    return x - 1.0, y - 1.0


def get_roman_canvas_geometry(pixscale=0.11):
    """
    Geometry of the full Roman/WFI canvas, computed only from pysiaf.
    The canvas lives in the ideal frame of WFI_CEN, with pixel size pixscale (arcsec).
    Returns the (x, y) ideal-frame origin of the canvas edge and its (ny, nx) shape.
    """
    center, detectors = _roman_apertures()
    edge = np.array([0.5, ROMAN_DETECTOR_SIZE + 0.5])
    xs, ys = [], []
    for detector in detectors:
        corner_x, corner_y = np.meshgrid(edge, edge)
        v2, v3 = detector.sci_to_tel(corner_x.ravel(), corner_y.ravel())
        xidl, yidl = center.tel_to_idl(v2, v3)
        xs.append(xidl)
        ys.append(yidl)
    xs, ys = np.concatenate(xs), np.concatenate(ys)
    origin = (xs.min(), ys.min())
    shape = (int(np.ceil((ys.max() - origin[1]) / pixscale)),
             int(np.ceil((xs.max() - origin[0]) / pixscale)))
    return origin, shape


def _sample_image(image, x, y):
    # Bilinear sampling at 0-based pixel positions. NaN outside the image.
    from scipy.ndimage import map_coordinates
    ny, nx = image.shape
    inside = (x >= -0.5) & (x <= nx - 0.5) & (y >= -0.5) & (y <= ny - 0.5)
    values = map_coordinates(image, [np.clip(y, 0, ny - 1), np.clip(x, 0, nx - 1)],
                             order=1, mode="nearest")
    return np.where(inside, values, np.nan)


def _fit_tan_wcs(pix_x, pix_y, ra, dec, sip_degree=None, array_shape=None):
    from astropy.wcs.utils import fit_wcs_from_points
    world = SkyCoord(ra, dec, unit="deg", frame="icrs")
    wcs = fit_wcs_from_points((pix_x, pix_y), world, proj_point="center",
                              projection="TAN", sip_degree=sip_degree)
    if array_shape is not None:
        wcs.array_shape = tuple(array_shape)
    return wcs


def get_roman_canvas_wcs(wcs_list, pixscale=0.11):
    """
    WCS and shape of the full Roman/WFI canvas from the 18 detector WCS, using only pysiaf
    for the layout (no reproject). Returns (canvas_wcs, shape) with shape = (ny, nx).
    """
    if len(wcs_list) != ROMAN_N_DETECTORS:
        raise ValueError("Expected %d detector WCS" % ROMAN_N_DETECTORS)
    center, detectors = _roman_apertures()
    origin, shape = get_roman_canvas_geometry(pixscale)
    grid = np.linspace(0, ROMAN_DETECTOR_SIZE - 1, 17)
    gx, gy = np.meshgrid(grid, grid)
    gx, gy = gx.ravel(), gy.ravel()
    cxs, cys, ras, decs = [], [], [], []
    for wcs, detector in zip(wcs_list, detectors):
        pcx, pcy = _detector_to_canvas(center, detector, gx, gy, pixscale, origin)
        ra, dec = wcs.all_pix2world(gx, gy, 0)
        cxs.append(pcx); cys.append(pcy); ras.append(ra); decs.append(dec)
    canvas_wcs = _fit_tan_wcs(np.concatenate(cxs), np.concatenate(cys),
                              np.concatenate(ras), np.concatenate(decs), array_shape=shape)
    return canvas_wcs, shape


def make_roman_canvas_from_detectors(data_list, wcs_list, pixscale=0.11):
    """
    Function 1. Combine the 18 Roman/WFI detector arrays (4088x4088, ordered WFI01 to WFI18)
    and their astropy WCS into a mosaic of the full field of view, using only pysiaf for the
    detector layout (no reproject). Gaps between detectors are NaN.

    Returns (canvas, canvas_wcs).
    """
    if len(data_list) != ROMAN_N_DETECTORS or len(wcs_list) != ROMAN_N_DETECTORS:
        raise ValueError("Expected %d detector arrays and WCS" % ROMAN_N_DETECTORS)

    center, detectors = _roman_apertures()
    origin, shape = get_roman_canvas_geometry(pixscale)
    canvas = np.full(shape, np.nan, dtype=np.float32)

    for data, wcs, detector in zip(data_list, wcs_list, detectors):
        data = np.asarray(data, dtype=np.float32)
        if data.shape != (ROMAN_DETECTOR_SIZE, ROMAN_DETECTOR_SIZE):
            raise ValueError("Each detector array must be %dx%d" % (ROMAN_DETECTOR_SIZE, ROMAN_DETECTOR_SIZE))

        # Canvas pixels covered by this detector: its footprint bounding box.
        cx_edge, cy_edge = _detector_to_canvas(center, detector, [0, ROMAN_DETECTOR_SIZE - 1] * 2,
                                               [0, 0, ROMAN_DETECTOR_SIZE - 1, ROMAN_DETECTOR_SIZE - 1],
                                               pixscale, origin)
        x0 = max(int(np.floor(cx_edge.min())) - 1, 0)
        x1 = min(int(np.ceil(cx_edge.max())) + 2, shape[1])
        y0 = max(int(np.floor(cy_edge.min())) - 1, 0)
        y1 = min(int(np.ceil(cy_edge.max())) + 2, shape[0])

        cx, cy = np.meshgrid(np.arange(x0, x1), np.arange(y0, y1))
        dx, dy = _canvas_to_detector(center, detector, cx, cy, pixscale, origin)
        values = _sample_image(data, dx, dy)
        region = canvas[y0:y1, x0:x1]
        fill = np.isnan(region) & ~np.isnan(values)
        region[fill] = values[fill]

    canvas_wcs, _ = get_roman_canvas_wcs(wcs_list, pixscale)
    return canvas, canvas_wcs


def recover_roman_detectors_from_canvas(canvas, canvas_wcs, pixscale=0.11, sip_degree=3, detector_wcs=None):
    """
    Function 2. Inverse of make_roman_canvas_from_detectors. Extract the 18 Roman/WFI detector
    arrays (4088x4088, WFI01 to WFI18) from a canvas built with the same pixscale, using pysiaf
    for the detector layout and the canvas WCS for the sky positions.

    Returns (data_list, wcs_list). Each detector WCS is a TAN fit with SIP distortion terms
    (sip_degree) to the canvas WCS evaluated at the detector pixels.

    If detector_wcs (18 astropy WCS) is given, the canvas can be on any grid (e.g. a swarp
    mosaic): each detector pixel is mapped to the sky with its WCS and then to the canvas with
    canvas_wcs. Those WCS are returned unchanged.
    """
    center, detectors = _roman_apertures()
    canvas = np.asarray(canvas, dtype=np.float32)
    if detector_wcs is not None:
        pix = np.arange(ROMAN_DETECTOR_SIZE)
        px, py = np.meshgrid(pix, pix)
        data_list = []
        for wcs in detector_wcs:
            ra, dec = wcs.all_pix2world(px, py, 0)
            cx, cy = canvas_wcs.all_world2pix(ra, dec, 0)
            data_list.append(_sample_image(canvas, cx, cy).astype(np.float32))
        return data_list, list(detector_wcs)

    origin, shape = get_roman_canvas_geometry(pixscale)
    if canvas.shape != shape:
        raise ValueError("Canvas shape %s does not match the expected %s for pixscale=%s"
                         % (canvas.shape, shape, pixscale))

    pix = np.arange(ROMAN_DETECTOR_SIZE)
    px, py = np.meshgrid(pix, pix)
    grid = np.linspace(0, ROMAN_DETECTOR_SIZE - 1, 17)
    gx, gy = np.meshgrid(grid, grid)
    gx, gy = gx.ravel(), gy.ravel()

    data_list, wcs_list = [], []
    for detector in detectors:
        cx, cy = _detector_to_canvas(center, detector, px, py, pixscale, origin)
        data_list.append(_sample_image(canvas, cx, cy).astype(np.float32))

        pcx, pcy = _detector_to_canvas(center, detector, gx, gy, pixscale, origin)
        ra, dec = canvas_wcs.all_pix2world(pcx, pcy, 0)
        wcs_list.append(_fit_tan_wcs(gx, gy, ra, dec, sip_degree=sip_degree,
                                     array_shape=(ROMAN_DETECTOR_SIZE, ROMAN_DETECTOR_SIZE)))

    return data_list, wcs_list


def scale_and_subtract_stars(input_name, ext, exposure_identity, g_mag_max=False, clean=False, verbose=False):
    ## TODO: Compute the scale factor for the object at (x,y)=(53,69) for
    ## the PSF (psf.fits). Compute it in the ring 20-30 pixels.
    # $ astscript-psf-scale-factor image.fits --mode=img \
    #      --center=53,69 --normradii=20,30 --psf=psf.fits

    ## Iterate over a catalog with RA,Dec positions of stars that are in
    ## the input image to compute their scale factors.
    #$ asttable catalog.fits | while read -r ra dec mag; do \
    #    astscript-psf-scale-factor image.fits \
    #        --mode=wcs \
    #        --psf=psf.fits \
    #        --center=$ra,$dec --quiet \
    # --normradii=20,30 > scale-"$ra"-"$dec".txt; done

    # Exposure identity properties
    lambda_ref = exposure_identity["FILTER_IDENTITY"]["filter_lambda_ref"]
    # Adding cache and style files 
    try: 
        psf_archive = os.environ["ROSALIACACHE"] + "/CORE/PSF_ARCHIVE/"
    except NameError:
        print("ROSALIACACHE environment variable not set. Please check your .bashrc or .zshrc file and define it. Example: export ROSALIACACHE=/path/to/rosaliacache/")
        print("rosalia.psf.scale_and_subtract_stars will not work without ROSALIACACHE/CORE/PSF_ARCHIVE defined.")

    psf_name = psf_archive + exposure_identity["FILTER"].lower() + "00_tinytim.fits"
    MJD = exposure_identity["EXPSTART"]

    # List of input images: Checking that it is a list.
    if isinstance(input_name, (list, pd.core.series.Series, np.ndarray)):
        input_list = input_name
        if verbose:
            print("Executing for images:")
            print(input_list)

    else:
        if verbose:
            print("Executing for single image " + input_name)
        input_list = [input_name]

    # Now we check if the extension is just one, or a list.
    if isinstance(ext, (list, pd.core.series.Series, np.ndarray)):
        ext_list = ext
    else:
        ext_list = [ext]

    ###########################
    # Fundamental estimations #
    ###########################
    data, wcs = rs.utils.get_data_and_wcs(input_name=input_name, ext=ext)

    star_search_radius = rs.utils.find_max_angular_size_of_image(data, wcs)

    # Now we execute the scale_and_subtract_stars_single_image program
    # File by file, and extension by extension.
    ###########################
    # Check that PSF is in the right format
    ##################

    rs.utils.check_number_of_extensions(psf_name, expected_next=2, verbose=verbose)

    subtracted_results = []
    output_name_list = []
    for input_name_i in tqdm(input_list):

        # We open a copy of the input_name_i that we are modifying.
        input_fits_i = fits.open(input_name_i)
        output_name_i =  input_name_i.replace(".fits", "_nostars.fits")
        for ext_i in ext_list:

            matched_stars_in_detector = find_stars_inside_detector(input_name=input_name_i,
                                                                   ext=ext_i,
                                                                   lambda_ref=lambda_ref,
                                                                   radius=star_search_radius,
                                                                   g_mag_max=g_mag_max,
                                                                   MJD=MJD, clean=clean,
                                                                   verbose=verbose)
            if verbose: print(matched_stars_in_detector)

            nstars_in_detector = len(matched_stars_in_detector["matched_catalog_image"]["RA_GAIA"])
            # A flat sky background correction is required to fit and subtract the PSF from stars
            # Otherwise it is impossible to scale correctly the PSF to the ring around the identified stars.
            # However, we do not want to do such correction with gradients just yet:
            # Ideally, we want to leave the gradient correction to the NDI, and then the residuals with Gnuastro.

            #rs.sky.correct_flat_sky(input_name=input_name_i, ext=ext_i, clean=clean)

            # If the telescope is HST ACS/WFC, then the Data Quality extension holds
            # a mask that marks which pixels are saturated.

            bool_is_saturated = np.zeros((nstars_in_detector), dtype="bool")
            if (exposure_identity["TELESCOP"] == "HST") & (exposure_identity["INSTRUME"] == "ACS"):
                # Saturated keyword in HST/ACS == 256. Lucas et al. 2018
                print(input_fits_i[ext+2].data)
                saturation_mask = input_fits_i[ext+2].data == 256

            # ------------------------------------------- #
            # Add here more telescopes to be supported.
            # ------------------------------------------- #

            ra_list  = np.zeros((nstars_in_detector))*np.nan
            dec_list = np.zeros((nstars_in_detector))*np.nan

            for i in range(len(ra_list)):
                x_ima = matched_stars_in_detector["matched_catalog_image"]["X_IMA"][i]
                y_ima = matched_stars_in_detector["matched_catalog_image"]["Y_IMA"][i]
                #x_ima, y_ima = rs.utils.radec_to_xy(ra=ra_list[i], dec=dec_list[i], fits_name=input_name_i, ext=ext_i)

                bool_is_saturated = saturation_mask[int(y_ima), int(x_ima)]
                if bool_is_saturated:
                    ra_list[i]  = matched_stars_in_detector["matched_catalog_image"]["RA_GAIA"][i]
                    dec_list[i] = matched_stars_in_detector["matched_catalog_image"]["DEC_GAIA"][i]
                    if verbose: print("Saturated Star! RA:" + str(ra_list[i]) + " DEC: " + str(dec_list[i]))
                else:
                    ra_list[i]  = matched_stars_in_detector["matched_catalog_image"]["RA_IMA"][i]
                    dec_list[i] = matched_stars_in_detector["matched_catalog_image"]["DEC_IMA"][i]



            #ra_list = matched_stars_in_detector["matched_catalog_image"]["RA_IMA"]
            #dec_list = matched_stars_in_detector["matched_catalog_image"]["DEC_IMA"]

            ds9_region_filename = input_name_i.replace(".fits", "_ext" + str(ext_i) + "_instars.reg")
            rs.utils.make_ds9_region(ra=ra_list, dec=dec_list,
                                     output=ds9_region_filename)

            # (input_name, ext, ra, dec, psf_name, clean=True):
            dictionary_with_results_from_star_subtraction = scale_and_subtract_stars_single_image(input_name = input_name_i, ext = ext_i,
                                                                                                  ra = ra_list, dec = dec_list,
                                                                                                  psf_name = psf_name, clean=False, verbose=verbose)
            subtracted_results.append(dictionary_with_results_from_star_subtraction)
            # Now we modify the extension with the residual image
            nostars_fits = fits.open(dictionary_with_results_from_star_subtraction["residuals"])
            input_fits_i[ext_i].data = nostars_fits[0].data
            nostars_fits.close()

        input_fits_i.verify("silentfix")
        input_fits_i.writeto(output_name_i, overwrite=True)
        input_fits_i.close()
        output_name_list.append(output_name_i)

        if clean:
            #for i in dictionary_with_results_from_star_subtraction["residuals"]:
            #    rs.utils.execute_cmd("rm " + i)

            rs.utils.execute_cmd("rm " + input_name_i.replace(".fits", "*detected.fits"))
            rs.utils.execute_cmd("rm " + input_name_i.replace(".fits", "*instars.reg"))
            rs.utils.execute_cmd("rm " + input_name_i.replace(".fits", "*_matched_stars_inside_detector.fits"))
            rs.utils.execute_cmd("rm " + input_name_i.replace(".fits", "*_gaia_stars_inside_detector.csv"))
            rs.utils.execute_cmd("rm " + input_name_i.replace(".fits", "*_gaia.dat"))
            rs.utils.execute_cmd("rm " + input_name_i.replace(".fits", "*_detected_segmented_cat.fits"))

    subtracted_results_dataframe = pd.DataFrame(subtracted_results)

    return(output_name_list)

############################

def astscript_psf_scale_factor(input_name, ra, dec, psf_name, clean=False, verbose=False):
    # Python wrapper of astscript-psf-scale-factor
    #
    # Based on https://www.gnu.org/software/gnuastro/manual/html_node/Invoking-astscript_002dpsf_002dscale_002dfactor.html
    #
    # Shell based example
    # cmd = "astscript-psf-scale-factor mastDownload/HST/j9en0pnbq/j9en0pnbq_flc_ext1.fits --mode=wcs --center=329.03478737,1.4031257 --normradii=10,20 --psf=PSF/psf_hst_test.fits  | tail -1"
    # scaling_factor = os.popen(cmd).read()
    # print(scaling_factor)



    #print("WARNING:")
    if verbose:
        print("Estimating norm_radiis:")
        print("RA:" + str(ra))
        print("DEC:" + str(dec))

    max_size_of_a_star = 600 # pixels
    # Lets gather the max detection radiis for all the stars

    measure_maxradii_result = measure_maxradii(input_name=input_name, ra=ra, dec=dec, rmax=max_size_of_a_star, clean=clean, verbose=verbose)

    rmax = measure_maxradii_result["rlim"]

    norm_radii_min = int(rmax/3.)
    norm_radii_max = int(2*rmax/3.)

    if verbose:
        print("norm_radii_min:" + str(norm_radii_min))
        print("norm_radii_max:" + str(norm_radii_max))


    # We construct the cmd to send to shell.
    cmd = "astscript-psf-scale-factor " + input_name +\
           " --mode=wcs --quiet --center=" + str(ra) + "," + str(dec) +\
           " --normradii=" + str(norm_radii_min) +"," + str(norm_radii_max) +\
           " --psf=" + psf_name + "  | tail -1"

    # We execute the command and get the output, but only the last line ( blahblah .... | tail -1)

    if verbose:
        print(cmd)

    scaling_factor = float(rs.utils.execute_cmd(cmd))
    if verbose: print("Scaling factor: " + str(scaling_factor))
    return(scaling_factor)

############################

def astscript_radial_profile(input_name, ra, dec, rmax, clean=False, verbose=False):
    # This is a wrapper for Gnuastro astcript_radial_profile
    output_name = input_name.replace(".fits", "_profile_ra" + str(ra) + "_dec" + str(dec) + ".fits")

    # We construct the cmd to send to shell.
    cmd = "astscript-radial-profile " + input_name +\
           " --mode=wcs --measure=median --center=" + str(ra) + "," + str(dec) +\
           " --rmax=" + str(rmax) +\
           " -o " + output_name

    # We execute the command and get the output, but only the last line ( blahblah .... | tail -1)
    if not clean:
        if verbose:
            print(cmd)

    rs.utils.execute_cmd(cmd)

    if not os.path.exists(output_name):
        if verbose:
            print("Caution: Profile not generated: " + input_name + " ra: " + str(ra) + " dec " + str(dec))
            print("Most likely, the target is not in the field of view.")
        return(False)

    if verbose:
        print("Profile for " + input_name + " ra: " + str(ra) + " dec " + str(dec))
        print("saved in " + output_name)



    return({"profile_name": output_name})

############################

def measure_maxradii(input_name, ra, dec, rmax, saturation_level=None, clean=False, verbose=False, clean_profile=False):
    star_profile_dict = astscript_radial_profile(input_name=input_name, ra=ra, dec=dec, rmax=rmax, verbose=verbose)

    # If the target is not in the FOV, then astscript_radial_profile will return None.
    # Propagate so it can be detected by other programs and the next line does not raise an error.

    if star_profile_dict is False:
        print("Caution: Star profile not found! Returning None")
        return(False)

    star_profile = fits.open(star_profile_dict["profile_name"])

    r = star_profile[1].data["RADIUS"]
    sky = star_profile[1].data["MEDIAN"][int(rmax/2):int(rmax-1)]
    sb = star_profile[1].data["MEDIAN"] - bn.nanmedian(sky) #+ 0.5# - np.nanmin(db[1].data["MEAN"])

    # We make a first approximation to the rlim
    noise_level = bn.nanstd(sky)

    if verbose:
        plt.axhline(noise_level, label="Noise level 1")

    rlim = star_profile[1].data["RADIUS"][np.where(((sb - 3*noise_level) < 0) & (r > 2))[0][0]] # We put the (r>2) to avoid the saturated core
    if verbose:
        print("Rlim 1: " + str(rlim))
        plt.axvline(rlim, color="purple", label="Rlim 1")
    # Refine the measurement using the s1_down limit (sigma-clipping)
    r = star_profile[1].data["RADIUS"]
    sky = star_profile[1].data["MEDIAN"][int(rlim):]
    sb_nosky_cor = star_profile[1].data["MEDIAN"]
    sb = sb_nosky_cor - bn.nanmedian(sky) #+ 0.5# - np.nanmin(db[1].data["MEAN"])
    noise_level = bn.nanstd(sky)

    #rlim = star_profile[1].data["RADIUS"][np.where((s1_down < 0) & (r > 2))[0][0]] # We put the (r>2) to avoid the saturated core
    rlim = star_profile[1].data["RADIUS"][np.where(((sb - 3*noise_level) < 0) & (r > 2))[0][0]] # We put the (r>2) to avoid the saturated core

    if saturation_level is not None:
        where_is_saturated = np.where(sb_nosky_cor > saturation_level)[0]
        if verbose:
            print("Where is saturated?")
            print(where_is_saturated)
        if len(where_is_saturated) > 1:
            rsat = star_profile[1].data["RADIUS"][where_is_saturated[-1]]
        elif len(where_is_saturated) == 1:
            rsat = star_profile[1].data["RADIUS"][where_is_saturated[0]]
        else:
            rsat = np.nan
    else:
        rsat = np.nan


    if verbose:
        print("Rlim 2: " + str(rlim))
        plt.scatter(r, sb, label="Star profile")
        plt.yscale("log")
        plt.xscale("log")
        if saturation_level is not None:
            plt.axvline(rsat, color="firebrick", label="Saturation level")

        plt.axhline(noise_level, color="orange", label="Noise level 2")
        plt.axvline(rlim, color="red", label="Rlim 2")
        plt.xlabel("R (pixels)")
        plt.ylabel("Intensity (flux)")
        plt.legend()
        plt.show(block=False)

    if clean_profile: rs.utils.execute_cmd("rm " + star_profile_dict["profile_name"])

    return({"rlim": rlim, "rsat": rsat})

############################

def astscript_psf_subtract(input_name, ext, ra, dec, scaling_factor, psf_name, clean=True):
    # Python wrapper of astscript-psf-subtract
    #
    # Based on https://www.gnu.org/software/gnuastro/manual/html_node/Invoking-astscript_002dpsf_002dsubtract.html
    #
    # Shell based example
    # execute_cmd("astscript-psf-subtract mastDownload/HST/j9en0pnaq/j9en0pnaq_flc_ext1.fits \
    #       --mode=wcs \
    #       --psf=PSF/psf_hst_test.fits \
    #       --scale=1.56063e+06 \
    #       --center=329.03478737,1.4031257 \
    #       --output=mastDownload/HST/j9en0pnaq/j9en0pnaq_flc_ext1_sub.fits")

    # First we generate an image with the main WCS, otherwise Noisechisel cannot read the header
    #image_wcs_cleaned_extension_name = rs.mast.make_wcs_clean_extension(input_name, ext)
    output_name = input_name.replace(".fits", "_" + str(ra) + "_" + str(dec) + "_sub.fits")

    # We construct the cmd to send to shell.
    cmd = "astscript-psf-subtract " + input_name +\
           " --mode=wcs --psf="+ psf_name +\
           " --scale="+ str(scaling_factor) + " --center=" + str(ra)+","+str(dec)+\
           " --output=" + output_name

    if not clean:
        print(cmd)
    output_from_shell = rs.utils.execute_cmd(cmd)
    return(output_name)

############################

def scale_and_subtract_stars_single_image(input_name, ext, ra, dec, psf_name, clean=True, extrapolate=False, verbose=False):

    # For now, we will set a fixed norm_radii_min, and norm_radii_max.
    # But this is a function of the saturation level of the star.
    # We need to measure this an scale accordingly.


    input_fits = fits.open(input_name)
    star_model = np.zeros(input_fits[ext].data.shape)
    # First we generate an image with the main WCS, otherwise Noisechisel cannot read the header
    image_wcs_cleaned_extension_name = rs.gnu.make_wcs_clean_extension(input_name, ext)
    # We correct the sky of the extension
    rs.sky.correct_flat_sky(input_name=image_wcs_cleaned_extension_name, ext=1, clean=False)

    # We go star by star executing the wrappers.
    for ra_i, dec_i in tqdm(zip(ra, dec), total=len(ra), position=0, leave=True):

        # First we find the scaling factor
        scaling_factor = astscript_psf_scale_factor(input_name=image_wcs_cleaned_extension_name,
                                                    ra=ra_i, dec=dec_i,
                                                    psf_name=psf_name, clean=clean, verbose=verbose)

        # Then we subtract the star
        subtracted_image_name = astscript_psf_subtract(input_name=image_wcs_cleaned_extension_name, ext=1,
                                                       ra=ra_i, dec=dec_i,
                                                       scaling_factor=scaling_factor, psf_name=psf_name, clean=clean)

        # We open the resulting subtracted image
        subtracted_fits = fits.open(subtracted_image_name)

        # And we add what we have subtracted to the star_model
        image_wcs_cleaned_extension = fits.open(image_wcs_cleaned_extension_name)
        star_model = star_model + image_wcs_cleaned_extension[1].data - subtracted_fits[1].data

        # Clean the temporary file
        rs.utils.execute_cmd("rm " + subtracted_image_name)

    # Finally, we save the star model in a new fits.
    output_starmodel_name = input_name.replace(".fits", "_ext" + str(ext) + "_stars.fits")
    output_residual_name = input_name.replace(".fits", "_ext" + str(ext) + "_nostars.fits")

    rs.utils.save_fits(array=star_model, name=output_starmodel_name, header=input_fits[ext].header, extname='STAR_MODEL')
    subtracted_image = input_fits[ext].data - star_model
    rs.utils.save_fits(array=subtracted_image, name=output_residual_name, header=input_fits[ext].header, extname='STAR_SUBTRACTED')

    return({"star_model":output_starmodel_name, "residuals":output_residual_name})

############################

def gaia_find_stars_in_and_out(input_name, ext, lambda_ref, ra=None, dec=None, radius=0.1, g_mag_max=False,
                               MJD=None, clean=True, verbose=False, gaia_query=False, match_inside_stars=False):
    # First find the stars using Gaia
    if verbose: print("> Finding Gaia/2MASS/WISE stars around image...")

    if gaia_query is False:
        gaia_query_dict = find_gaia_stars_around_image(lambda_ref=lambda_ref, input_name=input_name, ext=ext,
                                                   ra=ra, dec=dec, MJD=MJD, radius=radius, g_mag_max=g_mag_max, verbose=verbose)
        gaia_query = gaia_query_dict["gaia_query"]



    # Now identify which ones are actually inside the detector
    if verbose: print("> Identifying which stars are inside the FOV and which are outside...")
    input_fits = fits.open(input_name)
    infield_stars = identify_stars_in_out_field(data_shape=input_fits[ext].data.shape(),
                                                wcs=astropy_wcs.WCS(input_fits[ext].header, input_fits),
                                                catalog=gaia_query,
                                                verbose=verbose) # identify_stars_in_out_field(input_name=input_name, ext=ext, catalog=gaia_query, ra=ra, dec=dec, verbose=verbose)
    # Lets refine the position of the stars
    # Make a catalog of the objects inside the image
    bool_isIn_gaia_stars = infield_stars["bool_isIn"]

    if (ra is None) and (dec is None):
        if verbose: print("> Generating a catalog of objects in the image...")
        if match_inside_stars:
            object_catalog_in_detector = rs.gnu.astmkcatalog(input_name=input_name, ext=ext, clean=clean, verbose=verbose)
        else:
            object_catalog_in_detector = [False]

        if verbose: print("> Catalog done.")

        # Cross-match the coordinates, and find the real centroids of the stars.
        # We will need to define a magnitude limit, possibly. And a limit in distance.

        # Save a catalog only with the Gaia stars inside the detector
        in_query_filename = input_name.replace(".fits", "_ext" + str(ext) + "_gaia_stars_inside_detector.csv")
        out_query_filename = input_name.replace(".fits", "_ext" + str(ext) + "_gaia_stars_outside_detector.csv")

        gaia_query_in_detector = gaia_query[bool_isIn_gaia_stars]
        gaia_query_in_detector.to_csv(in_query_filename)

        gaia_query_out_detector = gaia_query[np.invert(bool_isIn_gaia_stars)]
        gaia_query_out_detector.to_csv(out_query_filename)
        return({"gaia_query_in_detector": gaia_query_in_detector, "in_query_filename": in_query_filename,
                "gaia_query_out_detector": gaia_query_out_detector, "out_query_filename": out_query_filename,
                "object_catalog_in_detector": object_catalog_in_detector})

    else:
        out_query_filename = str(ra) + "_" + str(dec) + "_gaia_stars_outside_detector.csv"
        gaia_query_out_detector = gaia_query[np.invert(bool_isIn_gaia_stars)]
        gaia_query_out_detector.to_csv(out_query_filename)
        return({"gaia_query_out_detector": gaia_query_out_detector, "out_query_filename": out_query_filename})

############################

def find_stars_inside_detector(input_name, ext, lambda_ref, radius, g_mag_max=False,  MJD=None, clean=True, verbose=False):
    """
    Queries Gaia archive and finds the stars inside a given detector, provided as a FITS file.

    """


    gaia_search_db = gaia_find_stars_in_and_out(input_name=input_name, ext=ext, lambda_ref=lambda_ref,
                                                ra=None, dec=None, radius=radius, g_mag_max=g_mag_max,
                                                MJD=MJD, clean=clean, verbose=verbose, gaia_query=False,
                                                match_inside_stars=True)
    if verbose: print(gaia_search_db)
    object_catalog_in_detector = gaia_search_db["object_catalog_in_detector"]

    matched_catalog_name =  input_name.replace(".fits", "_ext" + str(ext) + "_matched_stars_inside_detector.fits")

    # Reformat the catalog to a format that Gnuastro can read
    inside_star_catalog = pd.read_csv(gaia_search_db["in_query_filename"])
    inside_star_catalog_for_gnuastro = inside_star_catalog[['ra', 'dec', 'mag_lambda']].copy()
    inside_star_catalog_name = gaia_search_db["in_query_filename"].replace(".csv", "_gnuastro.csv")
    inside_star_catalog_for_gnuastro.to_csv(inside_star_catalog_name,  header=False, index=False)


    # Run astmatch and generate a new catalog of matched objects.
    cmd = "astmatch -K  " + inside_star_catalog_name + " -h2 " + object_catalog_in_detector["output_name"] + " --aperture=1/360 --ccol1=1,2 --ccol2=ra,dec --output " + matched_catalog_name + " --outcols=a1,a2,b1,b2,b3,b4,a3"
    if verbose: print(cmd)
    rs.utils.execute_cmd(cmd)

    # Write the correct table names.
    cmd = "astfits -h1  " + matched_catalog_name +\
           " --update=TTYPE1,RA_GAIA    --update=TUNIT1,deg "+\
           " --update=TTYPE2,DEC_GAIA   --update=TUNIT2,deg "+\
           " --update=TTYPE3,RA_IMA     --update=TUNIT3,deg "+\
           " --update=TTYPE4,DEC_IMA    --update=TUNIT4,deg "+\
           " --update=TTYPE5,X_IMA      --update=TUNIT5,pixel "+\
           " --update=TTYPE6,Y_IMA      --update=TUNIT6,pixel "+\
           " --update=TTYPE7,MAG_LAMBDA --update=TUNIT7,mag "

    if verbose: print(cmd)
    rs.utils.execute_cmd(cmd)

    matched = fits.open(matched_catalog_name)
    matched_catalog_image = matched[1].data


    return({"matched_catalog_image": matched_catalog_image,
            "matched_catalog_name": matched_catalog_name})
############################

# Function to append n empty rows to a DataFrame
def append_empty_rows(dataframe, n):
    for _ in range(n):
        dataframe.loc[len(dataframe)] = pd.Series(dtype='float64')

############################

def identify_stars_in_out_field(data_shape, wcs, catalog, ra=None, dec=None, verbose=False):
    # This program identifies the objects that are inside or outside an image in a fits file.
    # input_name: Name of the fits file
    # ext: Extension of the fits file where we are looking for the image.
    #      i.e., For ACS/WFC in HST, the science frames are stored in ext 1 and 4.
    #      i.e., For WFC/IR, the science frame is in ext 1.
    # catalog: A table containing at least two columns called "ra" and "dec".
    if (ra is not None) and (dec is not None):
        # Plot only stars brighter than N magnitudes
        plt.figure(figsize=(10,10))
        rs.plots.plot_stars_around(catalog)
        print('RA' + str(ra))
        print('DEC' + str(dec))
        plt.scatter(ra, dec, marker="+", color="navy", alpha=0.5, s=500)


        plt.gca().invert_xaxis()
        plt.show(block=False)
        return({"bool_isIn": np.zeros(len(catalog["ra"])).astype("bool"),
                "corners_pix": None,
                "corners_world": None,
                "catalog_inside": None})


    detector_corners = rs.detectors.get_detector_corners(wcs)
    corners_pix = detector_corners["corners_pix"]
    corners_world = detector_corners["corners_world"]

    ######### Check if something is infield or outfield #########
    if verbose: print(" Detector corners world: " + str(corners_world))
    if verbose: print(" Detector corners pixel: " + str(corners_pix))

    vertices_detector_dec_ra = np.zeros(corners_world.shape)
    vertices_detector_dec_ra[:,1] = corners_world[:,0]
    vertices_detector_dec_ra[:,0] = corners_world[:,1]

    from sphericalpolygon import Sphericalpolygon
    #print("DEMO WARNING: This can be done way faster by prefiltering those stars at high angular distances.")
    region_detector = Sphericalpolygon.from_array(vertices_detector_dec_ra)


    distance_detector_corners = rs.utils.angular_distance(ra1=np.array(corners_world[:,0]),
                                                          dec1=np.array(corners_world[:,1]),
                                                          ra2=np.array(corners_world[:,0]),
                                                          dec2=np.array(corners_world[:,1]))

    aprox_size_detector = bn.nanmax(distance_detector_corners)


    distance_detector_to_stars = bn.nanmax(rs.utils.angular_distance(ra1=catalog["ra"],
                                                                     dec1=catalog["dec"],
                                                                     ra2=np.array(corners_world[:,0]),
                                                                     dec2=np.array(corners_world[:,1])), axis=1)

    
    #if isinstance(catalog["ra"], (pd.Series, np.ndarray)):
    if len(catalog["ra"]) > 1:
        isIn = np.zeros(len(catalog["ra"])).astype("bool")
        far_away_stars = distance_detector_to_stars > 5*aprox_size_detector
        isIn[far_away_stars] = False
        close_stars_around_list_dec_ra = np.array([[j, i] for i, j in zip(catalog["ra"], catalog["dec"])])
        isIn[~far_away_stars] = region_detector.contains_points(close_stars_around_list_dec_ra[~far_away_stars])
        catalog_inside = catalog[isIn]
        if verbose: print(" Stars inside detector: " + str(catalog_inside))
        if verbose: print(" All stars catalog: " + str(catalog))
        return({"bool_isIn": isIn, "corners_pix": corners_pix, "corners_world": corners_world, "catalog_inside": catalog_inside})
    
    #elif isinstance(catalog["ra"], (np.floating, float)):
    else:
        isIn = region_detector.contains_points([[catalog["dec"].iloc[0], catalog["ra"].iloc[0]]])
        if isIn:
            return({"bool_isIn": True, 
                    "corners_pix": corners_pix, 
                    "corners_world": corners_world, 
                    "catalog_inside": catalog})
        else:
            return({"bool_isIn": False, 
                    "corners_pix": corners_pix, 
                    "corners_world": corners_world, 
                    "catalog_inside": catalog})

############################

def get_hybrid_catalog(ra, dec, radius, lambda_ref, MJD, observer, g_mag_max=False, verbose=False, query_filename="default_query.dat"):
    from scipy import interpolate


    # C is the central coordinates of the pointing.
    c = SkyCoord(ra, dec, frame='icrs', unit="deg")

    # First get the high resolution catalog.
    if verbose: print("Running query_gaia_2mass_wise query...")
    high_resolution_catalog = rs.gaia.query_gaia_2mass_wise(ra=ra, dec=dec,
                                                    radius=radius,
                                                    g_mag_max=g_mag_max,
                                                    verbose=verbose,
                                                    query_filename=query_filename)
    # Load the low resolution catalog.
    low_resolution_catalog = pd.read_csv(os.environ["ROSALIACACHE"] + "/CORE/gaia_2mass_wise_healpix_7_cat.csv")

    # Make a new low_resolution_catalog where the columns will fit.
    # Maybe stick this in the low_res
    # ----------------------------------------------- #

    # Remove the super stars in the high res region
    ipix_highres = sorted(set(high_resolution_catalog["healpix_lvl"]))
    if verbose: print("Healpix cells to remove from low resolution:" + str(ipix_highres))

    row_healpix_to_drop_list = []
    for healpix_covered_in_highres_area in ipix_highres:
        if verbose: print(healpix_covered_in_highres_area)
        row_healpix_to_drop = np.where(low_resolution_catalog["healpix_lvl"] == healpix_covered_in_highres_area)[0]
        if verbose: print(row_healpix_to_drop)
        if verbose: print(low_resolution_catalog["healpix_lvl"].iloc[row_healpix_to_drop])
        row_healpix_to_drop_list.append(row_healpix_to_drop[0])
        # Before dropping, check if the sum of fluxes is the same.
        #print("High res g flux: " + str(bn.nansum(high_resolution_catalog["g"].iloc[np.where(high_resolution_catalog["healpix_5"] == healpix_covered_in_highres_area)])))
        #print("Low res g flux: " + str(bn.nansum(low_resolution_catalog["g"].iloc[row_healpix_to_drop])))
        #print("High res j flux: " + str(bn.nansum(high_resolution_catalog["j"].iloc[np.where(high_resolution_catalog["healpix_5"] == healpix_covered_in_highres_area)])))
        #print("Low res j flux: " + str(bn.nansum(low_resolution_catalog["j"].iloc[row_healpix_to_drop])))
        #print("High res w1 flux: " + str(bn.nansum(high_resolution_catalog["w1"].iloc[np.where(high_resolution_catalog["healpix_5"] == healpix_covered_in_highres_area)])))
        #print("Low res w1 flux: " + str(bn.nansum(low_resolution_catalog["w1"].iloc[row_healpix_to_drop])))

    low_resolution_catalog = low_resolution_catalog.drop(row_healpix_to_drop_list)

    #high_resolution_catalog.append(low_resolution_catalog)
    high_resolution_catalog = pd.concat([high_resolution_catalog, low_resolution_catalog], axis=0)

    ######################################################################
    #
    # Now complement the query with the 230 naked eye stars from Hipparcos
    #
    ######################################################################
    gaia_naked_eye_stars = pd.read_csv(os.environ["ROSALIACACHE"] + "/CORE/gaia_naked_eye_stars.csv", sep=" ")

    # Find which ones are in the query radius.
    ra_stars = np.array(gaia_naked_eye_stars["ra"])*u.deg
    dec_stars = np.array(gaia_naked_eye_stars["dec"])*u.deg
    gaia_coords = SkyCoord(ra_stars, dec_stars, frame='icrs')
    sep = c.separation(gaia_coords)

    n_bright_stars = len(gaia_naked_eye_stars["ra"])
    n_highres_stars = len(high_resolution_catalog)

    print("\n[DEMO WARNING:] Naked eye stars do not have photometry besides g band")
    #append_empty_rows(high_resolution_catalog, n_bright_stars)
    if verbose: print("> Adding bright stars from Hipparcos survey...")

    source_id_bright_star = []
    ra_bright_star = []
    dec_bright_star = []
    bp_bright_star = []
    g_bright_star = []
    rp_bright_star = []
    j_bright_star = []
    h_bright_star = []
    ks_bright_star = []
    w1_bright_star = []
    w2_bright_star = []
    w3_bright_star = []
    w4_bright_star = []

    phot_bp_mean_mag_AB_bright_star = []
    phot_g_mean_mag_AB_bright_star = []
    phot_rp_mean_mag_AB_bright_star = []
    phot_j_mean_mag_AB_bright_star = []
    phot_h_mean_mag_AB_bright_star = []
    phot_ks_mean_mag_AB_bright_star = []
    phot_w1_mean_mag_AB_bright_star = []
    phot_w2_mean_mag_AB_bright_star = []
    phot_w3_mean_mag_AB_bright_star = []
    phot_w4_mean_mag_AB_bright_star = []

    for i in tqdm(range(len(gaia_naked_eye_stars)), disable=not verbose, position=0, leave=True):
        source_id_bright_star.append("Hipparcos bright star")
        ra_bright_star.append(gaia_naked_eye_stars["ra"][i])
        dec_bright_star.append(gaia_naked_eye_stars["dec"][i])
        bp_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_bp_ab"][i])))
        g_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_g_ab"][i])))
        rp_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_rp_ab"][i])))
        j_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_j_ab"][i])))
        h_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_h_ab"][i])))
        ks_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_ks_ab"][i])))
        w1_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_w1_ab"][i])))
        w2_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_w2_ab"][i])))
        w3_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_w3_ab"][i])))
        w4_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_w4_ab"][i])))

        phot_bp_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_bp_ab"][i])
        phot_g_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_g_ab"][i])
        phot_rp_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_rp_ab"][i])
        phot_j_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_j_ab"][i])
        phot_h_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_h_ab"][i])
        phot_ks_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_ks_ab"][i])
        phot_w1_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_w1_ab"][i])
        phot_w2_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_w2_ab"][i])
        phot_w3_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_w3_ab"][i])
        phot_w4_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_w4_ab"][i])

    bright_star_catalog = pd.DataFrame({"source_id": source_id_bright_star,
                           "ra": ra_bright_star,
                           "dec": dec_bright_star,
                           "bp": bp_bright_star, "phot_bp_mean_mag_AB": phot_bp_mean_mag_AB_bright_star,
                           "g": g_bright_star, "phot_g_mean_mag_AB": phot_g_mean_mag_AB_bright_star,
                           "rp": rp_bright_star, "phot_rp_mean_mag_AB": phot_rp_mean_mag_AB_bright_star,
                           "j": j_bright_star, "phot_j_mean_mag_AB": phot_j_mean_mag_AB_bright_star,
                           "h": h_bright_star, "phot_h_mean_mag_AB": phot_h_mean_mag_AB_bright_star,
                           "ks": ks_bright_star, "phot_ks_mean_mag_AB": phot_ks_mean_mag_AB_bright_star,
                           "w1": w1_bright_star, "phot_w1_mean_mag_AB": phot_w1_mean_mag_AB_bright_star,
                           "w2": w2_bright_star, "phot_w2_mean_mag_AB": phot_w2_mean_mag_AB_bright_star,
                           "w3": w3_bright_star, "phot_w3_mean_mag_AB": phot_w3_mean_mag_AB_bright_star,
                           "w4": w4_bright_star, "phot_w4_mean_mag_AB": phot_w4_mean_mag_AB_bright_star})

    high_resolution_catalog = pd.concat([high_resolution_catalog, bright_star_catalog], axis=0)

    # > Here we generate a new column with the custom user magnitude
    mag_lambda = np.zeros(len(high_resolution_catalog))
    filters_list_lambda = np.array([rs.telescopes.Passbands.Gaia.bp.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.Gaia.g.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.Gaia.rp.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.TwoMASS.J.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.TwoMASS.H.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.TwoMASS.Ks.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.WISE.W1.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.WISE.W2.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.WISE.W3.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.WISE.W4.lambda_ref.to("AA").value])

    magnitude_array = np.array([high_resolution_catalog["phot_bp_mean_mag_AB"],
                                high_resolution_catalog["phot_g_mean_mag_AB"],
                                high_resolution_catalog["phot_rp_mean_mag_AB"],
                                high_resolution_catalog["phot_j_mean_mag_AB"],
                                high_resolution_catalog["phot_h_mean_mag_AB"],
                                high_resolution_catalog["phot_ks_mean_mag_AB"],
                                high_resolution_catalog["phot_w1_mean_mag_AB"],
                                high_resolution_catalog["phot_w2_mean_mag_AB"],
                                high_resolution_catalog["phot_w3_mean_mag_AB"],
                                high_resolution_catalog["phot_w4_mean_mag_AB"]])

    # Transform the magnitudes to the observing filter.
    #print("Demo warning: Using Gaia g filter magnitudes instead of the correct filter")
    #print("Query Gaia/2MASS/WISE https://notebook.community/bosscha/alma-calibrator/notebooks/2mass/13_environment-with_analysis")

    for i in tqdm(range(len(high_resolution_catalog)), disable=not verbose, position=0, leave=True):
        good = np.where(np.isfinite(magnitude_array[:,i]))
        #print(magnitude_array[:,i])
        magnitude_interpolator = interpolate.interp1d(x=filters_list_lambda[good], y=magnitude_array[good,i], kind="linear", fill_value="extrapolate")
        mag_lambda[i] = magnitude_interpolator(lambda_ref.to("AA").value)[0]

    high_resolution_catalog["mag_lambda"] = mag_lambda


    # magnitude_interpolator = interpolate.interp1d(x=mag_sun_db["lambda_mum"]*10**4, y=mag_sun_db["AB_App"], kind="linear")
    # mag_ref = magnitude_interpolator(lambda_ref.to("AA"))


    # Here, calculate the positions of the planets in L2.
    # For the filter considered, add the magnitude of the planets.
    planets_catalog = rs.sso.get_SS0s_loc_magnitude(observer=observer, MJD=MJD, lambda_ref=lambda_ref, verbose=verbose) # Awesome it is working
    if verbose: print(planets_catalog)
    planets_catalog_to_concat = pd.DataFrame({"source_id": planets_catalog["source_id"],
                                              #"source_id_2": planets_catalog["source_id_2"],
                                              "ra": planets_catalog["RA"],
                                              "dec": planets_catalog["DEC"],
                                              "mag_lambda": planets_catalog["mag_custom"]})
    high_resolution_catalog = pd.concat([high_resolution_catalog, planets_catalog_to_concat], axis=0)
    #print("DEMO WARNING: ADD THE MOON AND REST OF SSOs to psf.get_hybrid_catalog")
    #print("https://arxiv.org/pdf/1808.01973.pdf")

    # Find the distance to the pointing.
    ra_stars = np.array(high_resolution_catalog["ra"])*u.deg
    dec_stars = np.array(high_resolution_catalog["dec"])*u.deg
    star_coords = SkyCoord(ra_stars, dec_stars, frame='icrs')
    sep = c.separation(star_coords)
    high_resolution_catalog["dist"] = sep.deg
    high_resolution_catalog = high_resolution_catalog.sort_values(by=["mag_lambda"], ascending=True)
    high_resolution_catalog["cat_id"] = np.linspace(1, len(high_resolution_catalog), len(high_resolution_catalog), dtype="int64")
    high_resolution_catalog.to_csv(query_filename)

    return(high_resolution_catalog)

############################

def psf_harvester_single_exposure(input_name, clean=True, verbose=False):
    # if True:
    import astropy.units as u
    from tqdm import tqdm

    if verbose: print("Scanning input exposure: " + input_name)
    exposure_identity = rs.utils.exposure_inspector(input_name=input_name)
    if verbose: print("Done.")


    max_size_of_a_star = 501
    g_mag_max = 19


    main_folder = os.getcwd() # os.path.dirname(input_name)

    sci_exts = exposure_identity["SCIEXTS"] # To be detected automatically with exposure inspector

    # If the header 0 contains PSFHARVS, then the image has already been analyzed, and we should jump it.
    input_fits = fits.open(input_name)
    try:
        print(input_fits[0].header["PSFHARVS"])
        print(input_name + " has been PSF Harvested. Jumping.")
        return()
    except:
        print("Analyzing " + input_name + "...")

    # For each SCIENCE extension we run the following loop.
    for ext in sci_exts:
        ###########################
        # Fundamental estimations #
        ###########################
        data, wcs = rs.utils.get_data_and_wcs(input_name=input_name, ext=ext)
        star_search_radius = rs.utils.find_max_angular_size_of_image(data, wcs)
        ###########################

        # Make a segmentation map to find out the saturated sources
        segmentation_map = rs.gnu.make_gnuastro_segmentation_map(input_name, ext, saturation_level=50000)


        # Identify the stars
        # At this point, we can have the approximate coordinates of the stars.
        # But in real images, the centers might differ.
        # find_stars_inside_detector cross-match the stars in the image and the ones in the Gaia catalog
        # to give a better position of the centroids.

        matched_stars_in_detector = rs.psf.find_stars_inside_detector(input_name=input_name,
                                                                      ext=ext,
                                                                      lambda_ref=exposure_identity["PHOTPLAM"]*u.AA,
                                                                      radius=star_search_radius,
                                                                      g_mag_max=g_mag_max,
                                                                      MJD=exposure_identity["EXPSTART"],
                                                                      clean=clean,
                                                                      verbose=verbose)
        if verbose: print(matched_stars_in_detector)
        ra_gaia_list = matched_stars_in_detector["matched_catalog_image"]["RA_GAIA"]
        dec_gaia_list = matched_stars_in_detector["matched_catalog_image"]["DEC_GAIA"]
        ra_ima_list = matched_stars_in_detector["matched_catalog_image"]["RA_IMA"]
        dec_ima_list = matched_stars_in_detector["matched_catalog_image"]["DEC_IMA"]
        x_ima_list = matched_stars_in_detector["matched_catalog_image"]["X_IMA"]
        y_ima_list = matched_stars_in_detector["matched_catalog_image"]["Y_IMA"]
        mag_list = matched_stars_in_detector["matched_catalog_image"]["MAG_LAMBDA"]

        sci_ima_name     = rs.gnu.make_wcs_clean_extension(input_name, ext)
        sky_corrected_db = rs.sky.correct_flat_sky(input_name=sci_ima_name, ext=1, clean=clean, verbose=verbose)

        basename = os.path.basename(sci_ima_name)
        psf_name = sci_ima_name.replace(".fits", "_psf.fits")

        ################################
        # Now we analyze star by star.
        ################################

        if not os.path.exists(main_folder + "/psf_harvester/"):
            stdout = os.mkdir(main_folder + "/psf_harvester/")
        if not os.path.exists(main_folder + "/psf_harvester/psf_level_0"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_0")
        if not os.path.exists(main_folder + "/psf_harvester/psf_level_0/stamps"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_0/stamps")

        if not os.path.exists(main_folder + "/psf_harvester/psf_level_1"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_1")
        if not os.path.exists(main_folder + "/psf_harvester/psf_level_1/stamps"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_1/stamps")

        if not os.path.exists(main_folder + "/psf_harvester/psf_level_2"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_2")
        if not os.path.exists(main_folder + "/psf_harvester/psf_level_2/stamps"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_2/stamps")

        if not os.path.exists(main_folder + "/psf_harvester/psf_level_3"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_3")
        if not os.path.exists(main_folder + "/psf_harvester/psf_level_3/stamps"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_3/stamps")

        outname_list = []
        rlim = np.zeros(len(ra_ima_list))
        rsat = np.zeros(len(ra_ima_list))
        norm_radii_min = np.zeros(len(ra_ima_list))
        norm_radii_max = np.zeros(len(ra_ima_list))
        psf_level = np.zeros(len(ra_ima_list))

        nstars = len(ra_ima_list)


        for i in tqdm(range(nstars)):
            ra = ra_ima_list[i]
            dec = dec_ima_list[i]
            x_ima = x_ima_list[i]
            y_ima = y_ima_list[i]
            mag = mag_list[i]

            # Measuring the normalization radius for each star
            measure_maxradii_result = rs.psf.measure_maxradii(input_name=sci_ima_name, ra=ra, dec=dec,
                                           rmax=max_size_of_a_star, clean=clean,
                                           verbose=verbose, saturation_level=50000, clean_profile=clean)

            rlim[i] = measure_maxradii_result["rlim"]
            rsat[i] = measure_maxradii_result["rsat"]


            if rlim[i] is None:
                rlim[i] = np.nan
                rsat[i] = np.nan
                norm_radii_min[i] = np.nan
                norm_radii_max[i] = np.nan
                outname_list.append(None)

            else: # If the star is found, make a stamp.
                # Three different saturation levels.
                # This has to remain flexible. We will add more saturation levels
                # if we find brigther stars.
                if rlim[i] < 10: # If we do not detect more than 10 pixels, then the star is not detected.
                    psf_level[i] = int(0)

                elif np.isnan(rsat[i]):
                    psf_level[i] = int(1)
                    norm_radii_min[i], norm_radii_max[i] = 5, 8

                elif (rsat[i] >= 0) and (rsat[i] < 5):
                    psf_level[i] = int(2)
                    norm_radii_min[i], norm_radii_max[i] = 10, 15

                elif (rsat[i] >= 5):
                    psf_level[i] = int(3)
                    norm_radii_min[i], norm_radii_max[i] = 20, 30



                outname = main_folder + "/psf_harvester/psf_level_" + str(int(psf_level[i])) + "/stamps/"+basename.replace(".fits","")+"_"+str(i).zfill(5)+"_mag"+str(np.round(mag,2))+".fits"
                outname_list.append(outname)

                cmd = "astscript-psf-stamp "+ sci_ima_name + " " +\
                      "--mode=img --normradii="+str(norm_radii_min[i])+","+str(norm_radii_max[i])+\
                      " --center="+str(x_ima)+","+str(y_ima)+" "+\
                      " --widthinpix="+str(max_size_of_a_star)+","+str(max_size_of_a_star)+" "+\
                      " --segment="+ segmentation_map + " " +\
                      " --output="+ outname
                stdout = rs.utils.execute_cmd(cmd)

                # Write the info on the header of the psf stamp
                cmd = "astfits -h1  " + outname+\
                      " --update=EXPOSURE,"+str(input_name)+\
                      " --update=EXT,"+str(ext)+\
                      " --update=RA_GAIA,"+str(ra_gaia_list[i])+\
                      " --update=DEC_GAIA,"+str(dec_gaia_list[i])+\
                      " --update=RA_IMA,"+str(ra_ima_list[i])+\
                      " --update=DEC_IMA,"+str(dec_ima_list[i])+\
                      " --update=X_IMA,"+str(x_ima_list[i])+\
                      " --update=Y_IMA,"+str(y_ima_list[i])+\
                      " --update=MAG_WAVE,"+str(mag_list[i])+\
                      " --update=RLIM,"+str(rlim[i])+\
                      " --update=RSAT,"+str(rsat[i])+\
                      " --update=PSF_LVL,"+str(psf_level[i])+\
                      " --update=NRM_RMIN,"+str(norm_radii_min[i])+\
                      " --update=NRM_RMAX,"+str(norm_radii_max[i])+\
                      " --update=FILTER,"+str(exposure_identity["FILTER"])+\
                      " --update=EXPTIME,"+str(exposure_identity["EXPTIME"])

                stdout = rs.utils.execute_cmd(cmd)

        if clean:
            input_directory = os.path.dirname(input_name)
            rs.utils.execute_cmd("rm " + input_directory + "/*_ext*.fits")
            rs.utils.execute_cmd("rm " + input_directory + "/*_gaia.dat")
            rs.utils.execute_cmd("rm " + input_directory + "/*_swarp_coadd.fits")

    # Save in the header a new keyword, to indicate that this image has been analyzed already.
    cmd = "astfits -h0  " + input_name + " --update=PSFHARVS,1"
    stdout = rs.utils.execute_cmd(cmd)

##########################

def find_stars_inside_detector(input_name, g_mag_max=15, verbose=False):

    import os
    radius = 1 # To find the stars inside the detector, we do not need more than 1 degree.

    # If the input name is an ASDF, transform it to a compiled FITS.
    if isinstance(input_name, (list, np.ndarray, pd.Series)):
        if input_name[0].split(".")[-1] == "asdf":
            fits_input_name = input_name[0].replace(".asdf",".fits")
            if verbose:
                print(input_name)
            fits_input_name = rs.utils.convert_ASDF_to_FITS(asdf_list = input_name, output=fits_input_name)

    else:
        fits_input_name = input_name


    # Get the image identity
    image_identity = rs.utils.exposure_inspector(input_name=fits_input_name, verbose=verbose, lite=True)


    # Get the star catalog:
    source_catalog_filename = fits_input_name.replace(".fits", "_source_catalog.csv")

    #if not os.path.exists(source_catalog_filename):
    hybrid_catalog = rs.psf.get_hybrid_catalog(ra=image_identity["RA_TARG"],
                                               dec=image_identity["DEC_TARG"],
                                               radius=radius,
                                               lambda_ref=image_identity["FILTER_IDENTITY"]["filter_lambda_ref"],
                                               MJD=image_identity["EXPSTART"],
                                               observer=image_identity["TELESCOP"],
                                               g_mag_max = g_mag_max,
                                               verbose=verbose,
                                               query_filename=source_catalog_filename)
    #else:
    #    hybrid_catalog = pd.read_csv(source_catalog_filename)

    from tqdm import tqdm
    names_of_bool_columns_if_star_is_inside = []
    for SCIEXT_i in tqdm(image_identity["SCIEXTS"]):
        if verbose: print("> Identifying which stars are inside the FOV and which are outside...")
        infield_stars = rs.psf.identify_stars_in_out_field(data_shape=image_identity["DATA_SHAPE"][SCIEXT_i-1],
                                                           wcs=image_identity["ASTROPYWCS"][SCIEXT_i-1],
                                                           catalog=hybrid_catalog,
                                                           verbose=verbose)

        name_column_is_star_inside_this_detector = "in_SCI" + str(SCIEXT_i)
        names_of_bool_columns_if_star_is_inside.append(name_column_is_star_inside_this_detector)

        hybrid_catalog[name_column_is_star_inside_this_detector] = infield_stars["bool_isIn"]
    # Once you are done checking if the stars are inside each detector,
    # find out which stars are outside ALL detectors.
    hybrid_catalog["is_inside_FPA"] = hybrid_catalog[names_of_bool_columns_if_star_is_inside].any(axis=1)

    return(hybrid_catalog)

##########################

def getWCS_galsim_dict_style(file_name):
    # This program takes the wcs in the galsim fashion, to be used by
    # galsim.roman.findSCA
    import galsim
    galsim_wcs_dict = {}
    for i in range(18):
        galsim_wcs_dict[i+1] = galsim.GSFitsWCS(file_name=file_name, hdu=i+1)
    return(galsim_wcs_dict)

##########################

def find_SCA_for_a_target(file_name, ra, dec, include_border=True):
    # This program returns the SCA ID for a target. If not in the detector, is nan.
    import coord
    import galsim
    import galsim.roman as galsim_roman

    wcs_dict = getWCS_galsim_dict_style(file_name=file_name)
    nstars = len(ra)
    IN_SCA = np.zeros(nstars)
    from tqdm import tqdm

    ra1 = np.array([wcs_dict[1].center.ra.deg])
    dec1 = np.array([wcs_dict[1].center.dec.deg])
    # Measure the distance to the target. If it is too far off, it is a NAN.
    distance = rs.utils.angular_distance(ra1=ra1, dec1=dec1, ra2=ra, dec2=dec)
    too_far_off = distance > 1.2
    too_far_off = too_far_off[0,:]
    print(too_far_off.shape)
    for i in tqdm(range(nstars)):
        if too_far_off[i]:
            IN_SCA[i] = np.nan
        else:
            star_pos = galsim.CelestialCoord(ra=ra[i]*coord.degrees,
                                         dec=dec[i]*coord.degrees)
            IN_SCA[i] = galsim_roman.findSCA(wcs_dict=wcs_dict,
                                         world_pos=star_pos,
                                         include_border=include_border)
    return(IN_SCA)
##########################

def generate_star_in_footprint(hybrid_catalog, sciexts, astropywcs, telescope='Roman', filter='F129', verbose=False):
    """
    Deprecated version of the generate_star_stamps function.
    This function generates star stamps for each SCIEXT in the input catalog.
    """
    canvas = None




    return canvas


def generate_star_stamps_deprecated(hybrid_catalog, telescope, filename, sciexts, astropywcs, filter, pa, verbose=False):
    """Generate individual PSF-based FITS stamps for catalog stars.

    Parameters
    ----------
    hybrid_catalog : pandas.DataFrame
        Catalog containing at least ``ra``, ``dec``, and ``mag_lambda`` columns.
    telescope : str
        Telescope name used to select the telescope-specific PSF implementation.
    filename : str
        Input FITS exposure used to determine each star's detector/SCA.
    sciexts : iterable of int
        Science-extension or detector identifiers to process.
    astropywcs : sequence of astropy.wcs.WCS
        WCS objects corresponding to the science extensions.
    filter : str
        Roman filter name used to generate the PSF and calculate stellar flux.
    pa : float
        Position angle in degrees written to each generated stamp's WCS.
    verbose : bool or int, optional
        Controls progress and diagnostic output. Values greater than 1 enable
        more detailed progress reporting.

    Returns
    -------
    list of list of str
        Output FITS filenames grouped by science extension.

    Notes
    -----
    This function is deprecated. It creates a directory named
    ``<filename without .fits>_star_stamps`` and generates one FITS PSF stamp
    per catalog star. Existing stamp files are reused.
    """

    # For each SCIEXT.
    foldername = filename.replace(".fits", "_star_stamps")
    
    # os.system("rm -r " + foldername)
    # os.system("mkdir " + foldername)

    stars_outnames = []
    from tqdm import tqdm
    #  Find the telescope class for the get_psf function
    telescope_class = rs.telescopes.telescope_class_finder(telescope=telescope)

    if verbose > 0: print("Generating stamps for infield stars...")

    print("We might need to get rid of stars very far away.")
    print(np.array(hybrid_catalog["ra"])) #
    # Find the right SCA for each star.
    IN_SCA = rs.psf.find_SCA_for_a_target(file_name=filename,
                                   ra=np.array(hybrid_catalog["ra"]),
                                   dec=np.array(hybrid_catalog["dec"]),
                                   include_border=True)



    rs.utils.execute_cmd("mkdir " + foldername)

    for SCIEXT_i in tqdm(sciexts):
        stars_sciext_outnames = []
        rs.utils.execute_cmd("mkdir " + foldername + "/SCIEXT_" + str(SCIEXT_i).zfill(2))
        ra_stars = np.array(hybrid_catalog["ra"][IN_SCA == SCIEXT_i])
        dec_stars =  np.array(hybrid_catalog["dec"][IN_SCA == SCIEXT_i])
        mag_stars =  np.array(hybrid_catalog["mag_lambda"][IN_SCA == SCIEXT_i])
        astropywcs_sca = astropywcs[SCIEXT_i-1]
        xcen, ycen = astropywcs_sca.wcs_world2pix(ra_stars, dec_stars, 0)

        # Get the flux of each star in electrons per second
        fe_stars = rs.roman.mag2fe(mag=mag_stars, bandpass=filter, sca=SCIEXT_i)


        for i in tqdm(range(len(ra_stars)), disable=(verbose < 2)):
            outname = foldername + "/SCIEXT_" + str(SCIEXT_i).zfill(2) + "/star_ext" + str(SCIEXT_i).zfill(2) + "_id_" + str(i).zfill(6) + ".fits"
            if os.path.exists(outname):
                if verbose > 0: print("Star stamp already exists: " + outname)
                stars_sciext_outnames.append(outname)
                continue
            
            psf_array = telescope_class.get_psf(detector_position=(xcen[i], ycen[i]),
                                                detector = SCIEXT_i,
                                                filter_name = filter)
            norm_psf = np.nansum(psf_array.data)

            pixscale = 0.11/60/60
            psf_array_shape = psf_array.data.shape
            header = rs.utils.create_custom_wcs(crpix=[int(psf_array_shape[1]/2), int(psf_array_shape[0]/2)],
                                                crval=[ra_stars[i], dec_stars[i]],
                                                cdelt=[-pixscale,pixscale],
                                                crota=[-pa,-pa])

            psf_array.data = (psf_array.data/norm_psf)*fe_stars[i]


            rs.utils.save_fits(psf_array.data.value, outname, header)

            stars_sciext_outnames.append(outname)

        stars_outnames.append(stars_sciext_outnames)
    return(stars_outnames)

##########################

def moffat_function_sb(r, mu0, alpha, beta):
    """
    Moffat surface brightness profile in magnitudes per square arcsecond.

    Parameters
    ----------
        r : float or np.ndarray
            Radial distance from the center.
        mu0 : float
            Central surface brightness.
        alpha : float
            Scale parameter.
        beta : float
            Shape parameter.

    Returns
    -------
        float or np.ndarray: 
            Surface brightness at radius r.
    """
    return mu0 + 2.5 * beta * np.log10(1 + (r/alpha)**2)

def galsim_roman_psf(position, SCA, bandpass, shape=None, SED=None, **kwargs):
    '''
    Generate a Roman WFI point-spread function using GalSim.

    Parameters
    ----------
        position : sequence of float
            Two-dimensional detector coordinates ``(x, y)`` in pixels.
        SCA : int
            Roman WFI SCA number, from 1 through 18.
        bandpass : str
            Roman WFI filter name, for example ``"F158"``.
        shape : tuple of int, optional
            The shape of the output PSF image in pixels, as ``(ny, nx)``. If
            ``None``, the default size used by GalSim's ``drawImage`` is used.
        SED : galsim.SED or None, optional
            The spectral energy distribution to use for the PSF. If ``None``, a plane SED is used.
        method : str, optional
            The method to use for PSF generation. Default is ``'stpsf'``.
        kwargs : dict, optional
            Additional keyword arguments to pass to the galsim.GSParams.

    Returns
    -------
        astropy.io.fits.PrimaryHDU
            The rendered PSF image in the primary HDU.
    '''
    from romanisim.bandpass import roman2galsim_bandpass
    import galsim
    import galsim.roman

    # Transforms bandpass formating to galsim format and to the banspass galsim object
    bandpass = roman2galsim_bandpass[bandpass]
    bandpass_galsim = galsim.roman.getBandpasses()[bandpass]

    # Set up the GalSim parameters for high-accuracy PSF rendering and override defaults with any additional kwargs
    gsparams_kwargs = dict(maximum_fft_size=2**16, maxk_threshold=1.e-10, folding_threshold=1.e-6, **kwargs)
    gsparams = galsim.GSParams(**gsparams_kwargs)

    # Set the SCA position for the PSF rendering
    SCA_pos = galsim.PositionD(position[0], position[1])

    # Call galsim to generate the PSF for the specified SCA, bandpass, and position
    psf = galsim.roman.getPSF(SCA, bandpass, SCA_pos=SCA_pos, pupil_bin=1, wcs=None,
                                n_waves=10, extra_aberrations=None,
                                gsparams=gsparams, logger=None,
                                high_accuracy=None, approximate_struts=None)

    # Define the default flat SED if none is provided
    if SED is None:
        SED = galsim.SED(spec=lambda wave: 1.0, wave_type='nm', flux_type='fnu')
    elif not isinstance(SED, galsim.SED):
        raise TypeError("SED must be an instance of galsim.SED")

    # Convolve the PSF with the SED using a delta function to create the final image profile
    psf = galsim.Convolve([psf, SED*galsim.DeltaFunction()])

    # Draw the final image using the convolved PSF and the specified bandpass
    img = psf.drawImage(nx=shape[1], ny=shape[0], bandpass=bandpass_galsim) if shape is not None else psf.drawImage(bandpass=bandpass_galsim)
    output = fits.PrimaryHDU(data=img.array, header=img.wcs.writeToFitsHeader(fits.Header(), img.bounds))
    
    return output

def stpsf_roman_psf(position, SCA, bandpass, shape=None, SED=None,
                    magnitude_weights=None, wave_weights=None, oversample_factor=1,
                    flux_weights=None, zp=23.9,
                    pixel_scale=0.11, **kwargs):
    """Render a flux-calibrated, chromatic Roman/WFI PSF with STPSF.

    Monochromatic STPSF PSFs are interpolated over the requested filter and
    combined using the supplied spectrum. Each image pixel is in nJy/pixel;
    the image sum is conserved to the requested integrated flux, including
    when a rectangular output is cropped.

    Parameters
    ----------
    position : sequence of float
        Detector position ``(x, y)`` in pixels.
    SCA : int
        Roman WFI SCA number (1--18).
    bandpass : str
        Roman WFI filter name, for example ``"F158"``.
    shape : tuple of int, optional
        Output shape ``(ny, nx)``. Defaults to 91 by 91 pixels.
    SED : galsim.SED, optional
        Optional SED whose values are interpreted as spectral ``f_nu`` in
        cgs units. If supplied without explicit weights, its band-averaged
        flux is converted to nJy and distributed across the filter.
    magnitude_weights : array-like, optional
        AB magnitudes at each wavelength in ``wave_weights``. These are
        converted to nJy/pixel using ``zp`` and ``pixel_scale``. Prefer
        ``flux_weights`` when flux contributions are already available.
    wave_weights : array-like or astropy.units.Quantity, optional
        Wavelengths corresponding to ``magnitude_weights`` or
        ``flux_weights``. Unitless values are interpreted as nanometers.
    flux_weights : array-like, optional
        Per-wavelength flux contributions in nJy/pixel. They are summed
        directly and should pair with ``wave_weights``. This is the preferred
        input for :meth:`Star.get_linear_weights` output.
    zp : float, optional
        AB zero point for converting magnitudes to nJy. Defaults to 23.9.
    pixel_scale : float, optional
        Pixel scale in arcsec/pixel, used to convert surface flux density to
        nJy/pixel for magnitude and SED inputs. Defaults to 0.11.
    **kwargs
        Parameters passed to ``galsim.GSParams`` for PSF interpolation setup.

    Returns
    -------
    astropy.io.fits.PrimaryHDU
        Chromatic PSF in nJy/pixel. Header includes ``BUNIT``, ``MAGZP``,
        ``PIXSCALE``, and ``ABMAG`` (the integrated surface-brightness
        magnitude corresponding to the pixel sum).

    Notes
    -----
    STPSF expects monochromatic wavelengths in meters; filter curves and
    ``wave_weights`` are handled in nanometers internally. 
    """

    import stpsf
    from scipy.interpolate import interp1d

    if magnitude_weights is not None and flux_weights is not None:
        raise ValueError("Provide either magnitude_weights or flux_weights, not both.")
    if (magnitude_weights is not None or flux_weights is not None) and wave_weights is None:
        raise ValueError("wave_weights must be provided with magnitude_weights or flux_weights.")
    if SED is not None and (magnitude_weights is not None or flux_weights is not None):
        raise ValueError("SED cannot be combined with explicit magnitude/flux weights.")
    
    # Get the filter response for the specified bandpass
    bp = rs.telescopes.Roman.get_filter(instrument="WFI", filter_name=bandpass)
    filter_wave_nm = np.asarray(bp["wavelength_bins"].to(u.nm).value, dtype=float)
    filter_throughput = np.asarray(bp["transmission_bins"], dtype=float)
    order = np.argsort(filter_wave_nm)
    filter_wave_nm = filter_wave_nm[order]
    filter_throughput = filter_throughput[order]
    valid_filter = np.isfinite(filter_wave_nm) & np.isfinite(filter_throughput)
    filter_wave_nm = filter_wave_nm[valid_filter]
    filter_throughput = np.clip(filter_throughput[valid_filter], 0.0, None)

    # Ensure the filter has a valid wavelength grid
    if len(filter_wave_nm) < 2 or np.any(np.diff(filter_wave_nm) <= 0):
        raise ValueError(f"Filter {bandpass!r} has an invalid wavelength grid.")
    
    filter_response = filter_throughput * np.gradient(filter_wave_nm)
    if not np.isfinite(filter_response.sum()) or filter_response.sum() <= 0:
        raise ValueError(f"Filter {bandpass!r} has no positive transmission.")

    # Measure the weights for the PSF based on the provided wave_weights, magnitude_weights, flux_weights, or SED
    if wave_weights is not None:
        if isinstance(wave_weights, u.Quantity):
            sample_wave_nm = np.asarray(wave_weights.to(u.nm).value, dtype=float)
        else:
            sample_wave_nm = np.asarray(wave_weights, dtype=float)
        
        if sample_wave_nm.ndim != 1:
            raise ValueError("wave_weights must be a one-dimensional wavelength array.")
        
        if magnitude_weights is not None:
            sample_magnitudes = np.asarray(magnitude_weights, dtype=float)
            
            # Check that the magnitude array matches the wavelength array in shape
            if sample_magnitudes.shape != sample_wave_nm.shape:
                raise ValueError("magnitude_weights and wave_weights must have matching shapes.")
            
            throughput = np.interp(sample_wave_nm, filter_wave_nm, filter_throughput, left=0, right=0)
            response = throughput * np.abs(np.gradient(sample_wave_nm))
            sample_flux = 10**(0.4 * (zp - sample_magnitudes)) * pixel_scale**2
            weights = sample_flux * response / response.sum()
        else:
            weights = np.asarray(flux_weights, dtype=float)
            if weights.shape != sample_wave_nm.shape:
                raise ValueError("flux_weights and wave_weights must have matching shapes.")
    
    elif SED is not None:
        sample_wave_nm = filter_wave_nm
        spectral_fnu = np.asarray(SED(sample_wave_nm), dtype=float)
        
        if spectral_fnu.shape != sample_wave_nm.shape:
            raise ValueError("SED must return one f_nu value for each sampled wavelength.")
        
        # GalSim f_nu is in erg/s/cm^2/Hz; 1 nJy = 1e-32 in these units.
        sample_flux = spectral_fnu / 1e-32
        response = filter_response / filter_response.sum()
        weights = sample_flux * response * pixel_scale**2
    else:
        # A flat 1 nJy spectrum is a useful unit-flux default.
        sample_wave_nm = filter_wave_nm
        weights = filter_response / filter_response.sum()

    # Filter out invalid or non-positive samples
    valid_samples = np.isfinite(sample_wave_nm) & np.isfinite(weights) & (weights >= 0)
    sample_wave_nm = sample_wave_nm[valid_samples]
    weights = weights[valid_samples]
    if len(sample_wave_nm) == 0 or weights.sum() <= 0:
        raise ValueError("The wavelength/flux weights contain no positive finite samples.")
    if np.any(sample_wave_nm < filter_wave_nm[0]) or np.any(sample_wave_nm > filter_wave_nm[-1]):
        raise ValueError(f"wave_weights must lie within the {bandpass} filter curve.")

    # Determine the final PSF stamp shape, field of view, and position
    if shape is None:
        ny = nx = 91
    else:
        if len(shape) != 2 or min(shape) <= 0:
            raise ValueError("shape must be a pair of positive dimensions (ny, nx).")
        ny, nx = map(int, shape)
    fov = max(ny, nx)

    # Get the PSF from the STPSF Roman WFI model
    wfi = stpsf.roman.WFI()
    wfi.filter = bandpass
    wfi.detector = f"SCA{int(SCA):02d}"
    wfi.detector_position = np.clip(position, 0, 4095) 

    # Define the wavelength anchors for the chromatic PSF integration
    wave_anchors_nm = np.linspace(filter_wave_nm[0], filter_wave_nm[-1], 10)
    psf_cube = np.empty((len(wave_anchors_nm), fov, fov), dtype=float)
    
    for i, wave_nm in enumerate(wave_anchors_nm):
        # STPSF's monochromatic parameter is in meters (not nanometers).
        psf_hdul = wfi.calc_psf(fov_pixels=fov, monochromatic=wave_nm * 1e-9, oversample=1)
        monochromatic_psf = np.asarray(psf_hdul[0].data, dtype=float)
        psf_sum = np.nansum(monochromatic_psf)
        
        if not np.isfinite(psf_sum) or psf_sum <= 0:
            raise RuntimeError(f"STPSF returned an invalid PSF at {wave_nm:.3f} nm.")
        psf_cube[i] = np.nan_to_num(monochromatic_psf / psf_sum)

    # Create an interpolator for the chromatic PSF based on the wavelength anchors
    psf_interpolator = interp1d(wave_anchors_nm, psf_cube, axis=0, kind="linear",
                                bounds_error=False, fill_value="extrapolate")

    # Combbine the monochromatic PSFs weighted by the flux at each wavelength
    chromatic_psf = np.zeros((fov, fov), dtype=float)
    for wave_nm, flux in zip(sample_wave_nm, weights):
        chromatic_psf += psf_interpolator(wave_nm) * flux

    
    # Normalize the chromatic PSF to match the requested total flux
    requested_flux = float(np.sum(weights))
    image_flux = float(np.sum(chromatic_psf))
    if image_flux <= 0 or not np.isfinite(image_flux):
        raise RuntimeError("Chromatic PSF integration produced no finite positive flux.")
    chromatic_psf *= requested_flux / image_flux

    # Save header information and create the FITS HDU for the chromatic PSF to match ROSALIA
    hdu = fits.PrimaryHDU(data=chromatic_psf)
    hdu.header["TELESCOP"] = "ROMAN"
    hdu.header["INSTRUME"] = "WFI"
    hdu.header["DETECTOR"] = wfi.detector
    hdu.header["FILTER"] = bandpass
    hdu.header["X_POS"] = position[0]
    hdu.header["Y_POS"] = position[1]
    hdu.header["BUNIT"] = "nJy/pixel"
    hdu.header["MAGZP"] = (float(zp), "AB zero point for flux values in nJy")
    hdu.header["PIXSCALE"] = (float(pixel_scale), "Pixel scale in arcsec/pixel")
    hdu.header["TOTFLUX"] = (requested_flux, "Integrated PSF flux in nJy")
    hdu.header["ABMAG"] = (float(zp - 2.5 * np.log10(requested_flux / pixel_scale**2)),
                            "Integrated AB magnitude per square arcsec")
    return hdu



def superback_roman_psf(position, SCA, bandpass, shape=None, SED=None,
                    magnitude_weights=None, wave_weights=None, oversample_factor=1,
                    flux_weights=None, zp=23.9,
                    pixel_scale=0.11, **kwargs):
    """Render a flux-calibrated, chromatic Roman/WFI PSF with SUPERBACK analytical model.

    Monochromatic STPSF PSFs are interpolated over the requested filter and
    combined using the supplied spectrum. Each image pixel is in nJy/pixel;
    the image sum is conserved to the requested integrated flux, including
    when a rectangular output is cropped.

    Parameters
    ----------
    position : sequence of float
        Detector position ``(x, y)`` in pixels.
    SCA : int
        Roman WFI SCA number (1--18).
    bandpass : str
        Roman WFI filter name, for example ``"F158"``.
    shape : tuple of int, optional
        Output shape ``(ny, nx)``. Defaults to 91 by 91 pixels.
    SED : galsim.SED, optional
        Optional SED whose values are interpreted as spectral ``f_nu`` in
        cgs units. If supplied without explicit weights, its band-averaged
        flux is converted to nJy and distributed across the filter.
    magnitude_weights : array-like, optional
        AB magnitudes at each wavelength in ``wave_weights``. These are
        converted to nJy/pixel using ``zp`` and ``pixel_scale``. Prefer
        ``flux_weights`` when flux contributions are already available.
    wave_weights : array-like or astropy.units.Quantity, optional
        Wavelengths corresponding to ``magnitude_weights`` or
        ``flux_weights``. Unitless values are interpreted as nanometers.
    flux_weights : array-like, optional
        Per-wavelength flux contributions in nJy/pixel. They are summed
        directly and should pair with ``wave_weights``. This is the preferred
        input for :meth:`Star.get_linear_weights` output.
    zp : float, optional
        AB zero point for converting magnitudes to nJy. Defaults to 23.9.
    pixel_scale : float, optional
        Pixel scale in arcsec/pixel, used to convert surface flux density to
        nJy/pixel for magnitude and SED inputs. Defaults to 0.11.
    **kwargs
        Parameters passed to ``galsim.GSParams`` for PSF interpolation setup.

    Returns
    -------
    astropy.io.fits.PrimaryHDU
        Chromatic PSF in nJy/pixel. Header includes ``BUNIT``, ``MAGZP``,
        ``PIXSCALE``, and ``ABMAG`` (the integrated surface-brightness
        magnitude corresponding to the pixel sum).

    Notes
    -----
    STPSF expects monochromatic wavelengths in meters; filter curves and
    ``wave_weights`` are handled in nanometers internally. 
    """

    import rosalia.superback as superback
    from scipy.interpolate import interp1d

    if magnitude_weights is not None and flux_weights is not None:
        raise ValueError("Provide either magnitude_weights or flux_weights, not both.")
    if (magnitude_weights is not None or flux_weights is not None) and wave_weights is None:
        raise ValueError("wave_weights must be provided with magnitude_weights or flux_weights.")
    if SED is not None and (magnitude_weights is not None or flux_weights is not None):
        raise ValueError("SED cannot be combined with explicit magnitude/flux weights.")
    
    # Get the filter response for the specified bandpass
    bp = rs.telescopes.Roman.get_filter(instrument="WFI", filter_name=bandpass)
    filter_wave_nm = np.asarray(bp["wavelength_bins"].to(u.nm).value, dtype=float)
    filter_throughput = np.asarray(bp["transmission_bins"], dtype=float)
    order = np.argsort(filter_wave_nm)
    filter_wave_nm = filter_wave_nm[order]
    filter_throughput = filter_throughput[order]
    valid_filter = np.isfinite(filter_wave_nm) & np.isfinite(filter_throughput)
    filter_wave_nm = filter_wave_nm[valid_filter]
    filter_throughput = np.clip(filter_throughput[valid_filter], 0.0, None)

    # Ensure the filter has a valid wavelength grid
    if len(filter_wave_nm) < 2 or np.any(np.diff(filter_wave_nm) <= 0):
        raise ValueError(f"Filter {bandpass!r} has an invalid wavelength grid.")
    
    filter_response = filter_throughput * np.gradient(filter_wave_nm)
    if not np.isfinite(filter_response.sum()) or filter_response.sum() <= 0:
        raise ValueError(f"Filter {bandpass!r} has no positive transmission.")

    # Measure the weights for the PSF based on the provided wave_weights, magnitude_weights, flux_weights, or SED
    if wave_weights is not None:
        if isinstance(wave_weights, u.Quantity):
            sample_wave_nm = np.asarray(wave_weights.to(u.nm).value, dtype=float)
        else:
            sample_wave_nm = np.asarray(wave_weights, dtype=float)
        
        if sample_wave_nm.ndim != 1:
            raise ValueError("wave_weights must be a one-dimensional wavelength array.")
        
        if magnitude_weights is not None:
            sample_magnitudes = np.asarray(magnitude_weights, dtype=float)
            
            # Check that the magnitude array matches the wavelength array in shape
            if sample_magnitudes.shape != sample_wave_nm.shape:
                raise ValueError("magnitude_weights and wave_weights must have matching shapes.")
            
            throughput = np.interp(sample_wave_nm, filter_wave_nm, filter_throughput, left=0, right=0)
            response = throughput * np.abs(np.gradient(sample_wave_nm))
            sample_flux = 10**(0.4 * (zp - sample_magnitudes)) * pixel_scale**2
            weights = sample_flux * response / response.sum()
        else:
            weights = np.asarray(flux_weights, dtype=float)
            if weights.shape != sample_wave_nm.shape:
                raise ValueError("flux_weights and wave_weights must have matching shapes.")
    
    elif SED is not None:
        sample_wave_nm = filter_wave_nm
        spectral_fnu = np.asarray(SED(sample_wave_nm), dtype=float)
        
        if spectral_fnu.shape != sample_wave_nm.shape:
            raise ValueError("SED must return one f_nu value for each sampled wavelength.")
        
        # GalSim f_nu is in erg/s/cm^2/Hz; 1 nJy = 1e-32 in these units.
        sample_flux = spectral_fnu / 1e-32
        response = filter_response / filter_response.sum()
        weights = sample_flux * response * pixel_scale**2
    else:
        # A flat 1 nJy spectrum is a useful unit-flux default.
        sample_wave_nm = filter_wave_nm
        weights = filter_response / filter_response.sum()

    # Filter out invalid or non-positive samples
    valid_samples = np.isfinite(sample_wave_nm) & np.isfinite(weights) & (weights >= 0)
    sample_wave_nm = sample_wave_nm[valid_samples]
    weights = weights[valid_samples]
    if len(sample_wave_nm) == 0 or weights.sum() <= 0:
        raise ValueError("The wavelength/flux weights contain no positive finite samples.")
    if np.any(sample_wave_nm < filter_wave_nm[0]) or np.any(sample_wave_nm > filter_wave_nm[-1]):
        raise ValueError(f"wave_weights must lie within the {bandpass} filter curve.")

    # Determine the final PSF stamp shape and field of view
    if shape is None:
        ny = nx = 91
    else:
        if len(shape) != 2 or min(shape) <= 0:
            raise ValueError("shape must be a pair of positive dimensions (ny, nx).")
        ny, nx = map(int, shape)
    fov = max(ny, nx)

    # Get the thetax, thetay coordinates for the source position
    thetax, thetay = get_detector_position_telescope_frame(SCA, idl_x=position[0], idl_y=position[1])
    x_src, y_src = create_psf_grid(fov, scale_arcsec=0.108, oversample_factor=oversample_factor)

    # Get the PSF from the STPSF Roman WFI model
    bard = superback.getbardict(thetax, thetay)

    # Define the wavelength anchors for the chromatic PSF integration
    wave_anchors_nm = np.linspace(filter_wave_nm[0], filter_wave_nm[-1], 10)
    psf_cube = np.empty((len(wave_anchors_nm), fov, fov), dtype=float)
    
    for i, wave_nm in enumerate(wave_anchors_nm):
        # STPSF's monochromatic parameter is in meters (not nanometers).
        psf_hdul = np.abs(superback.romanpsf(x_src, y_src, wave_nm * 1e-9, barparam=bard))**2
        monochromatic_psf = np.asarray(psf_hdul, dtype=float)
        psf_sum = np.nansum(monochromatic_psf)
        
        if not np.isfinite(psf_sum) or psf_sum <= 0:
            raise RuntimeError(f"STPSF returned an invalid PSF at {wave_nm:.3f} nm.")
        psf_cube[i] = np.nan_to_num(monochromatic_psf / psf_sum)

    # Create an interpolator for the chromatic PSF based on the wavelength anchors
    psf_interpolator = interp1d(wave_anchors_nm, psf_cube, axis=0, kind="linear",
                                bounds_error=False, fill_value="extrapolate")

    # Combbine the monochromatic PSFs weighted by the flux at each wavelength
    chromatic_psf = np.zeros((fov, fov), dtype=float)
    for wave_nm, flux in zip(sample_wave_nm, weights):
        chromatic_psf += psf_interpolator(wave_nm) * flux

    
    # Normalize the chromatic PSF to match the requested total flux
    requested_flux = float(np.sum(weights))
    image_flux = float(np.sum(chromatic_psf))
    if image_flux <= 0 or not np.isfinite(image_flux):
        raise RuntimeError("Chromatic PSF integration produced no finite positive flux.")
    chromatic_psf *= requested_flux / image_flux

    # Save header information and create the FITS HDU for the chromatic PSF to match ROSALIA
    hdu = fits.PrimaryHDU(data=chromatic_psf)
    hdu.header["TELESCOP"] = "ROMAN"
    hdu.header["INSTRUME"] = "WFI"
    hdu.header["DETECTOR"] = f'WFI{SCA:02d}_FULL'
    hdu.header["FILTER"] = bandpass
    hdu.header["X_POS"] = position[0]
    hdu.header["Y_POS"] = position[1]
    hdu.header["BUNIT"] = "nJy/pixel"
    hdu.header["MAGZP"] = (float(zp), "AB zero point for flux values in nJy")
    hdu.header["PIXSCALE"] = (float(pixel_scale), "Pixel scale in arcsec/pixel")
    hdu.header["TOTFLUX"] = (requested_flux, "Integrated PSF flux in nJy")
    hdu.header["ABMAG"] = (float(zp - 2.5 * np.log10(requested_flux / pixel_scale**2)),
                            "Integrated AB magnitude per square arcsec")
    return hdu



def moffat_2d(shape, mu0, alpha, beta):
    """
    Generate a 2D Moffat profile.

    Parameters
    ----------
        shape : tuple of int
            Shape of the 2D array to generate.
        mu0 : float
            Central surface brightness.
        alpha : float
            Scale parameter.
        beta : float
            Shape parameter.

    Returns
    -------
        np.ndarray
            2D array representing the surface brightness at each point.
    """
    # Create grid of given shape
    y, x = np.indices(shape)
    y_cen = (shape[0] - 1) / 2
    x_cen = (shape[1] - 1) / 2
    r = np.sqrt((x - x_cen)**2 + (y - y_cen)**2)
    return mu0 + 2.5 * beta * np.log10(1 + (r/alpha)**2)

def get_psf_extension_moffat(mag, depth=30, alpha=0.08, beta=1.52):
    """
    Analytically calculate the PSF extension based on a standard Moffat profile at a given depth.
        
    Parameters
    ----------
        mag: float or array-like
            Magnitude(s) of the star.
        depth: float, optional
            Depth at the extension in mag arcsecond^-2. Default is 30.
        alpha : float, optional
            Scale parameter of the Moffat profile. Default is 0.08.
        beta : float, optional
            Shape parameter of the Moffat profile. Default is 1.52.
        
    Returns
    -------
        float or np.ndarray
            Exact PSF extension in arcseconds.
    """
    mag = np.asarray(mag)

    # Measure the mu0 from integrating the moffat function analytically
    term = (np.pi * alpha**2) / (beta - 1)
    mu0 = mag + 2.5 * np.log10(term)
    
    # Invert the Moffat surface brightness equation to solve for r
    exponent = (depth - mu0) / (2.5 * beta)
    term = 10**exponent - 1
    
    # Prevent np.sqrt from encountering negative numbers if depth < mag
    term = np.maximum(term, 0) 
    
    extension = alpha * np.sqrt(term)
    
    return extension.item() if extension.ndim == 0 else extension

def build_psf_canvas(stars_catalog, shape_out, bandpass, depth=31, n_workers=None,
                     method="superback"):
    '''
    Build a PSF canvas from a catalog of stars in parallel.

    Parameters
    ----------
        stars_catalog: pandas.DataFrame
            Catalog of stars containing at least columns:
                `ra`, `dec`, `x`, `y`, `detector_id`, `x_det`, `y_det`
        shape_out: tuple
            Shape of the output PSF canvas (ny, nx).
        bandpass: str
            Bandpass of the observation in Roman WFI filter names (e.g., F062, F087, F106, F129, F158, F184, F213).
        depth: float, optional
            Depth at which to build the PSF in mag arcsecond^-2. Default is 31.
        n_workers: int, optional
            Maximum number of parallel workers. Defaults to at most two to limit
            the memory used by concurrently generated PSF models.
        method: str, optional
            PSF model to use. Defaults to ``"superback"``.

    Returns
    -------
        np.ndarray
            PSF canvas of shape `shape_out`.

    Notes
    -----
        The positions of each detector follow the SIAF convention in pixel coordinates.
        
    '''

    if not isinstance(stars_catalog, pd.DataFrame):
        raise TypeError("stars_catalog must be a pandas DataFrame.")
    if len(shape_out) != 2 or any(int(size) <= 0 for size in shape_out):
        raise ValueError("shape_out must contain two positive dimensions (ny, nx).")
    shape_out = tuple(int(size) for size in shape_out)
    if not np.isfinite(depth):
        raise ValueError("depth must be finite.")
    if not isinstance(bandpass, str) or not bandpass:
        raise ValueError("bandpass must be a non-empty filter name.")

    required_columns = {
        "ra", "dec", "source_id", "detector_id", "x_det", "y_det", "x", "y",
        "phot_bp_mean_mag_AB", "phot_g_mean_mag_AB", "phot_rp_mean_mag_AB",
        "phot_w1_mean_mag_AB", "phot_w2_mean_mag_AB", "phot_w3_mean_mag_AB",
        "phot_w4_mean_mag_AB", "phot_j_mean_mag_AB", "phot_h_mean_mag_AB",
        "phot_ks_mean_mag_AB",
    }
    missing_columns = required_columns.difference(stars_catalog.columns)
    if missing_columns:
        raise ValueError(
            "stars_catalog is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    canvas = np.zeros(shape_out, dtype=float)
    if stars_catalog.empty:
        return canvas

    if n_workers is None:
        n_workers = min(2, os.cpu_count() or 1)
    if isinstance(n_workers, bool) or not isinstance(n_workers, (int, np.integer)) or n_workers < 1:
        raise ValueError("n_workers must be a positive integer or None.")
    n_workers = min(int(n_workers), len(stars_catalog))

    from concurrent.futures import ThreadPoolExecutor
    from threading import Lock

    canvas_lock = Lock()
    with tqdm(
        total=len(stars_catalog),
        desc="Injecting star PSFs",
        unit="star",
        position=0,
        leave=True,
    ) as progress:
        worker_catalogs = [
            stars_catalog.iloc[offset::n_workers]
            for offset in range(n_workers)
        ]
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            futures = [
                executor.submit(
                    _build_psf_for_worker,
                    worker_catalog,
                    depth,
                    bandpass,
                    canvas,
                    canvas_lock,
                    method,
                    progress,
                )
                for worker_catalog in worker_catalogs
            ]
            for future in futures:
                future.result()

    return canvas


def _build_psf_for_worker(stars_catalog, depth, bandpass, canvas, canvas_lock,
                          method="superback", progress=None, pixel_scale=0.11):
    '''
    Generate and inject the PSFs assigned to a worker into the shared canvas.

    Parameters
    ----------
        stars_catalog: pandas.DataFrame
            Catalog of stars assigned to this worker.
        depth: float
            Depth at which to build the PSF in mag arcsecond^-2.
        bandpass: str
            Roman WFI bandpass used to compute the stellar spectrum.
        canvas: np.ndarray
            Shared output canvas. Updates are protected by ``canvas_lock``.
        canvas_lock: threading.Lock
            Lock protecting overlapping PSF additions to the canvas.
        progress: tqdm.tqdm, optional
            Shared progress bar advanced after each PSF is injected.
    '''
    for _, row in stars_catalog.iterrows():
        x, y = float(row["x"]), float(row["y"])
        x_det, y_det = float(row["x_det"]), float(row["y_det"])
        detector_id = int(row["detector_id"])
        if not np.all(np.isfinite([x, y, x_det, y_det])):
            raise ValueError(f"Star {row['source_id']!r} has non-finite pixel coordinates.")
        if detector_id < 1 or detector_id > 18 or detector_id != row["detector_id"]:
            raise ValueError(f"Star {row['source_id']!r} has invalid detector_id {row['detector_id']!r}.")

        star = Star(row)
        magnitude = star.get_magnitudes(bandpass)
        extension = get_psf_extension_moffat(magnitude, depth=depth)
        if not np.isfinite(extension):
            raise ValueError(f"Could not determine a finite PSF size for star {row['source_id']!r}.")

        radius_pixels = max(1, int(np.ceil(extension / pixel_scale)))
        stamp_size = 2 * radius_pixels + 1
        wave_weights, flux_weights = star.get_linear_weights(bandpass)
        psf_model = rs.telescopes.Roman.get_psf(
            position=(x_det, y_det),
            SCA=detector_id,
            bandpass=bandpass,
            shape=(stamp_size, stamp_size),
            method=method,
            wave_weights=wave_weights,
            flux_weights=flux_weights,
        )
        psf_data = np.asarray(psf_model.data, dtype=float)

        center_x, center_y = int(np.rint(x)), int(np.rint(y))
        x_start = center_x - psf_data.shape[1] // 2
        y_start = center_y - psf_data.shape[0] // 2
        x_end = x_start + psf_data.shape[1]
        y_end = y_start + psf_data.shape[0]
        canvas_x_start, canvas_y_start = max(0, x_start), max(0, y_start)
        canvas_x_end = min(canvas.shape[1], x_end)
        canvas_y_end = min(canvas.shape[0], y_end)

        if canvas_x_start < canvas_x_end and canvas_y_start < canvas_y_end:
            psf_x_start = canvas_x_start - x_start
            psf_y_start = canvas_y_start - y_start
            psf_x_end = psf_x_start + canvas_x_end - canvas_x_start
            psf_y_end = psf_y_start + canvas_y_end - canvas_y_start
            with canvas_lock:
                canvas[canvas_y_start:canvas_y_end, canvas_x_start:canvas_x_end] += (
                    psf_data[psf_y_start:psf_y_end, psf_x_start:psf_x_end]
                )

        del psf_data, psf_model, wave_weights, flux_weights, star
        if progress is not None:
            progress.update(1)


def get_detector_position_telescope_frame(detector_id, idl_x=0, idl_y=0):
    """
    Convert detector position to telescope frame coordinates.
    
    Parameters
    ----------
    detector_id : str
        Detector identifier (e.g., 'WFI01')
    idl_x : float
        X position in detector ideal frame (pixels)
    idl_y : float
        Y position in detector ideal frame (pixels)
        
    Returns
    -------
    thetax : float
        V2 position in telescope frame (radians)
    thetay : float
        V3 position in telescope frame (radians)
    """
    import stpsf

    siaf = stpsf.stpsf_core.get_siaf_with_caching('roman')
    aperture_name = f'WFI{detector_id:02d}_FULL'
    
    if aperture_name not in siaf.apertures:
        raise ValueError(f"Aperture {aperture_name} not found in SIAF")
    
    aperture = siaf.apertures[aperture_name]
    thetax_deg, thetay_deg = aperture.idl_to_tel(idl_x, idl_y)
    
    # Convert from degrees to radians
    thetax_rad = thetax_deg * np.pi / 180 / 60 / 60
    thetay_rad = thetay_deg * np.pi / 180 / 60 / 60
    
    return thetax_rad, thetay_rad


def create_psf_grid(size_pixels, scale_arcsec, oversample_factor=1):
    """
    Create coordinate grids for PSF calculation.
    
    Parameters
    ----------
    size_pixels : int
        PSF array size (will be size_pixels x size_pixels)
    scale_arcsec : float
        Pixel scale in arcsec/pixel
    oversample_factor : float
        Sub-pixel sampling factor
        
    Returns
    -------
    x : ndarray
        X coordinate grid in arcsec
    y : ndarray
        Y coordinate grid in arcsec
    """
    # Create fine grid with oversampling
    fine_size = size_pixels * oversample_factor
    fine_scale = scale_arcsec / oversample_factor
    
    # Centered at (0, 0), extend ±half the field
    half_field = (size_pixels * scale_arcsec) / 2
    coords = np.linspace(-half_field + fine_scale/2, half_field - fine_scale/2, fine_size)
    
    x, y = np.meshgrid(coords, coords)
    return x, y