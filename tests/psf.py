# This notebook demonstrates how to use ROSALIA to analyze psf in a simple image. 

import rosalia as rs
import numpy as np


STAR_ROW = {
    "source_id": 4413833440790780672,
    "ra": 229.7654312757245,
    "dec": -3.851853628629195,
    "phot_g_mean_mag_AB": 9.666927337646484,
    "phot_bp_mean_mag_AB": 10.327932357788086,
    "phot_rp_mean_mag_AB": 9.064289093017578,
    "phot_j_mean_mag_AB": 8.369778633117676,
    "phot_h_mean_mag_AB": 8.155250549316406,
    "phot_ks_mean_mag_AB": 8.408172607421875,
    "phot_w1_mean_mag_AB": 9.191,
    "phot_w2_mean_mag_AB": 9.955,
    "phot_w3_mean_mag_AB": 11.686,
    "phot_w4_mean_mag_AB": 13.049,
    "detector_id": 18,
}

_STAR_CATALOG = {'Unnamed: 0.2': {200425: 218,
  200409: 168,
  200405: 184,
  200398: 140,
  200384: 17},
 'Unnamed: 0': {200425: 1607.0,
  200409: 1478.0,
  200405: 1552.0,
  200398: 1329.0,
  200384: 272.0},
 'healpix_lvl': {200425: 167054.0,
  200409: 167076.0,
  200405: 167076.0,
  200398: 167076.0,
  200384: 167055.0},
 'source_id': {200425: '5877722672321796992',
  200409: '5878467728852386816',
  200405: '5878472642295175296',
  200398: '5878482503541385472',
  200384: '5877727551405514368'},
 'ra': {200425: 219.74011998126645,
  200409: 219.2932760228554,
  200405: 219.33920282717625,
  200398: 219.087138667354,
  200384: 220.51380856169752},
 'dec': {200425: -60.95482487005107,
  200409: -61.00864921257298,
  200405: -60.85332787840406,
  200398: -60.91373367889644,
  200384: -61.070269275939786},
 'wise_ra': {200425: 219.7402086,
  200409: 219.2932471,
  200405: 219.3393476,
  200398: 219.0871592,
  200384: 220.5138304},
 'wise_dec': {200425: -60.9548361,
  200409: -61.0086339,
  200405: -60.8533105,
  200398: -60.9137398,
  200384: -61.0702066},
 'tmass_ra': {200425: 219.740157,
  200409: 219.293353,
  200405: 219.339272,
  200398: 219.087193,
  200384: 220.513877},
 'tmass_dec': {200425: -60.954781,
  200409: -61.008644,
  200405: -60.853336,
  200398: -60.913727,
  200384: -61.070244},
 'phot_bp_mean_mag': {200425: 13.626878,
  200409: 16.235336,
  200405: 16.249563,
  200398: 15.44119,
  200384: 16.099075},
 'phot_g_mean_mag': {200425: 13.359948,
  200409: 15.06591,
  200405: 15.08639,
  200398: 14.497751,
  200384: 14.997941},
 'phot_rp_mean_mag': {200425: 13.062532,
  200409: 14.210691,
  200405: 14.229944,
  200398: 13.725222,
  200384: 14.16255},
 'parallax': {200425: 1.622817956341234,
  200409: 0.2863621251911822,
  200405: 0.3513024725172187,
  200398: 0.288384734336787,
  200384: 0.3232128572314011},
 'w1': {200425: 0.0066926790693584,
  200409: 0.0146554784095591,
  200405: 0.0168888469479091,
  200398: 0.0104039917705135,
  200384: 0.0136269919699486},
 'w1_error': {200425: 0.0, 200409: 0.0, 200405: 0.0, 200398: 0.0, 200384: 0.0},
 'w2': {200425: 0.0028549598169004,
  200409: 0.0064923212314321,
  200405: 0.0075787922091604,
  200398: 0.0052287755531344,
  200384: 0.0061659500186148},
 'w2_error': {200425: 0.0, 200409: 0.0, 200405: 0.0, 200398: 0.0, 200384: 0.0},
 'w3': {200425: 0.0020758692704951,
  200409: 0.0014655478409559,
  200405: 0.0017733726245425,
  200398: 0.0008464474480977,
  200384: 0.0011303163247386},
 'w3_error': {200425: 0.0, 200409: 0.0, 200405: 0.0, 200398: 0.0, 200384: 0.0},
 'w4': {200425: 0.0021737019751571,
  200409: 0.0020970075785446,
  200405: 0.0029593749953433,
  200398: 0.0064923212314321,
  200384: 0.0021066869257394},
 'w4_error': {200425: 0.0, 200409: 0.0, 200405: 0.0, 200398: 0.0, 200384: 0.0},
 'j': {200425: 0.025473498,
  200409: 0.024215225,
  200405: 0.024506915,
  200398: 0.025709163,
  200384: 0.025008546},
 'j_error': {200425: 0.0, 200409: 0.0, 200405: 0.0, 200398: 0.0, 200384: 0.0},
 'h': {200425: 0.020393858,
  200409: 0.034857772,
  200405: 0.031996336,
  200398: 0.023157917,
  200384: 0.033197075},
 'h_error': {200425: 0.0, 200409: 0.0, 200405: 0.0, 200398: 0.0, 200384: 0.0},
 'ks': {200425: 0.013649926,
  200409: 0.029589005,
  200405: 0.027210124,
  200398: 0.018498233,
  200384: 0.026590755},
 'ks_error': {200425: 0.0, 200409: 0.0, 200405: 0.0, 200398: 0.0, 200384: 0.0},
 'g': {200425: 0.016444497,
  200409: 0.0034169264,
  200405: 0.003353081,
  200398: 0.0057663275,
  200384: 0.0036376712},
 'bp': {200425: 0.012860231,
  200409: 0.001163765,
  200405: 0.0011486149,
  200398: 0.0024183767,
  200384: 0.0013193793},
 'rp': {200425: 0.021626528,
  200409: 0.007511439,
  200405: 0.0073794196,
  200398: 0.011746571,
  200384: 0.00785199},
 'phot_g_mean_mag_AB': {200425: 13.359948,
  200409: 15.06591,
  200405: 15.08639,
  200398: 14.497751,
  200384: 14.997941},
 'phot_bp_mean_mag_AB': {200425: 13.626878,
  200409: 16.235336,
  200405: 16.249563,
  200398: 15.44119,
  200384: 16.099075},
 'phot_rp_mean_mag_AB': {200425: 13.062532,
  200409: 14.210691,
  200405: 14.229944,
  200398: 13.725222,
  200384: 14.16255},
 'phot_j_mean_mag_AB': {200425: 12.884778,
  200409: 12.939778,
  200405: 12.926778,
  200398: 12.87478,
  200384: 12.9047785},
 'phot_h_mean_mag_AB': {200425: 13.126251,
  200409: 12.5442505,
  200405: 12.637249,
  200398: 12.988251,
  200384: 12.59725},
 'phot_ks_mean_mag_AB': {200425: 13.562174,
  200409: 12.722174,
  200405: 12.813173,
  200398: 13.232174,
  200384: 12.838173},
 'phot_w1_mean_mag_AB': {200425: 14.336,
  200409: 13.485,
  200405: 13.331,
  200398: 13.857,
  200384: 13.564},
 'phot_w2_mean_mag_AB': {200425: 15.261,
  200409: 14.369,
  200405: 14.201,
  200398: 14.604,
  200384: 14.425},
 'phot_w3_mean_mag_AB': {200425: 15.607,
  200409: 15.985,
  200405: 15.778,
  200398: 16.581,
  200384: 16.267},
 'phot_w4_mean_mag_AB': {200425: 15.557,
  200409: 15.596,
  200405: 15.222,
  200398: 14.369,
  200384: 15.591},
 'dist': {200425: 0.1483676732499288,
  200409: 0.350709163058619,
  200405: 0.2837378044905734,
  200398: 0.4131144632066638,
  200384: 0.3718731633926681},
 'Unnamed: 0.1': {200425: np.nan,
  200409: np.nan,
  200405: np.nan,
  200398: np.nan,
  200384: np.nan},
 'source_id_2': {200425: np.nan,
  200409: np.nan,
  200405: np.nan,
  200398: np.nan,
  200384: np.nan},
 'mag_lambda': {200425: 12.90796320809442,
  200409: 12.90180113552834,
  200405: 12.89897865964158,
  200398: 12.885675001709018,
  200384: 12.875250924380929},
 'cat_id': {200425: 200426,
  200409: 200410,
  200405: 200406,
  200398: 200399,
  200384: 200385}}


def test_superback_roman_psf_integrated_magnitude_matches_star():
    """The flux-calibrated chromatic PSF preserves the Star band magnitude."""
    star = rs.psf.Star(STAR_ROW)
    bandpass = "F158"
    wavelength_nm, flux_weights = star.get_linear_weights(bandpass=bandpass)
    star_magnitude = star.get_magnitudes(bandpass=bandpass)

    psf = rs.psf.superback_roman_psf(
        position=(200, 1300),
        SCA=3,
        bandpass=bandpass,
        shape=(512, 512),
        wave_weights=wavelength_nm,
        flux_weights=flux_weights,
    )

    total_flux_njy_per_pixel = np.sum(psf.data)
    integrated_magnitude = (
        23.9 - 2.5 * np.log10(total_flux_njy_per_pixel)
    )

    assert psf.data.shape == (512, 512)
    assert psf.header["BUNIT"] == "nJy/pixel"
    assert np.isfinite(total_flux_njy_per_pixel)
    assert total_flux_njy_per_pixel > 0
    assert np.isclose(integrated_magnitude, star_magnitude, atol=1e-4)

def test_stpsf_roman_psf_integrated_magnitude_matches_star():
    """The flux-calibrated chromatic PSF preserves the Star band magnitude."""
    star = rs.psf.Star(STAR_ROW)
    bandpass = "F158"
    wavelength_nm, flux_weights = star.get_linear_weights(bandpass=bandpass)
    star_magnitude = star.get_magnitudes(bandpass=bandpass)

    psf = rs.psf.superback_roman_psf(
        position=(200, 1300),
        SCA=3,
        bandpass=bandpass,
        shape=(512, 512),
        wave_weights=wavelength_nm,
        flux_weights=flux_weights,
    )

    total_flux_njy_per_pixel = np.sum(psf.data)
    integrated_magnitude = (
        23.9 - 2.5 * np.log10(total_flux_njy_per_pixel)
    )

    assert psf.data.shape == (512, 512)
    assert psf.header["BUNIT"] == "nJy/pixel"
    assert np.isfinite(total_flux_njy_per_pixel)
    assert total_flux_njy_per_pixel > 0
    assert np.isclose(integrated_magnitude, star_magnitude, atol=1e-4)


def test_exposure_psf_creation():
    # This notebook demonstrates how to use ROSALIA to analyze psf in a simple image. 

    import rosalia as rs
    import pandas as pd
    from astropy.time import Time
    from astropy.coordinates import SkyCoord
    from astropy.wcs import WCS

    # Define the target
    target = 'alpha Centauri'
    tgt = SkyCoord.from_name(target)
    ra =  219.92041  # Right ascension, in degrees. 
    dec = -60.835148 # Declination, in degrees.
    PA = 30 # Position angle, in degrees.
    date = Time("2027-01-01T00:00:00") # Date of the observation, in Astropy Time YYYY-MM-DDTHH:MM:SS format.
    bandpass = "F129"
    exptime = 600 # Exposure time, in seconds.

    observer={"TELESCOP": "Roman/WFI", 
            "pointing": [ra, dec], 
            "FILTER":bandpass, 
            "PA_Y": PA, 
            "EXPSTART": date.mjd, 
            "EXPTIME": exptime}

    prefix = "alphacen"
    custom_roman_exposure = rs.core.exposure(observer=observer, prefix=prefix) 

    psf_out = custom_roman_exposure.psf(catalog=_STAR_CATALOG, n_stars=4, depth=27, n_workers=4)


    ## Test outputs are correct
    assert psf_out is not None
    assert isinstance(psf_out, tuple)
    assert len(psf_out) == 3
    assert psf_out[0] is not None  # PSF HDUList
    for i in range(18):
        assert isinstance(psf_out[0][i], np.ndarray)  # Check that each PSF detector array is indeed a numpy arrayassert isinstance(psf_out[0][0], np.ndarray)  # PSF HDUList
    assert np.sum(psf_out[0]) > 0   # Not empty arrays
    assert len(psf_out[0]) == 18  # 
    assert psf_out[1] is not None  # WCS list
    assert isinstance(psf_out[1], list)  # WCS list
    assert len(psf_out[1]) == 18  # There should be 18 WCS objects, one for each detector.
    for i in range(18):
        assert isinstance(psf_out[1][i], WCS)  # Check that each WCS object is indeed a WCS instance
    assert psf_out[2] is not None  # Canvas geometry
    assert isinstance(psf_out[2], np.ndarray)  # Canvas geometry
    assert np.sum(psf_out[2]) > 0  # Canvas is not empty
