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
        23.9 - 2.5 * np.log10(total_flux_njy_per_pixel) + 5 * np.log10(0.11)
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
        23.9 - 2.5 * np.log10(total_flux_njy_per_pixel) + 5 * np.log10(0.11)
    )

    assert psf.data.shape == (512, 512)
    assert psf.header["BUNIT"] == "nJy/pixel"
    assert np.isfinite(total_flux_njy_per_pixel)
    assert total_flux_njy_per_pixel > 0
    assert np.isclose(integrated_magnitude, star_magnitude, atol=1e-4)