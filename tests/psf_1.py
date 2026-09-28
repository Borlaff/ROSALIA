import rosalia as rs
import numpy as np

import warnings
warnings.filterwarnings('ignore')

row = {'Unnamed: 0.2': 823,
 'healpix_lvl': 125448.0,
 'source_id': 4413833440790780672,
 'ra': 229.7654312757245,
 'dec': -3.851853628629195,
 'phot_g_mean_mag_AB': 9.666927337646484,
 'phot_bp_mean_mag_AB': 10.327932357788086,
 'phot_rp_mean_mag_AB': 9.064289093017578,
 'phot_j_mean_mag_AB': 8.369778633117676,
 'phot_h_mean_mag_AB': 8.155250549316406,
 'phot_ks_mean_mag_AB': 8.408172607421875,
 'phot_w1_mean_mag_AB': 9.191,
 'phot_w2_mean_mag_AB': 9.955,
 'phot_w3_mean_mag_AB': 11.686,
 'phot_w4_mean_mag_AB': 13.049,
 'dist': 0.4290934237250933,
 'mag_lambda': 8.207610060848406,
 'cat_id': 181975,
 'detector_id': 18}

star = rs.psf.Star(row)

wave_weights, fnu_weights = star.get_linear_weights(bandpass='F158')
mag = star.get_magnitudes(bandpass='F158')

psf = rs.psf.stpsf_roman_psf((200,1300), 3, 'F158', shape=(300,300), magnitude_weights=wave_weights, wave_weights=fnu_weights)


tot_flux = np.sum(psf.data)
tot_mag = 23.9 -2.5 * np.log10(tot_flux) + 5*np.log10(0.11)

print('Total PSF flux:', tot_flux)
print('Total PSF magnitude:', tot_mag)
print('Star magnitude:', mag)
