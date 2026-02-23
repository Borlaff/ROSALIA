
import astropy.units as u

############################
# Transformation constants #
############################

MJysr_to_Jyarcsec2 = 2.350443E-5 # 1 MJy sr-1 = 2.350443 × 10-5 Jy arcsec-2.

mapp_sun_V_AB = -26.77
mapp_moon_V_AB = -12.74
flux_moon_jy_V = 10**(-0.4*(mapp_moon_V_AB-8.9))
flux_sun_jy_V = 10**(-0.4*(mapp_sun_V_AB-8.9))
moon_sun_ratio = flux_moon_jy_V/flux_sun_jy_V
thermal_emissivity = 0.9

sigma1 = 0.682689492137086
sigma2 = 0.954499736103642
sigma3 = 0.997300203936740

s1_down_q = (1-sigma1)/2
s1_up_q = 1 - s1_down_q
s2_down_q = (1-sigma2)/2
s2_up_q = 1 - s2_down_q
s3_down_q = (1-sigma3)/2
s3_up_q = 1 - s3_down_q


# Reference guide for the Gaia photometric zeropoints
# https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_calibr_extern.html
# 25.8010 ± 0.0028 	25.3540 ±  0.0023 	25.1040 ± 0.0016
# Updated Zeropoints in DR3
# https://gea.esac.esa.int/archive/documentation/GEDR3/Data_processing/chap_cu5pho/cu5pho_sec_photProc/cu5pho_ssec_photCal.html#SSS3.P1
gaia_g_AB_zp =  25.8010 # -2.5*log10(e/s) + zp = magnitude
gaia_bp_AB_zp = 25.3540 # -2.5*log10(e/s) + zp = magnitude
gaia_rp_AB_zp = 25.1040  # -2.5*log10(e/s) + zp = magnitude

gaia_g_Vega_zp = 25.6874 # -2.5*log10(e/s) + zp = magnitude
gaia_bp_Vega_zp = 25.3385 # -2.5*log10(e/s) + zp = magnitude
gaia_rp_Vega_zp = 24.7479 # -2.5*log10(e/s) + zp = magnitude

# Conversion from WISE magnitudes to AB magnitudes. https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html
#Table 3 - Conversion to the AB system
#(mAB = mVega + Δm)
#Band	magnitude offset (Δm)
#W1	2.699
#W2	3.339
#W3	5.174
#W4	6.620

WISE_W1_delta_mag_Vega_to_AB = 2.699
WISE_W2_delta_mag_Vega_to_AB = 3.339
WISE_W3_delta_mag_Vega_to_AB = 5.174
WISE_W4_delta_mag_Vega_to_AB = 6.620

## Conversion from 2MASS magnitudes: https://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
#
# Table 1 - 2MASS Isophotal Bandpasses and Fluxes-for-0-magnitude from Cohen et al. (2003)
# Band 	Lambda (µm) 	Bandwidth (µm) 	Fnu - 0 mag (Jy) 	Flambda - 0 mag (W cm-2 µm-1)
# J 	1.235 ± 0.006 	0.162 ± 0.001 	1594 ± 27.8 	3.129E-13 ± 5.464E-15
# H  	1.662 ± 0.009 	0.251 ± 0.002 	1024 ± 20.0 	1.133E-13 ± 2.212E-15
# Ks 	2.159 ± 0.011 	0.262 ± 0.002 	666.7 ± 12.6 	4.283E-14 ± 8.053E-16

# mag_2mass = -2.5*log10(Flux_jy/Fnu0)
# flux_Jy = fnu0 * 10 ** (-0.4*mag_2mass)
# mag_2mass_AB = -2.5*log10(fluxJy) + 8.9
TwoMASS_fnu0_J = 1594
TwoMASS_fnu0_H = 1024
TwoMASS_fnu0_Ks = 666.7


##########################
# Solar system constants #

r_earth =  6378E+3 * u.m # m
