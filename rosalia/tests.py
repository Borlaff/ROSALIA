import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import bottleneck as bn
import pandas as pd
import pytest
from numpy.testing import assert_allclose

# Custom modules
#sys.path.append("/Users/aborlaff/NASA/STRAYCOR/")

import rosalia as rs

def test_fe2mu_mu2fe():
    input_flux = 40
    print("Input flux:" + str(input_flux))
    A = rs.detectors.fe2mu(fe=40, instrument="ACS", filter_name="F606W", telescope="HST", verbose=False)
    print(A)
    B = rs.detectors.mu2fe(mu=A, instrument="ACS", filter_name="F606W", telescope="HST", verbose=False)
    print("Output flux:" + str(B))
    assert_allclose(B.value, input_flux, rtol=0.1)
    print("> test_fe2mu_mu2fe: PASS")
    return(True)

def test_magnitude_conversion_gaia():
    from astroquery.gaia import Gaia

    # Reference guide for the Gaia photometric zeropoints
    # https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_calibr_extern.html

    # Updated Zeropoints in DR3
    # https://gea.esac.esa.int/archive/documentation/GEDR3/Data_processing/chap_cu5pho/cu5pho_sec_photProc/cu5pho_ssec_photCal.html#SSS3.P1
    #gaia_g_AB_zp = 25.7934 # -2.5*log10(e/s) + zp = magnitude
    #gaia_bp_AB_zp = 25.3806 # -2.5*log10(e/s) + zp = magnitude
    #gaia_rp_AB_zp = 25.1161 # -2.5*log10(e/s) + zp = magnitude

    #gaia_g_Vega_zp = 25.6874 # -2.5*log10(e/s) + zp = magnitude
    #gaia_bp_Vega_zp = 25.3385 # -2.5*log10(e/s) + zp = magnitude
    #gaia_rp_Vega_zp = 24.7479 # -2.5*log10(e/s) + zp = magnitude
    # Checking GAIA zeropoint

    ra = 199.5 # deg
    dec =  -22.0 # deg
    radius = 0.5 # deg

    cmd = "SELECT g.ra AS ra,\
       g.dec AS dec,\
       phot_g_mean_flux,\
       phot_bp_mean_flux,\
       phot_rp_mean_flux,\
       phot_g_mean_mag,\
       phot_bp_mean_mag,\
       phot_rp_mean_mag,\
       POWER(10,0.4*(8.9-j_m-0.8937792073497661)) AS flux_j,\
       POWER(10,0.4*(8.9-h_m-1.37425010840047)) AS flux_h,\
       POWER(10,0.4*(8.9-ks_m-1.8401738621860915)) AS flux_ks,\
       POWER(10,0.4*(8.9-w1mpro-2.699)) AS flux_w1,\
       POWER(10,0.4*(8.9-w2mpro-3.339)) AS flux_w2,\
       POWER(10,0.4*(8.9-w3mpro-5.174)) AS flux_w3,\
       POWER(10,0.4*(8.9-w4mpro-6.620)) AS flux_w4\
       FROM gaiadr3.gaia_source AS g\
       JOIN gaiaedr3.tmass_psc_xsc_best_neighbour AS xmatch USING (source_id)\
       JOIN gaiaedr3.tmass_psc_xsc_join AS xjoin USING (clean_tmass_psc_xsc_oid)\
       JOIN gaiadr1.tmass_original_valid AS tmass ON xjoin.original_psc_source_id = tmass.designation\
       INNER JOIN gaiaedr3.allwise_best_neighbour AS wise_x USING (source_id)\
       INNER JOIN gaiadr1.allwise_original_valid AS wise USING(allwise_oid)\
       WHERE g.ra > 199\
       AND g.ra < 201\
       AND g.dec > -23\
       AND g.dec < -22"
    print(cmd)

    job = Gaia.launch_job_async(cmd, dump_to_file=False, verbose=False)
    print("> Async query to Gaia/2MASS/WISE done.")
    filtered_table = job.get_results().to_pandas()
    print(filtered_table)

    g_mag = -2.5*np.log10(filtered_table["phot_g_mean_flux"]) + rs.constants.gaia_g_Vega_zp
    relative_magnitude_g = np.array(g_mag/filtered_table["phot_g_mean_mag"])
    print(relative_magnitude_g)

    bp_mag = -2.5*np.log10(filtered_table["phot_bp_mean_flux"]) + rs.constants.gaia_bp_Vega_zp
    relative_magnitude_bp = np.array(bp_mag/filtered_table["phot_bp_mean_mag"])
    print(relative_magnitude_bp)

    rp_mag = -2.5*np.log10(filtered_table["phot_rp_mean_flux"]) + rs.constants.gaia_rp_Vega_zp
    relative_magnitude_rp = np.array(rp_mag/filtered_table["phot_rp_mean_mag"])
    print(relative_magnitude_rp)

    assert_allclose(bn.nanmedian(relative_magnitude_g), 1, rtol=0.01)
    assert_allclose(bn.nanmedian(relative_magnitude_bp), 1, rtol=0.01)
    assert_allclose(bn.nanmedian(relative_magnitude_rp), 1, rtol=0.01)

    print("> test_magnitude_conversion_gaia: PASS")
    return(True)

def test_photometry_superstars():
    from astroquery.gaia import Gaia

    # First we make a catalog on high-resolution for all the stars inside a region.
    obj_ra = 199.5 # deg
    obj_dec =  -22.0 # deg
    cone_radius = 0.5 # deg
    query_result = rs.gaia.query_gaia_2mass_wise(ra=obj_ra, dec=obj_dec, radius=cone_radius, verbose=False)

    # Then we make a query averaing in healpix cells that contains that area

    cmd = "SELECT GAIA_HEALPIX_INDEX(5, source_id) AS healpix_5,\
       SUM(g.ra*phot_g_mean_flux)/SUM(phot_g_mean_flux) AS geo_RA,\
       SUM(g.dec*phot_g_mean_flux)/SUM(phot_g_mean_flux) AS geo_DEC,\
       AVG(g.ra) AS mean_RA,\
       AVG(g.dec) AS mean_DEC,\
       SUM(phot_g_mean_flux) AS healpix_phot_g_sum_flux,\
       SUM(phot_bp_mean_flux) AS healpix_phot_bp_sum_flux,\
       SUM(phot_rp_mean_flux) AS healpix_phot_rp_sum_flux,\
       SUM(POWER(10,0.4*(8.9-j_m-0.8937792073497661))) AS sum_j,\
       SUM(POWER(10,0.4*(8.9-h_m-1.37425010840047))) AS sum_h,\
       SUM(POWER(10,0.4*(8.9-ks_m-1.8401738621860915))) AS sum_ks,\
       SUM(POWER(10,0.4*(8.9-w1mpro-2.699))) AS sum_w1,\
       SUM(POWER(10,0.4*(8.9-w2mpro-3.339))) AS sum_w2,\
       SUM(POWER(10,0.4*(8.9-w3mpro-5.174))) AS sum_w3,\
       SUM(POWER(10,0.4*(8.9-w4mpro-6.620))) AS sum_w4\
       FROM gaiadr3.gaia_source AS g\
       JOIN gaiaedr3.tmass_psc_xsc_best_neighbour AS xmatch USING (source_id)\
       JOIN gaiaedr3.tmass_psc_xsc_join AS xjoin USING (clean_tmass_psc_xsc_oid)\
       JOIN gaiadr1.tmass_original_valid AS tmass ON xjoin.original_psc_source_id = tmass.designation\
       INNER JOIN gaiaedr3.allwise_best_neighbour AS wise_x USING (source_id)\
       INNER JOIN gaiadr1.allwise_original_valid AS wise USING(allwise_oid)\
       WHERE g.ra > 195\
       AND g.ra < 205\
       AND g.dec > -25\
       AND g.dec < -18\
       GROUP BY healpix_5"

    print(cmd)

    job = Gaia.launch_job_async(cmd, dump_to_file=False, verbose=False)
    print("> Async query to Gaia/2MASS/WISE done.")
    filtered_table = job.get_results().to_pandas()
    print(filtered_table)

    plt.scatter(filtered_table["mean_RA"], filtered_table["mean_DEC"], marker="o", s=filtered_table["healpix_phot_g_sum_flux"]*1E-7, alpha=0.5)
    plt.scatter(filtered_table["geo_RA"], filtered_table["geo_DEC"], marker="o", s=filtered_table["healpix_phot_g_sum_flux"]*1E-7, alpha=0.5)

    plt.show(block=False)

    bool_test_cell_11006 = query_result["healpix_lvl"] == 11006

    # Transform Gaia fluxes to AB magnitudes.

    # We fix the units in the Gaia table
    g_mag_AB =  -2.5*np.log10(filtered_table["healpix_phot_g_sum_flux"]) + rs.constants.gaia_g_AB_zp
    bp_mag_AB =  -2.5*np.log10(filtered_table["healpix_phot_bp_sum_flux"]) + rs.constants.gaia_bp_AB_zp
    rp_mag_AB =  -2.5*np.log10(filtered_table["healpix_phot_rp_sum_flux"]) + rs.constants.gaia_rp_AB_zp

    filtered_table["healpix_phot_g_sum_flux_AB"] = 10**(0.4*(8.9 - g_mag_AB))
    filtered_table["healpix_phot_bp_sum_flux_AB"] = 10**(0.4*(8.9 - bp_mag_AB))
    filtered_table["healpix_phot_rp_sum_flux_AB"] = 10**(0.4*(8.9 - rp_mag_AB))

    # Gaia G band
    flux_superstar_11006 = filtered_table["healpix_phot_g_sum_flux_AB"][filtered_table["healpix_5"]==11006]
    combined_flux_cell_11006_highres = np.sum(query_result["g"][bool_test_cell_11006])
    relative_flux = flux_superstar_11006/combined_flux_cell_11006_highres
    print(relative_flux)
    assert_allclose(relative_flux, 1, rtol=0.1)
    print("> Test healpix_phot_g_sum_flux: PASS")
    # Gaia Bp band
    flux_superstar_11006 = filtered_table["healpix_phot_bp_sum_flux_AB"][filtered_table["healpix_5"]==11006]
    combined_flux_cell_11006_highres = np.sum(query_result["bp"][bool_test_cell_11006])
    relative_flux = flux_superstar_11006/combined_flux_cell_11006_highres
    print(relative_flux)
    assert_allclose(relative_flux, 1, rtol=0.1)
    print("> Test healpix_phot_bp_sum_flux: PASS")
    # Gaia Rp band
    flux_superstar_11006 = filtered_table["healpix_phot_rp_sum_flux_AB"][filtered_table["healpix_5"]==11006]
    combined_flux_cell_11006_highres = np.sum(query_result["rp"][bool_test_cell_11006])
    relative_flux = flux_superstar_11006/combined_flux_cell_11006_highres
    print(relative_flux)
    assert_allclose(relative_flux, 1, rtol=0.1)
    print("> Test healpix_phot_rp_sum_flux: PASS")
    # 2MASS J band
    flux_superstar_11006 = filtered_table["sum_j"][filtered_table["healpix_5"]==11006]
    query_result["phot_j_mean_flux_AB"] = 10**(0.4*(8.9-query_result["phot_j_mean_mag_AB"]))
    combined_flux_cell_11006_highres = np.sum(query_result["j"][bool_test_cell_11006])
    relative_flux = flux_superstar_11006/combined_flux_cell_11006_highres
    print(relative_flux)
    assert_allclose(relative_flux, 1, rtol=0.1)
    print("> Test 2MASS j: PASS")
    # 2MASS H band
    flux_superstar_11006 = filtered_table["sum_h"][filtered_table["healpix_5"]==11006]
    query_result["phot_h_mean_flux_AB"] = 10**(0.4*(8.9-query_result["phot_h_mean_mag_AB"]))
    combined_flux_cell_11006_highres = np.sum(query_result["h"][bool_test_cell_11006])
    relative_flux = flux_superstar_11006/combined_flux_cell_11006_highres
    print(relative_flux)
    assert_allclose(relative_flux, 1, rtol=0.1)
    print("> Test 2MASS h: PASS")
    # 2MASS Ks band
    flux_superstar_11006 = filtered_table["sum_ks"][filtered_table["healpix_5"]==11006]
    query_result["phot_ks_mean_flux_AB"] = 10**(0.4*(8.9-query_result["phot_ks_mean_mag_AB"]))
    combined_flux_cell_11006_highres = np.sum(query_result["ks"][bool_test_cell_11006])
    relative_flux = flux_superstar_11006/combined_flux_cell_11006_highres
    print(relative_flux)
    assert_allclose(relative_flux, 1, rtol=0.1)
    print("> Test 2MASS ks: PASS")
    # WISE W1 band
    flux_superstar_11006 = filtered_table["sum_w1"][filtered_table["healpix_5"]==11006]
    query_result["phot_w1_mean_flux_AB"] = 10**(0.4*(8.9-query_result["phot_w1_mean_mag_AB"]))
    combined_flux_cell_11006_highres = np.sum(query_result["w1"][bool_test_cell_11006])
    relative_flux = flux_superstar_11006/combined_flux_cell_11006_highres
    print(relative_flux)
    assert_allclose(relative_flux, 1, rtol=0.1)
    print("> Test WISE W1: PASS")
    # WISE W2 band
    flux_superstar_11006 = filtered_table["sum_w2"][filtered_table["healpix_5"]==11006]
    query_result["phot_w2_mean_flux_AB"] = 10**(0.4*(8.9-query_result["phot_w2_mean_mag_AB"]))
    combined_flux_cell_11006_highres = np.sum(query_result["w2"][bool_test_cell_11006])
    relative_flux = flux_superstar_11006/combined_flux_cell_11006_highres
    print(relative_flux)
    assert_allclose(relative_flux, 1, rtol=0.1)
    print("> Test WISE W2: PASS")
    # WISE W3 band
    flux_superstar_11006 = filtered_table["sum_w3"][filtered_table["healpix_5"]==11006]
    query_result["phot_w3_mean_flux_AB"] = 10**(0.4*(8.9-query_result["phot_w3_mean_mag_AB"]))
    combined_flux_cell_11006_highres = np.sum(query_result["w3"][bool_test_cell_11006])
    relative_flux = flux_superstar_11006/combined_flux_cell_11006_highres
    print(relative_flux)
    assert_allclose(relative_flux, 1, rtol=0.1)
    print("> Test WISE W3: PASS")
    # WISE W4 band
    flux_superstar_11006 = filtered_table["sum_w4"][filtered_table["healpix_5"]==11006]
    query_result["phot_w4_mean_flux_AB"] = 10**(0.4*(8.9-query_result["phot_w4_mean_mag_AB"]))
    combined_flux_cell_11006_highres = np.sum(query_result["w4"][bool_test_cell_11006])
    relative_flux = flux_superstar_11006/combined_flux_cell_11006_highres
    print(relative_flux)
    assert_allclose(relative_flux, 1, rtol=0.1)
    print("> Test WISE W4: PASS")

    return(True)
