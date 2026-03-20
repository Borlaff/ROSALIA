# Alejandro S. Borlaff. NASA Ames Research Center. a.s.borlaff@nasa.gov / asborlaff@gmail.com
# June 14, 2024s.
#
# ROSALIA/GAIA module
# This module will hold all the programs related to the queries to the Gaia / 2MASS / WISE database
#
# Version log:
# v.1.0 - 14 June 2024. First loading of programs inherited from former psf.pys
#
##########################################################

############################
import os
import sys
import pandas as pd
import numpy as np
import healpy as hp
from tqdm import tqdm
import bottleneck as bn
from astropy.io import fits
import astropy.wcs as wcs
import matplotlib.pyplot as plt
from astropy.io import ascii
import astropy.units as u
from astropy.coordinates import match_coordinates_sky, ICRS, Angle, SkyCoord
from astropy.table import Table, join
import astropy_healpix as astro_hp
import rosalia as rs
from rosalia.plots import style
import rosalia.constants as rs_constants

# Suppress warnings. Comment this out if you wish to see the warning messages
import warnings
warnings.filterwarnings('ignore')


##############
def query_gaia_2mass_wise(ra, dec, radius, g_mag_max=False, verbose=False, query_filename="gaia_query.dat"):

    # C is the central coordinates of the pointing.
    c = SkyCoord(ra, dec, frame='icrs', unit="deg")
    from astroquery.gaia import Gaia

    # Initialize the tables
    if verbose: print("Loading Gaia databases...")
    tables = Gaia.load_tables(only_names=True)


    # Healpix cells in level 5 containing the stars:
    hp = astro_hp.HEALPix(nside=128, order='nested', frame=ICRS())
    healpix_ids_query_all = hp.cone_search_skycoord(c, radius=radius*u.deg)


    # Check if the healpix cells have already been queried.
    import os
    cache_gaia_path = os.environ["ROSALIACACHE"]  + '/cache/gaia_queries/'

    healpix_ids_query = []
    healpix_cached_pds = []
    for cell_i in healpix_ids_query_all:
        potential_saved_cell_i_name = cache_gaia_path + 'hp_' + str(cell_i).zfill(12) + '_gmax_' + str(g_mag_max).zfill(3) + '.csv'
        if os.path.exists(potential_saved_cell_i_name):
            # print("cell exists!")
            saved_cell_i_db = pd.read_csv(potential_saved_cell_i_name)
            healpix_cached_pds.append(saved_cell_i_db)
        else:
            healpix_ids_query.append(cell_i)

    # If all the cells have been queried, then open and merge the saved tables:
    if len(healpix_ids_query) == 0:
        return(pd.concat(healpix_cached_pds))

    healpix_radec_query = hp.healpix_to_skycoord(healpix_ids_query)
    healpix_level_lvl_sql_line = "WHERE ("
    for healpix_cell in healpix_ids_query:
        healpix_level_lvl_sql_line = healpix_level_lvl_sql_line + "sub_level_2.healpix_lvl=" + str(int(healpix_cell)) + " OR "
    healpix_level_lvl_sql_line = healpix_level_lvl_sql_line[:-4]
    healpix_level_lvl_sql_line = healpix_level_lvl_sql_line + ") "

    if verbose:
        print(healpix_level_lvl_sql_line)

    if verbose: print("Finding the RA DEC constraints for the Gaia server query...")
    bounding_circle_in_sql_subquery_level_2 = find_ra_dec_constraints(ra=ra, dec=dec, radius=radius)["sql_search_string"]

    # Here we design another filter for the Gaia query. If the user does not want to include all the stars
    # set g_mag_max to a lower value, so it only includes stars brighter than g_mag_max. This is,
    # only stars with g mag < g_mag_max.

    if g_mag_max is not False:
        if verbose: print("> Limiting the Gaia query to stars brighter than G=" + str(g_mag_max) + " mag")
        g_mag_max_string = " AND sub_level_2.phot_g_mean_mag <= " + str(g_mag_max) + " "
    else:
        g_mag_max_string = ""

    cmd = "SELECT sub.healpix_lvl,\
       sub.source_id, wise_x.source_id, sub.ra, sub.dec,\
       wise.ra as wise_ra, wise.dec as wise_dec,\
       tmass.ra as tmass_ra, tmass.dec as tmass_dec,\
       sub.phot_bp_mean_mag, sub.phot_g_mean_mag,\
       sub.phot_rp_mean_mag, sub.parallax, w1mpro, w1mpro_error,\
       w2mpro, w2mpro_error, w3mpro, w3mpro_error, w4mpro, w4mpro_error,\
       j_m, j_msigcom, h_m, h_msigcom, ks_m, ks_msigcom\
    FROM (\
	    SELECT sub_level_2.healpix_lvl,\
		sub_level_2.source_id,\
		sub_level_2.ra,\
		sub_level_2.dec,\
		sub_level_2.phot_bp_mean_mag,\
		sub_level_2.phot_g_mean_mag,\
		sub_level_2.phot_rp_mean_mag,\
		sub_level_2.parallax\
	    FROM (\
			SELECT GAIA_HEALPIX_INDEX(7, source_id) AS healpix_lvl, source_id, ra, dec,    phot_bp_mean_mag,    phot_g_mean_mag,    phot_rp_mean_mag, parallax\
		    FROM gaiadr3.gaia_source_lite AS g\
            " + bounding_circle_in_sql_subquery_level_2 + "\
        OFFSET 0) AS sub_level_2\
		" + healpix_level_lvl_sql_line + "\
        " + g_mag_max_string + "\
    OFFSET 0) AS sub\
    JOIN gaiaedr3.tmass_psc_xsc_best_neighbour AS xmatch USING (source_id)\
    JOIN gaiaedr3.tmass_psc_xsc_join AS xjoin USING (clean_tmass_psc_xsc_oid)\
    JOIN gaiadr1.tmass_original_valid AS tmass\
       ON xjoin.original_psc_source_id = tmass.designation\
    INNER JOIN gaiaedr3.allwise_best_neighbour AS wise_x USING (source_id)\
    INNER JOIN gaiadr1.allwise_original_valid AS wise USING(allwise_oid)"

    if verbose > 3:
        print("Running Gaia query!")
        print(cmd)

    #job = Gaia.launch_job_async(cmd, dump_to_file=False, verbose=verbose)
    print(style.BLUE)
    if verbose < 2:
        verbose_gaia = False
    else:
        verbose_gaia = True
    job = Gaia.launch_job_async(cmd, dump_to_file=False, verbose=verbose_gaia)
    print(style.RESET)
    if verbose: print("> Query to Gaia/2MASS/WISE databases finished.")
    filtered_table = job.get_results().to_pandas()

    # filtered_table.to_csv(query_filename.replace(".csv", "_temp.csv")) # Debugging

    # Correct the magnitudes to the AB system.

    # Gaia

    filtered_table["phot_g_mean_flux"] = 10**(0.4*(rs_constants.gaia_g_Vega_zp - filtered_table["phot_g_mean_mag"]))
    filtered_table["phot_bp_mean_flux"] = 10**(0.4*(rs_constants.gaia_bp_Vega_zp - filtered_table["phot_bp_mean_mag"]))
    filtered_table["phot_rp_mean_flux"] = 10**(0.4*(rs_constants.gaia_rp_Vega_zp - filtered_table["phot_rp_mean_mag"]))

    filtered_table["phot_g_mean_mag_AB"] = -2.5*np.log10(filtered_table["phot_g_mean_flux"]) + rs_constants.gaia_g_AB_zp
    filtered_table["phot_bp_mean_mag_AB"] = -2.5*np.log10(filtered_table["phot_bp_mean_flux"]) + rs_constants.gaia_bp_AB_zp
    filtered_table["phot_rp_mean_mag_AB"] = -2.5*np.log10(filtered_table["phot_rp_mean_flux"]) + rs_constants.gaia_rp_AB_zp

    # 2MASS
    #flux_Jy_J = TwoMASS_fnu0_J*10**(-0.4*filtered_table["j_m"])
    #flux_Jy_H = TwoMASS_fnu0_H*10**(-0.4*filtered_table["h_m"])
    #flux_Jy_Ks = TwoMASS_fnu0_Ks*10 **(-0.4*filtered_table["ks_m"])

    filtered_table["phot_j_mean_mag_AB"] = filtered_table["j_m"] + 8.9 - 2.5*np.log10(rs_constants.TwoMASS_fnu0_J) #-2.5*np.log10(flux_Jy_J) + 8.9
    filtered_table["phot_h_mean_mag_AB"] = filtered_table["h_m"] + 8.9 - 2.5*np.log10(rs_constants.TwoMASS_fnu0_H) #-2.5*np.log10(flux_Jy_H) + 8.9
    filtered_table["phot_ks_mean_mag_AB"] = filtered_table["ks_m"] + 8.9 - 2.5*np.log10(rs_constants.TwoMASS_fnu0_Ks) #-2.5*np.log10(flux_Jy_Ks) + 8.9

    # WISE

    filtered_table["phot_w1_mean_mag_AB"] = filtered_table["w1mpro"] + rs_constants.WISE_W1_delta_mag_Vega_to_AB
    filtered_table["phot_w2_mean_mag_AB"] = filtered_table["w2mpro"] + rs_constants.WISE_W2_delta_mag_Vega_to_AB
    filtered_table["phot_w3_mean_mag_AB"] = filtered_table["w3mpro"] + rs_constants.WISE_W3_delta_mag_Vega_to_AB
    filtered_table["phot_w4_mean_mag_AB"] = filtered_table["w4mpro"] + rs_constants.WISE_W4_delta_mag_Vega_to_AB



    final_query_table = pd.DataFrame({"healpix_lvl": filtered_table["healpix_lvl"],
                                     "source_id": filtered_table["source_id"],
                                     #"source_id_2":  filtered_table["source_id_2"],
                                  "ra": filtered_table["ra"],
                                  "dec": filtered_table["dec"],
                                  #"wise_ra": filtered_table["wise_ra"],
                                  #"wise_dec": filtered_table["wise_dec"],
                                  #"tmass_ra": filtered_table["tmass_ra"],
                                  #"tmass_dec": filtered_table["tmass_dec"],
                                  #"phot_bp_mean_mag": filtered_table["phot_bp_mean_mag_AB"],
                                  #"phot_g_mean_mag":  filtered_table["phot_g_mean_mag_AB"],
                                  #"phot_rp_mean_mag":  filtered_table["phot_rp_mean_mag_AB"],
                                  #"parallax": filtered_table["parallax"],
                                  #"w1": 10**(0.4*(8.9 - filtered_table["phot_w1_mean_mag_AB"])),
                                  #"w1_error": np.zeros(len(filtered_table)),
                                  #"w2": 10**(0.4*(8.9 - filtered_table["phot_w2_mean_mag_AB"])),
                                  #"w2_error": np.zeros(len(filtered_table)),
                                  #"w3": 10**(0.4*(8.9 - filtered_table["phot_w3_mean_mag_AB"])),
                                  #"w3_error": np.zeros(len(filtered_table)),
                                  #"w4": 10**(0.4*(8.9 - filtered_table["phot_w4_mean_mag_AB"])),
                                  #"w4_error": np.zeros(len(filtered_table)),
                                  #"j":  10**(0.4*(8.9 - filtered_table["phot_j_mean_mag_AB"])),
                                  #"j_error": np.zeros(len(filtered_table)),
                                  #"h":  10**(0.4*(8.9 - filtered_table["phot_h_mean_mag_AB"])),
                                  #"h_error": np.zeros(len(filtered_table)),
                                  #"ks":  10**(0.4*(8.9 - filtered_table["phot_ks_mean_mag_AB"])),
                                  #"ks_error": np.zeros(len(filtered_table)),
                                  #"g": 10**(0.4*(8.9 - filtered_table["phot_g_mean_mag_AB"])),
                                  #"bp": 10**(0.4*(8.9 - filtered_table["phot_bp_mean_mag_AB"])),
                                  #"rp": 10**(0.4*(8.9 - filtered_table["phot_rp_mean_mag_AB"])),
                                  "phot_g_mean_mag_AB":   filtered_table["phot_g_mean_mag_AB"],
                                  "phot_bp_mean_mag_AB":  filtered_table["phot_bp_mean_mag_AB"],
                                  "phot_rp_mean_mag_AB":  filtered_table["phot_rp_mean_mag_AB"],
                                  "phot_j_mean_mag_AB":   filtered_table["phot_j_mean_mag_AB"],
                                  "phot_h_mean_mag_AB":   filtered_table["phot_h_mean_mag_AB"],
                                  "phot_ks_mean_mag_AB":  filtered_table["phot_ks_mean_mag_AB"],
                                  "phot_w1_mean_mag_AB":  filtered_table["phot_w1_mean_mag_AB"],
                                  "phot_w2_mean_mag_AB":  filtered_table["phot_w2_mean_mag_AB"],
                                  "phot_w3_mean_mag_AB":  filtered_table["phot_w3_mean_mag_AB"],
                                  "phot_w4_mean_mag_AB":  filtered_table["phot_w4_mean_mag_AB"]})



    all_stars_coords = SkyCoord(np.array(final_query_table["ra"])*u.deg, np.array(final_query_table["dec"])*u.deg, frame='icrs')
    final_query_table["dist"] = c.separation(all_stars_coords).value

    # Merging the new table with the cached tables
    if len(healpix_cached_pds) > 0:
        cached_tables = pd.concat(healpix_cached_pds)
        final_query_table = pd.concat([cached_tables, final_query_table])

    final_query_table.to_csv(query_filename)

    if verbose: print(" Gaia Star by star query result:")
    if verbose: print(final_query_table)

    list_of_healpix_cells_queried = list(set(final_query_table["healpix_lvl"]))

    # If the directory does not exist, create it.
    import os
    cache_gaia_path =  os.environ["ROSALIACACHE"]  + '/cache/gaia_queries/'
    if not os.path.isdir(os.environ["ROSALIACACHE"]  + '/cache'): os.system('mkdir ' + os.environ["ROSALIACACHE"] + '/cache')
    if not os.path.isdir(cache_gaia_path): os.system('mkdir ' + cache_gaia_path)

    for cell_i in list_of_healpix_cells_queried:
        single_cell_db = final_query_table[final_query_table["healpix_lvl"] == cell_i]
        single_cell_db.to_csv(cache_gaia_path + 'hp_' + str(cell_i).zfill(12) + '_gmax_' + str(g_mag_max).zfill(3) + '.csv')

    return(final_query_table)

############################



def query_healpix_ra_slices(healpix_lvl=7, verbose=False):
    # This program is only run in case that the low resolution map is needed.
    # It can take one or two days, depending on the healpix mapping grid
    #healpix_lvl = 7
    from astroquery.gaia import Gaia

    filtered_table_list = []
    ra_list = np.linspace(0, 359, 360)
    for ra_i in tqdm(ra_list):
        #print(ra_i)
        ra_min = ra_i
        ra_max = ra_i + 1
        if os.path.exists("gaia_slice_ra" + str(ra_min) + "_" + str(ra_max) + ".csv"):
            #print("gaia_slice_ra" + str(ra_min) + "_" + str(ra_max) + ".csv found!")
            filtered_table = pd.read_csv("gaia_slice_ra" + str(ra_min) + "_" + str(ra_max) + ".csv")
            filtered_table_list.append(filtered_table)
            continue

        # Gaia fluxes come in photoelectrons / s
        # WISE and 2MASS are in Jy
        cmd = "SELECT GAIA_HEALPIX_INDEX("+str(healpix_lvl)+", source_id) AS healpix_lvl,\
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
        WHERE g.ra >= "+ str(ra_min) + " \
        AND g.ra < "+ str(ra_max) + "\
        GROUP BY healpix_lvl"

        print(cmd)

        job = Gaia.launch_job_async(cmd, dump_to_file=False, verbose=False)
        print("> Async query to Gaia/2MASS/WISE done.")
        filtered_table = job.get_results().to_pandas()
        filtered_table_list.append(filtered_table)
        filtered_table.to_csv("gaia_slice_ra" + str(ra_min) + "_" + str(ra_max) + ".csv")
        if verbose:
            plt.scatter(filtered_table["mean_RA"], filtered_table["mean_DEC"], marker="o", s=filtered_table["healpix_phot_g_sum_flux"]*1E-7, alpha=0.5)
            plt.scatter(filtered_table["geo_RA"], filtered_table["geo_DEC"], marker="o", s=filtered_table["healpix_phot_g_sum_flux"]*1E-7, alpha=0.5)
            plt.show()



    # We combine the list of tables into a single table
    g2w_raw_table = pd.concat(filtered_table_list, axis=0)

    # We fix the units in the Gaia table
    g_mag_AB =  -2.5*np.log10(g2w_raw_table["healpix_phot_g_sum_flux"]) + rs_constants.gaia_g_AB_zp
    bp_mag_AB =  -2.5*np.log10(g2w_raw_table["healpix_phot_bp_sum_flux"]) + rs_constants.gaia_bp_AB_zp
    rp_mag_AB =  -2.5*np.log10(g2w_raw_table["healpix_phot_rp_sum_flux"]) + rs_constants.gaia_rp_AB_zp

    g2w_raw_table["healpix_phot_g_sum_flux"] = 10**(0.4*(8.9 - g_mag_AB))
    g2w_raw_table["healpix_phot_bp_sum_flux"] = 10**(0.4*(8.9 - bp_mag_AB))
    g2w_raw_table["healpix_phot_rp_sum_flux"] = 10**(0.4*(8.9 - rp_mag_AB))
    # The resulting table has subslices in healpix cells. We need to combine them.
    # Basic Healpix operations
    ncells = astro_hp.nside_to_npix(astro_hp.level_to_nside(healpix_lvl))

    healpix_lvl_list = []
    geo_RA_list = []
    geo_DEC_list = []
    mean_RA_list = []
    mean_DEC_list = []
    healpix_phot_g_sum_flux_list = []
    healpix_phot_bp_sum_flux_list = []
    healpix_phot_rp_sum_flux_list = []
    sum_j_list = []
    sum_h_list = []
    sum_ks_list = []
    sum_w1_list = []
    sum_w2_list = []
    sum_w3_list = []
    sum_w4_list = []

    # 'healpix_lvl', 'source_id', 'source_id_2', 'ra', 'dec', 'wise_ra',
    #   'wise_dec', 'tmass_ra', 'tmass_dec', 'phot_bp_mean_mag',
    #   'phot_g_mean_mag', 'phot_rp_mean_mag', 'parallax', 'w1mpro',
    #   'w1mpro_error', 'w2mpro', 'w2mpro_error', 'w3mpro', 'w3mpro_error',
    #   'w4mpro', 'w4mpro_error', 'j_m', 'j_msigcom', 'h_m', 'h_msigcom',
    #   'ks_m', 'ks_msigcom', 'phot_g_mean_flux', 'phot_bp_mean_flux',
    #   'phot_rp_mean_flux', 'phot_g_mean_mag_AB', 'phot_bp_mean_mag_AB',
    #   'phot_rp_mean_mag_AB', 'phot_j_mean_mag_AB', 'phot_h_mean_mag_AB',
    #   'phot_ks_mean_mag_AB', 'phot_w1_mean_mag_AB', 'phot_w2_mean_mag_AB',
    #   'phot_w3_mean_mag_AB', 'phot_w4_mean_mag_AB', 'dist'

    for i in tqdm(range(ncells)):
        subtable_cell_i = g2w_raw_table.iloc[np.where(g2w_raw_table["healpix_lvl"]==i)]
        healpix_lvl_list.append(i)
        geo_RA_list.append(bn.nansum(subtable_cell_i["geo_RA"]*subtable_cell_i["healpix_phot_g_sum_flux"])/bn.nansum(subtable_cell_i["healpix_phot_g_sum_flux"]))
        geo_DEC_list.append(bn.nansum(subtable_cell_i["geo_DEC"]*subtable_cell_i["healpix_phot_g_sum_flux"])/bn.nansum(subtable_cell_i["healpix_phot_g_sum_flux"]))
        mean_RA_list.append(bn.nanmean(subtable_cell_i["mean_RA"]))
        mean_DEC_list.append(bn.nanmean(subtable_cell_i["mean_DEC"]))
        healpix_phot_g_sum_flux_list.append(bn.nansum(subtable_cell_i["healpix_phot_g_sum_flux"]))
        healpix_phot_bp_sum_flux_list.append(bn.nansum(subtable_cell_i["healpix_phot_bp_sum_flux"]))
        healpix_phot_rp_sum_flux_list.append(bn.nansum(subtable_cell_i["healpix_phot_rp_sum_flux"]))
        sum_j_list.append(bn.nansum(subtable_cell_i["sum_j"]))
        sum_h_list.append(bn.nansum(subtable_cell_i["sum_h"]))
        sum_ks_list.append(bn.nansum(subtable_cell_i["sum_ks"]))
        sum_w1_list.append(bn.nansum(subtable_cell_i["sum_w1"]))
        sum_w2_list.append(bn.nansum(subtable_cell_i["sum_w2"]))
        sum_w3_list.append(bn.nansum(subtable_cell_i["sum_w3"]))
        sum_w4_list.append(bn.nansum(subtable_cell_i["sum_w4"]))



    healpix_table = pd.DataFrame({"healpix_lvl": healpix_lvl_list,
                                  "source_id": healpix_lvl_list,
                                  "source_id_2": "Healpix_superstar",
                                  "ra": np.float32(geo_RA_list),
                                  "dec": np.float32(geo_DEC_list),
                                  "wise_ra": np.float16(geo_RA_list),
                                  "wise_dec": np.float16(geo_DEC_list),
                                  "tmass_ra": np.float16(geo_RA_list),
                                  "tmass_dec": np.float16(geo_DEC_list),
                                  #"phot_bp_mean_mag": np.float16(-2.5*np.log10(healpix_phot_bp_sum_flux_list) + 8.9),
                                  #"phot_g_mean_mag":  np.float16(-2.5*np.log10(healpix_phot_g_sum_flux_list) + 8.9),
                                  #"phot_rp_mean_mag":  np.float16(-2.5*np.log10(healpix_phot_rp_sum_flux_list) + 8.9),
                                  #"parallax": np.zeros(len(healpix_lvl_list)),
                                  "w1": np.float16(sum_w1_list),
                                  #"w1_error": np.zeros(len(healpix_lvl_list)),
                                  "w2": np.float16(sum_w2_list),
                                  #"w2_error": np.zeros(len(healpix_lvl_list)),
                                  "w3": np.float16(sum_w3_list),
                                  #"w3_error": np.zeros(len(healpix_lvl_list)),
                                  "w4": np.float16(sum_w4_list),
                                  #"w4_error": np.zeros(len(healpix_lvl_list)),
                                  "j":  np.float16(sum_j_list),
                                  #"j_error": np.zeros(len(healpix_lvl_list)),
                                  "h":  np.float16(sum_h_list),
                                  #"h_error": np.zeros(len(healpix_lvl_list)),
                                  "ks":  np.float16(sum_ks_list),
                                  #"ks_error": np.zeros(len(healpix_lvl_list)),
                                  "g": np.float16(healpix_phot_g_sum_flux_list),
                                  "bp": np.float16(healpix_phot_bp_sum_flux_list),
                                  "rp": np.float16(healpix_phot_rp_sum_flux_list),
                                  "phot_g_mean_mag_AB": np.float16(-2.5*np.log10(healpix_phot_g_sum_flux_list) + 8.9),
                                  "phot_bp_mean_mag_AB":  np.float16(-2.5*np.log10(healpix_phot_bp_sum_flux_list) + 8.9),
                                  "phot_rp_mean_mag_AB":  np.float16(-2.5*np.log10(healpix_phot_rp_sum_flux_list) + 8.9),
                                  "phot_j_mean_mag_AB":  np.float16(-2.5*np.log10(sum_j_list) + 8.9),
                                  "phot_h_mean_mag_AB":  np.float16(-2.5*np.log10(sum_h_list) + 8.9),
                                  "phot_ks_mean_mag_AB":  np.float16(-2.5*np.log10(sum_ks_list) + 8.9),
                                  "phot_w1_mean_mag_AB": np.float16(-2.5*np.log10(sum_w1_list) + 8.9),
                                  "phot_w2_mean_mag_AB": np.float16(-2.5*np.log10(sum_w2_list) + 8.9),
                                  "phot_w3_mean_mag_AB": np.float16(-2.5*np.log10(sum_w3_list) + 8.9),
                                  "phot_w4_mean_mag_AB": np.float16(-2.5*np.log10(sum_w4_list) + 8.9)})

    healpix_table.to_csv("gaia_2mass_wise_healpix_" + str(healpix_lvl) +"_cat.csv")
    return(healpix_table)

############################
############################

def find_ra_dec_constraints(ra, dec, radius, verbose=False):

    coord = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    hp = astro_hp.HEALPix(nside=128, order='nested', frame=ICRS())
    healpix_ids_bounding_circle = hp.cone_search_skycoord(coord, radius=radius*5*u.deg)
    #print(healpix_ids_bounding_circle)

    healpix_radec_bounding_circle = hp.healpix_to_skycoord(healpix_ids_bounding_circle)
    #print(healpix_radec_bounding_circle)
    if verbose:
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
        ax.scatter(healpix_radec_bounding_circle.ra.deg, healpix_radec_bounding_circle.dec.deg)
        ax.set_xlim(0, 360)
        ax.set_ylim(-90,90)

    ra_min = bn.nanmin(healpix_radec_bounding_circle.ra.deg)
    ra_max = bn.nanmax(healpix_radec_bounding_circle.ra.deg)
    dec_min = bn.nanmin(healpix_radec_bounding_circle.dec.deg)
    dec_max = bn.nanmax(healpix_radec_bounding_circle.dec.deg)

    if verbose:
        print(ra_min)
        print(ra_max)
        print(dec_min)
        print(dec_max)

    # Is any of the poles contained in the bounding circle?
    from regions import CircleSkyRegion, make_example_dataset

    sky_radius = Angle(radius*2*u.deg)
    sky_region = CircleSkyRegion(coord, sky_radius)
    north_pole_coord = SkyCoord(1*u.deg, 90*u.deg, frame='icrs')
    south_pole_coord = SkyCoord(1*u.deg, -90*u.deg, frame='icrs')

    mock_dataset = make_example_dataset(data='simulated')
    wcs = mock_dataset.wcs
    bool_north_pole_in = sky_region.contains(north_pole_coord, wcs)
    if verbose:
        print("Bool North pole in:" + str(bool_north_pole_in))
    bool_south_pole_in = sky_region.contains(south_pole_coord, wcs)
    if verbose:
        print("Bool South pole in:" + str(bool_south_pole_in))

    # RA string
    if bool_north_pole_in or bool_south_pole_in:
        ra_sql_search_string = ""
        bounding_ra_min = 0
        bounding_ra_max = 360

    elif np.abs(ra_min - ra_max) > 180:
        # We are close to the edge of 0-360.
        ra_sql_search_string = "ra<=" + str(ra_min) + " AND ra>=" + str(ra_max)
        bounding_ra_min = ra_max
        bounding_ra_max = ra_min

    else:
        ra_sql_search_string = "ra>=" + str(ra_min) + " AND ra<= " + str(ra_max)
        bounding_ra_min = ra_min
        bounding_ra_max = ra_max

    # Dec string:
    if bool_north_pole_in:
        dec_sql_search_string = "dec>=" + str(dec_min)
        bounding_dec_min = dec_min
        bounding_dec_max = 90

    elif bool_south_pole_in:
        dec_sql_search_string = "dec<=" + str(dec_max)
        bounding_dec_min = -90
        bounding_dec_max = dec_max

    else:
        dec_sql_search_string = "dec>=" + str(dec_min) + " AND dec<=" + str(dec_max)
        bounding_dec_min = dec_min
        bounding_dec_max = dec_max

    sql_search_string = "WHERE " + ra_sql_search_string + " AND " + dec_sql_search_string + " "

    return({"sql_search_string": sql_search_string,
            "ra_min": ra_min, "ra_max": ra_max,
            "dec_min": dec_min, "dec_max": dec_max})
