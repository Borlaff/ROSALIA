# Alejandro S. Borlaff. NASA Ames Research Center. a.s.borlaff@nasa.gov / asborlaff@gmail.com
# January 20, 2023.
#
# STRAYCOR/IRSA module
# This module will hold all the general programs strictly derived from IRSA
# including query wrappers.
#
# Version log:
# v.1.0 - 26 Enero 2023. First loading of programs inherited from former
#                        Euclid project - irsa_query.py
#
##########################################################
import os
import pandas as pd
import numpy as np
import multiprocess
from scipy import interpolate
from tqdm import tqdm
from astropy.time import Time
import rosalia as rs

#########################

def launch_irsa_query(ra, dec, wavelength, year, day, obslocin=3):

    # Keywords for IRSA background query:
    # Check: https://irsa.ipac.caltech.edu/applications/BackgroundModel/docs/dustProgramInterface.html
    #
    #
    # locstr 	(locstr examples) 	NA 	If the input is a recognizable object name, it will be resolved into coordinates using NED or, if that fails, SIMBAD (required).
    # wavelength 	float 	0.5-1000.0 	Wavelength in microns (defaults to 2.0).
    # year 	char 	varies 	Year. Limited to 2018 to 2029 for L2 position. Defaults to 2019.
    # day 	char 	1-366 	Day. Limited to 2018 Day 274 to 2029 Day 120 for L2 position and ido_viewin=0. Defaults to 180.
    # obslocin 	char 	0 or 3 	Observing location. 0 is Earth-Sun L2 region; 3 is Earth (defaults to 0).
    # obsverin 	char 	1 or 4 	Code version (optional, defaults to 4).
    # ido_viewin 	char 	0 or 1 	0 = find zodiacal on Day; 1 = find median zodiacal over a likely viewing range (defaults to 1, see Help).


    ra_str = str(round(ra,4)).zfill(4)
    dec_str = str(round(dec,4)).zfill(4)
    wave_str = str(round(wavelength,2)).zfill(2)
    year_str = str(year)
    day_str = str(day).zfill(2)

    outfile="irsa_ra"+ra_str+"_dec"+dec_str+"_wave"+wave_str+"_y"+year_str+"_d"+day_str+".xml"

    # rs.utils.execute_cmd("rm " + outfile)

    website="https://irsa.ipac.caltech.edu/cgi-bin/BackgroundModel/nph-bgmodel?"
    cmd = 'curl -o ' + outfile + ' "' + website + 'locstr='+str(ra)+'+'+str(dec)+'+equ+j2000&wavelength='+str(wavelength)+'&year='+str(year)+'&day='+str(day)+'&obslocin='+str(obslocin)+'&ido_viewin=0"'

    while not (os.path.exists(outfile) and os.path.getsize(outfile)) > 0:
        rs.utils.execute_cmd(cmd)

    return(outfile)


def read_irsa_query(xml_query_name):
    import xmltodict
    #print("Reading " + infile)
    xml_data = open(xml_query_name, 'r').read()  # Read data
    xmlDict = xmltodict.parse(xml_data)  # Parse XML

    d = {}
    coords_string = xmlDict['results']["input"]["objname"].split(" ")
    d['ra'] = float(coords_string[0])
    d['dec'] = float(coords_string[1])
    d['wavelength'] = float(xmlDict['results']["input"]["wavelength"].replace("microns",""))
    d['year'] = float(xmlDict['results']["input"]["year"])
    d['day'] = float(xmlDict['results']["input"]["day"])
    d['zody'] = rs.utils.MJysr_to_jyarcsec2(float(xmlDict['results']['result']['statistics']['zody'].replace("(MJy/sr)","")))
    d['cib'] = rs.utils.MJysr_to_jyarcsec2(float(xmlDict['results']['result']['statistics']['cib'].replace("(MJy/sr)","")))
    d['stars'] = rs.utils.MJysr_to_jyarcsec2(float(xmlDict['results']['result']['statistics']['stars'].replace("(MJy/sr)","")))
    d['ism'] = rs.utils.MJysr_to_jyarcsec2(float(xmlDict['results']['result']['statistics']['ism'].replace("(MJy/sr)","")))
    d['totbg'] = rs.utils.MJysr_to_jyarcsec2(float(xmlDict['results']['result']['statistics']['totbg'].replace("(MJy/sr)","")))

    os.system("rm " + xml_query_name)

    return(d)

def irsa_query_single(ra, dec, wavelength, year, day, obslocin):
    ra_str = str(round(ra,4)).zfill(4)
    dec_str = str(round(dec,4)).zfill(4)
    wave_str = str(round(wavelength,2)).zfill(2)
    year_str = str(year)
    day_str = str(day).zfill(2)

    outfile="irsa_ra"+ra_str+"_dec"+dec_str+"_wave"+wave_str+"_y"+year_str+"_d"+day_str+".xml"
    irsa_query_file = launch_irsa_query(ra, dec, wavelength, year, day, obslocin)
    irsa_query_result = read_irsa_query(irsa_query_file)
    return(irsa_query_result)


def irsa_query(ra, dec, wavelength, year, day, obslocin):
    # Check that ra and dec arguments are in array
    #print("IRSA Query")
    #print(ra)
    #print(dec)
    #print(wavelength)
    #print(year)
    #print(day)
    #print(obslocin)

    if not isinstance(ra, (list, pd.core.series.Series, np.ndarray)):
        ra = np.array([ra])
    if not isinstance(dec, (list, pd.core.series.Series, np.ndarray)):
        dec = np.array([dec])

    # If wavelength, year, or day are not in array, copy their values into one
    # as large as ra, dec
    if not isinstance(wavelength, (list, pd.core.series.Series, np.ndarray)):
        wavelength = np.array([wavelength]*len(ra))

    if not isinstance(year, (list, pd.core.series.Series, np.ndarray)):
        year = np.array([year]*len(ra))

    if not isinstance(day, (list, pd.core.series.Series, np.ndarray)):
        day = np.array([day]*len(ra))

    if not isinstance(obslocin, (list, pd.core.series.Series, np.ndarray)):
        obslocin = np.array([obslocin]*len(ra))

    #print("Launching parallel IRSA query:")
    #print("wavelength: ")
    #print(wavelength)
    arguments = zip(ra, dec, wavelength, year, day, obslocin)
    num_cores = np.max(np.array([multiprocess.cpu_count(), 1]))
    pool = multiprocess.Pool(processes=num_cores)
    query = pool.starmap(irsa_query_single, tqdm(arguments, total=len(ra)))
    pool.terminate()
    pool.close()

    return(pd.DataFrame(query))
