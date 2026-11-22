import os
import numpy as np
from astropy.time import Time
import pandas as pd
from astroquery.jplhorizons import Horizons
import time
from tqdm import tqdm 
from pathlib import Path
from rosalia.horizons import get_jpl_observer_name
from rosalia.utils import save_dict, load_dict
from scipy.interpolate import interp1d
from scipy import interpolate

# From https://ssd-api.jpl.nasa.gov/doc/horizons.html#encode
# To successfully transmit such a list via URL to the API,
# the spaces (or commas) and single quotes in the TLIST setting should be URL-encoded using
# %20 for spaces,
# %2C for commas,
# %27 for single-quotes.
# The resulting URL for successful API handling of the example would look like this:
# TLIST=%272455339.95748%27%20%272455354.92142%27%20%272033-Jan-17%2012:10:25.1%27
def sso_interpolators(observer="roman", start_epoch='2026-08-31', end_epoch=None, step='1h'):
    """
      _load_sso_interpolators: Generates a set of interpolators for the Solar System Objects from the POV of your observatory. 
    
      input: 
      location: JPL/Horizons observatory code. Defaults to 500@-680 (Roman Space Telescope)
      start_epoch: Start of epheremis. Defaults to August 31st, 2026, launch date for Roman Space Telescope
      end_epoch: End of ephemeris. Defaults to today + 1 day. 
      step: Time resolution of the ephemeris: Default to 1h.  
    """
    ori_observer = observer
    observer = get_jpl_observer_name(observer)
    observer = "500@" + observer

    file_path = Path(os.environ["ROSALIACACHE"] + observer + "_ephem_interpolator.pkl")

    # Check if the file exists and was modified more than 24 hours ago
    if file_path.exists():
        file_age = time.time() - file_path.stat().st_mtime
        if file_age > 86400: # one day in seconds 
            print(f"The ephemeris is older than 1 day! Reloading...")
        else:
            return(load_dict(file_path))

    SSO_list = ["10", "301", "399", "199", "299", "499", "599", "699", "799", "899", "999"]
    SSO_list_plain = ["Sun", "Moon", "Earth", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]

    if end_epoch is None: end_epoch = Time.now().isot

    ephemeris_interpolators = {} 
    
    for i in tqdm(range(len(SSO_list)), desc="Loading " + ori_observer + " ephemeris: "):
        sso = SSO_list[i] 
        sso_name = SSO_list_plain[i]
        obj = Horizons(id=sso, location=observer, epochs={'start': start_epoch, 'stop': end_epoch, 'step': step})
        eph = obj.ephemerides()
        mjd_list = Time(eph["datetime_jd"], format="jd").mjd
        ra_interp = interp1d(x=mjd_list, y=eph["RA"], kind="quadratic", fill_value="extrapolate")
        dec_interp = interp1d(x=mjd_list, y=eph["DEC"], kind="quadratic", fill_value="extrapolate")
        V_interp = interp1d(x=mjd_list, y=eph["V"], kind="quadratic", fill_value="extrapolate")
        ang_width_interp = interp1d(x=mjd_list, y=eph["ang_width"], kind="quadratic", fill_value="extrapolate")

        ephemeris_interpolators[sso_name] = {"ra_interp": ra_interp, 
                                             "dec_interp": dec_interp,
                                             "V_interp": V_interp, 
                                              "ang_width_interp": ang_width_interp}

    save_dict(ephemeris_interpolators, file_path)
    return(ephemeris_interpolators)


def get_ephemeris(observer, object, MJD):
    # If observer is a common keyword, transform
    observer = get_jpl_observer_name(observer)

    astropy_time = Time(MJD, format="mjd")

    obj = Horizons(id=object,
                   location=observer,
                   epochs=astropy_time.jd)

    eph = obj.ephemerides()
    return(eph)


def get_SSOs_ephemeris(observer, MJD):
    # If observer is a common keyword, transform
    observer = get_jpl_observer_name(observer)

    #          0 -> 1    SOLAR SYSTEM BARYCENTER -> MERCURY BARYCENTER
    #          0 -> 2    SOLAR SYSTEM BARYCENTER -> VENUS BARYCENTER
    #          0 -> 3    SOLAR SYSTEM BARYCENTER -> EARTH BARYCENTER
    #          0 -> 4    SOLAR SYSTEM BARYCENTER -> MARS BARYCENTER
    #          0 -> 5    SOLAR SYSTEM BARYCENTER -> JUPITER BARYCENTER
    #          0 -> 6    SOLAR SYSTEM BARYCENTER -> SATURN BARYCENTER
    #          0 -> 7    SOLAR SYSTEM BARYCENTER -> URANUS BARYCENTER
    #          0 -> 8    SOLAR SYSTEM BARYCENTER -> NEPTUNE BARYCENTER
    #          0 -> 9    SOLAR SYSTEM BARYCENTER -> PLUTO BARYCENTER
    #          0 -> 10   SOLAR SYSTEM BARYCENTER -> SUN
    #          3 -> 301  EARTH BARYCENTER -> MOON
    #          3 -> 399  EARTH BARYCENTER -> EARTH
    #          1 -> 199  MERCURY BARYCENTER -> MERCURY
    #          2 -> 299  VENUS BARYCENTER -> VENUS
    #          4 -> 499  MARS BARYCENTER -> MARS
    #          5 -> 599  Jupiter BARYCENTER -> Jupiter
    #          6 -> 699  Saturn BARYCENTER -> Saturn
    #          7 -> 799  Uranus BARYCENTER -> Uranus
    #          8 -> 899  Neptune BARYCENTER -> Neptune
    #          9 -> 999  Pluto BARYCENTER -> Pluto
    #.         10        Sun
    SSO_list = ["10", "301", "399", "199", "299", "499", "599", "699", "799", "899", "999"]
    SSO_list_plain = ["Sun", "Moon", "Earth", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
    SSO_ephemeris_list = []
    for object in SSO_list:
        SSO_ephemeris = get_ephemeris(observer, object, MJD).to_pandas().dropna(axis=1)
        SSO_ephemeris_list.append(SSO_ephemeris)

    SSO_ephemeris_db = pd.concat(SSO_ephemeris_list)

    # https://oro.open.ac.uk/54834/1/593910.pdf
    # Secondary correction to the V-band from Hipparcos.
    # The Johnson UBVRI system Horizons (and asteroid work in general) uses is
    # defined so that alpha-Lyrae (Vega) has V= 0.03, with U-B = B-V = 0.
    #
    # Correcting the magnitude by -0.03 give us Vega magnitudes.
    # The difference between Johnson V in Vega and AB is minimal, Vega - AB = -0.01.
    # So AB mag of the Sun is = V_Hozirons - 0.02
    # --------------------------------------------------------------------------#

    SSO_ephemeris_db["V_ab"] = SSO_ephemeris_db["V"] - 0.03 + 0.01
    SSO_ephemeris_db["source_id"] = SSO_list_plain
    SSO_ephemeris_db["source_id_2"] = np.array(SSO_list, dtype="int")

    return(SSO_ephemeris_db)


def get_SS0s_loc_magnitude(observer, MJD, lambda_ref, verbose=False):
    # JPL/Horizons does not support queries in the future.
    # If the exposure time is > than current, then take the "now" value as proxy.
    # This has to be upgraded in the future.

    observer = get_jpl_observer_name(observer)
    # Make the Horizons query
    SSO_ephemeris_db = get_SSOs_ephemeris(observer, MJD)
    #ephem = sso_interpolators(observer=observer)
    
    # Load the spectra of the Sun
    mag_sun_db = pd.read_csv(os.environ["ROSALIACACHE"] + "/CORE/magnitude_sun.txt", sep="\t")

    if verbose: print("Wavelengths to estimate apparent Magnitude of the Sun in UV / Optical / IR:")
    if verbose: print(mag_sun_db["lambda_mum"]*10**4)
    if verbose: print("Apparent Magnitude of the Sun in UV / Optical / IR:")
    if verbose: print(mag_sun_db["AB_App"])

    magnitude_interpolator = interpolate.interp1d(x=mag_sun_db["lambda_mum"]*10**4, y=mag_sun_db["AB_App"])
    mag_ref = magnitude_interpolator(lambda_ref.to("AA"))
    if verbose: print("Magnitude of the Sun in the selected filter - mag_ref:")
    if verbose: print(mag_ref)

    # We need to estimate the magnitude of the solar system bodies in different bands,
    # but the JPL / Horizons only provides V band magnitudes.
    # To solve the problem, the estimate a color factor for the magnitude of the Sun,
    # and then compensate for the color term in the magnitude of the rest of the bodies.
    sun_color = mag_ref - SSO_ephemeris_db[SSO_ephemeris_db["targetname"]=="Sun (10)"]["V_ab"]
    SSO_ephemeris_db["mag_custom"] = sun_color + SSO_ephemeris_db["V_ab"]
    return(SSO_ephemeris_db)
