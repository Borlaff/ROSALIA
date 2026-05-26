import time
import numpy as np
import pandas as pd
import ephessos as ep
from tqdm import tqdm
import astropy.units as u
from astropy.coordinates import SkyCoord
from rosalia.utils import divide_array_in_chunks

def get_jpl_observer_name(observer_name):
    obs_center = ep.core.jpl_name_translation(observer_name)
    return(obs_center)



def get_mpc_observer_name(observer_name):
    from astroquery.mpc import MPC

    if observer_name == "HST" or observer_name == "Hubble": return("250")
    if observer_name == "RST" or observer_name == "Roman": return("289")
    if observer_name == "Euclid": return("273")

    try:
        obs = MPC.get_observatory_codes()
        return(obs[obs["Name"] == observer_name]["Code"][0])
    except:
        print("MPC Code not found!")
    


# ROLO lunar model
# https://iopscience.iop.org/article/10.1086/430185/pdf
def horizons_query(mjd, source="301", location="@jwst", chunk_size=50):
    import pandas as pd
    if not isinstance(mjd, (list, pd.core.series.Series, np.ndarray)):
        return(single_horizons_query(mjd, source=source, location=location))

    else:
        mjd_chunks = divide_array_in_chunks(mjd, chunk_size = chunk_size)

        RA_all = []
        DEC_all = []
        V_all = []
        ang_width_all = []

        datetime_jd_all = []
        for mjd_chunk_i in mjd_chunks:
            target_query_chunk_i = single_horizons_query(mjd_chunk_i, source=source, location=location)
            #print(target_query_chunk_i)
            RA_all = RA_all + list(target_query_chunk_i["RA"])
            DEC_all = DEC_all + list(target_query_chunk_i["DEC"])
            V_all = V_all + list(target_query_chunk_i["V"])
            ang_width_all = ang_width_all + list(target_query_chunk_i["ang_width"])

            datetime_jd_all = datetime_jd_all + list(target_query_chunk_i["datetime_jd"])

        return({"RA": np.array(RA_all), "DEC": np.array(DEC_all),
            "V": np.array(V_all),  "ang_width": np.array(ang_width_all), "datetime_jd": np.array(datetime_jd_all)})


def single_horizons_query(mjd, source="301", location="@hst"):
    # ----------------------------- #
    # ROLO lunar model
    # https://iopscience.iop.org/article/10.1086/430185/pdf
    # Moon id - 301
    # From HST - location = @hst
    # ----------------------------- #
    from astroquery.jplhorizons import Horizons

    jpl_horizons_query = Horizons(id=source, location=location, epochs=np.array(mjd) + 2400000.5)
    #print("Epochs:")
    #print(np.array(mjd) + 2400000.5)

    tolerance_time = 1/60./60.
    residual_time = tolerance_time + 1

    while np.abs(np.nanmax(residual_time))>tolerance_time: # If there is a mismatch higher than 1 sec, repeat.
        try:
            ephemeris_source_from_location = jpl_horizons_query.ephemerides()
        except Exception as error:
            # Handle the exception
            print("An exception occurred:", error) # An exception occurred: division by zero
            dummy = np.array([np.nan]*len(mjd))
            return({"RA": dummy, "DEC": dummy, "V": dummy, "datetime_jd": dummy})

        # RA_moon[i] = bn.nanmedian(moon_from_hst["RA"])
        # DEC_moon[i] = bn.nanmedian(moon_from_hst["DEC"])
        #print(ephemeris_source_from_location)
        RA_source = ephemeris_source_from_location["RA"].value.data
        DEC_source = ephemeris_source_from_location["DEC"].value.data
        V_source = ephemeris_source_from_location["V"].value.data
        datetime_jd = ephemeris_source_from_location["datetime_jd"].value.data
        residual_time = datetime_jd - np.array(mjd) - 2400000.5
        ang_width = ephemeris_source_from_location["ang_width"]/60/60

        #print(residual_time)

        if np.abs(np.nanmax(residual_time))>tolerance_time:
            print("\nPotential mismatch in input/output MJD.")
            print("Epoch intro = " + str(np.array(mjd) + 2400000.5))
            print("Epoch output = " + str(np.array(datetime_jd)))
            print("Delta time = " + str(residual_time))
            print("Sleeping for 5 seconds...")
            time.sleep(5)

    return({"RA": RA_source, "DEC": DEC_source,
            "V": V_source, "ang_width": ang_width,
            "datetime_jd": datetime_jd,
            "jpl_horizons_query": jpl_horizons_query})



def find_moon_and_jupiter_in_HST_history(obs_history, chunksize=50):

    # Order the input by t_min if it wasn't

    n_obs = len(obs_history)
    RA_moon = np.empty(len(obs_history)) + np.nan
    DEC_moon = np.empty(len(obs_history)) + np.nan
    obs_history["t_med"] = np.array((obs_history["t_max"] + obs_history["t_min"])/2.)
    t_med = np.copy(np.array(obs_history["t_med"]))
    obs_history = obs_history.sort_values(by=['t_med'])

    RA_moon = np.empty(n_obs) * np.nan
    DEC_moon = np.empty(n_obs) * np.nan
    V_moon = np.empty(n_obs) * np.nan
    datetime_jd_moon = np.empty(n_obs) * np.nan
    diff_datetime_moon = np.empty(n_obs) * np.nan

    RA_jupiter = np.empty(n_obs) * np.nan
    DEC_jupiter = np.empty(n_obs) * np.nan
    V_jupiter = np.empty(n_obs) * np.nan
    datetime_jd_jupiter = np.empty(n_obs) * np.nan
    diff_datetime_jupiter = np.empty(n_obs) * np.nan

    counter = 0

    pbar = tqdm(total=n_obs, position=0, leave=True)
    while counter < n_obs:

        min_i = counter
        max_i = counter+chunksize
        obs_history_chunk = obs_history.iloc[min_i:max_i]#.sort_values(by=['t_med'])
        #print(len(obs_history_chunk))
        #print(min_i)
        #print(max_i)
        t_med_chunk = obs_history_chunk["t_med"]
        #print("Len tmedchunck")
        #print(len(t_med_chunk))
        #print(t_med_chunk)

        RA_DEC_moon_chunck = horizons_query(t_med_chunk, source="301", location="@hst")
        RA_DEC_moon_chunck = pd.DataFrame(RA_DEC_moon_chunck)
        RA_DEC_moon_chunck = RA_DEC_moon_chunck.sort_values(by=['datetime_jd'])

        RA_DEC_jupiter_chunck = horizons_query(t_med_chunk, source="599", location="@hst")
        RA_DEC_jupiter_chunck = pd.DataFrame(RA_DEC_jupiter_chunck)
        RA_DEC_jupiter_chunck = RA_DEC_jupiter_chunck.sort_values(by=['datetime_jd'])
        #print(RA_DEC_moon_chunck)

        RA_moon[min_i:max_i]       = RA_DEC_moon_chunck["RA"]
        DEC_moon[min_i:max_i]      = RA_DEC_moon_chunck["DEC"]
        V_moon[min_i:max_i]        = RA_DEC_moon_chunck["V"]
        datetime_jd_moon[min_i:max_i]   = np.array(RA_DEC_moon_chunck["datetime_jd"])
        diff_datetime_moon[min_i:max_i] = np.array(RA_DEC_moon_chunck["datetime_jd"]) - np.array(obs_history_chunk["t_med"]) - 2400000.5


        RA_jupiter[min_i:max_i]       = RA_DEC_jupiter_chunck["RA"]
        DEC_jupiter[min_i:max_i]      = RA_DEC_jupiter_chunck["DEC"]
        V_jupiter[min_i:max_i]        = RA_DEC_jupiter_chunck["V"]
        datetime_jd_jupiter[min_i:max_i]   = np.array(RA_DEC_jupiter_chunck["datetime_jd"])
        diff_datetime_jupiter[min_i:max_i] = np.array(RA_DEC_jupiter_chunck["datetime_jd"]) - np.array(obs_history_chunk["t_med"]) - 2400000.5

        # House-keeping procedures
        counter = counter + chunksize
        pbar.update(chunksize)

    obs_history["RA_moon"] = RA_moon
    obs_history["DEC_moon"] = DEC_moon
    obs_history["V_moon"] = V_moon
    obs_history["t_med"] = t_med
    obs_history["datetime_jd_moon"] = datetime_jd_moon
    obs_history["diff_datetime_moon"] = diff_datetime_moon

    obs_history["RA_jupiter"] = RA_jupiter
    obs_history["DEC_jupiter"] = DEC_jupiter
    obs_history["V_jupiter"] = V_jupiter
    obs_history["t_med"] = t_med
    obs_history["datetime_jd_jupiter"] = datetime_jd_jupiter
    obs_history["diff_datetime_jupiter"] = diff_datetime_jupiter

    hst_pointing = SkyCoord(np.array(obs_history["s_ra"])*u.deg, np.array(obs_history["s_dec"])*u.deg, frame='icrs')
    moon = SkyCoord(np.array(obs_history["RA_moon"])*u.deg, np.array(obs_history["DEC_moon"])*u.deg, frame='icrs')
    jupiter = SkyCoord(np.array(obs_history["RA_jupiter"])*u.deg, np.array(obs_history["DEC_jupiter"])*u.deg, frame='icrs')
    sep_moon = hst_pointing.separation(moon)
    sep_jupiter = hst_pointing.separation(jupiter)
    obs_history["Moon_sep"] = sep_moon.value
    obs_history["Jupiter_sep"] = sep_jupiter.value

    return(obs_history)
