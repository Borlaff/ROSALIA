import requests
import pandas as pd
import numpy as np
from tqdm import tqdm

def load_skysurf():
    acs_filters = ["F435W", "F475W", "F555W", "F606W",
                "F625W", "F775W", "F814W", "F850LP"]

    wfc3ir_filters = ["F098M", "F105W", "F110W", "F125W",
                   "F140W", "F160W"]

    wfc3uvis_filters = ["F225W", "F275W", "F300X", "F336W",
                    "F390W", "F438W", "F475W", "F555W",
                    "F606W", "F625W", "F775W", "F850LP",
                    "F814W"]

    skysurf_db = {}
    print("Loading ACS")
    instrument = "acswfc"
    skysurf_db[instrument] = {}
    for filter_hst in tqdm(acs_filters):
        try:
            skysurf_db[instrument][filter_hst] = load_skysurf_filter(instrument, filter_hst)
        except:
            print("WARNING: " +  instrument + " - " + filter_hst + ": Not found")

    print("Loading WFC3/IR")
    instrument = "wfc3ir"
    skysurf_db[instrument] = {}
    for filter_hst in tqdm(wfc3ir_filters):
        try:
            skysurf_db[instrument][filter_hst] = load_skysurf_filter(instrument, filter_hst)
        except:
            print("WARNING: " +  instrument + " - " + filter_hst + ": Not found")

    print("Loading WFC3/UVIS")
    instrument = "wfc3uvis"
    skysurf_db[instrument] = {}
    for filter_hst in tqdm(wfc3uvis_filters):
        try:
            skysurf_db[instrument][filter_hst] = load_skysurf_filter(instrument, filter_hst)
        except:
            print("WARNING: " +  instrument + " - " + filter_hst + ": Not found")
    return(skysurf_db)


def load_skysurf_filter(instrument, filter_hst):
    # Load the SKYSURF database
    # http://skysurf.asu.edu/sky_measurements/SKYSURF_profoundmed_acswfc_F435W.csv
    skysurf_server = "http://skysurf.asu.edu/sky_measurements/"

    csv_filename = "SKYSURF_profoundmed_"+instrument+"_"+filter_hst+".csv"
    url_table = "http://skysurf.asu.edu/sky_measurements/"+csv_filename
    #print(url_table)
    #skysurf_db = pd.read_html(requests.get(url_table).content)[-1].to_csv(csv_filename)
    r = requests.get(url_table) # create HTTP response object

    with open(csv_filename,'wb') as f:
        f.write(r.content)

    skysurf_db = pd.read_csv(csv_filename)

    return(skysurf_db)


def find_obsid_in_skysurf(obsid):

    # Find the background of a certain image in the skysurf db
    skysurf_instruments = ["acswfc", "wfc3ir", "wfc3uvis"]
    skysurf_db = load_skysurf()
    matched_filename_list = []
    matched_instrument_list = []
    matched_filter_list = []
    obsid_list = []

    if isinstance(obsid, (list,pd.core.series.Series,np.ndarray)):
        print("Obsids to check:" + str(obsid))
    else:
        obsid = [obsid]
        print("Obsids to check:" + str(obsid))


    for instrument in skysurf_instruments:
        instrument_filters = skysurf_db[instrument].keys()

        for instrument_filter in instrument_filters:
            filename_list = np.array(skysurf_db[instrument][instrument_filter]["File"])
            for i in range(len(filename_list)):

                for obsid_i in obsid:
                    if obsid_i in filename_list[i]:
                        #print(skysurf_db[instrument][instrument_filter].iloc[i])
                        matched_filename_list.append(i)
                        matched_instrument_list.append(instrument)
                        matched_filter_list.append(instrument_filter)
                        obsid_list.append(obsid_i)

    #print(skysurf_db[matched_instrument_list][matched_filter_list].iloc[matched_filename_list])
    #print(matched_instrument_list)
    #print(matched_filter_list)
    #print(matched_filename_list)

    if len(matched_filename_list) > 0:
        row_list = []
        for i in range(len(matched_filename_list)):
            row_list.append(skysurf_db[matched_instrument_list[i]][matched_filter_list[i]].iloc[matched_filename_list[i]].to_dict())

        matched_files_db = pd.DataFrame(row_list)
        matched_files_db["instrument"] = matched_instrument_list
        matched_files_db["filter"] = matched_filter_list
        matched_files_db["obsid"] = obsid_list
        return(matched_files_db)
