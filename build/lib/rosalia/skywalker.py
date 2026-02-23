import matplotlib.pyplot as plt
import numpy as np
import pyvo as vo
from tqdm import tqdm
import pandas as pd
import os
import urllib.request
import rosalia as rs
#for service in vo.regsearch(datamodel="obscore"):
#  print(service['ivoid'])


def skywalker(instrument, extension, calib_level, collection):
    service = vo.dal.TAPService("https://mast.stsci.edu/vo-tap/api/v0.1/caom")

    unfiltered_table_list = []
    ra_list = np.linspace(0, 359, 360)
    for ra_i in tqdm(ra_list):
        #print(ra_i)
        ra_min = ra_i
        ra_max = ra_i + 1

        ex_query = """
        SELECT *
        FROM ivoa.obscore o
        WHERE
        obs_collection = '"""+ collection +"""'
        and dataproduct_type = 'image'
        and instrument_name = '""" + instrument + """'
        and access_format = 'application/fits'
        and calib_level = """ + str(calib_level) + """
        and o.s_ra >= """ + str(ra_min) + """
        and o.s_ra <  """ + str(ra_max) + """
        """
        # print(ex_query)
        result = service.search(ex_query)

        # print(len(result))
        if len(result) == 100000: raise Exception("Query truncated! Length >> 100000")

        unfiltered_table_list.append(result.to_table().to_pandas())

    unfiltered_table = pd.concat(unfiltered_table_list, axis=0)
    filtered_table = unfiltered_table[unfiltered_table['access_url'].str.contains(extension)]

    return(filtered_table)


def download_url(url):
    if not os.path.exists(url.split("/")[-1]):
        urllib.request.urlretrieve(url, url.split("/")[-1])
    else:
        print("File found, jumping!")
    return(url.split("/")[-1])

def download_query(csv_file, parallel=True):
    import pandas as pd
    import multiprocessing
    from tqdm import tqdm

    query_table = pd.read_csv(csv_file)

    if not parallel:
        for i in tqdm(range(len(query_table))):
            access_url = query_table["access_url"].iloc[i]
            try:
                downloaded_image = download_url(access_url)
                downloaded_images.append(downloaded_image)
            except:
                print("Oops. This file could not be downloaded.")
                print("Maybe it has access restricted or there is a query error.")
                print(access_url)
                print("Skipping")


    else:
        url_list = np.array(query_table["access_url"])
        nprocs = multiprocessing.cpu_count() - 2
        ############
        if multiprocessing.cpu_count() > 2:
            num_cores = multiprocessing.cpu_count() - 2
        else:
            num_cores = 1
        pool = multiprocessing.Pool(processes=num_cores)
        #############
        downloaded_images = []
        #     results = pool.starmap(my_function, tqdm.tqdm(inputs, total=len(param1)))
        for result in tqdm(pool.starmap(download_url, zip(url_list)), total=len(url_list)):
            downloaded_images.append(result)
        ##########
        pool.terminate()
    return(downloaded_images)

    
#####################################################################

def sort_hst_flcs_by_filter(filelist):
    #print(filelist)
    filters_list_keywords = rs.utils.get_keys_from_header(filelist, ["FILTER1", "FILTER2"], ext=0)
    filters = []
    for i in tqdm(range(len(filelist))):
        if not "CLEAR" in filters_list_keywords[0][i]:
            exposure_filter = filters_list_keywords[0][i]
        elif not "CLEAR" in filters_list_keywords[1][i]:
            exposure_filter = filters_list_keywords[1][i]
        else:
            exposure_filter = "None"
        filters.append(exposure_filter)

    list_of_filters = np.array(list(set(filters)))
    filters = np.array(filters)

    for filter_name in tqdm(list_of_filters):
        os.system("mkdir " + filter_name)
        exposures_with_that_filter = np.array(filelist)[np.where(filters == filter_name)[0]]

        if filter_name != "None":
            for selected_exposure in exposures_with_that_filter:
                os.system("mv " + selected_exposure + " " + filter_name)


    return(list_of_filters)
