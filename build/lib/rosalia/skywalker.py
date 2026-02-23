import matplotlib.pyplot as plt

import pyvo as vo
from tqdm import tqdm
import pandas as pd

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
        print(ex_query)
        result = service.search(ex_query)

        print(len(result))
        if len(result) == 100000: raise Exception("Query truncated! Length >> 100000")

        unfiltered_table_list.append(result.to_table().to_pandas())

    unfiltered_table = pd.concat(unfiltered_table_list, axis=0)
    filtered_table = unfiltered_table[unfiltered_table['access_url'].str.contains(extension)]

    return(filtered_table)

def download_query(csv_file):
    import urllib.request
    import pandas as pd
    query_table = pd.read_csv(csv_file)

    for i in tqdm(range(len(query_table))):
        try:
            exposure_name = query_table["access_url"].iloc[i]
            urllib.request.urlretrieve(exposure_name, exposure_name.split("/")[-1])

        except:
            print("Oops. This file could not be downloaded.")
            print("Maybe it has access restricted or there is a query error.")
            print(exposure_name)
            print("Skipping")
