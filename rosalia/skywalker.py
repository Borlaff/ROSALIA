# Alejandro S. Borlaff. NASA Ames Research Center. a.s.borlaff@nasa.gov / asborlaff@gmail.com
# January 20, 2023.
#
# ROSALIA/GNU module
# This module will hold all the general programs strictly derived from MAST
# including query wrappers, HST cleaning software, astrodrizle tasks.
#
# Version log:
# v.1.0 - 20 Enero 2023. First loading of programs inherited from former monolithic straycor.py
#
##########################################################

import os
import numpy as np
import pandas as pd
import pyvo as vo
from tqdm import tqdm
from astropy.io import fits
from astropy import coordinates
from astroquery.mast import Observations
import urllib.request
import rosalia as rs

#for service in vo.regsearch(datamodel="obscore"):
#  print(service['ivoid'])


# Importing modules

### Setting the HST calibration directories ###
#os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
#os.environ['CRDS_PATH'] = os.path.abspath(os.path.join('.', 'reference_files'))



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


############################

def download_mast(ra, dec, radius, filters, extension, instrument_name, project, obs_collection="HST"):
    # This program automatically downloads HST MAST images from the archive with positions as input

    c = coordinates.SkyCoord(ra, dec, frame='icrs', unit="deg")
    obs_table = Observations.query_criteria(coordinates=c, radius=radius, obs_collection=obs_collection,
                                            filters=filters, instrument_name=instrument_name)
    products_table = Observations.get_product_list(obs_table)
    filtered_products = Observations.filter_products(products_table, extension=extension, project=project)
    print(filtered_products)

    print("Downloading files...")
    output = Observations.download_products(filtered_products)
    output = np.array(list(set(np.array(output["Local Path"]))))

    print("Checking file integrity...")
    rs.utils.check_file_integrity(output)

    print("Downloading best CRDS reference files...")
    bestrefs(output, clean=True)

    print(output)

    absolute_path_output = []
    for i in tqdm(range(len(output))):
        absolute_path_output.append(os.path.abspath(output[i]))

    print("Done")
    return(absolute_path_output)

############################

def bestrefs(filename, clean=True):
    # If we input a list of files, then check them all
    if isinstance(filename, (list,pd.core.series.Series,np.ndarray)):
        filename_to_check = ''

        for i in range(len(filename)):
            filename_to_check = filename_to_check + " " + filename[i]

    else:
        filename_to_check = filename

    print(filename_to_check)

    cmd = 'crds bestrefs  --files ' + filename_to_check + ' --sync-references=1 --update-bestrefs'
    if not clean:
        print(cmd)
    # subprocess.run('crds bestrefs  --files ' + filename_to_check + ' --sync-references=1 --update-bestrefs', shell=True, capture_output=True,  flush=True).stdout
    rs.utils.execute_cmd(cmd)
    os.environ['iref'] = os.path.abspath(os.path.join('.', 'reference_files', 'references', 'hst', 'wfc3')) + os.path.sep
    os.environ['jref'] = os.path.abspath(os.path.join('.', 'reference_files', 'references', 'hst', 'acs')) + os.path.sep

    print("Bestrefs done.")

############################

############################

def clean_wcs(input_name, ext):
    # This program removes secondary WCS from Hubble images.

    keywords_to_be_cleaned = ["WCSNAMEA", "WCSAXESA", "CRPIX1A", "CRPIX2A",
                              "CUNIT1A" , "CUNIT2A", "CTYPE1A", "CTYPE2A",
                              "CRVAL1A", "CRVAL2A", "RADESYSA", "CD1_1A",
                              "CD1_2A", "CD2_1A", "CD2_2A", "WCSTYPEA",
                              "TDD_CTA", "TDD_CYA", "TDD_CXA"]



    for keyword in keywords_to_be_cleaned:
        for wcs_extension in ["A", "B", "C", "D", "E", "F", "O"]:
            keyword = keyword[:-1] + wcs_extension
            #print(keyword)
            modify_keyword_hdr(input_name=input_name, ext=ext, keyword=keyword, mode="delete", verbose=False)
    return(input_name)



def ensure_list(input_list_or_string):
    # This function checks if the input is a string or a list, and always return a list.

    # List of input images: Checking that it is a list.
    if isinstance(input_list_or_string, (list, pd.core.series.Series, np.ndarray)):
        #print("Executing for images:")
        input_list = input_list_or_string
        print(input_list)

    else:
        #print("Executing for single image " + input_list_or_string)
        input_list = [input_list_or_string]

    return(input_list)


def heal_ACS_DQ_type(input_list):
    # Let's ensure that this is a list
    input_list = ensure_list(input_list)

    # For ACS, there is a weird bug that modifies the type of the first DQ extension [ext 3]
    # from int16 to uint16. That causes the files to be unusable for future astrodrizzle runs.
    # This program ensures that the DQ have the right type.

    for input_name in tqdm(input_list):
        input_fits = fits.open(input_name)
        input_fits[3].data = input_fits[3].data.astype(np.int16)
        input_fits.verify("silentfix")
        input_fits.writeto(input_name, overwrite=True)

    return(input_list)
