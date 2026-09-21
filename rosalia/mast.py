import os 
import roman_datamodels as rdm
from astroquery.mast import MastMissions # pip install git+https://github.com/snbianco/astroquery.git@mast-beta-release 


def stream_roman_mast(products, row=0):
    missions = MastMissions(mission='roman')
    token = os.getenv("MAST_API_TOKEN")    
    missions.login(token=token)
    af = missions.read_product(products[row]['filename'])
    dm = rdm.open(af)
    return(dm)

####################


def pack_exposures(products):
    # if True:
    # print(products["filename"])
    root_id = []
    for i in range(len(products["filename"])):
        split_filename = products["filename"][i].split("_")
        root_id.append(split_filename[0] + "_" + split_filename[1])
    unique_roots = list(set(root_id))
    products["exposure_id"] = root_id
    return(products, unique_roots)

def roman_query(criteria=None, coordinates=None, radius=None, filter=None, detector=None, file_suffix=None):

    """
    # Create a dictionary of search criteria
    # Example APT 1020, observation 5, exposure 4 pass 2 SCA 6 (F158)
    search = {'program': 1020,
            'observation': 5,
            'pass': 2, 
            #'exposure_start_time': '2026-09-15T02:01:47.4080000',
            'product_type': 'l2', 
            'detector': 'WFI06',
            'optical_element': 'F158',
            'exposure_type': "WFI_IMAGE"}
    """

    # Create MastMissions object and assign mission to 'roman'
    missions = MastMissions(mission='roman')

    # Login to search and retrieve Roman data
    token = os.getenv("MAST_API_TOKEN")    
    missions.login(token=token)
                
    print(f'Mission: {missions.mission}')
    print(f'Service: {missions.service}')

    #---------- Connection established ----------------# 

    # Make a list of column names to return in the search results. The results will not be ordered
    # in this order, so we will re-order later. The fileSetName, which we need for retrieval,
    # is not in this list, but will still be returned.
    col_list = ['program', 'execution_plan', 'pass', 'segment', 'visit', 'observation',
                'optical_element', 'exposure_type', 'instrument_name', 'detector', 'productLevel', 
                'product_type', 'exposure_time', 'exposure_start_time', 'exposure_end_time', 'fileSetName']

    if criteria is not None:
        # Query with column criteria
        results = missions.query_criteria(**criteria, select_cols=col_list)

    if coordinates is not None:
        results = missions.query_region(coordinates=coordinates, radius=radius)

    if filter is not None: results = results[results["optical_element"] == filter]
    if detector is not None: results = results[results["detector"] == detector]
    
    # Re-order the column names 
    # results = results[col_list]
    products = missions.get_unique_product_list(results)

    if file_suffix is not None:
        products = missions.filter_products(products, file_suffix='_cal')

    products, unique_roots = pack_exposures(products)

    return(results, products, unique_roots)

