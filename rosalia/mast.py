import os 
import roman_datamodels as rdm
from astroquery.mast import MastMissions # pip install git+https://github.com/snbianco/astroquery.git@mast-beta-release 
import xml.etree.ElementTree as ET
import astropy.units as u
import logging
from astropy.time import Time
import rosalia as rs
from tqdm import tqdm 

# Supress info from MAST
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)

class APT():
    def update_target_coordinates(xml_path, output_path, new_coords):
        """
        Updates RA/Dec coordinates for each FixedTarget in the XML.

        Parameters:
            xml_path (str): Path to the input XML file.
            output_path (str): Path to write the updated XML file.
            new_coords (list of tuples): List of (RA, Dec) strings.
                Example: [("13 47 3.30", "+49 01 40.49"), ...]
                Must be same length as number of FixedTarget entries.
        """

        tree = ET.parse(xml_path)
        root = tree.getroot()

        # XML uses namespaces; must extract them for searching.
        ns = {'ns': root.tag.split('}')[0].strip('{')}

        # Find all FixedTarget entries
        targets = root.findall('.//ns:FixedTarget', ns)

        if len(new_coords) != len(targets):
            raise ValueError(
                f"Provided {len(new_coords)} coordinates but XML contains {len(targets)} FixedTargets."
            )

        for (target, coord) in zip(targets, new_coords):
            ra = coord.ra.to_string(unit=u.hour, sep=' ', precision=4, pad=True)
            dec = coord.dec.to_string(unit=u.degree, sep=' ', precision=2, pad=True, alwayssign=True)
            # Find EquatorialCoordinates node
            eq = target.find('ns:EquatorialCoordinates', ns)
            if eq is not None:
                # Format must match XML's "Value" attribute: "RA Dec"
                eq.set("Value", f"{ra} {dec}")
                # print("Value", f"{ra} {dec}")
            else:
                print(f"Warning: FixedTarget missing EquatorialCoordinates element.")

        # Save modified XML
        tree.write(output_path, encoding="UTF-8", xml_declaration=True)
        print(f"Updated file written to: {output_path}")



    def sync_passplan_numbers(xml_input, xml_output):
        """
        Updates each <PassPlan Number="X"> so that X matches its TargetSelection Fixed target ID.
        Example:
            <TargetSelection>Fixed: 13</TargetSelection>
            → PassPlan Number becomes "13".
        """

        tree = ET.parse(xml_input)
        root = tree.getroot()

        # Namespace used by Roman APT XML
        ns = {'ns': root.tag.split('}')[0].strip('{')}

        # Iterate through all PassPlan entries
        for passplan in root.findall('.//ns:PassPlan', ns):
            ts = passplan.find('ns:TargetSelection', ns)
            if ts is not None and "Fixed:" in ts.text:
                # Extract target number from "Fixed: N"
                target_num = ts.text.split("Fixed:")[1].strip()
                passplan.set("Number", target_num)

        # Save output
        tree.write(xml_output, encoding="UTF-8", xml_declaration=True)
        print(f"Updated XML saved to {xml_output}")



    def sort_surveyplan_steps(xml_input, xml_output):
        """
        Sorts all <SurveyPlanStep> entries by their <PassPlan> value in increasing order.
        """

        tree = ET.parse(xml_input)
        root = tree.getroot()

        # Extract namespace automatically
        ns = {'ns': root.tag.split('}')[0].strip('{')}

        # Locate the <SurveyPlan> container
        survey_plan = root.find('.//ns:SurveyPlan', ns)
        if survey_plan is None:
            raise RuntimeError("Could not find <SurveyPlan> in XML.")

        # Extract all SurveyPlanStep elements
        steps = survey_plan.findall('ns:SurveyPlanStep', ns)

        # Sort by numeric PassPlan value
        def get_passplan_number(step):
            pp = step.find('ns:PassPlan', ns)
            return int(pp.text.strip()) if pp is not None else 999999999

        steps_sorted = sorted(steps, key=get_passplan_number)

        # Clear existing order
        for step in steps:
            survey_plan.remove(step)

        # Reinsert in sorted order
        for step in steps_sorted:
            survey_plan.append(step)

        # Save output
        tree.write(xml_output, encoding="UTF-8", xml_declaration=True)
        print(f"SurveyPlan sorted and saved to {xml_output}")



    def update_orient_ranges(xml_input, xml_output, orient_min_list, orient_max_list):
        """
        Updates OrientRange OrientMin/OrientMax in each SurveyPlanStep using
        user-provided arrays.
        
        orient_min_list and orient_max_list must have the same length as the
        number of SurveyPlanStep entries.
        """

        tree = ET.parse(xml_input)
        root = tree.getroot()

        # Extract namespace automatically
        ns = {'ns': root.tag.split('}')[0].strip('{')}

        # Find all SurveyPlanStep entries
        steps = root.findall('.//ns:SurveyPlan/ns:SurveyPlanStep', ns)

        if len(steps) != len(orient_min_list) or len(steps) != len(orient_max_list):
            raise ValueError("ERROR: Input arrays must match number of SurveyPlanStep entries.")

        # Update OrientMin and OrientMax for each step
        for step, new_min, new_max in zip(steps, orient_min_list, orient_max_list):
            orient_range = step.find('ns:SpecialRequirements/ns:OrientRange', ns)
            if orient_range is not None:
                orient_range.set("OrientMin", f"{new_min} Degrees")
                orient_range.set("OrientMax", f"{new_max} Degrees")
            else:
                print(f"Warning: No OrientRange found in SurveyPlanStep uid={step.get('uid')}")

        # Save updated XML
        tree.write(xml_output, encoding="UTF-8", xml_declaration=True)
        print(f"Updated XML saved to {xml_output}")


    # Example usage:
    # new_min = [142.1, 143.2, 144.3, ...]
    # new_max = [142.1, 143.2, 144.3, ...]

class query():
    ################
    def stream_roman_mast(products, row=0):
        missions = MastMissions(mission='roman')
        token = os.getenv("MAST_API_TOKEN")    
        missions.login(token=token)
        af = missions.read_product(products[row]['filename'])
        dm = rdm.open(af)
        return(dm)

    ####################

    from astropy.time import Time

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


    def roman_query(criteria=None, coordinates=None, radius=None, filter=None, mindate=None, detector=None, file_suffix=None, exposure_id=None, download_products=False, verbose=False):

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

        if criteria is not None:  # Query with column criteria
            results = missions.query_criteria(**criteria, select_cols=col_list)
            if verbose: print(results)

        if coordinates is not None:
            results = missions.query_region(coordinates=coordinates, radius=radius)
            if verbose: print(results)

        if filter is not None: 
            results = results[results["optical_element"] == filter]
            if verbose: print(results)

        if detector is not None: 
            results = results[results["detector"] == detector]
            if verbose: print(results)
        
        if mindate is not None: 
            results = results[Time(results["exposure_start_time"]) > mindate]
            if verbose: print(results)
        
        # Re-order the column names 
        # results = results[col_list]
        products = missions.get_unique_product_list(results)

        if file_suffix is not None:
            products = missions.filter_products(products, file_suffix='_cal')

        if exposure_id is not None:
            found_index = []
            for i in range(len(products["dataset"])):
                if exposure_id in products["dataset"][i]: found_index.append(i)
            products = products[found_index]

        products, unique_roots = rs.mast.query.pack_exposures(products)

        if download_products:
            downloaded_files = rs.mast.query.download_results(products=products, exposures=unique_roots)
            return(results, products, unique_roots, downloaded_files)
        else:
            return(results, products, unique_roots, None)


    def download_results(products, exposures, verbose=False):
        # Get the filters and prepare the folders 
        list_of_unique_filters = list(set(products["filters"]))
        for bandpass in list_of_unique_filters:
            bandpass_id = "".join(bandpass)
            print(bandpass_id)
            rs.utils.execute_cmd("mkdir " + bandpass_id)
        # Folders ready. 

        # Start download 
        pbar = tqdm(range(len(exposures)))
        
        downloaded_files = []
        for i in pbar:
            outname = exposures[i] + ".fits"
            pbar.set_description(f"Downloading {outname}")
            if verbose: print(outname)
            if not os.path.exists(outname): 
                exposure_products = products[products["exposure_id"] == exposures[i]]
                roman_exposure = rs.core.exposure(exposure_products)
                bandpass = roman_exposure.FILTER
                roman_exposure.save_flc(outname=bandpass+ "/" + outname)
                downloaded_files.append(bandpass+ "/" + outname)
        return(downloaded_files)