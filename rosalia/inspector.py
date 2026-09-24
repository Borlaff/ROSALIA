import os
import glob
import asdf
import s3fs
from tqdm import tqdm
import rosalia as rs
import numpy as np
from astropy.io import fits
import roman_datamodels.datamodels._datamodels
import astropy.table
from collections import defaultdict
from types import SimpleNamespace


############################
def exposure_inspector(input_name, verbose=False, lite=False):
    """
    exposure_inspector: 
    Inspects and returns critical information about the contents of an exposure telescope.
    :param input_name: Exposure file to be inspected. It can be a string with the name of the file, a list of files, or a pattern with *.
    :type input_name: str, list
    :return: Exposure identity - A dictionary with critical information about the exposure, including pointing (right ascension and declination), position angle, telescope, instrument, detector, filter, and a WCS (Astropy and GWCS).
    :rtype: dict, pd.DataFrame
    """

    # If the input is a pattern with *, then we will use glob to find the files 
    # that match the pattern, then we will run convert_ASDF_to_FITS, 
    # and finally run exposure_inspector of the final product. 

    #list_of_files = glob.glob(input_name)
    #if "*" in input_name and input_name.endswith(".asdf"):
    #    output_name = input_name.replace("*", "_").replace(".asdf", ".fits")
    #    input_name = rs.utils.convert_ASDF_to_FITS(asdf_list=list_of_files, output=output_name)
    #    exposure_identity = exposure_inspector_single(input_name, verbose=verbose, lite=lite)
    #    return(exposure_identity)



    # If it is just a string without * wildcard, then we will run exposure_inspector_single directly.
    if isinstance(input_name, (str,)):
        # If we provide a list of files, then we will run exposure_inspector_single for each file,
        # and return a dataframe with the results.
        if "*" in input_name:
            if "s3://" in input_name: # Then is a S3 bucket path. 
                fs = s3fs.S3FileSystem(anon=True)
                s3_files = fs.glob(input_name)
                s3_dir = os.path.dirname(input_name)
                s3_paths = []

                for s3_file in s3_files:
                    s3_paths.append(s3_dir + "/" + os.path.basename(s3_file))
                input_name = s3_paths
                # print(input_name)
            else:
                input_name = glob.glob(input_name)

        else:
            exposure_identity = exposure_inspector_single(input_name, verbose=verbose, lite=lite)
            return(exposure_identity)


    # If the input is a MAST stream, bypass the following. 
    if isinstance(input_name, roman_datamodels.datamodels._datamodels.ImageModel):
        exposure_identity = exposure_inspector_single(input_name, verbose=verbose, lite=lite)
        # If the input is a stream, it needs a local basename. 
        exposure_identity["FILENAME"] = input_name["meta"]["filename"]
        return(exposure_identity)


    # isinstance("Hello", (float, int, str, list, dict, tuple))     
    if isinstance(input_name, (list, astropy.table.Table)):
        """
        DATA = []
        ASTROPYWCS = []
        DATA_SHAPE = []
        SCIEXTS = []
        SCAS = []
        # PHOTOMETRY = []
        conversion_megajanskys = []
        conversion_megajanskys_uncertainty = []
        pixel_area = []
        # {'conversion_megajanskys': 0.6901055870968528,
        # 'conversion_megajanskys_uncertainty': 0.02684817117080749,
        # 'pixel_area': 2.6492509187216547e-13}


        for i in tqdm(range(len(input_name))):
            if isinstance(input_name, (list,)): 
                exposure_identity = rs.inspector.exposure_inspector(input_name[i], lite=lite)
            elif isinstance(input_name, (astropy.table.table.Table,)): 
                data_stream = rs.mast.stream_roman_mast(products=input_name, row=i)
                exposure_identity = rs.inspector.exposure_inspector(data_stream, lite=lite)
            # exposure_identities.append(exposure_identity)
            DATA.append(exposure_identity["DATA"][0])
            DATA_SHAPE.append(exposure_identity["DATA_SHAPE"][0])
            ASTROPYWCS.append(exposure_identity["ASTROPYWCS"][0])
            SCIEXTS.append(exposure_identity["SCIEXTS"][0])
            SCAS.append(exposure_identity["SCA"])
            conversion_megajanskys.append(exposure_identity["photometry"]["conversion_megajanskys"])
            conversion_megajanskys_uncertainty.append(exposure_identity["photometry"]["conversion_megajanskys_uncertainty"])
            pixel_area.append(exposure_identity["photometry"]["pixel_area"])
            # print()
            #PHOTOMETRY.append(exposure_identity["photometry"])


        exposure_identity["DATA"] = DATA
        exposure_identity["DATA_SHAPE"] = DATA_SHAPE
        exposure_identity["ASTROPYWCS"] = ASTROPYWCS
        exposure_identity["SCIEXTS"] = SCIEXTS
        exposure_identity["SCAS"] = SCAS
        #exposure_identity["PHOTOMETRY"] = PHOTOMETRY
        exposure_identity["conversion_megajanskys"] = conversion_megajanskys
        exposure_identity["conversion_megajanskys_uncertainty"] = conversion_megajanskys_uncertainty
        exposure_identity["pixel_area"] = pixel_area

        """
        exposure_identities = []
        for i in tqdm(range(len(input_name))):
            if isinstance(input_name, (list,)): 
                exposure_identity = rs.inspector.exposure_inspector(input_name[i], lite=lite)
                
            elif isinstance(input_name, (astropy.table.table.Table,)): 
                data_stream = rs.mast.stream_roman_mast(products=input_name, row=i)
                exposure_identity = rs.inspector.exposure_inspector(data_stream, lite=lite)
                
            exposure_identities.append(exposure_identity)

        exposure_identities = combine_dicts(exposure_identities)


        return(exposure_identities)


def exposure_inspector_single(input_name, verbose=False, lite=False):
    """
    Inspects and returns critical information about the contents of a telescope
    exposure, stored as a FITS or ADSF file.

    :param input_name: Exposure file to be inspected
    :type cmd: str
    :return: Exposure identity - A dictionary with critical information about the exposure, including pointing (right ascension and declination), position angle, telescope, instrument, detector, filter, and a WCS (Astropy and GWCS).
    :rtype: dict
    """

    # If the input is a MAST stream, bypass the following. 
    import roman_datamodels.datamodels._datamodels
    if isinstance(input_name, roman_datamodels.datamodels._datamodels.ImageModel):
        exposure_identity = exposure_inspector_asdf(input_name, verbose=verbose, lite=lite)
        return(exposure_identity)
    
    # Extract input file extension
    filename, file_extension = os.path.splitext(input_name)



    # If the input image is a FITS file, then use astropy.io.fits.open
    if file_extension == ".fits":
        exposure_identity = exposure_inspector_fits(input_name, verbose=verbose, lite=lite)

    # If the input image is a ASDF file, then use asdf.open
    if file_extension == ".asdf":
        exposure_identity = exposure_inspector_asdf(input_name, verbose=verbose, lite=lite)

    # Add the position of the telescope to the identity.
    try:
        telescope = rs.telescopes.telescope_class_finder(exposure_identity["TELESCOP"])
        exposure_identity["XYZ_HELIO_POS"] = telescope.get_location(mjd=exposure_identity["EXPSTART"])
    except:
        telescope = None
        print("Telescope " + exposure_identity["TELESCOP"] + " Class not found!")
        print("Some functions might not work")

    return(exposure_identity)


def exposure_inspector_asdf(input_name, telescope="roman", verbose=False, lite=False):
    # print(input_name)
    import roman_datamodels.datamodels._datamodels
    if isinstance(input_name, roman_datamodels.datamodels._datamodels.ImageModel):
        input_asdf = input_name

    if isinstance(input_name, str):
        if "s3://" in input_name: 
            if verbose: print("Nexus S3 bucket file detected")

            fs = s3fs.S3FileSystem(anon=True)
            input_asdf = asdf.open(fs.open(input_name, 'rb'))[telescope]

        else:
            if verbose: print("Local ASDF file detected")
            input_asdf = asdf.open(input_name)[telescope]

    # Setting up keywords to store the info from the file
    exposure_identity = {}
    exposure_identity["FILENAME"] = input_name

    #keywords = ["TELESCOP", "INSTRUME", "DETECTOR", "RA_TARG", "DEC_TARG", "SUNANGLE", "BUNIT",
    #            "EXPSTART", "EXPEND", "EXPTIME", "MOONANGL", "DRIZCORR",
    #            "PHOTCORR", "PHOTFLAM", "PHOTPLAM"]

    exposure_identity["INSTRUME"] = input_asdf["meta"]["instrument"]["name"]

    exposure_identity["TELESCOP"] = input_asdf["meta"]["telescope"]
    if exposure_identity["TELESCOP"] == "ROMAN":
        # telescope_class = rs.telescopes.Roman
        detector_svo = "WFI"

    print('Location is stored in data_stream["meta"]["ephemeris"]! ')

    exposure_identity["DETECTOR"]   = input_asdf["meta"]["instrument"]["detector"]
    exposure_identity["FILTER"]     = input_asdf["meta"]["instrument"]["optical_element"]

    exposure_identity["RA_TARG"]    = input_asdf["meta"]["pointing"]["target_ra"]
    exposure_identity["DEC_TARG"]   = input_asdf["meta"]["pointing"]["target_dec"]
    exposure_identity["photometry"] = input_asdf["meta"]["photometry"]# ["conversion_megajanskys"]
    # exposure_identity["photometry"] = input_asdf["meta"]["photometry"]["conversion_megajanskys_uncertainty"]
    # exposure_identity["pixel_area"] = input_asdf["meta"]["photometry"]["pixel_area"]
    # data_stream["meta"]["photometry"]["conversion_megajanskys_uncertainty"]
    #exposure_identity["SUNANGLE"]  = input_asdf["roman"]["meta"]["ephemeris"]["sun_angle"]
    #exposure_identity["MOONANGL"]  = input_asdf["roman"]["meta"]["ephemeris"]["moon_angle"]
    #exposure_identity["EARTHANGL"] = input_asdf["roman"]["meta"]["ephemeris"]["earth_angle"]
    exposure_identity["EXPSTART"]   = input_asdf["meta"]["exposure"]["start_time"].mjd
    exposure_identity["EXPEND"]     = input_asdf["meta"]["exposure"]["end_time"].mjd
    exposure_identity["EXPTIME"]    = input_asdf["meta"]["exposure"]["exposure_time"]
    exposure_identity["conversion_megajanskys"] = input_asdf["meta"]["photometry"]["conversion_megajanskys"]
    exposure_identity["conversion_megajanskys_uncertainty"] = input_asdf["meta"]["photometry"]["conversion_megajanskys_uncertainty"]
    exposure_identity["pixel_area"] = input_asdf["meta"]["photometry"]["pixel_area"]*((180/np.pi)*60*60)**2
    exposure_identity["EXPSTART_ISOT"] = input_asdf["meta"]["exposure"]["start_time"].isot
    exposure_identity["SCA"]        = int(input_asdf["meta"]["instrument"]["detector"].replace("WFI",""))
    exposure_identity["BUNIT"]      = "DN/s" # input_asdf["roman"]["meta"]["photometry"]["flux_unit"]

    ############## Get the filter identity #######################
    try: 
        exposure_identity["FILTER_IDENTITY"] = rs.telescopes.find_filter_in_svo(wavelength=exposure_identity["FILTER"],
                                                                                telescope=exposure_identity["TELESCOP"],
                                                                                instrument=exposure_identity["INSTRUME"],
                                                                                detector=detector_svo,
                                                                                verbose=verbose)
        
    except:
        print("The filter " + exposure_identity["FILTER"] + "/" + exposure_identity["TELESCOP"] +  "/" +  exposure_identity["INSTRUME"] + "/" + detector_svo + " was not found. Photometric calculations can be compromised.")
        print("Check FILTER keyword in the Science Headers and / or proceed with caution.")
        exposure_identity["FILTER_IDENTITY"] = None
    ##############################################################

    # exposure_identity["PHYSPIX"] = telescope_class.get_physical_pixelsize(instrument=exposure_identity["INSTRUME"])
    exposure_identity["PIXSCALE"] = np.sqrt(exposure_identity["pixel_area"])

    # ------------- #
    data = []
    gwcs = []
    astropywcs = []
    from astropy.wcs import WCS as astropy_wcs
    #  nSCAs = 1 # Right now (October 2024) exposure inspector only accepts Roman/WFI images with one SCA per ASDF file.
    data.append(np.array(input_asdf["data"]))
    gwcs.append(input_asdf["meta"]["wcs"])
    astropywcs.append(astropy_wcs(input_asdf["meta"]["wcs"].to_fits()[0]))

    exposure_identity["DATA"] = data[0]
    exposure_identity['DATA_SHAPE'] = data[0].shape
    exposure_identity["GWCS"] = gwcs[0]
    exposure_identity["ASTROPYWCS"] = astropywcs[0]

    exposure_identity["RA_PNT"] =  input_asdf["meta"]['pointing']['ra_v1']
    exposure_identity["DEC_PNT"] =  input_asdf["meta"]['pointing']['dec_v1']
    # exposure_identity["PA"] =  # input_asdf["roman"]["meta"]['pointing']["pa_v3"] - 60  # input_asdf["roman"]["meta"]['pointing']["pa_aperture"] # input_asdf["roman"]["meta"]['pointing']['pa_v3']
    exposure_identity["PA"] = input_asdf["meta"]['pointing']["pa_aperture"]
    exposure_identity["FILETYPE"] = "ASDF"
    exposure_identity["SCIEXTS"] = int(input_asdf["meta"]["instrument"]["detector"].replace("WFI",""))
    return(exposure_identity)



def exposure_inspector_fits(input_name, verbose=False, lite=False):
    from astropy.wcs import WCS

    exposure_identity = {}

    # Setting up keywords to store the info from the file
    exposure_identity["FILENAME"] = input_name

    keywords_ext0 = ["TELESCOP", "INSTRUME", "DETECTOR", "RA_TARG", "DEC_TARG", "SUNANGLE",
                     "SUN_ALT", "BUNIT",
                     "EXPSTART", "EXPEND", "EXPTIME", "MOONANGL", "DRIZCORR",
                     "PHOTCORR"]

    keywords_ext1 = ["BUNIT", "PHOTFLAM", "PHOTPLAM"]

    input_fits = fits.open(input_name)

    for keyword_i in keywords_ext0:
        try:
            exposure_identity[keyword_i] = input_fits[0].header[keyword_i]
        except:
            if verbose: print(rs.plots.style.YELLOW + "WARNING: KEYWORD " + keyword_i + " not found in ext 0" + rs.plots.style.RESET)
            if verbose: print(rs.plots.style.YELLOW + "Setting to None" + rs.plots.style.RESET)

    for keyword_i in keywords_ext1:
        try:
            exposure_identity[keyword_i] = input_fits[1].header[keyword_i]
        except:
            if verbose: print(rs.plots.style.YELLOW + "WARNING: KEYWORD " + keyword_i + " not found in ext 1" + rs.plots.style.RESET)
            if verbose: print(rs.plots.style.YELLOW + "Setting to None" + rs.plots.style.RESET)

    # Get the exposure time EXPSTART
    from astropy.time import Time
    t = Time(exposure_identity["EXPSTART"], format='mjd', scale='utc')
    exposure_identity["EXPSTART_ISOT"] = t.isot

    if exposure_identity["TELESCOP"] == "HST" or exposure_identity["TELESCOP"]=="Hubble":
        telescope_class = rs.telescopes.Hubble
        try: 
            exposure_identity["PA"] = input_fits[1].header["PA_APER"]
        except:
            print("Warning: HST image - PA_APER not found.")


    if exposure_identity["TELESCOP"] == "Roman" or exposure_identity["TELESCOP"] == "ROMAN" or exposure_identity["TELESCOP"]=="RST" or exposure_identity["TELESCOP"]=="NGRST":
        telescope_class = rs.telescopes.Roman
        exposure_identity["PA"] = input_fits[0].header["PA"]

        SCAs = []
        for i in range(len(input_fits)-1):
            SCAs.append(input_fits[i+1].header["SCA"])
        exposure_identity["SCA"] = SCAs

    if exposure_identity["TELESCOP"] == "CSST":
        telescope_class = rs.telescopes.CSST
        exposure_identity["PA"] = input_fits[0].header["PA"]

    if exposure_identity["TELESCOP"] == "SPHEREx":
        telescope_class = rs.telescopes.SPHEREx
        exposure_identity["PA"] = input_fits[0].header["PA"]

    if exposure_identity["TELESCOP"] == "ARRAKIHS":
        telescope_class = rs.telescopes.ARRAKIHS
        exposure_identity["PA"] = input_fits[0].header["PA"]

    # If this an flc image ?
    exposure_identity["HST_TYPE"] = "GEN"
    try:
        bunit = exposure_identity["BUNIT"]
        if (exposure_identity["DRIZCORR"] == "PERFORM") & (exposure_identity["BUNIT"] == "ELECTRONS"):
            exposure_identity["HST_TYPE"] = "FLT"

        if (exposure_identity["DRIZCORR"] == "COMPLETE") & (exposure_identity["BUNIT"] == "ELECTRONS/S"):
            exposure_identity["HST_TYPE"] = "DRZ"

    except:
        if verbose>1: print(rs.plots.style.YELLOW + "BUNIT not found!" + rs.plots.style.RESET)
        if (exposure_identity["HST_TYPE"] != "DRZ") or (exposure_identity["HST_TYPE"] != "FLT"):
            exposure_identity["HST_TYPE"] = "GEN"




    # Get the right filter. 
    if "FILTER1" in input_fits[0].header: 
        if "CLEAR" not in input_fits[0].header["FILTER1"]:
            filter = input_fits[0].header["FILTER1"]
    if "FILTER2" in input_fits[0].header: 
        if "CLEAR" not in input_fits[0].header["FILTER2"]:
            filter = input_fits[0].header["FILTER2"]
    else:
        try: 
            filter = input_fits[0].header["FILTER"]
        except: 
            print("ERROR! The FITS does not have a FILTER keyword in the header of EXT 0. ")
    exposure_identity["FILTER"] = filter

    ############## Get the filter identity #######################
    try: 
        exposure_identity["FILTER_IDENTITY"] = rs.telescopes.find_filter_in_svo(wavelength=exposure_identity["FILTER"],
                                                                            telescope=exposure_identity["TELESCOP"],
                                                                            instrument=exposure_identity["INSTRUME"],
                                                                            detector=exposure_identity["DETECTOR"],
                                                                            verbose=False)
        
    except:
        print("The filter " + exposure_identity["FILTER"] + "/" + exposure_identity["TELESCOP"] +  "/" +  exposure_identity["INSTRUME"] + "/" + exposure_identity["DETECTOR"] + " was not found. Photometric calculations can be compromised.")
        print("Check FILTER keyword in the Science Headers and / or proceed with caution.")
        exposure_identity["FILTER_IDENTITY"] = None
    ##############################################################

    # Get pixel size
    # exposure_identity["PHYSPIX"] = telescope_class.get_physical_pixelsize(instrument=exposure_identity["INSTRUME"])
    # exposure_identity["PIXSCALE"] = telescope_class.get_pixscale(instrument=exposure_identity["INSTRUME"])

    # Find which extensions are SCI
    exposure_identity["SCIEXTS"] = detect_sci_extensions(input_name)

    # Emancipate this as a separate program.
    # Find the central coordinates of the Multiextension fits.

    data = []
    data_shape = []
    astropywcs = []
    science_extension_hdulist = []
    header_list = []

    # If LITE, fill this anyways.
    for sci_ext_i in exposure_identity["SCIEXTS"]:
        data_shape.append(input_fits[sci_ext_i].data.shape)
        header_i = input_fits[sci_ext_i].header
        header_i["EXTNAME"] = "SCI"
        header_i["SCA"] = sci_ext_i      
        # print("Hey!")
        # print(input_fits[sci_ext_i].header) 
        header_list.append(input_fits[sci_ext_i].header)
        astropywcs_i = WCS(header=input_fits[sci_ext_i], fobj=input_fits, naxis=2)
        astropywcs.append(astropywcs_i)
    exposure_identity["HEADERS"] = header_list
    exposure_identity["DATA_SHAPE"] = data_shape
    exposure_identity["ASTROPYWCS"] = astropywcs
    exposure_identity["PIXSCALE"] = np.abs(astropywcs[0].proj_plane_pixel_scales()[0])

    # If not lite, do one more loop with the data to make a swarp coadd.
    if not "RA_TARG" in exposure_identity:
        lite = False
        print("RA_TARG and DEC_TARG keywords not present in header")
        print("Enabling mosaic generator to find the central coordinates of the image")

    if not "PA" in exposure_identity:
        wcs = exposure_identity["ASTROPYWCS"][0]
        exposure_identity["PA"] = np.arctan2(-wcs.wcs.cd[1, 0], wcs.wcs.cd[1, 1])*180.0/np.pi


    if not lite: # Avoid generating the mosaic header.
        for sci_ext_i in exposure_identity["SCIEXTS"]:
            #data_shape.append(input_fits[sci_ext_i].data.shape)
            science_extension_hdulist.append(input_fits[sci_ext_i])
            data.append(input_fits[sci_ext_i].data)
            #astropywcs.append(astropy_wcs.WCS(input_fits[sci_ext_i].header, input_fits))

        exposure_identity["DATA"] = data
        #exposure_identity["DATA_SHAPE"] = data_shape
        #exposure_identity["ASTROPYWCS"] = astropywcs

        try:
            from reproject.mosaicking import find_optimal_celestial_wcs
            wcs_out, shape_out = find_optimal_celestial_wcs(input_data=input_fits, hdu_in=exposure_identity["SCIEXTS"])
            reference_header = wcs_out.to_header()

        except:
            #print("WARNING: The input WCS has distortion. find_optimal_celestial_wcs does not currently support this mode.")
            #print("ROSALIA will apply a temporary solution with SWARP // Jul 29 2024.")

            os.system("swarp -dd > swarp.conf")
            swarp_cmd_str = ""

            for sci_ext_i in exposure_identity["SCIEXTS"]:
                swarp_cmd_str = swarp_cmd_str + '"' + input_name+"["+str(sci_ext_i)  + ']" '

            outname, outname_scaled = run_swarp(pattern=swarp_cmd_str, 
                                                outname=input_name.replace(".fits", "_drz.fits"),
                                                scale=0.11,
                                                coveredfrac=1)


            swarp_coadd = fits.open(outname)
            reference_header = swarp_coadd[0].header


        exposure_identity["RA_PNT"]  = reference_header["CRVAL1"]
        exposure_identity["DEC_PNT"] = reference_header["CRVAL2"]
        exposure_identity["RA_TARG"]  = reference_header["CRVAL1"]
        exposure_identity["DEC_TARG"] = reference_header["CRVAL2"]
        exposure_identity["X_PNT"]  = reference_header["CRPIX1"]
        exposure_identity["Y_PNT"]  = reference_header["CRPIX2"]

    
    ################
    exposure_identity["FILETYPE"] = "FITS"

    return(exposure_identity)




######################################################
######################################################
######################################################


def combine_dicts(items: list[dict]) -> dict:
    """Combines a list of dictionaries sharing the same keys

    into a single dictionary where each key maps to a list of values.
    """
    # Guard against an empty list input
    if not items:
        return {}

    # Iterate over the keys of the first dictionary and gather values from all dictionaries
    return {key: [d[key] for d in items] for key in items[0]}

import numpy as np


def values_are_equal(v1, v2) -> bool:
    """Recursively checks if two values are equal, safely handling

    NumPy arrays, PyTorch/TensorFlow tensors, Pandas structures, and nested containers.
    """
    # 1. Check if both are the exact same instance in memory (fast-path)
    if v1 is v2:
        return True

    # 2. Check type matching (allow numeric scalar types to compare cleanly)
    if type(v1) is not type(v2) and not (
        isinstance(v1, (int, float, complex))
        and isinstance(v2, (int, float, complex))
    ):
        return False

    # 3. Handle NumPy Arrays
    if isinstance(v1, np.ndarray):
        return np.array_equal(v1, v2, equal_nan=True)

    # 4. Handle PyTorch Tensors (if torch is installed)
    if "torch" in type(v1).__module__:
        import torch

        if isinstance(v1, torch.Tensor):
            return (
                v1.shape == v2.shape
                and bool(torch.equal(v1, v2))
                or bool(torch.all(torch.isnan(v1) == torch.isnan(v2)))
            )

    # 5. Handle Pandas DataFrame / Series (if pandas is installed)
    if "pandas" in type(v1).__module__:
        import pandas as pd

        if isinstance(v1, (pd.DataFrame, pd.Series)):
            return v1.equals(v2)

    # 6. Handle nested Dictionaries
    if isinstance(v1, dict):
        if v1.keys() != v2.keys():
            return False
        return all(values_are_equal(v1[k], v2[k]) for k in v1)

    # 7. Handle nested Lists / Tuples
    if isinstance(v1, (list, tuple)):
        if len(v1) != len(v2):
            return False
        return all(values_are_equal(item1, item2) for item1, item2 in zip(v1, v2))

    # 8. Handle Sets
    if isinstance(v1, (set, frozenset)):
        return v1 == v2

    # 9. Fallback comparison wrapped with array check
    try:
        result = v1 == v2
        # If the result itself is an array or collection of bools
        if isinstance(result, np.ndarray):
            return bool(result.all())
        return bool(result)
    except (ValueError, TypeError):
        return False


def are_all_dicts_equal(dict_list: list[dict]) -> bool:
    """Compares all dictionaries in a list against the first element recursively."""
    if not dict_list:
        return True

    ref_dict = dict_list[0]
    return all(values_are_equal(d, ref_dict) for d in dict_list)


def longest_common_substring(strings: list[str]) -> str:
    """Finds the longest contiguous substring common to all strings."""
    if not strings:
        return ""

    shortest_str = min(strings, key=len)
    substr = shortest_str

    # Narrow down the substring until all items contain it
    for length in range(len(shortest_str), 0, -1):
        for start in range(len(shortest_str) - length + 1):
            candidate = shortest_str[start : start + length]
            if all(candidate in s for s in strings):
                return candidate
    return ""
