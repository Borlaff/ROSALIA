# Alejandro S. Borlaff. NASA Ames Research Center. a.s.borlaff@nasa.gov / asborlaff@gmail.com
# January 20, 2023.
#
# ROSALIA/UTILS module
# This module will hold all the programs not related any other task, mainly system tools.
#
# Version log:
# v.1.0 - 20 Enero 2023. First loading of programs inherited from former monolithic straycor.py
# v.2.0 - March 4, 2026. Reorganization of the code into a more modular structure. 
#  This module utils.py is created to hold all the programs not related any other task,
#  mainly system tools.
##########################################################
import os
import glob
import numpy as np
from astropy.io import fits
import pandas as pd
import astropy.units as u
import astropy.wcs as astropy_wcs
from astropy.coordinates import SkyCoord  # High-level coordinates
from tqdm import tqdm
import bottleneck as bn
import rosalia as rs

# print("CHECK THE FLUX OF THE ZODIACAL LIGHT IN A NORMAL IMAGE!!!!!")


###############################
# MULTIORDER HEALPIX ROUTINES #
###############################

#################################

def hp_resol2nside(resolution):
    """
    This program converts a resolution in degrees to the nside parameter of healpix.
    :param resolution: Resolution in degrees (astropy units).
    :type resolution: astropy.units.Unit
    :return: nside parameter of healpix
    :rtype: int
    """
    # Input: resolution in degrees (astropy units).
    # Output: int: nside healpix parameter

    import healpy as hp
    from astropy import units as u
    i = 0
    hp_resolution = np.inf*u.degree
    while hp_resolution > resolution:
        nside = 2**i
        hp_resolution = np.degrees(hp.pixelfunc.nside2resol(nside=nside)*u.radian)
        i = i + 1
    return(nside)

#################################

""" DEPRECATED
def make_allsky_MOC(min_resolution):
    import mhealpy as hmap
    from astropy import units as u
    nside = hp_resol2nside(min_resolution)
    MOC = hmap.HealpixMap(nside=nside, density=True)
    return({"MOC": MOC, "nside": nside})
"""
#################################

""" DEPRECATED
def make_polygon_MOC(ra_vert, dec_vert, min_resolution):
    import healpy as hp
    import mhealpy as hmap
    from astropy import units as u

    #min_resolution = 10*u.arcsec # arcsec
    #ra_vert = np.array([-0.1, -0.1, 0.1, 0.1])
    #dec_vert = np.array([-0.1, 0.1, 0.1, -0.1])

    nside = hp_resol2nside(min_resolution)
    mEq = hmap.HealpixBase(nside = nside)

    vec = hp.ang2vec(ra_vert[::-1], dec_vert[::-1], lonlat=True)
    polygon_pix = mEq.query_polygon(vec)

    #print("Making MOC from pixels 1")
    MOC = hmap.HealpixMap.moc_from_pixels(mEq.nside, polygon_pix, density=True)

    return({"MOC": MOC, "nside": nside})
"""

###############################

def interpolate_location_in_fits(fits_name, ext, ra, dec):
    """
    This program interpolates the value of a FITS file at the location of the input coordinates (ra, dec). The FITS file requires to have a WCS in the correct EXT, and the coordinates are transformed to pixel coordinates using the WCS. The interpolation is done using a linear interpolation in the grid of the FITS file. The output is the interpolated value at the input coordinates.
    
    :param fits_name: Name of the FITS file to be interpolated
    :type fits_name: str
    :param ext: Extension of the FITS file to be interpolated
    :type ext: int
    :param ra: Right ascension of the location to be interpolated
    :type ra: float
    :param dec: Declination of the location to be interpolated
    :type dec: float
    :return: Interpolated value at the input coordinates
    :rtype: float

    """

    # Open the fits file
    fits_file = fits.open(fits_name)


    # Let's interpolate the NDI in the grid.
    import scipy.interpolate as spint
    LGI = spint.LinearNDInterpolator

    x = np.linspace(0, fits_file[ext].header["NAXIS1"]-1, fits_file[ext].header["NAXIS1"]) #  or  0.5*np.arange(3.) works too
    y = np.linspace(0, fits_file[ext].header["NAXIS2"]-1, fits_file[ext].header["NAXIS2"]) #  or  0.5*np.arange(3.) works too

    # populate the 3D array of values (re-using x because lazy)
    X, Y = np.meshgrid(x, y, indexing='ij')
    Z = fits_file[ext].data

    # Make the interpolator, (list of 1D axes, values at all points)
    lgi = LGI(points=list(zip(X.flatten(), Y.flatten())), values=Z.flatten())  # can also be [x]*3 or (x,)*3

    # Find the pixel locations that we want to interpolate of the input coordinates (world) using the WCS.
    x2interp, y2interp = rs.utils.radec_to_xy(ra, dec, fits_name, ext)
    # xy2interp = np.hstack([x2interp[:, np.newaxis], y2interp[:, np.newaxis]])
    # print(x2interp, y2interp)
    Z_interp = lgi(x2interp, y2interp)
    return(Z_interp)


############################

def load_dict(input_name, verbose=False):
    """
    Load a python object using pickle.
        :param str input_name: Name of the file to be loaded
        :returns: Python object stored in the file
        :rtype: Python object
    """
    import dill as pickle
    with open(input_name, 'rb') as f:
        if verbose: print("Loading " + input_name)
        dictionary = pickle.load(f)
    return(dictionary)

def save_dict(dictionary, input_name, verbose=False):
    """
    Save a python object using pickle.
        :param Python object dictionary: Any type of python object
        :returns: Input python object
        :rtype: Python object
    """
    import dill as pickle

    with open(input_name, 'wb') as f:
        if verbose: print("Saving " + input_name)
        pickle.dump(dictionary, f)
    return(dictionary)
############################

def download_file(url):
    """
    Downloads the file described as an url address.
        :param str url: Url of the file to be downloaded
        :returns: Local path of the downloaded file
        :rtype: str
    """
    import requests
    # From: https://stackoverflow.com/questions/16694907/download-large-file-in-python-with-requests
    local_filename = url.split('/')[-1]
    # NOTE the stream=True parameter below
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                # If you have chunk encoded response uncomment if
                # and set chunk_size parameter to None.
                #if chunk:
                f.write(chunk)
    return(local_filename)

############################
def execute_cmd(cmd, verbose=False):
    """
    Executes a shell command line.
    :param cmd: Command line to be executed
    :type cmd: str
    :return: Shell output
    :rtype: str
    """
    import subprocess
    log_error = open('subprocess.out', "w")
    # output = subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    #result = subprocess.run(cmd, capture_output=True, text=True)
    # result = subprocess.Popen([cmd], shell=True, text=True)
    if verbose:
        print(cmd)
    try:
        print(rs.plots.style.YELLOW)
        result = subprocess.check_output([cmd], shell=True, text=True, stderr=subprocess.STDOUT)
        print(rs.plots.style.RESET)
        return(result)

    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)
        return(None)

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
    list_of_files = glob.glob(input_name)
    if "*" in input_name and input_name.endswith(".asdf"):
        output_name = input_name.replace("*", "_").replace(".asdf", ".fits")
        input_name = rs.utils.convert_ASDF_to_FITS(asdf_list=list_of_files, output=output_name)
        exposure_identity = exposure_inspector_single(input_name, verbose=verbose, lite=lite)
        return(exposure_identity)

    # If it is just a string without * wildcard, then we will run exposure_inspector_single directly.
    if isinstance(input_name, (str,)):
        exposure_identity = exposure_inspector_single(input_name, verbose=verbose, lite=lite)
        return(exposure_identity)

    # If we provide a list of files, then we will run exposure_inspector_single for each file,
    # and return a dataframe with the results.
    if isinstance(input_name, (list,)):
        exposure_identities = []
        for i in tqdm(range(len(input_name))):
            exposure_identity = rs.utils.exposure_inspector(input_name[i], lite=lite)
            exposure_identities.append(exposure_identity)

        if exposure_identity["TELESCOP"] == "ROMAN":
            print(rs.plots.style.YELLOW)
            print("WARNING: Nancy Grace Roman Space Telescope has not been launched yet.")
            print("Approximating the location by using JWST's location in JPL/Horizons.")
            print(rs.plots.style.RESET)

        return(pd.DataFrame(exposure_identities))


def exposure_inspector_single(input_name, verbose=False, lite=False):
    """
    Inspects and returns critical information about the contents of a telescope
    exposure, stored as a FITS or ADSF file.

    :param input_name: Exposure file to be inspected
    :type cmd: str
    :return: Exposure identity - A dictionary with critical information about the exposure, including pointing (right ascension and declination), position angle, telescope, instrument, detector, filter, and a WCS (Astropy and GWCS).
    :rtype: dict
    """
    # Extract input file extension
    filename, file_extension = os.path.splitext(input_name)

    # If the input image is a FITS file, then use astropy.io.fits.open
    if file_extension == ".fits":
        exposure_identity = exposure_inspector_fits(input_name, verbose=verbose, lite=lite)
        exposure_identity["FILETYPE"] = "FITS"

    # If the input image is a ASDF file, then use asdf.open
    if file_extension == ".asdf":
        exposure_identity = exposure_inspector_asdf(input_name, verbose=verbose, lite=lite)
        exposure_identity["FILETYPE"] = "ASDF"
        exposure_identity["SCIEXTS"] = np.array([0])

    # Add the position of the telescope to the identity.
    telescope = rs.telescopes.telescope_class_finder(exposure_identity["TELESCOP"])
    exposure_identity["XYZ_HELIO_POS"] = telescope.get_location(mjd=exposure_identity["EXPSTART"])

    return(exposure_identity)


def exposure_inspector_asdf(input_name, verbose=False, lite=False):
    import asdf
    input_asdf = asdf.open(input_name)

    # Setting up keywords to store the info from the file
    exposure_identity = {}
    exposure_identity["FILENAME"] = input_name

    keywords = ["TELESCOP", "INSTRUME", "DETECTOR", "RA_TARG", "DEC_TARG", "SUNANGLE", "BUNIT",
                "EXPSTART", "EXPEND", "EXPTIME", "MOONANGL", "DRIZCORR",
                "PHOTCORR", "PHOTFLAM", "PHOTPLAM"]

    exposure_identity["INSTRUME"] = input_asdf["roman"]["meta"]["instrument"]["name"]

    exposure_identity["TELESCOP"] = input_asdf["roman"]["meta"]["telescope"]
    if exposure_identity["TELESCOP"] == "ROMAN":
        telescope_class = rs.telescopes.Roman
        detector_svo = "WFI"

    exposure_identity["DETECTOR"] = input_asdf["roman"]["meta"]["instrument"]["detector"]
    exposure_identity["FILTER"] = input_asdf["roman"]["meta"]["instrument"]["optical_element"]

    exposure_identity["RA_TARG"] = input_asdf["roman"]["meta"]["pointing"]["target_ra"]
    exposure_identity["DEC_TARG"] = input_asdf["roman"]["meta"]["pointing"]["target_dec"]
    #exposure_identity["SUNANGLE"] = input_asdf["roman"]["meta"]["ephemeris"]["sun_angle"]
    #exposure_identity["MOONANGL"] = input_asdf["roman"]["meta"]["ephemeris"]["moon_angle"]
    #exposure_identity["EARTHANGL"] = input_asdf["roman"]["meta"]["ephemeris"]["earth_angle"]
    exposure_identity["EXPSTART"] = input_asdf["roman"]["meta"]["exposure"]["start_time"].mjd
    exposure_identity["EXPEND"] = input_asdf["roman"]["meta"]["exposure"]["end_time"].mjd
    exposure_identity["EXPTIME"] = input_asdf["roman"]["meta"]["exposure"]["exposure_time"]
    exposure_identity["conversion_megajanskys"] = input_asdf["roman"]["meta"]["photometry"]["conversion_megajanskys"]
    exposure_identity["pixel_area"] = input_asdf["roman"]["meta"]["photometry"]["pixel_area"]*((180/np.pi)*60*60)**2
    exposure_identity["EXPSTART_ISOT"] = input_asdf["roman"]["meta"]["exposure"]["start_time"].isot
    exposure_identity["SCA"] = int(input_asdf["roman"]["meta"]["instrument"]["detector"].replace("WFI",""))


    ############## Get the filter identity #######################

    exposure_identity["FILTER_IDENTITY"] = rs.telescopes.find_filter_in_svo(wavelength=exposure_identity["FILTER"],
                                                                            telescope=exposure_identity["TELESCOP"],
                                                                            instrument=exposure_identity["INSTRUME"],
                                                                            detector=detector_svo,
                                                                            verbose=False)

    ##############################################################

    exposure_identity["PHYSPIX"] = telescope_class.get_physical_pixelsize(instrument=exposure_identity["INSTRUME"])
    exposure_identity["PIXSCALE"] = np.sqrt(exposure_identity["pixel_area"])

    # ------------- #
    data = []
    gwcs = []
    astropywcs = []
    from astropy.wcs import WCS as astropy_wcs
    nSCAs = 1 # Right now (October 2024) exposure inspector only accepts Roman/WFI images with one SCA per ASDF file.
    for i in range(nSCAs):
        data.append(np.array(input_asdf["roman"]["data"]))
        gwcs.append(input_asdf["roman"]["meta"]["wcs"])
        astropywcs.append(astropy_wcs(input_asdf["roman"]["meta"]["wcs"].to_fits()[0]))

    exposure_identity["DATA"] = data
    exposure_identity["GWCS"] = gwcs
    exposure_identity["ASTROPYWCS"] = astropywcs

    exposure_identity["RA_PNT"] =  input_asdf["roman"]["meta"]['pointing']['ra_v1']
    exposure_identity["DEC_PNT"] =  input_asdf["roman"]["meta"]['pointing']['dec_v1']
    exposure_identity["PA"] = input_asdf["roman"]["meta"]['pointing']["pa_aperture"] # input_asdf["roman"]["meta"]['pointing']['pa_v3']



    return(exposure_identity)


def exposure_inspector_fits(input_name, verbose=False, lite=False):
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
        exposure_identity["PA"] = input_fits[1].header["PA_APER"]


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




    # Get the right filter
    filter_1 = input_fits[0].header["FILTER1"]
    filter_2 = input_fits[0].header["FILTER2"]

    if "CLEAR" in filter_1:
        exposure_identity["FILTER"] = filter_2
    else:
        exposure_identity["FILTER"] = filter_1

    # Get the filter identity
    exposure_identity["FILTER_IDENTITY"] = rs.telescopes.find_filter_in_svo(wavelength=exposure_identity["FILTER"],
                                                              telescope=exposure_identity["TELESCOP"],
                                                              instrument=exposure_identity["INSTRUME"],
                                                              detector=exposure_identity["DETECTOR"],
                                                              verbose=False)

    # Get pixel size
    exposure_identity["PHYSPIX"] = telescope_class.get_physical_pixelsize(instrument=exposure_identity["INSTRUME"])
    exposure_identity["PIXSCALE"] = telescope_class.get_pixscale(instrument=exposure_identity["INSTRUME"])

    # Find which extensions are SCI
    exposure_identity["SCIEXTS"] = detect_sci_extensions(input_name)

    # Emancipate this as a separate program.
    # Find the central coordinates of the Multiextension fits.

    data = []
    data_shape = []
    astropywcs = []
    science_extension_hdulist = []

    # If LITE, fill this anyways.
    for sci_ext_i in exposure_identity["SCIEXTS"]:
        data_shape.append(input_fits[sci_ext_i].data.shape)
        astropywcs.append(astropy_wcs.WCS(input_fits[sci_ext_i].header, input_fits))
    exposure_identity["DATA_SHAPE"] = data_shape
    exposure_identity["ASTROPYWCS"] = astropywcs

    # If not lite, do one more loop with the data to make a swarp coadd.
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
                                                outname=input_name.replace(".fits", "_swarp_coadd.fits"),
                                                coveredfrac=1)
            #swarp_coadd_name =
            #swarp_cmd_str = swarp_cmd_str + ' -SUBTRACT_BACK N -BLANK_BADPIXELS Y -VERBOSE_TYPE QUIET -IMAGEOUT_NAME "' + input_name.replace(".fits", "_swarp_coadd.fits") + '"'

            #if verbose: print(swarp_cmd_str)
            #execute_cmd(swarp_cmd_str)

            swarp_coadd = fits.open(outname)
            reference_header = swarp_coadd[0].header


        exposure_identity["RA_PNT"]  = reference_header["CRVAL1"]
        exposure_identity["DEC_PNT"] = reference_header["CRVAL2"]
        exposure_identity["X_PNT"]  = reference_header["CRPIX1"]
        exposure_identity["Y_PNT"]  = reference_header["CRPIX2"]
    ################

    return(exposure_identity)



def run_swarp(pattern, outname, coveredfrac=1, verbose=False):
    rs.utils.execute_cmd("swarp -d > swarp.conf", verbose=verbose) # Generate a default config file for swarp
    rs.utils.execute_cmd("swarp -c swarp.conf -SUBTRACT_BACK N -BLANK_BADPIXELS Y -VERBOSE_TYPE QUIET " + pattern, verbose=verbose) # Run swarp on all the SCAs
    rs.utils.execute_cmd("mv coadd.fits " + outname, verbose=verbose) # Make a compressed version, for easiest visualization.
    outname_warped = outname.replace(".fits", "_warped.fits")
    rs.utils.execute_cmd("astwarp " + outname +" -h0 --coveredfrac=" + str(coveredfrac) + " --scale=0.1,0.1 --output=" + outname_warped, verbose=verbose) # Make a compressed version, for easiest visualization.
    outname_scaled = outname.replace(".fits", "_scaled.fits")
    rs.utils.execute_cmd("astarithmetic -h1 " + outname_warped + " 0.01 x --output=" + outname_scaled, verbose=verbose) # Make a compressed version, for easiest visualization.
    rs.utils.execute_cmd("rm " + outname_warped, verbose=verbose)
    print("Result in " + outname + " & " + outname_scaled)
    return([outname, outname_scaled])


#####################################################################

############################
def check_number_of_extensions(input_name, expected_next, verbose=False):
    # This program will check that the input fits
    # contains a number of extensions at least as large as the expected_next

    input_fits = fits.open(input_name)
    nextensions = len(input_fits)
    if nextensions < expected_next:
        print(r'\n\n\n\n\n')
        print("ERROR: ")
        print(input_name + " does not contain the expected " + str(expected_next) + " number of extensions")
        print("Please verify")
        print(r'\n\n\n\n\n')
    else:
        if verbose:
            print(input_name + " " + str(nextensions) + " # extensions OK")


def get_parameters_list(fits_list, index, ext=0):
    # Check that fits_list are in array
    if not isinstance(fits_list, (list, pd.core.series.Series, np.ndarray)):
        fits_list = np.array([fits_list])

    # If index not in array, copy their values into one
    # as large asfits_list
    if not isinstance(index, (list, pd.core.series.Series, np.ndarray)):
        index = np.array([index]*len(fits_list))

    PARAM = []
    for j in range(len(index)):
        PARAM.append([])
    for raw_name in fits_list:
        print(raw_name)
        raw_fits = fits.open(raw_name)
        for j in range(len(index)):
            try:
                PARAM[j].append(raw_fits[ext].header[index[j]])
            except KeyError:
                print("KeyError: Header keyword not found")
                PARAM[j].append("NONE")
    return(list(PARAM))


#######################################

############################

def find_subposition_EXTNAME(EXTNAME):
    import re
    splitted_string = re.findall(r"[-+]?(?:\d*\.*\d+)", EXTNAME)
    subpos1 = splitted_string[-2]
    subpos2 = splitted_string[-1]
    return([float(subpos1), float(subpos2)])

############################

def find_SCA_EXTNAME(EXTNAME):
    import re
    splitted_string = re.findall(r"[-+]?(?:\d*\.*\d+)", EXTNAME)
    sca = splitted_string[-3]
    return(int(sca))


############################
def make_ds9_region(ra, dec, output="ds9.reg"):
    # This program generates a DS9 region file for DS9.

    with open(output, 'w') as the_file:
        the_file.write('# Region file format: DS9 version 4.1\n')
        the_file.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        the_file.write('fk5\n')

        counter = 0
        for ra_i, dec_i in zip(ra, dec):
            label = str(counter).zfill(4)
            the_file.write('circle(' + str(ra_i) + ',' + str(dec_i) + ',1")  # color=red width=4 font="helvetica 16 normal roman" text={'+ label +'}\n')
            counter = counter + 1
    return(output)
############################

############################
def check_file_integrity(filenames):
    # This program ensures that the downloaded files are OK.
    was_any_file_removed = 0

    for filename in filenames:
        try:
            filesize = os.stat(filename).st_size
            if filesize > 0:
                print(f'{filename} has a size greater than 0 KB')
            else:
                print(f'{filename} has a size of 0 KB')
                execute_cmd('rm ' + filename)
                was_any_file_removed = 1

        except:
            print(rs.plots.style.YELLOW + "WARNING: File " + filename + " not found." + rs.plots.style.RESET)
            was_any_file_removed = 1

    if was_any_file_removed:
        print("CAUTION!!!! At least one file was removed")
        print("Please, re download the files, and check their integrity.")
        print("If this error persists, check the download query location.")
        print("or check the disk space. Downloads might be failing.")

    else:
        print("All files verified. Good to go.")
############################

############################

def save_fits(array, name, header=None, extname=None, overwrite=True, output_verify='silentfix'):

    if isinstance(array, (np.ndarray)):
        hdu = fits.PrimaryHDU(header=header, data=array)
        hdul = fits.HDUList([hdu])
        hdul.writeto(name, overwrite=overwrite, output_verify=output_verify)

    if isinstance(array, (list)):
        #print("Multiextension fits")
        hdu1 = fits.PrimaryHDU()
        ext_list = [hdu1]
        #print(ext_list)

        for i in range(len(array)):
            #print("ext " + str(i))
            if isinstance(header, (list)):
                header_ext = header[i]
            else:
                header_ext = None

            ext_list.append(fits.ImageHDU(data = array[i], header = header_ext))

        new_hdul = fits.HDUList(ext_list)
        new_hdul.writeto(name, overwrite=overwrite, output_verify=output_verify)


    if extname is not None:
        output_fits = fits.open(name)

        if isinstance(extname, (list, np.ndarray)):
            for i in range(len(extname)):
                output_fits[i+1].header["EXTNAME"] = extname[i]

        if isinstance(extname, (str)):
            output_fits[0].header["EXTNAME"] = extname

        output_fits.verify("silentfix")
        output_fits.writeto(name, overwrite=True)
        output_fits.close()
    return(name)
############################

############################
def modify_keyword_hdr(input_name, ext, keyword, mode, verbose=True):
    # This program modifies header keywords using Gnuastro astfits.
    # mode = Options ["update", "delete"]
    cmd = "astfits -h" + str(ext) + " " + input_name +  " --" + mode + "=" + keyword
    if verbose == False:
        cmd = cmd + " >/dev/null 2>&1"
    execute_cmd(cmd)

###############################

def flambda_to_fnu(flambda, wavelength):
    # https://en.wikipedia.org/wiki/AB_magnitude#Expression_in_terms_of_f%CE%BB
    # Input:
    # flambda: In erg cm-2 s-1 A-1
    # wavelength: in A
    #
    # Output: Jy
    fnu = (flambda*u.erg / u.cm**2 / u.s / u.AA).to(u.Jy, equivalencies=u.spectral_density(wavelength * u.AA)).value # (3.34E+4)*(wavelength)**2*flambda
    return(fnu)

def fnu_to_flambda(fnu, wavelength):
    # https://en.wikipedia.org/wiki/AB_magnitude#Expression_in_terms_of_f%CE%BB
    # Input:
    # flambda: In Jy
    # wavelength: in A
    flambda = (fnu*u.Jy).to(u.erg / u.cm**2 / u.s / u.AA, equivalencies=u.spectral_density(wavelength * u.AA)).value
    return(flambda)

########################################


def MJysr_to_jyarcsec2(flux_mjy_sr):
    return(flux_mjy_sr*(10**6)/(4.25e10))



###############################################

def get_pixscale(fits_name, ext):
    # We look for the pixsize in the ext 0, if it doesnt work, go to ext 1.
    input_fits = fits.open(fits_name)

    try:
        pixsize = np.abs(input_fits[ext].header["CDELT2"])*60*60
    except:
        pixsize = np.abs(input_fits[ext].header["CD2_2"])*60*60

    return(pixsize)

#################################

def radec_to_xy(ra, dec, fits_name, ext):
    input_fits = fits.open(fits_name)
    w = astropy_wcs.WCS(header=input_fits[ext].header, fobj=input_fits, naxis=2)
    # coord = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    xcen, ycen = w.wcs_world2pix(ra, dec, 0)
    # pixscale = get_pixscale(image_input)  # arcsec/pix
    # print("WCS available: xcen: " + str(xcen) + " ycen: " + str(ycen))
    return([xcen, ycen])


def xy_to_radec(x, y, fits_name, ext):
    input_fits = fits.open(fits_name)
    w = astropy_wcs.WCS(header=input_fits[ext].header, fobj=input_fits, naxis=2)
    ra, dec = w.wcs_pix2world(x, y, 0)
    return([ra, dec])


def get_pixscale(fits_name, ext):
    # We look for the pixsize in the ext 0, if it doesnt work, go to ext 1.
    input_fits = fits.open(fits_name)

    try:
        pixsize = input_fits[ext].header["IDCSCALE"]
    except:
        pixsize = np.abs(input_fits[ext].header["CD2_2"])*60*60
    #except:
    #    pixsize = np.abs(input_fits[ext].header["CDELT2"])*60*60

    return(pixsize)

###########################################

def create_radial_mask(xsize, ysize, q=1, theta=0, center=None, radius=None):
    if center is None: # use the middle of the image
        center = (int(xsize/2), int(ysize/2))

    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], xsize-center[0], ysize-center[1])
    X, Y = np.ogrid[:xsize, :ysize]

    np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    x_image = np.linspace(np.min(X - center[0]),np.max(X - center[0]), xsize)
    y_image = np.linspace(np.min(Y - center[1]),np.max(Y - center[1]), ysize)

    X_image, Y_image = np.meshgrid(x_image,y_image) # Images with the observer X and Y coordinates

    X_gal = (X_image*np.cos(np.radians(-theta))-Y_image*np.sin(np.radians(-theta)))/q # X in galactic frame
    Y_gal = (X_image*np.sin(np.radians(-theta))+Y_image*np.cos(np.radians(-theta)))   # Y in galactic frame

    radial_array = np.sqrt(X_gal**2 + Y_gal**2)
    return(np.flip(radial_array,axis=1))
    #return(radial_array)

###########################################

def create_angle_mask(xsize, ysize, q, theta, center=None, radius=None, pitch=90):
    if center is None: # use the middle of the image
        center = [int(xsize/2), int(ysize/2)]

    center[0] = xsize - center[0]

    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], xsize-center[0], ysize-center[1])
    X, Y = np.ogrid[:xsize, :ysize]

    radial_grid = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    x_image = np.linspace(np.min(X - center[0]),np.max(X - center[0]), xsize)
    y_image = np.linspace(np.min(Y - center[1]),np.max(Y - center[1]), ysize)

    X_image, Y_image = np.meshgrid(x_image,y_image) # Images with the observer X and Y coordinates

    X_gal = (X_image*np.cos(np.radians(theta))-Y_image*np.sin(np.radians(theta)))/q # X in galactic frame
    Y_gal = (X_image*np.sin(np.radians(theta))+Y_image*np.cos(np.radians(theta)))   # Y in galactic frame

    # Calculate pitch angle in the galactic plane coordinate frame (Radial direction-90 + pitch)
    rotation = pitch - 90
    U_gal_pitch = X_gal*np.cos(np.radians(rotation)) - Y_gal*np.sin(np.radians(rotation))
    V_gal_pitch = X_gal*np.sin(np.radians(rotation)) + Y_gal*np.cos(np.radians(rotation))

    # Blot back the angles to the observer frame
    U_ima_pitch = (U_gal_pitch*q*np.cos(np.radians(-theta))-V_gal_pitch*np.sin(np.radians(-theta))) #
    V_ima_pitch = (U_gal_pitch*q*np.sin(np.radians(-theta))+V_gal_pitch*np.cos(np.radians(-theta))) #

    #U_ima = (U_gal_pitch*q*np.cos(np.radians(theta))-V_gal_pitch*np.sin(np.radians(theta)))
    #V_ima = (U_gal_pitch*q*np.sin(np.radians(theta))+V_gal_pitch*np.cos(np.radians(theta)))
    angle_array = np.degrees(np.arctan2(U_ima_pitch,V_ima_pitch))
    # angle_array[angle_array < 0] = angle_array[angle_array < 0] + 360
    return(np.flip(angle_array,axis=1))

#####################################


def make_profile(image, radial_mask=None, ext=None, ra_cen=None, dec_cen=None, xcen=None, ycen=None, rbins=None, pixscale = 1, nbins=100, nsimul=100, q=1, theta=0):

    if isinstance(image, str):
        print("Input image:" + image)
        image_fits = fits.open(image)
        image_name = image
        image = image_fits[ext].data
        header = image_fits[ext].header

        try:
            w = astropy_wcs.WCS(header)
            pixscale = astropy_wcs.utils.proj_plane_pixel_scales(w)[0]
            xcen, ycen = radec_to_xy(ra=ra_cen, dec=dec_cen, fits_name=image_name, ext=ext)
            print("RA" + str(ra_cen) + " - DEC: " + str(dec_cen))
        except:
            print("No valid WCS found. Pixscale set to 1")
            print("Using xcen and ycen supplied by user")
            pixscale = 1


        # Find the coordinates of the target in the image
        if ext is None:
            print("WARNING: No ext was specified")
        print("X" + str(xcen) + " - Y: " + str(ycen))


    shape_image = image.shape

    if xcen is None:
        xcen = int(shape_image[0]/2)

    if ycen is None:
        ycen = int(shape_image[1]/2)


    if radial_mask is None:
        radial_mask = create_radial_mask(xsize=shape_image[1], ysize=shape_image[0],
                                               q=q, theta=theta, center=[xcen,ycen], radius=None)

    if rbins is None:
        print("No rbins provided.")
        min_rbins = np.min(radial_mask)
        max_rbins = np.max(radial_mask)
        print("Assuming " + str(min_rbins) + " - " + str(max_rbins))
        rbins = np.linspace(min_rbins, max_rbins, nbins)

    profile = np.zeros((len(rbins)-1,5))
    rad = np.zeros((len(rbins)-1,3))

    # Make a model at the same time
    model = np.zeros(image.shape)

    for i in tqdm(range(len(rbins)-1)):
        pixels_to_bin = np.where((radial_mask >= rbins[i]) & (radial_mask < rbins[i+1]) & (~np.isnan(radial_mask)))
        if len(pixels_to_bin[0]) > 0:
            rad[i,0] = np.median(radial_mask[pixels_to_bin])*pixscale
            rad[i,1] = np.max(radial_mask[pixels_to_bin])*pixscale
            rad[i,2] = np.min(radial_mask[pixels_to_bin])*pixscale
            # boot_g = bm.bootmedian(image[pixels_to_bin], nsimul=nsimul)
            profile[i,0] = bn.nanmedian(image[pixels_to_bin]) # boot_g["median"]
            profile[i,1] = profile[i,0]+bn.nanstd(image[pixels_to_bin]) # boot_g["s1_up"]
            profile[i,2] = profile[i,0]-bn.nanstd(image[pixels_to_bin]) #boot_g["s1_down"]
            profile[i,3] = profile[i,0]+2*bn.nanstd(image[pixels_to_bin]) #boot_g["s2_up"]
            profile[i,4] = profile[i,0]-2*bn.nanstd(image[pixels_to_bin]) #boot_g["s2_down"]

        else:
            profile[i,0] = np.nan
            profile[i,1] = np.nan
            profile[i,2] = np.nan
            profile[i,3] = np.nan
            profile[i,4] = np.nan

        model[pixels_to_bin] = profile[i,0]


    df = pd.DataFrame(data=np.array([rad[:,0],
                                     rad[:,1],
                                     rad[:,2],
                                     profile[:,0],
                                     profile[:,1],
                                     profile[:,2],
                                     profile[:,3],
                                     profile[:,4]]).T,
                      columns=["r", "r_s1up", "r_s1down", "int", "int_s1up", "int_s1down", "int_s2up", "int_s2down"])


    return({"profile": df, "model": model})


###########################################
def run_basic_noisechisel(file_name, ext):
    os.system("astnoisechisel --tilesize=5,5 -K -h" + str(ext) + " " + file_name)
    return(file_name.replace(".fits", "_detected.fits"))


def run_basic_astrodrizzle(file_name):
    was_a_single_image = 0
    if not isinstance(file_name, (list, pd.core.series.Series, np.ndarray)):
        file_name = np.array([file_name])
        was_a_single_image = 1

    for file_name_i in file_name:
        rs.mast.heal_ACS_DQ_type(input_list=[file_name_i])
        fits_file = fits.open(file_name_i)
        #print(fits_file[0].header)
        print(fits_file[0].header["INSTRUME"])
        #print(fits_file[0].header["DETECTOR"])
        print(fits_file[0].header["FILTER1"])
        print(fits_file[0].header["FILTER2"])
        print(fits_file[0].header["RA_TARG"])
        print(fits_file[0].header["DEC_TARG"])

    #xy_to_radec(x, y, fits_name, ext)
    from drizzlepac import astrodrizzle

    astrodrizzle.AstroDrizzle(input=file_name, static=False,
                          skymethod='globalmin+match', skysub=0,
                          driz_sep_bits='64, 32', final_rot=0,
                          combine_type='imedian', updatewcs=False,
                          driz_combine=1)

    if was_a_single_image:
        if (file_name[-9:0]=="_flt.fits"):
            drz_name = os.path.dirname(file_name) + "/" + os.path.basename(file_name).split("_")[0] + "_drz_sci.fits"

        if (file_name[-9:0]=="_flc.fits"):
            drz_name = os.path.dirname(file_name) + "/" + os.path.basename(file_name).split("_")[0] + "_drc_sci.fits"

        else:
            drz_name = os.path.dirname(file_name) + "/" + os.path.basename(file_name)[:-5] + "_drz_sci.fits"
    else:
        drz_name = os.path.dirname(file_name) + "/" + "final_drz_sci.fits"

    return(drz_name)

#################

def mask_sources(file_name, ext):
    input_fits = fits.open(file_name)
    detected_name = run_basic_noisechisel(file_name, ext)
    detected_fits = fits.open(detected_name)
    masked_image = np.copy(input_fits[ext].data)
    masked_image[detected_fits["DETECTIONS"].data==1] = np.nan
    output_name = file_name.replace(".fits","_masked.fits")
    save_fits(array=masked_image, name=output_name, header=input_fits[ext].header, extname=["MASKED"])
    return(output_name)

###########################################


def create_custom_wcs(crpix, crval, cdelt, crota=[0,0], projection="TAN"):
   # WCS code ripped from
   # http://docs.astropy.org/en/latest/wcs/index.html

    w = astropy_wcs.WCS(naxis=2)

    # What is the center pixel of the XY grid.
    w.wcs.crpix = crpix

    # what is the galactic coordinate of that pixel.
    w.wcs.crval = crval

    # what is the pixel scale in lon, lat.
    # This essentially modifies how apparently large is the galaxy. Tweak it at your desire.
    w.wcs.cdelt = cdelt
    w.wcs.crota = crota

    # you would have to determine if this is in fact a tangential projection.
    w.wcs.ctype = ["RA---"+projection, "DEC--"+projection]

    # write the HDU object WITH THE HEADER
    header = w.to_header()
    return(header)

def create_dummy_image_with_wcs(ra_cen, dec_cen, ra_size, dec_size, pixscale, outname="default_dummy.fits"):
    # create_dummy_image_with_wcs.
    #
    # Input:
    # ra_cen, dec_cen: float (degrees): Central coordinates of the image.
    # ra_range, dec_range: tuple, 2 float (degrees): Min and max values of ra, dec. They need to be in order.
    # pixscale: float (degrees/px): Pixel scale in degree per pixel.
    # outname: string: Name of the output image.


    # Calculate how large should be the image
    x_size = int(np.ceil(1.01*ra_size/pixscale))
    y_size = int(np.ceil(1.01*dec_size/pixscale))

    # Make the dummy image
    dummy_array = np.zeros((x_size, y_size))

    header = create_custom_wcs(crpix=[x_size/2, y_size/2], crval=[ra_cen, dec_cen], cdelt=[pixscale, pixscale])

    rs.utils.save_fits(array=dummy_array, name=outname, header=header)

    return(outname)

def detect_sci_extensions(input_name):
    # detect_sci_extensions: It finds which extensions have the SCI id in
    # the FITS header.
    # Input: input_name. String. The name of the fits file to analyze.

    input_fits = fits.open(input_name)

    sci_ext_list = []
    for ext_i in range(len(input_fits)):
        try:
            if input_fits[ext_i].header["EXTNAME"] == "SCI":
                sci_ext_list.append(ext_i)
        except:
            continue

    return(sci_ext_list)

#########################
def find_nearest_index(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return(idx)


# Set the WCS information manually by setting properties of the WCS
# object.

def create_dummy_exposure(telescope, instrument, detector, exposure_params, binning=1, dummy_name="dummy.fits"):

    # Find telescope class
    telescope_class = rs.telescopes.telescope_class_finder(telescope)

    # Create a new WCS object.  The number of axes must be set from the start
    w = astropy_wcs.WCS(naxis=2)

    # Setting some dummy image
    dummy_shape = telescope_class.get_canvas_shape(instrument=instrument)
    pixscale = telescope_class.get_pixscale(instrument=instrument)
    image_data = np.random.normal(loc=0, scale=1, size=[int(dummy_shape[0]/binning),int(dummy_shape[1]/binning)])

    # Vector properties may be set with Python lists, or Numpy arrays
    w.wcs.crpix = [dummy_shape[1]/binning/2., dummy_shape[0]/binning/2.]
    w.wcs.cdelt = np.array([-binning*pixscale.to("arcsec").value/3600., binning*pixscale.to("arcsec").value/3600.])
    w.wcs.crval = [exposure_params["ra"], exposure_params["dec"]]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    #w.wcs.set_pv([(2, 1, exposure_params["pa"])])
    #w.wcs.crota([30, 30])

    # Now, write out the WCS object as a FITS header
    header = w.to_header()
    header["CROTA1"] = -exposure_params["pa"]
    header["CROTA2"] = -exposure_params["pa"]
    header["EXTNAME"] = "SCI"
    header["RA_TARG"] = exposure_params["ra"]
    header["DEC_TARG"] = exposure_params["dec"]
    header["EXPSTART"] = exposure_params["mjd_start"]
    header["EXPEND"] = exposure_params["mjd_end"]
    header["PA"] = exposure_params["pa"]
    header["TELESCOP"] = telescope
    header["INSTRUME"] = instrument
    header["DETECTOR"] = detector


    # header is an astropy.io.fits.Header object.  We can use it to create a new
    # PrimaryHDU and write it to a file.
    hdu = fits.PrimaryHDU(header=header)
    save_fits(image_data, dummy_name, header)
    return(dummy_name)

########################################################

def delta_angular_separation(ra, dec):
    npoints = len(ra)
    coords = SkyCoord(ra, dec, unit="degree", frame="icrs")
    delta_separation = np.zeros((npoints))

    for i in range(npoints-1):
        #print(coords[i+1].separation(coords[i]).value)
        delta_separation[i+1] = coords[i+1].separation(coords[i]).value
    delta_separation[0] = delta_separation[1]
    return(delta_separation)

######################################################

def great_circle_ra_dec_shift(theta, phi):
    """
    Estimates the new coordinates (ra2, dec2)
    after moving an angular distance theta (radial direction) in the phi direction
    (0 deg, North, counterclockwise) along a Great Circle from the Equator (ra=0, dec=0)
    Adapted from Great Circle navigation: https://en.wikipedia.org/wiki/Great-circle_navigation

    :param theta:
    :type theta: :class:`np.ndarray` or float in degrees
    :param phi:
    :type phi: :class:`np.ndarray` or float in degrees
    :return: radec_2 :dict: Contains the new coordinates "ra", "dec" in degrees.
    """
    import numpy as np
    y1 = np.cos(-np.radians(phi))*np.sin(np.radians(theta))
    y2 = np.sqrt( (np.cos(np.radians(theta)))**2 + ((np.sin(-np.radians(phi)))**2 * (np.sin(np.radians(theta))**2) ))
    dec_2 = 180*np.arctan2(y1, y2)/np.pi

    y1 = np.sin(-np.radians(phi))*np.sin(np.radians(theta))
    y2 = np.cos(np.radians(theta))
    ra_2 = 180*np.arctan2(y1, y2)/np.pi

    print(ra_2, dec_2, np.sqrt(ra_2**2 + dec_2**2))

    return({"ra": ra_2, "dec": dec_2})

def sphere_dist(ra1, dec1, ra2, dec2):
    """
    A faster separation calculator than astropy solution.
    # https://stackoverflow.com/questions/71506875/faster-alternative-to-skycoord-separation

    Haversine formula for angular distance on a sphere: more stable at poles.
    This version uses arctan instead of arcsin and thus does better with sign
    conventions.  This uses numexpr to speed expression evaluation by a factor
    of 2 to 3.

    :param ra1: first RA (deg)
    :param dec1: first Dec (deg)
    :param ra2: second RA (deg)
    :param dec2: second Dec (deg)

    :returns: angular separation distance (deg)
    """

    ra1 = np.radians(ra1).astype(np.float64)
    ra2 = np.radians(ra2).astype(np.float64)
    dec1 = np.radians(dec1).astype(np.float64)
    dec2 = np.radians(dec2).astype(np.float64)
    import numexpr
    numerator = numexpr.evaluate('sin((dec2 - dec1) / 2) ** 2 + '
                                 'cos(dec1) * cos(dec2) * sin((ra2 - ra1) / 2) ** 2')

    dists = numexpr.evaluate('2 * arctan2(numerator ** 0.5, (1 - numerator) ** 0.5)')
    return np.degrees(dists)

#########################
def angular_distance(ra1, dec1, ra2, dec2):
    """
    Einsum version of the angular separation function by Astropy einsum_angular_distance_harversine
    v1: October 17, 2024 - Fastest version compared against numba, or astropy.

    """
    import numexpr

    pi = np.pi
    ra1_rad  = numexpr.evaluate('ra1*pi/180')
    ra2_rad  = numexpr.evaluate('ra2*pi/180')
    dec1_rad = numexpr.evaluate('dec1*pi/180')
    dec2_rad = numexpr.evaluate('dec2*pi/180')
    cdec2    = numexpr.evaluate('cos(dec2_rad)')
    cdec1    = numexpr.evaluate('cos(dec1_rad)')

    ra2_m_ra1     = np.add.outer(-ra1_rad, ra2_rad)
    dec2_m_dec1   = np.add.outer(-dec1_rad, dec2_rad)
    cdec1_x_cdec2 = np.einsum('i,j', cdec1, cdec2, optimize='optimal')
    distance      = numexpr.evaluate('180/pi*2*arcsin(sqrt((1 - cos(dec2_m_dec1) + cdec1_x_cdec2 * (1 - cos(ra2_m_ra1)) )/2))')
    return(distance)


######### POSITION ANGLE ##############
def position_angle(ra1, dec1, ra2, dec2):
    """
    Einsum version of the position angle function by Astropy einsum_angular_distance_harversine
    v1: October 17, 2024
    """
    #radec1_to_radians = np.radians([ra1, dec1])
    #radec2_to_radians = np.radians([ra2, dec2])
    pi = np.pi
    import numexpr

    ra1_rad = numexpr.evaluate('ra1*pi/180')
    ra2_rad = numexpr.evaluate('ra2*pi/180')
    dec1_rad = numexpr.evaluate('dec1*pi/180')
    dec2_rad = numexpr.evaluate('dec2*pi/180')

    cdec2 = numexpr.evaluate('cos(dec2_rad)')
    sdec2 = numexpr.evaluate('sin(dec2_rad)')
    cdec1 = numexpr.evaluate('cos(dec1_rad)')
    sdec1 = numexpr.evaluate('sin(dec1_rad)')

    deltalon   = np.add.outer(-ra1_rad, ra2_rad)
    colat = np.einsum('i,j', np.ones(len(dec1)), cdec2, optimize='optimal')
    sdec2_x_cdec1 = np.einsum('i,j', cdec1, sdec2, optimize='optimal')
    cdec2_x_sdec1 = np.einsum('i,j', sdec1, cdec2, optimize='optimal')

    position_angle = numexpr.evaluate('180/pi*arctan2(sin(deltalon) * colat, sdec2_x_cdec1 - cdec2_x_sdec1 * cos(deltalon))')
    return(position_angle)

def separation_and_position_angle(ra1, dec1, ra2, dec2):

    """
    Einsum version of the angular separation function by Astropy einsum_angular_distance_harversine
    v1: October 17, 2024 - Fastest version compared against numba, or astropy.

    """
    pi = np.pi
    import numexpr
    ra1_rad = numexpr.evaluate('ra1*pi/180')
    ra2_rad = numexpr.evaluate('ra2*pi/180')
    dec1_rad = numexpr.evaluate('dec1*pi/180')
    dec2_rad = numexpr.evaluate('dec2*pi/180')
    cdec2 = numexpr.evaluate('cos(dec2_rad)')
    cdec1 = numexpr.evaluate('cos(dec1_rad)')
    sdec2 = numexpr.evaluate('sin(dec2_rad)')
    sdec1 = numexpr.evaluate('sin(dec1_rad)')

    #if len(ra1.shape) > 0:
    if not isinstance(ra1, (np.floating, float)):
        deltalon   = np.add.outer(-ra1_rad, ra2_rad)
        deltalat = np.add.outer(-dec1_rad, dec2_rad)
        colat = np.einsum('i,j', np.ones(len(dec1)), cdec2, optimize='optimal')
        sdec2_x_cdec1 = np.einsum('i,j', cdec1, sdec2, optimize='optimal')
        cdec2_x_sdec1 = np.einsum('i,j', sdec1, cdec2, optimize='optimal')
        cdec1_x_cdec2 = np.einsum('i,j', cdec1, cdec2, optimize='optimal')
    else:
        deltalon   = ra2_rad-ra1_rad
        deltalat = dec2_rad-dec1_rad
        colat = cdec2
        sdec2_x_cdec1 = cdec1*sdec2
        cdec2_x_sdec1 = sdec1*cdec2
        cdec1_x_cdec2 = cdec1*cdec2

    distance = numexpr.evaluate('180/pi*2*arcsin(sqrt((1 - cos(deltalat) + cdec1_x_cdec2 * (1 - cos(deltalon)) )/2))')
    position_angle = numexpr.evaluate('180/pi*arctan2(sin(deltalon) * colat, sdec2_x_cdec1 - cdec2_x_sdec1 * cos(deltalon))')
    return({"distance": distance, "position_angle": position_angle})


####################
def circle_around_position(ra, dec, radius, npoints=200):
    import astropy.units as u
    from astropy.coordinates import ICRS, Angle, SkyCoord
    central_coords = SkyCoord(ra, dec, frame="icrs", unit="deg")

    position_angle = np.linspace(0, 360, npoints)*u.degree
    circle = central_coords.directional_offset_by(position_angle, separation=radius*u.degree)
    return(circle)



####################
def thetaphi_2_radec(theta, phi):
    #phi = np.array(phi)
    #theta = np.array(theta)
    pi = np.pi
    phi_radians = numexpr.evaluate('phi*pi/180')
    theta_radians = numexpr.evaluate('theta*pi/180')
    dec_in_NDI = numexpr.evaluate('180/pi*arctan2(cos(-phi_radians)*sin(theta_radians),  sqrt((cos(theta_radians))**2 + ((sin(-phi_radians))**2 * (sin(theta_radians)**2) )))')
    ra_in_NDI = numexpr.evaluate('-180/pi*arctan2(sin(-phi_radians)*sin(theta_radians), cos(theta_radians))')

    return([ra_in_NDI, dec_in_NDI])


####################

def generate_image_interpolator(data, x=None, y=None):
    from scipy.interpolate import LinearNDInterpolator

    if x is None:
        x = np.linspace(0, data.shape[0], data.shape[0])-data.shape[0]/2.
    if y is None:
        y = np.linspace(0, data.shape[1], data.shape[1])-data.shape[0]/2.

    XX, YY = np.meshgrid(x, y)
    XY_points = np.c_[XX.flatten(), YY.flatten()]
    Z_points = data.flatten()

    image_interpolator = LinearNDInterpolator(XY_points, Z_points, fill_value=0)
    return(image_interpolator)


########################

def convert_ASDF_to_FITS(asdf_list, output):
    # We sort alphabetically the list.
    asdf_list.sort()

    array_list = []
    header_list = []
    extname_list = []

    for asdf_name in tqdm(asdf_list):
        exposure_identity = rs.utils.exposure_inspector(asdf_name)
        array_list.append(exposure_identity["DATA"][0])
        header = exposure_identity["ASTROPYWCS"][0].to_header()

        keywords_to_paste = ['INSTRUME','TELESCOP','DETECTOR','FILTER',
                             'RA_TARG','DEC_TARG', #'SUNANGLE','MOONANGL','EARTHANGL',
                             'EXPSTART','EXPEND','EXPTIME','conversion_megajanskys',
                             'pixel_area','EXPSTART_ISOT','SCA',
                             'RA_PNT','DEC_PNT', 'PA','FILETYPE']

        for keyword_to_paste in keywords_to_paste:
            #print(keyword_to_paste)
            #print(exposure_identity[keyword_to_paste])
            if isinstance(exposure_identity[keyword_to_paste], (float, int, str, bool,)):
                header[keyword_to_paste] = exposure_identity[keyword_to_paste]
            else:
                header[keyword_to_paste] = (exposure_identity[keyword_to_paste].value, str(exposure_identity[keyword_to_paste].unit))

        header["FILTER1"] = exposure_identity["FILTER"]
        header_list.append(header)
        extname_list.append("SCI")

    rs.utils.save_fits(array=array_list, name=output, header=header_list, extname=extname_list)

    # Adding some required keywords in the ext 0
    output_fits = fits.open(output)
    output_fits[0].header["FILTER"] = output_fits[1].header["FILTER"]
    output_fits[0].header["FILTER1"]= output_fits[1].header["FILTER1"]
    output_fits[0].header["FILTER2"]= "CLEAR"
    output_fits[0].header["INSTRUME"] = output_fits[1].header["INSTRUME"]
    output_fits[0].header["TELESCOP"]= output_fits[1].header["TELESCOP"]
    output_fits[0].header["EXPSTART"]= output_fits[1].header["EXPSTART"]
    output_fits[0].header["EXPEND"]= output_fits[1].header["EXPEND"]
    output_fits[0].header["EXPTIME"]= output_fits[1].header["EXPTIME"]

    if "ROMAN" in output_fits[1].header["TELESCOP"]:
        output_fits[0].header["DETECTOR"] = "WFI"
    else:
        output_fits[0].header["DETECTOR"] = output_fits[1].header["DETECTOR"]

    output_fits[0].header["RA_PNT"]= output_fits[1].header["RA_PNT"]
    output_fits[0].header["DEC_PNT"]= output_fits[1].header["DEC_PNT"]
    output_fits[0].header["RA_TARG"]= output_fits[1].header["RA_TARG"]
    output_fits[0].header["DEC_TARG"]= output_fits[1].header["DEC_TARG"]
    output_fits[0].header["PA"]= output_fits[1].header["PA"]
    output_fits.verify("silentfix")
    output_fits.writeto(output, overwrite=True)
    return(output)

def read_ds9reg(fname):
    #######################s

    # A program to transform ds9 polygon vertice files to fits masks
    # Alejandro Serrano Borlaff
    # PhD Student - Instituto de Astrofisica de Canarias
    # v.1.0 - First version (31/07/2017)
    #####################################################################

    # The main idea was taken from
    # http://xingxinghuang.blogspot.fi/2014/05/imagebad-pixel-mask-using-ds9-region.html
    # And https://stackoverflow.com/questions/3654289/scipy-create-2d-polygon-mask

    with open(fname) as f:
        content = f.readlines()
    # you may also want to remove whitespace characters like `\n` at the end of each line
    content = [x.strip() for x in content]
    content = content[3:]
    polygons = []
    for i in range(len(content)):
        polygon_unit = []
        polygon = content[i].replace("polygon(","").replace(")","").split(",")
        #print(len(polygon))
        #print(np.array(np.linspace(start=0,stop=len(polygon)-2,num=len(polygon)/2),dtype="int"))
        for j in np.array(np.linspace(start=0,stop=len(polygon)-2,num=int(len(polygon)/2)),dtype="int"):
            x = int(np.round(np.array(float(polygon[j]))))
            y = int(np.round(np.array(float(polygon[j+1]))))
            polygon_unit.append((x,y))
        polygons.append(polygon_unit)
    return(polygons)

#####################################################################

def ds9tomask(fname, nx, ny, outname):
    polygons = read_ds9reg(fname)
    #nx, ny = 2000, 2000
    mask = np.zeros(shape=(ny,nx))
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))
    x, y = x.flatten(), y.flatten()
    points = np.vstack((x,y)).T

    from matplotlib.path import Path

    for i in tqdm(range(len(polygons))):
        path = Path(polygons[i])
        grid = path.contains_points(points)
        grid = grid.reshape((ny,nx))
        mask[grid] = 1

    hdu = fits.PrimaryHDU(mask)
    stout = execute_cmd("rm "+outname)
    hdu.writeto(outname)
    print("File saved:"+ outname)
    print("Note: 1 is masked. 0 is not masked")
    return(outname)

#####################################################################
def erwinspacing(minimum, number, rate=1.03, mode="int"):

    # A tool to make log-like spacing arrays.
    # Inspired on Peter Erwin's
    # THE OUTER DISKS OF EARLY-TYPE GALAXIES. I. SURFACE-BRIGHTNESS PROFILES OF BARRED GALAXIES
    # https://iopscience.iop.org/article/10.1088/0004-6256/135/1/20/pdf
    # Erwin et al. 2008. ApJ
    ###########################################
    bins = []
    bins.append(minimum)

    last = minimum
    for i in range(number-1):
        new = rate*last

        if mode=="int":
            new = np.round(new,0)
            if last == new:
                new = new + 1

        bins.append(new)

        last = new

    if mode=="int":
        bins = np.array(bins, dtype="int")
    else:
        bins = np.array(bins)

    return(bins)

#####################################################################

def get_data_and_wcs(input_name, ext):
    """
    Opens a fits file and returns the data and wcs at a certain extension.

    :param input_name:
    :type input_name: :class:`str`
    :param ext:
    :type ext: :class:`int`
    :return: :list: A two element list. [data (numpy.ndarray), wcs (astropy.wcs.wcs.WCS)]
    """

    input_fits = fits.open(input_name)
    wcs = astropy_wcs.WCS(input_fits[ext].header, input_fits)
    return([input_fits[ext].data, wcs])

#####################################################################

def get_keys_from_header(fits_list, index, ext=0):
    PARAM = []
    for j in range(len(index)):
        PARAM.append([])
    for raw_name in fits_list:
        # print(raw_name)
        raw_fits = fits.open(raw_name)
        for j in range(len(index)):
            try:
                PARAM[j].append(raw_fits[ext].header[index[j]])
            except KeyError:
                print("KeyError: Header keyword not found")
                PARAM[j].append("NONE")
    return(list(PARAM))

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
#####################################################################

def find_max_angular_size_of_image(data, wcs):
    """
    Returns the maximum angular extension of an image.

    :param data:
    :type data: :class:`numpy.ndarray`
    :param wcs:
    :type wcs: :class:`astropy.wcs.wcs.WCS`
    :return: :float: The maximum angular extension of the image in sky coordinates in degrees.
    """

    data_shape = wcs.array_shape
    ra_cen, dec_cen = wcs.wcs_pix2world(data_shape[0]/2, data_shape[1]/2, 0)
    corners = rs.detectors.get_detector_corners(data=data, wcs=wcs)
    distance = rs.utils.sphere_dist(ra1=ra_cen, dec1=dec_cen, ra2=corners["corners_world"][:,0], dec2=corners["corners_world"][:,1])
    return(np.nanmax(distance))



def check_fits_integrity(input):
    rs.utils.execute_cmd("mkdir CORRUPT")
    if isinstance(input, (str,)):
        input_list = [input]
    else:
        input_list = input

    clean_list = []
    for input_name in tqdm(input_list):
        try:
            dummy_open = fits.open(input_name)
            dummy_open.close()
            clean_list.append(input_name)
        except:
            print("FITS file " + input_name + " seems to be corrupted. Moving to CORRUPT folder")
            rs.utils.execute_cmd("mv " + input_name + " CORRUPT/")

    return(clean_list)


def circular_hist(ax, x, bins=16, density=True, offset=0, gaps=True, color="dodgerblue", edgecolor='C0'):
    """
    Produce a circular histogram of angles on ax.

    Parameters
    ----------
    ax : matplotlib.axes._subplots.PolarAxesSubplot
        axis instance created with subplot_kw=dict(projection='polar').

    x : array
        Angles to plot, expected in units of radians.

    bins : int, optional
        Defines the number of equal-width bins in the range. The default is 16.

    density : bool, optional
        If True plot frequency proportional to area. If False plot frequency
        proportional to radius. The default is True.

    offset : float, optional
        Sets the offset for the location of the 0 direction in units of
        radians. The default is 0.

    gaps : bool, optional
        Whether to allow gaps between bins. When gaps = False the bins are
        forced to partition the entire [-pi, pi] range. The default is True.

    Returns
    -------
    n : array or list of arrays
        The number of values in each bin.

    bins : array
        The edges of the bins.

    patches : `.BarContainer` or list of a single `.Polygon`
        Container of individual artists used to create the histogram
        or list of such containers if there are multiple input datasets.
    """
    # Wrap angles to [-pi, pi)
    x = (x+np.pi) % (2*np.pi) - np.pi

    # Force bins to partition entire circle
    if not gaps:
        bins = np.linspace(-np.pi, np.pi, num=bins+1)

    # Bin data and record counts
    n, bins = np.histogram(x, bins=bins)

    # Compute width of each bin
    widths = np.diff(bins)

    # By default plot frequency proportional to area
    if density:
        # Area to assign each bin
        area = n / x.size
        # Calculate corresponding bin radius
        radius = (area/np.pi) ** .5
    # Otherwise plot frequency proportional to radius
    else:
        radius = n

    # Plot data on ax
    patches = ax.bar(bins[:-1], radius, zorder=1, align='edge', width=widths, color=color,
                     edgecolor=edgecolor, fill=True, alpha=0.75, linewidth=1)

    # Set the direction of the zero angle
    ax.set_theta_offset(offset)

    # Remove ylabels for area plots (they are mostly obstructive)
    if density:
        ax.set_yticks([])

    return n, bins, patches


def divide_array_in_chunks(array, chunk_size):
    tail_size = len(array) % chunk_size
    n_chunks = int(np.floor(len(array)/chunk_size) + 1)
    counter = 0
    chunks = []
    start_i = 0
    end_i = 0
    while counter < (n_chunks - 1):
        end_i = end_i + chunk_size
        chunks.append(array[start_i:end_i])
        counter = counter + 1
        start_i = end_i

    if tail_size > 0:
        chunks.append(array[end_i:])

    return(chunks)
