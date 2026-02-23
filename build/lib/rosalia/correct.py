# Alejandro S. Borlaff. NASA Ames Research Center. a.s.borlaff@nasa.gov / asborlaff@gmail.com
# June 27, 2023.
#
# STRAYCOR main module
# This module will hold all the general programs to be interacting with the user
#
# Version log:
# v.1.0 - 27 June 2023. First loading of programs inherited from former straycor.
#
##########################################################
#  Required Packages
#  pip install sep
#
#
##########################################################

# General modules
import os
import glob
import pandas as pd
import numpy as np
import bottleneck as bn
import astropy.units as u
from tqdm import tqdm
from astropy import coordinates
from astroquery.mast import Observations
from astropy.io import fits
import astropy.wcs as wcs
from astropy.time import Time
from scipy import interpolate
import matplotlib.pyplot as plt
from astropy import constants as const
import rosalia as rs

# from rosalia.utils import exposure_inspector
# from rosalia.utils import convert_ASDF_to_FITS

###########################

def rosalia_stray(input_name, output_name="rosalia_stray_output.fits", radius=1, g_mag_max=15, sun_block=False, verbose=False, catalog=None):
    main_offender_db = main_offender(input_name=input_name,
                                     radius=radius,
                                     g_mag_max=g_mag_max,
                                     verbose=verbose,
                                     input_catalog=catalog,
                                     match_inside_stars=False,
                                     sun_block = sun_block,
                                     output_name = output_name)
    return(main_offender_db)

###########################

def rosalia_zody(input_name, output_name="rosalia_zody_output.fits"):
    exposure_identity = rs.utils.exposure_inspector(input_name)
    zodiacal_background_list = []

    for i in tqdm(range(len(exposure_identity["SCIEXTS"]))):
        zodiacal_background = rs.sky.get_zodiacal_background(input_name=input_name,
                                                      ext=exposure_identity["SCIEXTS"][i],
                                                      wavelength=exposure_identity["FILTER"],
                                                      telescope=exposure_identity["TELESCOP"],
                                                      instrument=exposure_identity["INSTRUME"],
                                                      detector=exposure_identity["DETECTOR"],
                                                      expstart=exposure_identity["EXPSTART"],
                                                      step=1000, zody_mode="zodipy",
                                                      nbins_wavelength=20, obslocin=3,
                                                      grid_method="random", output_units=None, verbose=False)
        zodiacal_background_list.append(zodiacal_background)

    ########################################
    # Save the results to a fits file.
    ########################################
    data_output = []
    header_output = []
    for SCIEXT_i, zodiacal_background_i in tqdm(zip(exposure_identity["SCIEXTS"], zodiacal_background_list)):
        data_output.append(zodiacal_background_i)
        header_output.append(exposure_identity["ASTROPYWCS"][SCIEXT_i-1].to_header())

    rs.utils.save_fits(array=data_output, name=output_name, header=header_output,
                       extname=None, overwrite=True, output_verify='silentfix')

    print("Output saved in: " + output_name)

    return({"image_identity":exposure_identity,
            "zodi_list": zodiacal_background_i,
            "output_name": output_name})


###########################


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
    rs.mast.bestrefs(output, clean=True)

    print(output)

    absolute_path_output = []
    for i in tqdm(range(len(output))):
        absolute_path_output.append(os.path.abspath(output[i]))

    print("Done")
    return(absolute_path_output)

############################

def do_astrometry(input_name, refcat=None, refxcol=None, refycol=None, force_tweakreg=False, correct_crs=False):
    output = rs.astrometry.astrometry_list(fitslist=input_name, refcat=refcat, refxcol=refxcol, refycol=refycol,
                                        force_tweakreg=force_tweakreg, correct_crs=correct_crs)
    return(output)

############################

def subtract_stars(input_name, clean=True, verbose=False):

    output_list = []
    for input_name_i in tqdm(input_name):
        # Run HST_inspector
        input_inspection = rs.utils.exposure_inspector(input_name_i)

        # Select the right order of extensions
        if (input_inspection["INSTRUME"] == "ACS"):
            sci_exts = [1, 4]
        if (input_inspection["INSTRUME"] == "WFC3IR"):
            sci_exts = [1]

        # Select the right PSF
        #
        print("WARNING: For testing purposes the PSF is fixed")
        psf_name = os.path.dirname(rs.utils.__file__) + "/../PSF_ARCHIVE/f814w00_arith.fits"
        print("PSF: " + psf_name)

        output = rs.psf.scale_and_subtract_stars(input_name=input_name_i, ext=sci_exts,
                                              psf_name=psf_name, clean=clean, verbose=verbose)

        output_list.append(output[0])

        if clean:
            rs.utils.execute_cmd("rm " + input_name_i.replace(".fits", "*_ext*.fits"))

    return(output_list)


##########################################

def correct_zody(input_name, verbose=False):
    output_list = []
    for input_name_i in tqdm(input_name):
        # Run HST_inspector
        if verbose:
            print(input_name_i)
        input_inspection = rs.utils.exposure_inspector(input_name_i, verbose)

        # Select the right order of extensions
        if (input_inspection["INSTRUME"] == "ACS"):
            corrected_image = rs.sky.remove_zodiacal_light_acs(input_name_i)
        if (input_inspection["INSTRUME"] == "WFC3IR"):
            print("Zodiacal light correction not yet implemented for WFC3!")
            raise(Exception)

        output_list.append(corrected_image)
    return(output_list)




def main_offender(input_name=None, ext=None, ra=None, dec=None, phi=0,
                  filter_name=None, instrument=None, telescope=None, detector=None,
                  exptime=None, expstart=None, radius=1, g_mag_max=False,
                  step=200, grid_method="random", verbose=False, ndi_mode="legacy",
                  input_catalog=None, match_inside_stars=False, sun_block=False,
                  output_name="main_offender_default_output.fits"):

    #######################################
    # main_offender: Alejandro S. Borlaff. NASA/Ames STA. a.s.borlaff@nasa.gov
    # -------------------------------
    # The objective of this program is to identify and estimate the straylight from stars outside the field of view for Optical and NIR observations.
    # As input, the user can provide either a fits file with a WCS, or a pair of coordinates (ra, dec, ICRS).
    # The output is an image with the estimated straylight per pixel if an image is provided.
    # If the input is a coordinate, then the estimation is only performed in the center of the coordinates.
    # --------------------------------
    # History:
    # v1 - 29 Feb 2024. First working version.
    # v2 - 14 June 2024. Adding g_mag_max to limit the number for Gaia stars.
    #                    Adding Multiextension FITS functionality.
    # v3 - 18 June 2024. Reshaping the code to estimate first which stars are inside which detector before straylight estimations.
    # v4 - 22 October 2024. Adapting main_offender to accept Roman ASDF files.
    #      Nov 5 2024 - Roman Single SCA ASDF files now working in main_offender
    #
    #######################################

    # Gathering info about the scene
    if input_name is not None:
        if verbose: print("> Image mode:" + input_name)
        bool_image_mode = True

        # Get the exposure time
        image_identity = rs.utils.exposure_inspector(input_name, verbose)
        exptime = image_identity["EXPTIME"]
        MJD     = image_identity["EXPSTART"]
        SCIEXTS = image_identity["SCIEXTS"]

        if (len(SCIEXTS) == 0) and (ext is None):
            raise(Exception("No SCI extensions in FITS file, and no extensions (ext) defined in input. Please flag the science extensions in the header of the FITS file, or specify them in your call to main_offender."))

        # Find closest stars with Gaia.
        filter_name = image_identity["FILTER"]
        instrument  = image_identity["INSTRUME"]
        telescope   = image_identity["TELESCOP"]
        lambda_ref  = image_identity["FILTER_IDENTITY"]["filter_lambda_ref"]
        RA_PNT      = image_identity["RA_PNT"]
        DEC_PNT     = image_identity["DEC_PNT"]
        PA_PNT      = image_identity["PA"]

    # If the user does not input an image, just a set of coordinates, then bool_image_mode is False
    else:
        if verbose: print("> Single location mode: RA: " + str(ra) + " DEC: " + str(dec))
        bool_image_mode = False
        filter_db  = rs.telescopes.find_filter_in_svo(wavelength=filter_name, telescope=telescope,
                                                      instrument=instrument, detector=detector, verbose=False)
        lambda_ref = filter_db["filter_lambda_ref"]
        MJD        = expstart
        RA_PNT     = ra
        DEC_PNT    = dec
    # Find the central coordinates of the image - If there are more than one science extensions,
    # then find the center of them all. <--- This goes into exposure_inspector


    # Find the stars around the central coordinate of the scene.
    # Save the stellar catalog into a catalog object, and run main_offender as if input_catalog was set by the User.

    #print("Demo warning: ADD FEATURE planets: https://github.com/skyfielders/python-skyfield/tree/master")

    if input_catalog is None:
        if verbose: print("> Querying stars in the surroundings using ESA/Gaia Archive")

        if radius > 0.5:
            print("INFO: radius parameter (minimum distance to search for individual stars) is > 0.5 degrees.")
            print("Gaia/2MASS/WISE query database can take several minutes to process. Please be patient.")

        loader = rs.plots.Loader("Querying Gaia/2MASS/WISE/JPL Horizons databases. This might take a few minutes...",
                                 "All-sky source map constructed.", 0.05).start()
        gaia_query = rs.psf.find_gaia_stars_around_image(lambda_ref=lambda_ref,
                                                         observer=telescope,
                                                         input_name=None,
                                                         ext=None,
                                                         ra=RA_PNT, dec=DEC_PNT, MJD=MJD,
                                                         radius=radius,
                                                         g_mag_max=g_mag_max, verbose=verbose)["gaia_query"]
        loader.stop()
    else:
        gaia_query = input_catalog

    # If sun_block is True, then remove the Sun from the catalog.
    if sun_block:
        gaia_query = gaia_query[~(gaia_query["source_id"] == "Sun")]
        gaia_query = gaia_query[~(gaia_query["source_id"] == "Earth")]

    # Find where each star lands (detector ID or outside FOV)
    names_of_bool_columns_if_star_is_inside = []
    detector_square_list = []

    """
    If the input file is a multi-extension fits, then exposure_inspector will scan for extensions with EXTNAME = SCI.
    The extension ID in the FITS file will be stored in SCIEXTS = image_identity["SCIEXTS"].

    In that case, gaia_query, the catalog of stars, will have a set of N columns called in_SCI[i] (boolean), where the catalog
    stores if that particular star is inside each detector or not.

    """

    for SCIEXT_i in tqdm(SCIEXTS):
        if verbose: print("> Identifying which stars are inside the FOV and which are outside...")
        infield_stars = rs.psf.identify_stars_in_out_field(data=image_identity["DATA"][SCIEXT_i-1],
                                                           wcs=image_identity["ASTROPYWCS"][SCIEXT_i-1],
                                                           catalog=gaia_query,
                                                           verbose=verbose)

        name_column_is_star_inside_this_detector = "in_SCI" + str(SCIEXT_i)
        names_of_bool_columns_if_star_is_inside.append(name_column_is_star_inside_this_detector)

        gaia_query[name_column_is_star_inside_this_detector] = infield_stars["bool_isIn"]

        # If verbose, make a plot of the stars with the footprint.
        detector_corners = rs.detectors.get_detector_corners(data=image_identity["DATA"][SCIEXT_i-1],
                                                             wcs=image_identity["ASTROPYWCS"][SCIEXT_i-1])
        detector_square_list.append(np.concatenate([detector_corners["corners_world"], detector_corners["corners_world"]]))
    # Once you are done checking if the stars are inside each detector,
    # find out which stars are outside ALL detectors.
    gaia_query["is_inside_FPA"] = gaia_query[names_of_bool_columns_if_star_is_inside].any(axis=1)

    #### TODO: INDEPENDIZE THIS INTO STRAYCOR.PLOTS ##########
    ############# IF VERBOSE, MAKE AN INFIELD - OUTFIELD PLOT ###################

    if verbose:

        plt.figure(figsize=(1.618*10,10))
        plot_size = rs.plots.plot_stars_around(catalog=gaia_query, max_plot_size=100, min_plot_size=5, alpha=0.2)
        plot_radec_limits = rs.gaia.find_ra_dec_constraints(ra=RA_PNT, dec=DEC_PNT, radius=2*radius)
        for detector_square in detector_square_list:
            plt.plot(detector_square[:,0], detector_square[:,1], alpha=0.5, color="red")
        plt.xlim((plot_radec_limits["ra_max"], plot_radec_limits["ra_min"]))
        plt.ylim((plot_radec_limits["dec_min"], plot_radec_limits["dec_max"]))
        plt.show()

    ########################################################

    straylevel_list = []
    main_offender_list = []

    ## Prepare the coordinates of the stars that do not fall inside the Focal Plane Array ##
    ## This step is common for all SCI extensions #
    ra_stars_outside      = np.array(gaia_query["ra"][~gaia_query["is_inside_FPA"]])
    dec_stars_outside     = np.array(gaia_query["dec"][~gaia_query["is_inside_FPA"]])
    source_id_outside     = np.array(gaia_query["source_id"][~gaia_query["is_inside_FPA"]])
    cat_id_outside        = np.array(gaia_query["cat_id"][~gaia_query["is_inside_FPA"]])
    synthetic_mag_outside = np.array(gaia_query["mag_lambda"][~gaia_query["is_inside_FPA"]])
    stars_world_location = coordinates.SkyCoord(ra_stars_outside, dec_stars_outside, frame='icrs', unit="deg")

    irradiance_stars = const.c*(image_identity["FILTER_IDENTITY"]["filter_lambda_max"]-image_identity["FILTER_IDENTITY"]["filter_lambda_min"])/(image_identity["FILTER_IDENTITY"]["filter_lambda_ref"]**2)*((10**(-0.4*(synthetic_mag_outside+56.1)))*u.W/u.meter**2/u.Hz)

    ########################################
    # Here we estimate the stray-light
    ########################################
    for SCIEXT_i, name_column_is_star_inside_this_detector in tqdm(zip(SCIEXTS, names_of_bool_columns_if_star_is_inside)):


        if verbose: print("> Estimating stray-light in detector positions")

        if image_identity["TELESCOP"] == "Roman" or image_identity["TELESCOP"] == "ROMAN" or image_identity["TELESCOP"]=="RST" or image_identity["TELESCOP"]=="NGRST":
            #print("DEMO WARNING: The straylight estimator must be a method of rs.telescopes.CLASS")
            if verbose:
                print("rs.roman.roman_estimate_straylight_SCA")
                print("RA DEC irradiance_stars")
                print(ra_stars_outside)
                print(dec_stars_outside)
                print(irradiance_stars)



            straylevel_image_db = rs.roman.roman_estimate_straylight_SCA(data=image_identity["DATA"][SCIEXT_i-1],
                                                                        wcs=image_identity["ASTROPYWCS"][SCIEXT_i-1],
                                                                        SCA=image_identity["SCA"][SCIEXT_i-1],
                                                                        filter_identity=image_identity["FILTER_IDENTITY"],
                                                                        ra_stars=ra_stars_outside,
                                                                        dec_stars=dec_stars_outside,
                                                                        cat_id = cat_id_outside,
                                                                        source_id=source_id_outside,
                                                                        irradiance_stars=irradiance_stars,
                                                                        ra_point=RA_PNT, dec_point=DEC_PNT,
                                                                        pa_point=PA_PNT, verbose=verbose)
            straylevel_image_i = straylevel_image_db["straylight_SCA"]
            main_offender_image_i = straylevel_image_db["main_offender_SCA"]

            #     return({"straylight_SCA": straylight_SCA, "main_offender_SCA": main_offender_SCA})


        else:
            if verbose: print("estimate_straylight_in_detector_locations")
            straylevel_image_i = estimate_straylight_in_detector_locations(input_name=input_name,
                                                                    ext=SCIEXT_i,
                                                                    ra=RA_PNT,
                                                                    dec=DEC_PNT,
                                                                    stars_world_location=stars_world_location,
                                                                    synthetic_mag=synthetic_mag_outside,
                                                                    image_identity=image_identity,
                                                                    step=step,
                                                                    grid_method=grid_method,
                                                                    ndi_mode=ndi_mode,
                                                                    bool_image_mode=bool_image_mode,
                                                                    verbose=verbose)

        straylevel_list.append(straylevel_image_i)
        main_offender_list.append(main_offender_image_i)

    ########################################
    # Save the results to a fits file.
    ########################################

    # Stray-light
    input_fits = fits.open(input_name)
    data_output = []
    header_output = []
    for SCIEXT_i, straylevel_image_i in tqdm(zip(SCIEXTS, straylevel_list)):
        data_output.append(straylevel_image_i)
        header_output.append(image_identity["ASTROPYWCS"][SCIEXT_i-1].to_header())

    rs.utils.save_fits(array=data_output, name=output_name, header=header_output,
                       extname=None, overwrite=True, output_verify='silentfix')

    # Main-offender
    input_fits = fits.open(input_name)
    data_output = []
    header_output = []
    for SCIEXT_i, main_offender_i in tqdm(zip(SCIEXTS, main_offender_list)):
        data_output.append(main_offender_i)
        header_output.append(image_identity["ASTROPYWCS"][SCIEXT_i-1].to_header())

    main_offender_output_name = output_name.replace(".fits", "_main_off.fits")
    rs.utils.save_fits(array=data_output, name=main_offender_output_name, header=header_output,
                       extname=None, overwrite=True, output_verify='silentfix')



    print("Output saved in: " + output_name)

    return({"image_identity":image_identity,
            "straylevel_list": straylevel_list,
            "star_catalog": gaia_query,
            "output_name": output_name,
            "main_offender_output": main_offender_output_name,
            "detector_square": detector_square_list})

########################################
########################################
########################################

def estimate_straylight_in_detector_locations(input_name, ext, ra, dec,
                                              stars_world_location, synthetic_mag,
                                              image_identity, step, grid_method,
                                              ndi_mode, bool_image_mode, verbose):

    if bool_image_mode:
        # If the input is an image - Then calculate the straylight flux per position in the detector grid.
        # Else, only do it for the central location.
        # Make a grid of points in the detector.
        if verbose: print("> Generating a WCS grid in the detector...")
        detector_grid = rs.detectors.make_detector_grid(input_name=input_name, ext=ext, step=step, mode=grid_method)
        if verbose: print("> Detector WCS grid done.")

        in_detector_world_locations = detector_grid["grid_world"]
        in_detector_xy_locations = detector_grid["grid_xy"]

    else:
        # If the mode is RADEC - Then calculate the straylight flux in a single location.
        in_detector_world_locations = np.array([[ra], [dec]]) # Central world coordinates.
        in_detector_xy_locations = np.array([[0,0]]) # Arbitrary center of no-detector


    ######################
    ## Here we calculate the straylight for each pixel.
    ######################
    sum_of_straylevel_list_on_locations = []
    for i_reference_location in tqdm(range(len(in_detector_world_locations[0])), disable=not verbose, position=0, leave=True):
        #if verbose: print("> Organizing coordinates...")
        ra_pixel = in_detector_world_locations[0][i_reference_location]
        dec_pixel = in_detector_world_locations[1][i_reference_location]
        #if verbose: print("> Done")

        #if verbose: print("> Setting up coordinates of the pixels in the sky")
        pixel_world_location = coordinates.SkyCoord(ra_pixel, dec_pixel, frame='icrs', unit="deg")
        #if verbose: print("> Done")


        #if verbose: print("> Finding the distance of the stars to the pixel...")
        distance_from_detector_center_to_star = pixel_world_location.separation(stars_world_location) #catalog_close_stars["gaia_query_out_detector"]["dist"]
        #if verbose: print("> Done.")

        #if verbose: print("Demo warning: Position angle of detector set to 0 for NDI estimation")
        nstars = len(synthetic_mag)
        position_angle_detector_star = np.zeros(nstars)
        theta = distance_from_detector_center_to_star
        phi = position_angle_detector_star

        #if verbose: print("> Estimating the straylight...")

        # Find closest stars with Gaia.
        filter_name = image_identity["FILTER"]
        instrument = image_identity["INSTRUME"]
        telescope = image_identity["TELESCOP"]
        lambda_ref = image_identity["FILTER_IDENTITY"]["filter_lambda_ref"]
        exptime = image_identity["EXPTIME"]
        MJD = image_identity["EXPSTART"]

        straylevel = rs.ndi.straylight_flux(mag=synthetic_mag,
                                         theta=theta,
                                         phi=phi,
                                         filter_name=filter_name,
                                         instrument=instrument,
                                         telescope=telescope,
                                         exptime=exptime,
                                         mu_mode=False,
                                         ndi_mode=ndi_mode)
        #if verbose: print("> Done.")

        #if verbose: print("> Appending to list")
        sum_of_straylevel_list_on_locations.append(bn.nansum(straylevel))
        #if verbose: print("> Done.")
        # TODO: Organize the return in a coherent way.
        # return: A dict with
        # central_straylight = Straylight level at the center of the detector if provided. If a ra,dec is the input, then at this location.
        # straylight_image = Estimation of the straylight gradient across the FOV, if provided.
        # star_catalog =

    if bool_image_mode:
        # Reconstruct the straylight map
        x_pixel_list = detector_grid["grid_xy"][:,0]
        y_pixel_list = detector_grid["grid_xy"][:,1]
        n_pixels = len(x_pixel_list)
        xsize = int(np.max(x_pixel_list))
        ysize = int(np.max(y_pixel_list))
        x_lin = np.linspace(0,ysize,ysize)
        y_lin = np.linspace(0,xsize,xsize)
        xv, yv = np.meshgrid(y_lin, x_lin, indexing='ij')
        points = [(x_pixel_list[i], y_pixel_list[i]) for i in range(n_pixels)]
        straylight_image_interp =  interpolate.griddata(points, sum_of_straylevel_list_on_locations, (xv, yv), method="cubic").T
        straylight_output_name = input_name.replace(".fits", "_" + str(ext) + "_ofs.fits")
        input_fits = fits.open(input_name)
        rs.utils.save_fits(array=straylight_image_interp, name=straylight_output_name, header=input_fits[ext].header, extname="OFS")

    return(straylight_image_interp)
