# Alejandro S. Borlaff. NASA Ames Research Center. a.s.borlaff@nasa.gov / asborlaff@gmail.com
# June 27, 2023.
#
# ROSALIA main module
# This module will hold all the general programs to be interacting with the user
#
# Version log:
# v.1.0 - 27 June 2023. First loading of programs inherited from former straycor.
# v.2.0 - 29 Feb 2024. First working version in a public release. The main_offender program is now fully working and tested with Roman Dummy images. The main program to run is rosalia_stray, which calls main_offender internally.
##########################################################

# General modules
import os
import glob
import numpy as np
import pandas as pd
import rosalia as rs
import bottleneck as bn
import astropy.units as u

from tqdm import tqdm
from scipy import interpolate
import matplotlib.pyplot as plt
from astropy import coordinates
from astroquery.mast import Observations
from astropy.io import fits
import astropy.wcs as wcs
from astropy.time import Time
from astropy import constants as const
from astropy.coordinates import ICRS, Angle, SkyCoord

###########################

def rosalia_stray(ra, dec, PA, date, bandpass, exptime, prefix="", input_fits=None, radius=1,
                  g_mag_max=15, sun_block=False, verbose=False, catalog=None, 
                  figsize=(10,7), mu_vmin=None, mu_vmax=None):
    """
    Identify and estimate straylight from stars outside the field of view.

    Estimates the straylight from stars outside the field of view for Roman Space
    Telescope Wide Field Instrument observations.

    Args:
        ra (float): 
            Right ascension of the pointing, in degrees.

        dec (float): 
            Declination of the pointing, in degrees.

        PA (float): 
            Position angle of the observation, in degrees.
        
        date (astropy.time.Time): 
            Date of the observation in YYYY-MM-DDTHH:MM:SS format.
        
        bandpass (str): 
            Bandpass of the observation in Roman WFI filter names (e.g., F062, F087, F106, F129, F158, F184, F213).
        
        exptime (float): 
            Exposure time of the observation, in seconds.
        
        radius (float, optional): 
            Radius around the pointing to search for stars, in degrees. Default is 1.
        
        g_mag_max (float, optional): 
            Maximum g magnitude of the stars to consider in the straylight estimation. Default is 15.
        
        sun_block (bool, optional): 
            If True, the Sun will be removed from the star catalog. Default is False.
        
        verbose (bool, optional): 
            If True, print more information about progress. Default is False.
        
        catalog (pandas.DataFrame, optional): 
            User-provided catalog of stars. Must contain columns: "ra", "dec", "source_id", "cat_id", "mag_lambda".
            If provided, skips querying Gaia/2MASS/WISE catalogs. Default is None.

    Returns:
        pandas.DataFrame: Estimated straylight from each star outside the field
            of view with metadata. Columns: "source_id", "cat_id", "ra", "dec",
            "mag_lambda", "straylevel", "main_offender", "mosaic_name",
            "mosaic_name_scaled".

    History:
        v1 - 29 Feb 2024. First working version in a public release.

    Author:
        Alejandro S. Borlaff (NASA Ames Research Center, a.s.borlaff@nasa.gov)
    """
    from tqdm import tqdm
    import logging
    logger = logging.getLogger()
    logger.setLevel(logging.CRITICAL)

    if catalog is not None:
        if len(catalog)==1:
            # If the catalog has only one star, add a mag = inf at the same location to avoid problems with the plotting code. 
            dummy_catalog = catalog.copy()
            dummy_catalog["mag_lambda"]=np.inf
            dummy_catalog["cat_id"]=999
            dummy_catalog["source_id"]="Dummy star"
            catalog = pd.concat([catalog, dummy_catalog])


    if input_fits is None:
        # Make the Roman Dummy image
        roman_dummy_name = os.getcwd() + "/" + prefix + "WFI_" + bandpass +\
                                     "_RA_" + '{:07.3f}'.format(ra) +\
                                     "_DEC_" + '{:07.3f}'.format(dec) +\
                                     "_MJD_" + '{:07.5f}'.format(date.mjd) +\
                                     "_PA_" + '{:06.2f}'.format(PA) + ".fits"
        output_name = roman_dummy_name.replace(".fits","_stray.fits")
        central_coords = SkyCoord(ra, dec, frame="icrs", unit="deg")

        input_name = rs.roman.create_roman_dummy(point=central_coords, date=date,
                                               band=bandpass, PA=PA, exptime=exptime,
                                               output=roman_dummy_name)
        print(input_name)
    else: 
        input_name = input_fits
        output_name = input_name.replace(".fits","_stray.fits")

    # Get the image identity
    exposure_dict = rs.utils.exposure_inspector(input_name=input_name, verbose=verbose, lite=True)


    ra = exposure_dict["RA_TARG"]
    dec = exposure_dict["DEC_TARG"]
    PA = exposure_dict["PA"]
    date = exposure_dict["EXPSTART"]
    bandpass = exposure_dict["FILTER"]
    exptime = exposure_dict["EXPTIME"]
 
    if input_name is None:
        input_name = rs.roman.create_roman_dummy(point=exposure_dict["point"],
                                                 date=date,
                                                 band=bandpass,
                                                 PA=PA,
                                                 exptime=exptime,
                                                 output=output_name.replace(".fits", "_dummy.fits"))

    # If the input name is an ASDF, transform it to a compiled FITS.
    if isinstance(input_name, (list, np.ndarray, pd.Series)):
        if input_name[0].split(".")[-1] == "asdf":
            fits_input_name = input_name[0].replace(".asdf",".fits")
            if verbose:
                print(input_name)
            rs.utils.convert_ASDF_to_FITS(asdf_list = input_name, output=fits_input_name)

    else:
        fits_input_name = input_name

    main_offender_db = main_offender(input_name=fits_input_name,
                                     radius=radius,
                                     g_mag_max=g_mag_max,
                                     verbose=verbose,
                                     input_catalog=catalog,
                                     match_inside_stars=False,
                                     sun_block = sun_block,
                                     output_name = output_name)


    if verbose > 0: print("Generating ROSALIA summary report...")

    drz_name, scaled_drz_name = rs.utils.run_swarp(pattern=main_offender_db["output_name"], 
                                                   outname=main_offender_db["output_name"].replace(".fits","_drz.fits"), scale=0.1)
    
    #mainoff_name, scaled_mainoff_name = rs.utils.run_swarp(pattern=main_offender_db["main_offender_output"],
    #                                                       outname=main_offender_db["main_offender_output"].replace(".fits","_drz.fits"), scale=0.1)

    data, header = rs.utils.reproject_roman_wfi_fits(input_name=scaled_drz_name, 
                                                     input_ext=1,
                                                     reference_name=main_offender_db["main_offender_output"], 
                                                     reference_ext=rs.telescopes.Roman.WFI_SCAs)
    
    scaled_main_off = main_offender_db["main_offender_output"].replace(".fits", "_scaled.fits")

    catalog_name = main_offender_db["star_catalog"]

    rs.utils.save_fits(array=data, name=scaled_main_off, header=header, overwrite=True)

    # Writing necessary keywords in the output mosaics. 
    keywords = ["RA_TARG", "DEC_TARG", "EXPSTART", "EXPTIME", "FILTER", "WAVEREF", "WAVEMIN", "WAVEMAX", "TELESCOP", "INSTRUME", "DETECTOR"]
    key_values = rs.utils.get_keys_from_header([fits_input_name], keywords, ext=0)
    rs.utils.write_parameters_list([drz_name], keywords, key_values, ext=0)
    rs.utils.write_parameters_list([scaled_drz_name], keywords, key_values, ext=0)
    rs.utils.write_parameters_list([scaled_drz_name], ["PIXSCALE"], [[1]], ext=0)
    rs.utils.write_parameters_list([scaled_drz_name], ["REBINNED"], [[10]], ext=0)

    rs.utils.write_parameters_list([main_offender_db["main_offender_output"]], keywords, key_values, ext=0)
    rs.utils.write_parameters_list([scaled_main_off], keywords, key_values, ext=0)
    rs.utils.write_parameters_list([scaled_main_off], ["PIXSCALE"], [[1]], ext=0)
    rs.utils.write_parameters_list([scaled_main_off], ["REBINNED"], [[10]], ext=0)

    # Make the plots. 
    if verbose > 0: print("Plot: Stray-light surface brightness magnitude.")
    fe2mu_png = rs.plots.make_stray_plot(input_name=scaled_drz_name,
                                        ext=1, mode="fe2mu",
               color_label = "Surface brightness (mag arcsec$^{-2}$)", figsize=figsize,mu_vmin=mu_vmin, mu_vmax=mu_vmax)

    if verbose > 0: print("Plot: Stray-light surface brightness flux.")
    fe_png = rs.plots.make_stray_plot(input_name=scaled_drz_name, ext=1, mode="fe", 
                color_label = "Flux (e/s/px)",cmap="RdYlBu_r", figsize=figsize)

    if verbose > 0: print("Plot: Main offender map.")
    main_offender_png = rs.plots.make_stray_plot(input_name=scaled_main_off, ext=0, 
                                                 mode="main_offender", 
                catalog_name=catalog_name, 
                color_label = "Main offending source (ID)", figsize=figsize)

    
    if verbose > 0: print("Plot: Source environment map.")
    try:
        stars_around_png = rs.plots.make_stars_around_plot(fits_input_name, catalog_name, output_name=None, figsize=figsize)
    except:
        print(rs.plots.style.YELLOW + "WARNING: Could not make stars around plot. Check if the catalog has the right columns." + rs.plots.style.RESET)
        stars_around_png = ""


    if verbose > 0: print("Plot: Location of Main offenders on NDI map.")
    try:
        ndi_mainoff_name = rs.plots.plot_ndi_main_offenders(input_name=fits_input_name, scaled_main_off=scaled_main_off, 
                                                        catalog_name=catalog_name, ndi_level=2, figsize=figsize)
    except:
        print(rs.plots.style.YELLOW + "WARNING: Could not make stars Main offender plot. Check if the catalog has the right columns." + rs.plots.style.RESET)
        ndi_mainoff_name = ""

    rs.utils.execute_cmd("magick -adjoin " +\
                        fe2mu_png + " " +\
                        fe_png + " "+\
                        main_offender_png + " " +\
                        stars_around_png + " "+\
                        ndi_mainoff_name + " "+\
                        fits_input_name.replace(".fits", ".pdf"))


    return(main_offender_db)

###########################

def rosalia_zody(ra, dec, PA, date, bandpass, exptime, verbose=False, output_name=None, output_units="e/s"):

    from tqdm import tqdm
    import logging
    logger = logging.getLogger()
    logger.setLevel(logging.CRITICAL)

    # Make the Roman Dummy image
    roman_dummy_name = os.getcwd() + "/WFI_" + bandpass +\
                                     "_RA_" + '{:07.3f}'.format(ra) +\
                                     "_DEC_" + '{:07.3f}'.format(dec) +\
                                     "_MJD_" + '{:07.5f}'.format(date.mjd) +\
                                     "_PA_" + '{:06.2f}'.format(PA) + ".fits"

    central_coords = SkyCoord(ra, dec, frame="icrs", unit="deg")

    roman_dummy_name = rs.roman.create_roman_dummy(point=central_coords, date=date,
                                               band=bandpass, PA=PA, exptime=exptime,
                                               output=roman_dummy_name)
    print(roman_dummy_name)

    # Get the image identity
    exposure_identity = rs.utils.exposure_inspector(input_name=roman_dummy_name, verbose=verbose, lite=True)



    zodiacal_background_list = []
    zodiacal_background_unit_list = []
    if output_name is None:
        output_name = roman_dummy_name.replace(".fits", "_zody.fits")

    for i in tqdm(range(len(exposure_identity["SCIEXTS"]))):
        zodiacal_background = rs.sky.get_zodiacal_background(input_name=roman_dummy_name,
                                                      ext=exposure_identity["SCIEXTS"][i],
                                                      wavelength=exposure_identity["FILTER"],
                                                      telescope=exposure_identity["TELESCOP"],
                                                      instrument=exposure_identity["INSTRUME"],
                                                      detector=exposure_identity["DETECTOR"],
                                                      expstart=exposure_identity["EXPSTART"],
                                                      step=1000, zody_mode="zodipy",
                                                      nbins_wavelength=20, obslocin=3,
                                                      grid_method="random", output_units=output_units,
                                                      verbose=False)
        # print(zodiacal_background)
        zodiacal_background_list.append(zodiacal_background.value)
        zodiacal_background_unit_list.append(zodiacal_background.unit.to_string())

    ########################################
    # Save the results to a fits file.
    ########################################
    data_output = []
    header_output = []
    for i, SCIEXT_i, zodiacal_background_i in tqdm(zip(range(len(zodiacal_background_list)), exposure_identity["SCIEXTS"], zodiacal_background_list)):
        data_output.append(zodiacal_background_i)
        temp_header = exposure_identity["ASTROPYWCS"][i].to_header()
        temp_header["UNITS"] = zodiacal_background_unit_list[i]
        header_output.append(temp_header)

    rs.utils.save_fits(array=data_output, name=output_name, header=header_output,
                       extname=None, overwrite=True, output_verify='silentfix')

    # Make the scaled model and the summary plot.
    drz_name, scaled_drz_name = rs.utils.run_swarp(pattern=output_name, 
                                                   outname=output_name.replace(".fits","_drz.fits"), scale=0.1)
    
    ext = 1
    rs.plots.make_stray_plot(input_name=scaled_drz_name, ext=ext, mode="fe2mu", catalog_name=None, 
                    vmin=None, vmax=None, 
                    color_label = 'Surface brightness (mag arcsec$^{-2}$)',
                    cmap="RdYlBu", output_name=None, figsize=(10,7), mu_vmin=None, mu_vmax=None)

    print("Output saved in: " + output_name)

    return({"image_identity":exposure_identity,
            "zodi_list": zodiacal_background_i,
            "output_name": output_name})


###########################


def rosalia_psf(ra, dec, PA, g_mag_max, date, bandpass, exptime, input_catalog=None, verbose=False):
    #######################################
    # rosalia_psf: Alejandro S. Borlaff. NASA/Ames STA. a.s.borlaff@nasa.gov
    # -------------------------------
    # The objective of this program is to make a model of the stars inside a Roman WFI image
    # --------------------------------
    # History:
    # v1 - 22 January 2026. First working version.
    #
    #######################################

    '''
    rosalia_psf: Alejandro S. Borlaff. NASA Ames Research Center.
    Model the stars inside a Roman WFI image. This is useful for estimating the straylight from stars 
    inside the field of view, and for subtracting the stars from the image.

    Args:
        ra (float): 
            Right ascension of the pointing, in degrees.
        
        dec (float): 
            Declination of the pointing, in degrees.
        
        PA (float): 
            Position angle of the observation, in degrees.
        
        g_mag_max (float):
            Maximum g magnitude of the stars to consider in the model.
        
        date (astropy.time.Time): 
            Date of the observation in YYYY-MM-DDTHH:MM:SS format.
        
        bandpass (str): 
            Bandpass of the observation in Roman WFI filter names (e.g., F062, F087, F106, F129, F158, F184, F213).
        
        exptime (float): 
            Exposure time of the observation, in seconds.

        input_catalog (pandas.DataFrame, optional): 
            User-provided catalog of stars. Must contain columns: "ra", "dec", "source_id", "cat_id", "mag_lambda".
        
        verbose (bool, optional): 
            If True, print more information about progress. Default is False.

    '''
    from tqdm import tqdm
    import logging
    logger = logging.getLogger()
    logger.setLevel(logging.CRITICAL)

    # Make the Roman Dummy image
    roman_dummy_name = os.getcwd() + "/WFI_" + bandpass +\
                                     "_RA_" + '{:07.3f}'.format(ra) +\
                                     "_DEC_" + '{:07.3f}'.format(dec) +\
                                     "_MJD_" + '{:07.5f}'.format(date.mjd) +\
                                     "_PA_" + '{:06.2f}'.format(PA) + ".fits"

    central_coords = SkyCoord(ra, dec, frame="icrs", unit="deg")

    roman_dummy_name = rs.roman.create_roman_dummy(point=central_coords, date=date,
                                                   band=bandpass, PA=PA, exptime=exptime,
                                                   output=roman_dummy_name)
    print(roman_dummy_name)

    # Get the image identity
    image_identity = rs.utils.exposure_inspector(input_name=roman_dummy_name, verbose=verbose, lite=True)


    # Get the catalog of the stars around the FOV
    # hybrid_catalog = rs.psf.find_stars_inside_detector(input_name=roman_dummy_name, g_mag_max=g_mag_max, verbose=verbose)
    #hybrid_catalog = rs.psf.find_gaia_stars_around_image(lambda_ref=lambda_ref, input_name=input_name, ext=ext,
                                                          # ra=ra, dec=dec, MJD=MJD, radius=radius, g_mag_max=g_mag_max, verbose=verbose)
    #hybrid_catalog = gaia_query_dict["gaia_query"]

    source_catalog_filename = roman_dummy_name.replace(".fits", "_source_catalog.csv")
    if input_catalog is not None:
        hybrid_catalog = input_catalog
    
    else:
        hybrid_catalog = rs.psf.get_hybrid_catalog(ra=ra,
                                               dec=dec,
                                               radius=1,
                                               lambda_ref=image_identity["FILTER_IDENTITY"]["filter_lambda_ref"],
                                               MJD=date.mjd,
                                               observer=image_identity["TELESCOP"],
                                               g_mag_max=g_mag_max, verbose=verbose, 
                                               query_filename=source_catalog_filename)

    # def rosalia_stray(input_name, output_name="rosalia_stray_output.fits", radius=1, g_mag_max=15, sun_block=False, verbose=False, catalog=None):

    # input_name = "/Users/aborlaff/NASA/ROSALIA/notebooks/DEVEL/default_roman_dummy.fits"
    # verbose=1
    # g_mag_max=15
    #


    # Generate the star stamps (PSFs)
    print("TO DO: Make stamps with a more reasonable size. Dim stars can have smaller PSFs.")
    print("To do this, make a profile of the Roman / PSF, and find out when would it be essentially 0.")
    star_stamps = rs.psf.generate_star_stamps(hybrid_catalog=hybrid_catalog, image_identity=image_identity)


    # Now combine all the stamps in the mosaiced frame and blot back to the single SCAs.
    # This is more efficient than reprojecting each star into all SCAs.
    # Flattening the list of lists.
    star_stamps_flat = []
    for i in range(len(star_stamps)):
        star_stamps_flat = star_stamps_flat + star_stamps[i]

    # Making the combined frame.
    os.system("swarp -dd > swarp.conf")
    swarp_cmd_str = ""
    for star_stamp in star_stamps_flat:
        swarp_cmd_str = swarp_cmd_str + '"' + star_stamp +'" '

    if verbose > 1: print("Combining star stamps into WCS frame...")
    cmd = "swarp -c swarp.conf -SUBTRACT_BACK N -BLANK_BADPIXELS Y -COMBINE_TYPE SUM -VERBOSE_TYPE QUIET " + swarp_cmd_str
    if verbose > 2: print(cmd)
    rs.utils.execute_cmd(cmd)  # Run swarp on all the SCAs
    star_swarp_name = roman_dummy_name.replace(".fits", "_stars_drz.fits")
    rs.utils.execute_cmd("mv coadd.fits " + star_swarp_name) # Make a compressed version, for easiest visualization.

    # Now blot back to the dummy SCA per SCA frame
    from reproject import reproject_interp
    # Let's make a dummy copy to reproject the stars into
    roman_dummy = fits.open(roman_dummy_name)
    star_model = fits.open(star_swarp_name)

    if verbose > 1: print("Storing stars in each SCA")
    for SCIEXT_i in tqdm(image_identity["SCIEXTS"]):
        # Open the star fits
        star_reprojected, footprint = reproject_interp(star_model[0], roman_dummy[SCIEXT_i].header, parallel=True)
        star_reprojected[np.isnan(star_reprojected)] = 0
        roman_dummy[SCIEXT_i].data = roman_dummy[SCIEXT_i].data + star_reprojected

    roman_dummy.verify("silentfix")
    star_output_name = roman_dummy_name.replace(".fits", "_stars.fits")
    roman_dummy.writeto(star_output_name, overwrite=True)
    if verbose > 1: print("In-field stray-light model completed: " + star_output_name)
    return(star_output_name)


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
        if verbose > 0: print(rs.plots.styles.YELLOW + "WARNING: For testing purposes the PSF is fixed" + rs.plots.styles.RESET)
        psf_name = os.environ["ROSALIACACHE"] + "/CORE/PSF_ARCHIVE/f814w00_arith.fits"
        if verbose > 0: print("PSF: " + psf_name)

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
########################################

def main_offender(input_name=None, ext=None, ra=None, dec=None, phi=0,
                  filter_name=None, instrument=None, telescope=None, detector=None,
                  exptime=None, expstart=None, radius=1, g_mag_max=False,
                  step=200, grid_method="random", verbose=False, ndi_mode="legacy",
                  input_catalog=None, match_inside_stars=False, sun_block=False,
                  output_name="main_offender_default_output.fits"):
    #print("input_catalog")

    #print(input_catalog)
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
    # Resetting the Roman / WFI Loading bar
    rs.plots.ascii_progress_focal_plane.canvas = np.copy(rs.plots.ascii_progress_focal_plane.canvas_zero)

    if input_name is not None: # If the input is an image, use exposure inspector to get the arguments to the following programs.
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
        # Find the central coordinates of the image - If there are more than one science extensions,
        # then find the center of them all. <--- exposure_inspector does this for us
        RA_PNT      = image_identity["RA_TARG"]
        DEC_PNT     = image_identity["DEC_TARG"]
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


    #print("Demo warning: ADD FEATURE planets: https://github.com/skyfielders/python-skyfield/tree/master")

    #print("input_catalog")
    #print(input_catalog)

    if input_name is None:
        source_catalog_filename = str(RA_PNT) + "_" + str(DEC_PNT) + "_source_catalog.csv"
    else:
        source_catalog_filename = input_name.replace(".fits", "_source_catalog.csv")

    if input_catalog is None:
        #print("input catalog is none!")
        existing_catalog_name = output_name.replace(".fits", "_catalog.csv")
        if os.path.exists(existing_catalog_name):
            print("WARNING: Loading existing catalog! Remove " + existing_catalog_name + " if this is a mistake.")
            input_catalog = pd.read_csv(existing_catalog_name)

        else:
            if verbose: print("> Querying stars in the surroundings using ESA/Gaia Archive")
 
            if radius > 0.5:
                print("INFO: radius parameter (minimum distance to search for individual stars) is > 0.5 degrees.")
                print("Gaia/2MASS/WISE query database can take several minutes to process. Please be patient.")

            # Find the stars around the central coordinate of the scene.
            # Save the stellar catalog into a catalog object, and run main_offender as if input_catalog was set by the User.
            loader = rs.plots.Loader("Querying Gaia/2MASS/WISE/JPL Horizons databases. This might take a few minutes...",
                                     "All-sky source map constructed.", 0.05).start()

            hybrid_catalog = rs.psf.get_hybrid_catalog(ra=RA_PNT, dec=DEC_PNT,
                                                   radius=radius,
                                                   lambda_ref=lambda_ref,
                                                   MJD=MJD,
                                                   observer=telescope,
                                                   g_mag_max = g_mag_max,
                                                   verbose=verbose,
                                                   query_filename=source_catalog_filename)
            loader.stop()
    else:
        hybrid_catalog = input_catalog

    hybrid_catalog.to_csv(source_catalog_filename)


    # If sun_block is True, then remove the Sun from the catalog.
    if sun_block:
        hybrid_catalog = hybrid_catalog[~(hybrid_catalog["source_id"] == "Sun")]
        hybrid_catalog = hybrid_catalog[~(hybrid_catalog["source_id"] == "Earth")]

    # Find where each star lands (detector ID or outside FOV)
    names_of_bool_columns_if_star_is_inside = []
    detector_square_list = []

    """
    If the input file is a multi-extension fits, then exposure_inspector will scan for extensions with EXTNAME = SCI.
    The extension ID in the FITS file will be stored in SCIEXTS = image_identity["SCIEXTS"].

    In that case, hybrid_catalog, the catalog of stars, will have a set of N columns called in_SCI[i] (boolean), where the catalog
    stores if that particular star is inside each detector or not.

    """

    for SCIEXT_i in tqdm(SCIEXTS):
        if verbose: print("> Identifying which stars are inside the FOV and which are outside...")
        infield_stars = rs.psf.identify_stars_in_out_field(data_shape=image_identity["DATA"][SCIEXT_i-1].shape,
                                                           wcs=image_identity["ASTROPYWCS"][SCIEXT_i-1],
                                                           catalog=hybrid_catalog,
                                                           verbose=verbose)

        name_column_is_star_inside_this_detector = "in_SCI" + str(SCIEXT_i)
        names_of_bool_columns_if_star_is_inside.append(name_column_is_star_inside_this_detector)

        hybrid_catalog[name_column_is_star_inside_this_detector] = infield_stars["bool_isIn"]

        # If verbose, make a plot of the stars with the footprint.
        detector_corners = rs.detectors.get_detector_corners(wcs=image_identity["ASTROPYWCS"][SCIEXT_i-1])
        detector_square_list.append(np.concatenate([detector_corners["corners_world"], detector_corners["corners_world"]]))
    # Once you are done checking if the stars are inside each detector,
    # find out which stars are outside ALL detectors.
    hybrid_catalog["is_inside_FPA"] = hybrid_catalog[names_of_bool_columns_if_star_is_inside].any(axis=1)

    #### TODO: INDEPENDIZE THIS INTO STRAYCOR.PLOTS ##########
    ############# IF VERBOSE, MAKE AN INFIELD - OUTFIELD PLOT ###################

    if verbose:

        plt.figure(figsize=(1.618*10,10))
        plot_size = rs.plots.plot_stars_around(catalog=hybrid_catalog, max_plot_size=100, min_plot_size=5, alpha=0.2)
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
    ra_stars_outside      = np.array(hybrid_catalog["ra"][~hybrid_catalog["is_inside_FPA"]])
    dec_stars_outside     = np.array(hybrid_catalog["dec"][~hybrid_catalog["is_inside_FPA"]])
    source_id_outside     = np.array(hybrid_catalog["source_id"][~hybrid_catalog["is_inside_FPA"]])
    cat_id_outside        = np.array(hybrid_catalog["cat_id"][~hybrid_catalog["is_inside_FPA"]])
    synthetic_mag_outside = np.array(hybrid_catalog["mag_lambda"][~hybrid_catalog["is_inside_FPA"]])
    stars_world_location = coordinates.SkyCoord(ra_stars_outside, dec_stars_outside, frame='icrs', unit="deg")

    irradiance_stars = const.c*(image_identity["FILTER_IDENTITY"]["filter_lambda_max"]-image_identity["FILTER_IDENTITY"]["filter_lambda_min"])/(image_identity["FILTER_IDENTITY"]["filter_lambda_ref"]**2)*((10**(-0.4*(synthetic_mag_outside+56.1)))*u.W/u.meter**2/u.Hz)

    ########################################
    # Here we estimate the stray-light
    ########################################
    # Reset the Roman / WFI loading bar:

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


            #
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
    rs.utils.save_fits(array=data_output, 
                       name=main_offender_output_name, 
                       header=header_output,
                       extname=None, 
                       overwrite=True, 
                       output_verify='silentfix')


    # Let's do one more step to include the needed metadata from the dummy file. 
    stray_image = fits.open(output_name)
    main_offender_image = fits.open(main_offender_output_name)
    exposure_identity_keys_to_copy = ["TELESCOP", "INSTRUME", "DETECTOR", "FILTER", "RA_TARG", "DEC_TARG", 
                                      "RA_PNT", "DEC_PNT", "X_PNT", "Y_PNT", "PA",
                                      "EXPTIME", "EXPSTART", "EXPSTART_ISOT"]
    for key in exposure_identity_keys_to_copy:
        stray_image[0].header[key] = image_identity[key]
        main_offender_image[0].header[key] = image_identity[key]
    
    # Additional keys: 
    stray_image[0].header["WAVEREF"] = image_identity["FILTER_IDENTITY"]["filter_lambda_ref"].to("nm").value
    main_offender_image[0].header["WAVEREF"] = image_identity["FILTER_IDENTITY"]["filter_lambda_ref"].to("nm").value
    stray_image[0].header["WAVEMIN"] = image_identity["FILTER_IDENTITY"]["filter_lambda_min"].to("Angstrom").value
    main_offender_image[0].header["WAVEMIN"] = image_identity["FILTER_IDENTITY"]["filter_lambda_min"].to("nm").value
    stray_image[0].header["WAVEMAX"] = image_identity["FILTER_IDENTITY"]["filter_lambda_max"].to("nm").value
    main_offender_image[0].header["WAVEMAX"] = image_identity["FILTER_IDENTITY"]["filter_lambda_max"].to("nm").value
    
    # Verify, save and close
    stray_image.verify("silentfix")
    main_offender_image.verify("silentfix")
    stray_image.writeto(output_name, overwrite=True)
    main_offender_image.writeto(main_offender_output_name, overwrite=True)
    

    print("Output saved in: " + output_name)

    return({"image_identity":image_identity,
            "straylevel_list": straylevel_list,
            "star_catalog": source_catalog_filename,
            "output_name": output_name,
            "main_offender_output": main_offender_output_name,
            "detector_square": detector_square_list})

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
########################################
