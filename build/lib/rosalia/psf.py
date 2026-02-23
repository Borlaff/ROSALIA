# Alejandro S. Borlaff. NASA Ames Research Center. a.s.borlaff@nasa.gov / asborlaff@gmail.com
# January 20, 2023.
#
# STRAYCOR/PSF module
# This module will hold all the programs related to the modelling and removal
# of the PSF.
#
# Version log:
# v.1.0 - 20 Enero 2023. First loading of programs inherited from former monolithic straycor.py
#
##########################################################

############################
import os
import sys
import pandas as pd
import numpy as np
import bottleneck as bn
from tqdm import tqdm
from astropy.io import fits
import astropy.wcs as astropy_wcs
import matplotlib.pyplot as plt
import astropy.units as u
from astroquery.mast import Observations
from astropy.coordinates import match_coordinates_sky, ICRS, Angle, SkyCoord
import rosalia as rs

# Suppress warnings. Comment this out if you wish to see the warning messages
import warnings
warnings.filterwarnings('ignore')

# Custom modules
#sys.path.append("/Users/aborlaff/NASA/STRAYCOR/")


# STRAYCOR modules


psf_archive = os.path.dirname(rs.utils.__file__) + "/../PSF_ARCHIVE/"
plt.style.use(os.path.dirname(rs.utils.__file__) + "/CORE/presi_style.mplstyle")


def scale_and_subtract_stars(input_name, ext, lambda_ref, psf_name, g_mag_max=False, MJD=None, clean=False, verbose=False):
    ## TODO: Compute the scale factor for the object at (x,y)=(53,69) for
    ## the PSF (psf.fits). Compute it in the ring 20-30 pixels.
    # $ astscript-psf-scale-factor image.fits --mode=img \
    #      --center=53,69 --normradii=20,30 --psf=psf.fits

    ## Iterate over a catalog with RA,Dec positions of stars that are in
    ## the input image to compute their scale factors.
    #$ asttable catalog.fits | while read -r ra dec mag; do \
    #    astscript-psf-scale-factor image.fits \
    #        --mode=wcs \
    #        --psf=psf.fits \
    #        --center=$ra,$dec --quiet \
    # --normradii=20,30 > scale-"$ra"-"$dec".txt; done

    # List of input images: Checking that it is a list.
    if isinstance(input_name, (list, pd.core.series.Series, np.ndarray)):
        input_list = input_name
        if verbose:
            print("Executing for images:")
            print(input_list)

    else:
        if verbose:
            print("Executing for single image " + input_name)
        input_list = [input_name]

    # Now we check if the extension is just one, or a list.
    if isinstance(ext, (list, pd.core.series.Series, np.ndarray)):
        ext_list = ext
    else:
        ext_list = [ext]

    ###########################
    # Fundamental estimations #
    ###########################
    data, wcs = rs.utils.get_data_and_wcs(input_name=input_name, ext=ext)

    star_search_radius = rs.utils.find_max_angular_size_of_image(data, wcs)

    # Now we execute the scale_and_subtract_stars_single_image program
    # File by file, and extension by extension.
    ###########################
    # Check that PSF is in the right format
    ##################

    rs.utils.check_number_of_extensions(psf_name, expected_next=2, verbose=verbose)

    subtracted_results = []
    output_name_list = []
    for input_name_i in tqdm(input_list):

        # We open a copy of the input_name_i that we are modifying.
        input_fits_i = fits.open(input_name_i)
        output_name_i =  input_name_i.replace(".fits", "_nostars.fits")
        for ext_i in ext_list:

            matched_stars_in_detector = find_stars_inside_detector(input_name=input_name_i,
                                                                   ext=ext_i,
                                                                   lambda_ref=lambda_ref,
                                                                   radius=star_search_radius,
                                                                   g_mag_max=g_mag_max,
                                                                   MJD=MJD, clean=False,
                                                                   verbose=verbose)
            if verbose: print(matched_stars_in_detector)
            # A flat sky background correction is required to fit and subtract the PSF from stars
            # Otherwise it is impossible to scale correctly the PSF to the ring around the identified stars.
            # However, we do not want to do such correction with gradients just yet:
            # Ideally, we want to leave the gradient correction to the NDI, and then the residuals with Gnuastro.

            #rs.sky.correct_flat_sky(input_name=input_name_i, ext=ext_i, clean=clean)

            ra_list = matched_stars_in_detector["matched_catalog_image"]["RA"]
            dec_list = matched_stars_in_detector["matched_catalog_image"]["DEC"]

            ds9_region_filename = input_name_i.replace(".fits", "_ext" + str(ext_i) + "_instars.reg")
            rs.utils.make_ds9_region(ra=ra_list, dec=dec_list,
                                     output=ds9_region_filename)

            # (input_name, ext, ra, dec, psf_name, clean=True):
            dictionary_with_results_from_star_subtraction = scale_and_subtract_stars_single_image(input_name = input_name_i, ext = ext_i,
                                                                               ra = ra_list, dec = dec_list,
                                                                               psf_name = psf_name, clean=False, verbose=verbose)
            subtracted_results.append(dictionary_with_results_from_star_subtraction)
            # Now we modify the extension with the residual image
            nostars_fits = fits.open(dictionary_with_results_from_star_subtraction["residuals"])
            input_fits_i[ext_i].data = nostars_fits[0].data
            nostars_fits.close()

        input_fits_i.verify("silentfix")
        input_fits_i.writeto(output_name_i, overwrite=True)
        input_fits_i.close()
        output_name_list.append(output_name_i)

        if clean:
            #for i in dictionary_with_results_from_star_subtraction["residuals"]:
            #    rs.utils.execute_cmd("rm " + i)

            rs.utils.execute_cmd("rm " + input_name_i.replace(".fits", "*detected.fits"))
            rs.utils.execute_cmd("rm " + input_name_i.replace(".fits", "*instars.reg"))
            rs.utils.execute_cmd("rm " + input_name_i.replace(".fits", "*_matched_stars_inside_detector.fits"))
            rs.utils.execute_cmd("rm " + input_name_i.replace(".fits", "*_gaia_stars_inside_detector.csv"))
            rs.utils.execute_cmd("rm " + input_name_i.replace(".fits", "*_gaia.dat"))
            rs.utils.execute_cmd("rm " + input_name_i.replace(".fits", "*_detected_segmented_cat.fits"))

    subtracted_results_dataframe = pd.DataFrame(subtracted_results)

    return(output_name_list)
############################


############################
def astscript_psf_scale_factor(input_name, ra, dec, psf_name, clean=False, verbose=False):
    # Python wrapper of astscript-psf-scale-factor
    #
    # Based on https://www.gnu.org/software/gnuastro/manual/html_node/Invoking-astscript_002dpsf_002dscale_002dfactor.html
    #
    # Shell based example
    # cmd = "astscript-psf-scale-factor mastDownload/HST/j9en0pnbq/j9en0pnbq_flc_ext1.fits --mode=wcs --center=329.03478737,1.4031257 --normradii=10,20 --psf=PSF/psf_hst_test.fits  | tail -1"
    # scaling_factor = os.popen(cmd).read()
    # print(scaling_factor)



    #print("WARNING:")
    if verbose:
        print("Estimating norm_radiis:")
        print("RA:" + str(ra))
        print("DEC:" + str(dec))

    max_size_of_a_star = 600 # pixels
    # Lets gather the max detection radiis for all the stars

    rmax = measure_maxradii(input_name=input_name, ra=ra, dec=dec, rmax=max_size_of_a_star, clean=clean, verbose=verbose)
    norm_radii_min = int(rmax/3.)
    norm_radii_max = int(2*rmax/3.)

    if verbose:
        print("norm_radii_min:" + str(norm_radii_min))
        print("norm_radii_max:" + str(norm_radii_max))


    # We construct the cmd to send to shell.
    cmd = "astscript-psf-scale-factor " + input_name +\
           " --mode=wcs --center=" + str(ra) + "," + str(dec) +\
           " --normradii=" + str(norm_radii_min) +"," + str(norm_radii_max) +\
           " --psf=" + psf_name + "  | tail -1"

    # We execute the command and get the output, but only the last line ( blahblah .... | tail -1)

    if verbose:
        print(cmd)

    scaling_factor = float(rs.utils.execute_cmd(cmd))
    if verbose: print("Scaling factor: " + str(scaling_factor))
    return(scaling_factor)
############################


def astscript_radial_profile(input_name, ra, dec, rmax, clean=False, verbose=False):
    # This is a wrapper for Gnuastro astcript_radial_profile
    output_name = input_name.replace(".fits", "_profile_ra" + str(ra) + "_dec" + str(dec) + ".fits")

    # We construct the cmd to send to shell.
    cmd = "astscript-radial-profile " + input_name +\
           " --mode=wcs --measure=median --center=" + str(ra) + "," + str(dec) +\
           " --rmax=" + str(rmax) +\
           " -o " + output_name

    # We execute the command and get the output, but only the last line ( blahblah .... | tail -1)
    if not clean:
        if verbose:
            print(cmd)

    rs.utils.execute_cmd(cmd)

    if not os.path.exists(output_name):
        if verbose:
            print("Caution: Profile not generated: " + input_name + " ra: " + str(ra) + " dec " + str(dec))
            print("Most likely, the target is not in the field of view.")
        return(False)

    if verbose:
        print("Profile for " + input_name + " ra: " + str(ra) + " dec " + str(dec))
        print("saved in " + output_name)



    return({"profile_name": output_name})

############################



def measure_maxradii(input_name, ra, dec, rmax, saturation_level=None, clean=False, verbose=False, clean_profile=False):
    star_profile_dict = astscript_radial_profile(input_name=input_name, ra=ra, dec=dec, rmax=rmax, verbose=verbose)

    # If the target is not in the FOV, then astscript_radial_profile will return None.
    # Propagate so it can be detected by other programs and the next line does not raise an error.

    if star_profile_dict is False:
        print("Caution: Star profile not found! Returning None")
        return(False)

    star_profile = fits.open(star_profile_dict["profile_name"])

    r = star_profile[1].data["RADIUS"]
    sky = star_profile[1].data["MEDIAN"][int(rmax/2):int(rmax-1)]
    sb = star_profile[1].data["MEDIAN"] - bn.nanmedian(sky) #+ 0.5# - np.nanmin(db[1].data["MEAN"])

    # We make a first approximation to the rlim
    noise_level = bn.nanstd(sky)

    if verbose:
        plt.axhline(noise_level, label="Noise level 1")

    rlim = star_profile[1].data["RADIUS"][np.where(((sb - 3*noise_level) < 0) & (r > 2))[0][0]] # We put the (r>2) to avoid the saturated core
    if verbose:
        print("Rlim 1: " + str(rlim))
        plt.axvline(rlim, color="purple", label="Rlim 1")
    # Refine the measurement using the s1_down limit (sigma-clipping)
    r = star_profile[1].data["RADIUS"]
    sky = star_profile[1].data["MEDIAN"][int(rlim):]
    sb_nosky_cor = star_profile[1].data["MEDIAN"]
    sb = sb_nosky_cor - bn.nanmedian(sky) #+ 0.5# - np.nanmin(db[1].data["MEAN"])
    noise_level = bn.nanstd(sky)

    #rlim = star_profile[1].data["RADIUS"][np.where((s1_down < 0) & (r > 2))[0][0]] # We put the (r>2) to avoid the saturated core
    rlim = star_profile[1].data["RADIUS"][np.where(((sb - 3*noise_level) < 0) & (r > 2))[0][0]] # We put the (r>2) to avoid the saturated core

    if saturation_level is not None:
        where_is_saturated = np.where(sb_nosky_cor > saturation_level)[0]
        if verbose:
            print("Where is saturated?")
            print(where_is_saturated)
        if len(where_is_saturated) > 1:
            rsat = star_profile[1].data["RADIUS"][where_is_saturated[-1]]
        elif len(where_is_saturated) == 1:
            rsat = star_profile[1].data["RADIUS"][where_is_saturated[0]]
        else:
            rsat = np.nan
    else:
        rsat = np.nan


    if verbose:
        print("Rlim 2: " + str(rlim))
        plt.scatter(r, sb, label="Star profile")
        plt.yscale("log")
        plt.xscale("log")
        if saturation_level is not None:
            plt.axvline(rsat, color="firebrick", label="Saturation level")

        plt.axhline(noise_level, color="orange", label="Noise level 2")
        plt.axvline(rlim, color="red", label="Rlim 2")
        plt.xlabel("R (pixels)")
        plt.ylabel("Intensity (flux)")
        plt.legend()
        plt.show(block=False)

    if clean_profile: rs.utils.execute_cmd("rm " + star_profile_dict["profile_name"])

    return({"rlim": rlim, "rsat": rsat})

############################
def astscript_psf_subtract(input_name, ext, ra, dec, scaling_factor, psf_name, clean=True):
    # Python wrapper of astscript-psf-subtract
    #
    # Based on https://www.gnu.org/software/gnuastro/manual/html_node/Invoking-astscript_002dpsf_002dsubtract.html
    #
    # Shell based example
    # execute_cmd("astscript-psf-subtract mastDownload/HST/j9en0pnaq/j9en0pnaq_flc_ext1.fits \
    #       --mode=wcs \
    #       --psf=PSF/psf_hst_test.fits \
    #       --scale=1.56063e+06 \
    #       --center=329.03478737,1.4031257 \
    #       --output=mastDownload/HST/j9en0pnaq/j9en0pnaq_flc_ext1_sub.fits")

    # First we generate an image with the main WCS, otherwise Noisechisel cannot read the header
    #image_wcs_cleaned_extension_name = rs.mast.make_wcs_clean_extension(input_name, ext)
    output_name = input_name.replace(".fits", "_" + str(ra) + "_" + str(dec) + "_sub.fits")

    # We construct the cmd to send to shell.
    cmd = "astscript-psf-subtract " + input_name +\
           " --mode=wcs --psf="+ psf_name +\
           " --scale="+ str(scaling_factor) + " --center=" + str(ra)+","+str(dec)+\
           " --output=" + output_name

    if not clean:
        print(cmd)
    output_from_shell = rs.utils.execute_cmd(cmd)
    return(output_name)



def scale_and_subtract_stars_single_image(input_name, ext, ra, dec, psf_name, clean=True, extrapolate=False, verbose=False):

    # For now, we will set a fixed norm_radii_min, and norm_radii_max.
    # But this is a function of the saturation level of the star.
    # We need to measure this an scale accordingly.


    input_fits = fits.open(input_name)
    star_model = np.zeros(input_fits[ext].data.shape)
    # First we generate an image with the main WCS, otherwise Noisechisel cannot read the header
    image_wcs_cleaned_extension_name = rs.mast.make_wcs_clean_extension(input_name, ext)
    # We correct the sky of the extension
    rs.sky.correct_flat_sky(input_name=image_wcs_cleaned_extension_name, ext=1, clean=False)

    # We go star by star executing the wrappers.
    for ra_i, dec_i in tqdm(zip(ra, dec), total=len(ra), position=0, leave=True):

        # First we find the scaling factor
        scaling_factor = astscript_psf_scale_factor(input_name=image_wcs_cleaned_extension_name,
                                                    ra=ra_i, dec=dec_i,
                                                    psf_name=psf_name, clean=clean, verbose=verbose)

        # Then we subtract the star
        subtracted_image_name = astscript_psf_subtract(input_name=image_wcs_cleaned_extension_name, ext=1,
                                                       ra=ra_i, dec=dec_i,
                                                       scaling_factor=scaling_factor, psf_name=psf_name, clean=clean)

        # We open the resulting subtracted image
        subtracted_fits = fits.open(subtracted_image_name)

        # And we add what we have subtracted to the star_model
        image_wcs_cleaned_extension = fits.open(image_wcs_cleaned_extension_name)
        star_model = star_model + image_wcs_cleaned_extension[1].data - subtracted_fits[1].data

        # Clean the temporary file
        rs.utils.execute_cmd("rm " + subtracted_image_name)

    # Finally, we save the star model in a new fits.
    output_starmodel_name = input_name.replace(".fits", "_ext" + str(ext) + "_stars.fits")
    output_residual_name = input_name.replace(".fits", "_ext" + str(ext) + "_nostars.fits")

    rs.utils.save_fits(array=star_model, name=output_starmodel_name, header=input_fits[ext].header, extname='STAR_MODEL')
    subtracted_image = input_fits[ext].data - star_model
    rs.utils.save_fits(array=subtracted_image, name=output_residual_name, header=input_fits[ext].header, extname='STAR_SUBTRACTED')

    return({"star_model":output_starmodel_name, "residuals":output_residual_name})
############################


############################
def find_gaia_stars_around_image(lambda_ref, observer=None, input_name=None,
                                 ext=None, ra=None, dec=None, MJD=None, radius=1,
                                 g_mag_max = False,
                                 verbose=False):

    if input_name is not None:
        # First identify where is the image pointing.
        if verbose: print("Reading input image...")
        input_fits = fits.open(input_name)
        w = astropy_wcs.WCS(header=input_fits[ext].header, fobj=input_fits, naxis=2)
        image_shape = input_fits[ext].data.shape
        center = w.all_pix2world(image_shape[0]/2., image_shape[1]/2., 0)
        center_ra = center[0]
        center_dec = center[1]
        gaia_query_filename = input_name.replace(".fits","_gaia.dat")
        if verbose: print("Running exposure inspector...")
        db_inspection_exposure = rs.utils.exposure_inspector(input_name, verbose)
        MJD = db_inspection_exposure["EXPSTART"]
        if observer==None:
            observer = db_inspection_exposure["TELESCOP"]

    else:
        center_ra = ra
        center_dec = dec
        gaia_query_filename = str(ra) + "_" + str(dec) + "_gaia.dat"

        if observer==None:
            print("> Observer frame has not been set. Default is Hubble")
            observer = "Hubble"

    if verbose: print("Running get_hybrid_catalog...")
    gaia_query = get_hybrid_catalog(ra=center_ra, dec=center_dec, radius=radius,
                                    lambda_ref=lambda_ref, MJD=MJD, observer=observer,
                                    g_mag_max = g_mag_max,
                                    verbose=verbose, query_filename=gaia_query_filename)

    gaia_query = gaia_query.sort_values(by=["mag_lambda"], ascending=False)
    gaia_query["cat_id"] = np.linspace(0, len(gaia_query)-1, len(gaia_query), dtype="int64")

    if verbose:
        print("Query saved in " + gaia_query_filename)


    return({"gaia_query": gaia_query, "gaia_query_filename": gaia_query_filename})


############################

def gaia_find_stars_in_and_out(input_name, ext, lambda_ref, ra=None, dec=None, radius=0.1, g_mag_max=False,
                               MJD=None, clean=True, verbose=False, gaia_query=False, match_inside_stars=False):
    # First find the stars using Gaia
    if verbose: print("> Finding Gaia/2MASS/WISE stars around image...")

    if gaia_query is False:
        gaia_query_dict = find_gaia_stars_around_image(lambda_ref=lambda_ref, input_name=input_name, ext=ext,
                                                   ra=ra, dec=dec, MJD=MJD, radius=radius, g_mag_max=g_mag_max, verbose=verbose)
        gaia_query = gaia_query_dict["gaia_query"]



    # Now identify which ones are actually inside the detector
    if verbose: print("> Identifying which stars are inside the FOV and which are outside...")
    input_fits = fits.open(input_name)
    infield_stars = identify_stars_in_out_field(data=input_fits[ext].data,
                                                wcs=astropy_wcs.WCS(input_fits[ext].header, input_fits),
                                                catalog=gaia_query,
                                                verbose=verbose) # identify_stars_in_out_field(input_name=input_name, ext=ext, catalog=gaia_query, ra=ra, dec=dec, verbose=verbose)
    # Lets refine the position of the stars
    # Make a catalog of the objects inside the image
    bool_isIn_gaia_stars = infield_stars["bool_isIn"]

    if (ra is None) and (dec is None):
        if verbose: print("> Generating a catalog of objects in the image...")
        if match_inside_stars:
            object_catalog_in_detector = rs.gnu.astmkcatalog(input_name=input_name, ext=ext, clean=clean, verbose=verbose)
        else:
            object_catalog_in_detector = [False]

        if verbose: print("> Catalog done.")

        # Cross-match the coordinates, and find the real centroids of the stars.
        # We will need to define a magnitude limit, possibly. And a limit in distance.

        # Save a catalog only with the Gaia stars inside the detector
        in_query_filename = input_name.replace(".fits", "_ext" + str(ext) + "_gaia_stars_inside_detector.csv")
        out_query_filename = input_name.replace(".fits", "_ext" + str(ext) + "_gaia_stars_outside_detector.csv")

        gaia_query_in_detector = gaia_query[bool_isIn_gaia_stars]
        gaia_query_in_detector.to_csv(in_query_filename)

        gaia_query_out_detector = gaia_query[np.invert(bool_isIn_gaia_stars)]
        gaia_query_out_detector.to_csv(out_query_filename)
        return({"gaia_query_in_detector": gaia_query_in_detector, "in_query_filename": in_query_filename,
            "gaia_query_out_detector": gaia_query_out_detector, "out_query_filename": out_query_filename,
            "object_catalog_in_detector": object_catalog_in_detector})

    else:
        out_query_filename = str(ra) + "_" + str(dec) + "_gaia_stars_outside_detector.csv"
        gaia_query_out_detector = gaia_query[np.invert(bool_isIn_gaia_stars)]
        gaia_query_out_detector.to_csv(out_query_filename)
        return({"gaia_query_out_detector": gaia_query_out_detector, "out_query_filename": out_query_filename})


############################

def find_stars_inside_detector(input_name, ext, lambda_ref, radius, g_mag_max=False,  MJD=None, clean=True, verbose=False):
    """
    Queries Gaia archive and finds the stars inside a given detector, provided as a FITS file.

    """


    gaia_search_db = gaia_find_stars_in_and_out(input_name=input_name, ext=ext, lambda_ref=lambda_ref,
                                                ra=None, dec=None, radius=radius, g_mag_max=g_mag_max,
                                                MJD=MJD, clean=clean, verbose=verbose, gaia_query=False,
                                                match_inside_stars=True)
    if verbose: print(gaia_search_db)
    object_catalog_in_detector = gaia_search_db["object_catalog_in_detector"]

    matched_catalog_name =  input_name.replace(".fits", "_ext" + str(ext) + "_matched_stars_inside_detector.fits")

    # Reformat the catalog to a format that Gnuastro can read
    inside_star_catalog = pd.read_csv(gaia_search_db["in_query_filename"])
    inside_star_catalog_for_gnuastro = inside_star_catalog[['ra', 'dec', 'mag_lambda']].copy()
    inside_star_catalog_name = gaia_search_db["in_query_filename"].replace(".csv", "_gnuastro.csv")
    inside_star_catalog_for_gnuastro.to_csv(inside_star_catalog_name,  header=False, index=False)


    # Run astmatch and generate a new catalog of matched objects.
    cmd = "astmatch -K  " + inside_star_catalog_name + " -h2 " + object_catalog_in_detector["output_name"] + " --aperture=1/3600 --ccol1=1,2 --ccol2=ra,dec --output " + matched_catalog_name + " --outcols=a1,a2,b1,b2,b3,b4,a3"
    if verbose: print(cmd)
    rs.utils.execute_cmd(cmd)

    # Write the correct table names.
    cmd = "astfits -h1  " + matched_catalog_name +\
           " --update=TTYPE1,RA_GAIA    --update=TUNIT1,deg "+\
           " --update=TTYPE2,DEC_GAIA   --update=TUNIT2,deg "+\
           " --update=TTYPE3,RA_IMA     --update=TUNIT3,deg "+\
           " --update=TTYPE4,DEC_IMA    --update=TUNIT4,deg "+\
           " --update=TTYPE5,X_IMA      --update=TUNIT5,pixel "+\
           " --update=TTYPE6,Y_IMA      --update=TUNIT6,pixel "+\
           " --update=TTYPE7,MAG_LAMBDA --update=TUNIT7,mag "

    if verbose: print(cmd)
    rs.utils.execute_cmd(cmd)

    matched = fits.open(matched_catalog_name)
    matched_catalog_image = matched[1].data


    return({"matched_catalog_image": matched_catalog_image,
            "matched_catalog_name": matched_catalog_name})
############################


# Function to append n empty rows to a DataFrame
def append_empty_rows(dataframe, n):
    for _ in range(n):
        dataframe.loc[len(dataframe)] = pd.Series(dtype='float64')





def identify_stars_in_out_field(data, wcs, catalog, ra=None, dec=None, verbose=False):
    # This program identifies the objects that are inside or outside an image in a fits file.
    # input_name: Name of the fits file
    # ext: Extension of the fits file where we are looking for the image.
    #      i.e., For ACS/WFC in HST, the science frames are stored in ext 1 and 4.
    #      i.e., For WFC/IR, the science frame is in ext 1.
    # catalog: A table containing at least two columns called "ra" and "dec".
    if (ra is not None) and (dec is not None):
        # Plot only stars brighter than N magnitudes
        plt.figure(figsize=(10,10))
        rs.plots.plot_stars_around(catalog)
        print('RA' + str(ra))
        print('DEC' + str(dec))
        plt.scatter(ra, dec, marker="+", color="navy", alpha=0.5, s=500)


        plt.gca().invert_xaxis()
        plt.show(block=False)
        return({"bool_isIn": np.zeros(len(catalog["ra"])).astype("bool"),
                "corners_pix": None,
                "corners_world": None,
                "catalog_inside": None})


    detector_corners = rs.detectors.get_detector_corners(data, wcs)
    corners_pix = detector_corners["corners_pix"]
    corners_world = detector_corners["corners_world"]

    ######### Check if something is infield or outfield #########
    if verbose: print(" Detector corners world: " + str(corners_world))
    if verbose: print(" Detector corners pixel: " + str(corners_pix))

    vertices_detector_dec_ra = np.zeros(corners_world.shape)
    vertices_detector_dec_ra[:,1] = corners_world[:,0]
    vertices_detector_dec_ra[:,0] = corners_world[:,1]

    from sphericalpolygon import Sphericalpolygon
    #print("DEMO WARNING: This can be done way faster by prefiltering those stars at high angular distances.")
    region_detector = Sphericalpolygon.from_array(vertices_detector_dec_ra)

    distance_detector_corners = rs.utils.angular_distance(ra1=np.array(corners_world[:,0]),
                                                          dec1=np.array(corners_world[:,1]),
                                                          ra2=np.array(corners_world[:,0]),
                                                          dec2=np.array(corners_world[:,1]))

    aprox_size_detector = bn.nanmax(distance_detector_corners)

    distance_detector_to_stars = bn.nanmax(rs.utils.angular_distance(ra1=np.array(corners_world[:,0]),
                                                                     dec1=np.array(corners_world[:,1]),
                                                                     ra2=catalog["ra"],
                                                                     dec2=catalog["dec"]), axis=0)

    isIn = np.zeros(len(catalog["ra"])).astype("bool")
    far_away_stars = distance_detector_to_stars > 5*aprox_size_detector
    isIn[far_away_stars] = False

    close_stars_around_list_dec_ra = np.array([[j, i] for i, j in zip(catalog["ra"], catalog["dec"])])

    isIn[~far_away_stars] = region_detector.contains_points(close_stars_around_list_dec_ra[~far_away_stars])

    catalog_inside = catalog[isIn]
    if verbose: print(" Stars inside detector: " + str(catalog_inside))
    if verbose: print(" All stars catalog: " + str(catalog))


    return({"bool_isIn": isIn, "corners_pix": corners_pix, "corners_world": corners_world, "catalog_inside": catalog_inside})


def get_hybrid_catalog(ra, dec, radius, lambda_ref, MJD, observer, g_mag_max=False, verbose=False, query_filename="default_query.dat"):
    from scipy import interpolate

    # C is the central coordinates of the pointing.
    c = SkyCoord(ra, dec, frame='icrs', unit="deg")

    # First get the high resolution catalog.
    if verbose: print("Running query_gaia_2mass_wise query...")
    high_resolution_catalog = rs.gaia.query_gaia_2mass_wise(ra=ra, dec=dec,
                                                    radius=radius,
                                                    g_mag_max=g_mag_max,
                                                    verbose=verbose,
                                                    query_filename=query_filename)
    # Load the low resolution catalog.
    low_resolution_catalog = pd.read_csv(os.path.dirname(rs.utils.__file__) + "/CORE/gaia_2mass_wise_healpix_7_cat.csv")

    # Make a new low_resolution_catalog where the columns will fit.
    # Maybe stick this in the low_res
    # ----------------------------------------------- #

    # Remove the super stars in the high res region
    ipix_highres = sorted(set(high_resolution_catalog["healpix_lvl"]))
    if verbose: print("Healpix cells to remove from low resolution:" + str(ipix_highres))

    row_healpix_to_drop_list = []
    for healpix_covered_in_highres_area in ipix_highres:
        if verbose: print(healpix_covered_in_highres_area)
        row_healpix_to_drop = np.where(low_resolution_catalog["healpix_lvl"] == healpix_covered_in_highres_area)[0]
        if verbose: print(row_healpix_to_drop)
        if verbose: print(low_resolution_catalog["healpix_lvl"].iloc[row_healpix_to_drop])
        row_healpix_to_drop_list.append(row_healpix_to_drop[0])
        # Before dropping, check if the sum of fluxes is the same.
        #print("High res g flux: " + str(bn.nansum(high_resolution_catalog["g"].iloc[np.where(high_resolution_catalog["healpix_5"] == healpix_covered_in_highres_area)])))
        #print("Low res g flux: " + str(bn.nansum(low_resolution_catalog["g"].iloc[row_healpix_to_drop])))
        #print("High res j flux: " + str(bn.nansum(high_resolution_catalog["j"].iloc[np.where(high_resolution_catalog["healpix_5"] == healpix_covered_in_highres_area)])))
        #print("Low res j flux: " + str(bn.nansum(low_resolution_catalog["j"].iloc[row_healpix_to_drop])))
        #print("High res w1 flux: " + str(bn.nansum(high_resolution_catalog["w1"].iloc[np.where(high_resolution_catalog["healpix_5"] == healpix_covered_in_highres_area)])))
        #print("Low res w1 flux: " + str(bn.nansum(low_resolution_catalog["w1"].iloc[row_healpix_to_drop])))

    low_resolution_catalog = low_resolution_catalog.drop(row_healpix_to_drop_list)

    #high_resolution_catalog.append(low_resolution_catalog)
    high_resolution_catalog = pd.concat([high_resolution_catalog, low_resolution_catalog], axis=0)

    ######################################################################
    #
    # Now complement the query with the 230 naked eye stars from Hipparcos
    #
    ######################################################################
    gaia_naked_eye_stars = pd.read_csv(os.path.dirname(rs.utils.__file__) + "/CORE/gaia_naked_eye_stars.csv", delim_whitespace=True)

    # Find which ones are in the query radius.
    ra_stars = np.array(gaia_naked_eye_stars["ra"])*u.deg
    dec_stars = np.array(gaia_naked_eye_stars["dec"])*u.deg
    gaia_coords = SkyCoord(ra_stars, dec_stars, frame='icrs')
    sep = c.separation(gaia_coords)

    n_bright_stars = len(gaia_naked_eye_stars["ra"])
    n_highres_stars = len(high_resolution_catalog)

    print("\n[DEMO WARNING:] Naked eye stars do not have photometry besides g band")
    #append_empty_rows(high_resolution_catalog, n_bright_stars)
    if verbose: print("> Adding bright stars from Hipparcos survey...")

    source_id_bright_star = []
    ra_bright_star = []
    dec_bright_star = []
    bp_bright_star = []
    g_bright_star = []
    rp_bright_star = []
    j_bright_star = []
    h_bright_star = []
    ks_bright_star = []
    w1_bright_star = []
    w2_bright_star = []
    w3_bright_star = []
    w4_bright_star = []

    phot_bp_mean_mag_AB_bright_star = []
    phot_g_mean_mag_AB_bright_star = []
    phot_rp_mean_mag_AB_bright_star = []
    phot_j_mean_mag_AB_bright_star = []
    phot_h_mean_mag_AB_bright_star = []
    phot_ks_mean_mag_AB_bright_star = []
    phot_w1_mean_mag_AB_bright_star = []
    phot_w2_mean_mag_AB_bright_star = []
    phot_w3_mean_mag_AB_bright_star = []
    phot_w4_mean_mag_AB_bright_star = []

    for i in tqdm(range(len(gaia_naked_eye_stars)), disable=not verbose, position=0, leave=True):
        source_id_bright_star.append("Hipparcos bright star")
        ra_bright_star.append(gaia_naked_eye_stars["ra"][i])
        dec_bright_star.append(gaia_naked_eye_stars["dec"][i])
        bp_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_bp_ab"][i])))
        g_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_g_ab"][i])))
        rp_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_rp_ab"][i])))
        j_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_j_ab"][i])))
        h_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_h_ab"][i])))
        ks_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_ks_ab"][i])))
        w1_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_w1_ab"][i])))
        w2_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_w2_ab"][i])))
        w3_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_w3_ab"][i])))
        w4_bright_star.append(10**(0.4*(8.9-gaia_naked_eye_stars["mag_w4_ab"][i])))

        phot_bp_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_bp_ab"][i])
        phot_g_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_g_ab"][i])
        phot_rp_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_rp_ab"][i])
        phot_j_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_j_ab"][i])
        phot_h_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_h_ab"][i])
        phot_ks_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_ks_ab"][i])
        phot_w1_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_w1_ab"][i])
        phot_w2_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_w2_ab"][i])
        phot_w3_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_w3_ab"][i])
        phot_w4_mean_mag_AB_bright_star.append(gaia_naked_eye_stars["mag_w4_ab"][i])

    bright_star_catalog = pd.DataFrame({"source_id": source_id_bright_star,
                           "ra": ra_bright_star,
                           "dec": dec_bright_star,
                           "bp": bp_bright_star, "phot_bp_mean_mag_AB": phot_bp_mean_mag_AB_bright_star,
                           "g": g_bright_star, "phot_g_mean_mag_AB": phot_g_mean_mag_AB_bright_star,
                           "rp": rp_bright_star, "phot_rp_mean_mag_AB": phot_rp_mean_mag_AB_bright_star,
                           "j": j_bright_star, "phot_j_mean_mag_AB": phot_j_mean_mag_AB_bright_star,
                           "h": h_bright_star, "phot_h_mean_mag_AB": phot_h_mean_mag_AB_bright_star,
                           "ks": ks_bright_star, "phot_ks_mean_mag_AB": phot_ks_mean_mag_AB_bright_star,
                           "w1": w1_bright_star, "phot_w1_mean_mag_AB": phot_w1_mean_mag_AB_bright_star,
                           "w2": w2_bright_star, "phot_w2_mean_mag_AB": phot_w2_mean_mag_AB_bright_star,
                           "w3": w3_bright_star, "phot_w3_mean_mag_AB": phot_w3_mean_mag_AB_bright_star,
                           "w4": w4_bright_star, "phot_w4_mean_mag_AB": phot_w4_mean_mag_AB_bright_star})

    high_resolution_catalog = pd.concat([high_resolution_catalog, bright_star_catalog], axis=0)

    # > Here we generate a new column with the custom user magnitude
    mag_lambda = np.zeros(len(high_resolution_catalog))
    filters_list_lambda = np.array([rs.telescopes.Passbands.Gaia.bp.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.Gaia.g.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.Gaia.rp.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.TwoMASS.J.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.TwoMASS.H.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.TwoMASS.Ks.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.WISE.W1.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.WISE.W2.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.WISE.W3.lambda_ref.to("AA").value,
                                    rs.telescopes.Passbands.WISE.W4.lambda_ref.to("AA").value])

    magnitude_array = np.array([high_resolution_catalog["phot_bp_mean_mag_AB"],
                                high_resolution_catalog["phot_g_mean_mag_AB"],
                                high_resolution_catalog["phot_rp_mean_mag_AB"],
                                high_resolution_catalog["phot_j_mean_mag_AB"],
                                high_resolution_catalog["phot_h_mean_mag_AB"],
                                high_resolution_catalog["phot_ks_mean_mag_AB"],
                                high_resolution_catalog["phot_w1_mean_mag_AB"],
                                high_resolution_catalog["phot_w2_mean_mag_AB"],
                                high_resolution_catalog["phot_w3_mean_mag_AB"],
                                high_resolution_catalog["phot_w4_mean_mag_AB"]])

    # Transform the magnitudes to the observing filter.
    #print("Demo warning: Using Gaia g filter magnitudes instead of the correct filter")
    #print("Query Gaia/2MASS/WISE https://notebook.community/bosscha/alma-calibrator/notebooks/2mass/13_environment-with_analysis")

    for i in tqdm(range(len(high_resolution_catalog)), disable=not verbose, position=0, leave=True):
        good = np.where(np.isfinite(magnitude_array[:,i]))
        #print(magnitude_array[:,i])
        magnitude_interpolator = interpolate.interp1d(x=filters_list_lambda[good], y=magnitude_array[good,i], kind="nearest", fill_value="extrapolate")
        mag_lambda[i] = magnitude_interpolator(lambda_ref.to("AA").value)

    high_resolution_catalog["mag_lambda"] = mag_lambda


    # magnitude_interpolator = interpolate.interp1d(x=mag_sun_db["lambda_mum"]*10**4, y=mag_sun_db["AB_App"])
    # mag_ref = magnitude_interpolator(lambda_ref.to("AA"))


    # Here, calculate the positions of the planets in L2.
    # For the filter considered, add the magnitude of the planets.
    planets_catalog = rs.sso.get_SS0s_loc_magnitude(observer=observer, MJD=MJD, lambda_ref=lambda_ref, verbose=verbose) # Awesome it is working
    if verbose: print(planets_catalog)
    planets_catalog_to_concat = pd.DataFrame({"source_id": planets_catalog["source_id"],
                                              #"source_id_2": planets_catalog["source_id_2"],
                                              "ra": planets_catalog["RA"],
                                              "dec": planets_catalog["DEC"],
                                              "mag_lambda": planets_catalog["mag_custom"]})
    high_resolution_catalog = pd.concat([high_resolution_catalog, planets_catalog_to_concat], axis=0)
    #print("DEMO WARNING: ADD THE MOON AND REST OF SSOs to psf.get_hybrid_catalog")
    #print("https://arxiv.org/pdf/1808.01973.pdf")

    # Find the distance to the pointing.
    ra_stars = np.array(high_resolution_catalog["ra"])*u.deg
    dec_stars = np.array(high_resolution_catalog["dec"])*u.deg
    star_coords = SkyCoord(ra_stars, dec_stars, frame='icrs')
    sep = c.separation(star_coords)
    high_resolution_catalog["dist"] = sep.deg

    high_resolution_catalog.to_csv(query_filename)

    return(high_resolution_catalog)





def f_hst_attenuation(theta):
    return(10**f_hst_attenuation_interpolator(np.log10(theta)))


##############



def psf_harvester_single_exposure(input_name, clean=True, verbose=False):
    # if True:
    import astropy.units as u
    from tqdm import tqdm

    if verbose: print("Scanning input exposure: " + input_name)
    exposure_identity = rs.utils.exposure_inspector(input_name=input_name)
    if verbose: print("Done.")


    max_size_of_a_star = 501
    g_mag_max = 30


    main_folder = os.getcwd() # os.path.dirname(input_name)

    sci_exts = exposure_identity["SCIEXTS"] # To be detected automatically with exposure inspector

    # For each SCIENCE extension we run the following loop.


    for ext in sci_exts:
        ###########################
        # Fundamental estimations #
        ###########################
        data, wcs = rs.utils.get_data_and_wcs(input_name=input_name, ext=ext)
        star_search_radius = rs.utils.find_max_angular_size_of_image(data, wcs)
        ###########################

        # Make a segmentation map to find out the saturated sources
        segmentation_map = rs.gnu.make_gnuastro_segmentation_map(input_name, ext, saturation_level=50000)


        # Identify the stars
        # At this point, we can have the approximate coordinates of the stars.
        # But in real images, the centers might differ.
        # find_stars_inside_detector cross-match the stars in the image and the ones in the Gaia catalog
        # to give a better position of the centroids.

        matched_stars_in_detector = rs.psf.find_stars_inside_detector(input_name=input_name,
                                                                      ext=ext,
                                                                      lambda_ref=exposure_identity["PHOTPLAM"]*u.AA,
                                                                      radius=star_search_radius,
                                                                      g_mag_max=g_mag_max,
                                                                      MJD=exposure_identity["EXPSTART"],
                                                                      clean=clean,
                                                                      verbose=verbose)
        if verbose: print(matched_stars_in_detector)
        ra_gaia_list = matched_stars_in_detector["matched_catalog_image"]["RA_GAIA"]
        dec_gaia_list = matched_stars_in_detector["matched_catalog_image"]["DEC_GAIA"]
        ra_ima_list = matched_stars_in_detector["matched_catalog_image"]["RA_IMA"]
        dec_ima_list = matched_stars_in_detector["matched_catalog_image"]["DEC_IMA"]
        x_ima_list = matched_stars_in_detector["matched_catalog_image"]["X_IMA"]
        y_ima_list = matched_stars_in_detector["matched_catalog_image"]["Y_IMA"]
        mag_list = matched_stars_in_detector["matched_catalog_image"]["MAG_LAMBDA"]

        sci_ima_name     = rs.gnu.make_wcs_clean_extension(input_name, ext)
        sky_corrected_db = rs.sky.correct_flat_sky(input_name=sci_ima_name, ext=1, clean=clean, verbose=verbose)

        basename = os.path.basename(sci_ima_name)
        psf_name = sci_ima_name.replace(".fits", "_psf.fits")

        ################################
        # Now we analyze star by star.
        ################################

        if not os.path.exists(main_folder + "/psf_harvester/"):
            stdout = os.mkdir(main_folder + "/psf_harvester/")
        if not os.path.exists(main_folder + "/psf_harvester/psf_level_0"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_0")
        if not os.path.exists(main_folder + "/psf_harvester/psf_level_0/stamps"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_0/stamps")

        if not os.path.exists(main_folder + "/psf_harvester/psf_level_1"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_1")
        if not os.path.exists(main_folder + "/psf_harvester/psf_level_1/stamps"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_1/stamps")

        if not os.path.exists(main_folder + "/psf_harvester/psf_level_2"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_2")
        if not os.path.exists(main_folder + "/psf_harvester/psf_level_2/stamps"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_2/stamps")

        if not os.path.exists(main_folder + "/psf_harvester/psf_level_3"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_3")
        if not os.path.exists(main_folder + "/psf_harvester/psf_level_3/stamps"):
            stdout = os.mkdir(main_folder + "/psf_harvester/psf_level_3/stamps")

        outname_list = []
        rlim = np.zeros(len(ra_ima_list))
        rsat = np.zeros(len(ra_ima_list))
        norm_radii_min = np.zeros(len(ra_ima_list))
        norm_radii_max = np.zeros(len(ra_ima_list))
        psf_level = np.zeros(len(ra_ima_list))

        nstars = len(ra_ima_list)


        for i in tqdm(range(nstars)):
            ra = ra_ima_list[i]
            dec = dec_ima_list[i]
            x_ima = x_ima_list[i]
            y_ima = y_ima_list[i]
            mag = mag_list[i]

            # Measuring the normalization radius for each star
            measure_maxradii_result = rs.psf.measure_maxradii(input_name=sci_ima_name, ra=ra, dec=dec,
                                           rmax=max_size_of_a_star, clean=clean,
                                           verbose=verbose, saturation_level=50000, clean_profile=clean)

            rlim[i] = measure_maxradii_result["rlim"]
            rsat[i] = measure_maxradii_result["rsat"]


            if rlim[i] is None:
                rlim[i] = np.nan
                rsat[i] = np.nan
                norm_radii_min[i] = np.nan
                norm_radii_max[i] = np.nan
                outname_list.append(None)

            else: # If the star is found, make a stamp.
                # Three different saturation levels.
                # This has to remain flexible. We will add more saturation levels
                # if we find brigther stars.
                if rlim[i] < 10: # If we do not detect more than 10 pixels, then the star is not detected.
                    psf_level[i] = int(0)

                elif np.isnan(rsat[i]):
                    psf_level[i] = int(1)
                    norm_radii_min[i], norm_radii_max[i] = 5, 8

                elif (rsat[i] >= 0) and (rsat[i] < 5):
                    psf_level[i] = int(2)
                    norm_radii_min[i], norm_radii_max[i] = 10, 15

                elif (rsat[i] >= 5):
                    psf_level[i] = int(3)
                    norm_radii_min[i], norm_radii_max[i] = 20, 30



                outname = main_folder + "/psf_harvester/psf_level_" + str(int(psf_level[i])) + "/stamps/"+basename.replace(".fits","")+"_"+str(i).zfill(5)+"_mag"+str(np.round(mag,2))+".fits"
                outname_list.append(outname)

                cmd = "astscript-psf-stamp "+ sci_ima_name + " " +\
                      "--mode=img --normradii="+str(norm_radii_min[i])+","+str(norm_radii_max[i])+\
                      " --center="+str(x_ima)+","+str(y_ima)+" "+\
                      " --widthinpix="+str(max_size_of_a_star)+","+str(max_size_of_a_star)+" "+\
                      " --segment="+ segmentation_map + " " +\
                      " --output="+ outname
                stdout = rs.utils.execute_cmd(cmd)

                # Write the info on the header of the psf stamp
                cmd = "astfits -h1  " + outname+\
                      " --update=EXPOSURE,"+str(input_name)+\
                      " --update=EXT,"+str(ext)+\
                      " --update=RA_GAIA,"+str(ra_gaia_list[i])+\
                      " --update=DEC_GAIA,"+str(dec_gaia_list[i])+\
                      " --update=RA_IMA,"+str(ra_ima_list[i])+\
                      " --update=DEC_IMA,"+str(dec_ima_list[i])+\
                      " --update=X_IMA,"+str(x_ima_list[i])+\
                      " --update=Y_IMA,"+str(y_ima_list[i])+\
                      " --update=MAG_WAVE,"+str(mag_list[i])+\
                      " --update=RLIM,"+str(rlim[i])+\
                      " --update=RSAT,"+str(rsat[i])+\
                      " --update=PSF_LVL,"+str(psf_level[i])+\
                      " --update=NRM_RMIN,"+str(norm_radii_min[i])+\
                      " --update=NRM_RMAX,"+str(norm_radii_max[i])+\
                      " --update=FILTER,"+str(exposure_identity["FILTER"])+\
                      " --update=EXPTIME,"+str(exposure_identity["EXPTIME"])

                stdout = rs.utils.execute_cmd(cmd)

        if clean:
            input_directory = os.path.dirname(input_name)
            rs.utils.execute_cmd("rm " + input_directory + "/*_ext*.fits")
            rs.utils.execute_cmd("rm " + input_directory + "/*_gaia.dat")
            rs.utils.execute_cmd("rm " + input_directory + "/*_swarp_coadd.fits")
