import os
#import sep Problem with numpy 2 .
from tqdm import tqdm
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import rosalia as rs

def clean_crs(input_list):
    # Use astrodrizzle for cleaning the Cosmic Rays.
    # This works OK for HST ACS.
    # TODO: We need to find a way to generalize this.
    #
    # Advice from Mohammad Akhlaghi
    # hdu=1
    # in=j9en0pnaq-h1-det.fits
    # out=j9en0pnaq-h1-det-nocs.fits
    # astarithmetic $in -h$hdu set-i i i 0 eq nan where set-i \
    #                5 5 i filter-median set-f i i f - f / 5 gt \
    #                2 dilate nan where --output=$out
    ##################################################
    from drizzlepac import astrodrizzle
    from drizzlepac import tweakreg
#
    astrodrizzle.AstroDrizzle(input=list(input_list), skysub=1, skymethod='globalmin+match',
                                                  static=0,
                                                  in_memory=1, preserve=1, restore=0,
                                                  clean=1,
                                                  driz_sep_pixfrac=1.,
                                                  num_cores=1,
                                                  driz_sep_wcs=1,
                                                  driz_cr_corr=1,
                                                  resetbits=128,
                                                  driz_cr_snr='4.0 3.0',
                                                  driz_cr_scale='1.5 0.7',
                                                  driz_combine=0, updatewcs=0,
                                                  combine_type="imedian")

    # Get the crclean files
    crclean_list = [i.replace("_flc.fits", "_crclean.fits") for i in input_list]
    print(crclean_list)
    return(crclean_list)


def source_extractor(input_fits, ext, sn=3):
    # This program runs SEP and generates a catalog of sources.
    data = fits.open(input_fits)[ext].data.byteswap(False).newbyteorder()
    bkg = sep.Background(data)
    # Subtract the background
    data_sub = data - bkg
    # Run SEP on the subtracted data
    objects = sep.extract(data_sub, sn, err=bkg.globalrms)
    catalog = pd.DataFrame(objects)
    ra, dec = rs.utils.xy_to_radec(x=catalog["x"], y=catalog["y"], fits_name=input_fits, ext=ext)

    # We define a new ASCII Astropy Table
    data = Table()
    data['x'] = catalog["x"]
    data['y'] = catalog["y"]
    data['ra'] = ra
    data['dec'] = dec
    data['flux'] = catalog["flux"]


    catfile_name = input_fits.replace(".fits", "_ext" + str(ext) + "_cat.csv")
    #catalog.to_csv(catfile_name, sep=' ')
    ascii.write(data, catfile_name, overwrite=True, format="commented_header")

    x_col_number = 1
    y_col_number = 2
    ra_col_number = 3
    dec_col_number = 4
    flux_col_number = 5

    return({"CATALOG": data, "CATFILE_NAME": catfile_name,
            "X_COL": x_col_number, "Y_COL": y_col_number,
            "RA_COL": ra_col_number, "DEC_COL": dec_col_number,
            "FLUX_COL": flux_col_number})



def generate_tweakreg_catfile(input_list, catfile_name = "catfile_default.txt"):

    if not isinstance(input_list, (list, pd.core.series.Series, np.ndarray)):
        input_list = np.array([input_list])

    with open(catfile_name, 'w') as catfile:
        for input_name in tqdm(input_list):
            # Make the catfile of an HST image
            input_inspection = rs.mast.HST_inspector(input_name)
            if (input_inspection["INSTRUME"] == "ACS"):
                sci_exts = [1, 4]
            if (input_inspection["INSTRUME"] == "WFC3IR"):
                sci_exts = [1]

            catfile_line = input_name
            for sci_ext_i in sci_exts:
                input_name_ext_i_SEP = source_extractor(input_fits=input_name, ext=sci_ext_i)
                refcat_i = input_name_ext_i_SEP["CATFILE_NAME"]
                catfile_line = catfile_line + " " + refcat_i
                refxcol = input_name_ext_i_SEP["X_COL"]
                refycol = input_name_ext_i_SEP["Y_COL"]
                refracol = input_name_ext_i_SEP["RA_COL"]
                refdeccol = input_name_ext_i_SEP["DEC_COL"]
                reffluxcol = input_name_ext_i_SEP["FLUX_COL"]
            catfile.write(catfile_line + '\n')

    return({"CATFILE_NAME": catfile_name,
            "REFXCOL": refxcol, "REFYCOL": refycol,
            "REFRACOL": refracol, "REFDECCOL": refdeccol, "FLUX_COL": reffluxcol})



def get_parameters_list(fits_list, index, ext=0):
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


def separate_visits(fitslist):
    PARAMS = get_parameters_list(fitslist, ["ASN_ID"])
    unique_visits = np.array(list(set(PARAMS[0])))
    flocks = []
    fits_array = np.array(fitslist)
    for i in unique_visits:
        print(i)
        flock = fits_array[np.array(PARAMS[0]) == i]
        print(flock)
        flocks.append(list(flock))
    return(flocks)



def astrometry_list(fitslist, refcat=None, refxcol=None, refycol=None, force_tweakreg=False, correct_crs=False):
    # This program runs on two stages.
    # First - Removing the CRs.
    # Second - Aligning the images if necesary.
    ##################################

    # First we separate the observations in visits.
    # That ensures that the differences in astrometry will not be TOO high
    # before CR correction.

    for fitsimages in fitslist:
        rootname = os.path.basename(fitsimages).split("_")[0]
        os.system("rm " + rootname + "_crclean.fits")
    chunks_final = separate_visits(fitslist)

    print(chunks_final)
    #####################
    # 1 - This is where we remove the Cosmic Rays.
    # Required for running Tweakreg.

    crclean_list = []
    for subset in chunks_final:
        if len(subset) > 20:
            step1, step2 = split_list(subset)
            steps_list = [step1, step2]
        else:
            steps_list = [subset]

        for step_i in steps_list:
            crclean_step_i = clean_crs(step_i)
            for crclean_image in crclean_step_i:
                crclean_list.append(crclean_image)

    #####################
    # 2 - This is where we generate the master catalog for astrometry.
    #     If it is not provided.

    do_tweakreg = True

    if refcat is None:
        # Find the WCSNAME of the input files.
        # Select one of those with GAIA WCS name.
        # If it does not exists, select the one with the largest EXPTIME.
        WCSNAME_list = get_parameters_list(fits_list=fitslist, index=["WCSNAME"], ext=1)[0]
        bool_GAIA_WCSLIST = []
        for WCSNAME_i in WCSNAME_list:
            bool_GAIA_WCSLIST.append("GAIA" in WCSNAME_i)

        if all(bool_GAIA_WCSLIST):
            print("All images contain a GAIA MAST revised WCS")
            print("Astrometric calibration not required. Jumping")
            do_tweakreg = False # Change this to avoid forcing the astrometry

        if force_tweakreg:
            print("force_tweakreg set to True")
            print("Forcing astrometric calibration with Tweakreg")
            do_tweakreg = True # Change this to avoid forcing the astrometry


        #else:
            # Find a good image to act as master frame.
            # TO BE DONE: Now select the first one
            reference_frame = crclean_list[0]
            print(reference_frame)
            reference_frame_drz = rs.utils.run_basic_astrodrizzle(reference_frame)
            reference_frame_SEP = source_extractor(input_fits=reference_frame_drz, ext=0)
            refcat = reference_frame_SEP["CATFILE_NAME"]
            refxcol = reference_frame_SEP["RA_COL"]
            refycol = reference_frame_SEP["DEC_COL"]
            rfluxcol = reference_frame_SEP["FLUX_COL"]
    else:
        refcat = refcat
        refxcol = refxcol
        refycol = refycol

    if do_tweakreg:
        print("\n\n\n\n\n Tweakreg stage! Now aligning \n\n\n\n")
        print(refcat)

        print("Generating SEP catalogs for automatic alignment")
        catfile_dict = generate_tweakreg_catfile(input_list=crclean_list, catfile_name = "catfile_default.txt")
        catfile_name = catfile_dict["CATFILE_NAME"]
        xcol = catfile_dict["REFXCOL"]
        ycol = catfile_dict["REFYCOL"]
        fluxcol = catfile_dict["FLUX_COL"]
        from drizzlepac import tweakreg
        tweakreg.TweakReg(files=crclean_list,
                          updatewcs=False, updatehdr=True,
                          verbose=True, searchrad=1, searchunits="arcseconds", headerlet=True, clobber=True,
                          catfile=catfile_name, xcol=xcol, ycol=ycol, xyunits="pixels", fluxcol=fluxcol, fluxunits="cps",
                          refimage=reference_frame_drz, refcat=refcat, refxcol=refxcol, refycol=refycol,
                          refxyunits="degrees",  rfluxcol=rfluxcol, rfluxunits = "cps",
                          interactive=0, use2dhist=True, see2dplot=True, fitgeometry="general")

    print("Copying astormetry to original files...")
    for file_to_write, reference_file in zip(fitslist, crclean_list):
        print("Header " + reference_file + " -> " + file_to_write)
        copy_header(file_to_write, reference_file, copy_data=correct_crs)
        print("Done")

    return(fitslist)


def copy_header(file_to_write, reference_file, copy_data=False):
    file_to_write_fits = fits.open(file_to_write)
    reference_file_fits = fits.open(reference_file)

    for i in range(len(file_to_write_fits)):
        file_to_write_fits[i].header = reference_file_fits[i].header
        if copy_data:
            file_to_write_fits[i].data = reference_file_fits[i].data

    file_to_write_fits.verify("silentfix")
    file_to_write_fits.writeto(file_to_write, overwrite=True)
    file_to_write_fits.close()
    reference_file_fits.close()
