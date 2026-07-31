def f_hst_attenuation(theta):
    return(10**f_hst_attenuation_interpolator(np.log10(theta)))

##############


def sort_hst_flcs_by_filter(filelist):
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


def measure_sky_level_HST_ACS(exposure_name, verbose=False):
    # Open the fits file
    from astropy.io import fits

    exposure_fits = fits.open(exposure_name) # We open the fits file with astropy

    # Try to retrieve a previous analysis that was performed in this exposure
    try:
        return({"zody": exposure_fits[1].header["ZODY"],
                "sky":  exposure_fits[1].header["SKYLVL"]})

    except:
        if verbose: print("No previous sky-zody analysis available in this header.")


    # Get the basic parameters
    PHOTPLAM =  exposure_fits[1].header["PHOTPLAM"] # / Pivot wavelength (Angstroms)   <---- This is the observation wavelength. We will use this to estimate the Zodiacal light.
    PHOTFLAM =  exposure_fits[1].header["PHOTFLAM"] # / inverse sensitivity, ergs/cm2/Ang/electron   <---- This is the transformation between electrons to flux. We will use this to turn the e/s/px to Jy/arcsec2.
    RA_TARG  = exposure_fits[0].header["RA_TARG"]    #  Right ascension of the observation target .
    DEC_TARG = exposure_fits[0].header["DEC_TARG"]  #  Declination of the observation target
    EXPSTART = exposure_fits[0].header["EXPSTART"]  #  Time of the observation, in Modified Julian Date
    # Get the right filter
    filter_1 = exposure_fits[0].header["FILTER1"]
    filter_2 = exposure_fits[0].header["FILTER2"]

    if "CLEAR" in filter_1:
        filter_name = filter_2
    else:
        filter_name = filter_1


    # Lets import some more astropy packages to deal with the time, units, and coordinates.

    # -------------------------------------- #
    # MODEL ZODIACAL LIGHT                   #
    # -------------------------------------- #
    # Initialize a zodiacal light model at a wavelength/frequency or over a bandpass
    zody_background_sci1 = rs.sky.get_zodiacal_background(input_name=exposure_name, ext=1,
                                                          wavelength=filter_name,
                                                          telescope="HST",
                                                          instrument="ACS",
                                                          detector="WFC",
                                                          expstart=EXPSTART,
                                                          step=4000, zody_mode="zodipy",
                                                          nbins_wavelength=20, obslocin=3,
                                                          grid_method="random", verbose=False, interpolate=False)
    import bottleneck as bn
    median_zody_jy_arcsec2 = bn.nanmedian(zody_background_sci1)

    sky_sci1 = rs.sky.correct_flat_sky(input_name=exposure_name, ext=1, overwrite=False, clean=True, verbose=False)
    median_sky_es = sky_sci1["skylvl"]
    # print(median_sky_es)
    median_sky_jy_arcsec2 = rs.detectors.HST_ACS_counts_to_jy(flux_ACS=median_sky_es, photflam=PHOTFLAM, photplam=PHOTPLAM)

    if not np.isnan(median_sky_jy_arcsec2):    exposure_fits[1].header["SKYLVL"] = median_sky_jy_arcsec2
    if not np.isnan(median_zody_jy_arcsec2):    exposure_fits[1].header["ZODY"]   = median_zody_jy_arcsec2

    exposure_fits.verify("silentfix")
    exposure_fits.writeto(exposure_name, overwrite=True)
    exposure_fits.close()

    return({"zody": median_zody_jy_arcsec2, "sky": median_sky_jy_arcsec2})
