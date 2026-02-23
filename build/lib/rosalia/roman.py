import os
import numpy as np
from tqdm import tqdm
import rosalia as rs
import coord

from astropy.time import Time
from astropy.io import fits
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy import constants as const
from datetime import datetime


# --------------------- #
def create_roman_dummy(point, date, band, PA=None, exptime=500, output="default_roman_dummy.fits"):
    # create_roman_dummy:
    # Uses GalSim subroutines to generate a dummy empty Roman/WFI FITS file with the right WCS and array sizes.
    # Input:
    # pointing = Astropy coord object with ra, dec.
    #            Example c = SkyCoord(1, -30, frame="icrs", unit="deg")
    #
    # date = Astropy time object with the exposure date.
    #            Example t = Time(56000, format='mjd', scale='utc')
    #            Example 2: t = Time('1999-01-01T00:00:00.123456789', format='isot', scale='utc')
    # band = String. Name of the filter - http://svo2.cab.inta-csic.es/theory/fps/
    #
    # PA = Position angle of the telescope Y axis in degrees. Note that WFI is rotated 30 degrees from the Y direction,
    #      so PA = 60 will give an image approximately aligned with the North axis. If left None, Galsim will find the optimal
    #      PA for the telescope on that date.
    #
    # output = Name of the output image.
    # --------------------------------------------- #
    import galsim
    import galsim.roman as roman

    # Setting up the parameters in GalSim format
    ra_targ = coord.Angle(point.ra.degree*coord.degrees)
    dec_targ = coord.Angle(point.dec.degree*coord.degrees)
    targ_pos = galsim.CelestialCoord(ra=ra_targ, dec=dec_targ)

    if PA is not None:
        WFI_PA = galsim.Angle(PA, coord.degrees)
    else:
        WFI_PA = None

    # Starting the loop for the 18 SCAs

    temp_out_sca_filename_list = []
    temp_array_list = []
    temp_header_list = []
    temp_sca_name_list = []

    for SCA_i in range(18):

        SCA_name = str(SCA_i+1).zfill(3)
        #print(SCA_name)

        temp_sca_name_list.append(SCA_name)
        out_sca_filename = output.replace(".fits", "_sca" + SCA_name + ".fits")
        temp_out_sca_filename_list.append(out_sca_filename)

        wcs_dict = roman.getWCS(world_pos=targ_pos, PA=WFI_PA, SCAs=SCA_i+1, date=date.tt.datetime) # world_pos, PA=None, date=None, SCAs=None, PA_is_FPA=False)
        wcs = wcs_dict[SCA_i+1]

        # Set up the full image for the galaxies
        full_image = galsim.ImageF(roman.n_pix, roman.n_pix, wcs=wcs)
        full_image.write(out_sca_filename)

        # Open the temporary SCA image
        temp_sca = fits.open(out_sca_filename)

        temp_array_list.append(temp_sca[0].data)

        temp_sca[0].header["EXTNAME"] = "SCA" + SCA_name


        os.remove(out_sca_filename)

        ### ADD THE HEADER KEYWORDS THAT DRIZZLEPAC NEEDS HERE ###
        temp_sca[0].header["TELESCOP"] = "RST"
        temp_sca[0].header["INSTRUME"] = 'WFI   '
        temp_sca[0].header["DETECTOR"] = 'WFI'
        temp_sca[0].header["RA_TARG"] = point.ra.degree
        temp_sca[0].header["DEC_TARG"] = point.dec.degree
        temp_sca[0].header["PA_V3"] = PA
        temp_sca[0].header["EXPSTART"] = date.mjd
        temp_sca[0].header["EXTNAME"] = 'SCI     '
        temp_sca[0].header["EXTVER"] = SCA_i + 1
        temp_sca[0].header["EXPTIME"] = exptime
        temp_sca[0].header["BUNIT"] = "ELECTRONS"
        temp_sca[0].header["FILTER1"] = band
        temp_sca[0].header["FILTER2"] = "CLEAR"
        temp_sca[0].header["NGOODPIX"] = np.sum(np.isfinite(temp_sca[0].data))


        temp_header_list.append(temp_sca[0].header)

    utils.save_fits(array=temp_array_list, name=output, header=temp_header_list)

    # We might need to add info on the first extension. Load the multifits, modify and save again.

    final_fits = fits.open(output)
    final_fits[0].header["RA_TARG"] = final_fits[1].header["RA_TARG"]
    final_fits[0].header["DEC_TARG"] = final_fits[1].header["DEC_TARG"]
    final_fits[0].header["TELESCOP"] = final_fits[1].header["TELESCOP"]
    final_fits[0].header["INSTRUME"] = final_fits[1].header["INSTRUME"]
    final_fits[0].header["DETECTOR"] = final_fits[1].header["DETECTOR"]
    final_fits[0].header["EXPTIME"] = final_fits[1].header["EXPTIME"]
    final_fits[0].header["EXPSTART"] = final_fits[1].header["EXPSTART"]
    final_fits[0].header["FILTER1"] = final_fits[1].header["FILTER1"]
    final_fits[0].header["FILTER2"] = final_fits[1].header["FILTER2"]

    final_fits[0].header["PA_V3"] = final_fits[1].header["PA_V3"]
    #final_fits[0].header["TELESCOP"] = final_fits[1].header["TELESCOP"]

    final_fits.verify("silentfix")
    final_fits.writeto(output, overwrite=True)

    #    if i>0:
    #        roman_dummy[i].header["NGOODPIX"] = np.sum(np.isfinite(roman_dummy[i].data))
    #        print(i)
    #        print(roman_dummy[i].header["NGOODPIX"])


    return(output)

#################################################

def get_subarray_locations(SCA, verbose=False):
    """
    This program estimates the locations on the SCA detector of the different subarrays used for NDI estimation.
    Since some SCAS are flipped, and the number of subarrays used for might be variable, it is better to estimate
    it case by case.

    # Scott R. : SCA 3, 6, 9, 12, 15, 18 are flipped 180 deg.
    # 17.5, -17.5 mm correspond to the upper right corner points in each SCA.

    """

    ##############################################################################
    # Each SCA is divided in 8 x 8 subarrays.                                    #
    # Each subarray has its own NDI, so we need to run this process per subarray.#
    # First, we need to generate the limits of each subarray.                    #
    ##############################################################################
    if verbose: print(datetime.now().isoformat() + ": Setting up the grid of subarrays to estimate the NDI...")
    n_subarray_grid = 8
    data_shape = [4088, 4088]
    subarray_indexes = np.linspace(0, n_subarray_grid**2-1, n_subarray_grid**2).astype("int")
    subarray_step = int(data_shape[0]/n_subarray_grid)
    x_subarray_limits = np.linspace(0, data_shape[0], n_subarray_grid+1).astype("int")
    y_subarray_limits = np.linspace(0, data_shape[1], n_subarray_grid+1).astype("int")

    xy_subarray_limits = np.zeros(((n_subarray_grid+1)**2, 2))

    xmin = np.zeros((n_subarray_grid**2)).astype("int")
    xmax = np.zeros((n_subarray_grid**2)).astype("int")
    ymin = np.zeros((n_subarray_grid**2)).astype("int")
    ymax = np.zeros((n_subarray_grid**2)).astype("int")
    xmid = np.zeros((n_subarray_grid**2)).astype("int")
    ymid = np.zeros((n_subarray_grid**2)).astype("int")
    xlabel = np.zeros((n_subarray_grid**2))
    ylabel = np.zeros((n_subarray_grid**2))

    counter = 0
    # Scott R. : SCA 3, 6, 9, 12, 15, 18 are flipped 180 deg.
    # 17.5, -17.5 mm correspond to the upper right corner points in each SCA.

    #if (SCA == 3) or (SCA == 6) or (SCA == 9) or (SCA == 12) or (SCA == 15) or (SCA == 18):
    #    subarray_NDI_labels = np.array([17.5, 12.5, 7.5, 2.5, -2.5, -7.5, -12.5, -17.5])
    #else:
    #    subarray_NDI_labels = np.array([-17.5, -12.5, -7.5, -2.5, 2.5, 7.5, 12.5, 17.5])
    # I think they are already flipped in the NDI files
    #subarray_NDI_labels = np.array([-17.5, -12.5, -7.5, -2.5, 2.5, 7.5, 12.5, 17.5])

    subarray_NDI_labels = np.array([-17.92, -12.8, -7.68, -2.56, 2.56, 7.68, 12.8, 17.92])
    NDI_labels = []
    mask = np.zeros(data_shape)
    subarray_index = []
    for m in range(n_subarray_grid):
        for n in range(n_subarray_grid):
            xlabel[counter] = subarray_NDI_labels[m]
            ylabel[counter] = subarray_NDI_labels[n]

            xmin[counter] = x_subarray_limits[m]
            xmax[counter] = x_subarray_limits[m+1]
            ymin[counter] = y_subarray_limits[n]
            ymax[counter] = y_subarray_limits[n+1]
            xmid[counter] = int((x_subarray_limits[m] + x_subarray_limits[m+1])/2)
            ymid[counter] = int((y_subarray_limits[n] + y_subarray_limits[n+1])/2)
            mask[xmin[counter]:xmax[counter], ymin[counter]:ymax[counter]] = counter
            NDI_labels.append("X" + str(subarray_NDI_labels[m]) + "_Y" + str(subarray_NDI_labels[n]))
            subarray_index.append(counter)
            counter = counter + 1


    if verbose:
        fig, ax = plt.subplots(figsize=(16,16))
        im = ax.imshow(mask, origin="lower")
        for i in range(len(NDI_labels)):
            ax.text(xmid[i], ymid[i], NDI_labels[i])
        fig.colorbar(im, ax=ax)
    if verbose: print(datetime.now().isoformat() + ": Done")

    return({"subarray_index": np.array(subarray_index), "mask": mask, "xlabel": xlabel, "ylabel": ylabel, "xmin": xmin, "xmax": xmax, "xmid": xmid, "ymin": ymin, "ymax": ymax, "ymid": ymid, "NDI_labels": NDI_labels, "subarray_step": subarray_step})


#################################################

def _test_run_rosalia_asdf(catalog, exposure_name_list, outname="TEST_rosalia_", dry_mode=False, verbose=False, clean=True):
    from scipy.interpolate import RegularGridInterpolator

    for exposure_name in tqdm(exposure_name_list, position=0, leave=True):
        exposure_identity = rs.utils.exposure_inspector(exposure_name)

        irradiance_stars = const.c*(exposure_identity["FILTER_IDENTITY"]["filter_lambda_max"]-exposure_identity["FILTER_IDENTITY"]["filter_lambda_min"])/(exposure_identity["FILTER_IDENTITY"]["filter_lambda_ref"]**2)*((10**(-0.4*(catalog["mag"]+56.1)))*u.W/u.meter**2/u.Hz)

        test_stray = roman_estimate_straylight_SCA(data=exposure_identity["DATA"][0], wcs=exposure_identity["ASTROPYWCS"][0], SCA=exposure_identity["SCA"],
                                                   filter_identity=exposure_identity["FILTER_IDENTITY"],
                                                   ra_point=exposure_identity["RA_TARG"],
                                                   dec_point=exposure_identity["DEC_TARG"], pa_point=exposure_identity["PA"],
                                                   ra_stars=catalog["ra"], dec_stars=catalog["dec"], source_id = catalog["source_id"],
                                                   irradiance_stars=irradiance_stars,
                                                   verbose=verbose, dry_mode=dry_mode)

        rs.utils.save_fits(exposure_identity["DATA"][0], outname + str(exposure_identity["SCA"]).zfill(2) + ".fits", exposure_identity["ASTROPYWCS"][0].to_header())
        rs.utils.save_fits(test_stray, outname + str(exposure_identity["SCA"]).zfill(2) + "_stray.fits", exposure_identity["ASTROPYWCS"][0].to_header())

    #rs.utils.execute_cmd("swarp -d > swarp.conf") # Generate a default config file for swarp
    #rs.utils.execute_cmd("swarp -c swarp.conf -SUBTRACT_BACK N -BLANK_BADPIXELS Y -VERBOSE_TYPE QUIET " + outname +"*_stray.fits") # Run swarp on all the SCAs
    #rs.utils.execute_cmd("astwarp coadd.fits -h0 --scale=0.1,0.1") # Make a compressed version, for easiest visualization.
    #print("Result in coadd.fits and coadd_scaled.fits")

    final_mosaic = rs.utils.run_swarp(pattern=outname +"*_stray.fits", outname=outname +"_swarp_stray.fits")

    if clean:
        rs.utils.execute_cmd("rm " + outname + "*.fits") # Remove temporary files
    print("Done")
    return(final_mosaic)
#################################################


def roman_WFI_NDI_estimator_direct(ra_stars, dec_stars, ra_point, dec_point, pa_point, SCA, level, X_label=None, Y_label=None, verbose=True):
    import astropy.wcs as astropy_wcs
    from scipy.interpolate import RegularGridInterpolator
    from astropy.io import fits

    # Identify NDI calibration file

    if (level == 1) or (level == 2):
        ndi_name = os.path.dirname(rs.__file__) + "/CORE/NDI/RST/ndi_lvl" + str(level) + "/" + "lvl" + str(level) + "_SCA_" + str(SCA)+ "_SUB_X" + str(X_label) + "_Y" + str(Y_label) + "_TAN.fits"
        #if (level == 1): ndi_name = os.path.dirname(rs.__file__) + "/CORE/NDI/RST/ndi_lvl" + str(level) + "/" + "lvl" + str(level) + "_SCA_" + str(SCA)+ "_SUB_X" + str(X_label) + "_Y" + str(Y_label) + "_TAN.fits"
        #if (level == 2): ndi_name = os.path.dirname(rs.__file__) + "/CORE/NDI/RST/ndi_lvl" + str(level) + "/" + "lvl" + str(level) + "_SCA_" + str(SCA)+ "_TAN.fits"

        ndi_fits = fits.open(ndi_name)

        w = astropy_wcs.WCS(header=ndi_fits[0].header, fobj=ndi_fits, naxis=2)
        w.wcs.crota = -pa_point,-pa_point
        w.wcs.crval = ra_point,dec_point

        if verbose:
            ndi_fits[0].header = w.to_header()
            ndi_fits.verify("silentfix")
            ndi_fits.writeto(ndi_name.replace(".fits", "_temp_aligned.fits"), overwrite=True)
        #w = astropy_wcs.WCS(header=ndi_fits[0].header, fobj=ndi_fits, naxis=2)
        x_stars, y_stars = w.wcs_world2pix(ra_stars, dec_stars, 0)

        x_grid = np.linspace(0, ndi_fits[0].header["NAXIS1"]-1, ndi_fits[0].header["NAXIS1"])
        y_grid = np.linspace(0, ndi_fits[0].header["NAXIS2"]-1, ndi_fits[0].header["NAXIS2"])

        f = RegularGridInterpolator((x_grid, y_grid), np.flip(ndi_fits[0].data, axis=1).T, bounds_error=False, fill_value=0)

        xy_points = np.c_[x_stars, y_stars]
        return(f(xy_points))



    if (level==3):
        import healpy as hp

        ndi_name = os.path.dirname(rs.__file__) + "/CORE/NDI/RST/ndi_lvl" + str(level) + "/" + "lvl" + str(level) + "_SCA_" + str(SCA) + "_HP.fits"
        NDI_map_hp = hp.read_map(ndi_name)

        nside = hp.pixelfunc.get_nside(NDI_map_hp)
        longitude = ra_point * u.deg
        latitude = dec_point * u.deg
        position_angle = -pa_point*u.deg

        rot_custom = hp.Rotator(rot=[longitude.to_value(u.deg), latitude.to_value(u.deg), position_angle.to_value(u.deg)], inv=True)
        NDI_map_hp_rotated = rot_custom.rotate_map_pixel(NDI_map_hp) # The use of rotate_map_pixel from healpy introduces artifacts at low SNR.
        # This is relatively contradictory to what the documentation claims. https://healpy.readthedocs.io/en/latest/generated/healpy.rotator.Rotator.rotate_map_pixel.html#healpy.rotator.Rotator.rotate_map_pixel

        if verbose:
            hp.write_map(filename=ndi_name.replace(".fits", "_temp_aligned.fits"), m=NDI_map_hp_rotated, overwrite=True)

        pix = hp.ang2pix(nside=hp.pixelfunc.get_nside(NDI_map_hp), theta=ra_stars, phi=dec_stars, lonlat=True)

        return(NDI_map_hp_rotated[pix])

#################################################


#################################################

def roman_estimate_straylight_SCA(data, wcs, SCA, filter_identity, ra_stars, dec_stars, irradiance_stars, cat_id, source_id, ra_point, dec_point, pa_point, verbose=False, dry_mode=False):

    from tqdm import tqdm
    import numexpr
    import bottleneck as bn

    # Level 1 stars are the closest. R to the center of the SCA of 1 degree.
    # Level 2 stars have a distance between 1 degree and 10 degrees to the center of the SCA.
    # Level 3 stars have a distance of 10 degrees or more to the SCA.

    level_1_critical_distance_to_star = 1 # degree
    level_2_critical_distance_to_star = 10 # degree

    data_shape = data.shape
    pixsize = rs.telescopes.Roman.get_physical_pixelsize(instrument="WFI") # Physical pixel size, in meters

    # Find the coordinates of each pixel in the sky.
    if verbose: print(datetime.now().isoformat() + ": Finding coordinates of SCA pixels in the Sky...")
    SCA_pixel_radec = wcs.pixel_to_world(int(data_shape[0])/2,int(data_shape[1])/2)
    ra_SCA = SCA_pixel_radec.ra.value   # Right ascension of the center of the SCA
    dec_SCA = SCA_pixel_radec.dec.value # Declination of the center of the SCA
    if verbose: print(datetime.now().isoformat() + ": Done")

    radec_stars =  SkyCoord(ra_stars, dec_stars, frame="icrs", unit="deg")
    ra_stars     = radec_stars.ra.value.astype("float32")
    dec_stars    = radec_stars.dec.value.astype("float32")


    # Identify those at distances > 1 degrees (~100 times the size of the subarray).
    # For stars at those distances, we will only estimate one NDI and assume the straylight is flat across
    # the subarray (we will still have gradients!).
    if verbose: print(datetime.now().isoformat() + ": Estimating the relative coordinates of the SCA to the stars...")
    theta_phi_center_FPA = rs.utils.angular_distance(ra1=np.array([ra_point]), dec1=np.array([dec_point]),
                                                     ra2=ra_stars, dec2=dec_stars)
    if verbose: print(datetime.now().isoformat() + ": Done")

    if verbose:
        print("Stars - RA: " + str(radec_stars.ra) + " - DEC: " + str(radec_stars.dec))
        print("Irradiance: " + str(irradiance_stars))
        print("SCA - RA: " + str(ra_point) + " - DEC: " + str(dec_point))
    ##############################################################################

    # Here we split the stars into the different boundaries, depending on the distance to the center of the SCA.
    bool_is_the_star_level_1 = np.array(theta_phi_center_FPA[0] < level_1_critical_distance_to_star)
    bool_is_the_star_level_2 = np.array((theta_phi_center_FPA[0] >= level_1_critical_distance_to_star) & (theta_phi_center_FPA[0] < level_2_critical_distance_to_star))
    bool_is_the_star_level_3 = np.array(theta_phi_center_FPA[0] >= level_2_critical_distance_to_star)

    radec_level_1_stars = radec_stars[bool_is_the_star_level_1]
    radec_level_2_stars = radec_stars[bool_is_the_star_level_2]
    radec_level_3_stars = radec_stars[bool_is_the_star_level_3]

    id_stars = np.linspace(0, len(ra_stars)-1, len(ra_stars), dtype="int64") # np.array(cat_id)

    ra_level_1_stars   = radec_level_1_stars.ra.value
    dec_level_1_stars  = radec_level_1_stars.dec.value
    ra_level_2_stars   = radec_level_2_stars.ra.value
    dec_level_2_stars  = radec_level_2_stars.dec.value
    ra_level_3_stars   = radec_level_3_stars.ra.value
    dec_level_3_stars  = radec_level_3_stars.dec.value

    irradiance_level_1_stars = irradiance_stars[bool_is_the_star_level_1]

    source_id_level_1 = source_id[bool_is_the_star_level_1]
    id_level_1_stars = id_stars[bool_is_the_star_level_1]
    #print(id_level_1_stars)

    irradiance_level_2_stars = irradiance_stars[bool_is_the_star_level_2]
    source_id_level_2 = source_id[bool_is_the_star_level_2]
    id_level_2_stars = id_stars[bool_is_the_star_level_2]
    #print(id_level_2_stars)
    irradiance_level_3_stars = irradiance_stars[bool_is_the_star_level_3]
    source_id_level_3 = source_id[bool_is_the_star_level_3]
    id_level_3_stars = id_stars[bool_is_the_star_level_3]
    #print(id_level_3_stars)

    ##############################################################################


    n_level_1_stars = len(ra_level_1_stars)
    n_level_2_stars = len(ra_level_2_stars)
    n_level_3_stars = len(ra_level_3_stars)


    if verbose:
        print(datetime.now().isoformat() + ": Sources identified.")
        print("Number of stars at R < " + str(level_1_critical_distance_to_star) + " degrees - " + str(n_level_1_stars))
        print("Number of stars at " + str(level_1_critical_distance_to_star) + " < R < " + str(level_2_critical_distance_to_star) + " degrees - " + str(n_level_2_stars))
        print("Number of stars at R > " + str(level_2_critical_distance_to_star) + " degrees - " + str(n_level_3_stars))

    # This is the canvas array where we will store all the straylight level.
    straylight_SCA = np.zeros(data_shape).astype(np.float32)
    # This is the canvas array where we will store the ID of the largest stray-light contributor
    main_offender_SCA = np.zeros(data_shape).astype(np.float32)

    if True:
        subarray_locations_db = rs.roman.get_subarray_locations(SCA=SCA, verbose=False)
        xmid = subarray_locations_db["xmid"]
        ymid = subarray_locations_db["ymid"]

        # Estimating the coordinates of the center of the subarrays
        RA_mid_subarrays =       SCA_pixel_radec.ra.value#[ymid, xmid]
        DEC_mid_subarrays =      SCA_pixel_radec.dec.value#[ymid, xmid]
        mid_subarrays_Skycoord = SkyCoord(RA_mid_subarrays, DEC_mid_subarrays, frame="icrs", unit="deg")

        if verbose: print(datetime.now().isoformat() + ": Done")
        if verbose: print(datetime.now().isoformat() + ": Finding the NDI at the estimated angles...")


        ##############################################################################
        #  Level 3 stars only have one subarray.
        # So there is no need to make this pass through the subarray loop.
        #
        ##############################################################################


        # Level 3
        if len(ra_level_3_stars) > 0:
            NDI_level_3 = roman_WFI_NDI_estimator_direct(ra_stars=ra_level_3_stars, dec_stars=dec_level_3_stars,
                                                     ra_point=ra_point, dec_point=dec_point, pa_point=pa_point,
                                                     SCA=SCA, level=3, verbose=True)
            straylight_level_3 = (NDI_level_3*(pixsize**2)*filter_identity["filter_transmission_ref"]*filter_identity["filter_lambda_ref"]*irradiance_level_3_stars/const.c/const.h).decompose()
            where_max_stray = np.where(straylight_level_3.value == bn.nanmax(straylight_level_3))[0][0]
            id_main_offender_level_3 = id_level_3_stars[where_max_stray]
            source_id_main_offender_level_3 = source_id_level_3[where_max_stray]
            stray_main_offender_level_3 = bn.nanmax(straylight_level_3)
            max_NDI_level_3 = NDI_level_3[where_max_stray]

        else:
            straylight_level_3 = 0/u.s#
            id_main_offender_level_3 = 0
            stray_main_offender_level_3 = 0
            max_NDI_level_3 = 0
            source_id_main_offender_level_3 = None
        straylight_SCA = straylight_SCA + bn.nansum(straylight_level_3)


        ##############################################################################
        # Level 1 - The closest to the focal plane array.
        # Here we go subarray per subarray.
        # First - We measure the average Straylight at the center of the SCA for stars at large distances
        #.        Theta > 1 degrees
        # Second - We measure the straylight at the center of the subarray for stars at closer distances.
        #
        ################################################################################
        for i in range(len(subarray_locations_db["NDI_labels"])):
            X_label = subarray_locations_db["xlabel"][i]
            Y_label = subarray_locations_db["ylabel"][i]
            xmin    = subarray_locations_db["xmin"][i]
            xmax    = subarray_locations_db["xmax"][i]
            ymin    = subarray_locations_db["ymin"][i]
            ymax    = subarray_locations_db["ymax"][i]


            # -------------------------------------------------------------------------- #
            # 1 - Estimate the stray-light at the SCA level for close-distance stars     #
            # -------------------------------------------------------------------------- #
            if len(ra_level_1_stars) > 0:
                NDI_level_1 = roman_WFI_NDI_estimator_direct(ra_stars=ra_level_1_stars, dec_stars=dec_level_1_stars,
                                                           ra_point=ra_point, dec_point=dec_point, pa_point=pa_point,
                                                           SCA=SCA, X_label=X_label, Y_label=Y_label, level=1, verbose=True)
                straylight_level_1 = (NDI_level_1*(pixsize**2)*filter_identity["filter_transmission_ref"]*filter_identity["filter_lambda_ref"]*irradiance_level_1_stars/const.c/const.h).decompose()
                #id_main_offender_level_1 = id_level_1_stars[np.where(straylight_level_1.value == bn.nanmax(straylight_level_1))[0][0]]
                where_max_stray = np.where(straylight_level_1.value == bn.nanmax(straylight_level_1))[0][0]
                id_main_offender_level_1 = id_level_1_stars[where_max_stray]
                source_id_main_offender_level_1 = source_id_level_1[where_max_stray]
                max_NDI_level_1 = NDI_level_1[where_max_stray]

            else:
                stray_main_offender_level_1 = 0
                straylight_level_1 = 0/u.s
                id_main_offender_level_1 = 0
                max_NDI_level_1 = 0
                source_id_main_offender_level_1 = None
            # -------------------------------------------------------------------------- #
            # 3 - Combine all the values                                                 #
            # -------------------------------------------------------------------------- #
            stray_main_offender_level_1 = bn.nanmax(straylight_level_1)
            straylight_SCA[ymin:ymax, xmin:xmax] = straylight_SCA[ymin:ymax, xmin:xmax] + bn.nansum(straylight_level_1)

            # Level 2
            if len(ra_level_2_stars) > 0:
                NDI_level_2 = roman_WFI_NDI_estimator_direct(ra_stars=ra_level_2_stars, dec_stars=dec_level_2_stars,
                                                     ra_point=ra_point, dec_point=dec_point, pa_point=pa_point,
                                                     SCA=SCA, X_label=X_label, Y_label=Y_label, level=2, verbose=True)
                straylight_level_2 = (NDI_level_2*(pixsize**2)*filter_identity["filter_transmission_ref"]*filter_identity["filter_lambda_ref"]*irradiance_level_2_stars/const.c/const.h).decompose()  #
                where_max_stray = np.where(straylight_level_2.value == bn.nanmax(straylight_level_2))[0][0]
                id_main_offender_level_2 = id_level_2_stars[where_max_stray]
                source_id_main_offender_level_2 = source_id_level_2[where_max_stray]
                max_NDI_level_2 = NDI_level_2[where_max_stray]

            else:
                straylight_level_2 = 0/u.s
                id_main_offender_level_2 = 0
                stray_main_offender_level_2 = 0
                max_NDI_level_2 = 0
                source_id_main_offender_level_2 = None

            stray_main_offender_level_2 = bn.nanmax(straylight_level_2)
            straylight_SCA[ymin:ymax, xmin:xmax] = straylight_SCA[ymin:ymax, xmin:xmax] + bn.nansum(straylight_level_2)


            main_offender_level_123 = np.array([stray_main_offender_level_1,
                                                stray_main_offender_level_2,
                                                stray_main_offender_level_3])

            id_main_offender_level_123 = np.array([id_main_offender_level_1,
                                                   id_main_offender_level_2,
                                                   id_main_offender_level_3])

            source_name_main_offence_level_123 = np.array([source_id_main_offender_level_1, source_id_main_offender_level_2,source_id_main_offender_level_3])


            if verbose:
                print("Main offender - RA :" +  str(ra_stars[id_main_offender_level_123]) + " DEC: " + str(dec_stars[id_main_offender_level_123]))
                print("ID :" +  str(id_main_offender_level_123))
                print("Source name :" +  str(source_name_main_offence_level_123))
                print("Irradiance " + str(irradiance_stars[id_main_offender_level_123]))
                print("Stray " + str(main_offender_level_123))
                print("NDI " + str([max_NDI_level_1, max_NDI_level_2, max_NDI_level_3]))


            main_offender_SCA[ymin:ymax, xmin:xmax] = id_main_offender_level_123[main_offender_level_123 == bn.nanmax(main_offender_level_123)]
            #######################################################################


        # -------------------------------------------------------------------------- #
        # 4 - TO-DO - Interpolate the values across the SCA to smooth it out         #
        # -------------------------------------------------------------------------- #

        # from scipy.interpolate import RegularGridInterpolator
        #f = RegularGridInterpolator((x_grid, y_grid), np.flip(ndi_fits[0].data, axis=1).T, bounds_error=False, fill_value=0)

    if verbose: print(datetime.now().isoformat() + ": Done")
    return({"straylight_SCA": straylight_SCA, "main_offender_SCA": main_offender_SCA})

#################################################
