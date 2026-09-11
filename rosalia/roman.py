import os
import numpy as np
from tqdm import tqdm
import rosalia as rs
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord  # High-level coordinates
import astropy.units as u
from astropy import constants as const
from datetime import datetime
import logging
import healpy as hp
import pickle
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)

romanisim_bandpasses_names = ["R062", "Z087", "Y106", "J129", "H158", "W146", "F184", "K213"]
rosalia_bandpasses_names   = ["F062", "F087", "F106", "F129", "F158", "F146", "F184", "F213"]

# Build both mapping directions
bandpasses_names_to_rosalia = dict(zip(romanisim_bandpasses_names, rosalia_bandpasses_names))
bandpasses_names_to_romanisim = dict(zip(rosalia_bandpasses_names, romanisim_bandpasses_names))

def translate_bandpass_name(name):
    if name in bandpasses_names_to_rosalia:
        return bandpasses_names_to_rosalia[name]
    elif name in bandpasses_names_to_romanisim:
        return bandpasses_names_to_romanisim[name]
    else:
        return None


# --------------------- #
def create_roman_dummy(point, date, band, PA=None, exptime=500, output="default_roman_dummy.fits"):
    # create_roman_dummy:
    # Alejandro S. Borlaff - NASA/STA a.s.borlaff@nasa.gov
    #
    # Uses GalSim subroutines to generate a dummy empty Roman/WFI FITS file with the right WCS and array sizes.
    #
    # Update: Jan 15, 2026. Switching to use romanisim.get_wcs(...). Galsim WCS uses different conventions.
    # https://github.com/spacetelescope/romanisim/issues/302#issuecomment-3756578230
    #
    #
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
    #
    # --------------------------------------------- #

    #import galsim
    #import galsim.roman as roman

    # Setting up the parameters in GalSim format
    #ra_targ = coord.Angle(point.ra.degree*coord.degrees)
    #dec_targ = coord.Angle(point.dec.degree*coord.degrees)
    #targ_pos = galsim.CelestialCoord(ra=ra_targ, dec=dec_targ)

    #if PA is not None:
    #    WFI_PA = galsim.Angle(PA, coord.degrees)
    #else:
    #    WFI_PA = None
    from romanisim.ris_make_utils import set_metadata
    from romanisim.models import wcs as ris_wcs

    # import roman_datamodels

    # Starting the loop for the 18 SCAs
    temp_out_sca_filename_list = []
    temp_array_list = []
    temp_header_list = []
    temp_sca_name_list = []

    for SCA_i in range(18):
        SCA_name = str(SCA_i+1).zfill(3)

        temp_sca_name_list.append(SCA_name)
        out_sca_filename = output.replace(".fits", "_sca" + SCA_name + ".fits")
        temp_out_sca_filename_list.append(out_sca_filename)

        # Romanisim procedure>
        metadata = set_metadata(date=date,
                            bandpass=band,
                            sca=SCA_i+1,
                            ma_table_number=4,
                            usecrds=True,
                            truncate=None,
                            scale_factor=1)

        # Add this to include the distortion: 
        #distortion_file = parameters.reference_data["distortion"]
        #dist_model = roman_datamodels.datamodels.DistortionRefModel(distortion_file)
        #if distortion_file is not None:
        #    distortion = dist_model.coordinate_distortion_transform
        #else:
        #    distortion = None

        ris_wcs.fill_in_parameters(metadata, point, boresight=False, pa_aper=PA)
        sca_gwcs = ris_wcs.get_wcs(image=metadata, usecrds=True, distortion=None)
        sca_astropy_wcs = sca_gwcs.wcs.to_fits()[0]
        sca_dummy = np.zeros((4088, 4088))
        rs.utils.save_fits(array=sca_dummy, name=out_sca_filename, header=sca_astropy_wcs)

        # Open the temporary SCA image
        temp_sca = fits.open(out_sca_filename)
        temp_array_list.append(temp_sca[0].data)

        temp_sca[0].header["EXTNAME"] = "SCA" + SCA_name
        temp_sca[0].header["SCA"] = SCA_i+1


        os.remove(out_sca_filename)

        ### ADD THE HEADER KEYWORDS THAT DRIZZLEPAC NEEDS HERE ###
        temp_sca[0].header["TELESCOP"] = "RST"
        temp_sca[0].header["INSTRUME"] = 'WFI   '
        temp_sca[0].header["DETECTOR"] = 'WFI'
        temp_sca[0].header["RA_TARG"] = point.ra.degree
        temp_sca[0].header["DEC_TARG"] = point.dec.degree
        temp_sca[0].header["PA_V3"] = PA
        temp_sca[0].header["PA"] = PA
        temp_sca[0].header["EXPSTART"] = date.mjd
        temp_sca[0].header["EXTNAME"] = 'SCI     '
        temp_sca[0].header["EXTVER"] = SCA_i + 1
        temp_sca[0].header["EXPTIME"] = exptime
        temp_sca[0].header["BUNIT"] = "ELECTRONS"
        temp_sca[0].header["FILTER1"] = band
        temp_sca[0].header["FILTER2"] = "CLEAR"
        temp_sca[0].header["NGOODPIX"] = np.sum(np.isfinite(temp_sca[0].data))


        temp_header_list.append(temp_sca[0].header)

    rs.utils.save_fits(array=temp_array_list, name=output, header=temp_header_list)

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
    final_fits[0].header["PA"] = final_fits[1].header["PA_V3"]

    final_fits.verify("silentfix")
    final_fits.writeto(output, overwrite=True)

    return(output)

#################################################

def get_subarray_locations(SCA, verbose=False):
    """
    This program estimates the locations on the SCA detector of the different subarrays used for NDI estimation.
    Since some SCAS are flipped, and the number of subarrays used for might be variable, it is better to estimate
    it case by case.
    """
    # Scott R. : SCA 3, 6, 9, 12, 15, 18 are flipped 180 deg.
    # 17.5, -17.5 mm correspond to the upper right corner points in each SCA.

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

        # Reset the Roman / WFI loading bar:
        rs.plots.ascii_progress_focal_plane.canvas = rs.plots.ascii_progress_focal_plane.canvas_zero
        test_stray = roman_estimate_straylight_SCA(data=exposure_identity["DATA"][0],
                                                   wcs=exposure_identity["ASTROPYWCS"][0],
                                                   SCA=exposure_identity["SCA"],
                                                   filter_identity=exposure_identity["FILTER_IDENTITY"],
                                                   ra_point=exposure_identity["RA_TARG"],
                                                   dec_point=exposure_identity["DEC_TARG"],
                                                   pa_point=exposure_identity["PA"],
                                                   ra_stars=catalog["ra"],
                                                   dec_stars=catalog["dec"],
                                                   source_id = catalog["source_id"],
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


############################
def mag2fe(mag, bandpass, sca):
    from romanisim.bandpass import compute_count_rate

    back_single = 0
    if not isinstance(mag, (list, pd.core.series.Series, np.ndarray)):
        mag = np.array([mag])
        back_single = 1

    if isinstance(mag, (list)): mag = np.array(mag)


    # Zeropoint flux AB -> 0 = -2.5 log10(fzpAB) + zp
    # zp = 2.5*np.log10(fzpAB)
    flux = 10**(-0.4*(mag + 48.60))*u.erg/u.s/u.hertz/u.cm**2
    flux_e_s = np.zeros(len(mag))
    romanisim_bandpass_name = translate_bandpass_name(bandpass)
    for i in range(len(mag)):
        flux_e_s[i] = compute_count_rate(flux[i],romanisim_bandpass_name,sca,
                                         filename=None, effarea=None, wavedist=None)
    if back_single: flux_e_s = flux_e_s[0]
    return(flux_e_s/u.s)


def fe2mag(fe, bandpass, sca):
    back_single = 0
    if not isinstance(fe, (list, pd.core.series.Series, np.ndarray)):
        fe = np.array([fe])
        back_single = 1

    if isinstance(fe, (list)): fe = np.array(fe)

    # Zeropoint flux AB -> 0 = -2.5 log10(fzpAB) + zp
    # zp = 2.5*np.log10(fzpAB)

    zp = 2.5*np.log10(compute_abflux(sca, effarea=None)["SCA"+str(sca).zfill(2)][bandpass])

    mag = -2.5*np.log10(fe.value) + zp
    if back_single: mag = mag[0]
    return(mag)

#######################################


def roman_WFI_NDI_estimator_direct(ra_stars, dec_stars, ra_point, dec_point, pa_point, 
                                   level, ndi_name=None, ndi_wcs=None, ndi_grid_interpolator=None, rot_custom=None, verbose=True):

    # Identify NDI calibration file
    if (level == 1) or (level == 2):
        w = ndi_wcs
        w.wcs.crota = -pa_point,-pa_point
        w.wcs.crval = ra_point,dec_point
        x_stars, y_stars = w.wcs_world2pix(ra_stars, dec_stars, 0)
        xy_points = np.c_[x_stars, y_stars]
        return(ndi_grid_interpolator(xy_points))

    if (level==3):
        NDI_map_hp = hp.read_map(ndi_name)
        # rot_custom.rotate_map_pixel is relatively slow. We must find an alternative solution.
        # Perhaps there is a way to not rotate the map at all.
        # This is relatively contradictory to what the documentation claims. https://healpy.readthedocs.io/en/latest/generated/healpy.rotator.Rotator.rotate_map_pixel.html#healpy.rotator.Rotator.rotate_map_pixel
        NDI_map_hp_rotated = rot_custom.rotate_map_pixel(NDI_map_hp) # The use of rotate_map_pixel from healpy introduces artifacts at low SNR.
        pix = hp.ang2pix(nside=hp.pixelfunc.get_nside(NDI_map_hp), theta=ra_stars, phi=dec_stars, lonlat=True)
        return(NDI_map_hp_rotated[pix])

#################################################
#################################################

def roman_estimate_straylight_SCA(data_shape, wcs, SCA, filter_identity, ra_stars,
                                  dec_stars, irradiance_stars, cat_id,
                                  source_id, ra_point, dec_point, pa_point,
                                  verbose=False, loadingbar=False):

    import bottleneck as bn
    from IPython.display import clear_output

    # Level 1 stars are the closest. R to the center of the SCA of 1 degree.
    # Level 2 stars have a distance between 1 degree and 10 degrees to the center of the SCA.
    # Level 3 stars have a distance of 10 degrees or more to the SCA.

    level_1_critical_distance_to_star = 1 # degree
    level_2_critical_distance_to_star = 10 # degree

    pixsize = rs.telescopes.Roman.get_physical_pixelsize(instrument="WFI") # Physical pixel size, in meters

    # Find the coordinates of each pixel in the sky.
    if verbose > 1: print(datetime.now().isoformat() + ": Finding coordinates of SCA pixels in the Sky...")
    # SCA_pixel_radec = wcs.pixel_to_world(int(data_shape[0])/2,int(data_shape[1])/2)
    if verbose > 1: print(datetime.now().isoformat() + ": Done")

    radec_stars =  SkyCoord(ra_stars, dec_stars, frame="icrs", unit="deg")
    ra_stars     = radec_stars.ra.value.astype("float32")
    dec_stars    = radec_stars.dec.value.astype("float32")
    n_stars = len(ra_stars)

    # Identify those at distances > 1 degrees (~100 times the size of the subarray).
    # For stars at those distances, we will only estimate one NDI and assume the straylight is flat across
    # the subarray (we will still have gradients!).
    if verbose > 1: print(datetime.now().isoformat() + ": Estimating the relative coordinates of the SCA to the stars...")
    theta_phi_center_FPA = rs.utils.angular_distance(ra1=np.array([ra_point]), 
                                                     dec1=np.array([dec_point]),
                                                     ra2=ra_stars, 
                                                     dec2=dec_stars)
    
    if verbose > 1: print(datetime.now().isoformat() + ": Done")

    if verbose > 1:
        print("Stars - RA: " + str(radec_stars.ra) + " - DEC: " + str(radec_stars.dec))
        print("Irradiance: " + str(irradiance_stars))
        print("SCA - RA: " + str(ra_point) + " - DEC: " + str(dec_point))
    ##############################################################################

    # Here we split the stars into the different boundaries, depending on the distance to the center of the SCA.
    bool_is_the_star_level_1 = np.array(theta_phi_center_FPA[0] < level_1_critical_distance_to_star)
    bool_is_the_star_level_2 = np.array((theta_phi_center_FPA[0] >= level_1_critical_distance_to_star) & (theta_phi_center_FPA[0] < level_2_critical_distance_to_star))
    bool_is_the_star_level_3 = np.array(theta_phi_center_FPA[0] >= level_2_critical_distance_to_star)

    star_level = np.zeros(n_stars)
    star_level[bool_is_the_star_level_1] = 1
    star_level[bool_is_the_star_level_2] = 2
    star_level[bool_is_the_star_level_3] = 3

    radec_level_1_stars = radec_stars[bool_is_the_star_level_1]
    radec_level_2_stars = radec_stars[bool_is_the_star_level_2]
    radec_level_3_stars = radec_stars[bool_is_the_star_level_3]

    id_stars = cat_id # 

    ra_level_1_stars   = radec_level_1_stars.ra.value
    dec_level_1_stars  = radec_level_1_stars.dec.value
    ra_level_2_stars   = radec_level_2_stars.ra.value
    dec_level_2_stars  = radec_level_2_stars.dec.value
    ra_level_3_stars   = radec_level_3_stars.ra.value
    dec_level_3_stars  = radec_level_3_stars.dec.value

    # Level 1
    irradiance_level_1_stars = irradiance_stars[bool_is_the_star_level_1]
    source_id_level_1 = source_id[bool_is_the_star_level_1]
    id_level_1_stars = id_stars[bool_is_the_star_level_1]
    # Level 2
    irradiance_level_2_stars = irradiance_stars[bool_is_the_star_level_2]
    source_id_level_2 = source_id[bool_is_the_star_level_2]
    id_level_2_stars = id_stars[bool_is_the_star_level_2]
    # Level 3
    irradiance_level_3_stars = irradiance_stars[bool_is_the_star_level_3]
    source_id_level_3 = source_id[bool_is_the_star_level_3]
    id_level_3_stars = id_stars[bool_is_the_star_level_3]

    ##############################################################################


    n_level_1_stars = len(ra_level_1_stars)
    n_level_2_stars = len(ra_level_2_stars)
    n_level_3_stars = len(ra_level_3_stars)

    # 

    if verbose > 1:
        print(datetime.now().isoformat() + ": Sources identified.")
        print("Number of stars at R < " + str(level_1_critical_distance_to_star) + " degrees - " + str(n_level_1_stars))
        print("Number of stars at " + str(level_1_critical_distance_to_star) + " < R < " + str(level_2_critical_distance_to_star) + " degrees - " + str(n_level_2_stars))
        print("Number of stars at R > " + str(level_2_critical_distance_to_star) + " degrees - " + str(n_level_3_stars))

    # The transfer NDI maps units of 1/(superpixel size in mm2). 
    # We need to correct by that factor. 
    superpixel_mm2_area = (rs.telescopes.Roman.get_physical_pixelsize("WFI").to("mm").value * 512)**2

    # Get the locations of the subarrays in the selected SCA.
    subarray_locations_db = rs.roman.get_subarray_locations(SCA=SCA, verbose=False)

    # Initialize a minimal storage array for the straylight. 
    # The minimal storage array has 64 x 18, which is 64 locations, for each one of the 18 SCAs.
    n_subarrays = len(subarray_locations_db["NDI_labels"])

    # Initialize the minimal storage arrays: straylight, mainoffender, ra, dec
    straylight_total   = np.zeros((n_subarrays,)).astype(np.float32) 
    mainoffender_total = np.zeros((n_subarrays,)).astype(np.float32) 

    xmid = subarray_locations_db["xmid"]
    ymid = subarray_locations_db["ymid"]
    ramid, decmid = wcs.wcs_pix2world(xmid, ymid, 0) # subarray_locations_db["xmid"]
    

    # Estimating the coordinates of the center of the subarrays
    # RA_mid_subarrays =       SCA_pixel_radec.ra.value #[ymid, xmid]
    # DEC_mid_subarrays =      SCA_pixel_radec.dec.value #[ymid, xmid]
    #mid_subarrays_Skycoord = SkyCoord(RA_mid_subarrays, DEC_mid_subarrays, frame="icrs", unit="deg")

    if verbose > 1: print(datetime.now().isoformat() + ": Done")
    if verbose > 1: print(datetime.now().isoformat() + ": Finding the NDI at the estimated angles...")

    ##############################################################################
    # Initialize storage arrays for xmid, ymid, NDI, stray-light, main_offender
    ##############################################################################
    #### 
    # Prepare the healpix rotator for level 3 stars.
    longitude = ra_point * u.deg
    latitude = dec_point * u.deg
    position_angle = -pa_point*u.deg
    rot_custom = hp.Rotator(rot=[longitude.to_value(u.deg), latitude.to_value(u.deg), position_angle.to_value(u.deg)], inv=True)


    ############################################
    # Load the NDI maps for the level 1 and 2  #
    ############################################
    ndi_name_1 = os.environ['ROSALIACACHE'] + "/CORE/NDI/RST/ndi_lvl1/lvl1_SCA_" + str(SCA) + ".pkl"
    ndi_name_2 = os.environ['ROSALIACACHE'] + "/CORE/NDI/RST/ndi_lvl2/lvl2_SCA_" + str(SCA) + ".pkl"
    with open(ndi_name_1, "rb") as f:
        ndi_lvl1 = pickle.load(f)
    with open(ndi_name_2, "rb") as f:
        ndi_lvl2 = pickle.load(f)

    ###############################################
    # Precalculate constant factors for straylight estimation                                            
    ###############################################
    NDI_conversion_factor = ((pixsize**2)*filter_identity["filter_transmission_ref"]*filter_identity["filter_lambda_ref"]/superpixel_mm2_area/const.c/const.h)
    NDI_conversion_factor_1 = (NDI_conversion_factor*irradiance_level_1_stars).decompose()
    NDI_conversion_factor_2 = (NDI_conversion_factor*irradiance_level_2_stars).decompose()
    NDI_conversion_factor_3 = (NDI_conversion_factor*irradiance_level_3_stars).decompose()

    ##############################################################################
    # Here we go subarray per subarray.
    # First - We measure the average Straylight at the center of the SCA for stars at large distances
    #.        Theta > 1 degrees
    # Second - We measure the straylight at the center of the subarray for stars at closer distances.
    ################################################################################

    for i in range(n_subarrays): # Replaced by Roman / WFI loading bar
        X_label = subarray_locations_db["xlabel"][i]
        Y_label = subarray_locations_db["ylabel"][i]
        xmid    = subarray_locations_db["xmid"][i]
        ymid    = subarray_locations_db["ymid"][i]


        # -------------------------------------------------------------------------- #
        # 1 - Estimate the stray-light at the SCA level for close-distance stars     #
        # -------------------------------------------------------------------------- #
        if len(ra_level_1_stars) > 0:
            level=1
            transfer_level_1 = roman_WFI_NDI_estimator_direct(ra_stars=ra_level_1_stars, dec_stars=dec_level_1_stars,
                                                              ra_point=ra_point, dec_point=dec_point, pa_point=pa_point,
                                                                ndi_wcs=ndi_lvl1["wcs"][i], 
                                                                ndi_grid_interpolator=ndi_lvl1["ndi_interpolator"][i],
                                                                level=level, verbose=True)

            straylight_level_1 = transfer_level_1*NDI_conversion_factor_1  # (NDI_level_1*(pixsize**2)*filter_identity["filter_transmission_ref"]*filter_identity["filter_lambda_ref"]*irradiance_level_1_stars/const.c/const.h).decompose()
            
            # Once we have the straylight level, we can estimate: 
            # What is the star that produces the maximum straylight level on its own?
            # A.K.A. the main offender.
            where_max_stray = np.where(straylight_level_1.value == bn.nanmax(straylight_level_1))[0][0]
            id_main_offender_level_1        = id_level_1_stars[where_max_stray]
            source_id_main_offender_level_1 = source_id_level_1[where_max_stray]

        else:
            stray_main_offender_level_1 = 0
            straylight_level_1 = np.array([0/u.s])
            id_main_offender_level_1 = 0
            where_max_stray = 0
            source_id_main_offender_level_1 = None

        # -------------------------------------------------------------------------- #
        # 3 - Combine all the values                                                 #
        # -------------------------------------------------------------------------- #
        stray_main_offender_level_1 = straylight_level_1[where_max_stray]
        total_straylight_level_1 = bn.nansum(straylight_level_1)
        # straylight_SCA[ymin:ymax, xmin:xmax] = straylight_SCA[ymin:ymax, xmin:xmax] + total_straylight_level_1
        # straylight_minimal_storage_array[i] = straylight_minimal_storage_array[i] + total_straylight_level_1

        # -------------------------------------------------------------------------- #
        # 2 - Estimate the stray-light at the SCA level for mid-distance stars     #
        # -------------------------------------------------------------------------- #
        if len(ra_level_2_stars) > 0:
            level = 2
            transfer_level_2 = roman_WFI_NDI_estimator_direct(ra_stars=ra_level_2_stars, dec_stars=dec_level_2_stars,
                                                                ra_point=ra_point, dec_point=dec_point, pa_point=pa_point,
                                                                level=level,  ndi_wcs=ndi_lvl2["wcs"][i], 
                                                                ndi_grid_interpolator=ndi_lvl2["ndi_interpolator"][i], verbose=True)
            #NDI_level_2 = transfer_level_2/superpixel_mm2_area

            #straylight_level_2 = (NDI_level_2*(pixsize**2)*filter_identity["filter_transmission_ref"]*filter_identity["filter_lambda_ref"]*irradiance_level_2_stars/const.c/const.h).decompose()  #
            straylight_level_2 = transfer_level_2*NDI_conversion_factor_2

            where_max_stray = np.where(straylight_level_2.value == bn.nanmax(straylight_level_2))[0][0]
            id_main_offender_level_2 = id_level_2_stars[where_max_stray]
            source_id_main_offender_level_2 = source_id_level_2[where_max_stray]
            # max_NDI_level_2 = NDI_level_2[where_max_stray]

        else:
            straylight_level_2 = np.array([0/u.s])
            id_main_offender_level_2 = 0
            stray_main_offender_level_2 = 0
            where_max_stray = 0
            # max_NDI_level_2 = 0
            source_id_main_offender_level_2 = None

        stray_main_offender_level_2 = straylight_level_2[where_max_stray] # bn.nanmax(straylight_level_2)
        total_straylight_level_2 = bn.nansum(straylight_level_2)
        #straylight_SCA[ymin:ymax, xmin:xmax] = straylight_SCA[ymin:ymax, xmin:xmax] + total_straylight_level_2
        # straylight_minimal_storage_array[i] = straylight_minimal_storage_array[i] + total_straylight_level_2

        # -------------------------------------------------------------------------- #
        # 3 - Estimate the stray-light at the SCA level for long-distance stars      #
        # -------------------------------------------------------------------------- #

        if len(ra_level_3_stars) > 0:
            level = 3
            ndi_name_3 = os.environ['ROSALIACACHE'] + "/CORE/NDI/RST/ndi_lvl" + str(level) + "/" + "lvl" + str(level) + "_SCA_" + str(SCA) + "_SUB_X" + str(X_label) + "_Y" + str(Y_label) + "_HP.fits"

            transfer_level_3 = roman_WFI_NDI_estimator_direct(ra_stars=ra_level_3_stars, dec_stars=dec_level_3_stars,
                                                                ra_point=ra_point, dec_point=dec_point, pa_point=pa_point, 
                                                                ndi_name=ndi_name_3, level=level, rot_custom=rot_custom, verbose=True)
            
            #NDI_level_3 = transfer_level_3/superpixel_mm2_area
            #straylight_level_3 =     (NDI_level_3*(pixsize**2)*filter_identity["filter_transmission_ref"]*filter_identity["filter_lambda_ref"]*irradiance_level_3_stars/const.c/const.h).decompose()
            straylight_level_3 = transfer_level_3*NDI_conversion_factor_3
            where_max_stray = np.where(straylight_level_3.value == bn.nanmax(straylight_level_3))[0][0]
            id_main_offender_level_3 = id_level_3_stars[where_max_stray]
            source_id_main_offender_level_3 = source_id_level_3[where_max_stray]
            # stray_main_offender_level_3 = straylight_level_3[where_max_stray]
            # max_NDI_level_3 = NDI_level_3[where_max_stray]

        else:
            straylight_level_3 = np.array([0/u.s])#
            id_main_offender_level_3 = 0
            stray_main_offender_level_3 = 0
            where_max_stray = 0
            # max_NDI_level_3 = 0
            source_id_main_offender_level_3 = None



        stray_main_offender_level_3 = straylight_level_3[where_max_stray]
        total_straylight_level_3 = bn.nansum(straylight_level_3)
        #straylight_SCA[ymin:ymax, xmin:xmax] = straylight_SCA[ymin:ymax, xmin:xmax] + total_straylight_level_3
        # straylight_minimal_storage_array[i] = straylight_minimal_storage_array[i] + total_straylight_level_3


        # Combine the three levels of straylight and identify the main offender across all levels.
        #print(stray_main_offender_level_1)
        #print(stray_main_offender_level_2)
        #print(stray_main_offender_level_3)
        main_offender_level_123 = np.array([stray_main_offender_level_1.value,
                                            stray_main_offender_level_2.value,
                                            stray_main_offender_level_3.value])

        id_main_offender_level_123 = np.array([id_main_offender_level_1,
                                                id_main_offender_level_2,
                                                id_main_offender_level_3])

        source_name_main_offence_level_123 = np.array([source_id_main_offender_level_1,
                                                        source_id_main_offender_level_2,
                                                        source_id_main_offender_level_3])


        if verbose > 2:
            print("Main offender - RA :" +  str(ra_stars[id_main_offender_level_123]) + " DEC: " + str(dec_stars[id_main_offender_level_123]))
            print("ID :" +  str(id_main_offender_level_123))
            print("Source name :" +  str(source_name_main_offence_level_123))
            print("Irradiance " + str(irradiance_stars[id_main_offender_level_123]))
            print("Stray " + str(main_offender_level_123))
            # print("NDI " + str([max_NDI_level_1, max_NDI_level_2, max_NDI_level_3]))


        #straylight_SCA[ymin:ymax, xmin:xmax] = straylight_SCA[ymin:ymax, xmin:xmax] + total_straylight_level_1 + total_straylight_level_2 + total_straylight_level_3
        #main_offender_SCA[ymin:ymax, xmin:xmax] = id_main_offender_level_123[main_offender_level_123 == bn.nanmax(main_offender_level_123)][0]
        straylight_total[i] = total_straylight_level_1 + total_straylight_level_2 + total_straylight_level_3 # straylight_SCA[ymid, xmid]
        mainoffender_total[i] = id_main_offender_level_123[main_offender_level_123 == bn.nanmax(main_offender_level_123)][0] #main_offender_SCA[ymid, xmid]

        # Plot the Roman/WFI loading bar
        # os.system('clear')
        if loadingbar: clear_output(wait=True)
        if loadingbar: print(rs.plots.style.CYAN + "ROSALIA Stray-light Mapper: Scanning ... " + rs.plots.style.RESET)
        if loadingbar: rs.plots.print_ascii_focal_plane(x=X_label, y=Y_label, SCA=SCA)
        if loadingbar: print(rs.plots.style.CYAN + 'SCA ' + str(SCA) + " - Subarray X=" + str(X_label) + " - Subarray Y=" + str(Y_label))
        if loadingbar: print('SCA ' + str(SCA) + " out of 18: " + str(np.round(100*i/n_subarrays,2)) + '% completed' + rs.plots.style.RESET)
        #######################################################################


        #main_offender_db = pd.DataFrame({"xmid": xmid, "ymid": ymid,
        #                                 "ramid": ramid, "decmid": decmid,
        #                                 "xmin": xmin, "ymin": ymin,
        #                                 "xmax": xmax, "ymid": ymax,
        #                                 "RA_mid": RA_mid_subarrays, "DEC_mid": DEC_mid_subarrays,
        #                                 "stray": col_stray, "main_off_id": col_main_off_id})
        # -------------------------------------------------------------------------- #
        # 4 - TO-DO - Interpolate the values across the SCA to smooth it out         #
        # -------------------------------------------------------------------------- #

        # from scipy.interpolate import RegularGridInterpolator
        #f = RegularGridInterpolator((x_grid, y_grid), np.flip(ndi_fits[0].data, axis=1).T, bounds_error=False, fill_value=0)

    main_offender_db = pd.DataFrame({"SCA": len(straylight_total)*[SCA], 
                                     "xmid": subarray_locations_db["xmid"], "ymid": subarray_locations_db["ymid"],
                                     "ramid": ramid, "decmid": decmid,
                                     "xmin": subarray_locations_db["xmin"], "ymin": subarray_locations_db["ymin"],
                                     "xmax": subarray_locations_db["xmax"], "ymax": subarray_locations_db["ymax"],
                                     # "RA_mid": RA_mid_subarrays, "DEC_mid": DEC_mid_subarrays,
                                     "straylight_total": straylight_total, 
                                     "mainoffender_total": mainoffender_total})
    #######################################################################
    ##########################################
    # OUT OF THE FOR LOOP                    #
    ##########################################
    #######################################################################

    if verbose > 1: print(datetime.now().isoformat() + ": Done")
    #return({"straylight_SCA": straylight_SCA, 
    #        "main_offender_SCA": main_offender_SCA, 
    #        "main_offender_db": main_offender_db, 
    #        #"straylight_minimal_storage_array": straylight_minimal_storage_array,
    #        })
    return(main_offender_db)
#################################################


def make_average_ndi_maps():
    import glob
    ndi_lvl_1_list = glob.glob(os.environ["ROSALIACACHE"] + "/CORE/NDI/RST/ndi_lvl1/lvl1_SCA_*_SUB_X*_Y*_TAN.fits")
    ndi_lvl_2_list = glob.glob(os.environ["ROSALIACACHE"] + "/CORE/NDI/RST/ndi_lvl2/lvl2_SCA_*_SUB_X*_Y*_TAN.fits")
    ndi_lvl_3_list = glob.glob(os.environ["ROSALIACACHE"] + "/CORE/NDI/RST/ndi_lvl3/lvl3_SCA_*_SUB_X*_Y*_HP.fits")


    product_name_list = []
    for i, ndi_list in zip([1,2,3], [ndi_lvl_1_list, ndi_lvl_2_list, ndi_lvl_3_list]):

        # Check if the average maps are already there.
        name_average_ndi = os.path.dirname(os.environ["ROSALIACACHE"]) + "/CORE/NDI/RST/ndi_lvl" + str(i) + "_mean.fits"

        if os.path.exists(name_average_ndi):
            print("Average map found! " + name_average_ndi + " - Skipping")
            continue

        n_frames = len(ndi_list)
        if (i < 3):
            ndi_fits = fits.open(ndi_list[0])

        else:
            import healpy as hp
            ndi_fits = hp.read_map(ndi_list[0])

        canvas = np.zeros((ndi_fits[0].data.shape))
        for j in tqdm(range(n_frames)):
            if (i < 3):
                ndi_fits = fits.open(ndi_list[j])
                canvas  = canvas  + ndi_fits[0].data/n_frames
            else:
                ndi_fits = hp.read_map(ndi_list[0])
                canvas  = canvas  + ndi_fits/n_frames


        if (i < 3):
            rs.utils.save_fits(array = canvas, header=ndi_fits[0].header, name = name_average_ndi)

        else:
            hp.write_map(filename=name_average_ndi, m=canvas, overwrite=True)

        product_name_list.append(name_average_ndi)

    print("Average NDI maps completed! ")
    print(product_name_list)


###############

def make_romanisim_dummy(name, ra, dec, pa, bandpass, date, catalog=None):
    import glob
    previous_run = glob.glob(name + "*.asdf")

    if len(previous_run) == 0:
        os.system("romanisim-make-image --radec " + str(ra) + " " + str(dec) + " " + name + "{}.asdf --roll " + str(pa) + " --date " + date + "  --sca -1 --bandpass " + bandpass + " --level 2 --nobj 0 --usecrds --")
    else:
        print("ASDF exposures with same name already in place! Check local files")
        print(previous_run)

    current_run = glob.glob(name + "*.asdf")
    return(current_run)
