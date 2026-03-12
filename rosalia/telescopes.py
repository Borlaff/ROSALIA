import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
import bottleneck as bn
import astropy.units as u
from astropy import constants as const
from astropy.io import fits
from scipy import interpolate
from astropy.time import Time
import pandas as pd
import rosalia as rs
import skyfield as sf
from skyfield import api as sf_api
from astropy.convolution import Gaussian2DKernel

#####################
# TELESCOPE OBJECTS #
#####################

def ndi_estimator(theta, phi=None, wavelength=None, mode="legacy"):
    if mode == "legacy":
        #print("DEMO WARNING: Using legacy values of Bely 2003 for the HST NDI")
        ndi_db = pd.read_csv(os.environ["ROSALIACACHE"] + "/CORE/ndi_HST_legacy_bely2003.csv")
        ndi_theta_dummy = ndi_db["theta"]
        ndi_dummy = ndi_db["NDI"]
    if mode == "euclid":
        #print("DEMO WARNING: Using Euclid NDI instead of HST!")
        ndi_theta_dummy = np.linspace(0,180, 1000)
        ndi_dummy = rs.ndi.ndi_euclid.ndi_envelope(x=0, y=0, theta=ndi_theta_dummy, n_lambda=wavelength.to("nm").value)

    f_NDI_interpolator = interpolate.interp1d(x=ndi_theta_dummy, y=ndi_dummy, fill_value='extrapolate')

    return(f_NDI_interpolator(theta))

###########################
# TELESCOPE FINDER OBJECT #
###########################

def telescope_class_finder(telescope):

    # Supported telescopes
    supported_telescopes = {"Hubble": Hubble,
                            "Roman": Roman,
                            "CSST": CSST,
                            "SPHEREx": SPHEREx,
                            "ARRAKIHS": ARRAKIHS,
                            "Euclid": Euclid}

    if ("Hubble" in telescope) or ("HST" in telescope) or ("hst" in telescope):
        telescope_classname = "Hubble"

    elif ("Roman" in telescope) or ("ROMAN" in telescope) or ("RST" in telescope) or ("NGRST" in telescope):
        telescope_classname = "Roman"

    elif ("CSST" in telescope):
        telescope_classname = "CSST"

    elif ("SPHEREx" in telescope):
        telescope_classname = "SPHEREx"

    elif ("ARRAKIHS" in telescope):
        telescope_classname = "ARRAKIHS"

    elif ("Euclid" in telescope):
        telescope_classname = "Euclid"
    else:
        print("Telescope " + telescope + " not supported.")
        print("Try 'Hubble', 'Roman', 'SPHEREx', 'ARRAKIHS' or 'CSST' instead")
        telescope = None

    return(supported_telescopes[telescope_classname])
#################################################################

#############################
# TELESCOPE FILTER ROUTINES #
#############################

def find_filter_in_svo(wavelength, telescope, instrument, detector, verbose=False):
    from astroquery.svo_fps import SvoFps

    if verbose:
        print("Input filter name " + wavelength + "")
        print("Looking for " + wavelength + " in SVO library")

    # Handle exceptions with naming
    # ----------------------------------- #
    if telescope== "Hubble":
        telescope = "HST"
    if instrument== "ACS/WFC" or instrument== "acs" or instrument== "ACS/HRC" or instrument== "ACS/SBC":
        instrument = "ACS"
    if telescope== "RST":
        telescope = "Roman"

    #if True:
    try:
        filter_list = SvoFps.get_filter_list(facility=telescope, instrument=instrument).to_pandas()

    except:
        print("Filter " + wavelength +  " not found.")
        print("Please use a wavelength between 0.5 - 2000 microns")
        print("Or a STRAYCOR compatible filter name")
        print("No filter found. Check input data")
        print("Wavelength: " + wavelength)
        print("Telescope: " + telescope)
        print("Instrument: " + instrument)
        print("Detector: " + detector)

    #return(matched_filter)


    for i in range(len(filter_list)):
        if wavelength in filter_list["filterID"][i] and detector in filter_list["filterID"][i]:
            filter_id = filter_list["filterID"][i]
            filter_transmission_curve = SvoFps.get_transmission_data(filter_id)
            filter_lambda_ref= filter_list["WavelengthEff"][i]
            filter_lambda_width= filter_list["WidthEff"][i]

            good = np.where(np.isfinite(filter_transmission_curve["Wavelength"]) & np.isfinite(filter_transmission_curve["Transmission"]))
            wavelength_good = filter_transmission_curve["Wavelength"][good]
            transmission_good = filter_transmission_curve["Transmission"][good]
            transmission_interpolator = interpolate.interp1d(x=wavelength_good,
                                                             y=transmission_good,
                                                             kind="linear", fill_value="extrapolate")
            transmission_ref = transmission_interpolator(filter_lambda_ref)

            filter_lambda_min = wavelength_good[np.where(transmission_good > 0.05)[0][0]]
            filter_lambda_max = wavelength_good[np.where(transmission_good > 0.05)[0][-1]]

    if True:
        return({"filter_id": filter_id,
            "filter_transmission_curve": filter_transmission_curve,
            "filter_lambda_ref": filter_lambda_ref*u.AA,
            "filter_lambda_min": filter_lambda_min*u.AA,
            "filter_lambda_max": filter_lambda_max*u.AA,
            "filter_lambda_width": filter_lambda_width*u.AA,
            "filter_transmission_ref": transmission_ref,
            "wavelength": wavelength,
            "telescope": telescope,
            "instrument": instrument,
            "detector": detector})

    if False:
        print("No filter found. Check input data")
        print("Wavelength: " + wavelength)
        print("Telescope: " + telescope)
        print("Instrument: " + instrument)
        print("Detector: " + detector)



        #return({"instrument": instrument, "filter_name": filter_name,
        #        "wavelength_bins": wavelength_bins, "transmission_bins": transmission_bins, "average_wave_filter": average_wave_filter,
        #        "lambda_min": lambda_min, "lambda_max": lambda_max, "lambda_ref": lambda_ref, "transmission_ref": transmission_ref})
#################################################################

def get_filter(telescope, instrument, detector, filter_name, verbose=False):
    svo_filter_identity = rs.telescopes.find_filter_in_svo(wavelength=filter_name, telescope=telescope, instrument=instrument, detector=detector, verbose=False)
    if verbose: print(svo_filter_identity)
    wavelength_bins = np.array(svo_filter_identity["filter_transmission_curve"]["Wavelength"])*u.AA
    transmission_bins = svo_filter_identity["filter_transmission_curve"]["Transmission"]

    lambda_min = wavelength_bins[np.where(transmission_bins > 0.05)[0][0]]
    lambda_max = wavelength_bins[np.where(transmission_bins > 0.05)[0][-1]]
    lambda_ref = (lambda_min + lambda_max)/2.
    transmission_interpolator = interpolate.interp1d(x=wavelength_bins, y=transmission_bins)
    transmission_ref = transmission_interpolator(lambda_ref)

    # Now the median would be:
    average_wave_filter = np.sum(wavelength_bins*transmission_bins)/np.sum(transmission_bins)

    if verbose:
        import matplotlib.pyplot as plt
        plt.plot(wavelength_bins, transmission_bins, label=instrument +","+ filter_name)
        plt.axvline(average_wave_filter, color="black", linestyle = ":", linewidth=3, label="Median wavelength")
        plt.legend()
        plt.show()

    telescope_class = rs.telescopes.telescope_class_finder(telescope)
    telescope_area = telescope_class.mirror_radius**2 * np.pi

    if transmission_bins.unit == u.meter**2:
        print("SVO reported the transmission in m2. Correcting by area of the telescope (" + str(telescope_area) + " )")
        transmission_bins = np.array(transmission_bins)/telescope_area.value

    return({"instrument": instrument, "filter_name": filter_name,
            "wavelength_bins": wavelength_bins, "transmission_bins": transmission_bins, "average_wave_filter": average_wave_filter,
            "lambda_min": lambda_min, "lambda_max": lambda_max, "lambda_ref": lambda_ref, "transmission_ref": transmission_ref})
#################################################################

##########################
# FILTER BAND PROPERTIES #
##########################

class Passbands:
    class Gaia:
        class bp:
            lambda_ref = 5124.20*u.AA
            filter_width = 2333.06*u.AA

        class g:
            lambda_ref = 6251.50*u.AA
            filter_width = 4203.61*u.AA

        class rp:
            lambda_ref = 7829.65*u.AA
            filter_width = 2842.11*u.AA

    class TwoMASS:
        class J:
            lambda_ref = 12350.00*u.AA
            filter_width = 1624.32*u.AA

        class H:
            lambda_ref = 16620.00*u.AA
            filter_width = 2509.40*u.AA

        class Ks:
            lambda_ref = 21590.00*u.AA
            filter_width = 2618.87*u.AA

    class WISE:
        class W1:
            lambda_ref = 33526.00*u.AA
            filter_width = 6626.42*u.AA

        class W2:
            lambda_ref = 46028.00*u.AA
            filter_width = 10422.66*u.AA

        class W3:
            lambda_ref = 115608.00*u.AA
            filter_width = 55055.23*u.AA

        class W4:
            lambda_ref = 220883.00*u.AA
            filter_width = 41016.80*u.AA
#################################################################






#################################################################
################### TELESCOPE OBJECTS  ##########################
#################################################################

def find_closest_TLE(epoch, TLE_history):
    # This has to go under class telescopes. Same function but different list of TLEs
    ts = sf_api.load.timescale()
    TLE_mjd_list = np.zeros((len(TLE_history),))
    for i in range(len(TLE_history)):
        TLE_mjd_list[i] = TLE_history[i].epoch.to_astropy().mjd

    closest_TLE_index = rs.utils.find_nearest_index(array=TLE_mjd_list, value=epoch)
    closest_TLE = TLE_history[closest_TLE_index]

    return(closest_TLE)


#################################################################
#################################################################
#################################################################


class Hubble:
    """Hubble Space Telescope properties class"""
    TELESCOP = "HST"

    mirror_radius = 1.2*u.meter

    def f(self):
        return 'ROSALIA Hubble Space Telescope class'

    ndi_estimator = ndi_estimator


    def load_HST_TLEs(verbose=False):
        """ Loads the history of TLEs for Hubble Space Telescope.
            It return a list of Skyfield EarthSatellite objects TLE.
            Current version:   HST_TLE_history_sat000020580.txt - From launch (1990) - up to June 19, 2024.
        """
        ts = sf_api.load.timescale()

        with sf_api.load.open(os.environ["ROSALIACACHE"] + "/CORE/HST_TLE_history_sat000020580.txt") as f:
            satellites = list(sf.iokit.parse_tle_file(f, ts))
        if verbose: print('Loaded', len(satellites), ' TLEs of Hubble')
        return(satellites)


    def TLE_exposure(epoch):
        return(find_closest_TLE(epoch=epoch, TLE_history=Hubble.load_HST_TLEs()))


    def get_physical_pixelsize(instrument):
        if ("acs" in instrument) or ("ACS" in instrument):
            pixelsize = 15E-6*u.meter

        if instrument == "WFC3/UVIS":
            pixelsize = 15E-6*u.meter

        if instrument == "WFC3/IR":
            pixelsize = 18E-6*u.meter

        if instrument == "WFPC/PC":
            pixelsize = 15E-6*u.meter

        if instrument == "WFPC/WF":
            pixelsize = 15E-6*u.meter

        return(pixelsize)


    def get_canvas_shape(instrument):
        if instrument == "ACS":
            canvas_shape = [4144, 4144]

        return(canvas_shape)



    def get_pixscale(instrument):

        if instrument == "ACS":
            pixscale = 0.05*u.arcsec

        if instrument == "WFPC/WF":
            pixscale = 0.100*u.arcsec

        if instrument == "WFPC/PC":
            pixscale = 0.043*u.arcsec

        if instrument == "WFPC2/WF":
            pixscale = 0.100*u.arcsec

        if instrument == "WFPC2/PC":
            pixscale = 0.046*u.arcsec


        return(pixscale)


    def get_filter(instrument, filter_name, verbose=False):
        return(rs.telescopes.get_filter(telescope="HST", instrument=instrument, detector=instrument, filter_name=filter_name, verbose=False))

    def make_dummy_exposure(ra, dec, pa, outname):
        exp_name = "/Users/aborlaff/NASA/SPARKLES/notebooks/SATELLITES/HST_ACS_mock_es.fits"

        exp = fits.open(exp_name)
        canvas = exp[0].data
        canvas_shape = canvas.shape
        cdelt = [-0.05/60/60, 0.05/60/60]

        header = rs.utils.create_custom_wcs(crpix=[int(canvas_shape[1]/2), int(canvas_shape[0]/2)],
                                            crval=[ra,dec],
                                            cdelt=cdelt,
                                            crota=[-pa,-pa])
        rs.utils.save_fits(canvas, outname, header)
        return(outname)

    def mu2fe(mu):
        return(rs.detectors.mu2fe(mu=mu, instrument="ACS", filter_name="F814W", telescope="Hubble", verbose=False))

    def fe2mu(fe):
        return(rs.detectors.fe2mu(fe=fe, instrument="ACS", filter_name="F814W", telescope="Hubble",verbose=False))

    def get_HST_ACS_psf(wavelength, chip="WFC1", xy=[2000, 1000], output=None, verbose=1):
        # from focus_diverse_epsfs import interp_epsf, psf_retriever
        from acstools.focus_diverse_epsfs import psf_retriever, interp_epsf
        # provide an existing file location for download
        download_location = "."
        # call the psf_retriever function with observation rootname

        print("DEMO WARNING: HST PSF currently only uses ACS/WFC F775W. This will be improved in future ROSALIA versions. Use with caution")

        retrieved_download = psf_retriever("jezz03zbq", download_location)
        # open the file with astropy . io and grab the image data
        retrieved_filepath = (retrieved_download)
        psf_hdu = fits.open(retrieved_filepath)
        ePSFs = fits.getdata(retrieved_filepath, ext=0)
        interpolated_epsf = interp_epsf(ePSFs, x=xy[0], y=xy[1], chip="WFC2", pixel_space=True)
        psf_hdu[0].data = interpolated_epsf

        if output is None:
            output = "acs_" + chip + "_" + wavelength + "_x" + str(xy[0]) + "_y" + str(xy[1]) + "_ovsamp1.fits"

        psf_hdu.verify("silentfix")
        psf_hdu.writeto(output, overwrite=True, output_verify='ignore')
        if verbose: print("HST ACS PSF model saved in: " + output)

        return([output, psf_hdu])

    def get_PSF():
        #output, psf_hdu = rs.telescopes.Hubble.get_HST_ACS_psf(wavelength="F814W", chip="WFC1", xy=[2000, 1000], output=None, verbose=1)
        psf_hdu = fits.open("/Users/aborlaff/NASA/SPARKLES/notebooks/SATELLITES/psf_tinytim00_psf_scaled_cropped.fits")
        return(psf_hdu[1].data)

    def get_location(mjd, pov="@sun"):
        query = rs.horizons.horizons_query(mjd=mjd, source="@hst", location=pov)
        vectors = query["jpl_horizons_query"].vectors()
        return(np.array([vectors["x"].value, vectors["y"].value, vectors["z"].value])*u.AU)

#################################################################

class Roman:
    """Nancy Grace Space Telescope properties class"""
    TELESCOP = "Roman"
    WFI_SCAs = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18])

    mirror_radius = 1.2*u.meter


    def f(self):
        return 'ROSALIA Roman Space Telescope class'

    ndi_estimator = ndi_estimator

    def get_physical_pixelsize(instrument):
        if instrument == "WFI":
            pixelsize = 10E-6*u.meter

        return(pixelsize)

    def get_pixscale(instrument):
        if instrument == "WFI":
            pixscale = 0.11*u.arcsec

        return(pixscale)


    def get_filter(instrument, filter_name, verbose=False):
        return(rs.telescopes.get_filter(telescope="Roman", instrument=instrument, detector=instrument, filter_name=filter_name, verbose=False))

    def get_location(mjd, pov="@sun"):
        query = rs.horizons.horizons_query(mjd=mjd, source="@jwst", location=pov)
        vectors = query["jpl_horizons_query"].vectors()
        return(np.array([vectors["x"].value, vectors["y"].value, vectors["z"].value])*u.AU)

    def stpsf_get_psf(detector_position, detector, filter_name, fov_pixels=1024, oversample=1):
        # This is just a wrapper for the stpsf wfi module to generate Roman's PSF.
        # import stpsf
        wfi = stpsf.roman.WFI()
        wfi.filter = filter_name
        wfi.detector = 'WFI' + str(detector).zfill(2)
        wfi.detector_position = detector_position
        psf_wfi = wfi.calc_psf(fov_pixels = fov_pixels, oversample=oversample)
        return(psf_wfi[0])

    def get_bestPA(ra, dec, mjd):
        """
        get_bestPA calculates the best position angle for Roman's WFI given the target coordinates and the time of observation. It uses the galsim.roman.bestPA function to compute the optimal position angle that gets the Sun on the sunshield. The function takes in the right ascension (ra), declination (dec), and modified Julian date (mjd) as inputs and returns the best position angle in degrees.
       
        Input parameters:
        - ra: Right ascension of the target in degrees.
        - dec: Declination of the target in degrees.
        - mjd: Modified Julian Date of the observation.

        Output:
        - WFI_PA: The best position angle for Roman's WFI in degrees.

        """

        import galsim
        import coord
        import galsim.roman as galsim_roman
        from astropy.time import Time

        # Setting up the parameters in GalSim format
        ra_targ = coord.Angle(ra*coord.degrees)
        dec_targ = coord.Angle(dec*coord.degrees)
        targ_pos = galsim.CelestialCoord(ra=ra_targ, dec=dec_targ)
        t = Time(mjd, format='mjd')
        dt = t.to_datetime()
    
        WFI_PA = galsim_roman.bestPA(targ_pos, dt)
        return(WFI_PA.deg)



    def get_psf(detector_position, detector, filter_name):
        SCA = detector
        from romanisim.bandpass import roman2galsim_bandpass
        import galsim
        import galsim.roman as galsim_roman

        # Get the filter identity
        FILTER_IDENTITY = rs.telescopes.find_filter_in_svo(wavelength=filter_name,
                                                           telescope="Roman",
                                                           instrument="WFI",
                                                           detector="WFI",
                                                           verbose=False)

        bandpass = roman2galsim_bandpass[filter_name]
        gsparams = galsim.GSParams(maximum_fft_size=2**16, maxk_threshold=1.e-10, folding_threshold=1.e-6)

        class SCA_pos:
            x = detector_position[0]
            y = detector_position[1]

        psf = galsim_roman.getPSF(SCA, bandpass, SCA_pos=SCA_pos, pupil_bin=1, wcs=None,
                                  n_waves=10, extra_aberrations=None,
                                  wavelength=FILTER_IDENTITY["filter_lambda_ref"].to("nm").value,
                                  gsparams=gsparams, logger=None,
                                  high_accuracy=None, approximate_struts=None)


        img = psf.drawImage()
        outname = 'temp_psf_' + filter_name + '_SCA' + str(SCA).zfill(2) + '.fits'
        img.write(outname)
        psf_wfi = fits.open(outname)
        os.system("rm " + outname)
        return(psf_wfi[0])

    def find_wfi_center_for_offset_target(ra_target, dec_target, mjd, dX, dY, PA_wfi=None, verbose=False):
        # dX dY in degrees
    
        import pysiaf
        # Define the Roman Space Telescope Frame 
        rsiaf = pysiaf.Siaf('Roman')
        wfi_cen = rsiaf['WFI_CEN']

        # Find the optimal position angle of the observatory at that RA, Dec, and date
        # if PA_wfi is None:
        # This is the position angle of the observatory
        # WFI focal plane Y direction is -60 degrees offset from this. 
        PA_v3 = rs.telescopes.Roman.get_bestPA(ra=ra_target, dec=dec_target, mjd=mjd)
            #PA_wfi = PA_v3 - 60*u.degree
        # else:
            #PA_v3 = PA_wfi + 60*u.degree
        
        if verbose:
            print("RA: " + str(ra_target) + " - DEC: " + str(dec_target) + " PA_v3: " + str(PA_v3))
    
        # First, we find v2, v3 for the stray light point we want to hit
        # these are the coordinates in the telescope (V2, V3) plane, 
        # with origin in WFI_CEN, the center of the focal plane array.
        v2, v3 = wfi_cen.idl_to_tel(dX*60*60, dY*60*60, method="spherical", 
                                    input_coordinates="polar", output_coordinates="polar")


        # Then we define the attitude matrix required to place that v2 and v3 on the star
        attmat = pysiaf.utils.rotations.attitude_matrix(nu2=v2, nu3=v3,
                                                    ra=ra_target, dec=dec_target, 
                                                    pa=PA_v3)

        # We apply the attitude matrix to the observatory
        wfi_cen.set_attitude_matrix(attmat)
        # compute the position angle at that aperture
        V2Ref = rsiaf.apertures['WFI_CEN'].V2Ref
        V3Ref = rsiaf.apertures['WFI_CEN'].V3Ref
        pa = pysiaf.utils.rotations.posangle(attmat, V2Ref, V3Ref)
    
        # Compute the sky coordinates of the WFI_CEN aperture reference position
        wfi_ra, wfi_dec = wfi_cen.idl_to_sky(0, 0)
        if verbose: print(f'' + name.iloc[i] + f' - WFI_CEN: RA = {wfi_ra:.5f} deg, Dec = {wfi_dec:.5f} deg')
        return({"ra_wficen": wfi_ra, "dec_wficen": wfi_dec, "pa_wfi": pa - 60})
    
#################################################################

class CSST:
    """Chinese Space Station Telescope properties class"""
    TELESCOP = "CSST"

    #########
    # https://iopscience.iop.org/article/10.3847/1538-4365/acfc3a/pdf
    # CSST wide-field
    # multiband imaging survey is set to image approximately 17,500
    # deg2 of the sky using the near-ultraviolet (NUV), u, g, r, i, z,
    ##########

    # CSST Properties: https://arxiv.org/pdf/2304.02196
    filter_names = ["MCI.F275W", "MCI.F336W", "MCI.F375M",
                    "MCI.F450M", "MCI.F500M", "MCI.F630M",
                    "MCI.F763M", "MCI.F845M", "MCI.F960M"]
    mirror_radius = 1.0*u.meter

    def f(self):
        return 'STRAYCOR Chinese Space Station Telescope class'

    ndi_estimator = ndi_estimator

    def load_example_TLE_CSST():
        from skyfield.api import EarthSatellite
        ts = sf_api.load.timescale()
        # CSS (TIANHE) 9 Feb 2026
        # 1 48274U 21035A   26040.41958422  .00018345  00000+0  22932-3 0  9994
        # 2 48274  41.4684  20.0103 0004697  39.0740 321.0437 15.59384169273192
        #
        #
        line1 = '1 48274U 21035A   26040.41958422  .00018345  00000+0  22932-3 0  9994'
        line2 = '2 48274  41.4684  20.0103 0004697  39.0740 321.0437 15.59384169273192'
        satellite = EarthSatellite(line1, line2, 'CSST', ts)
        return(satellite)

    def TLE_exposure(epoch):
        return(rs.telescopes.CSST.load_example_TLE_CSST())

    def get_physical_pixelsize(instrument):
        # https://www.teledyneimaging.com/en/aerospace-and-defense/products/sensors-overview/ccd/ccd290-99/
        if ("CSST" in instrument) or ("csst" in instrument) or ("MCI" in instrument) or ("mci" in instrument):
            focal_length = 28*u.meter
            pixelsize = 9.9E-6*u.meter # http://www.waloszek.de/astro_foto_sensor_e.php
        return(pixelsize)

    def get_pixscale(instrument):

        if ("CSST" in instrument) or ("csst" in instrument) or ("MCI" in instrument) or ("mci" in instrument):
            pixscale = 0.073*u.arcsec

        return(pixscale)

    def get_canvas_shape(instrument):
        # Detector type for CSST
        # 30 CCD290-99
        # https://indico.cern.ch/event/1026001/contributions/4385928/attachments/2264178/3873429/Detector_Effects_Simulation_For_CSST_MXM_20210712.pdf
        if instrument == "CSST":
            canvas_shape = np.array(np.array([1.1, 1.1])*60*60/0.073*u.arcsec).astype("int")
        return(canvas_shape)

    def old_make_dummy_exposure(ra, dec, pa, outname, binning=8):
        xCCD = 9216/binning
        yCCD = 9232/binning
        gap_x = 100/binning
        gap_y = 400/binning
        canvas_shape = np.array([int(yCCD*5+gap_y*4), int(xCCD*6+gap_x*5)]) # rs.telescopes.CSST.get_canvas_shape(instrument="CSST")

        canvas = np.zeros(canvas_shape)

        # gap_x
        canvas[:, int(xCCD):int(xCCD+gap_x)] = np.nan
        canvas[:, int(2*xCCD+gap_x):int(2*xCCD+2*gap_x)] = np.nan
        canvas[:, int(3*xCCD+2*gap_x):int(3*xCCD+3*gap_x)] = np.nan
        canvas[:, int(4*xCCD+3*gap_x):int(4*xCCD+4*gap_x)] = np.nan
        canvas[:, int(5*xCCD+4*gap_x):int(5*xCCD+5*gap_x)] = np.nan

        # gap_y
        canvas[int(yCCD):int(yCCD+gap_y), :] = np.nan
        canvas[int(2*yCCD+gap_y):int(2*yCCD+2*gap_y), :] = np.nan
        canvas[int(3*yCCD+2*gap_y):int(3*yCCD+3*gap_y), :] = np.nan
        canvas[int(4*yCCD+3*gap_y):int(4*yCCD+4*gap_y), :] = np.nan
        canvas[int(5*yCCD+4*gap_y):int(5*yCCD+5*gap_y), :] = np.nan


        #canvas[:,(2*canvas_shape[0]+gap_width):2*(canvas_shape[0]+gap_width)] = np.nan
        header = rs.utils.create_custom_wcs(crpix=[int(canvas_shape[0]/2), int(canvas_shape[1]/2)],
                                            crval=[ra,dec],
                                            cdelt=binning*np.array([-rs.telescopes.CSST.get_pixscale(instrument="CSST").to("degree").value, rs.telescopes.CSST.get_pixscale(instrument="CSST").to("degree").value]),
                                            crota=[-pa,-pa])

        header["TELESCOP"] =  "CSST"
        header["INSTRUME"] =  "MCI"
        header["DETECTOR"] =  "MCI"
        header["FILTER"] = "MCI.F630M"
        header["FILTER1"] = "MCI.F630M"
        header["FILTER2"] = "None"
        header["PA"] = pa
        header["EXTNAME"] = "SCI"
        header["RA_TARG"] = ra
        header["DEC_TARG"] = dec


        rs.utils.save_fits(canvas, outname, header)
        return(outname)

    def get_PSF():
    #if True:
        ### CSST PSF ###
        # https://academic.oup.com/mnras/article/528/2/2770/7529202
        # The PSF of CSST grating is assumed to be a 2D Gaussian function
        # with a radius of REE80 ≲ 0.3 enclosing 80 per cent energy on average.
        sigma_CSST_psf = 1/1.282 * 0.3/0.0730000000000008 / 3 # We divide by 3 because the CSST images are undersampled by a factor 3, for computational purposes.
        gaussian_2D_kernel = Gaussian2DKernel(sigma_CSST_psf, x_size=31, y_size=31)
        return(gaussian_2D_kernel.array)


    def make_dummy_exposure(ra, dec, pa, outname):
        exp_name = "/Users/aborlaff/NASA/SPARKLES/notebooks/SATELLITES/CSST_mock_es.fits"

        exp = fits.open(exp_name)
        canvas = exp[0].data
        canvas_shape = canvas.shape

        cdelt = [-0.0730000000000008/60/60*3, 0.0730000000000008/60/60*3]
        header = rs.utils.create_custom_wcs(crpix=[int(canvas_shape[1]/2), int(canvas_shape[0]/2)],
                                            crval=[ra,dec],
                                            #cdelt=[-rs.telescopes.CSST.get_pixscale(instrument="CSST").to("degree").value,
                                            #        rs.telescopes.CSST.get_pixscale(instrument="CSST").to("degree").value],
                                            cdelt=cdelt,
                                            crota=[-pa,-pa])
        rs.utils.save_fits(canvas, outname, header)
        return(outname)


    def get_filter(instrument, filter_name, verbose=False):
        return(rs.telescopes.get_filter(telescope="CSST", instrument=instrument, detector="MCI", filter_name=filter_name, verbose=False))


    def mu2fe(mu, instrument="MCI", filter_name="MCI.F630M", telescope="CSST", verbose=False):
        return(rs.detectors.mu2fe(mu=mu,
                                  instrument=instrument,
                                  filter_name=filter_name,
                                  telescope=telescope,
                                  verbose=False))

    def fe2mu(fe, instrument="MCI", filter_name="MCI.F630M", telescope="CSST", verbose=False):
        return(rs.detectors.fe2mu(fe=fe,
                                  instrument=instrument,
                                  filter_name=filter_name,
                                  telescope=telescope,
                                  verbose=False))


#################################################################

class SPHEREx:
    """ SPHEREx Telescope properties class"""
    TELESCOP = "SPHEREx"
    mirror_radius = 0.1*u.meter

    def TLE_exposure(epoch):
        #return(rs.satellites.find_closest_TLE(epoch=epoch, TLE_history=CSST.load_CSST_TLE()))
        return(rs.telescopes.SPHEREx.load_example_TLE())

    def load_example_TLE():
        from skyfield.api import EarthSatellite
        ts = sf_api.load.timescale()
        #line1 = '1 44874U 19092B   24260.76087593  .00001534  00000+0  32292-3 0  9990'
        #line2 = '2 44874  98.1547  87.1562 0009664 302.2867  57.7397 14.60891889252646'
        #line1 = '1 99999U          24153.00000000  .00000000  00000-0  00000-0 0 00003'
        #line2 = '2 99999 097.9960 339.9709 0011528 246.2443 113.7557 14.75451568000015'
        line1 = '1 63182U 25047E   26040.49042002  .00001639  00000+0  25830-3 0  9990'
        line2 = '2 63182  97.9549 227.6915 0009366 239.5726 120.4558 14.73753633 49230'
        satellite = EarthSatellite(line1, line2, 'SPHEREx', ts)
        return(satellite)

    def get_canvas_shape(instrument):
        if instrument == "SPHEREx":
            canvas_shape = np.array(np.array([3.5, 11.3])*60*60/6.2*u.arcsec).astype("int")
        return(canvas_shape)


    def get_pixscale(instrument):
        if ("SPHEREx" in instrument) or ("spherex" in instrument):
            pixscale = 6.2*u.arcsec
        return(pixscale)

    def old_make_dummy_exposure(ra, dec, pa, outname):
        canvas_shape = rs.telescopes.SPHEREx.get_canvas_shape(instrument="SPHEREx")
        gap_width = int((canvas_shape[1] - 3*canvas_shape[0])/2)
        canvas = np.zeros(canvas_shape)
        canvas[:,canvas_shape[0]:canvas_shape[0]+gap_width] = np.nan
        canvas[:,(2*canvas_shape[0]+gap_width):2*(canvas_shape[0]+gap_width)] = np.nan
        header = rs.utils.create_custom_wcs(crpix=[int(canvas_shape[1]/2), int(canvas_shape[0]/2)],
                                            crval=[ra,dec],
                                            cdelt=[-rs.telescopes.SPHEREx.get_pixscale(instrument="SPHEREx").to("degree").value, rs.telescopes.SPHEREx.get_pixscale(instrument="SPHEREx").to("degree").value],
                                            crota=[-pa,-pa])
        rs.utils.save_fits(canvas, outname, header)
        return(outname)

    def make_dummy_exposure(ra, dec, pa, outname):
        exp1_name = "/Users/aborlaff/NASA/SPARKLES/notebooks/SATELLITES/SPHEREx_simdata_AAS2025/level2_2025W17_2A_0089_4D1_simu_HighLat.fits"
        exp2_name = "/Users/aborlaff/NASA/SPARKLES/notebooks/SATELLITES/SPHEREx_simdata_AAS2025/level2_2025W17_2A_0089_4D2_simu_HighLat.fits"
        exp3_name = "/Users/aborlaff/NASA/SPARKLES/notebooks/SATELLITES/SPHEREx_simdata_AAS2025/level2_2025W17_2A_0089_4D3_simu_HighLat.fits"

        exp1 = fits.open(exp1_name)
        exp2 = fits.open(exp2_name)
        exp3 = fits.open(exp3_name)

        # Convert the MJy/sr to e/s
        exp1[1].data = exp1[1].data * exp1[1].header["CAL"]
        exp2[1].data = exp2[1].data * exp2[1].header["CAL"]
        exp3[1].data = exp3[1].data * exp3[1].header["CAL"]

        spherex_image_size = 2040
        gap_width = 230 # pixels
        canvas = np.zeros((spherex_image_size, spherex_image_size*3 + gap_width*2))
        canvas[:, 0:spherex_image_size] = exp1[1].data
        canvas[:, spherex_image_size+gap_width:2*spherex_image_size+gap_width] = exp2[1].data
        canvas[:, 2*spherex_image_size+2*gap_width:3*spherex_image_size+2*gap_width] = exp3[1].data
        canvas_shape = canvas.shape
        canvas[canvas == 0] = np.nan

        # Transform from MJy sr-1 to Jy/arcsec2
        # canvas = canvas * rs.constants.MJysr_to_Jyarcsec2

        header = rs.utils.create_custom_wcs(crpix=[int(canvas_shape[1]/2), int(canvas_shape[0]/2)],
                                            crval=[ra,dec],
                                            cdelt=[-rs.telescopes.SPHEREx.get_pixscale(instrument="SPHEREx").to("degree").value, rs.telescopes.SPHEREx.get_pixscale(instrument="SPHEREx").to("degree").value],
                                            crota=[-pa,-pa])
        rs.utils.save_fits(canvas, outname, header)
        return(outname)



    # SPHEREx mu2fe and fe2mu

    def mu2fe(mu):
        jy_arcsec2 = 10**(-0.4*(mu-8.9))
        MJysr = jy_arcsec2 / rs.constants.MJysr_to_Jyarcsec2
        cal = 8 # e/s sr/MJy
        fe = MJysr*cal
        return(fe/u.second)


    def fe2mu(fe):
        cal = 8 # e/s sr/MJy
        MJysr = fe/cal
        jy_arcsec2 = MJysr * rs.constants.MJysr_to_Jyarcsec2
        mu = -2.5*np.log10(jy_arcsec2) + 8.9
        return(mu)


    def get_PSF():
        # /Users/aborlaff/NASA/ROSALIA/rosalia/CORE/PSF_ARCHIVE/SPHEREx_psf.fits
        psf = fits.open(os.environ["ROSALIACACHE"] + "/CORE/PSF_ARCHIVE/SPHEREx_psf.fits")
        return(psf[0].data)


#################################################################

class ARRAKIHS:
    """ ARRAKIHS Telescope properties class"""
    TELESCOP = "ARRAKIHS"
    mirror_radius = 0.075*u.meter

    def TLE_exposure(epoch):
        #return(rs.satellites.find_closest_TLE(epoch=epoch, TLE_history=CSST.load_CSST_TLE()))
        return(rs.telescopes.ARRAKIHS.load_example_TLE())

    def load_example_TLE():
        from skyfield.api import EarthSatellite
        ts = sf_api.load.timescale()
        # 1 99999U          25001.00000000  .00000000  00000-0  00000-0 0 00006
        # 2 99999 098.6066 190.9024 0019445 243.2152 354.9930 14.25484048052155

        line1 = '1 99999U          25001.00000000  .00000000  00000-0  00000-0 0 00006'
        line2 = '2 99999 098.6066 190.9024 0019445 243.2152 354.9930 14.25484048052155'
        satellite = EarthSatellite(line1, line2, 'ARRAKIHS', ts)
        return(satellite)

    def get_canvas_shape(instrument):
        if instrument == "ARRAKIHS":
            canvas_shape = np.array(np.array([1.18, 1.18])*60*60/1.5*u.arcsec).astype("int")
        return(canvas_shape)


    def get_pixscale(instrument):
        if ("ARRAKIHS" in instrument) or ("ARRAKIHS" in instrument):
            pixscale = 1.5*u.arcsec
        return(pixscale)



    def get_filter(instrument, filter_name, verbose=False):
        return(rs.telescopes.get_filter(telescope="Euclid", instrument="VIS",
                                        detector="VIS", filter_name="VIS", verbose=False))

    def old_make_dummy_exposure(ra, dec, pa, outname):
        canvas_shape = rs.telescopes.ARRAKIHS.get_canvas_shape(instrument="ARRAKIHS")
        canvas = np.zeros(canvas_shape)
        header = rs.utils.create_custom_wcs(crpix=[int(canvas_shape[1]/2), int(canvas_shape[0]/2)],
                                            crval=[ra,dec],
                                            cdelt=[-rs.telescopes.ARRAKIHS.get_pixscale(instrument="ARRAKIHS").to("degree").value,
                                                    rs.telescopes.ARRAKIHS.get_pixscale(instrument="ARRAKIHS").to("degree").value],
                                            crota=[-pa,-pa])
        rs.utils.save_fits(canvas, outname, header)
        return(outname)

    def modified_make_dummy_exposure(ra, dec, pa, outname):
        canvas_shape = [4001, 4001] # rs.telescopes.ARRAKIHS.get_canvas_shape(instrument="ARRAKIHS")
        canvas = np.zeros(canvas_shape)
        header = rs.utils.create_custom_wcs(crpix=[int(canvas_shape[1]/2), int(canvas_shape[0]/2)],
                                            crval=[ra,dec],
                                            cdelt=[-rs.telescopes.ARRAKIHS.get_pixscale(instrument="ARRAKIHS").to("degree").value,
                                                    rs.telescopes.ARRAKIHS.get_pixscale(instrument="ARRAKIHS").to("degree").value],
                                            crota=[-pa,-pa])
        rs.utils.save_fits(canvas, outname, header)
        return(outname)


    def make_dummy_exposure(ra, dec, pa, outname):
        exp_name = "/Users/aborlaff/NASA/SPARKLES/notebooks/SATELLITES/ARRAKIHS_mock_es.fits"

        exp = fits.open(exp_name)
        canvas = exp[0].data
        canvas_shape = canvas.shape
        cdelt = [-1.5/60/60, 1.5/60/60]

        header = rs.utils.create_custom_wcs(crpix=[int(canvas_shape[1]/2), int(canvas_shape[0]/2)],
                                            crval=[ra,dec],
                                            cdelt=cdelt,
                                            crota=[-pa,-pa])
        rs.utils.save_fits(canvas, outname, header)
        return(outname)

    def get_PSF():
        sigma_ARRAKIHS_psf = 0.756/1.5
        print("ARRAKIHS Moffat 2D PSF")
        #from astropy.convolution import Gaussian2DKernel
        #gaussian_2D_kernel = Gaussian2DKernel(sigma_ARRAKIHS_psf, x_size=101, y_size=101)
        from astropy.modeling.functional_models import Moffat2D

        a = Moffat2D(amplitude=1, gamma=1)
        x = np.linspace(-20,20,41)
        y = np.linspace(-20,20,41)
        X, Y = np.meshgrid(y,x)
        psf = a.evaluate(amplitude=1, gamma=0.446, x_0=0, y_0=0, alpha=2.5, x=X, y=Y)
        return(psf)
        #return(gaussian_2D_kernel.array)

    def mu2fe(mu):
        return(rs.detectors.mu2fe(mu=mu, instrument="ARRAKIHS", filter_name="VIS", telescope="ARRAKIHS", verbose=False))

    def fe2mu(fe):
        return(rs.detectors.fe2mu(fe=fe, instrument="ARRAKIHS", filter_name="VIS", telescope="ARRAKIHS", verbose=False))

#################################################################


class Euclid:
    """ Euclid Telescope properties class"""
    TELESCOP = "Euclid"
    mirror_radius = 0.6*u.meter
