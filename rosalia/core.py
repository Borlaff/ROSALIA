# Telescope definition object 
# telescope = {"pointing": [RA_TARG, DEC_TARG], "PA_Y": PA_Y, "detector_shape": [NAXIS1, NAXIS2], "pixscale": pixscale,
#              "R_mirror": R_mirror, "OBSLOC": OBSLOC, "EXPSTART": EXPSTART, "EXPTIME": EXPTIME}

from astropy.time import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import rosalia as rs
from astropy.coordinates import SkyCoord
from rosalia.correct import rosalia_stray
import os
from datetime import datetime
import pandas as pd

class exposure():
    """
    Exposure class. This class holds all the information about a particular image, including the pointing orientation,  the filter, the position angle, the date of the observation, and many other parameters. The ``exposure`` class also contains methods to estimate the stray-light signal in a given exposure, and to correct for it.

    Exposure can be initialized in two ways: 
    - By providing a filename of a fits file with the necessary WCS and header information (filename). In this case, the class will read the fits file and extract the necessary information from the header. This is the most straightforward way to initialize the class, as long as the input fits file has the necessary information in the header.

    - By providing the custom observer parameters in a dictionary (observer). In this case, the class will use the observer parameters to create a dummy fits file with the necessary WCS and header information. This is useful when the user does not have a fits file with the necessary information, but has the observer parameters. If the observer keyword "TELESCOP" is set to "Roman/WFI", then the class will create a dummy Roman/WFI fits file with the necessary information. In that case, only needs to contain "TELESCOP", "pointing", "FILTER", "PA_Y", "EXPSTART", and "EXPTIME". For example, the observer dictionary can be defined as follows:

    
    observer={"TELESCOP": "Roman/WFI", 
            "pointing": [ra, dec], 
            "FILTER":bandpass, 
            "PA_Y": PA, 
            "EXPSTART": date.mjd, 
            "EXPTIME": exptime}

    prefix = "Pleiades"
    custom_roman_exposure = rs.core.exposure(observer=observer, prefix=prefix) 
    straylight_out = custom_roman_exposure.straylight() 
     
     
    If the observer keyword "TELESCOP" is set to any other value, then the class will create a generic dummy fits file with the necessary information. The observer dictionary should contain the following parameters:

    TELESCOP = "MYTELESCOPE"
    INSTRUME = "MYTELESCOPE"
    DETECTOR = "MYTELESCOPE"
    RA_TARG = 100
    DEC_TARG = 23
    PA_Y = 32
    NAXIS1=3048
    NAXIS2=4096
    PIXSCALE=8/60/60
    R_mirror=1
    OBSLOC=273
    EXPSTART=61000
    EXPTIME=10
    FILTER_PARAMS={"NAME": "F606W", "TELESCOPE": "HST", "INSTRUMENT": "ACS", "DETECTOR": "WFC"}

    observer = {"TELESCOP":TELESCOP, 
                "INSTRUME":INSTRUME, 
                "DETECTOR": DETECTOR,
                "pointing": [RA_TARG, DEC_TARG], 
                "PA_Y": PA_Y, 
                "DATA_SHAPE": [NAXIS1, NAXIS2], 
                "PIXSCALE": PIXSCALE,
                "R_mirror": R_mirror, 
                "FILTER_PARAMS":FILTER_PARAMS, 
                "OBSLOC": OBSLOC, 
                "EXPSTART": EXPSTART, 
                "EXPTIME": EXPTIME}

    custom_exposure = rs.core.exposure(observer=observer)

    ---------- 
    
    init parameters:
    - filename: str. Path to the input fits file. If None, then the user needs to provide the observer and telescope parameters.
    - prefix: str. Prefix to add to the output filename. Default is "".
    - observer: dict. Dictionary with the observer parameters. If None, then the user needs to provide the filename parameter. The dictionary should contain the following
        - "TELESCOP": str. For predefined telescopes, this replaces the telescope dictionary, simplifying the input. Name of the telescope. For example, "Roman/WFI". 
        - "pointing": list. List with the RA and DEC of the pointing, in degrees. For example, [RA_TARG, DEC_TARG].
        - "FILTER": str. Name of the filter. For example, "F129".
        - "PA_Y": float. Position angle of the Y axis of the detector, in degrees. For example, PA_Y = 0 means that the Y axis of the detector is aligned with the North direction.
        - "EXPSTART": float. Start time of the exposure, in MJD. For example, EXPSTART = 59300.0.
        - "EXPTIME": float. Exposure time, in seconds. For example, EXPTIME = 600.0.
        
    """


    def __init__(self, filename=None, prefix="", observer=None, telescope=None):
        import os 
        from astropy.wcs import WCS   

        if observer is not None:
            if observer["TELESCOP"] == "Roman/WFI":
                self = self.roman_wfi_exposure(observer, prefix=prefix)

            if telescope is not None: 
                self.TELESCOP = telescope['TELESCOP']
                self.INSTRUME = telescope['INSTRUME']
                self.DETECTOR = telescope['DETECTOR']
                self.DATA_SHAPE = telescope['DATA_SHAPE']
                self.RA_TARG = telescope["pointing"][0]
                self.DEC_TARG = telescope["pointing"][1]
                self.PA = telescope['PA_Y']
                self.EXPSTART = telescope['EXPSTART']
                self.EXPTIME = telescope['EXPTIME']
                self.EXPMID = (Time(self.EXPSTART, format="mjd") + self.EXPTIME*u.s/2).mjd
                self.EXPEND = (Time(self.EXPSTART, format="mjd") + self.EXPTIME*u.s).mjd
                # self.PHYSPIX = telescope['PHYSPIX']
                self.PIXSCALE = telescope['PIXSCALE']
                # self.EXPSTART_ISOT = 
                self.MPC_OBSLOC = self.get_mpc_observer_location()
                self.JPL_OBSLOC = self.get_jpl_observer_location()
                self.FILTER_IDENTITY = rs.telescopes.find_filter_in_svo(wavelength=telescope["FILTER_PARAMS"]["NAME"],
                                                                        telescope=telescope["FILTER_PARAMS"]["TELESCOPE"],
                                                                        instrument=telescope["FILTER_PARAMS"]["INSTRUMENT"],
                                                                        detector=telescope["FILTER_PARAMS"]["DETECTOR"], verbose=False)
                # self.XYZ_HELIO_POS = exposure_identity['XYZ_HELIO_POS']


                self.SCIEXTS = [0]
                header = rs.utils.create_custom_wcs(crpix=np.array(self.DATA_SHAPE)/2, 
                                                            crval=[self.RA_TARG, self.DEC_TARG], 
                                                            cdelt=[-self.PIXSCALE,self.PIXSCALE], 
                                                            crota=[-self.PA,-self.PA], 
                                                            projection="TAN")
            
                header["NAXIS1"] = self.DATA_SHAPE[0]
                header["NAXIS2"] = self.DATA_SHAPE[1]

                self.ASTROPYWCS = [WCS(header)]
                self.FPA_NEAR_RADIUS = self.get_max_angular_size()

                if prefix != "":
                    prefix = prefix + "_"
                self.FILENAME = os.getcwd() + "/" + prefix + self.TELESCOP + "_RA_" + '{:07.3f}'.format(self.RA_TARG) +\
                                                "_DEC_" + '{:07.3f}'.format(self.DEC_TARG) +\
                                                "_MJD_" + '{:07.5f}'.format(self.EXPSTART) +\
                                                "_PA_" + '{:06.2f}'.format(self.PA) + ".fits"

        if filename is not None:
            exposure_identity = rs.utils.exposure_inspector(filename, lite=False)
            self.DATA = exposure_identity['DATA']
            self.FILENAME = exposure_identity['FILENAME']
            self.TELESCOP = exposure_identity['TELESCOP']
            self.INSTRUME = exposure_identity['INSTRUME']
            self.DETECTOR = exposure_identity['DETECTOR']
            self.RA_TARG = exposure_identity['RA_TARG']
            self.DEC_TARG = exposure_identity['DEC_TARG']
            self.EXPSTART = exposure_identity['EXPSTART']
            self.EXPTIME = exposure_identity['EXPTIME']
            self.EXPMID = (Time(self.EXPSTART, format="mjd") + self.EXPTIME*u.s/2).mjd
            self.EXPEND = (Time(self.EXPSTART, format="mjd") + self.EXPTIME*u.s).mjd
            self.BUNIT = exposure_identity['BUNIT']
            self.EXPSTART_ISOT = exposure_identity['EXPSTART_ISOT']
            self.PA = exposure_identity['PA']
            # self.SCA = exposure_identity['SCA']
            # self.HST_TYPE = exposure_identity['HST_TYPE']
            self.FILTER = exposure_identity['FILTER']
            self.FILTER_IDENTITY = exposure_identity['FILTER_IDENTITY']
            # self.PHYSPIX = exposure_identity['PHYSPIX']
            self.PIXSCALE = exposure_identity['PIXSCALE']
            self.SCIEXTS = exposure_identity['SCIEXTS']
            self.DATA_SHAPE = exposure_identity['DATA_SHAPE']
            self.ASTROPYWCS = exposure_identity['ASTROPYWCS']
            self.FILETYPE = exposure_identity['FILETYPE']
            self.MPC_OBSLOC = rs.horizons.get_mpc_observer_name(self.TELESCOP)
            self.JPL_OBSLOC = rs.horizons.get_jpl_observer_name(self.TELESCOP)
            # self.XYZ_HELIO_POS = exposure_identity['XYZ_HELIO_POS']
            self.FPA_NEAR_RADIUS = self.get_max_angular_size()

    def roman_wfi_exposure(self, observer, prefix=""):
        print("> roman_wfi_exposure")
        # Here we expect observer={"TELESCOP": "Roman/WFI", "pointing": [RA_TARG, DEC_TARG], "FILTER":FILTER, "PA_Y": PA_Y, "EXPSTART": EXPSTART, "EXPTIME": EXPTIME}
        # Fixed parameters for Roman/WFI 
        observer['TELESCOP'] = "Roman"
        observer['INSTRUME'] = "WFI"
        observer['DETECTOR'] = "WFI"
        observer['DATA_SHAPE'] = [4088, 4088]
        observer['PIXSCALE'] = rs.telescopes.Roman.get_pixscale(instrument=observer['INSTRUME']).to("degree").value
        observer['FILTER_PARAMS'] = {"NAME": observer["FILTER"], "TELESCOPE": "RST", "INSTRUMENT": "WFI", "DETECTOR": "WFI"}
        
        # Generic derived definitions. 
        self.TELESCOP = observer['TELESCOP']
        self.INSTRUME = observer['INSTRUME']
        self.DETECTOR = observer['DETECTOR']
        self.DATA_SHAPE = observer['DATA_SHAPE']
        self.RA_TARG = observer["pointing"][0]
        self.DEC_TARG = observer["pointing"][1]
        self.PA = observer['PA_Y']
        self.EXPSTART = observer['EXPSTART']
        self.EXPTIME = observer['EXPTIME']
        self.EXPTIME = observer['EXPTIME']
        self.EXPMID = (Time(self.EXPSTART, format="mjd") + self.EXPTIME*u.s/2).mjd
        self.EXPEND = (Time(self.EXPSTART, format="mjd") + self.EXPTIME*u.s).mjd
        self.PHYSPIX = rs.telescopes.Roman.get_physical_pixelsize(instrument=observer['INSTRUME'])
        self.PIXSCALE = rs.telescopes.Roman.get_pixscale(instrument=observer['INSTRUME'])
        self.EXPSTART_ASTROPY = Time(self.EXPSTART, format="mjd")
        self.EXPSTART_ISOT = self.EXPSTART_ASTROPY.isot
        self.MPC_OBSLOC = rs.horizons.get_mpc_observer_name(observer['TELESCOP'])
        self.JPL_OBSLOC = rs.horizons.get_jpl_observer_name(observer['TELESCOP'])

        self.FILTER_IDENTITY = rs.telescopes.find_filter_in_svo(wavelength=observer["FILTER_PARAMS"]["NAME"],
                                                                telescope=observer["FILTER_PARAMS"]["TELESCOPE"],
                                                                instrument=observer["FILTER_PARAMS"]["INSTRUMENT"],
                                                                detector=observer["FILTER_PARAMS"]["DETECTOR"], verbose=False)
        # self.XYZ_HELIO_POS = exposure_identity['XYZ_HELIO_POS']

        if "FILENAME" not in observer:
            # If the user did not define an output filename, do it for them
            if prefix != "":
                prefix = prefix + "_"
            observer["FILENAME"] = os.getcwd() + "/" + prefix + self.TELESCOP +\
                                    "_RA_" + '{:07.3f}'.format(self.RA_TARG) +\
                                    "_DEC_" + '{:07.3f}'.format(self.DEC_TARG) +\
                                    "_MJD_" + '{:07.5f}'.format(self.EXPSTART) +\
                                    "_PA_" + '{:06.2f}'.format(self.PA) + ".fits"
            self.FILENAME = observer["FILENAME"]

        central_coords = SkyCoord(self.RA_TARG, self.DEC_TARG, frame="icrs", unit="deg")

        observer["FILENAME"] = rs.roman.create_roman_dummy(point=central_coords, date=self.EXPSTART_ASTROPY,
                                                            band=observer["FILTER_PARAMS"]["NAME"],
                                                            PA=self.PA, exptime=self.EXPTIME,
                                                            output=observer["FILENAME"])
        
        self.SCIEXTS = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
        astropywcs_info = rs.utils.get_astropywcs_info_from_sciexts(filename=self.FILENAME, sciexts=self.SCIEXTS)
        self.ASTROPYWCS = astropywcs_info["ASTROPYWCS"]
        self.DATA_SHAPE = astropywcs_info["DATA_SHAPE"]
        self.PIXSCALE = astropywcs_info["PIXSCALE"]
        self.FPA_NEAR_RADIUS = self.get_max_angular_size()
        #    exposure_identity["HEADERS"] = header_list
        #    exposure_identity["DATA_SHAPE"] = data_shape
        #    exposure_identity["ASTROPYWCS"] = astropywcs
        #    exposure_identity["PIXSCALE"] = np.abs(astropywcs[0].proj_plane_pixel_scales()[0])


        return(self)
    
    def get_detector_corners(self):
        detector_corners = []
        for i in range(len(self.SCIEXTS)):
            detector_corners.append(rs.detectors.get_detector_corners(self.ASTROPYWCS[i]))
        return(detector_corners)
    

    def get_max_angular_size(self):
        return(rs.utils.find_max_angular_size_of_image(wcs=self.ASTROPYWCS, ra_cen=self.RA_TARG, dec_cen=self.DEC_TARG))
    

    def plot_footprint(self, figsize=(10,10), verbose=True, ax=None, color='red', label=None):
        print("Hey! plot_footprint")
        if verbose: print("Finding ra dec constraints")
        # ra_dec_constraints = rs.gaia.find_ra_dec_constraints(self.RA_TARG, self.DEC_TARG, radius=self.FPA_NEAR_RADIUS, verbose=verbose)
        if verbose: print("Getting detector corners")
        exp_corners = self.get_detector_corners()
        if verbose: print("Plotting...")
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        for i in range(len(self.SCIEXTS)):
            x = np.array(exp_corners[i]["corners_world"][:,0].tolist() + [exp_corners[i]["corners_world"][0,0]])
            y = np.array(exp_corners[i]["corners_world"][:,1].tolist() + [exp_corners[i]["corners_world"][0,1]])
            ax.plot(x, y, color=color, label=label)
            # ax.set_xlim(ra_dec_constraints["ra_max"], ra_dec_constraints["ra_min"])
            # ax.set_ylim(ra_dec_constraints["dec_min"], ra_dec_constraints["dec_max"])
            ax.set_xlabel("Right Ascension (degree)")
            ax.set_ylabel("Declination (degree)")
        return(ax)
    

    def find_nearby_ssos(self):
        import ephessos as ep 
        sso_cone_search = ep.core.cone_search(ra=self.RA_TARG, dec=self.DEC_TARG, mjd=self.EXPSTART, search_radius=self.FPA_NEAR_RADIUS, verbose=False)
        return(sso_cone_search)


    def get_nearby_sources(self, g_mag_max=15, verbose=False):
        hybrid_catalog = rs.psf.get_hybrid_catalog(ra=self.RA_TARG, dec=self.DEC_TARG,
                                                   radius=1,
                                                   lambda_ref=self.FILTER_IDENTITY["filter_lambda_ref"],
                                                   MJD=self.EXPSTART,
                                                   observer=self.TELESCOP,
                                                   g_mag_max = g_mag_max,
                                                   verbose=verbose,
                                                   query_filename=self.FILENAME.replace(".fits", ".csv"))
        return(hybrid_catalog)
        

    def get_source_catalog(self, g_mag_max=15, verbose=False):
        import pandas as pd

        self.source_catalog_filename = os.path.splitext(os.path.basename(self.FILENAME))[0] + "_source_catalog.csv" #
        
        search_radius = 1.5*self.get_max_angular_size()

        if os.path.exists(self.source_catalog_filename):
            print("WARNING: Loading existing catalog! Remove " + self.source_catalog_filename + " if this is a mistake.")
            hybrid_catalog = pd.read_csv(self.source_catalog_filename)

        else:
            if verbose: print("> Querying stars in the surroundings using ESA/Gaia Archive")

            if search_radius > 0.5:
                print("INFO: radius parameter (minimum distance to search for individual stars) is > 0.5 degrees.")
                print("Gaia/2MASS/WISE query database can take several minutes to process. Please be patient.")

            # Find the stars around the central coordinate of the scene.
            # Save the stellar catalog into a catalog object, and run main_offender as if input_catalog was set by the User.
            loader = rs.plots.Loader("Querying Gaia/2MASS/WISE/JPL Horizons databases. This might take a few minutes...",
                                    "All-sky source map constructed.", 0.05).start()

            hybrid_catalog = rs.psf.get_hybrid_catalog(ra=self.RA_TARG, dec=self.DEC_TARG,
                                                    radius=search_radius,
                                                    lambda_ref=self.FILTER_IDENTITY["filter_lambda_ref"],
                                                    MJD=self.EXPMID,
                                                    observer=self.TELESCOP,
                                                    g_mag_max = g_mag_max,
                                                    verbose=verbose,
                                                    query_filename=self.source_catalog_filename)
            loader.stop()
        return(hybrid_catalog)



    def find_which_stars_are_inside_each_detector(self, verbose=False):
        # Find where each star lands (detector ID or outside FOV)

        from tqdm import tqdm
        names_of_bool_columns_if_star_is_inside = []
        detector_square_list = []

        """
        If the input file is a multi-extension fits, then exposure_inspector will scan for extensions with EXTNAME = SCI.
        The extension ID in the FITS file will be stored in SCIEXTS = image_identity["SCIEXTS"].

        In that case, hybrid_catalog, the catalog of stars, will have a set of N columns called in_SCI[i] (boolean), where the catalog
        stores if that particular star is inside each detector or not.

        """
        
        fpa_detector_corners = self.get_detector_corners()
        for SCIEXT_i, ASTROPYWCS_i, detector_corners in tqdm(zip(self.SCIEXTS, self.ASTROPYWCS, fpa_detector_corners)):
            if verbose: print("> Identifying which stars are inside the FOV and which are outside...")
            infield_stars = rs.psf.identify_stars_in_out_field(data_shape=self.DATA_SHAPE,
                                                            wcs=ASTROPYWCS_i,
                                                            catalog=self.source_catalog,
                                                            verbose=verbose)

            name_column_is_star_inside_this_detector = "in_SCI" + str(SCIEXT_i)
            names_of_bool_columns_if_star_is_inside.append(name_column_is_star_inside_this_detector)

            self.source_catalog[name_column_is_star_inside_this_detector] = infield_stars["bool_isIn"]
            # detector_corners = fpa_detector_corners[SCIEXT_i-1]
            # If verbose, make a plot of the stars with the footprint.
            # rs.detectors.get_detector_corners(wcs=self.ASTROPYWCS[SCIEXT_i-1])
            detector_square_list.append(np.concatenate([detector_corners["corners_world"], detector_corners["corners_world"]]))
        # Once you are done checking if the stars are inside each detector,
        # find out which stars are outside ALL detectors.
        self.source_catalog["is_inside_FPA"] = self.source_catalog[names_of_bool_columns_if_star_is_inside].any(axis=1)
        return(self.source_catalog)


    def get_nearby_ssos(self, ra=None, dec=None, radius=None, mjd=None, verbose=False, time_step="30s"):
        import ephessos as ep
        if ra==None: ra = self.RA_TARG
        if dec==None: dec = self.DEC_TARG
        if radius==None: radius = self.FPA_NEAR_RADIUS*60*60
        if mjd==None: mjd = self.EXPSTART

        cone_search = ep.core.cone_search(ra=ra, dec=dec, mjd=mjd, 
                                          search_radius=radius, observatory=self.MPC_OBSLOC, verbose=verbose)

        if len(cone_search) == 0:
            print("No SSOs found!")
            return()
        else:        
            ephessos_df = ep.core.ephessos(sso_search=cone_search, mjd_start=self.EXPSTART, mjd_end=self.EXPEND, 
                                           obs_center=self.JPL_OBSLOC, step_size=time_step, verbose=verbose)
        return(cone_search, ephessos_df)
    

    def make_generic_dummy(self, binning=1):
        data = []
        headers = []
        for SCIEXT in self.SCIEXTS:
            # Create a new WCS object.  The number of axes must be set from the start
            w = self.ASTROPYWCS[SCIEXT]
    
            # Setting some dummy image
            image_data = np.random.normal(loc=0, scale=1, 
                                          size=[int(self.DATA_SHAPE[0]/binning),
                                                int(self.DATA_SHAPE[1]/binning)])

            # Now, write out the WCS object as a FITS header
            header = w.to_header()
            header["EXTNAME"] = "SCI"
            header["RA_TARG"] = self.RA_TARG
            header["DEC_TARG"] = self.DEC_TARG
            header["EXPSTART"] = self.EXPSTART
            header["EXPEND"] = self.EXPEND
            header["PA"] = self.PA
            header["TELESCOP"] = self.TELESCOP
            header["INSTRUME"] = self.INSTRUME
            header["DETECTOR"] = self.DETECTOR

            # header is an astropy.io.fits.Header object.  We can use it to create a new
            # PrimaryHDU and write it to a file.
            data.append(image_data)
            headers.append(header)
        
        rs.utils.save_fits(data, self.FILENAME, headers)
        return(self.FILENAME)



    from concurrent.futures import ProcessPoolExecutor
    from tqdm import tqdm

    @staticmethod
    def _parallel_worker(args):
        data_shape, wcs, SCIEXT_i, filter_identity, ra_stars_outside, dec_stars_outside, cat_id_outside, source_id_outside, irradiance_stars, ra_point, dec_point, pa_point, verbose = args
        return rs.roman.roman_estimate_straylight_SCA(
            data_shape=data_shape,
            wcs=wcs,
            SCA=SCIEXT_i,
            filter_identity=filter_identity,
            ra_stars=ra_stars_outside,
            dec_stars=dec_stars_outside,
            cat_id=cat_id_outside,
            source_id=source_id_outside,
            irradiance_stars=irradiance_stars,
            ra_point=ra_point,
            dec_point=dec_point,
            pa_point=pa_point,
            verbose=verbose
        )

    def straylight(self, catalog=None, g_mag_max=15, sun_block=False, verbose=False):
        from astropy import constants as const
        from tqdm import tqdm 
        from astropy.io import fits
        
        #######################################
        # straylight: AKA. main_offender: Alejandro S. Borlaff. NASA/Ames STA. a.s.borlaff@nasa.gov
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
        # v5 - 4 June 2026 - Adapted into a method of the exposure class. Now it can be used for any exposure, but the straylight estimation is currently only available for Roman/WFI exposures.
        #######################################

    
        if self.TELESCOP != "Roman" and self.TELESCOP != "RST" and self.TELESCOP != "ROMAN":
            print("Straylight modeling is currently only available for Roman/WFI exposures.")
            return(None)
        #    stray_db = rosalia_stray(ra=self.RA_TARG, dec=self.DEC_TARG, PA=self.PA,
        #                date=self.EXPSTART, bandpass=self.FILTER_IDENTITY["NAME"], 
        #                exptime=self.EXPTIME, prefix=prefix, input_fits=input_fits, radius=radius,
        #                g_mag_max=g_mag_max, sun_block=sun_block, verbose=verbose, catalog=catalog, 
        #                figsize=figsize, mu_vmin=mu_vmin, mu_vmax=mu_vmax)

        # Find out which stars belong to each detector. 
        print(datetime.now().isoformat() + " > Fetching catalog")

        if catalog is None:
            self.source_catalog = self.get_source_catalog(g_mag_max=g_mag_max, verbose=False)
        else:
            self.source_catalog = rs.utils.fix_custom_catalog(catalog)

        self.source_catalog = self.find_which_stars_are_inside_each_detector(verbose=False)
        print(datetime.now().isoformat() + " > Done.")

        # If sun_block is True, then remove the Sun from the catalog.
        if sun_block:
            self.source_catalog = self.source_catalog[~(self.source_catalog["source_id"] == "Sun")]
            self.source_catalog = self.source_catalog[~(self.source_catalog["source_id"] == "Earth")]


        straylevel_list = []
        main_offender_list = []

        ## Prepare the coordinates of the stars that do not fall inside the Focal Plane Array ##
        ## This step is common for all SCI extensions #
        ra_stars_outside      = np.array(self.source_catalog["ra"][~self.source_catalog["is_inside_FPA"]])
        dec_stars_outside     = np.array(self.source_catalog["dec"][~self.source_catalog["is_inside_FPA"]])
        source_id_outside     = np.array(self.source_catalog["source_id"][~self.source_catalog["is_inside_FPA"]])
        cat_id_outside        = np.array(self.source_catalog["cat_id"][~self.source_catalog["is_inside_FPA"]])
        synthetic_mag_outside = np.array(self.source_catalog["mag_lambda"][~self.source_catalog["is_inside_FPA"]])
        # stars_world_location = coordinates.SkyCoord(ra_stars_outside, dec_stars_outside, frame='icrs', unit="deg")

        lambda_max = self.FILTER_IDENTITY["filter_lambda_max"]
        lambda_min = self.FILTER_IDENTITY["filter_lambda_min"]
        lambda_ref = self.FILTER_IDENTITY["filter_lambda_ref"]
        irradiance_stars = const.c*(lambda_max-lambda_min)/(lambda_ref**2)*((10**(-0.4*(synthetic_mag_outside+56.1)))*u.W/u.meter**2/u.Hz)

        ##########################################
        # Define the output tables and filenames #
        ##########################################

        # Save the stray-light full scale map
        if "s3://" in self.FILENAME:
            # 's3://stpubdata/roman/nexus/soc_simulations/tutorial_data/roman-2026.1/r0003201001001001004_0001_wfi01_f106_cal.asdf'
            self.output_name = self.FILENAME.split("/")[-1].replace(".asdf", "_stray.fits")
        else: 
            self.output_name = self.FILENAME.replace(".fits", "_stray.fits")

        self.main_offender_output_name = self.output_name.replace(".fits", "_main_off.fits")
        self.straylight_db_output_name = self.output_name.replace(".fits", "_db.csv")

 

        ########################################
        # Here we estimate the stray-light
        ########################################
        # Reset the Roman / WFI loading bar:
        # rs.plots.ascii_progress_focal_plane.canvas = np.copy(rs.plots.ascii_progress_focal_plane.canvas_zero)
        straylevel_all_SCAS = []


        # Running parallel computation # 
        straylevel_all_SCAS = []
        t = datetime.now()
        print(datetime.now().isoformat() + " > Starting Stray-light scan: ")
        from concurrent.futures import ProcessPoolExecutor
        from tqdm import tqdm

        with ProcessPoolExecutor() as executor:
            inputs = [
                (
                    DATA_SHAPE_i,
                    ASTROPYWCS_i,
                    SCIEXT_i,
                    self.FILTER_IDENTITY,
                    ra_stars_outside,
                    dec_stars_outside,
                    cat_id_outside,
                    source_id_outside,
                    irradiance_stars,
                    self.RA_TARG,
                    self.DEC_TARG,
                    self.PA,
                    verbose,
                )
                for SCIEXT_i, DATA_SHAPE_i, ASTROPYWCS_i in zip(self.SCIEXTS, self.DATA_SHAPE, self.ASTROPYWCS)
            ]
            results = list(tqdm(executor.map(self._parallel_worker, inputs),total=len(inputs),))
            straylevel_all_SCAS.extend(results)

        print(datetime.now().isoformat() + " > Done : " + str(datetime.now() - t) + " elapsed.")

        ###########
        # Reconstruct the stray-light maps and main offender maps from the straylevel_all_SCAS database. 
        ###########
        NSCAs = len(self.SCIEXTS)# NSCAs 
        print(" > Reconstructing the Stray-light / Main offender map: ")
        straylevel_list = [] 
        main_offender_list = [] 
        for SCA in range(len(self.SCIEXTS)):
                # This is the canvas array where we will store all the straylight level.
            straylight_SCA = np.zeros(self.DATA_SHAPE[0]).astype(np.float32)
            # This is the canvas array where we will store the ID of the largest stray-light contributor
            main_offender_SCA = np.zeros(self.DATA_SHAPE[0]).astype(np.float32)

            for subarray_i in range(len(straylevel_all_SCAS[SCA])):
                xmin = straylevel_all_SCAS[SCA]["xmin"].iloc[subarray_i]
                xmax = straylevel_all_SCAS[SCA]["xmax"].iloc[subarray_i]
                ymin = straylevel_all_SCAS[SCA]["ymin"].iloc[subarray_i]
                ymax = straylevel_all_SCAS[SCA]["ymax"].iloc[subarray_i]
                straylight_SCA[ymin:ymax, xmin:xmax] = straylevel_all_SCAS[SCA]["straylight_total"].iloc[subarray_i]
                main_offender_SCA[ymin:ymax, xmin:xmax] = straylevel_all_SCAS[SCA]["mainoffender_total"].iloc[subarray_i]

            
            straylevel_list.append(straylight_SCA)
            main_offender_list.append(main_offender_SCA)

        straylevel_db = pd.concat(straylevel_all_SCAS)
        straylevel_db.to_csv(self.straylight_db_output_name)

        print(datetime.now().isoformat() + " > Done : ")


        ########################################
        # Save the results to a fits file.
        ########################################

        # Stray-light
        data_output = []
        header_output = []
        for SCIEXT_i, straylevel_image_i, ASTROPYWCS_i in tqdm(zip(self.SCIEXTS, straylevel_list, self.ASTROPYWCS)):
            data_output.append(straylevel_image_i)
            header_output.append(ASTROPYWCS_i.to_header())
        
        rs.utils.save_fits(array=data_output, name=self.output_name, header=header_output,
                           extname=None, overwrite=True, output_verify='silentfix')

        # Main-offender
        data_output = []
        header_output = []
        for SCIEXT_i, main_offender_i, ASTROPYWCS_i in tqdm(zip(self.SCIEXTS, main_offender_list, self.ASTROPYWCS)):
            data_output.append(main_offender_i)
            header_output.append(ASTROPYWCS_i.to_header())


        rs.utils.save_fits(array=data_output, 
                        name=self.main_offender_output_name, 
                        header=header_output,
                        extname=None, 
                        overwrite=True, 
                        output_verify='silentfix')


        # Let's do one more step to include the needed metadata from the dummy file. 
        stray_image = fits.open(self.output_name)
        main_offender_image = fits.open(self.main_offender_output_name)
        #exposure_identity_keys_to_copy = ["TELESCOP", "INSTRUME", "DETECTOR", "FILTER", "RA_TARG", "DEC_TARG", 
        #                                  "RA_PNT", "DEC_PNT", "X_PNT", "Y_PNT", "PA",
        #                                  "EXPTIME", "EXPSTART", "EXPSTART_ISOT"]
        #for key in exposure_identity_keys_to_copy:
        
        stray_image[0].header["TELESCOP"] = self.TELESCOP
        stray_image[0].header["INSTRUME"] = self.INSTRUME
        stray_image[0].header["DETECTOR"] = self.DETECTOR
        stray_image[0].header["FILTER"] = self.FILTER_IDENTITY["wavelength"]
        stray_image[0].header["RA_TARG"] = self.RA_TARG
        stray_image[0].header["DEC_TARG"] = self.DEC_TARG
        stray_image[0].header["PA"] = self.PA
        stray_image[0].header["EXPTIME"] = self.EXPTIME
        stray_image[0].header["EXPSTART"] = self.EXPSTART
        stray_image[0].header["EXPSTART_ISOT"] = self.EXPSTART_ISOT

        stray_image[0].header["WAVEREF"] = self.FILTER_IDENTITY["filter_lambda_ref"].to("nm").value
        stray_image[0].header["WAVEMIN"] = self.FILTER_IDENTITY["filter_lambda_min"].to("Angstrom").value
        stray_image[0].header["WAVEMAX"] = self.FILTER_IDENTITY["filter_lambda_max"].to("nm").value

        
        # Add the keys to the main_offender image as well.: 
        main_offender_image[0].header = stray_image[0].header

        # Verify, save and close
        stray_image.verify("silentfix")
        main_offender_image.verify("silentfix")
        stray_image.writeto(self.output_name, overwrite=True)
        main_offender_image.writeto(self.main_offender_output_name, overwrite=True)
        
        # Generate the drizzled and scaled version of the images
        print(datetime.now().isoformat() + " > Drizzling maps... ")
        scaled_drz_names = rs.utils.generate_scaled_drz(stray_flc_name=self.output_name, straylevel_db=straylevel_db,
                                                        mainoff_flc_name=self.main_offender_output_name,
                                                        #input_ext=self.SCIEXTS,
                                                        verbose=verbose)
        self.stray_drz_name = scaled_drz_names["stray_drz_name"]
        self.scaled_stray_drz_name = scaled_drz_names["scaled_stray_drz_name"]
        self.scaled_main_off_name = scaled_drz_names["scaled_main_off_name"]

        # Writing necessary keywords in the output mosaics. 
        keywords = ["RA_TARG", "DEC_TARG", "EXPSTART", "EXPTIME", "FILTER", "WAVEREF", "WAVEMIN", "WAVEMAX", "TELESCOP", "INSTRUME", "DETECTOR"]
        key_values = rs.utils.get_keys_from_header([self.output_name], keywords, ext=0)
        rs.utils.write_parameters_list([self.stray_drz_name], keywords, key_values, ext=0)
        rs.utils.write_parameters_list([self.scaled_stray_drz_name], keywords, key_values, ext=0)
        rs.utils.write_parameters_list([self.scaled_stray_drz_name], ["PIXSCALE"], [[1]], ext=0)
        rs.utils.write_parameters_list([self.scaled_stray_drz_name], ["REBINNED"], [[10]], ext=0)

        rs.utils.write_parameters_list([self.main_offender_output_name], keywords, key_values, ext=0)
        rs.utils.write_parameters_list([self.scaled_main_off_name], keywords, key_values, ext=0)
        rs.utils.write_parameters_list([self.scaled_main_off_name], ["PIXSCALE"], [[1]], ext=0)
        rs.utils.write_parameters_list([self.scaled_main_off_name], ["REBINNED"], [[10]], ext=0)
        print(datetime.now().isoformat() + " > Done")


        ### Generate the straylight report pdf

        print(datetime.now().isoformat() + " > Summary plots... ")
        self.pdf_report_name = rs.plots.make_straylight_plots(RA_TARG=self.RA_TARG, 
                                       DEC_TARG=self.DEC_TARG, 
                                       PA=self.PA, 
                                       source_catalog=self.source_catalog, 
                                       ASTROPYWCS=self.ASTROPYWCS, 
                                       stray_flc_name=self.output_name, 
                                       scaled_stray_drz_name=self.scaled_stray_drz_name, 
                                       scaled_main_off_name=self.scaled_main_off_name, 
                                       figsize=(10,7), mu_vmin = 25, 
                                       mu_vmax = 35, verbose=1)   
        print(datetime.now().isoformat() + " > Done")

        print("Output saved in: " + self.output_name)
        print("Report saved in: " + self.pdf_report_name)

        return({"stray_flc_name": self.output_name,
                "mainoff_flc_name": self.scaled_main_off_name})
    

