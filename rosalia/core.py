# Telescope definition object 
# telescope = {"pointing": [RA_TARG, DEC_TARG], "PA_Y": PA_Y, "detector_shape": [NAXIS1, NAXIS2], "pixscale": pixscale,
#              "R_mirror": R_mirror, "OBSLOC": OBSLOC, "EXPSTART": EXPSTART, "EXPTIME": EXPTIME}

import os
import numpy as np
import pandas as pd
import rosalia as rs
from datetime import datetime
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord

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
                                                                        detector=telescope["FILTER_PARAMS"]["DETECTOR"], verbose=verbose)
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
            exposure_identity = rs.inspector.exposure_inspector(filename, lite=False)
            
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
            self.conversion_megajanskys = exposure_identity['conversion_megajanskys']
            self.conversion_megajanskys_uncertainty = exposure_identity['conversion_megajanskys_uncertainty']
            self.pixel_area = exposure_identity['pixel_area']
            self.EXPSTART_ISOT = exposure_identity['EXPSTART_ISOT']
            self.PA = exposure_identity['PA']
            self.SCA = exposure_identity['SCA']
            # self.HST_TYPE = exposure_identity['HST_TYPE']
            self.FILTER = exposure_identity['FILTER']
            self.FILTER_IDENTITY = exposure_identity['FILTER_IDENTITY']
            # self.PHYSPIX = exposure_identity['PHYSPIX']
            self.PIXSCALE = exposure_identity['PIXSCALE']
            self.SCIEXTS = exposure_identity['SCIEXTS']
            self.DATA_SHAPE = exposure_identity['DATA_SHAPE']
            self.ASTROPYWCS = exposure_identity['ASTROPYWCS']
            self.FILETYPE = exposure_identity['FILETYPE']

            # Simplify the attributes that are supposed to be constant across the exposure. 
            if len(list(set(self.TELESCOP))) == 1: self.TELESCOP = self.TELESCOP[0]
            if len(list(set(self.INSTRUME))) == 1: self.INSTRUME = self.INSTRUME[0]
            if len(list(set(self.DETECTOR))) == 1: self.DETECTOR = self.DETECTOR[0]
            if len(list(set(self.RA_TARG))) == 1: self.RA_TARG = self.RA_TARG[0]
            if len(list(set(self.DEC_TARG))) == 1: self.DEC_TARG = self.DEC_TARG[0]
            if len(list(set(self.FILTER))) == 1: self.FILTER = self.FILTER[0]
            # Check if all dictionaries in the list are identical copies
            are_all_FILTERS_same = not self.FILTER_IDENTITY or rs.inspector.are_all_dicts_equal(self.FILTER_IDENTITY)
            if are_all_FILTERS_same: self.FILTER_IDENTITY = self.FILTER_IDENTITY[0] 
            #print(self.FILTER_IDENTITY[0])
            if len(list(set(self.EXPSTART_ISOT))) == 1: self.EXPSTART_ISOT = self.EXPSTART_ISOT[0] 
            if len(list(set(self.EXPSTART))) == 1: self.EXPSTART = self.EXPSTART[0] 
            if len(list(set(self.EXPTIME))) == 1: self.EXPTIME = self.EXPTIME[0] 
            if len(list(set(self.EXPMID))) == 1: self.EXPMID = self.EXPMID[0] 
            if len(list(set(self.EXPEND))) == 1: self.EXPEND = self.EXPEND[0] 


            print("TEMP WARNING: RA_TARG is set to element 0 - This has to be fixed on MAST")
            print("TEMP WARNING: DEC_TARG is set to element 0 - This has to be fixed on MAST")
            print("TEMP WARNING: PA is set to element 0 - This has to be fixed on MAST")
            self.RA_TARG = self.RA_TARG[0]
            self.DEC_TARG = self.DEC_TARG[0]
            self.PA = self.PA[0]
            # --------------------------------------------------------- #
            self.ROOTNAME = rs.inspector.longest_common_substring(self.FILENAME)
            self.MPC_OBSLOC = rs.horizons.get_mpc_observer_name(self.TELESCOP)
            self.JPL_OBSLOC = rs.horizons.get_jpl_observer_name(self.TELESCOP)
            state_vectors = rs.horizons.get_heliocoords(self.EXPSTART, body=self.TELESCOP) # rs.horizons.interpolate_Roman_state_vectors(self.EXPSTART)
            self.XYZ_HELIO_POS = [state_vectors["x"].to("AU").value[0], state_vectors["y"].to("AU").value[0], state_vectors["z"].to("AU").value[0]]*u.AU # in AU.
            self.FPA_NEAR_RADIUS = 0.6 # Temporary fix until 'target.ra' 'target.dec' are fixed -> Switch to self.get_max_angular_size() when done (Borlaff - Sept 23, 2026)


        self.fits_keywords = {"TELESCOP": self.TELESCOP, "INSTRUME": self.INSTRUME, "DETECTOR": self.DETECTOR, 
            "FILTER": self.FILTER_IDENTITY["wavelength"], 
            "RA_TARG": self.RA_TARG, "DEC_TARG": self.DEC_TARG, 
            "PA": self.PA, "EXPTIME": self.EXPTIME, 
            "EXPSTART": self.EXPSTART, 
            "EXPSTART_ISOT": self.EXPSTART_ISOT, 
            "WAVEREF": self.FILTER_IDENTITY["filter_lambda_ref"].to("nm").value, 
            "WAVEMIN": self.FILTER_IDENTITY["filter_lambda_min"].to("Angstrom").value, 
            "WAVEMAX": self.FILTER_IDENTITY["filter_lambda_max"].to("nm").value}

    def roman_wfi_exposure(self, observer, prefix=""):
        print("> Synthetic Roman_wfi_exposure")
        # Here we expect observer={"TELESCOP": "Roman/WFI", "pointing": [RA_TARG, DEC_TARG], "FILTER":FILTER, "PA_Y": PA_Y, "EXPSTART": EXPSTART, "EXPTIME": EXPTIME}
        # Fixed parameters for Roman/WFI 
        observer['TELESCOP'] = "Roman"
        observer['INSTRUME'] = "WFI"
        observer['DETECTOR'] = "WFI"
        observer['DATA_SHAPE'] = [4088, 4088]
        observer['DATA'] = 18*[np.zeros(observer['DATA_SHAPE'])]
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
                                                                detector=observer["FILTER_PARAMS"]["DETECTOR"], verbose=verbose)
        state_vectors = rs.horizons.get_heliocoords(self.EXPSTART, body=self.TELESCOP) # rs.horizons.interpolate_Roman_state_vectors(self.EXPSTART)
        self.XYZ_HELIO_POS = [state_vectors["x"].to("AU").value[0], state_vectors["y"].to("AU").value[0], state_vectors["z"].to("AU").value[0]]*u.AU # in AU. 

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
        self.FPA_NEAR_RADIUS = 0.6 # Temporary fix until 'target.ra' 'target.dec' are fixed -> Switch to self.get_max_angular_size() when done (Borlaff - Sept 23, 2026)

        return(self)
    
    def get_detector_corners(self):
        detector_corners = []
        for i in range(len(self.SCIEXTS)):
            detector_corners.append(rs.detectors.get_detector_corners(self.ASTROPYWCS[i]))
        return(detector_corners)
    

    def get_max_angular_size(self):
        # print(self.ASTROPYWCS)
        return(rs.utils.find_max_angular_size_of_image(wcs=self.ASTROPYWCS, ra_cen=self.RA_TARG, dec_cen=self.DEC_TARG))
    

    def plot_footprint(self, figsize=(10,10), verbose=True, ax=None, color='red', label=None):
        # print("Hey! plot_footprint")
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

        self.source_catalog_filename = os.path.basename(self.ROOTNAME) + "_source_catalog.csv" #
        
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


    def get_nearby_ssos(self, ra=None, dec=None, radius=None, mjd=None, verbose=False, time_step="1m"):
        import ephessos as ep
        if ra==None: ra = self.RA_TARG
        if dec==None: dec = self.DEC_TARG
        if radius==None: radius = self.FPA_NEAR_RADIUS*60*60
        if mjd==None: mjd = self.EXPSTART

        cone_search = ep.core.cone_search(ra=ra, dec=dec, mjd=mjd, 
                                          search_radius=radius, observatory=self.MPC_OBSLOC, verbose=verbose)

        if verbose: print(cone_search)
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


    # def make_mosaic(self, outname=None):

    

    def HST_straylight(self):

        source_catalog = self.get_source_catalog()

        target_name = self.FILENAME

        straylight_images = []
        for i in range(len(self.SCIEXTS)):
            ra_detector, dec_detector, straylight_lvl  = rs.hst.hst_estimate_straylight_SCA(
                                                                wcs=self.ASTROPYWCS[i],
                                                                filter=self.FILTER_IDENTITY,
                                                                catalog=source_catalog,
                                                                instrument="ACS",
                                                                verbose=False)

            # Now interpolate the skypoints to the original image

            target_ext = self.SCIEXTS[i]
            straylight_image = rs.utils.interpolate_skypoints_to_image(target_name, target_ext, 
                                                                    ra=ra_detector, dec=dec_detector, 
                                                                    z=straylight_lvl,
                                                                    mask_original_nan=True,
                                                                    method="linear") 
            straylight_images.append(straylight_image)

        return(straylight_images)

    def make_close_stars_ds9_region(self, radius=5, output="nearby_stars.reg"):
        import numpy as np
        source_coords = SkyCoord(ra=np.array(self.source_catalog["ra"])*u.deg, dec=np.array(self.source_catalog["dec"])*u.deg, frame="icrs")
        pointing_coords = SkyCoord(ra=self.RA_TARG*u.deg, dec=self.DEC_TARG*u.deg, frame="icrs")
        sep = source_coords.separation(pointing_coords)
        close_catalog = self.source_catalog[sep < radius*u.deg]
        rs.utils.make_ds9_region(ra=np.array(close_catalog["ra"]), dec=np.array(close_catalog["dec"]), 
                                 radius=np.array([1]), label=np.array(close_catalog["source_id"]), output=output)
        print(output)


#####################################################################################
### Straylight ######################################################################
#####################################################################################

    def save_flc(self, outname, data_list=None, wcs_list=None, keywords=None, extname=None, overwrite=True):
        if data_list is None: data_list = self.DATA
        if wcs_list is None: wcs_list = self.ASTROPYWCS
        if keywords is None: keywords = self.fits_keywords

        headers_output = []
        for i in range(len(wcs_list)):
            ASTROPYWCS_i = wcs_list[i]
            header = ASTROPYWCS_i.to_header()

            if isinstance(keywords["DETECTOR"], (list,)):
                keywords["DETECTOR"] = rs.inspector.longest_common_substring(self.DETECTOR)

            for key in keywords.keys():
                header[key] = keywords[key]

            header["EXTNAME"] = "SCI"
            header["SCA"] = self.SCA[i]
            header["RA_TARG"] = self.RA_TARG
            header["DEC_TARG"] = self.DEC_TARG
            header["EXPSTART"] = self.EXPSTART
            header["EXPEND"] = self.EXPEND
            header["PA"] = self.PA
            header["TELESCOP"] = self.TELESCOP
            header["INSTRUME"] = self.INSTRUME

            headers_output.append(header)

        rs.utils.save_fits(array=data_list, 
                           name=outname, 
                           header=headers_output,
                           extname=extname, 
                           overwrite=overwrite, 
                           output_verify='silentfix')        

        return(outname)
        
    def save_drz(self, outname,  data_list=None, wcs_list=None, keywords=None, resolution=1, extname=None, overwrite=True):
        if data_list is None: data_list = self.DATA
        if wcs_list is None: wcs_list = self.ASTROPYWCS
        if keywords is None: keywords = self.fits_keywords

        if isinstance(keywords["DETECTOR"], (list,)):
            keywords["DETECTOR"] = rs.inspector.longest_common_substring(self.DETECTOR)
            
        headers_output = []
        for i in range(len(wcs_list)):
            ASTROPYWCS_i = wcs_list[i]
            header = ASTROPYWCS_i.to_header()

            header["EXTNAME"] = "DRZ"
            header["RA_TARG"] = self.RA_TARG
            header["DEC_TARG"] = self.DEC_TARG
            header["EXPSTART"] = self.EXPSTART
            header["EXPEND"] = self.EXPEND
            header["PA"] = self.PA
            header["TELESCOP"] = self.TELESCOP
            header["INSTRUME"] = self.INSTRUME
            header["BUNIT"] = "MJy sr-1"

            # Convert the output units to MJy/sr
            data_list[i] = data_list[i]*self.conversion_megajanskys[i]
            keywords["DETECTOR"] = rs.inspector.longest_common_substring(self.DETECTOR)
            for key in keywords.keys():
                header[key] = keywords[key]

            headers_output.append(header)

        drz_data, drz_wcs = rs.utils.generate_mosaic(data=data_list, astropywcs=headers_output, resolution=resolution)



        rs.utils.save_fits(array=np.float32(drz_data), name=outname, header=drz_wcs,
                           extname=extname, 
                           overwrite=overwrite, 
                           output_verify='silentfix')
        return(drz_data, drz_wcs, outname)



    def get_mosaic_wcs(self):
        from reproject.mosaicking import find_optimal_celestial_wcs
        resolution = self.PIXSCALE*u.arcsec
        input_data_for_reproject = list(zip(self.DATA, self.ASTROPYWCS))

        optimal_wcs = find_optimal_celestial_wcs(input_data=input_data_for_reproject, 
                                                resolution=resolution)

        return(optimal_wcs[0], optimal_wcs[1])



    from concurrent.futures import ProcessPoolExecutor
    from tqdm import tqdm

    @staticmethod
    def _parallel_roman_estimate_straylight_SCA(args):
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


    def straylight(self, catalog=None, g_mag_max=15, sun_block=False, resolution=1, verbose=False):
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


        if self.TELESCOP.lower() == "hst":
            straylight_images= self.HST_straylight()
            mos_data, mos_wcs = rs.utils.generate_mosaic(data=straylight_images, astropywcs=self.ASTROPYWCS, resolution=None)
            self.output_name = self.FILENAME.replace(".fits", "_stray.fits")
            rs.utils.save_fits(mos_data, self.output_name, mos_wcs)
            return(mos_data, mos_wcs)

        elif self.TELESCOP.lower() == "roman" or self.TELESCOP.lower() == "rst":
            # Find out which stars belong to each detector. 
            print(datetime.now().isoformat() + " > Fetching source catalog")

            if not hasattr(self, 'source_catalog'):
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
            elif "asdf" in self.FILENAME:
                self.output_name = self.ROOTNAME + "_stray.fits"
            else: 
                self.output_name = self.ROOTNAME + "_stray.fits"

            self.main_offender_output_name = self.output_name.replace(".fits", "_main_off.fits")
            self.straylight_db_output_name = self.output_name.replace(".fits", "_db.csv")

            if verbose: print(self.main_offender_output_name)
            if verbose: print(self.straylight_db_output_name)


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
                results = list(tqdm(executor.map(self._parallel_roman_estimate_straylight_SCA, inputs),total=len(inputs),))
                straylevel_all_SCAS.extend(results)

            print(datetime.now().isoformat() + " > Done : " + str(datetime.now() - t) + " elapsed.")

            ###########
            # Reconstruct the stray-light maps and main offender maps from the straylevel_all_SCAS database. 
            ###########
            NSCAs = len(self.SCIEXTS)# NSCAs 
            print(" > Reconstructing the Stray-light / Main offender map: ")
            straylevel_list = [] 
            main_offender_list = [] 
            for SCA in range(NSCAs):
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
            # Save the stray-light FLC file
            self.save_flc(outname=self.output_name, data_list=straylevel_list, wcs_list=self.ASTROPYWCS, keywords=self.fits_keywords)

            # Save the stray-light DRZ file 
            self.stray_drz_name = self.output_name.replace(".fits","_drz.fits")
            self.save_drz(outname=self.stray_drz_name, data_list=straylevel_list, wcs_list=self.ASTROPYWCS, keywords=self.fits_keywords, resolution=resolution)

            # Save the main-offender FLC file
            self.save_flc(outname=self.main_offender_output_name, data_list=main_offender_list, wcs_list=self.ASTROPYWCS, keywords=self.fits_keywords)

            # Save the main-offender DRZ file 
            self.mainoff_drz_name = self.main_offender_output_name.replace(".fits","_drz.fits")
            self.save_drz(outname=self.mainoff_drz_name, data_list=main_offender_list, wcs_list=self.ASTROPYWCS, keywords=self.fits_keywords, resolution=resolution)

            # --------------
            # We need to apply one extra correction to the Main_offender DRZ file
            # so no drizzling artifacts are included in the plots.
            # -------------- 
            main_off_map = fits.open(self.mainoff_drz_name)
            mask_main_off_ok = np.zeros(main_off_map[0].data.shape)
            main_offenders = list(set(straylevel_db["mainoffender_total"]))

            # Remove all those pixels that have values not included in the list of main offenders 
            for i in range(len(main_offenders)):
                main_off_id = main_offenders[i]
                main_offended_pixels = np.where(main_off_map[0].data == main_off_id)
                mask_main_off_ok[main_offended_pixels] = 1

            main_off_map[0].data[mask_main_off_ok == 0] = np.nan
            main_off_map.verify("silentfix")
            main_off_map.writeto(self.mainoff_drz_name, overwrite=True)
    

            ################################################################
            ############ Generate the straylight report pdf ################
            ################################################################

            print(datetime.now().isoformat() + " > Summary plots... ")

            self.pdf_report_name = rs.plots.make_straylight_plots(RA_TARG=self.RA_TARG, 
                                        DEC_TARG=self.DEC_TARG, 
                                        PA=self.PA, 
                                        source_catalog=self.source_catalog, 
                                        ASTROPYWCS=self.ASTROPYWCS, 
                                        stray_flc_name=self.output_name, 
                                        scaled_stray_drz_name=self.stray_drz_name, 
                                        scaled_main_off_name=self.mainoff_drz_name, 
                                        figsize=(10,7), mu_vmin = 25, 
                                        mu_vmax = 35, verbose=1) 
              
            print(datetime.now().isoformat() + " > Done")

            print("Output saved in: " + self.output_name)
            print("Report saved in: " + self.pdf_report_name)

            return({"straylevel_db": straylevel_db, 
                    "stray_flc_name": self.output_name,
                    "mainoff_flc_name": self.main_offender_output_name,
                    "stray_drz_name": self.stray_drz_name,
                    "mainoff_drz_name": self.mainoff_drz_name})

        else:
            print("Straylight modeling is currently only available for Roman/WFI exposures.")
            return(None)

#####################################################################################
### PSF         ########################################################################
#####################################################################################


    def psf(self, g_mag_max=15, catalog=None, verbose=False):
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
        from astropy.io import fits
        import logging
        logger = logging.getLogger()
        logger.setLevel(logging.CRITICAL)

        # Get the catalog of the stars around the FOV
        if not hasattr(self, 'source_catalog'):
            if catalog is None:
                self.source_catalog = self.get_source_catalog(g_mag_max=g_mag_max, verbose=False)
            else:
                self.source_catalog = rs.utils.fix_custom_catalog(catalog)


        # Generate the star stamps (PSFs)
        print("TO DO: Make stamps with a more reasonable size. Dim stars can have smaller PSFs.")
        print("To do this, make a profile of the Roman / PSF, and find out when would it be essentially 0.")
        star_stamps = rs.psf.generate_star_stamps(hybrid_catalog=self.source_catalog,
                                                # image_identity=image_identity),
                                                 telescope=self.TELESCOP, 
                                                 filename=self.FILENAME, 
                                                 sciexts=self.SCIEXTS, 
                                                 astropywcs=self.ASTROPYWCS,
                                                 filter=self.FILTER_IDENTITY["wavelength"],
                                                 pa=self.PA, verbose=verbose)
        # def generate_star_stamps(hybrid_catalog, telescope, filename, sciexts, astropywcs, filter, pa, verbose=False):


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
        star_swarp_name = self.FILENAME.replace(".fits", "_stars_drz.fits")
        rs.utils.execute_cmd("mv coadd.fits " + star_swarp_name) # Make a compressed version, for easiest visualization.

        # Now blot back to the dummy SCA per SCA frame
        from reproject import reproject_interp
        # Let's make a dummy copy to reproject the stars into
        roman_dummy = fits.open(self.FILENAME, memmap=True)
        star_model = fits.open(star_swarp_name, memmap=True)

        if verbose: print("Storing stars in each SCA")
        for SCIEXT_i in tqdm(self.SCIEXTS):
            # Open the star fits
            star_reprojected, footprint = reproject_interp(star_model[0],
                                                           roman_dummy[SCIEXT_i].header, 
                                                           parallel=True)
            star_reprojected[np.isnan(star_reprojected)] = 0
            roman_dummy[SCIEXT_i].data = roman_dummy[SCIEXT_i].data + star_reprojected

        roman_dummy.verify("silentfix")
        star_output_name = self.FILENAME.replace(".fits", "_stars.fits")
        roman_dummy.writeto(star_output_name, overwrite=True)

        # Now make again the drz, this time with the correct gaps.
        drz_name, scaled_drz_name = rs.utils.run_swarp(pattern=star_output_name, 
                                                       outname=star_output_name.replace(".fits","_drz.fits"), scale=0.11)


        if verbose: 
            print("In-field stray-light model completed. Level 2 multi-extension FITS: " + star_output_name)
            print("Mosaic image: " + drz_name)
            print("Scaled mosaic: " + scaled_drz_name)
        return(star_output_name)

#####################################################################################
### Zodiacal ########################################################################
#####################################################################################

    def zodiacal(self, zody_mode="zodipy", verbose=False, output_name=None, output_units="e/s", resolution=1):
        from tqdm import tqdm
        import logging
        logger = logging.getLogger()
        logger.setLevel(logging.CRITICAL)

        nSCIEXTS = len(self.SCIEXTS)
        # Make the Roman Dummy image
        roman_dummy_name = self.FILENAME 
        if output_name is None:
            output_name = self.FILENAME.replace(".fits", "_zody.fits").replace(".asdf", "_zody.fits")

        zodiacal_background_list = []
        zodiacal_background_unit_list = []

        print("Computing Zodiacal light...")
        for i in tqdm(range(nSCIEXTS)):
            zodiacal_background = rs.sky.get_zodiacal_background(astropywcs=self.ASTROPYWCS[i],
                                                                 wavelength=self.FILTER_IDENTITY,
                                                                 expstart=self.EXPSTART,
                                                                 step=1000, zody_mode=zody_mode,
                                                                 nbins_wavelength=10, obslocin=0,
                                                                 grid_method="random",
                                                                 obspos = self.XYZ_HELIO_POS,
                                                                 output_units=output_units,
                                                                 verbose=verbose)

            zodiacal_background_list.append(zodiacal_background.value)
            zodiacal_background_unit_list.append(zodiacal_background.unit.to_string())

        ########################################
        # Save the Zodiacal light FLC file
        ########################################
        self.drz_zody_name = output_name.replace(".fits","_drz.fits")
        self.save_flc(outname=output_name, data_list=zodiacal_background_list, wcs_list=self.ASTROPYWCS)

        # Save the Zodiacal light DRZ file 
        self.save_drz(outname=self.drz_zody_name, data_list=zodiacal_background_list, wcs_list=self.ASTROPYWCS, resolution=resolution)

        rs.plots.make_stray_plot(input_name=self.drz_zody_name, ext=0, mode="fe", 
                                 color_label = "Flux (e/s/px)", cmap="RdYlBu_r", figsize=(10,7))

        rs.plots.make_stray_plot(input_name=self.drz_zody_name, ext=0, mode="fe2mu", 
                                 vmin=None, vmax=None, 
                                 color_label = 'Surface brightness (mag arcsec$^{-2}$)',
                                 cmap="RdYlBu", output_name=None, figsize=(10,7), mu_vmin=None, mu_vmax=None)

        print("Output saved in: " + output_name)

        return({"zodi_list": zodiacal_background_list,
                "output_name": output_name})



