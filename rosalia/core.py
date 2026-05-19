# Telescope definition object 
# telescope = {"pointing": [RA_TARG, DEC_TARG], "PA_Y": PA_Y, "detector_shape": [NAXIS1, NAXIS2], "pixscale": pixscale,
#              "R_mirror": R_mirror, "OBSLOC": OBSLOC, "EXPSTART": EXPSTART, "EXPTIME": EXPTIME}

from astropy.time import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import rosalia as rs
import ephessos as ep

class exposure():
    
    def __init__(self, filename=None, telescope=None):

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

            from astropy.wcs import WCS
            self.ASTROPYWCS = [WCS(header)]
            self.FPA_NEAR_RADIUS = self.get_max_angular_size()
            import os 
            self.FILENAME = os.getcwd() + "/" + self.TELESCOP + "_RA_" + '{:07.3f}'.format(self.RA_TARG) +\
                                     "_DEC_" + '{:07.3f}'.format(self.DEC_TARG) +\
                                     "_MJD_" + '{:07.5f}'.format(self.EXPSTART) +\
                                     "_PA_" + '{:06.2f}'.format(self.PA) + ".fits"

        if filename is not None:
            exposure_identity = rs.utils.exposure_inspector(filename, lite=False)
        
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
            self.HST_TYPE = exposure_identity['HST_TYPE']
            self.FILTER = exposure_identity['FILTER']
            self.FILTER_IDENTITY = exposure_identity['FILTER_IDENTITY']
            # self.PHYSPIX = exposure_identity['PHYSPIX']
            self.PIXSCALE = exposure_identity['PIXSCALE']
            self.SCIEXTS = exposure_identity['SCIEXTS']
            self.DATA_SHAPE = exposure_identity['DATA_SHAPE']
            self.ASTROPYWCS = exposure_identity['ASTROPYWCS']
            self.FILETYPE = exposure_identity['FILETYPE']
            self.MPC_OBSLOC = self.get_mpc_observer_location()
            self.JPL_OBSLOC = self.get_jpl_observer_location()
            # self.XYZ_HELIO_POS = exposure_identity['XYZ_HELIO_POS']
            self.FPA_NEAR_RADIUS = self.get_max_angular_size()


    def get_jpl_observer_location(self):
        obs_center = ep.core.jpl_name_translation(self.TELESCOP)
        return(obs_center)


    def get_mpc_observer_location(self):
        from astroquery.mpc import MPC

        if self.TELESCOP == "HST" or self.TELESCOP == "Hubble": return("250")
        if self.TELESCOP == "RST" or self.TELESCOP == "Roman": return("289")
        if self.TELESCOP == "Euclid": return("273")

        try:
            obs = MPC.get_observatory_codes()
            return(obs[obs["Name"] == self.TELESCOP]["Code"][0])
        except:
            print("MPC Code not found!")
        


    def get_detector_corners(self):
        detector_corners = []
        for i in range(len(self.SCIEXTS)):
            detector_corners.append(rs.detectors.get_detector_corners(self.ASTROPYWCS[i]))
        return(detector_corners)
    

    def get_max_angular_size(self):
        return(rs.utils.find_max_angular_size_of_image(wcs=self.ASTROPYWCS, ra_cen=self.RA_TARG, dec_cen=self.DEC_TARG))
    

    def plot_footprint(self, verbose=True, ax=None, color='red', label=None):
        print("Hey! plot_footprint")
        if verbose: print("Finding ra dec constraints")
        # ra_dec_constraints = rs.gaia.find_ra_dec_constraints(self.RA_TARG, self.DEC_TARG, radius=self.FPA_NEAR_RADIUS, verbose=verbose)
        if verbose: print("Getting detector corners")
        exp_corners = self.get_detector_corners()
        if verbose: print("Plotting...")
        if ax is None:
            fig, ax = plt.subplots(figsize=(8,8))
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
        

    def get_nearby_ssos(self, ra=None, dec=None, radius=None, mjd=None, verbose=False, time_step="30s"):
        import ephessos as ep
        if ra==None: ra = self.RA_TARG
        if dec==None: dec = self.DEC_TARG
        if radius==None: radius = self.FPA_NEAR_RADIUS*60*60
        if mjd==None: mjd = self.EXPSTART

        cone_search = ep.core.cone_search(ra=ra, dec=dec, mjd=mjd, 
                                          search_radius=radius, observatory=self.MPC_OBSLOC, verbose=verbose)
        ephessos_df = ep.core.ephessos(sso_search=cone_search, mjd_start=self.EXPSTART, mjd_end=self.EXPEND, 
                                       obs_center=self.JPL_OBSLOC, step_size=time_step, verbose=verbose)
        return(cone_search, ephessos_df)
    

    def make_dummy(self, binning=1):
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

    
