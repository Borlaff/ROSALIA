# Alejandro S. Borlaff. NASA Ames Research Center. a.s.borlaff@nasa.gov / asborlaff@gmail.com
# January 20, 2023.
#
# STRAYCOR/SKY module
# This module will hold all the general programs related to the sky background
# not related to other tasks, or more general than other wrappers.
#
#
# Version log:
# v.1.0 - 20 Enero 2023. First loading of programs inherited from former monolithic straycor.py
# v.1.1 August 2026 integration with StSci zodiacal light calculator by SLN
##########################################################
import os
import glob
import numpy as np
import pandas as pd
from tqdm import tqdm
import bottleneck as bn
import astropy.units as u
from astropy.io import fits
from astropy.io import ascii
from astropy.time import Time
from astropy.coordinates import SkyCoord
from scipy import constants
import matplotlib.pyplot as plt
# import gunagala.sky as skies # Gunagala is not maintained. It requires to detach from astropy_helpers. https://github.com/astropy/astropy-helpers
import io
import copy
import healpy
from scipy.interpolate import interp1d
from astropy.coordinates import get_body
import warnings
from astropy.utils.exceptions import AstropyWarning
import urllib.request

import rosalia as rs

def remove_zodiacal_light_acs(input_name, zody_mode="stsci", verbose=False):
    zody_interp_sci1 = get_zodiacal_background(input_name=input_name, ext=1, zody_mode=zody_mode, verbose=verbose)
    zody_interp_sci2 = get_zodiacal_background(input_name=input_name, ext=4, zody_mode=zody_mode, verbose=verbose)

    input_fits = fits.open(input_name)
    input_fits[1].data = input_fits[1].data - zody_interp_sci1
    input_fits[4].data = input_fits[4].data - zody_interp_sci2

    input_fits.verify("silentfix")
    zody_cor_name = input_name.replace(".fits","_zodycor.fits")
    input_fits.writeto(zody_cor_name, overwrite=True)

    input_fits[1].data = zody_interp_sci1
    input_fits[4].data = zody_interp_sci2

    input_fits.verify("silentfix")
    zody_name = input_name.replace(".fits","_zody.fits")
    input_fits.writeto(zody_name, overwrite=True)
    return({"input": input_name, "zody": zody_name, "zody_cor": zody_cor_name})


def correct_flat_sky(input_name, ext, overwrite=True, clean=True, verbose=False):
    # So, lets write a program that does a flat sky background correction for us.
    # Run noisechisel with the default parameters
    stdout = rs.utils.execute_cmd("astnoisechisel --outliernumngb=5 -K -h" + str(ext) + " " + input_name)
    detected_name = input_name.replace(".fits", "_detected.fits")

    if not os.path.exists(detected_name):
        return({"skylvl": np.nan, "input_name": input_name, "ext": ext})

    # Open the detected fits file
    input_fits = fits.open(input_name)
    detected_fits = fits.open(detected_name)
    sky_sample = input_fits[ext].data[np.where(detected_fits["DETECTIONS"].data == 0)]
    sky_level = bn.nanmedian(sky_sample)

    if verbose:
        print("Flat sky level:" + str(sky_level))
        print("Subtracting from " + input_name + "[" + str(ext) + "]")
        print("Storing in header KEYWORD SKYLVL1")


    if overwrite:
        try:
            input_fits[ext].header["SKYLVL1"] = input_fits[ext].header["SKYLVL1"] + sky_level
        except:
            input_fits[ext].header["SKYLVL1"] = sky_level
        input_fits[ext].data = input_fits[ext].data - sky_level
        input_fits.verify("silentfix")
        input_fits.writeto(input_name, overwrite=True)

    input_fits.close()
    detected_fits.close()

    if clean:
        rs.utils.execute_cmd("rm " + detected_name)

    return({"skylvl": sky_level, "input_name": input_name, "ext": ext})


#####################################################

def rebin( a, newshape ):
        '''Rebin an array to a new shape.
        '''
        assert len(a.shape) == len(newshape)
        slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
        coordinates = np.mgrid[slices]
        indices = coordinates.astype('i')   #choose the biggest smaller integer index
        return a[tuple(indices)]



def rebin_transmission_curve(filter_transmission_curve, nbins, verbose=False):
    from scipy import interpolate

    # This program takes as input a table with two columns
    # 1) lambda_AA: Wavelength in Angstroms
    # 2) Transmission: Transmission (fraction)
    #
    # and a number of bins as nbins_wavelength
    #
    # Output:
    # A dictionary with:
    # wavelength: Original wavelength points
    # transmission: Original transmission points
    # rebinned_wavelength: Rebinned wavelength
    # rebinned_transmission: Rebinned transmission
    # rebinned_dlambda: Rebinned delta wavelength (spacing between each bin)

    # First we read the file with the curve
    # We might be able to do this automatically online using:
    # Response calculated using the stsynphot python package by STScI.
    # In particular, function stsynphot.band("acs,wfc1,f814w").
    # Please, take a look to stsynphot documentation and HST instrument documentation

    #transmission_curve = ascii.read(filter_curve_name)

    # We assume the first column contains lambda, and the second the transmission
    wavelength = filter_transmission_curve["Wavelength"]
    transmission = filter_transmission_curve["Transmission"]

    # We generate an interpolator function with the datapoints from the table
    good = np.where(np.isfinite(wavelength) & np.isfinite(transmission))

    if verbose:
        print("wavelength[good]")
        print(wavelength[good])
        print("transmission[good]")
        print(transmission[good])

    filter_interpolator = interpolate.interp1d(wavelength[good], transmission[good])

    # We rebin to a new grid, defined by nbins
    rebinned_wavelength = rebin(wavelength, (nbins,))

    rebinned_transmission = filter_interpolator(np.array(rebinned_wavelength))

    # Now we calculate the bin size in lambda
    dlambda = (np.max(rebinned_wavelength) - np.min(rebinned_wavelength))/(nbins+1)

    if verbose:
        # We plot the interpolated bins
        plt.plot(wavelength,transmission)
        plt.ylabel("Transmission")
        plt.xlabel(r'$\lambda$ ($\AA$)')
        plt.scatter(rebinned_wavelength, rebinned_transmission, marker="o", s=20, color="black")

        # And the dlambda sizes
        plt.plot(0, 0, linewidth=2, color="black", label = "Pseudo-spectra for background model")
        for i in range(len(rebinned_wavelength)):
             plt.plot(np.array([rebinned_wavelength[i]-dlambda/2., rebinned_wavelength[i]+dlambda/2.]),
             np.array([rebinned_transmission[i], rebinned_transmission[i]]), linewidth=2,
             color="black")
        plt.xlim(0.99*np.min(wavelength), 1.01*np.max(wavelength))
        plt.legend(frameon=False)
        plt.show()

    return({"wavelength": wavelength, "transmission": transmission,
            "rebinned_wavelength": rebinned_wavelength, "rebinned_transmission": rebinned_transmission,
            "rebinned_dlambda": dlambda})



#accepted zody_mode: zodipy and stsci; if something else is entered, defaults to stsci
def get_zodiacal_background(input_name=None, ext=None, exposure_identity=None, wavelength=None, telescope=None, instrument=None,
                            detector=None, expstart=None, step=1000, zody_mode="stsci",
                            nbins_wavelength=20, obslocin=3, grid_method="random",
                            output_units=None, verbose=False, interpolate=True):

    from scipy import interpolate
    import astropy.wcs as astropy_wcs

    # Keywords for IRSA background query:
    # Check: https://irsa.ipac.caltech.edu/applications/BackgroundModel/docs/dustProgramInterface.html
    #
    #
    # locstr 	(locstr examples) 	NA 	If the input is a recognizable object name, it will be resolved into coordinates using NED or, if that fails, SIMBAD (required).
    # wavelength 	float 	0.5-1000.0 	Wavelength in microns (defaults to 2.0).
    # year 	char 	varies 	Year. Limited to 2018 to 2029 for L2 position. Defaults to 2019.
    # day 	char 	1-366 	Day. Limited to 2018 Day 274 to 2029 Day 120 for L2 position and ido_viewin=0. Defaults to 180.
    # obslocin 	char 	0 or 3 	Observing location. 0 is Earth-Sun L2 region; 3 is Earth (defaults to 0).
    # obsverin 	char 	1 or 4 	Code version (optional, defaults to 4).
    # ido_viewin 	char 	0 or 1 	0 = find zodiacal on Day; 1 = find median zodiacal over a likely viewing range (defaults to 1, see Help).

    if exposure_identity is None:
        input_fits = fits.open(input_name)
        astropywcs = astropy_wcs.WCS(header=input_fits[ext].header, fobj=input_fits, naxis=2)

        exposure_identity = rs.utils.exposure_inspector(input_name, lite=True)

    # print(exposure_identity)
    #if output_units == None:
    #    output_units = exposure_identity["BUNIT"]

    # Step zero: Check the wavelength variable.
    # If it is string, it might be a filter name.
    # Look for it in the library FILTERS
    # Check that ra and dec arguments are in array
    if isinstance(wavelength, (str)):
        filter_curve_name = rs.telescopes.find_filter_in_svo(wavelength=exposure_identity["FILTER"],
                                                             telescope=exposure_identity["TELESCOP"], 
                                                             instrument=exposure_identity["INSTRUME"], 
                                                             detector=exposure_identity["DETECTOR"], 
                                                             verbose=verbose)
        
        rebinned_filter_curve = rebin_transmission_curve(filter_transmission_curve=filter_curve_name["filter_transmission_curve"],
                                                         nbins=nbins_wavelength, verbose=verbose)
        rebinned_transmission = rebinned_filter_curve["rebinned_transmission"]

        # Stored filter curves are in Angstrom
        rebinned_wavelength   = rebinned_filter_curve["rebinned_wavelength"]
        dlambda = rebinned_filter_curve["rebinned_dlambda"]

    # If the input is just a wavelength, then emulate the output of rebin_transmission_curve
    if isinstance(exposure_identity["FILTER"], (float)):
        if verbose:
            print("Input wavelength " + str(wavelength))
        rebinned_transmission = np.array([1])
        dlambda = np.array([1])
        rebinned_wavelength = np.array([wavelength])
    ###############################

    # We calculate the expstart
    t = Time(exposure_identity["EXPSTART"], format='mjd', scale='utc')
    year = t.yday.split(":")[0]
    day = t.yday.split(":")[1]


    ########################################

    # First we get the detector grid for interpolating the zodiacal measurements
    detector_grid = rs.detectors.make_detector_grid(w=astropywcs, step=step, mode=grid_method)
    npoints_grid = len(detector_grid["grid_world"][0])
    # Then we query IRSA to get the Zody surface brightness at those positions and time
    if verbose:
        print("Year: " + str(year))
        print("Day: " + str(day))
        print("RA: " + str(np.median(detector_grid["grid_world"][0])))
        print("DEC: " + str(np.median(detector_grid["grid_world"][1])))

    if zody_mode == "zodipy":
        if verbose:
            print("Estimating zodiacal light with Zodipy...")

        # Zodipy queries must go in um so we need to multiply by 1E+7 the m above
        obspos = np.array([exposure_identity["XYZ_HELIO_POS"][0][0].value,
                           exposure_identity["XYZ_HELIO_POS"][1][0].value,
                           exposure_identity["XYZ_HELIO_POS"][2][0].value])*u.AU

        if verbose:
            print("Heliocentric position of telescope:")
            print(obspos)
        zody_MJysr = zodipy_zody(ra=detector_grid["grid_world"][0],
                         dec=detector_grid["grid_world"][1],
                         wavelength=rebinned_wavelength.to("um").value,
                         weights=rebinned_transmission,
                         expstart=expstart,
                         obspos=obspos)
        # print('ra',detector_grid["grid_world"][0])
        # print('dec',detector_grid["grid_world"][1])
        print('wavelength',rebinned_wavelength.to("um").value)
        print('weights',rebinned_transmission)
        print('expstart',expstart)
        print('obspos',obspos)
        print(zody_MJysr)
    else:
        if not zody_mode == "stsci":
            print("Invalid zodiacal light mode entered. Using stsci model as default.")
        zody_MJysr = stsci_zody(ra=detector_grid["grid_world"][0],
                                 dec=detector_grid["grid_world"][1],
                                 wavelength=rebinned_wavelength.to("um").value,
                                 weights=rebinned_transmission,
                                 expstart=expstart,verbose=verbose)


    x = detector_grid["grid_xy"][:,0]
    y = detector_grid["grid_xy"][:,1]
    xsize = int(np.max(x))
    ysize = int(np.max(y))
    x_lin = np.linspace(0,ysize,ysize)
    y_lin = np.linspace(0,xsize,xsize)
    xv, yv = np.meshgrid(y_lin, x_lin, indexing='ij')
    points = [(x[i], y[i]) for i in range(len(x))]

    if interpolate:
        zody_interp =  interpolate.griddata(points, zody_MJysr, (xv, yv), method="cubic")
    else:
        zody_interp = zody_MJysr

    zody_interp = zody_interp*(u.MJy * u.steradian**-1) 
    # print(telescope)

    if output_units == "e/s":
        if verbose: print("Output units:" + output_units)
        #  HST_ACS_jy_to_counts(flux_ACS, photflam, photplam):
        if telescope == "Hubble" or telescope == "HST":
            zody_interp = rs.detectors.HST_ACS_jy_to_counts(flux_jy = zody_interp,
                                           photflam = exposure_identity["PHOTFLAM"],
                                           photplam = exposure_identity["PHOTPLAM"])

        if telescope == "Roman" or telescope == "RST":
            from romanisim import bandpass as ris_bandpass
            es_to_MJysr = ris_bandpass.etomjysr(bandpass=wavelength, sca=ext)*u.MJy * u.steradian**-1 * u.s # The factor F such that MJy / sr = F * DN/s
            zody_interp = zody_interp/es_to_MJysr
        if verbose: print("Output units: e/s")

    else:
        if verbose: print("Output units: Jy/arcsec2")

    return(zody_interp.T)

########################################################
"""
def find_filter_curve_file(wavelength):
    straycor_path = os.path.dirname(os.environ["ROSALIACACHE"])
    filters_path = straycor_path + "/FILTERS/"
    filters_list = glob.glob(filters_path + "*")

    filter_match = []
    for filter_name in filters_list:
        if wavelength in filter_name:
            filter_match.append(filter_name)

    if len(filter_match) > 1:
        print("More than one filter was found with " + wavelength + " name")
        print(filter_match)
        print("Check http://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse&gname=HST")
        raise Exception("ERROR: Undefined filter name. Please use a more specific name")

    else:
        filter_match = filter_match[0]
        print("Filter found: " + filter_match)

    return(filter_match)
"""

#########################################

def zodipy_zody(ra, dec, wavelength, weights, expstart, obspos="earth"):
    """
    This program makes use of the Zodipy Zodiacal model to estimate the
    Zodiacal light surface brightness in a given position, at a wavelength,
    day and year.

    Input:
    ra = Right ascension (degrees)
    dec = Declination (degrees)
    wavelength = Wavelength (micron)
    weights = Transmission at the proper wavelength
    expstart = Modified Julian Day of the observation (MJD)

    Output:
    mu_zody = surface brightness in Jy arcsec-2
    """

    from astropy.coordinates import SkyCoord
    import multiprocessing
    import pandas as pd
    import zodipy
    import astropy.units as u
    from astropy.time import Time

    if not isinstance(ra, (list, pd.core.series.Series, np.ndarray)):
        ra = np.array([ra])
    if not isinstance(dec, (list, pd.core.series.Series, np.ndarray)):
        dec = np.array([dec])

    n_pointings = len(ra)
    # If wavelength, year, or day are not in array, copy their values into one
    # as large as ra, dec
    if not isinstance(wavelength, (list, pd.core.series.Series, np.ndarray)):
        wavelength = np.array([wavelength]*n_pointings)

    if not isinstance(weights, (list, pd.core.series.Series, np.ndarray)):
        weights = np.array([weights]*n_pointings)

    # Use Astropy's `SkyCoord` object to specify coordinates
    skycoord = SkyCoord(ra, dec, unit="deg",  frame="icrs")
    # Note that we manually set the obstime attribute
    skycoord.obstime = Time(expstart, format='mjd', scale='utc')


    #print("Zodipy wave:")
    #print(wavelength)
    # Initialize a zodiacal light model at a wavelength/frequency or over a bandpass
    model = zodipy.Model(wavelength*u.micron, weights=weights, extrapolate=True)

    # Evaluate the zodiacal light model
    # Here the solution depends from obspos.
    # The exposure_inspector should provide the
    # heliocentric ecliptic cartesian position of the satellite.
    #
    # This argument accepts both a string representing
    # a body recognized by the solar system ephemeris,
    # or a heliocentric ecliptic cartesian position.
    emission = model.evaluate(skycoord, obspos=obspos) # , nprocesses=multiprocessing.cpu_count()

    return(emission)


#########################################
#implementing the STSci's zodiacal light model into rosalia
def stsci_zody(ra, dec, wavelength, weights, expstart,verbose=False):
    """
    This program makes use of the STSci's Zodiacal model for Roman to estimate the
    Zodiacal light surface brightness in a given position, at a wavelength,
    day and year.

    Input:
    ra = Right ascension (degrees)
    dec = Declination (degrees)
    wavelength = Wavelength (micron)
    weights = Transmission at the proper wavelength
    expstart = Modified Julian Day of the observation (MJD)

    Output:
    emisison = surface brightness in MJy/sr
    """

    from astropy.coordinates import SkyCoord
    import multiprocessing
    import pandas as pd
    import zodipy
    import astropy.units as u
    from astropy.time import Time

    if not isinstance(ra, (list, pd.core.series.Series, np.ndarray)):
        ra = np.array([ra])
    if not isinstance(dec, (list, pd.core.series.Series, np.ndarray)):
        dec = np.array([dec])

    n_pointings = len(ra)
    # If wavelength, year, or day are not in array, copy their values into one
    # as large as ra, dec
    if not isinstance(wavelength, (list, pd.core.series.Series, np.ndarray)):
        wavelength = np.array([wavelength]*n_pointings)

    if not isinstance(weights, (list, pd.core.series.Series, np.ndarray)):
        weights = np.array([weights]*n_pointings)

    #In the background calendar array, January 1 has index of 0
    #December 31 has index of 364 for non-leap years; index of 365 for leap years
    day = float(Time(expstart, format='mjd').yday.split(':')[1]) - 1

    emission = []
    cache_record = []
    emission_record = []

    for r,d in zip(ra,dec):

        #set nside = 128 as in background class
        cache_check = get_healpix(128,r,d)

        #check to see if zodiacal emission has already been
        # calculated for this cache file (i.e. a very close ra and dec has been called before)
        #if new cache file then calculate background form scratch
        #zodiacal light model is the same for every ra and dec for a given cache file
        #the available calendar days does explicitly rely on the exact ra and dec
        # but assumption is that for close enough ra and dec that use the same cache file
        # will not impact the calendar days
        if cache_check not in cache_record:

            #Call the background class lifted from StSci's rbt.py module
            bkg = background(r, d, verbose=verbose)

            #return the wavelength array and matching zodiacal light array if
            # day is in the available calendary days, otherwise return nans
            cal =  bkg.bkg_data['calendar']
            is_visible = (len(cal) > 0) and (day in cal)
            day_idx = np.where(cal == day)[0][0] if day in cal else -1
            if is_visible:
                wavelength_zodi = bkg.bkg_data['wave_array']
                background_zodi = bkg.bkg_data["zodi_bg"][day_idx]
            else:
                wavelength_zodi = [np.nan]
                background_zodi = [np.nan]

            #interpolate the background array onto the grid of weights to filter
            #and sum for entire filter
            if not np.isnan(wavelength_zodi[0]):
                background_zodi_interp=np.interp(wavelength,wavelength_zodi,background_zodi)
                emission.append(float(np.sum(weights*background_zodi_interp)/len(weights)))
            else:
                emission.append(np.nan)

            #if new cache file record it for re-use
            cache_file = bkg.cache_file
            if cache_file not in cache_record :
                cache_record.append(cache_file)
                emission_record.append(emission[-1])

        #re-use emisison value from this cache file
        else:
            cachei=np.where(cache_check==np.array(cache_record))
            emission.append(emission_record[cachei[0][0]])

    return emission

#########################################
#background and comput_visibility largely lifted from StSci rbt.py module with minimal modifications
#Source:
# https://github.com/spacetelescope/roman_notebooks/blob/main/notebooks/background_visualization_tool/rbt.py

#repeat of myfile_from_healpix, so that it can be called outside the background class
# for speed to quickly check the filename of a given ra and dec
def get_healpix(nside, ra, dec):
    """Map (RA, DEC) to the cache file path via healpix indexing."""
    healpix_idx = healpy.pixelfunc.ang2pix(nside, ra, dec, nest=False, lonlat=True)
    healpix_str_pad = str(healpix_idx).zfill(6)
    return f"{healpix_str_pad[0:4]}/sl_pix_{healpix_str_pad}.bin"

#########################################

class background:
    """
    Main background class. Loads background data for a specific (RA, DEC)
    directly from the online cache at STScI.

    Parameters
    ----------
    ra : float
        Right ascension in decimal degrees
    dec : float
        Declination in decimal degrees
    wavelength : float
        Wavelength (micron)
    """

    def __init__(self, ra, dec,verbose=False):
        # Remote source (no local caching)
        self.remote_dir = "https://archive.stsci.edu/missions/roman/simulations/straylight/sl_cache/"
        self.cache_version = "2025.5"


        # Static refdata (still read from local repo files)

        #sln
        # self.local_path = Path(__file__).parent / "refdata"
        self.local_path = 'https://raw.githubusercontent.com/spacetelescope/roman_notebooks/main/notebooks/background_visualization_tool/refdata' #sln

        self.wave_file = "std_spectrum_wavelengths.txt"
        self.thermal_file = "thermal_curve_roman_rryan_v1.0.csv"

        # Healpix details used by cache partitioning
        self.nside = 128

        # Load static spectral grids / thermal background
        self.abs_wave_array, self.thermal_wave_array, self.thermal_bg = self.read_static_data(verbose=verbose)
        self.sl_abs_nwave = self.abs_wave_array.size
        self.sl_thermal_nwave = self.thermal_wave_array.size

        # Inputs
        self.ra = ra
        self.dec = dec

        # Pick the remote filename for this sky position
        self.cache_file = self.myfile_from_healpix(ra, dec)

        # Load and prepare all per-position data
        self.bkg_data = self.read_bkg_data_from_url(self.cache_file,verbose=verbose)

    # ---------- File / data loading ----------

    def myfile_from_healpix(self, ra, dec):
        """Map (RA, DEC) to the cache file path via healpix indexing."""
        healpix_idx = healpy.pixelfunc.ang2pix(self.nside, ra, dec, nest=False, lonlat=True)
        healpix_str_pad = str(healpix_idx).zfill(6)
        return f"{healpix_str_pad[0:4]}/sl_pix_{healpix_str_pad}.bin"

    def read_static_data(self,verbose=False):
        """Load static wavelength grid and thermal curve from refdata."""
        if verbose: #sln
            print(f"Loading static wavelength grid from {self.local_path + '/' + self.wave_file}.") #sln
        abs_wave_array = np.loadtxt(self.local_path + '/' + self.wave_file)
        if verbose: #sln
            print(f"Loading thermal wavelength grid from {self.local_path + '/' + self.thermal_file}.")
        thermal = np.transpose(np.genfromtxt(self.local_path + '/' + self.thermal_file, delimiter=","))

        thermal_wave_array = thermal[0]
        thermal_flux = thermal[1]

        return abs_wave_array, thermal_wave_array, thermal_flux

    # ---------- Remote cache reading ----------

    def read_bkg_data_from_url(self, cache_file, verbose=False):
        """
        Read one Roman background file (.bin) directly from the STScI-hosted URL,
        parse it into arrays, and return the structured dict.
        """

        url = self.remote_dir.rstrip("/") + "/" + cache_file
        if verbose: #sln
            print(f"Loading background file from {url}")
        with urllib.request.urlopen(url) as response:
            file_data = response.read()

        buf = io.BytesIO(file_data)

        # Dtypes for cache reading
        nonzodi_pix_dtype = np.dtype(
            [
                ("pix_ra", "f8"),
                ("pix_dec", "f8"),
                ("upos", [("x", "f8"), ("y", "f8"), ("z", "f8")]),
                ("nonzodi_bg", ("f8", self.sl_abs_nwave)),
                ("iday_index", ("i4", 366)),
            ]
        )
        zodi_sl_dtype = np.dtype(
            [
                ("zodi_bg", ("f8", self.sl_abs_nwave)),
                ("stray_light_bg", ("f8", self.sl_abs_nwave)),  # unused for Roman (often -1)
            ]
        )

        # Parse directly from in-memory buffer
        nonzodi_bg = np.frombuffer(buf.getbuffer(), dtype=nonzodi_pix_dtype, count=1, offset=0)
        offset = nonzodi_pix_dtype.itemsize
        zodi_sl_bgs = np.frombuffer(buf.getbuffer(), dtype=zodi_sl_dtype, offset=offset)

        ra = nonzodi_bg["pix_ra"]
        dec = nonzodi_bg["pix_dec"]
        pos = nonzodi_bg["upos"]
        nonzodi_bg_flux = nonzodi_bg["nonzodi_bg"][0]
        date_map = nonzodi_bg["iday_index"][0]  # 366 entries: -1 for invalid, >=0 for valid mapping

        # Initial "calendar" are the day-of-year indices with data
        calendar = np.where(date_map >= 0)[0]

        if verbose:
            print("Valid days in file:", calendar.size, "of", len(date_map))

        # Apply Roman sun-angle constraint using tgt_vis (visibility mask)
        target = SkyCoord(self.ra * u.degree, self.dec * u.degree, frame="icrs")
        c = compute_visibility(
            target, report=True, fileout="vis_debug.txt",
            interval_sampling_days=1, interval_start_time=None, interval_duration_days=366,
        )
        c.get_good_angles()
        good_indices = np.where(c.df_results["good_angles"].values)[0]  # day-of-year indices (0..365)

        # Keep only days that are both in the file and visible
        calendar = calendar[np.isin(calendar, good_indices)]

        # Map from calendar (day-of-year) to the sequential index space in zodi_sl_bgs
        # The file packs only valid days sequentially; the mapping is via date_map
        packed_index = date_map[calendar]  # guaranteed >= 0 for valid days
        # extra safety: ensure indices are in range
        packed_index = packed_index[packed_index < len(zodi_sl_bgs)]

        # Extract zodi on the valid packed indices
        zodi_bgs_full = zodi_sl_bgs["zodi_bg"]
        zodi_bgs = zodi_bgs_full[packed_index]

        Ndays = len(packed_index)

        # Interpolate to the thermal wavelength grid
        zodi_bgs_int = np.zeros((Ndays, self.sl_thermal_nwave))
        for dd in range(Ndays):
            zodi_bgs_int[dd] = self.interpolate_spec(
                self.abs_wave_array, zodi_bgs[dd], self.thermal_wave_array, fill=0.0
            )
        nonzodi_bg_int = self.interpolate_spec(
            self.abs_wave_array, nonzodi_bg_flux, self.thermal_wave_array, fill=0.0
        )

        # Base total = static (nonzodi + thermal) + zodi(day)
        total_bg = np.tile(nonzodi_bg_int + self.thermal_bg, (Ndays, 1)) + zodi_bgs_int

        # Apply NIRCam-informed modification to total and components
        for dd in range(Ndays):
            if dd == 0:
                mod_wave_array, mod_total_bg_first = self.modify_background(total_bg[dd])
                mod_total_bg = np.zeros((Ndays, len(mod_total_bg_first)))
                mod_total_bg[dd] = mod_total_bg_first
            else:
                _, mod_total_bg[dd] = self.modify_background(total_bg[dd])

        mod_zodi_bgs_int = np.zeros((Ndays, self.sl_thermal_nwave))
        for dd in range(Ndays):
            _, mod_zodi_bgs_int[dd] = self.modify_background(zodi_bgs_int[dd])
        _, mod_nonzodi_flux = self.modify_background(nonzodi_bg_int)

        return {
            "calendar": np.array(calendar),      # day-of-year indices that survived visibility gating
            "ra": ra,
            "dec": dec,
            "pos": pos,
            "wave_array": self.thermal_wave_array,
            "nonzodi_bg": mod_nonzodi_flux,
            "thermal_bg": self.thermal_bg,
            "zodi_bg": mod_zodi_bgs_int,
            "total_bg": mod_total_bg,
        }

    # ---------- Spectral modification / interpolation ----------

    def modify_background(self, bg_flux):
        """
        E. Han's modification to incorporate NIRCam measurements.
        Returns (modified_wavelength_grid, modified_flux_interpolated_to_thermal_grid).
        """
        bg_wvl_to_mod = self.thermal_wave_array
        bg_flux_to_mod = copy.deepcopy(bg_flux)

        # Smooth the hard 0.5 μm cutoff to 0.4–0.5 μm
        first_pass = np.where((bg_wvl_to_mod >= 0.4) & (bg_wvl_to_mod <= 0.5))[0]
        if first_pass.size > 1:
            bg_flux_to_mod[first_pass] = np.interp(
                bg_wvl_to_mod[first_pass],
                [bg_wvl_to_mod[first_pass[0]], bg_wvl_to_mod[first_pass[-1]]],
                [bg_flux_to_mod[first_pass[0]], bg_flux_to_mod[first_pass[-1]]],
            )

        # NIRCam pivot wavelengths and ratios (technical report)
        pivot_wvl = [0.705, 0.902, 1.154, 1.501]
        measured_bg_ratio = [0.656, 0.680, 0.819, 0.982]

        bg_model = np.interp(pivot_wvl, bg_wvl_to_mod, bg_flux_to_mod)
        measured_bg_nircam = np.array(measured_bg_ratio) * bg_model

        # Quadratic fit of ratio for extrapolation 0.41–0.624 μm
        fit = np.polyfit(pivot_wvl, measured_bg_ratio, 2)
        wvl_fix_indices = np.where((bg_wvl_to_mod >= 0.41) & (bg_wvl_to_mod <= 0.624))[0]
        bg_wvl_to_fix = bg_wvl_to_mod[wvl_fix_indices]
        solution = fit[0] * bg_wvl_to_fix**2 + fit[1] * bg_wvl_to_fix + fit[2]
        bg_flux_polyfit = solution * bg_flux_to_mod[wvl_fix_indices]

        original_blue = np.where(bg_wvl_to_mod <= 0.41)[0]
        original_red = np.where(bg_wvl_to_mod >= 1.668)[0]

        bg_wvl_mod = np.concatenate((
            bg_wvl_to_mod[original_blue],
            bg_wvl_to_fix,
            np.array(pivot_wvl),
            bg_wvl_to_mod[original_red]
        ))
        bg_flux_mod = np.concatenate((
            bg_flux_to_mod[original_blue],
            bg_flux_polyfit,
            measured_bg_nircam,
            bg_flux_to_mod[original_red]
        ))

        # Interpolate back to the thermal grid
        bg_flux_mod_interp = self.interpolate_spec(
            bg_wvl_mod, bg_flux_mod, self.thermal_wave_array, fill=0.0
        )
        return bg_wvl_mod, bg_flux_mod_interp

    def interpolate_spec(self, wave, specin, new_wave, fill=np.nan):
        """Interpolate spectral data to a new wavelength grid."""
        f = interp1d(wave, specin, bounds_error=False, fill_value=fill)
        return f(new_wave)

#########################################

class compute_visibility():

    def __init__(self,targets_coordinates,fileout=None,report=False,interval_sampling_days=None,interval_start_time=None,interval_duration_days=None):

        if isinstance(targets_coordinates, list):
            self.targets_coordinates = targets_coordinates
        else:
            self.targets_coordinates = [targets_coordinates]

        self.fileout = fileout
        self.report = report
        self.interval = {'sampling_days':interval_sampling_days,
                         'start_time':interval_start_time,
                         'duration_days':interval_duration_days}
        self.sampled_times = self.set_sampled_times(self.interval)
        self.sun_coord = get_body('Sun',self.sampled_times)   # get coordinate object for the Sun for each day of the year
        self.min_sun_angle = (90. - 36.) * u.deg
        self.max_sun_angle = (90. + 36.) * u.deg
        self.df_results = self.initialize_dataframe()


    def initialize_dataframe(self):
        radec_string = ['({}, {})'.format(coords.ra.to_string(u.hour), coords.dec.to_string(u.degree, alwayssign=True)) for coords in self.targets_coordinates]
        time_string  = [self.format_time(time) for time in self.sampled_times]
        index_levels = [radec_string,time_string]
        index_names  = ['(RA, Dec)', 'DOY']
        multi_index  = pd.MultiIndex.from_product(index_levels, names=index_names)
        column_names = ['Sun_RA','Sun_Dec','separation','good_angles','nominal_roll','pa_obs_y','pa_fpa_local_x','pa_fpa_local_y','sunang_x','sunang_y','sunang_z']

        # The reindexing of dataframe below is necessary due to this pandas bug:
        # https://stackoverflow.com/questions/71837659/trying-to-sort-multiindex-index-using-categorical-index/73766126#73766126

        return (pd.DataFrame(index=multi_index,columns=column_names)).reindex(radec_string,level=0)

    def format_time(self,time_object):
        """Converts a datetime object's time to DOY.ddddd """
        datetime = time_object.datetime
        decimal_hours = datetime.hour + datetime.minute/ 60. + datetime.second/3600
        return datetime.strftime("%Y")+'-'+datetime.strftime("%j")+'.'+'{:7.5f}'.format(decimal_hours/24)[2:]


    def set_sampled_times(self,interval):

        '''
        define the array of Time object at which the visitbility will be sampled.
        If a date is in the future, this will generate a 'dubious date' warning message
        the reason is that it is unknown how many leap seconds will be needed in the future.
        the results will still be valid
        '''

        if interval['start_time'] is None:
            t_start_str = ['2024-01-01T00:00:00.0']
            t_start = Time(t_start_str,format='isot', scale='utc')
        elif isinstance(interval['start_time'], Time) == False:
            if isinstance(interval['start_time'], str):
                t_start = Time(interval['start_time'],format='isot', scale='utc')
            else:
                print('Start time needs to be a astropy.Time object or a string')
                assert False
        else:
            t_start = interval['start_time']

        if interval['duration_days'] is None:
            t_end = 365.
        else:
            t_end = interval['duration_days']

        if interval['sampling_days'] is None:
            t_step = 1.
        else:
            t_step = interval['sampling_days']

        if t_step > t_end:
            print('sampling interval cannot exceeed total duration')
            assert False

        return t_start+np.arange(0.,t_end,t_step)*u.d


    def get_good_angles(self):
        for i,target_coordinates in enumerate(self.targets_coordinates):

            with warnings.catch_warnings():
                warnings.simplefilter('ignore', AstropyWarning)
                sun_angle = self.sun_coord.separation(target_coordinates)
                good_angles = (sun_angle >= self.min_sun_angle) & (sun_angle <= self.max_sun_angle)
            self.df_results.loc[pd.IndexSlice[self.df_results.index.levels[0][i],:],'good_angles'] = good_angles
            self.df_results.loc[pd.IndexSlice[self.df_results.index.levels[0][i],:],'separation'] = sun_angle
            self.df_results.loc[pd.IndexSlice[self.df_results.index.levels[0][i],:],'Sun_RA'] = self.sun_coord.ra
            self.df_results.loc[pd.IndexSlice[self.df_results.index.levels[0][i],:],'Sun_Dec'] = self.sun_coord.dec
