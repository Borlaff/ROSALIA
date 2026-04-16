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
#
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

import rosalia as rs

def remove_zodiacal_light_acs(input_name, verbose=False):
    zody_interp_sci1 = get_zodiacal_background(input_name=input_name, ext=1, verbose=verbose)
    zody_interp_sci2 = get_zodiacal_background(input_name=input_name, ext=4, verbose=verbose)

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



def get_zodiacal_background(input_name, ext, wavelength=None, telescope=None, instrument=None,
                            detector=None, expstart=None, step=1000, zody_mode="zodipy",
                            nbins_wavelength=20, obslocin=3, grid_method="random",
                            output_units=None, verbose=False, interpolate=True):

    from scipy import interpolate
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

    switch_filter_curve = 0

    if (wavelength == None) or (expstart == None) or (output_units == None) or (zody_mode=="zodipy"):
        exposure_identity = rs.utils.exposure_inspector(input_name, lite=True)

    if wavelength == None:
        wavelength = exposure_identity["FILTER"]

    if telescope == None:
        telescope = exposure_identity["TELESCOP"]

    if instrument == None:
        instrument = exposure_identity["INSTRUME"]

    if detector == None:
        detector = exposure_identity["DETECTOR"]

    if expstart == None:
        expstart = exposure_identity["EXPSTART"]

    # print(exposure_identity)
    #if output_units == None:
    #    output_units = exposure_identity["BUNIT"]

    # Step zero: Check the wavelength variable.
    # If it is string, it might be a filter name.
    # Look for it in the library FILTERS
    # Check that ra and dec arguments are in array
    if isinstance(wavelength, (str)):
        filter_curve_name = rs.telescopes.find_filter_in_svo(wavelength, telescope, instrument, detector, verbose=verbose)
        rebinned_filter_curve = rebin_transmission_curve(filter_transmission_curve=filter_curve_name["filter_transmission_curve"],
                                                         nbins=nbins_wavelength, verbose=verbose)
        rebinned_transmission = rebinned_filter_curve["rebinned_transmission"]

        # Stored filter curves are in Angstrom
        rebinned_wavelength   = rebinned_filter_curve["rebinned_wavelength"]
        dlambda = rebinned_filter_curve["rebinned_dlambda"]

    # If the input is just a wavelength, then emulate the output of rebin_transmission_curve
    if isinstance(wavelength, (float)):
        if verbose:
            print("Input wavelength " + str(wavelength))
        rebinned_transmission = np.array([1])
        dlambda = np.array([1])
        rebinned_wavelength = np.array([wavelength])
    ###############################

    # We calculate the expstart
    t = Time(expstart, format='mjd', scale='utc')
    year = t.yday.split(":")[0]
    day = t.yday.split(":")[1]


    ########################################

    # First we get the detector grid for interpolating the zodiacal measurements
    detector_grid = rs.detectors.make_detector_grid(input_name=input_name, ext=ext, step=step, mode=grid_method)
    npoints_grid = len(detector_grid["grid_world"][0])
    # Then we query IRSA to get the Zody surface brightness at those positions and time
    if verbose:
        print("Year: " + str(year))
        print("Day: " + str(day))
        print("RA: " + str(np.median(detector_grid["grid_world"][0])))
        print("DEC: " + str(np.median(detector_grid["grid_world"][1])))

    # The anticipated results are len(detector_grid)*2*nbins
    db_irsa = np.zeros((npoints_grid, nbins_wavelength))*np.nan

    # For each wavelength bin, we do this query
    if zody_mode == "IRSA":
        print("Launching queries to IRSA/IPAC...")
        for i in tqdm(range(len(rebinned_wavelength))):
            # IRSA queries must go in um so we need to multiply by 1E+7 the m above
            db = rs.irsa.irsa_query(ra=detector_grid["grid_world"][0], dec=detector_grid["grid_world"][1],
                                    wavelength=rebinned_wavelength[i].to("um").value, year=year, day=day, obslocin=obslocin)

            db_irsa[:,i] = np.array(db["zody"])

    """
    if zody_mode == "gunagala":
        #gunagala_zody(ra, dec, wavelength, year, day)
        if verbose:
            print("Estimating zodiacal light with Gunagala...")
        for i in tqdm(range(len(rebinned_wavelength))):
            # IRSA queries must go in um so we need to multiply by 1E+7 the m above
            db = gunagala_zody(ra=detector_grid["grid_world"][0], dec=detector_grid["grid_world"][1],
                                    wavelength=rebinned_wavelength[i].to("AA").value, year=year, day=day)

            db_irsa[:,i] = np.array(db["zody"])


    if zody_mode == "IRSA" or zody_mode == "gunagala":
        # Now we calculate the integrated Zodiacal light
        # https://www.aanda.org/articles/aa/pdf/2022/01/aa41935-21.pdf
        # μVIS,AB = −2.5 log10( ∑ f (νi) (h νi)−1 e(νi) ∆νi/A ∑(h νi)−1 e(νi) ∆νi)+ 8.90.

        # T_array = np.loadtxt("VIS_full_transmission_IRSA_bins_v2.txt", delimiter=",")
        # evdv = (T_array[:,2]) * evp.cs*T_array[:,1]/(T_array[:,0]**2)
        # evdv = Transmission * c lambda / dlambda**2
        evdv = rebinned_transmission*rebinned_wavelength*constants.c/dlambda**2

        # We integrate the flux
        zodiacal_flux_convolved_by_filter = np.zeros(npoints_grid)*np.nan

        for i in range(npoints_grid):
            zodiacal_flux_convolved_by_filter[i] = np.sum(db_irsa[i,:]*evdv)/np.sum(evdv)
    """

    if zody_mode == "zodipy":
        #gunagala_zody(ra, dec, wavelength, year, day)
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

def gunagala_zody(ra, dec, wavelength, year, day):
    """
    This program makes use of the Gunagala Zodiacal model to estimate the
    Zodiacal light surface brightness in a given position, at a wavelength,
    day and year.

    Input:
    ra = Right ascension (degrees)
    dec = Declination (degrees)
    wavelength = wavelength (micron)
    year = year of the observation
    day = day of the observation


    Output:
    mu_zody = surface brightness in Jy arcsec-2
    """

    if not isinstance(ra, (list, pd.core.series.Series, np.ndarray)):
        ra = np.array([ra])
    if not isinstance(dec, (list, pd.core.series.Series, np.ndarray)):
        dec = np.array([dec])

    # If wavelength, year, or day are not in array, copy their values into one
    # as large as ra, dec
    if not isinstance(wavelength, (list, pd.core.series.Series, np.ndarray)):
        wavelength = np.array([wavelength]*len(ra))

    time = Time("" + str(year) + ":" + str(day) + ":00:00:00.000", scale='utc')

    if not isinstance(year, (list, pd.core.series.Series, np.ndarray)):
        year = np.array([year]*len(ra))

    if not isinstance(day, (list, pd.core.series.Series, np.ndarray)):
        day = np.array([day]*len(ra))


    zodi_gungala = skies.ZodiacalLight()

    skycoord = SkyCoord(ra, dec, unit='deg', frame='icrs')
    wavelength = wavelength * u.angstrom
    zodi_relative_to_poles = zodi_gungala.relative_brightness(skycoord, time)

    #print(wavelength)
    f_zodi_gunagala_at_the_ecliptic_poles = zodi_gungala.surface_brightness()
    #print(f_zodi_gunagala_at_the_ecliptic_poles(wavelength))
    zodi_intensity_poles_Mysr = f_zodi_gunagala_at_the_ecliptic_poles(wavelength).to(u.MJy * u.steradian**-1, equivalencies=u.spectral_density(wavelength))
    #print(zodi_relative_to_poles*zodi_intensity_poles_Mysr)
    zodi_intensity_poles_jyarcsec2 = f_zodi_gunagala_at_the_ecliptic_poles(wavelength).to(u.Jy * u.arcsec**-2, equivalencies=u.spectral_density(wavelength))
    #print(zodi_relative_to_poles*zodi_intensity_poles_jyarcsec2)
    zodi_intensity_pointing =  zodi_relative_to_poles*zodi_intensity_poles_jyarcsec2

    d = {}
    d['ra'] = ra
    d['dec'] = dec
    d['wavelength'] = wavelength
    d['year'] = year
    d['day'] = day
    d['zody'] = zodi_intensity_pointing


    return(d)

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

