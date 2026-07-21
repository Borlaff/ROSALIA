# This notebook demonstrates how to use ROSALIA to analyze psf in a simple image. 

import rosalia as rs
import pandas as pd

from astropy.time import Time
from astropy.coordinates import SkyCoord

# Define the target
target = 'alpha Centauri'
tgt = SkyCoord.from_name(target)
ra =  219.92041  # Right ascension, in degrees. 
dec = -60.835148 # Declination, in degrees.
PA = 30 # Position angle, in degrees.
date = Time("2027-01-01T00:00:00") # Date of the observation, in Astropy Time YYYY-MM-DDTHH:MM:SS format.
bandpass = "F129"
exptime = 600 # Exposure time, in seconds.

#catalog = {"ra": [123.01, 123.02, 123.03],
#           "dec": [23.01, 23.02, 23.03],
#           "mag_lambda": [14, 14.5, 15],
#           "source_id": [1, 2, 3],
#           "cat_id": [1, 2, 3]}
           
# input_catalog = pd.DataFrame(catalog)

rosalia_psf = rs.correct.rosalia_psf(ra=ra, dec=dec, PA=PA, date=date, 
                                     bandpass=bandpass, exptime=exptime, 
                                     #input_catalog=input_catalog,
                                     g_mag_max=15, verbose=False)


# Let's plot the result. 
import matplotlib.pyplot as plt
from astropy.io import fits

hdu = fits.open("Roman_Alpha_cenA.fits")
mu = rs.detectors.fe2mu(hdu[0].data, instrument="WFI", telescope="Roman", filter_name=bandpass)
plt.imshow(mu, origin='lower', cmap='RdYlBu_r', vmin=20, vmax=28)
plt.colorbar()
plt.savefig("Roman_Alpha_cenA.png", dpi=300)
plt.close()