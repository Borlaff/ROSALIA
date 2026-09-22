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

observer={"TELESCOP": "Roman/WFI", 
          "pointing": [ra, dec], 
          "FILTER":bandpass, 
          "PA_Y": PA, 
          "EXPSTART": date.mjd, 
          "EXPTIME": exptime}

prefix = "alphacen"
custom_roman_exposure = rs.core.exposure(observer=observer, prefix=prefix) 

# psf_out = custom_roman_exposure.psf()