"""
ROSALIA: ROman Sky Analyst for Low surface brightness Imaging & Astronomy

"""

__name__ = "ROSALIA"
__version__ = "0.9.6"
__author__ = "Alejandro S. Borlaff"
__author_email__ = "a.s.borlaff@nasa.gov"
__description__ = "A software to calibrate the sky background of Space Telescope images"
__url__ = "https://github.com/Borlaff/STRAYCOR"

import rosalia
import rosalia.utils
import rosalia.astrometry
import rosalia.bootima
import rosalia.correct
import rosalia.detectors
import rosalia.horizons
import rosalia.irsa
import rosalia.mast
import rosalia.ndi
import rosalia.psf
import rosalia.gaia
import rosalia.gnu
import rosalia.point
import rosalia.roman
import rosalia.sky
import rosalia.skysurf
import rosalia.sso
import rosalia.skywalker
import rosalia.telescopes
import rosalia.tests # No need to run the tests every time we load the module.
import rosalia.constants
import rosalia.plots
import rosalia.albedo
import rosalia.render
import rosalia.attitude
