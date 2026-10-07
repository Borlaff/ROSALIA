"""
ROSALIA: ROman Sky Analyst for Low surface brightness Imaging & Astronomy

Alejandro S. Borlaff
NASA Ames Research Center, Moffett Field, 94035, California, USA.
a.s.borlaff@nasa.gov
"""

from importlib.metadata import version, PackageNotFoundError
from pathlib import Path
import re

try:
    __version__ = version("rosalia-wfi")
except PackageNotFoundError:
   __version__ = "1.2.4.dev0"


__author__ = "Alejandro S. Borlaff"
__author_email__ = "a.s.borlaff@nasa.gov"
__description__ = "A software to calibrate the sky background of Space Telescope images"
__url__ = "https://github.com/Borlaff/ROSALIA"

from . import (
    albedo as albedo,
    astrometry as astrometry,
    attitude as attitude,
    bootima as bootima,
    constants as constants,
    core as core,
    correct as correct,
    detectors as detectors,
    gaia as gaia,
    gnu as gnu,
    horizons as horizons,
    hst as hst,
    irsa as irsa,
    ndi as ndi,
    plots as plots,
    psf as psf,
    render as render,
    roman as roman,
    sky as sky,
    skysurf as skysurf,
    skywalker as skywalker,
    sso as sso,
    telescopes as telescopes,
    tests as tests,
    utils as utils,
)

__all__ = [
    "utils",
    "astrometry",
    "bootima",
    "correct",
    "core",
    "detectors",
    "horizons",
    "irsa",
    "ndi",
    "psf",
    "gaia",
    "gnu",
    "roman",
    "sky",
    "skysurf",
    "sso",
    "skywalker",
    "telescopes",
    "tests",
    "constants",
    "plots",
    "albedo",
    "render",
    "attitude",
    "roman_estimate_straylight_SCA"
]