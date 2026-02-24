# ROSALIA: ROman Sky Analyst for Low surface brightness Imaging & Astronomy

[![Contributors](https://img.shields.io/github/contributors/Borlaff/ROSALIA.svg?style=for-the-badge)](https://github.com/Borlaff/ROSALIA/graphs/contributors)
[![Forks](https://img.shields.io/github/forks/Borlaff/ROSALIA.svg?style=for-the-badge)](https://github.com/Borlaff/ROSALIA/network/members)
[![Stargazers](https://img.shields.io/github/stars/Borlaff/ROSALIA.svg?style=for-the-badge)](https://github.com/Borlaff/ROSALIA/stargazers)
[![Issues](https://img.shields.io/github/issues/Borlaff/ROSALIA.svg?style=for-the-badge)](https://github.com/Borlaff/ROSALIA/issues)
[![License](https://img.shields.io/github/license/Borlaff/ROSALIA.svg?style=for-the-badge)](https://github.com/Borlaff/ROSALIA/blob/main/LICENSE.txt)

## About The Project

ROSALIA (Roman Sky Analyst for Low surface brightness Imaging & Astronomy) is a pipeline to model the sky background level on astronomical images obtained with [NASA/Nancy Grace Roman Space Telescope](https://roman.gsfc.nasa.gov) and its direct predecessor, the legendary [NASA/Hubble Space Telescope](https://science.nasa.gov/mission/hubble/).

In particular, ROSALIA is focused on the prediction and calibration of **stray-light** in the [Roman Wide Field Instrument](https://roman-docs.stsci.edu/roman-instruments-home/wfi-imaging-mode-user-guide/introduction-to-wfi), one of the main contaminants in ultra deep low surface brightness observations, and the main source of gradients of parasitic light for space telescopes.

ROSALIA combines the information from existing photometric catalogs (Gaia, 2MASS, WISE) with precise optical and payload ray-tracing models of the Roman Space Telescope, allowing to generate images of stray-light and other components of the sky-background for user-defined observational conditions.

ROSALIA is funded through a NASA Grant (D.14 Roman 2022), ROSES/Nancy Grace Roman Space Telescope Research and Support Participation Opportunities.

- **Sci-PI**: Alejandro S. Borlaff (NASA ARC)
- **Admin-PI**: Pamela M. Marcum (NASA ARC)

## Installation

### Managing Dependencies

ROSALIA is based on multiple packages, including [Astropy](https://www.astropy.org/), [Astroquery](https://astroquery.readthedocs.io/en/latest/), and [Romanisim](https://romanisim.readthedocs.io/en/latest/), [NumPy](https://numpy.org/), [SciPy](https://scipy.org/), and [Matplotlib](https://matplotlib.org/) among many others.

The easiest way to install all the dependencies is through a package manager like [Conda](https://anaconda.org/anaconda/conda) or [Mamba](https://github.com/mamba-org/mamba). If you have a Conda/Mamba package manager already installed in your system, skip to the next section. If you do not have a package manager, follow the Conda installation instructions at the [Space Telescope stenv environment webpage](https://stenv.readthedocs.io/en/latest/getting_started.html).

### Installing ROSALIA

#### For Developers

Create a clean environment for ROSALIA:

```bash
conda create -n rosalia python=3.12 conda-forge::astromatic-swarp
```

After the new environment is created, activate it:

```bash
conda activate rosalia
```

Once in a clean conda environment, download the ROSALIA package from GitHub.

#### For General Users

Create a clean environment for ROSALIA:

```bash
conda create -n rosalia python=3.12 conda-forge::astromatic-swarp
```

After the new environment is created, activate it:

```bash
conda activate rosalia
```

Once in a clean conda environment, install ROSALIA using pip:

```bash
pip install rosalia-wfi
```

That is it! You are ready to start analyzing Space Telescope images.

## Usage

### Minimal Use Example

ROSALIA estimates the amount of stray-light from Roman Space Telescope images. To do this, it calculates how many photons reach the focal plane array from secondary optical paths, based on a function called Normalized Detector Irradiance (NDI).

Those photons represent a source of contamination and typically must be modeled and removed before the images are ready for science. ROSALIA calculates the flux of photons for each pixel of the focal plane array. For Roman/WFI, that is a total of 300,811,392 pixels! (18 4088x4088 H4RG-10 detectors).

ROSALIA focuses on modeling three types of backgrounds:

1. **Stray-light from sources outside the field of view** (Normalized Detector Irradiance)
2. **Stray-light from sources inside the field of view** (Point Spread Function)
3. **Zodiacal light**

#### Out-of-field Stray-light

Open a Python terminal and type:

```python
import rosalia as rs 
from astropy.time import Time

# First define a Roman Space Telescope WFI exposure basic parameters.
ra = 123  # Right ascension at the center of the FOV, in degrees. 
dec = 23  # Declination at the center of the FOV, in degrees.
PA = 45   # Position angle, counter-clockwise from North, in degrees.
date = Time("2024-06-01T00:00:00")  # Date of the observation, in Astropy Time format.
bandpass = "F129"  # A string with the bandpass name for WFI. See https://roman.gsfc.nasa.gov/science/WFI_technical.html
exptime = 600  # Exposure time, in seconds.

rosalia_stray = rs.correct.rosalia_stray(ra=ra, dec=dec, PA=PA, date=date, 
                                         bandpass=bandpass, exptime=exptime, 
                                         radius=1, g_mag_max=15, 
                                         sun_block=False, verbose=False, 
                                         catalog=None)
```

While Nancy Grace Roman Space Telescope is scheduled to be launched no earlier than September 2026, you can start simulating the observations using romanisim.

#### Generating Mock Roman/WFI Observations

1. Install [romanisim](https://romanisim.readthedocs.io/en/latest/) and generate a Roman/WFI example image. For this experiment (to maximize the visualization of stray-light), simulate an exposure near Orion's Belt:

```bash
pip install romanisim
romanisim-make-image --radec 83.3419927 -1.9665163 RST_WFI_ROSALIA_test_Orion_Belt_SCA{}.asdf \
   --roll -45 --sca -1 --bandpass F158 --level 2 --usecrds
```

> **Note:** romanisim is a package in active development. Please visit the official webpage for more information on usage.

The result will be a series of 18 files (one ASDF file per Roman WFI detector, or SCA) in the local directory:

```bash
ls -lah

-rw-r--r--   1 user  staff   391M Nov 26 16:02 RST_WFI_ROSALIA_test_Orion_Belt_SCAWFI01.asdf
-rw-r--r--   1 user  staff   391M Nov 26 15:39 RST_WFI_ROSALIA_test_Orion_Belt_SCAWFI02.asdf
-rw-r--r--   1 user  staff   391M Nov 26 15:42 RST_WFI_ROSALIA_test_Orion_Belt_SCAWFI03.asdf
...
-rw-r--r--   1 user  staff   391M Nov 27 09:04 RST_WFI_ROSALIA_test_Orion_Belt_SCAWFI18.asdf
```

These 18 files represent a single Roman/WFI level 2 (calibrated, non-combined) exposure.

#### Analyzing Stray-light Level

2. Analyze the Roman/WFI example image with ROSALIA. The script `rosalia-stray` will extract all necessary information from the Roman/WFI exposure file metadata and generate a series of ASDF and FITS files with the pixel-to-pixel flux level expected for this particular exposure:

```bash
rosalia-sky RST_WFI_ROSALIA_test_Orion_Belt_SCAWFI01.asdf
```

For more examples, please refer to the [Documentation](https://rosalia.readthedocs.io/).

## ASDF Format

### What is ASDF? Where are the FITS Files?

[ASDF](https://asdf.readthedocs.io/en/latest/) is the successor of the [FITS](https://www.stsci.edu/hst/wfpc2/Wfpc2_dhb/intro_ch23.html) format and has been adopted by JWST. While GUI visualizers like [SAOImageDS9](https://ds9.si.edu/) are not yet compatible with ASDF, ROSALIA provides an easy way to extract most useful information from ASDF files through the *exposure-inspector* tool:

```bash
exposure-inspector RST_WFI_ROSALIA_test_Orion_Belt_SCAWFI01.asdf
```

`exposure-inspector` will print basic information from the ASDF tree, including:

- Telescope name, instrument, and detector
- Filter name
- Pointing information (right ascension and declination)
- Transmission curve of the filter
- World Coordinate System (WCS) of the header

## Roadmap

- [x] Automatic queries of catalogs of bright sources
  - [x] Gaia, 2MASS, WISE
  - [x] Horizons/JPL Solar System Objects
- [x] Retrieval of stray-light blocking efficiency from ray-tracing models
- [x] Ingestion of ASDF Roman/WFI simulated files
- [ ] Add diffraction modelling to Roman/WFI
- [ ] Add thermal emission model (internal stray-light)
- [ ] Complete support for Hubble Space Telescope ACS & WFC3/IR
- [ ] Automatic identification of SSOs in Roman/WFI observations

See the [open issues](https://github.com/Borlaff/ROSALIA/issues) for a full list of proposed features and known issues.

## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would improve ROSALIA, you have a few options:

1. Open an issue with the tag "enhancement"
2. Fork the repo and create a pull request
3. Email the project PIs with your ideas (a.s.borlaff@nasa.gov)

### Instructions for fork/pull contributions:

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

### Top Contributors

- [Alejandro S. Borlaff](https://github.com/Borlaff)
- [Pamela M. Marcum](https://github.com/pmarcum)
- [Scott Rohrbach](https://github.com/ScottRohrbach)

## License

ROSALIA © 2025 by Alejandro S. Borlaff is licensed under [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International](https://creativecommons.org/licenses/by-nc-sa/4.0/).

## Contact

- **Alejandro S. Borlaff** - [@asborlaff](https://bsky.app/profile/asborlaff.bsky.social) - a.s.borlaff@nasa.gov
- **Project Link**: https://github.com/Borlaff/ROSALIA

## Acknowledgments

- The Nancy Grace Roman Space Telescope Science Center
- Space Telescope Science Institute (STScI)
- NASA Ames Research Center
