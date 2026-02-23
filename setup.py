import os
from setuptools import setup
import setuptools
# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "ROSALIA",
    version = "0.9.4",
    author = "Alejandro S. Borlaff",
    author_email = "a.s.borlaff@nasa.gov",
    description = ("A software to calibrate the sky background of Space Telescope images"),
    license = "BSD",
    keywords = "Hubble / Roman LSB sky-background",
    url = "https://github.com/Borlaff/ROSALIA",
    packages=setuptools.find_packages(where=".", exclude=()),
    package_data={'': ['*.mplstyle', '*.csv', '*.txt']},
    install_requires=['multiprocess', 'bottleneck', 'xmltodict',
                      'romanisim @ git+https://github.com/spacetelescope/romanisim.git', 'asdf',
                      'cartopy',
                      'healpy', 'regions', 'reproject', 'LSSTDESC.Coord',  'numpy',#  'numpy==1.26.4',
                      'galsim', 'skyfield', 'tqdm', 'pybind11>=2.12', 'pandas',
                      'requests', 'numexpr', 'astroquery', 'scipy', 'matplotlib',
                      'astropy_healpix', 'pytest', 'celluloid', 'ipython', 'pysynphot',
                      'sphericalpolygon', 'psutil', 'zodipy'],
    long_description=read('README.md'),
    scripts=['bin/rosalia-sky', 'bin/run-rosalia-tests', 'bin/exposure-inspector', 'bin/rosalia-stray', 'bin/test-rosalia-parser'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
)
