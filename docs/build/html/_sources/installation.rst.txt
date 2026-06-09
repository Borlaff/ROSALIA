Installation
============

Preparing your system
-----------------------
..
  Example of link `CNN <http://cnn.com>`_

ROSALIA is based on multiple packages, including `Astropy <https://www.astropy.org/>`_, `Astroquery <https://astroquery.readthedocs.io/en/latest/>`_, `Romanisim <https://romanisim.readthedocs.io/en/latest/>`_, `NumPy <https://numpy.org/>`_, `SciPy <https://scipy.org>`_, `Matplotlib <https://matplotlib.org/>`_ among many others. The easiest way to install all the dependencies is through a package manager like `Conda <https://anaconda.org/anaconda/conda>`_ or `Mamba <https://github.com/mamba-org/mamba>`_. If you have a *Conda/Mamba* package manager already installed in your system, skip to the following section. If you do not have a package manager, follow the Conda installation instructions at the `Space Telescope *stenv* environment webpage <https://stenv.readthedocs.io/getting_started.html>`_. 

There are many versions of Conda managers. We recommend using free `Miniforge <https://github.com/conda-forge/miniforge>` with mamba as the package manager. You can install Miniforge by following the instructions at the `Miniforge GitHub repository <https://github.com/conda-forge/miniforge>`_. 

Installing ROSALIA
------------------
Create a clean environment for ROSALIA

.. code-block::

       conda create -n rosalia python=3.12 conda-forge::astromatic-swarp

After the new environment is created, we can activate it.

.. code-block::

       conda activate rosalia


Once in a clean conda environment, we can install ROSALIA. The preferred method to install it is through pip.

.. code-block::

       pip install rosalia

ROSALIA needs a series of cache files to work. If they are not installed, most of the functions will work, but you won't be able to estimate the surface brightness of the stray light, which is one of the main functionalities of ROSALIA. To install the cache files, download the following folder:

`ROSALIA Cache <https://zenodo.org/records/20335392/files/rosalia-cache.tar.gz>`_,

Then, extract the contents of the downloaded file and move the resulting folder to a location of your choice. Finally, set the environment variable ``ROSALIACACHE`` to point to the location of the extracted folder. For example, if you extracted the folder to ``/path/to/rosalia-cache``, you can set the environment variable as follows:
.. code-block:: sh

   export ROSALIACACHE=/path/to/rosalia-cache

CRDS installation
---------------------
ROSALIA relies on the `CRDS <https://roman-crds.stsci.edu/>`_ package to access the calibration reference files for the Roman Space Telescope. CRDS is installed by default as a ROSALIA dependency, but it needs to be configured to work properly. To configure CRDS, we need to set the environment variable ``CRDS_PATH`` to point to a directory where CRDS can store its cache files. For example, you can create a directory called ``crds_cache`` in your home directory and set the environment variable as follows:

.. code-block:: sh

   export CRDS_PATH=$HOME/crds_cache
   export CRDS_SERVER_URL=https://roman-crds.stsci.edu


That is it! We are ready to start analyzing Roman Space Telescope images.


