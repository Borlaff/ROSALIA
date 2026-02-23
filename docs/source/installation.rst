Installation
============

Preparing your system
-----------------------
..
  Example of link `CNN <http://cnn.com>`_

ROSALIA is based on multiple packages, including `Astropy <https://www.astropy.org/>`_, `Astroquery <https://astroquery.readthedocs.io/en/latest/>`_, `Romanisim <https://romanisim.readthedocs.io/en/latest/>`_, `NumPy <https://numpy.org/>`_, `SciPy <https://scipy.org>`_, `Matplotlib <https://matplotlib.org/>`_ among many others. The easiest way to install all the dependencies is through a package manager like `Conda <https://anaconda.org/anaconda/conda>`_ or `Mamba <https://github.com/mamba-org/mamba>`_. If you have a *Conda/Mamba* package manager already installed in your system, skip to the following section. If you do not have a package manager, follow the Conda installation instructions at the `Space Telescope *stenv* environment webpage <https://stenv.readthedocs.io/en/latest/getting_started.html>`_.

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

That is it! We are ready to start analyzing Space Telescope images.
