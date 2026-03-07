from pathlib import Path
from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Compiler import Options

Options.docstrings = True
Options.annotate = False

scripts = [str(s) for s in Path('bin/').iterdir()
           if s.is_file() and s.name != '__pycache__']

setup(scripts=scripts)