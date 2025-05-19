# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
from datetime import datetime
import tomllib
import os
import sys # Path manipulation

# Read information from pyproject.toml
with open("../../pyproject.toml", "rb") as _f:
    _config = tomllib.load(_f)

__project = _config['project']
__year = datetime.now().year

proj_name = __project['name']
proj_release = __project['version']
proj_license = __project['license']
proj_authors = ', '.join([d['name'] for d in __project['authors']])
proj_copyright = f"{__year}, {proj_authors}"
project = proj_name
copyright = proj_copyright
license = proj_license
author = proj_authors
release = proj_release

#
# project = 'acspype'
# copyright = '2025, Ian Black, 2025, Thor Black'
# license = 'MIT'
# author = 'Ian Black, Thor Black'
# release = '0.3.0'

sys.path.insert(0, os.path.abspath('../..')) # Adjust the path to point to the root of your project

extensions = [
    'sphinx.ext.autodoc', # Automatically generate documentation from docstrings
    'sphinx.ext.napoleon', # Google style docstrings
    'sphinx.ext.viewcode', # Add links to _source code
    'm2r2',  # m2r2 needs to be before myst_nb
    'myst_nb', # Jupyter notebook support, should also cover markdown
]

nb_execution_mode = 'off' # Disable execution of code cells in notebooks

templates_path = ['_templates']
exclude_patterns = []

source_suffix = [
    '.rst',
    '.md',
]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme' # testing new theme, I mean the package IS written in python
html_static_path = ['_static']
