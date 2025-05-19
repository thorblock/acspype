# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os 
import sys # Path manipulation


project = 'acspype'
copyright = '2025, Ian Black, Thor Black'
license = 'MIT'
author = 'Ian Black, Thor Black'
release = '0.2.9'

sys.path.insert(0, os.path.abspath('../..')) # Adjust the path to point to the root of your project

extensions = [
    'sphinx.ext.autodoc', # Automatically generate documentation from docstrings
    'sphinx.ext.napoleon', # Google style docstrings
    'sphinx.ext.viewcode', # Add links to _source code
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
