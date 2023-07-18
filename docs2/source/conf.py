# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

sys.path.insert(0, os.path.abspath('..'))

project = 'EUReCA'
copyright = '2023, Enrico Prataviera, Jacopo Vivian'
author = 'Enrico Prataviera, Jacopo Vivian'
release = '0.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
              # 'sphinx_rtd_theme',
              # 'sphinx.ext.duration',
              # 'sphinx.ext.doctest',
              'sphinx.ext.autodoc',
              # 'sphinx.ext.autosummary',
              # 'sphinx.ext.napoleon',
              #'myst_parser',
              ]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
