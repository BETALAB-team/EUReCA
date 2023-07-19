# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys


# def skip(app, what, name, obj, would_skip, options):
#     if name == "__init__":
#         return False
#     return would_skip
#
# def setup(app):
#     app.connect("autodoc-skip-member", skip)

sys.path.insert(0, os.path.abspath(os.path.join('..', '..')))

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
              'sphinx.ext.autosummary',
              'sphinx.ext.napoleon',
              'sphinx.ext.viewcode',
              'sphinx.ext.todo',
              'sphinx_mdinclude',
              #'myst_parser',
              ]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
# autoclass_content = 'both'
# autosummary_generate = True
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    # 'undoc-members': True,
    # 'exclude-members': '__weakref__'
}
