# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "openpharmmdflow"
copyright = "2024, OMSF"
author = "OMSF"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions: list[str] = ["nbsphinx", "nbsphinx_link", "sphinx.ext.autodoc"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
# see https://github.com/sphinx-doc/sphinx/issues/12300
# pickling environment... WARNING: cannot cache unpickable configuration value:
# 'nbsphinx_custom_formats' (because it contains a function, class, or module object)
# [config.cache]
suppress_warnings = ["config.cache"]

autodoc_mock_imports = ["openff", "pooch"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"
# html_static_path = ["_static"]
