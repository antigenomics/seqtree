"""Sphinx configuration for seqtree."""
import seqtree

project = "seqtree"
copyright = "2026, antigenomics"
author = "antigenomics"
# From the installed package, not a literal: these two read 0.6.1 while the package was 0.7.0.
release = seqtree.__version__
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
]

# The compiled extension is installed in the build environment, so no mocking.
autosummary_generate = False
autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_google_docstring = True
napoleon_numpy_docstring = False

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_title = "seqtree"
html_theme_options = {
    "github_url": "https://github.com/antigenomics/seqtree",
    "show_prev_next": False,
}
