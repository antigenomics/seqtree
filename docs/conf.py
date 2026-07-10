"""Sphinx configuration for seqtree."""

project = "seqtree"
copyright = "2026, antigenomics"
author = "antigenomics"
release = "0.3.0"
version = "0.3.0"

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
