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


def _teach_autodoc_about_nanobind_static_methods():
    """Make ``@staticmethod``-equivalent nanobind functions document as methods.

    nanobind binds instance methods as ``nb_method`` (a method descriptor, which
    ``inspect.isroutine`` accepts) but ``def_static`` members as ``nb_func``, which it does not.
    autodoc therefore classified every static method as a *data attribute* and rendered its
    ``repr()`` -- so ``Index.build``, ``TextIndex.build``, ``load``, and all seven
    ``SubstitutionMatrix`` factories appeared on the API page as
    ``build = <nanobind.nb_func object>``, with their signature and docstring dropped. Sixteen
    members, including the library's main entry point.

    ``MethodDocumenter.can_document_member`` asks ``sphinx.util.inspect.isroutine``, so widening
    that one predicate is enough.
    """
    from sphinx.util import inspect as sphinx_inspect

    try:
        from seqtree._core import Index
    except ImportError:  # pragma: no cover - only when the extension is not built
        return
    nb_func = type(Index.__dict__["build"])
    if nb_func.__name__ != "nb_func":  # nanobind changed; leave autodoc alone
        return

    original = sphinx_inspect.isroutine

    def isroutine(obj):
        return isinstance(obj, nb_func) or original(obj)

    sphinx_inspect.isroutine = isroutine


_teach_autodoc_about_nanobind_static_methods()

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_title = "seqtree"
html_theme_options = {
    "github_url": "https://github.com/antigenomics/seqtree",
    "show_prev_next": False,
    # Two levels of page headings in the right-hand "On this page", not one: most pages here
    # carry subsections that are the actual thing a reader is looking for.
    "show_toc_level": 2,
    "navigation_depth": 3,
}

# The left sidebar rendered as a bare "Section Navigation" heading with nothing under it: the
# toctree was flat, so every page was a top-level sibling in the navbar and no page had a
# subtree to show. The captioned groups above give it the whole site, grouped, on every page.
html_sidebars = {
    "**": ["site-nav"],
    "index": [],  # the landing page has its own card grid; a duplicate tree adds nothing
}
