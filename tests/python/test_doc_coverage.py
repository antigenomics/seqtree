"""Doc coverage: every public symbol is documented, exported, and reachable from the API page.

Sphinx's own coverage builder only sees what ``api.rst`` asks it to autodoc, so a module that is
never referenced looks perfectly covered. These checks come at it from the package side instead.
"""
import dataclasses
import importlib
import inspect
import pathlib
import pkgutil

import pytest

import seqtree

DOCS = pathlib.Path(__file__).resolve().parents[2] / "docs"
SUBMODULES = ["gapblock", "pairwise", "seeds", "evalue", "layout", "pmhc", "control",
              "pmhc_evalue"]

#: Modules with no public functions/classes of their own need no page of their own.
NO_PAGE_REQUIRED = {"control", "pmhc_evalue"}


def _public(module):
    """Public functions and classes *defined in* this module (not re-exports)."""
    names = getattr(module, "__all__", None)
    if names is None:
        names = [n for n in dir(module) if not n.startswith("_")]
    out = []
    for n in names:
        obj = getattr(module, n, None)
        if inspect.isfunction(obj) or inspect.isclass(obj):
            if getattr(obj, "__module__", "").startswith("seqtree"):
                out.append((n, obj))
    return out


def test_every_package_export_resolves():
    missing = [n for n in seqtree.__all__ if not hasattr(seqtree, n)]
    assert not missing, f"seqtree.__all__ names nothing: {missing}"


def test_no_public_symbol_is_missing_from_all():
    """A public name that is not in __all__ is invisible to `from seqtree import *` and to docs."""
    public = {n for n in dir(seqtree)
              if not n.startswith("_") and not inspect.ismodule(getattr(seqtree, n))}
    undeclared = sorted(public - set(seqtree.__all__))
    assert not undeclared, f"public but not exported in __all__: {undeclared}"


@pytest.mark.parametrize("mod", SUBMODULES)
def test_submodule_has_a_module_docstring(mod):
    m = importlib.import_module(f"seqtree.{mod}")
    assert (m.__doc__ or "").strip(), f"seqtree.{mod} has no module docstring"


@pytest.mark.parametrize("mod", SUBMODULES)
def test_every_public_symbol_is_documented(mod):
    m = importlib.import_module(f"seqtree.{mod}")
    bare = [n for n, obj in _public(m) if not (obj.__doc__ or "").strip()]
    assert not bare, f"undocumented public symbols in seqtree.{mod}: {bare}"


@pytest.mark.parametrize("mod", SUBMODULES)
def test_every_public_method_is_documented(mod):
    """``__init__`` is exempt: dataclasses generate it, and the class docstring documents it."""
    m = importlib.import_module(f"seqtree.{mod}")
    bare = []
    for cls_name, cls in _public(m):
        if not inspect.isclass(cls):
            continue
        assert (cls.__doc__ or "").strip(), f"seqtree.{mod}.{cls_name} has no class docstring"
        if dataclasses.is_dataclass(cls):
            continue
        for name, fn in inspect.getmembers(cls, inspect.isfunction):
            if name.startswith("_"):
                continue
            if fn.__qualname__.split(".")[0] != cls_name:
                continue        # inherited
            if not (fn.__doc__ or "").strip():
                bare.append(f"{cls_name}.{name}")
    assert not bare, f"undocumented public methods in seqtree.{mod}: {bare}"


def test_every_submodule_is_reachable_from_the_docs():
    """A module nobody automodules is a module nobody reads. Any page counts, not just api.rst."""
    pages = "\n".join(p.read_text() for p in DOCS.glob("*.rst"))
    absent = [m for m in SUBMODULES
              if m not in NO_PAGE_REQUIRED and f"seqtree.{m}" not in pages]
    assert not absent, f"no docs page references: {absent}"


def test_no_submodule_was_added_without_being_listed_here():
    """Guards this file itself: a new seqtree/*.py must join SUBMODULES and the docs."""
    found = {name for _, name, ispkg in pkgutil.iter_modules(seqtree.__path__)
             if not ispkg and not name.startswith("_")}
    assert found <= set(SUBMODULES), (
        f"new submodule(s) {sorted(found - set(SUBMODULES))}: add to SUBMODULES and docs/api.rst"
    )


def test_examples_compile_and_carry_a_docstring():
    """The examples are documentation. If one stops parsing, the docs are lying."""
    import py_compile

    root = pathlib.Path(__file__).resolve().parents[2] / "examples"
    scripts = sorted(root.glob("*.py"))
    assert len(scripts) >= 3, f"expected the three worked examples, found {scripts}"
    for s in scripts:
        py_compile.compile(str(s), doraise=True)
        body = "\n".join(ln for ln in s.read_text().split("\n") if not ln.startswith("#!"))
        assert body.lstrip().startswith('"""'), f"{s.name} has no module docstring"


def test_examples_page_references_every_script():
    page = (DOCS / "examples.rst").read_text()
    root = pathlib.Path(__file__).resolve().parents[2] / "examples"
    missing = [s.name for s in sorted(root.glob("*.py")) if s.name not in page]
    assert not missing, f"docs/examples.rst does not mention: {missing}"
