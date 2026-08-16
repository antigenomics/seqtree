"""The version a user sees must be the version that was packaged.

``seqtree.__version__`` is read from the installed distribution metadata rather than written as a
literal, so these are cheap. They exist because the failure they catch is silent: a release commit
that bumps ``pyproject.toml`` alone ships a package whose ``__version__`` still names the previous
release, and nothing anywhere raises. ``docs/conf.py`` had drifted exactly that way (0.6.1 against
a 0.7.0 package) before it was pointed at ``__version__`` too.
"""
import pathlib
from importlib.metadata import version as distribution_version

import pytest

import seqtree

PYPROJECT = pathlib.Path(__file__).resolve().parents[2] / "pyproject.toml"


def test_version_is_exposed_and_looks_like_a_version():
    assert isinstance(seqtree.__version__, str)
    assert seqtree.__version__.split(".")[0].isdigit(), seqtree.__version__


def test_version_matches_the_installed_distribution_metadata():
    assert seqtree.__version__ == distribution_version("seqtree")


def test_version_matches_pyproject():
    """The single source. Skipped off a source checkout (a wheel ships no pyproject.toml)."""
    tomllib = pytest.importorskip("tomllib")        # 3.11+; CI also runs 3.10
    if not PYPROJECT.exists():
        pytest.skip("not a source checkout")
    declared = tomllib.loads(PYPROJECT.read_text())["project"]["version"]
    assert seqtree.__version__ == declared, (
        f"seqtree.__version__ is {seqtree.__version__} but pyproject.toml declares {declared}; "
        "the installed package is stale -- reinstall, or the release bumped only one of them"
    )
