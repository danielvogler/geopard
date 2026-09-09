"""Shared fixtures.

Every path here is derived from this file's own location rather than from
``geopard.settings.PROJECT_ROOT``, so the suite keeps working if the package
layout moves. ``PROJECT_ROOT`` itself is asserted separately, in
``test_public_api.py``, where it is the thing under test rather than the way
the test finds its data.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from geopard import Geopard

REPO_ROOT = Path(__file__).resolve().parent.parent

# Relative tolerance for dynamic time warping distances.
#
# ``interpolate`` perturbs every coordinate by up to 1e-7 degrees to keep the
# spline solver away from duplicate points, and it draws that noise from the
# unseeded global numpy RNG. Six consecutive runs of the Sunnestube match
# spread over 0.1% of the mean, so 1% is loose enough never to flake and tight
# enough that a real numerical regression fails. Everything else the matcher
# returns -- flags, times, trackpoint coordinates -- is bit-stable, and is
# asserted exactly.
DTW_RTOL = 1e-2


@pytest.fixture(autouse=True)
def _quiet_logging():
    """Silence the library's INFO chatter, which is per-trackpoint and loud."""
    logging.disable(logging.CRITICAL)
    yield
    logging.disable(logging.NOTSET)


@pytest.fixture(scope="session")
def gp() -> Geopard:
    """The library entry point. It holds no state, so one instance will do."""
    return Geopard()


@pytest.fixture(scope="session")
def gpx_dir() -> Path:
    """Example GPX tracks shipped with the repository."""
    return REPO_ROOT / "data" / "gpx_files"


@pytest.fixture(scope="session")
def polygon_dir() -> Path:
    """Example start/finish region polygons shipped with the repository."""
    return REPO_ROOT / "data" / "csv_polygon_files"


@pytest.fixture(scope="session")
def fixture_dir() -> Path:
    """Small hand-written fixtures used only by the tests."""
    return Path(__file__).resolve().parent / "data"
