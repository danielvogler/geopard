"""The import surface, and the promises the README makes about it.

The README has always shown ``import geopard; geopard.Geopard()``. Until the
package grew an ``__init__`` that did not actually work, and nothing caught it.
These tests are what stop that happening again.
"""

from __future__ import annotations

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_the_documented_import_works():
    import geopard

    assert geopard.Geopard is not None


def test_public_names_are_exported():
    import geopard

    assert set(geopard.__all__) == {
        "Geopard",
        "GeopardException",
        "GeopardResponse",
        "MatchFlag",
        "__version__",
    }
    for name in geopard.__all__:
        assert hasattr(geopard, name), name


def test_the_legacy_import_path_still_works():
    """``from geopard.geopard import Geopard`` predates the package __init__.

    Existing scripts use it, so it stays -- and it must be the same class, not
    a copy.
    """
    import geopard
    from geopard.geopard import Geopard, GeopardException, GeopardResponse

    assert Geopard is geopard.Geopard
    assert GeopardException is geopard.GeopardException
    assert GeopardResponse is geopard.GeopardResponse


@pytest.mark.parametrize(
    "module",
    [
        "geopard.dtw",
        "geopard.exceptions",
        "geopard.geodesy",
        "geopard.gpx",
        "geopard.interpolation",
        "geopard.matching",
        "geopard.plotting",
        "geopard.regions",
        "geopard.response",
        "geopard.selection",
        "geopard.settings",
    ],
)
def test_every_stage_is_importable_on_its_own(module):
    __import__(module)


def test_version_is_reported_consistently():
    import geopard
    from geopard.version import __version__

    assert geopard.__version__ == __version__


def test_pyproject_takes_its_version_from_the_package():
    """One version number, in one place, or releases go out mislabelled."""
    import tomllib

    pyproject = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))

    assert "version" in pyproject["project"]["dynamic"]
    assert pyproject["tool"]["hatch"]["version"]["path"] == "src/geopard/version.py"


def test_the_installed_distribution_reports_the_same_version():
    """Catches a build that shipped a version the source does not claim."""
    from importlib.metadata import version

    import geopard

    assert version("geopard") == geopard.__version__


def test_match_flags_keep_their_documented_integer_values():
    """Callers compare against the raw integers, and always have."""
    from geopard import MatchFlag

    assert MatchFlag.MATCH == 2
    assert MatchFlag.GREY_ZONE == 1
    assert MatchFlag.NO_MATCH == -1
    assert MatchFlag.NOT_EVALUATED == -2


def test_success_is_exactly_the_positive_flags():
    from geopard import MatchFlag

    assert all(flag > 0 for flag in (MatchFlag.MATCH, MatchFlag.GREY_ZONE))
    assert all(flag < 0 for flag in (MatchFlag.NO_MATCH, MatchFlag.NOT_EVALUATED))


def test_project_root_points_at_the_repository():
    """The examples build data paths from it, so it has to be the real root."""
    from geopard.settings import PROJECT_ROOT, PROJECT_ROOT_PATH

    assert PROJECT_ROOT_PATH == REPO_ROOT
    assert str(REPO_ROOT) == PROJECT_ROOT


def test_the_shipped_data_directories_exist_where_settings_says():
    from geopard.settings import GPX_DATA_DIR, POLYGON_DATA_DIR

    assert (GPX_DATA_DIR / "tds_sunnestube_segment.gpx").is_file()
    assert (POLYGON_DATA_DIR / "example_start_region.csv").is_file()


def test_the_geopard_class_exposes_every_documented_method():
    """The README and AGENTS.md list these. Removing one is a breaking change."""
    from geopard import Geopard

    for name in (
        "dtw_match",
        "gpx_loading",
        "gpx_track_crop",
        "nearest_neighbours",
        "spheroid_point_distance",
        "acm",
        "dtw_computation",
        "interpolate",
        "create_polygon",
        "plot_track_comparison",
        "gpx_plot",
        "parse_response",
    ):
        assert callable(getattr(Geopard, name)), name


def test_readme_shows_the_import_that_works():
    """Guards the specific line that was wrong before the package had an __init__."""
    readme = (REPO_ROOT / "README.md").read_text(encoding="utf-8")

    assert "import geopard" in readme or "from geopard import" in readme


def test_the_package_ships_a_py_typed_marker():
    """The `Typing :: Typed` classifier is a promise downstream mypy checks.

    Without this file the annotations are invisible to anyone installing it.
    """
    import geopard

    assert (Path(geopard.__file__).parent / "py.typed").is_file()
