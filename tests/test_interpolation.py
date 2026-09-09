"""Resampling a track onto a fixed number of equidistant points.

Both curves must have the same point count before dynamic time warping can
compare them, which is the whole reason this step exists.
"""

from __future__ import annotations

import numpy as np
import pytest

INTERPOLATION_POINTS = 1000


def test_returns_a_fixed_number_of_two_dimensional_points(gp, gpx_dir):
    track = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))

    points = gp.interpolate(track)

    assert points.shape == (INTERPOLATION_POINTS, 2)


def test_output_is_latitude_then_longitude(gp, gpx_dir):
    """Same order as the input rows -- a swap here silently ruins every match."""
    track = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))

    points = gp.interpolate(track)

    assert points[:, 0].min() > 47.0
    assert points[:, 0].max() < 48.0
    assert points[:, 1].min() > 9.0
    assert points[:, 1].max() < 10.0


def test_a_short_and_a_long_track_come_out_the_same_length(gp, gpx_dir):
    """61 points and 1616 points both become 1000. That is the point."""
    short = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))
    long = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"))

    assert short.shape[1] != long.shape[1]
    assert gp.interpolate(short).shape == gp.interpolate(long).shape


def test_endpoints_stay_where_they_were(gp, gpx_dir):
    """The spline is fitted through the data, so it starts and ends on it."""
    track = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))

    points = gp.interpolate(track)

    assert points[0, 0] == pytest.approx(track[0][0], abs=1e-5)
    assert points[0, 1] == pytest.approx(track[1][0], abs=1e-5)
    assert points[-1, 0] == pytest.approx(track[0][-1], abs=1e-5)
    assert points[-1, 1] == pytest.approx(track[1][-1], abs=1e-5)


def test_stays_within_the_bounding_box_of_the_track(gp, gpx_dir):
    """A spline can overshoot; on these tracks it must not run off the map."""
    track = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))
    lat = np.array(track[0], dtype=float)
    lon = np.array(track[1], dtype=float)

    points = gp.interpolate(track)

    margin = 1e-3
    assert points[:, 0].min() >= lat.min() - margin
    assert points[:, 0].max() <= lat.max() + margin
    assert points[:, 1].min() >= lon.min() - margin
    assert points[:, 1].max() <= lon.max() + margin


def test_result_is_reproducible_under_a_seeded_rng(gp, gpx_dir):
    """Interpolation jitters each point by up to 1e-7 deg to avoid duplicates.

    That noise comes from the global numpy RNG, so seeding it makes the whole
    pipeline reproducible. This test is what documents that -- and what would
    fail if the jitter were ever moved to a private generator.
    """
    track = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))

    np.random.seed(1234)  # noqa: NPY002 - pins the library's own global-RNG path
    first = gp.interpolate(track)
    np.random.seed(1234)  # noqa: NPY002 - pins the library's own global-RNG path
    second = gp.interpolate(track)

    assert np.array_equal(first, second)


def test_seeded_output_matches_the_recorded_baseline(gp, gpx_dir):
    """Golden values recorded from the pre-refactor implementation."""
    track = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))

    np.random.seed(1234)  # noqa: NPY002 - pins the library's own global-RNG path
    points = gp.interpolate(track)

    assert points[0].tolist() == pytest.approx([47.14997993830391, 9.149930024421755], abs=1e-9)
    assert points[-1].tolist() == pytest.approx([47.1642400453317, 9.176630080017569], abs=1e-9)
    assert points[:, 0].sum() == pytest.approx(47156.34808628759, abs=1e-4)
    assert points[:, 1].sum() == pytest.approx(9163.825149270324, abs=1e-4)


def test_the_jitter_is_smaller_than_gps_precision(gp, gpx_dir):
    """1e-7 degrees is ~1 cm. Anything larger would change match outcomes."""
    track = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))

    runs = [gp.interpolate(track) for _ in range(3)]
    spread = max(np.abs(runs[i] - runs[0]).max() for i in (1, 2))

    assert spread < 1e-5


def test_an_explicit_seed_is_reproducible(gp, gpx_dir):
    """``seed=`` pins the jitter without touching the global RNG."""
    track = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))
    from geopard.interpolation import interpolate

    assert np.array_equal(interpolate(track, seed=7), interpolate(track, seed=7))


def test_different_seeds_give_different_jitter(gp, gpx_dir):
    track = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))
    from geopard.interpolation import interpolate

    assert not np.array_equal(interpolate(track, seed=7), interpolate(track, seed=8))


def test_an_explicit_seed_leaves_the_global_rng_alone(gp, gpx_dir):
    """Seeding geopard must not reach into a caller's own random stream."""
    track = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))
    from geopard.interpolation import interpolate

    np.random.seed(99)  # noqa: NPY002 - the state under test
    expected = np.random.random()  # noqa: NPY002

    np.random.seed(99)  # noqa: NPY002
    interpolate(track, seed=7)

    assert np.random.random() == expected  # noqa: NPY002


def test_a_custom_point_count_is_honoured(gp, gpx_dir):
    track = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))
    from geopard.interpolation import interpolate

    assert interpolate(track, points=250).shape == (250, 2)
