"""End-to-end segment matching -- the benchmark the refactor is measured against.

Every number here was recorded from the implementation as it stood before any
restructuring. Flags, elapsed times and trackpoint coordinates are asserted
exactly; only the DTW distance carries a tolerance, because interpolation
jitters coordinates from the unseeded global RNG (see ``DTW_RTOL`` in
``conftest.py``).
"""

from __future__ import annotations

import pytest

from conftest import DTW_RTOL

MATCH_WITHIN_THRESHOLD = 2
MATCH_IN_GREY_ZONE = 1
NO_MATCH = -1
COULD_NOT_BE_EVALUATED = -2


@pytest.mark.slow
def test_green_marathon(gp, gpx_dir):
    result = gp.dtw_match(
        str(gpx_dir / "green_marathon_segment.gpx"),
        str(gpx_dir / "green_marathon_activity_4_15_17.gpx"),
        radius=20,
    )

    assert result.match_flag == MATCH_WITHIN_THRESHOLD
    assert result.is_success()
    assert result.dtw == pytest.approx(0.134818, rel=DTW_RTOL)
    assert result.time.seconds == 15303
    assert result.start_point[0] == 47.376267
    assert result.start_point[1] == 8.535439
    assert result.end_point[0] == 47.375947
    assert result.end_point[1] == 8.535789
    assert result.error is None


@pytest.mark.slow
def test_sunnestube(gp, gpx_dir):
    result = gp.dtw_match(
        str(gpx_dir / "tds_sunnestube_segment.gpx"),
        str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"),
        radius=7,
    )

    assert result.match_flag == MATCH_WITHIN_THRESHOLD
    assert result.is_success()
    assert result.dtw == pytest.approx(0.096486, rel=DTW_RTOL)
    assert result.time.seconds == 1534
    assert result.start_point[0] == 47.150031
    assert result.start_point[1] == 9.149962
    assert result.end_point[0] == 47.164231
    assert result.end_point[1] == 9.176556
    assert result.error is None


@pytest.mark.slow
def test_sunnestube_with_start_and_finish_regions(gp, gpx_dir, polygon_dir):
    """Polygon regions pick different endpoints than a circle does.

    Both are valid matches of the same segment; the numbers differ because the
    regions admit a slightly different stretch of the activity, and this test
    exists to keep that difference intentional.
    """
    result = gp.dtw_match(
        str(gpx_dir / "tds_sunnestube_segment.gpx"),
        str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"),
        start_region=gp.create_polygon(str(polygon_dir / "example_start_region.csv")),
        finish_region=gp.create_polygon(str(polygon_dir / "example_finish_region.csv")),
    )

    assert result.match_flag == MATCH_WITHIN_THRESHOLD
    assert result.dtw == pytest.approx(0.100879, rel=DTW_RTOL)
    assert result.time.seconds == 1511
    assert result.start_point[0] == 47.150105
    assert result.start_point[1] == 9.150117
    assert result.end_point[0] == 47.16394
    assert result.end_point[1] == 9.176325


def test_radius_too_small_reports_why_rather_than_raising(gp, gpx_dir):
    """A failure to find candidates is a returned result, not an exception."""
    result = gp.dtw_match(
        str(gpx_dir / "tds_sunnestube_segment.gpx"),
        str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"),
        radius=0.001,
    )

    assert result.match_flag == COULD_NOT_BE_EVALUATED
    assert not result.is_success()
    assert result.error == "No trackpoints found near centroid"
    assert result.time is None
    assert result.dtw is None
    assert result.start_point is None
    assert result.end_point is None


def test_an_unrelated_activity_does_not_match(gp, gpx_dir):
    """A Zurich marathon segment against a Liechtenstein hill run."""
    result = gp.dtw_match(
        str(gpx_dir / "green_marathon_segment.gpx"),
        str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"),
        radius=20,
    )

    assert result.match_flag == COULD_NOT_BE_EVALUATED
    assert not result.is_success()
    assert result.error == "No trackpoints found near centroid"


@pytest.mark.slow
def test_min_trkps_can_rule_out_every_combination(gp, gpx_dir):
    """Candidates exist, but none are far enough apart, so nothing is tested.

    That is a different outcome from finding no candidates at all: flag -1 with
    no error, rather than flag -2 with one.
    """
    result = gp.dtw_match(
        str(gpx_dir / "tds_sunnestube_segment.gpx"),
        str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"),
        radius=7,
        min_trkps=10**9,
    )

    assert result.match_flag == NO_MATCH
    assert not result.is_success()
    assert result.error is None
    assert result.time is None


@pytest.mark.slow
def test_a_strict_threshold_downgrades_a_match_to_the_grey_zone(gp, gpx_dir):
    """Below ``dtw_threshold`` but within ``dtw_threshold * dtw_margin_range``.

    The shortest such candidate is kept and reported as flag 1 -- a match the
    caller is told to treat as provisional rather than one silently discarded.
    """
    result = gp.dtw_match(
        str(gpx_dir / "tds_sunnestube_segment.gpx"),
        str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"),
        radius=7,
        dtw_threshold=0.05,
        dtw_margin_range=4.0,
    )

    assert result.match_flag == MATCH_IN_GREY_ZONE
    assert result.is_success()
    assert 0.05 < result.dtw <= 0.2


@pytest.mark.slow
def test_a_loose_threshold_still_matches(gp, gpx_dir):
    result = gp.dtw_match(
        str(gpx_dir / "tds_sunnestube_segment.gpx"),
        str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"),
        radius=7,
        dtw_threshold=0.5,
    )

    assert result.match_flag == MATCH_WITHIN_THRESHOLD
    assert result.is_success()


@pytest.mark.slow
def test_the_matched_endpoints_come_from_the_activity(gp, gpx_dir):
    """Reported points must be real trackpoints, not interpolated ones."""
    activity = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"))
    result = gp.dtw_match(
        str(gpx_dir / "tds_sunnestube_segment.gpx"),
        str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"),
        radius=7,
    )

    latitudes = [float(value) for value in activity[0]]

    assert float(result.start_point[0]) in latitudes
    assert float(result.end_point[0]) in latitudes


@pytest.mark.slow
def test_the_elapsed_time_matches_the_two_endpoints(gp, gpx_dir):
    result = gp.dtw_match(
        str(gpx_dir / "tds_sunnestube_segment.gpx"),
        str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"),
        radius=7,
    )

    assert result.time == result.end_point[3] - result.start_point[3]


def test_a_corrupt_track_is_rejected_before_matching(gp, gpx_dir, tmp_path):
    """Missing latitude, longitude or time is fatal; missing elevation is not."""
    corrupt = tmp_path / "corrupt.gpx"
    corrupt.write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<gpx version="1.1" creator="test" '
        'xmlns="http://www.topografix.com/GPX/1/1">\n'
        " <trk><trkseg>\n"
        '  <trkpt lat="47.10" lon="9.10"><ele>1000.0</ele></trkpt>\n'
        '  <trkpt lat="47.20" lon="9.20"><ele>1001.0</ele></trkpt>\n'
        " </trkseg></trk>\n"
        "</gpx>\n",
        encoding="utf-8",
    )

    from geopard import GeopardException

    with pytest.raises(GeopardException, match="GPX file corrupted"):
        gp.dtw_match(str(gpx_dir / "tds_sunnestube_segment.gpx"), str(corrupt))
