"""Finding candidate start/finish trackpoints, and cropping to them.

These two steps decide which slices of an activity are even considered for
matching, so their indices are the hinge the whole result turns on.
"""

from __future__ import annotations

import pytest

from geopard import GeopardException


@pytest.fixture(scope="session")
def gold(gp, gpx_dir):
    return gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))


@pytest.fixture(scope="session")
def activity(gp, gpx_dir):
    return gp.gpx_loading(str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"))


# --- nearest neighbours, circular radius --------------------------------


def test_finds_candidates_around_the_gold_start(gp, gold, activity):
    _, indices = gp.nearest_neighbours(activity, gold[:4, 0], 7.0)

    assert indices == [36, 41]


def test_finds_candidates_around_the_gold_finish(gp, gold, activity):
    _, indices = gp.nearest_neighbours(activity, gold[:4, -1], 7.0)

    assert indices == [1574, 1579]


def test_returns_the_edges_of_each_run_not_every_point(gp, gold, activity):
    """Consecutive candidates collapse to the first and last of each run.

    Without that, a runner standing at the start line contributes dozens of
    near-identical start points and the combination count explodes.
    """
    _, indices = gp.nearest_neighbours(activity, gold[:4, 0], 7.0)

    assert len(indices) < 10


def test_a_larger_radius_finds_at_least_as_many(gp, gold, activity):
    _, narrow = gp.nearest_neighbours(activity, gold[:4, 0], 7.0)
    _, wide = gp.nearest_neighbours(activity, gold[:4, 0], 50.0)

    assert len(wide) >= len(narrow)


def test_returned_points_line_up_with_the_returned_indices(gp, gold, activity):
    points, indices = gp.nearest_neighbours(activity, gold[:4, 0], 7.0)

    assert points.shape == (4, len(indices))
    for column, index in enumerate(indices):
        assert points[0, column] == activity[0, index]
        assert points[1, column] == activity[1, index]


def test_every_candidate_really_is_within_the_radius(gp, gold, activity):
    radius = 7.0
    _, indices = gp.nearest_neighbours(activity, gold[:4, 0], radius)

    for index in indices:
        distance = gp.spheroid_point_distance([activity[0, index], activity[1, index]], gold[:4, 0][:2])
        assert distance < radius


def test_raises_when_nothing_is_near_the_centroid(gp, activity):
    with pytest.raises(GeopardException, match="No trackpoints found near centroid"):
        gp.nearest_neighbours(activity, [0.0, 0.0], 7.0)


def test_raises_when_the_radius_is_vanishingly_small(gp, gold, activity):
    with pytest.raises(GeopardException):
        gp.nearest_neighbours(activity, gold[:4, 0], 1e-6)


# --- nearest neighbours, polygon region ---------------------------------


def test_a_region_polygon_replaces_the_radius(gp, gold, activity, polygon_dir):
    region = gp.create_polygon(str(polygon_dir / "example_start_region.csv"))

    _, indices = gp.nearest_neighbours(activity, gold[:4, 0], None, region)

    assert indices == [0, 47]


def test_the_finish_region_selects_a_different_stretch(gp, gold, activity, polygon_dir):
    region = gp.create_polygon(str(polygon_dir / "example_finish_region.csv"))

    _, indices = gp.nearest_neighbours(activity, gold[:4, -1], None, region)

    assert indices == [1557, 1615]


def test_the_centroid_is_ignored_once_a_region_is_given(gp, activity, polygon_dir):
    """A region is absolute, so a nonsense centroid must not change the answer."""
    region = gp.create_polygon(str(polygon_dir / "example_start_region.csv"))

    _, with_real = gp.nearest_neighbours(activity, [47.15, 9.15], None, region)
    _, with_nonsense = gp.nearest_neighbours(activity, [0.0, 0.0], None, region)

    assert with_real == with_nonsense


def test_a_region_far_from_the_track_raises(gp, activity, gp_far_region):
    with pytest.raises(GeopardException):
        gp.nearest_neighbours(activity, [0.0, 0.0], None, gp_far_region)


@pytest.fixture(scope="session")
def gp_far_region():
    from shapely.geometry.polygon import Polygon

    return Polygon([(0.0, 0.0), (0.0, 0.1), (0.1, 0.1), (0.1, 0.0)])


# --- cropping ------------------------------------------------------------


def test_crop_narrows_the_activity_to_the_matched_stretch(gp, gold, activity):
    cropped = gp.gpx_track_crop(gold, activity, 7.0)

    assert cropped.shape == (4, 1544)
    assert cropped.shape[1] < activity.shape[1]


def test_crop_keeps_the_four_rows(gp, gold, activity):
    cropped = gp.gpx_track_crop(gold, activity, 7.0)

    assert cropped.shape[0] == 4


def test_crop_spans_the_recorded_time_window(gp, gold, activity):
    cropped = gp.gpx_track_crop(gold, activity, 7.0)

    assert str(cropped[3][0]) == "2021-01-23 11:00:36+00:00"
    assert str(cropped[3][-1]) == "2021-01-23 11:26:20+00:00"


def test_crop_starts_no_earlier_and_ends_no_later_than_the_activity(gp, gold, activity):
    cropped = gp.gpx_track_crop(gold, activity, 7.0)

    assert cropped[3][0] >= activity[3][0]
    assert cropped[3][-1] <= activity[3][-1]


def test_crop_is_contiguous_in_time(gp, gold, activity):
    """It selects a window, not a scatter -- times must be non-decreasing."""
    cropped = gp.gpx_track_crop(gold, activity, 7.0)

    times = list(cropped[3])
    assert times == sorted(times)


def test_crop_with_regions_selects_a_different_window(gp, gold, activity, polygon_dir):
    start = gp.create_polygon(str(polygon_dir / "example_start_region.csv"))
    finish = gp.create_polygon(str(polygon_dir / "example_finish_region.csv"))

    cropped = gp.gpx_track_crop(gold, activity, None, start, finish)

    assert cropped.shape[0] == 4
    assert cropped.shape[1] > 0
    assert cropped[3][0] <= cropped[3][-1]


def test_crop_propagates_the_no_candidate_error(gp, gold, activity):
    with pytest.raises(GeopardException):
        gp.gpx_track_crop(gold, activity, 1e-6)


def test_neither_a_region_nor_a_radius_is_an_error(gp, activity):
    """Silently returning nothing here would look like "no candidates found"."""
    with pytest.raises(GeopardException, match="region polygon"):
        gp.nearest_neighbours(activity, [47.15, 9.15], None, None)


def test_a_centroid_without_a_radius_is_an_error(gp, activity):
    with pytest.raises(GeopardException, match="region polygon"):
        gp.nearest_neighbours(activity, centroid=[47.15, 9.15])
