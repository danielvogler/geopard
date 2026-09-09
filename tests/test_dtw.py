"""The accumulated cost matrix and the distance read off its corner.

``acm`` is the only genuinely deterministic part of the pipeline -- no file
reading, no RNG -- so it is checked against hand-computed matrices rather than
against a recorded run.
"""

from __future__ import annotations

import numpy as np
import pytest

# A straight line and a line with one detour, small enough to compute by hand.
REFERENCE = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])
QUERY = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]])


def test_matrix_has_one_cell_per_pair_of_points(gp):
    matrix = gp.acm(REFERENCE, QUERY)

    assert matrix.shape == (len(REFERENCE), len(QUERY))


def test_first_cell_is_the_raw_distance_between_the_first_points(gp):
    matrix = gp.acm(REFERENCE, QUERY)

    assert matrix[0, 0] == 0.0


def test_matches_a_hand_computed_matrix(gp):
    """Müller (2007), Theorem 4.3 -- the recurrence this implements."""
    matrix = gp.acm(REFERENCE, QUERY, distance_metric="cityblock")

    expected = np.array(
        [
            [0.0, 2.0, 4.0],
            [1.0, 1.0, 2.0],
            [3.0, 3.0, 1.0],
            [6.0, 6.0, 2.0],
        ]
    )

    assert np.allclose(matrix, expected)


def test_euclidean_is_the_default_metric(gp):
    assert np.allclose(gp.acm(REFERENCE, QUERY), gp.acm(REFERENCE, QUERY, "euclidean"))


def test_identical_curves_cost_nothing(gp):
    matrix = gp.acm(REFERENCE, REFERENCE)

    assert matrix[-1, -1] == 0.0


def test_first_row_and_column_accumulate_without_choice(gp):
    """With one point on a side there is no path to choose; it is a running sum."""
    matrix = gp.acm(REFERENCE, QUERY, distance_metric="cityblock")

    assert list(matrix[:, 0]) == [0.0, 1.0, 3.0, 6.0]
    assert list(matrix[0, :]) == [0.0, 2.0, 4.0]


def test_every_cost_is_non_negative(gp):
    matrix = gp.acm(REFERENCE, QUERY)

    assert np.all(matrix >= 0.0)


def test_distance_is_symmetric_in_its_two_curves(gp):
    """The recurrence treats both sequences alike, so swapping them changes nothing."""
    forward = gp.acm(REFERENCE, QUERY)[-1, -1]
    backward = gp.acm(QUERY, REFERENCE)[-1, -1]

    assert forward == pytest.approx(backward)


def test_scaling_both_curves_scales_the_distance(gp):
    """Euclidean cost is homogeneous, so the metric carries no hidden constant."""
    plain = gp.acm(REFERENCE, QUERY)[-1, -1]
    scaled = gp.acm(REFERENCE * 3.0, QUERY * 3.0)[-1, -1]

    assert scaled == pytest.approx(3.0 * plain)


def test_a_detour_costs_more_than_a_straight_line(gp):
    straight = gp.acm(REFERENCE, np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]))
    detour = gp.acm(REFERENCE, QUERY)

    assert detour[-1, -1] > straight[-1, -1]


def test_dtw_computation_returns_distance_and_elapsed_time(gp, gpx_dir):
    gold = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))
    activity = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"))
    cropped = gp.gpx_track_crop(gold, activity, 7.0)

    distance, elapsed = gp.dtw_computation(cropped, gp.interpolate(gold))

    assert distance == pytest.approx(0.096, rel=5e-2)
    assert elapsed.total_seconds() == 1544.0


def test_elapsed_time_is_measured_from_the_activity_not_the_gold(gp, gpx_dir):
    """The gold segment here carries no timestamps at all, so it cannot be."""
    gold = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))
    activity = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"))
    cropped = gp.gpx_track_crop(gold, activity, 7.0)

    _, elapsed = gp.dtw_computation(cropped, gp.interpolate(gold))

    assert elapsed == cropped[3][-1] - cropped[3][0]
