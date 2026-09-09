"""Great-circle distance on the WGS-84 mean sphere.

Values here are checked against the closed form rather than against a previous
run, because the haversine formula has exact answers at the degenerate cases
and those are where an implementation actually goes wrong.
"""

from __future__ import annotations

import pytest

# Mean Earth radius used by the implementation, https://doi.org/10.1007/s001900050278
EARTH_RADIUS_M = 6371008.7714
HALF_CIRCUMFERENCE_M = 3.141592653589793 * EARTH_RADIUS_M


def test_returns_zero_for_identical_points(gp):
    assert gp.spheroid_point_distance([47.1, 9.1], [47.1, 9.1]) == 0.0


def test_one_degree_of_latitude_is_about_111_km(gp):
    distance = gp.spheroid_point_distance([0.0, 0.0], [1.0, 0.0])

    assert distance == pytest.approx(111195.0797, abs=1e-3)


def test_antipodal_points_are_half_the_circumference(gp):
    distance = gp.spheroid_point_distance([0.0, 0.0], [0.0, 180.0])

    assert distance == pytest.approx(HALF_CIRCUMFERENCE_M, rel=1e-12)


def test_pole_to_pole_matches_antipodal(gp):
    """Different inputs, same true answer -- catches a lat/lon transposition."""
    poles = gp.spheroid_point_distance([90.0, 0.0], [-90.0, 0.0])
    equator = gp.spheroid_point_distance([0.0, 0.0], [0.0, 180.0])

    assert poles == pytest.approx(equator, rel=1e-12)


def test_known_distance_nebraska_to_kansas(gp):
    """The case the original suite pinned, kept so the benchmark is a superset."""
    distance = gp.spheroid_point_distance([41.507483, -99.436554], [38.504048, -98.315949])

    assert round(distance / 1000, 1) == 347.3
    assert distance == pytest.approx(347328.826230, abs=1e-4)


def test_distance_is_symmetric(gp):
    zurich = [47.376267, 8.535439]
    liechtenstein = [47.150031, 9.149962]

    there = gp.spheroid_point_distance(zurich, liechtenstein)
    back = gp.spheroid_point_distance(liechtenstein, zurich)

    assert there == pytest.approx(back, rel=1e-12)


def test_short_distance_is_metres_not_kilometres(gp):
    """Two adjacent trackpoints in Zurich are ~44 m apart, not ~0.044."""
    distance = gp.spheroid_point_distance([47.376267, 8.535439], [47.375947, 8.535789])

    assert distance == pytest.approx(44.2796118, abs=1e-4)


def test_ignores_trailing_elements_of_a_trackpoint(gp):
    """Callers pass full [lat, lon, ele, time] rows; only the first two count."""
    plain = gp.spheroid_point_distance([47.1, 9.1], [47.2, 9.2])
    with_extras = gp.spheroid_point_distance([47.1, 9.1, 1300.0, "t"], [47.2, 9.2, 1400.0, "t"])

    assert plain == with_extras
