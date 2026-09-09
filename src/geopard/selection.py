"""Choosing which stretches of an activity are worth matching.

An activity is usually much longer than the segment inside it, and may cover
that segment more than once. These two steps narrow it down: find the
trackpoints near the segment's start and finish, then throw away everything
outside the time window they span.
"""

from __future__ import annotations

import logging
from itertools import pairwise

import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from geopard.exceptions import GeopardException
from geopard.geodesy import Coordinates, spheroid_point_distance

logger = logging.getLogger(__name__)


def nearest_neighbours(
    gpx_data: np.ndarray,
    centroid: Coordinates | None = None,
    radius: float | None = None,
    region: Polygon | None = None,
) -> tuple[np.ndarray, list[int]]:
    """Find the trackpoints that count as being at a given place.

    Two ways to say where that place is. Pass a ``region`` polygon and a
    trackpoint qualifies by falling inside it; otherwise a trackpoint qualifies
    by lying within ``radius`` metres of ``centroid``. A region takes
    precedence, and ``centroid`` is then ignored entirely.

    Consecutive qualifying trackpoints collapse to the first and last of each
    run. Someone waiting at a start line produces dozens of near-identical
    fixes, and every one of them would otherwise become a candidate start --
    multiplying the number of combinations to warp without adding a single
    distinct one.

    Args:
        gpx_data: A ``(4, n)`` track.
        centroid: ``[latitude, longitude]`` to measure from. Ignored when
            ``region`` is given.
        radius: Radius in metres. Ignored when ``region`` is given.
        region: A polygon in ``(longitude, latitude)``.

    Returns:
        The qualifying trackpoints as a ``(4, k)`` array, and their indices
        into ``gpx_data``.

    Raises:
        GeopardException: If nothing qualifies.
    """
    logger.info("Finding nearest neighbours ...")

    if region is not None:
        indices = within_region(gpx_data, region)
    else:
        indices = within_radius(gpx_data, centroid, radius)

    indices = run_boundaries(indices)

    if not indices:
        raise GeopardException("No trackpoints found near centroid")

    neighbours = np.asarray([[gpx_data[row, i] for i in indices] for row in range(4)])

    return neighbours, indices


def within_region(gpx_data: np.ndarray, region: Polygon) -> list[int]:
    """Indices of trackpoints inside the polygon.

    The polygon is ``(longitude, latitude)`` while the track is latitude-first,
    so the two are swapped here.
    """
    indices = [i for i in range(gpx_data.shape[1]) if region.contains(Point(gpx_data[1, i], gpx_data[0, i]))]

    logger.info("\tNearest neighbors: %s", len(indices))
    logger.info("\t\twithin region polygon: %s", region)

    return indices


def within_radius(gpx_data: np.ndarray, centroid: Coordinates | None, radius: float | None) -> list[int]:
    """Indices of trackpoints within ``radius`` metres of ``centroid``."""
    if centroid is None or radius is None:
        raise GeopardException(
            "nearest_neighbours needs either a region polygon, or both a centroid and a radius"
        )

    indices = [
        i
        for i in range(gpx_data.shape[1])
        if spheroid_point_distance((gpx_data[0, i], gpx_data[1, i]), centroid[:2]) < radius
    ]

    logger.info("\tNearest neighbors: %s", len(indices))
    logger.info("\tRadius: %.1fm", radius)
    logger.info("\tCentroid: %s", list(centroid[:2]))

    return indices


def run_boundaries(indices: list[int]) -> list[int]:
    """Reduce each run of consecutive indices to its first and last."""
    ordered = sorted(set(indices))
    gaps = [edge for a, b in pairwise(ordered) if a + 1 < b for edge in (a, b)]

    return ordered[:1] + gaps + ordered[-1:]


def crop_track(
    gold: np.ndarray,
    gpx_data: np.ndarray,
    radius: float | None = None,
    start_region: Polygon | None = None,
    finish_region: Polygon | None = None,
) -> np.ndarray:
    """Trim an activity to the time window that could contain the segment.

    The window runs from the earliest trackpoint near the segment's start to
    the latest one near its finish, which keeps every repetition of the segment
    while dropping the commute either side of them.

    Args:
        gold: The gold segment, ``(4, n)``. Only its first and last
            trackpoints are read.
        gpx_data: The activity to crop, ``(4, m)``.
        radius: Radius in metres, when not using regions.
        start_region: Polygon around the segment start.
        finish_region: Polygon around the segment finish.

    Returns:
        The cropped activity, ``(4, k)``.

    Raises:
        GeopardException: If no trackpoint is near the start, or none near the
            finish.
    """
    logger.info("Cropping GPX track.")

    near_start, _ = nearest_neighbours(gpx_data, gold[:4, 0], radius, start_region)
    near_finish, _ = nearest_neighbours(gpx_data, gold[:4, -1], radius, finish_region)

    earliest = min(near_start[3])
    latest = max(near_finish[3])

    inside = [i for i in range(len(gpx_data[3])) if earliest <= gpx_data[3, i] <= latest]

    return np.asarray([[gpx_data[row, i] for i in inside] for row in range(4)])
