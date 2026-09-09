"""Dynamic time warping between two tracks."""

from __future__ import annotations

import logging
from datetime import timedelta

import numpy as np
from scipy.spatial.distance import cdist

from geopard.interpolation import interpolate

logger = logging.getLogger(__name__)

DEFAULT_DISTANCE_METRIC = "euclidean"


def accumulated_cost_matrix(
    reference: np.ndarray,
    query: np.ndarray,
    distance_metric: str = DEFAULT_DISTANCE_METRIC,
) -> np.ndarray:
    """Accumulated cost matrix of two point sequences.

    Cell ``(n, m)`` holds the cost of the cheapest alignment of the first
    ``n + 1`` reference points with the first ``m + 1`` query points, so the
    bottom-right corner is the dynamic time warping distance between the two
    curves as a whole.

    Follows Müller, *Information Retrieval for Music and Motion*, Springer 2007,
    Theorem 4.3, https://doi.org/10.1007/978-3-540-74048-3.

    The recurrence is inherently sequential -- each cell reads the one to its
    left, which was computed in the same pass -- so this stays an explicit
    double loop. It is the slowest step in the library, and quadratic in
    :data:`~geopard.interpolation.INTERPOLATION_POINTS`.

    Args:
        reference: An ``(n, d)`` array of points.
        query: An ``(m, d)`` array of points.
        distance_metric: Any metric ``scipy.spatial.distance.cdist`` accepts.

    Returns:
        The ``(n, m)`` accumulated cost matrix.
    """
    logger.debug("Compute accumulated cost matrix")

    cost = cdist(reference, query, metric=distance_metric)
    length_n, length_m = cost.shape

    acm = np.zeros((length_n, length_m))
    acm[0, 0] = cost[0, 0]

    # First column and first row: with one point on the other side there is no
    # alignment to choose between, so each is a running sum.
    for n in range(1, length_n):
        acm[n, 0] = acm[n - 1, 0] + cost[n, 0]

    for m in range(1, length_m):
        acm[0, m] = acm[0, m - 1] + cost[0, m]

    # The interior: step diagonally, down, or right -- whichever is cheapest.
    for n in range(1, length_n):
        for m in range(1, length_m):
            acm[n, m] = cost[n, m] + min(acm[n - 1, m], acm[n, m - 1], acm[n - 1, m - 1])

    return acm


def dtw_computation(gpx_data: np.ndarray, gold: np.ndarray) -> tuple[float, timedelta]:
    """Compare one stretch of an activity against an interpolated gold segment.

    Args:
        gpx_data: A ``(4, n)`` slice of the activity, still at its own point
            count and still carrying timestamps.
        gold: The gold segment, already interpolated to ``(points, 2)``.

    Returns:
        The dynamic time warping distance, and how long the activity took to
        cover this stretch. The elapsed time is read from the activity, since
        a gold segment is often a drawn route with no clock at all.
    """
    logger.debug("Compute dynamic time warping between two gpx-track curves.")

    delta_time = gpx_data[3][-1] - gpx_data[3][0]
    distance = accumulated_cost_matrix(interpolate(gpx_data), gold)[-1, -1]

    return float(distance), delta_time
