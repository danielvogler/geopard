"""Distance between two points on the Earth."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from math import asin, cos, radians, sin, sqrt

import numpy as np

logger = logging.getLogger(__name__)

#: Mean Earth radius in metres, https://doi.org/10.1007/s001900050278
MEAN_EARTH_RADIUS_M = 6371008.7714

#: A point, as anything indexable by 0 and 1. Callers pass Python lists; the
#: matcher passes a column of a track, which is a numpy array. Narrowing this
#: to ``Sequence[float]`` would misdescribe half the call sites.
Coordinates = Sequence[float] | np.ndarray


def spheroid_point_distance(coordinates_1: Coordinates, coordinates_2: Coordinates) -> float:
    """Great-circle distance between two points, in metres.

    Uses the `haversine formula
    <https://en.wikipedia.org/wiki/Haversine_formula>`_, which stays accurate
    for the short distances this library cares about -- the radius around a
    segment's start is metres, where the naive spherical law of cosines loses
    precision badly.

    Args:
        coordinates_1: ``[latitude, longitude]`` in degrees. Further elements
            are ignored, so a whole ``[lat, lon, ele, time]`` trackpoint can be
            passed straight in.
        coordinates_2: The other point, same shape.

    Returns:
        Distance in metres.
    """
    logger.debug("Compute distance between two points on spheroid.")

    lat_1, lon_1 = radians(coordinates_1[0]), radians(coordinates_1[1])
    lat_2, lon_2 = radians(coordinates_2[0]), radians(coordinates_2[1])

    delta_lat = (lat_2 - lat_1) / 2
    delta_lon = (lon_2 - lon_1) / 2

    haversine = sin(delta_lat) ** 2 + cos(lat_1) * cos(lat_2) * sin(delta_lon) ** 2

    return 2 * MEAN_EARTH_RADIUS_M * asin(sqrt(haversine))
