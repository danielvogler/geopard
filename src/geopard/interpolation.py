"""Resampling a track onto a fixed number of points.

Dynamic time warping compares two curves point by point, so both have to be
the same length before it can run. Interpolating each track onto the same
count is what makes a 61-point drawn route comparable to a 1600-point
recording.
"""

from __future__ import annotations

import logging
from collections.abc import Callable

import numpy as np
from scipy.interpolate import splev, splprep

logger = logging.getLogger(__name__)

#: Points every track is resampled to. Both curves must agree on this, so it
#: is a module constant rather than an argument.
INTERPOLATION_POINTS = 1000

#: Half-width of the jitter added to each coordinate, in degrees. About a
#: centimetre -- far below GPS precision, and far above the spacing the spline
#: solver needs to tell two points apart.
COORDINATE_JITTER_DEG = 1e-7


def interpolate(
    gpx_data: np.ndarray,
    points: int = INTERPOLATION_POINTS,
    seed: int | None = None,
) -> np.ndarray:
    """Fit a spline through a track and resample it at equal intervals.

    Coordinates are jittered by up to :data:`COORDINATE_JITTER_DEG` first.
    ``splprep`` cannot fit through repeated points, and even after stationary
    trackpoints are dropped a track can still contain two consecutive fixes
    equal to the last decimal.

    That jitter is the one stochastic step in the library, and it is what makes
    a warping distance move in the fourth decimal between runs. Pass ``seed``
    to pin it, or seed numpy's global RNG with ``np.random.seed(...)`` to pin a
    whole pipeline. The seed is logged either way.

    Args:
        gpx_data: A track whose first two rows are latitude and longitude.
        points: How many points to return.
        seed: Seed for the jitter. ``None`` draws from numpy's global RNG,
            which is the default so that seeding globally still works.

    Returns:
        An ``(points, 2)`` array of latitude, longitude.
    """
    logger.debug("Interpolate gpx data along track (seed=%s).", seed)

    jitter = jitter_source(seed)

    jittered = np.array(
        [(lat + jitter(), lon + jitter()) for lat, lon in zip(gpx_data[0][:], gpx_data[1][:], strict=True)]
    )

    tck, u = splprep(jittered.T, u=None, s=0.0, per=0)

    latitude, longitude = splev(np.linspace(u.min(), u.max(), points), tck, der=0)

    return np.column_stack((latitude, longitude))


def jitter_source(seed: int | None) -> Callable[[], float]:
    """A function drawing one jitter value per call.

    With no seed this stays on numpy's global RNG, so the draw order -- and
    therefore the result under ``np.random.seed(...)`` -- is exactly what it
    has always been.
    """
    if seed is None:
        # NPY002 is suppressed deliberately. The legacy global RNG *is* this
        # branch: swapping in a Generator would silently stop
        # `np.random.seed(...)` from pinning a run, which is how the tests and
        # every existing caller make one reproducible.
        legacy_uniform = np.random.uniform

        return lambda: float(legacy_uniform(-COORDINATE_JITTER_DEG, COORDINATE_JITTER_DEG))

    generator = np.random.default_rng(seed)

    return lambda: float(generator.uniform(-COORDINATE_JITTER_DEG, COORDINATE_JITTER_DEG))
