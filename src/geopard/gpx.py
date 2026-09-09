"""Reading a GPX file into the array the rest of the library works on."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from gpxpy import parse

from geopard.exceptions import GeopardException

logger = logging.getLogger(__name__)

#: Row offsets into a loaded track. ``track[LATITUDE, i]`` is the i-th latitude.
LATITUDE, LONGITUDE, ELEVATION, TIME = 0, 1, 2, 3

#: Rows a track cannot do without. Elevation is allowed to be missing --
#: plenty of devices do not record it, and nothing here uses it.
REQUIRED_ROWS = (LATITUDE, LONGITUDE, TIME)

#: How many corrupt indices to name in the error message. Enough to see the
#: pattern, few enough to stay readable.
CORRUPT_INDICES_REPORTED = 3


def load_track(file_name: str) -> np.ndarray:
    """Read a GPX file into a ``(4, n)`` array of lat, lon, elevation, time.

    Every track and every segment in the file is concatenated in document
    order, and stationary repeats are dropped (see
    :func:`find_stationary_trackpoints`).

    Args:
        file_name: Path to a ``.gpx`` file.

    Returns:
        A ``(4, n)`` object array. Rows are latitude, longitude, elevation and
        timestamp, in that order.
    """
    logger.info("Load gpx file: %s", file_name)

    with Path(file_name).open(encoding="utf8") as gpx_file:
        gpx_data = parse(gpx_file)

    latitude, longitude, elevation, time = extract_trackpoints(gpx_data)
    stationary = find_stationary_trackpoints(latitude, longitude)

    logger.info("\tRemoving %s redundant trackpoints", len(stationary))

    return np.asarray(
        [
            np.delete(latitude, stationary, 0),
            np.delete(longitude, stationary, 0),
            np.delete(elevation, stationary, 0),
            np.delete(time, stationary, 0),
        ]
    )


def extract_trackpoints(
    gpx_data: object,
) -> tuple[list, list, list, list]:
    """Flatten every segment of every track into four parallel lists."""
    logger.info("\tExtract relevant gpx data")

    latitude: list[float] = []
    longitude: list[float] = []
    elevation: list[float | None] = []
    time: list[object] = []

    for track in gpx_data.tracks:  # type: ignore[attr-defined]
        for segment in track.segments:
            for point in segment.points:
                latitude.append(point.latitude)
                longitude.append(point.longitude)
                elevation.append(point.elevation)
                time.append(point.time)

    return latitude, longitude, elevation, time


def find_stationary_trackpoints(latitude: list[float], longitude: list[float]) -> list[int]:
    """Indices of trackpoints identical to both of their neighbours.

    A GPS device left standing emits the same fix over and over. The spline in
    :mod:`geopard.interpolation` cannot be fitted through duplicated points, so
    the middle of each such run is removed while its ends are kept -- which is
    what preserves the shape of the track either side of the pause.

    Args:
        latitude: Latitudes in track order.
        longitude: Longitudes in track order, same length.

    Returns:
        Sorted indices to drop.
    """
    logger.info("\tChecking for duplicate lat/lon trackpoints")

    def repeated(values: list[float]) -> set[int]:
        return {
            i for i in range(1, len(values) - 1) if values[i] == values[i - 1] and values[i] == values[i + 1]
        }

    # Only a point that repeats in *both* coordinates is stationary. A track
    # running due north repeats its longitude for hundreds of points and must
    # not be thinned.
    return sorted(repeated(latitude) & repeated(longitude))


def validate_track(track: np.ndarray) -> None:
    """Raise if latitude, longitude or time is missing anywhere.

    Missing elevation is accepted; nothing in the matching pipeline reads it.

    Args:
        track: A ``(4, n)`` array from :func:`load_track`.

    Raises:
        GeopardException: If any required row holds a ``None``.
    """
    for row in REQUIRED_ROWS:
        missing = [i for i in range(len(track[row])) if track[row, i] is None]

        if missing:
            first = sorted(set(missing))[:CORRUPT_INDICES_REPORTED]
            raise GeopardException(f"GPX file corrupted. Check GPX points (E.g. {first})")
