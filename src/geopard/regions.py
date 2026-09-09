"""Start and finish regions, read from CSV into polygons."""

from __future__ import annotations

import logging
from csv import DictReader
from pathlib import Path

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

logger = logging.getLogger(__name__)


def create_polygon(file_name: str) -> Polygon:
    """Build a region polygon from a CSV of latitude/longitude corners.

    The file needs ``Latitude`` and ``Longitude`` columns; any others (an
    ``ID``, say) are ignored. Rows are taken in file order as the polygon's
    corners, and the ring is closed for you.

    Note the axis order. The CSV reads latitude-first because that is how a
    person writes down a coordinate, while Shapely is ``(x, y)`` and therefore
    longitude-first. This function is where the two conventions meet, and the
    swap happens here so nothing downstream has to think about it.

    Args:
        file_name: Path to the CSV.

    Returns:
        The region, as a Shapely polygon in ``(longitude, latitude)``.
    """
    logger.debug("Create polygon from %s", file_name)

    with Path(file_name).open(newline="", encoding="utf-8") as csv_file:
        corners = [Point(float(row["Longitude"]), float(row["Latitude"])) for row in DictReader(csv_file)]

    return Polygon([(point.x, point.y) for point in corners])
