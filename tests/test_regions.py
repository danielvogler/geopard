"""Start/finish regions read from CSV into Shapely polygons.

The column order is the trap: the CSV is written ``Latitude,Longitude`` for a
human to read, and Shapely wants ``(x, y)`` -- longitude first. Every test here
exists to keep that swap in place.
"""

from __future__ import annotations

import pytest
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


def test_builds_the_documented_polygon(gp, fixture_dir):
    """The assertion the original suite made, kept verbatim."""
    polygon = gp.create_polygon(str(fixture_dir / "example_polygon.csv"))

    expected = Polygon(((1.0, 3.0), (2.0, 4.0), (4.0, 2.0), (3.0, 1.0), (1.0, 3.0)))

    assert polygon == expected


def test_longitude_becomes_x_and_latitude_becomes_y(gp, fixture_dir):
    """The CSV row ``1,3.0,1.0`` is lat 3.0, lon 1.0 -- so the point is (1, 3)."""
    polygon = gp.create_polygon(str(fixture_dir / "example_polygon.csv"))

    assert polygon.exterior.coords[0] == (1.0, 3.0)


def test_ring_is_closed(gp, fixture_dir):
    polygon = gp.create_polygon(str(fixture_dir / "example_polygon.csv"))

    assert polygon.exterior.coords[0] == polygon.exterior.coords[-1]


def test_polygon_is_valid_and_has_area(gp, fixture_dir):
    polygon = gp.create_polygon(str(fixture_dir / "example_polygon.csv"))

    assert polygon.is_valid
    assert polygon.area == pytest.approx(4.0)


@pytest.mark.parametrize(
    "name, wkt",
    [
        (
            "example_start_region.csv",
            "POLYGON ((9.1485 47.15, 9.1509 47.1486, 9.1514 47.1494, 9.149 47.15075, 9.1485 47.15))",
        ),
        (
            "example_finish_region.csv",
            "POLYGON ((9.1743 47.1648, 9.175 47.1655, 9.1792 47.1639, 9.1785 47.163, 9.1743 47.1648))",
        ),
    ],
)
def test_shipped_regions_have_stable_geometry(gp, polygon_dir, name, wkt):
    polygon = gp.create_polygon(str(polygon_dir / name))

    assert polygon.wkt == wkt


def test_shipped_start_region_contains_the_start_of_the_activity(gp, polygon_dir, gpx_dir):
    """The regions and the tracks have to agree, or the example stops working."""
    region = gp.create_polygon(str(polygon_dir / "example_start_region.csv"))
    track = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"))

    inside = [i for i in range(track.shape[1]) if region.contains(Point(track[1, i], track[0, i]))]

    assert inside, "no trackpoint of the shipped activity falls in the start region"


def test_missing_file_raises(gp, tmp_path):
    with pytest.raises(FileNotFoundError):
        gp.create_polygon(str(tmp_path / "nope.csv"))


def test_requires_latitude_and_longitude_columns(gp, tmp_path):
    path = tmp_path / "wrong-headers.csv"
    path.write_text("ID,lat,lon\n1,47.0,9.0\n", encoding="utf-8")

    with pytest.raises(KeyError):
        gp.create_polygon(str(path))
