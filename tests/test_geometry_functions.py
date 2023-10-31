"""Provide unit test cases."""
import logging
import unittest

from shapely.geometry.polygon import Polygon

from geopard.geopard import Geopard
from geopard.settings import PROJECT_ROOT

logging.basicConfig(encoding="utf-8", level=logging.INFO)


class TestGeometry(unittest.TestCase):
    """Test example for geopard."""

    def test_csv_polygon(self):
        """Test Polygon creation from CSV."""
        gp = Geopard()

        # create polygon from CSV file
        example_polygon_file = PROJECT_ROOT + "/tests/data/example_polygon.csv"
        example_polygon = gp.create_polygon(example_polygon_file)

        # polygon to test against
        coords = ((1.0, 3.0), (2.0, 4.0), (4.0, 2.0), (3.0, 1.0), (1.0, 3.0))
        test_polygon = Polygon(coords)

        self.assertEqual(example_polygon, test_polygon)

    def test_gpx_reading(self):
        """Test GPX track reading from *.gpx file."""
        gp = Geopard()

        # Read GPX track
        example_gpx_file = PROJECT_ROOT + "/tests/data/example_gpx_track.gpx"
        example_gpx = gp.gpx_loading(example_gpx_file)

        self.assertEqual(example_gpx[0][2], 47.13)
        self.assertEqual(example_gpx[1][1], 9.12)
        self.assertEqual(example_gpx[2][0], 1300.0)

    def test_spheroid_distance(self):
        """Test GPX track reading from *.gpx file."""
        gp = Geopard()

        # Example coordinates (NEBRASKA, USA, KANSAS, USA)
        coordinates_1 = [41.507483, -99.436554]
        coordinates_2 = [38.504048, -98.315949]

        # compute distance with haversine formula
        distance_in_m = gp.spheroid_point_distance(
            coordinates_1=coordinates_1, coordinates_2=coordinates_2
        )
        distance_in_km = round(distance_in_m / 1000, 1)

        self.assertEqual(distance_in_km, 347.3)
