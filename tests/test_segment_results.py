"""Provide unit test cases."""
import logging
import unittest

from geopard.geopard import Geopard
from geopard.settings import PROJECT_ROOT

logging.basicConfig(encoding="utf-8", level=logging.INFO)


class TestGeopard(unittest.TestCase):
    """Test example for geopard."""

    def test_green_marathon(self):
        """Green marathon test segment."""
        gold_name = PROJECT_ROOT + "/data/gpx_files/green_marathon_segment.gpx"
        activity_name = (
            PROJECT_ROOT
            + "/data/gpx_files/green_marathon_activity_4_15_17.gpx"
        )

        gp = Geopard()
        res = gp.dtw_match(gold_name, activity_name, radius=20)

        self.assertEqual(res.match_flag, 2)
        self.assertAlmostEqual(res.dtw, 0.134, places=2)
        self.assertEqual(res.time.seconds, 15303)
        self.assertEqual(res.start_point[0], 47.376267)
        self.assertEqual(res.start_point[1], 8.535439)
        self.assertEqual(res.end_point[0], 47.375947)
        self.assertEqual(res.end_point[1], 8.535789)

    def test_sunnestube(self):
        """Sunnestube test segment."""
        gold_name = PROJECT_ROOT + "/data/gpx_files/tds_sunnestube_segment.gpx"
        activity_name = (
            PROJECT_ROOT + "/data/gpx_files/tds_sunnestube_activity_25_25.gpx"
        )

        gp = Geopard()
        res = gp.dtw_match(gold_name, activity_name, radius=7)

        self.assertEqual(res.match_flag, 2)
        self.assertAlmostEqual(res.dtw, 0.096, places=2)
        self.assertEqual(res.time.seconds, 1534)
        self.assertEqual(res.start_point[0], 47.150031)
        self.assertEqual(res.start_point[1], 9.149962)
        self.assertEqual(res.end_point[0], 47.164231)
        self.assertEqual(res.end_point[1], 9.176556)


if __name__ == "__main__":
    unittest.main()
