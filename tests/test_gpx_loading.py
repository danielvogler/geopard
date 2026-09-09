"""Reading a GPX file into the [lat, lon, ele, time] array the rest uses."""

from __future__ import annotations

import numpy as np
import pytest

from geopard import GeopardException

DUPLICATE_TRACK = """<?xml version="1.0" encoding="UTF-8"?>
<gpx version="1.1" creator="test" xmlns="http://www.topografix.com/GPX/1/1">
 <trk><trkseg>
  <trkpt lat="47.10" lon="9.10"><ele>1000.0</ele></trkpt>
  <trkpt lat="47.20" lon="9.20"><ele>1001.0</ele></trkpt>
  <trkpt lat="47.20" lon="9.20"><ele>1002.0</ele></trkpt>
  <trkpt lat="47.20" lon="9.20"><ele>1003.0</ele></trkpt>
  <trkpt lat="47.30" lon="9.30"><ele>1004.0</ele></trkpt>
 </trkseg></trk>
</gpx>
"""

MULTI_SEGMENT_TRACK = """<?xml version="1.0" encoding="UTF-8"?>
<gpx version="1.1" creator="test" xmlns="http://www.topografix.com/GPX/1/1">
 <trk>
  <trkseg><trkpt lat="47.10" lon="9.10"/><trkpt lat="47.11" lon="9.11"/></trkseg>
  <trkseg><trkpt lat="47.12" lon="9.12"/></trkseg>
 </trk>
 <trk>
  <trkseg><trkpt lat="47.13" lon="9.13"/></trkseg>
 </trk>
</gpx>
"""


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return str(path)


def test_returns_four_rows_of_lat_lon_ele_time(gp, fixture_dir):
    track = gp.gpx_loading(str(fixture_dir / "example_gpx_track.gpx"))

    assert track.shape == (4, 3)


def test_reads_the_documented_values(gp, fixture_dir):
    """The assertions the original suite made, kept verbatim."""
    track = gp.gpx_loading(str(fixture_dir / "example_gpx_track.gpx"))

    assert track[0][2] == 47.13
    assert track[1][1] == 9.12
    assert track[2][0] == 1300.0


def test_rows_are_latitude_longitude_elevation_in_that_order(gp, fixture_dir):
    track = gp.gpx_loading(str(fixture_dir / "example_gpx_track.gpx"))

    assert list(track[0]) == [47.11, 47.12, 47.13]
    assert list(track[1]) == [9.11, 9.12, 9.13]
    assert list(track[2]) == [1300.0, 1400.0, 1500.0]


def test_drops_a_trackpoint_repeated_three_times(gp, tmp_path):
    """A stationary GPS fix repeats a coordinate; the spline cannot take it.

    Only the middle of a run of three identical points is removed, so the run
    above collapses from five points to four.
    """
    track = gp.gpx_loading(_write(tmp_path, "dup.gpx", DUPLICATE_TRACK))

    assert track.shape == (4, 4)
    assert list(track[0]) == [47.10, 47.20, 47.20, 47.30]
    # the dropped point takes its own elevation with it
    assert list(track[2]) == [1000.0, 1001.0, 1003.0, 1004.0]


def test_concatenates_every_segment_of_every_track(gp, tmp_path):
    track = gp.gpx_loading(_write(tmp_path, "multi.gpx", MULTI_SEGMENT_TRACK))

    assert list(track[0]) == [47.10, 47.11, 47.12, 47.13]


def test_missing_elevation_is_kept_as_none(gp, tmp_path):
    """Elevation is optional. Latitude, longitude and time are not."""
    track = gp.gpx_loading(_write(tmp_path, "multi.gpx", MULTI_SEGMENT_TRACK))

    assert all(value is None for value in track[2])


@pytest.mark.parametrize(
    "name, points, first_lat, last_lon",
    [
        ("tds_sunnestube_segment.gpx", 61, 47.14998, 9.17663),
        ("tds_sunnestube_activity_25_25.gpx", 1616, 47.14974, 9.176462),
        ("green_marathon_segment.gpx", 4802, 47.376189, 8.53579),
        ("green_marathon_activity_4_15_17.gpx", 15190, 47.376394, 8.535869),
    ],
)
def test_shipped_tracks_load_to_a_stable_shape(gp, gpx_dir, name, points, first_lat, last_lon):
    """Pins the deduplication result on the real tracks.

    A change to duplicate detection shifts every downstream index, so this is
    the cheapest place to catch one.
    """
    track = gp.gpx_loading(str(gpx_dir / name))

    assert track.shape == (4, points)
    assert track[0][0] == pytest.approx(first_lat, abs=1e-9)
    assert track[1][-1] == pytest.approx(last_lon, abs=1e-9)


def test_times_are_timezone_aware_and_increasing(gp, gpx_dir):
    track = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_activity_25_25.gpx"))

    assert track[3][0].tzinfo is not None
    assert track[3][-1] > track[3][0]


def test_a_gold_segment_may_carry_no_timestamps_at_all(gp, gpx_dir):
    """``tds_sunnestube_segment.gpx`` is a drawn route, not a recorded one.

    Only the activity's clock is ever read, so a gold segment without one is
    matched normally. Rejecting it would break the shipped example.
    """
    track = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))

    assert all(value is None for value in track[3])


def test_track_is_a_two_dimensional_numpy_array(gp, gpx_dir):
    """Downstream code slices it as ``track[:4, 0]``, which needs an ndarray."""
    track = gp.gpx_loading(str(gpx_dir / "tds_sunnestube_segment.gpx"))

    assert isinstance(track, np.ndarray)
    assert track.ndim == 2


def test_missing_file_raises(gp, tmp_path):
    with pytest.raises(FileNotFoundError):
        gp.gpx_loading(str(tmp_path / "does-not-exist.gpx"))


def test_geopard_exception_is_an_exception():
    """It is caught by name in dtw_match; it must stay catchable as one."""
    assert issubclass(GeopardException, Exception)
