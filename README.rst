geopard
=======

Matching of gpx segments with dynamic time warping.

Pre-processing and analysis of gpx tracks (activities) for comparison to an existing gpx track (gold standard, segment). Checked are the joint start and end points with a given tolerance to trim the activity. Both tracks are interpolated to allow for both curves to be compared with dynamic time warping. Dynamic time warping allows to assess whether the activity actually completed the gold standard segment, and what the shortest required time was (in case of multiple repetitions or many points of the activity within the allowed distance of start and end points).

Usage
=====
- Example GPX tracks are available at [https://github.com/danielvogler/geopard_tests](https://github.com/danielvogler/geopard_tests)
- Example files for construction of start/finish region in `<PROJECT_ROOT>/data/csv_polygon_files/`
- Example usage demonstrated in `<PROJECT_ROOT>/geopard_example`
- Example usage with circular start region around gold start/end points:

.. code:: python

  import geopard

  ### initialize
  gp = geopard.Geopard()

  ### dtw matching of example segments/activities
  response = gp.dtw_match(gold_name, activity_name)

  ### optional with:
  ### min_trkps - minimum number of trackpoints between start and finish
  ### radius - radius around start and finish of gold segment
  ### dtw_threshold - segment match quality
  response = gp.dtw_match(gold_name, activity_name, min_trkps = 100, radius=15, dtw_threshold=0.3)

  ### GeopardResponse

  # final time
  response.time

  # final dtw
  response.dtw

  # final start point
  response.start_point

  # final end point
  response.end_point

  # match flag
  response.match_flag

  # is_success and error
  response.is_success
  response.error


- Example usage with start/finish regions:

.. code:: python

  import geopard

  start_region = "<PROJECT_ROOT>/data/csv_polygon_files/example_start_region.csv"
  finish_region = "<PROJECT_ROOT>/data/csv_polygon_files/example_finish_region.csv"

  gold_name = '<PROJECT_ROOT>/data/gpx_files/tds_sunnestube_segment.gpx'
  activity_name = '<PROJECT_ROOT>/data/gpx_files/tds_sunnestube_activity_25_25.gpx'

  ### initialize
  gp = geopard.Geopard()

  response = gp.dtw_match(gold_name, activity_name, start_region=start_region, finish_region=finish_region)

  ### GeopardResponse

  # final time
  response.time

  # final dtw
  response.dtw

  # final start point
  response.start_point

  # final end point
  response.end_point

  # match flag
  response.match_flag

  # is_success and error
  response.is_success
  response.error


Dependencies
============

Install python dependencies with pip:

::

   pip install -r requirements.txt


Images
======

.. figure:: https://raw.githubusercontent.com/danielvogler/geopard/master/docs/images/example_track.png
  :alt: Example of gold segment, total activity and activity cropped to gold segment length.
Example of gold segment, total activity and activity cropped to gold segment length.

.. figure:: https://raw.githubusercontent.com/danielvogler/geopard/master/docs/images/example_track_start-finish.png
  :alt: Example of start and end points to crop gpx tracks and obtain pairs for dtw matching.
Example of start and end points to crop gpx tracks and obtain pairs for dtw matching.
