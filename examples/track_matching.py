"""Geopard example input file.

Example usage:

python ./examples/track_matching.py
    --GOLD_FILE_NAME=./data/gpx_files/tds_sunnestube_segment.gpx
    --ACTIVITY_FILE_NAME=./data/gpx_files/tds_sunnestube_segment_25_25.gpx
    --RADIUS=7
"""
import argparse
import logging
import sys

import matplotlib.pyplot as plt

from geopard.geopard import Geopard
from geopard.settings import PROJECT_ROOT

logging.basicConfig(encoding="utf-8", level=logging.INFO)

parser = argparse.ArgumentParser()
parser.add_argument("--GOLD_FILE_NAME", type=str, required=False)
parser.add_argument("--ACTIVITY_FILE_NAME", type=str, required=False)
parser.add_argument("--RADIUS", type=str, required=False)
parser.add_argument("--START_REGION_FILE_NAME", type=str, required=False)
parser.add_argument("--FINISH_REGION_FILE_NAME", type=str, required=False)
args, unknown = parser.parse_known_args()

# initialize
gp = Geopard()

# gold file, GPX track against which other tracks are compared
if args.GOLD_FILE_NAME:
    GOLD_FILE_NAME = args.GOLD_FILE_NAME
else:
    GOLD_FILE_NAME = (
        PROJECT_ROOT + "/data/gpx_files/tds_sunnestube_segment.gpx"
    )

# activity to evaluate
if args.ACTIVITY_FILE_NAME:
    ACTIVITY_FILE_NAME = args.ACTIVITY_FILE_NAME
else:
    ACTIVITY_FILE_NAME = (
        PROJECT_ROOT + "/data/gpx_files/tds_sunnestube_activity_25_25.gpx"
    )

# RADIUS (m) around start/end trackpoints
if args.RADIUS:
    RADIUS = args.RADIUS
else:
    RADIUS = 7

# start area geometry
# example:
# START_REGION_FILE_NAME = (
#        PROJECT_ROOT + "/data/csv_polygon_files/example_start_region.csv"
#    )
if args.START_REGION_FILE_NAME:
    START_REGION_FILE_NAME = args.START_REGION_FILE_NAME
else:
    START_REGION_FILE_NAME = None

# end area geometry:
# example:
# FINISH_REGION_FILE_NAME = (
#        PROJECT_ROOT + "/data/csv_polygon_files/example_finish_region.csv"
#    )
if args.FINISH_REGION_FILE_NAME:
    FINISH_REGION_FILE_NAME = args.FINISH_REGION_FILE_NAME
else:
    FINISH_REGION_FILE_NAME = None

logging.info("Track matching of example segments/activities")
geopard_response = gp.dtw_match(
    gold_name=GOLD_FILE_NAME, activity_name=ACTIVITY_FILE_NAME, radius=RADIUS
)

logging.info("Evaluate Geopard reponse")
if not geopard_response.is_success():
    logging.debug("\n----- Matching not successful -----")
    logging.debug("Error: %s", geopard_response.error)
    sys.exit(-1)

gp.parse_response(geopard_response)

logging.info("Plotting evaluated GPX tracks")
gp.plot_track_comparison(
    gold_file_name=GOLD_FILE_NAME,
    activity_file_name=ACTIVITY_FILE_NAME,
    radius=7,
)


if START_REGION_FILE_NAME or FINISH_REGION_FILE_NAME:
    logging.info("Plotting polygons")
    fig = plt.figure(
        num=None, figsize=(14, 10), dpi=80, facecolor="w", edgecolor="k"
    )
    FONT_SIZE = 30
    LW = 5
    plt.rcParams.update({"font.size": FONT_SIZE})

    if START_REGION_FILE_NAME:
        start_region_polygon = gp.create_polygon(
            file_name=START_REGION_FILE_NAME
        )
        plt.plot(
            *start_region_polygon.exterior.xy,
            c="k",
            label="Start",
            linewidth=LW,
        )

    if FINISH_REGION_FILE_NAME:
        finish_region_polygon = gp.create_polygon(
            file_name=FINISH_REGION_FILE_NAME
        )
        plt.plot(
            *finish_region_polygon.exterior.xy,
            c="b",
            label="Finish",
            linewidth=LW,
        )
    activity_gpx_track = gp.gpx_loading(ACTIVITY_FILE_NAME)
    gp.gpx_plot(fig, activity_gpx_track, ["Track", "o", "r"])
    plt.show()
