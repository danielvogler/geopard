"""Example to use polygons as start/end regions."""
import argparse
import logging

import matplotlib.pyplot as plt

from geopard.geopard import Geopard
from geopard.settings import PROJECT_ROOT

# initialize
gp = Geopard()

parser = argparse.ArgumentParser()
parser.add_argument("--start_region_file_name", type=str, required=False)
parser.add_argument("--finish_region_file_name", type=str, required=False)
parser.add_argument("--activity_name", type=str, required=False)
args, unknown = parser.parse_known_args()

if args.activity_name:
    ACTIVITY_NAME = args.activity_name
else:
    ACTIVITY_NAME = "tds_sunnestube_activity_25_25.gpx"

# start area geometry
if args.start_region_file_name:
    START_REGION_FILE_NAME = args.start_region_file_name
else:
    START_REGION_FILE_NAME = (
        PROJECT_ROOT + "/data/csv_polygon_files/example_start_region.csv"
    )

# end area geometry
if args.finish_region_file_name:
    FINISH_REGION_FILE_NAME = args.finish_region_file_name
else:
    FINISH_REGION_FILE_NAME = (
        PROJECT_ROOT + "/data/csv_polygon_files/example_finish_region.csv"
    )

logging.info("Track: %s", ACTIVITY_NAME)
logging.info("Start region: %s", START_REGION_FILE_NAME)
logging.info("Finish region: %s", FINISH_REGION_FILE_NAME)

# start and end area polygons
start_region_polygon = gp.create_polygon(START_REGION_FILE_NAME)
finish_region_polygon = gp.create_polygon(FINISH_REGION_FILE_NAME)

# load gold standard/baseline segment
track = gp.gpx_loading(file_name=ACTIVITY_NAME)

logging.info("Plotting")
fig = plt.figure(
    num=None, figsize=(12, 8), dpi=80, facecolor="w", edgecolor="k"
)
FONT_SIZE = 30
LW = 5
plt.rcParams.update({"font.size": FONT_SIZE})

# plot start/finish polygons
plt.plot(*start_region_polygon.exterior.xy, c="k", label="Start", linewidth=LW)
plt.plot(
    *finish_region_polygon.exterior.xy, c="b", label="Finish", linewidth=LW
)
gp.gpx_plot(fig, track, ["Track", "o", "r"])
plt.show()
