"""Compare two tracks.

track_comparison.py

run as:
python track_comparison.py --gold_file_name=gold_file.gpx
    --activity_file_name=activity_file.gpx

"""
import argparse
import logging

from matplotlib import pyplot as plt

from geopard.geopard import Geopard

# initialize
gp = Geopard()

parser = argparse.ArgumentParser()
parser.add_argument("--gold_file_name", type=str, required=False)
parser.add_argument("--activity_file_name", type=str, required=False)
parser.add_argument("--radius", type=str, required=False)
args, unknown = parser.parse_known_args()

if args.gold_file_name:
    GOLD_FILE_NAME = args.gold_file_name
else:
    GOLD_FILE_NAME = "tds_sunnestube_segment.gpx"

if args.activity_file_name:
    ACTIVITY_FILE_NAME = args.activity_file_name
else:
    ACTIVITY_FILE_NAME = "tds_sunnestube_activity_25_25.gpx"

# radius (m) around start/end trackpoints
if args.radius:
    RADIUS = args.radius
else:
    RADIUS = 20

# GOLD STANDARD
# load gold standard/baseline segment
gold_gpx_track = gp.gpx_loading(GOLD_FILE_NAME)
# interpolate gold data
gold_gpx_track_interpolated = gp.interpolate(gold_gpx_track)

# ACTIVITY
# load activity data to be edited
activity_gpx_track = gp.gpx_loading(ACTIVITY_FILE_NAME)
# crop activity data to segment length
activity_gpx_track_cropped = gp.gpx_track_crop(
    gold=gold_gpx_track, gpx_data=activity_gpx_track, radius=RADIUS
)

logging.info("Track plotting")

# PLOTTING
# plot gpx tracks
fig = plt.figure(
    num=None, figsize=(12, 8), dpi=80, facecolor="w", edgecolor="k"
)
gp.gpx_plot(fig, activity_gpx_track, ["Activity", ".", "k"])
gp.gpx_plot(fig, gold_gpx_track, ["Gold", ".", "r"])

# plot interpolated gpx tracks
fig = plt.figure(
    num=None, figsize=(12, 8), dpi=80, facecolor="w", edgecolor="k"
)
gpx_interpolated = gp.interpolate(activity_gpx_track_cropped)
gp.gpx_plot(fig, gpx_interpolated.T, ["Activity Interpolated", ".", "k"])
gp.gpx_plot(
    fig, gold_gpx_track_interpolated.T, ["Gold Interpolated", ".", "r"]
)

plt.show()
