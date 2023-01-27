"""Geopard example."""
import logging

from matplotlib import pyplot as plt

from geopard.geopard import Geopard
from geopard.settings import PROJECT_ROOT

### initialize
gp = Geopard()

folder_path = PROJECT_ROOT + "/gpx_files/"

logging.info("Example activities")

### example - one-way skimo
### dtw=0.09702, radius=7m, t=0:25:22
gold_name = "tds_sunnestube_segment.gpx"
activity_name = "tds_sunnestube_activity_25_25.gpx"  # 0:25:22
start_region = gp.create_polygon("../utils/example_start_region.csv")
finish_region = gp.create_polygon("../utils/example_finish_region.csv")

logging.info("Plotting data.")

radius = 7

### start and end area polygons
start_region_polygon = gp.create_polygon(start_region)
finish_region_polygon = gp.create_polygon(finish_region)

### load gold standard/baseline segment
gold = gp.gpx_loading(folder_path + gold_name)

### load activity data to be edited
trkps = gp.gpx_loading(folder_path + activity_name)
### crop activity data to segment length
gpx_cropped = gp.gpx_track_crop(gold, trkps, start_region=start_region, finish_region=finish_region)

### find potential start/end trackpoints - just for plotting
nn_start, nn_start_idx = gp.nearest_neighbours(gpx_cropped, region=start_region)
nn_finish, nn_finish_idx = gp.nearest_neighbours(gpx_cropped, region=finish_region)

logging.info("Track matching")

### dtw matching of example segments/activities
geopard_response = gp.dtw_match(
    folder_path + gold_name, folder_path + activity_name, start_region=start_region, finish_region=finish_region
)

if not geopard_response.is_success():
    logging.warning("\n----- Matching not successful -----")
    logging.warning("Error:", geopard_response.error)
    exit(-1)

final_time = geopard_response.time
final_dtw = geopard_response.dtw
match_flag = geopard_response.match_flag
final_start_point = geopard_response.start_point
final_end_point = geopard_response.end_point

logging.info("Track plotting")

### plot gpx tracks
fig = plt.figure(num=None, figsize=(200, 150), dpi=80, facecolor="w", edgecolor="k")
### plot start/finish polygons
plt.plot(*start_region_polygon.exterior.xy, c="b", label="Start region", linewidth=5)
plt.plot(*finish_region_polygon.exterior.xy, c="b", label="Finish region", linewidth=5, linestyle="dotted")
gp.gpx_plot(fig, nn_start, ["NN Start Cropped", "X", "b"], 1200)
gp.gpx_plot(fig, nn_finish, ["NN Finish Cropped", "P", "b"], 1200)
gp.gpx_plot(fig, trkps, ["Activity", ".", "k"])
gp.gpx_plot(fig, gpx_cropped, ["Activity Cropped", "o", "k"])
gp.gpx_plot(fig, gold, ["Gold", "o", "r"])

plt.show()
exit()
