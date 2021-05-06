"""
Daniel Vogler
track_comparison

run as:
python track_comparison.py gold_file.gpx activity_file.gpx

"""

from geopard.geopard import Geopard
from matplotlib import pyplot as plt
import sys

### initialize
gp = Geopard()

gold_name = sys.argv[1]
activity_name = sys.argv[2]

### radius (m) around start/end trackpoints
radius = 20

"""
Track plotting
"""

### load gold standard/baseline segment
gold = gp.gpx_loading(gold_name)
### interpolate gold data
gold_interpolated = gp.interpolate(gold)

### load activity data to be edited
trkps = gp.gpx_loading(activity_name)
### crop activity data to segment length
gpx_cropped = gp.gpx_track_crop(gold, trkps, radius)

### plot gpx tracks
fig = plt.figure(num=None, figsize=(200, 150), dpi=80, facecolor='w', edgecolor='k')
gp.gpx_plot(fig,trkps,["Activity",".","k"])
gp.gpx_plot(fig,gold,["Gold",".","r"])

### plot interpolated gpx tracks
fig = plt.figure(num=None, figsize=(200, 150), dpi=80, facecolor='w', edgecolor='k')
gpx_interpolated = gp.interpolate(gpx_cropped)
gp.gpx_plot(fig,gpx_interpolated.T,["Activity Interpolated",".","k"])
gp.gpx_plot(fig,gold_interpolated.T,["Gold Interpolated",".","r"])

plt.show()
exit()