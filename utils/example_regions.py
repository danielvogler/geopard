### daniel vogler
### example to create polygon

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import csv
import sys
import matplotlib.pyplot as plt
from geopard.geopard import Geopard

### initialize
gp = Geopard()

### start area geometry
try:
    file_start_region = str(sys.argv[1])
except: 
    file_start_region = 'example_start_region.csv'

### end area geometry
try:
    file_finish_region = str(sys.argv[2])
except: 
    file_finish_region = 'example_finish_region.csv'

### track example
### example - one-way skimo
### dtw=0.09702, radius=7m, t=0:25:22
try: 
    gpx_track = str(sys.argv[3])
except:
    gpx_track = "../gpx_files/tds_sunnestube_segment.gpx"

print('Track:',gpx_track)
print('Start region:',file_start_region)
print('Finish region:',file_finish_region)

def create_polygon(file_name):
    region = []
    ### read in lat/lon and conver to points
    with open(file_name, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            region.append(Point([float(row['Longitude']), float(row['Latitude'])]))

    ### convert points to polygon
    poly = Polygon([[p.x, p.y] for p in region])

    return poly

### start and end area polygons
start_region_polygon = create_polygon(file_start_region)
finish_region_polygon = create_polygon(file_finish_region)

### load gold standard/baseline segment
track = gp.gpx_loading(gpx_track)

###
### PLOT
###
fig = plt.figure(num=None, figsize=(200, 150), dpi=80, facecolor='w', edgecolor='k')
font_size = 30
lw=5
plt.rcParams.update({'font.size': font_size})

### plot start/finish polygons
plt.plot(*start_region_polygon.exterior.xy,c='k',label='Start',linewidth=lw)
plt.plot(*finish_region_polygon.exterior.xy,c='b',label='Finish',linewidth=lw)

gp.gpx_plot(fig,track,["Track","o","r"])

plt.show()