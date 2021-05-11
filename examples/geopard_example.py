"""
Daniel Vogler
geopard example
"""

from geopard.geopard import Geopard
from matplotlib import pyplot as plt

folder_path = "../gpx_files/"

"""
Example activities
"""

### example - cropping
### dtw = 0.28447 with radius = 40m
# gold_name = "tdh1_dv.gpx"
# activity_name = "tdh1_mg.gpx"

### example - one-way skimo
### dtw=0.09702, radius=7m, t=0:25:22
gold_name = "tds_sunnestube_segment.gpx"
activity_name = "tds_sunnestube_activity_25_25.gpx"          # 0:25:22
# activity_name = "tds_sunnestube_activity_25_55.gpx"        # 0:25:55
# activity_name = "tds_sunnestube_activity_25_39.gpx"        # 0:25:39
# gold_name = "strava.segments.25881647.TdU2_-Chez-le-Coiffeur.gpx"
# activity_name = "strava.activities.1108969469.Uetzgi-.gpx"
#gold_name = "strava-segments-26694751.tds1_sunnestube-605459c1ce1de130091189.gpx"
# activity_name = "4786250737-wixoox-6053d46bbc2ca546371668.gpx"
#activity_name = "4780044736-7ruftv-6053d4ee33ca3008909763.gpx"

### example - two loops
### dtw=0.06490, radius=34m, t=0:25:46
# gold_name = "nordicstar_weltcup_segment.gpx"
# activity_name = "nordicstar_weltcup_activity_25_52.gpx"    # 0:25:35

### example - one-way cross-country ski
### dtw=0.09671, radius=4m, t=0:44:26
# gold_name = "nordicstar_dischmatal_segment.gpx"
# activity_name = "nordicstar_dischmatal_activity_44_39.gpx" # 0:44:37

### example - green marathon zurich
### dtw = 0.12866, radius = 15m, time = 4:15:11
# gold_name = "green_marathon_segment.gpx"
# activity_name = "green_marathon_activity_4_15_17.gpx"        # 4:15:03

### example activity - no matching start/end points found
# gold_name = "tdh1_dv.gpx"
# activity_name = "tdu2a.gpx"

### example activity - gpx track jump during activity - tdh2
### dtw=0.08136, radius=7m, t=0:36:14
# gold_name = "tdh2.gpx"
# activity_name = "tdh2_error.gpx"

### example - intersecting tracks - tdu3
### dtw=0.01154, radius=7m, t=0:26:53
# gold_name = "tdu3_dv.gpx"
# activity_name = "tdu3_ls.gpx"

### example - match TH
# gold_name = "th1_gold.gpx"
# gold_name = "th2_gold.gpx"
# gold_name = "th3_gold.gpx"
# activity_name = "th1_ttb.gpx"

### radius (m) around start/end trackpoints
radius = 7

"""
Track matching
"""

### initialize
gp = Geopard()

### dtw matching of example segments/activities
geopard_response = gp.dtw_match(folder_path+gold_name, folder_path+activity_name,radius=radius)

if not geopard_response.is_success():
    print("\n----- Matching not successful -----")
    print("Error:" , geopard_response.error)
    exit(-1)

final_time = geopard_response.time
final_dtw = geopard_response.dtw
match_flag = geopard_response.match_flag
final_start_point = geopard_response.start_point
final_end_point = geopard_response.end_point

"""
Track plotting
"""

### load gold standard/baseline segment
gold = gp.gpx_loading(folder_path + gold_name)
### interpolate gold data
gold_interpolated = gp.interpolate(gold)

### load activity data to be edited
trkps = gp.gpx_loading(folder_path + activity_name)
### crop activity data to segment length
gpx_cropped = gp.gpx_track_crop(gold, trkps, radius)

### find potential start/end trackpoints - just for plotting
nn_start, nn_start_idx = gp.nearest_neighbours(gpx_cropped,gold[:4,0],radius)
nn_finish, nn_finish_idx = gp.nearest_neighbours(gpx_cropped,gold[:4,-1],radius)

### plot gpx tracks
fig = plt.figure(num=None, figsize=(200, 150), dpi=80, facecolor='w', edgecolor='k')
gp.gpx_plot(fig,nn_start,["NN Start Cropped","X","b"],1200)
gp.gpx_plot(fig,nn_finish,["NN Finish Cropped","P","b"],1200)
gp.gpx_plot(fig,trkps,["Activity",".","k"])
gp.gpx_plot(fig,gpx_cropped,["Activity Cropped","o","k"])
gp.gpx_plot(fig,gold,["Gold","o","r"])

### plot interpolated gpx tracks
fig = plt.figure(num=None, figsize=(200, 150), dpi=80, facecolor='w', edgecolor='k')
gp.gpx_plot(fig,gpx_cropped,["Activity Cropped","o","k"])
gpx_interpolated = gp.interpolate(gpx_cropped)
gp.gpx_plot(fig,gpx_interpolated.T,["Activity Interpolated",".","k"])
gp.gpx_plot(fig,gold,["Gold","o","r"])
gp.gpx_plot(fig,gold_interpolated.T,["Gold Interpolated",".","r"])

plt.show()
exit()
