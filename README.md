# geopard
Matching of gpx segments with dynamic time warping.

Pre-processing and analysis of gpx tracks (activities) for comparison to an existing gpx track (gold standard, segment). Checked are the joint start and end points with a given tolerance to trim the activity. Both tracks are interpolated to allow for both curves to be compared with dynamic time warping. Dynamic time warping allows to assess whether the activity actually completed the gold standard segment, and what the shortest required time was (in case of multiple repetitions or many points of the activity within the allowed distance of start and end points). 

### Usage 
- Example GPX tracks are available at `https://github.com/danielvogler/gpx_processing/tree/main/gpx_files`
- Example usage demonstrated in `./geopard_example`
- Example usage:

```
import geopard

### initialize
gp = geopard.Geopard()

### dtw matching of example segments/activities
final_time, final_dtw = gp.dtw_match(gold_name, activity_name)

### optional with:
### min_trkps - minimum number of trackpoints between start and finish
### radius - radius around start and finish of gold segment
### dtw_threshold - segment match quality
final_time, final_dtw = gp.dtw_match(gold_name, activity_name, min_trkps = 100, radius=15, dtw_threshold=0.3)
```

### Dependencies

Install python dependencies with pip:

`pip install -r requirements.txt`

![Example image](/images/example_track.png "Example of gpx crop")
Example of gold segment, total activity and activity cropped to gold segment length.


![Example image](/images/example_track_start-finish.png "Example start and end points")
Example of start and end points to crop gpx tracks and obtain pairs for dtw matching.
