# geopard
Matching of gpx segments with dynamic time warping.

Pre-processing and analysis of gpx tracks (activities) for comparison to an existing gpx track (gold standard, segment). Checked are the joint start and end points with a given tolerance to trim the activity. Both tracks are interpolated to allow for both curves to be compared with dynamic time warping. Dynamic time warping allows to assess whether the activity actually completed the gold standard segment, and what the shortest required time was (in case of multiple repetitions or many points of the activity within the allowed distance of start and end points). 

### Usage 
- Example GPX tracks are available at [https://github.com/danielvogler/geopard_tests](https://github.com/danielvogler/geopard_tests)
- Example files for construction of start/finish region in `./utils/`
- Example usage demonstrated in `./geopard_example`
- Example usage with circular start region around gold start/end points:

```
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
```


- Example usage with start/finish regions:

```
import geopard

start_region = "./utils/example_start_region.csv"
finish_region = "./utils/example_finish_region.csv"

### initialize
gp = geopard.Geopard()

geopard_response = gp.dtw_match(folder_path+gold_name, folder_path+activity_name, start_region=start_region, finish_region=finish_region)

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
```

### Dependencies

Install python dependencies with pip:

`pip install -r requirements.txt`

![Example image](/images/example_track.png "Example of gpx crop")
Example of gold segment, total activity and activity cropped to gold segment length.


![Example image](/images/example_track_start-finish.png "Example start and end points")
Example of start and end points to crop gpx tracks and obtain pairs for dtw matching.
