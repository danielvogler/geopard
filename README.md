# geopard
Matching of gpx segments with dynamic time warping.

[![python](https://img.shields.io/badge/Python-3.8-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/pre-commit/action/main.svg)](https://results.pre-commit.ci/latest/github/pre-commit/action/main)
[![linter](https://img.shields.io/badge/code%20linting-pylint-blue.svg)](https://github.com/PyCQA/pylint)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pydocstyle](https://img.shields.io/badge/pydocstyle-enabled-AD4CD3)](http://www.pydocstyle.org/en/stable/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Description

Pre-processing and analysis of gpx tracks (activities) for comparison to an existing gpx track (gold standard, segment). Checked are the joint start and end points with a given tolerance to trim the activity. Both tracks are interpolated to allow for both curves to be compared with dynamic time warping. Dynamic time warping allows to assess whether the activity actually completed the gold standard segment, and what the shortest required time was (in case of multiple repetitions or many points of the activity within the allowed distance of start and end points).

## Setup

### Poetry (recommended)

- Poetry is used for the virtual Python environment.
  - Install [python poetry](https://github.com/python-poetry/poetry):
    ```bash
    curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python
    ```
    or
    ```bash
    pip install --user poetry
    ```

  - Make sure that the virtual environment is installed in the root folder of this projects:
    ```bash
    poetry config virtualenvs.in-project true
    ```

  - Install dependencies:
    ```bash
    poetry install --no-root
    ```

  - Add packages:
    To add further packages, run:
    ```bash
    poetry add <package-name>
    ```

### pip

Alternatively, install python dependencies with pip:

`pip install -r requirements.txt`

### Code formatting
- Formatting via pre-commit hook:
  Python code is formatted using pre-commit hooks.
- Install pre-commit:
  ```bash
  pip install pre-commit
  ```
  or
  ```bash
  brew install pre-commit
  ```

- Install pre-commit hooks from the config file .pre-commit-config.yaml
  ```bash
  pre-commit install
  ```

- Run the pre-commit hooks:
  ```bash
  pre-commit run -a
  ```


## Usage

- Example GPX tracks are available at [https://github.com/danielvogler/geopard_tests](https://github.com/danielvogler/geopard_tests)
- Example files for construction of start/finish region in `./utils/`
- Example usage demonstrated in `./geopard_example`

### Example usage with circular start region around gold start/end points:

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


### Example usage with start/finish regions:

```
import geopard

### initialize
gp = geopard.Geopard()

start_region = gp.create_polygon("./utils/example_start_region.csv")
finish_region = gp.create_polygon("./utils/example_finish_region.csv")

gold_name = './gpx_files/tds_sunnestube_segment.gpx'
activity_name = './gpx_files/tds_sunnestube_activity_25_25.gpx'

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
```

## Examples

![Example image](/images/example_track.png "Example of gpx crop")
Example of gold segment, total activity and activity cropped to gold segment length.


![Example image](/images/example_track_start-finish.png "Example start and end points")
Example of start and end points to crop gpx tracks and obtain pairs for dtw matching.
