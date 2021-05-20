from setuptools import setup
from geopard.version import __version__

setup(name='geopard',
      version=__version__,
      description='Matching of gpx segments with dynamic time warping',
      long_description='Matching of gpx segments with dynamic time warping. Pre-processing and analysis of gpx tracks (activities) for comparison to an existing gpx track (gold standard, segment). Checked are the joint start and end points with a given tolerance to trim the activity. Both tracks are interpolated to allow for both curves to be compared with dynamic time warping. Dynamic time warping allows to assess whether the activity actually completed the gold standard segment, and what the shortest required time was (in case of multiple repetitions or many points of the activity within the allowed distance of start and end points).',
      # long_description_content_type ='text/markdown',
      url='https://github.com/danielvogler/geopard',
      author='Daniel Vogler, Sebastian de Castelberg',
      author_email='geopard.py@gmail.com',
      license='MIT',
      packages=['geopard'],
      install_requires=['numpy>=1.11',
                        'matplotlib>=1.5',
                        'similaritymeasures>=0.4.4',
                        'haversine>=2.3.0',
                        'gpxpy>=1.4.2',
                        'scipy>=1.5.4',
                        'shapely>=1.7.1'])
