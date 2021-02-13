from setuptools import setup
setup(name='geopard',
      version='0.1',
      description='Matching of gpx segments with dynamic time warping',
      url='https://github.com/geopard-py/geopard',
      author='Daniel Vogler',
      author_email='geopard.py@gmail.com',
      license='MIT',
      packages=['geopard'],
      install_requires=['numpy>=1.11',
                        'matplotlib>=1.5',
                        'similaritymeasures>=0.4.4',
                        'haversine>=2.3.0',
                        'gpxpy>=1.4.2',
                        'scipy>=1.5.4'])
