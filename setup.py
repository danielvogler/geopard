from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='geopard',
      version='0.1.10',
      description='Matching of gpx segments with dynamic time warping',
      long_description=readme(),
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
