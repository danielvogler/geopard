"""Setup file."""
from pathlib import Path

from setuptools import setup

from geopard.version import __version__

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="geopard",
    version=__version__,
    description="Matching of gpx segments with dynamic time warping",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/danielvogler/geopard",
    author="Daniel Vogler, Sebastian de Castelberg",
    author_email="geopard.py@gmail.com",
    license="MIT",
    packages=["geopard"],
    install_requires=[
        "numpy>=1.11",
        "matplotlib>=1.5",
        "gpxpy>=1.4.2",
        "scipy>=1.5.4",
        "shapely>=1.7.1",
    ],
)
