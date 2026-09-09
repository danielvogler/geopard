"""Matching of GPX segments with dynamic time warping.

The names below are the public API::

    import geopard

    gp = geopard.Geopard()
    response = gp.dtw_match("segment.gpx", "activity.gpx")

    if response.is_success():
        print(response.time, response.dtw)

Each stage is also importable on its own -- :mod:`geopard.gpx`,
:mod:`geopard.selection`, :mod:`geopard.interpolation`, :mod:`geopard.dtw`,
:mod:`geopard.geodesy`, :mod:`geopard.regions` and :mod:`geopard.plotting`.
"""

from geopard.exceptions import GeopardException
from geopard.matching import Geopard
from geopard.response import GeopardResponse, MatchFlag
from geopard.version import __version__

__all__ = [
    "Geopard",
    "GeopardException",
    "GeopardResponse",
    "MatchFlag",
    "__version__",
]
