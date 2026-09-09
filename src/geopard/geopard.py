"""Backwards-compatible import location.

``from geopard.geopard import Geopard`` was the only way to reach the class
before the package grew an ``__init__``, and it is what the shipped examples,
the tests and anyone's existing scripts use. It keeps working.

New code should prefer ``from geopard import Geopard``.
"""

from __future__ import annotations

from geopard.exceptions import GeopardException
from geopard.matching import Geopard
from geopard.response import GeopardResponse, MatchFlag

__all__ = ["Geopard", "GeopardException", "GeopardResponse", "MatchFlag"]
