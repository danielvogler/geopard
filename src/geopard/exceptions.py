"""Errors raised by geopard."""

from __future__ import annotations


class GeopardException(Exception):
    """A track could not be read, or holds nothing usable for a match.

    Raised by the loading and selection steps. ``Geopard.dtw_match`` catches it
    around the selection stage and turns it into a
    :class:`~geopard.response.GeopardResponse` carrying
    :data:`~geopard.response.MatchFlag.NOT_EVALUATED`, so a caller matching many
    activities never has to wrap the call in a ``try``.
    """
