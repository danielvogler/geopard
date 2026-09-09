"""What a match returns."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import timedelta
from enum import IntEnum
from typing import Any

logger = logging.getLogger(__name__)


class MatchFlag(IntEnum):
    """How a match ended.

    An ``IntEnum`` rather than a plain enum because these values have always
    been plain integers on the wire, and callers compare them as such. The
    sign is the contract: anything above zero is a match.
    """

    #: Distance is at or below ``dtw_threshold``.
    MATCH = 2
    #: Above the threshold but within ``dtw_threshold * dtw_margin_range``.
    #: The shortest such candidate is kept, and the caller decides.
    GREY_ZONE = 1
    #: Candidate start/finish pairs existed; none matched closely enough.
    NO_MATCH = -1
    #: The match never ran -- no usable trackpoints near start or finish.
    #: ``error`` says why.
    NOT_EVALUATED = -2


@dataclass(frozen=True)
class GeopardResponse:
    """The outcome of one segment match.

    Frozen: a response describes a match that already happened, and nothing
    downstream has any business editing it.

    Every field except ``match_flag`` is ``None`` when no match was found, so
    check :meth:`is_success` before reading them.
    """

    #: Elapsed time between the matched start and finish trackpoints.
    time: timedelta | None
    #: Dynamic time warping distance. Lower is a closer match.
    dtw: float | None
    #: The matched start trackpoint, as ``[lat, lon, ele, time]``.
    start_point: Any | None
    #: The matched finish trackpoint, as ``[lat, lon, ele, time]``.
    end_point: Any | None
    #: One of :class:`MatchFlag`.
    match_flag: int
    #: Why the match could not be evaluated, when it could not be.
    error: str | None = None

    def is_success(self) -> bool:
        """Whether a match was found, in the threshold or in the grey zone.

        A method rather than a property because it has always been one, and the
        shipped examples call it.
        """
        return self.match_flag > 0


def log_response(response: GeopardResponse) -> None:
    """Write a response to the log at INFO level.

    Tolerates the unmatched case, where every field but the flag is ``None``.
    """
    logger.info("Response time: %s", response.time)
    logger.info(
        "Response DTW: %s",
        "n/a" if response.dtw is None else f"{response.dtw:.3f}",
    )
    logger.info("Response start_point: %s", response.start_point)
    logger.info("Response end_point: %s", response.end_point)
    logger.info("Response match_flag: %s", response.match_flag)
    logger.info("Response success: %s", response.is_success())
    logger.info("Response error: %s", response.error)
