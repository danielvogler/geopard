"""Matching an activity against a gold-standard segment."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from datetime import datetime, timedelta
from typing import Any

import numpy as np
from shapely.geometry.polygon import Polygon

from geopard import dtw as dtw_module
from geopard import gpx, interpolation, plotting, regions, selection
from geopard.exceptions import GeopardException
from geopard.geodesy import Coordinates, spheroid_point_distance
from geopard.response import GeopardResponse, MatchFlag, log_response

logger = logging.getLogger(__name__)

#: Minimum trackpoints between a candidate start and finish. Guards against
#: matching a start against a finish a few seconds later, which happens
#: whenever a segment's start and finish regions overlap.
DEFAULT_MIN_TRACKPOINTS = 50

#: Radius in metres around the gold start/finish within which an activity
#: trackpoint counts as being there.
DEFAULT_RADIUS_M = 7.0

#: Warping distance at or below which a candidate is a match.
DEFAULT_DTW_THRESHOLD = 0.2

#: Multiplier on the threshold marking the grey zone. A candidate between the
#: threshold and ``threshold * margin`` is reported, flagged as provisional.
DEFAULT_DTW_MARGIN_RANGE = 1.5


class Geopard:
    """Entry point for segment matching.

    Holds no state, so one instance can be shared freely::

        from geopard import Geopard

        gp = Geopard()
        response = gp.dtw_match("segment.gpx", "activity.gpx")

    Every method here delegates to a module that can also be imported on its
    own -- :mod:`geopard.gpx`, :mod:`geopard.dtw`, :mod:`geopard.selection` and
    the rest. The class exists so that the whole library is reachable from one
    object, which is how it has always been used.
    """

    # --- matching --------------------------------------------------------

    def dtw_match(
        self,
        gold_name: str,
        activity_name: str,
        min_trkps: int = DEFAULT_MIN_TRACKPOINTS,
        radius: float = DEFAULT_RADIUS_M,
        dtw_threshold: float = DEFAULT_DTW_THRESHOLD,
        dtw_margin_range: float = DEFAULT_DTW_MARGIN_RANGE,
        start_region: Polygon | None = None,
        finish_region: Polygon | None = None,
    ) -> GeopardResponse:
        """Match an activity against a gold-standard segment.

        Finds every plausible pairing of a start trackpoint with a later finish
        trackpoint, then warps them against the segment shortest-first and stops
        at the first one close enough to count. Shortest-first is what makes the
        answer meaningful when an activity covers the segment several times: the
        result is the best time actually run on it, not merely one of them.

        Args:
            gold_name: Path to the gold-standard segment GPX file.
            activity_name: Path to the activity GPX file to evaluate.
            min_trkps: Minimum trackpoints between a candidate start and finish.
            radius: Radius in metres around the gold start and finish. Ignored
                where a region is given.
            dtw_threshold: Warping distance at or below which a candidate is a
                match. Lower is stricter.
            dtw_margin_range: Multiplier on the threshold defining the grey
                zone.
            start_region: Polygon to use instead of a radius around the start.
            finish_region: Polygon to use instead of a radius around the finish.

        Returns:
            A :class:`~geopard.response.GeopardResponse`. Check
            :meth:`~geopard.response.GeopardResponse.is_success` before reading
            the times and points.

        Raises:
            GeopardException: If the activity is missing latitude, longitude or
                time. Not finding a match is reported in the response instead.
        """
        logger.info("Starting DTW match")
        started_at = datetime.now()

        gold = gpx.load_track(gold_name)
        gold_interpolated = interpolation.interpolate(gold)

        activity = gpx.load_track(activity_name)
        gpx.validate_track(activity)

        try:
            cropped = selection.crop_track(gold, activity, radius, start_region, finish_region)
            _, start_indices = selection.nearest_neighbours(cropped, gold[:4, 0], radius, start_region)
            _, finish_indices = selection.nearest_neighbours(cropped, gold[:4, -1], radius, finish_region)
        except GeopardException as error:
            # Nothing near the start, or nothing near the finish. That is an
            # ordinary outcome when matching an activity against a segment it
            # never went near, so it is reported rather than raised -- a caller
            # sweeping a season of activities should not need a try block.
            return GeopardResponse(None, None, None, None, MatchFlag.NOT_EVALUATED, str(error))

        candidates = candidate_segments(cropped, start_indices, finish_indices, min_trkps)

        return self._warp_candidates(
            cropped,
            gold_interpolated,
            candidates,
            dtw_threshold,
            dtw_threshold * dtw_margin_range,
            started_at,
        )

    def _warp_candidates(
        self,
        cropped: np.ndarray,
        gold_interpolated: np.ndarray,
        candidates: list[tuple[timedelta, int, int]],
        dtw_threshold: float,
        grey_zone_threshold: float,
        started_at: datetime,
    ) -> GeopardResponse:
        """Warp candidates shortest-first, stopping at the first real match.

        Candidates arrive sorted by elapsed time, so the first one to clear
        ``dtw_threshold`` is also the fastest one that does, and the search can
        stop there. A candidate that only clears ``grey_zone_threshold`` is
        remembered but does not stop the search -- a later, slower candidate may
        still be a proper match, and a proper match always wins over a
        provisional one.
        """
        best = GeopardResponse(None, None, None, None, MatchFlag.NO_MATCH)
        best_dtw = float("inf")
        last_dtw: float | None = None
        tested = 0

        for _elapsed, start_idx, finish_idx in candidates:
            if best_dtw <= dtw_threshold:
                break

            tested += 1
            distance, delta_time = dtw_module.dtw_computation(
                cropped[:, start_idx : finish_idx + 1], gold_interpolated
            )
            last_dtw = distance

            logger.info(
                "DTW (y): %2.5f / T [s]: %s (%s/%s)",
                distance,
                delta_time,
                tested - 1,
                len(candidates),
            )

            if distance <= dtw_threshold:
                flag = MatchFlag.MATCH
            elif distance <= grey_zone_threshold and best.match_flag < 0:
                flag = MatchFlag.GREY_ZONE
            else:
                continue

            best_dtw = distance
            best = GeopardResponse(
                delta_time,
                distance,
                cropped[:, start_idx],
                cropped[:, finish_idx],
                flag,
            )

        log_summary(started_at, len(candidates), tested, best)

        if best.is_success():
            return best

        logger.warning("No segment match found")

        return GeopardResponse(None, last_dtw, None, None, MatchFlag.NO_MATCH)

    # --- the pipeline, individually ---------------------------------------

    def gpx_loading(self, file_name: str) -> np.ndarray:
        """Load a GPX file. See :func:`geopard.gpx.load_track`."""
        return gpx.load_track(file_name)

    def gpx_track_crop(
        self,
        gold: np.ndarray,
        gpx_data: np.ndarray,
        radius: float | None = None,
        start_region: Polygon | None = None,
        finish_region: Polygon | None = None,
    ) -> np.ndarray:
        """Crop an activity. See :func:`geopard.selection.crop_track`."""
        return selection.crop_track(gold, gpx_data, radius, start_region, finish_region)

    def nearest_neighbours(
        self,
        gpx_data: np.ndarray,
        centroid: Coordinates | None = None,
        radius: float | None = None,
        region: Polygon | None = None,
    ) -> tuple[np.ndarray, list[int]]:
        """Find candidate trackpoints.

        See :func:`geopard.selection.nearest_neighbours`.
        """
        return selection.nearest_neighbours(gpx_data, centroid, radius, region)

    def spheroid_point_distance(self, coordinates_1: Coordinates, coordinates_2: Coordinates) -> float:
        """Distance in metres. See :func:`geopard.geodesy.spheroid_point_distance`."""
        return spheroid_point_distance(coordinates_1, coordinates_2)

    def acm(
        self,
        reference: np.ndarray,
        query: np.ndarray,
        distance_metric: str = dtw_module.DEFAULT_DISTANCE_METRIC,
    ) -> np.ndarray:
        """Accumulated cost matrix.

        See :func:`geopard.dtw.accumulated_cost_matrix`.
        """
        return dtw_module.accumulated_cost_matrix(reference, query, distance_metric)

    def dtw_computation(self, gpx_data: np.ndarray, gold: np.ndarray) -> tuple[float, timedelta]:
        """Warp one stretch against the gold. See :func:`geopard.dtw.dtw_computation`."""
        return dtw_module.dtw_computation(gpx_data, gold)

    def interpolate(self, gpx_data: np.ndarray) -> np.ndarray:
        """Resample a track. See :func:`geopard.interpolation.interpolate`."""
        return interpolation.interpolate(gpx_data)

    def create_polygon(self, file_name: str) -> Polygon:
        """Read a region from CSV. See :func:`geopard.regions.create_polygon`."""
        return regions.create_polygon(file_name)

    # --- plotting ---------------------------------------------------------

    def plot_track_comparison(
        self,
        gold_file_name: str,
        activity_file_name: str,
        radius: float = DEFAULT_RADIUS_M,
        show: bool = True,
        save_to: str | None = None,
    ) -> Any:
        """Draw a gold segment against an activity.

        See :func:`geopard.plotting.plot_track_comparison`.
        """
        return plotting.plot_track_comparison(
            gold_file_name, activity_file_name, radius, show=show, save_to=save_to
        )

    def gpx_plot(
        self,
        fig: Any,
        gpx_data: np.ndarray,
        plot_info: Sequence[str],
        marker_size: int | None = None,
    ) -> Any:
        """Scatter a track onto a figure. See :func:`geopard.plotting.gpx_plot`."""
        return plotting.gpx_plot(fig, gpx_data, plot_info, marker_size)

    # --- reporting --------------------------------------------------------

    def parse_response(self, response: GeopardResponse) -> None:
        """Log a response. See :func:`geopard.response.log_response`."""
        log_response(response)


def candidate_segments(
    cropped: np.ndarray,
    start_indices: list[int],
    finish_indices: list[int],
    min_trkps: int,
) -> list[tuple[timedelta, int, int]]:
    """Every start/finish pairing worth warping, fastest first.

    A pairing qualifies when the finish comes after the start with at least
    ``min_trkps`` trackpoints in between.
    """
    candidates = [
        (cropped[3, finish] - cropped[3, start], start, finish)
        for start in start_indices
        for finish in finish_indices
        if start < finish - min_trkps
    ]

    return sorted(candidates, key=lambda candidate: candidate[0])


def log_summary(
    started_at: datetime,
    total: int,
    tested: int,
    best: GeopardResponse,
) -> None:
    """Write the end-of-match summary at INFO level."""
    elapsed = datetime.now() - started_at

    logger.info("----- Finished DTW segment match -----")
    logger.info("Total combinations to test: %s", total)
    logger.info("Total combinations tested: %s", tested)
    logger.info("Total execution time: %s", elapsed)

    if tested:
        logger.info("Execution time per combination: %s", elapsed / tested)

    logger.info("Match flag [-]: %s", int(best.match_flag))

    if best.is_success():
        logger.info("Final T [s]: %s", best.time)
        logger.info("Final DTW (y): %2.5f", best.dtw)
